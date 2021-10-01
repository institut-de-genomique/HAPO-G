#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"

#include "alipile.h"
#include "polished.h"
#include "hash.h"

const unsigned int MAX_COV = 50;
const unsigned int MIN_COV = 3;
const unsigned int BUFFER = 1000;

void usage();
int parse_bam(char*, char*, char*, char*, int);
void print_step(int, int);

int main(int argc, char* argv[]) {
  char *BAMfile = NULL, *outfile = NULL, *changefile = NULL, *FAfile = NULL;
  int silent = 0;
  int c;

  // Invokes member function `int operator ()(void);'
  while ((c = getopt(argc, argv, "f:b:c:o:hs")) != -1) {
    switch (c) {
    case 'b':
      BAMfile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'c':
      changefile = optarg;
      break;
    case 'f':
      FAfile = optarg;
      break;
      /*    case 'u':
      unmasked_len = atoi(optarg);
      break;
    case 'l':
      level = atoi(optarg);
      break;*/
    case 'h':
      usage();
    case 's':
      silent = 1;
      break;
    default :
      abort();
    }
  }
  
  if (BAMfile==NULL) error("Could not load bam: %s\n", BAMfile);
  if (FAfile==NULL) error("Could not load faidx: %s\n", FAfile);
  if (outfile==NULL) error("Could not open output fasta file: %s\n", outfile);
  if (changefile==NULL) error("Could not open output changes file: %s\n", changefile);
  return parse_bam(BAMfile, FAfile, outfile, changefile, silent);
}

int parse_bam(char* bam, char* fa, char* outfa, char *changefile, int silent) {
  setvbuf(stdout, NULL, _IONBF, 0);
  
  //open BAM for reading
  samFile *in = sam_open(bam, "r");
  if(in == NULL) error("Unable to open BAM/SAM file: %s\n", bam);
  
  // open fasta file and index
  faidx_t *ref = fai_load(fa);
  if (ref == NULL) error("Could not load faidx: %s\n", fa);
  
  // open output fasta file
  FILE* out = fopen(outfa,"w");
  if (out == NULL) error("Could not open output fasta file: %s\n", outfa);
  
  // open output changes file
  FILE* changes = fopen(changefile,"w");
  if (changes == NULL) error("Could not open output changes file: %s\n", changefile);
  
  //Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  //Initiate the alignment record
  bam1_t *aln = bam_init1();
  bam1_t *aln2 = bam_init1();
  int ret=0, i=0;
  int current_position = 0, previous_position = 0;
  alipile_t* ap = alipile_init(changes);
  alipile_t* allali = alipile_init(NULL);
  polished_t* s = polished_init();
  hash_t* readname = hash_init();
  int tot_read = 0, added_read = 0, too_short = 0;
  int current_seq = -1, progress = 0, nb_changes = 0, nb_hapB = 0;
  
  while ((ret = sam_read1(in, header, aln)) >= 0) { 
    bam_copy1(aln2, aln);
    // exclude unmapped reads
    if (aln->core.tid < 0 || (aln->core.flag&BAM_FUNMAP)) continue;
    //printf("Read: %s\t%s\t%i\n", bam_get_qname(aln), get_bamseq(aln), aln->core.l_qseq);
    
    if(current_seq != aln->core.tid) {
      progress = 0;
      if(current_seq != -1) {
	if(!silent) {
	  print_step(ap->len_seq, ap->len_seq);
	  printf("\n");
	}
	int i = 0;
	for(i = previous_position ; i < ap->len_seq ; i++) 
	  select_base(ap, allali, i, s, readname);
	print_seq(s, out);
	nb_changes += ap->nb_changes;
	nb_hapB += ap->nbhapB;
	alipile_free(ap);
	alipile_free(allali);
	polished_free(s);
	hash_free(readname);
	previous_position = 0;
	ap = alipile_init(changes);
	allali = alipile_init(NULL);
	s = polished_init();
	readname = hash_init();
      }
      if(!silent) printf("Reference sequence : %s\n", header->target_name[aln->core.tid]);
      current_seq = aln->core.tid;
      ap->name_seq = header->target_name[aln->core.tid];
      s->name_seq = header->target_name[aln->core.tid];
      ap->len_seq = header->target_len[aln->core.tid];
      int len = 0;
      ap->seq = faidx_fetch_seq(ref, ap->name_seq, 0, ap->len_seq, &len);
    }
    
    current_position = aln->core.pos;
    if(tot_read%10000 == 0 && !silent) print_step(current_position, ap->len_seq);
    
    //if(!strcmp(header->target_name[aln->core.tid], "HS_assemblymutate") && current_position > 96300 && current_position < 96400) ap->debug = 1;
    //else ap->debug = 0;
    //ap->debug = 1;
    int i;
    for(i = previous_position ; i < current_position ; i++)
      select_base(ap, allali, i, s, readname);
    
    // exclude too short alignments
    if(get_lseqali(aln) >= 31) {
      if(!hash_search(readname, bam_get_qname(aln))) {
	add(ap, aln); 
	added_read++;
      }
      else {
	//printf("======> hapremove %s\n", bam_get_qname(aln));
	hash_delete(readname, bam_get_qname(aln));
	bam_destroy1(aln);
      }
      add(allali, aln2);
      aln = bam_init1();
      aln2 = bam_init1();
      tot_read++;
    } else too_short++;
    previous_position = current_position;
  }
  for(i = previous_position ; i < ap->len_seq ; i++) 
    select_base(ap, allali, i, s, readname);
  
  if(!ap->len_seq) {
    print_seq(s, out);
    if(!silent) print_step(ap->len_seq, ap->len_seq);
  }
  nb_changes += ap->nb_changes;
  nb_hapB += ap->nbhapB;
  printf("\n\nNumber of reads               : %i\n", tot_read);
  printf("Number of too short alignment : %i\n", too_short);
  printf("Number of hapB reads found    : %i\n", nb_hapB*2);
  printf("Number of changes             : %i\n", nb_changes);

  fai_destroy(ref);
  bam_hdr_destroy(header);
  alipile_free(ap);
  alipile_free(allali);
  polished_free(s);
  hash_free(readname);
  bam_destroy1(aln);
  bam_destroy1(aln2);
  sam_close(in);
  fclose(changes);
  fclose(out);
  return 0;
}

void print_step(int pos, int len_seq) {	   
  int p = (int)((pos * 100)/ len_seq);
  if(p == 99) p = 100;
  int r = p, i = 50;
  printf("\r|");
  while(p >= 2) {
    printf("*");
    p -= 2;
    i--;
  }
  while(i > 0) { printf(" "); i--; }
  printf("| %i%%", r);
}

void usage() {
  fprintf(stderr, "--------------------------------------------------------------------------------------------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   polish_consensus -b <in.bam> -f <in.fasta> -o <out.fasta> -c <out.changes>\n\n");
  fprintf(stderr, "Parameters: -b FILE : input BAM file\n");
  fprintf(stderr, "            -f FILE : input fasta file\n");
  fprintf(stderr, "            -o FILE : output fasta file\n");
  fprintf(stderr, "            -c FILE : output file that describes the corrections made.\n");
  fprintf(stderr, "            -h      : this help\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "--------------------------------------------------------------------------------------------\n");
  exit(1);
}
