/// @file alipile.c
/// Pile of alignment and padded alignments

/*
################################################################################
# * Copyright Jean-Marc Aury / Genoscope / DRF / CEA
# *           <jmaury@genoscope.cns.fr>
# *
# * This software is governed by the CeCILL license under French law and
# * abiding by the rules of distribution of free software.  You can  use,
# * modify and/ or redistribute the software under the terms of the CeCILL
# * license as circulated by CEA, CNRS and INRIA at the following URL
# * "http://www.cecill.info".
# *
# * As a counterpart to the access to the source code and  rights to copy,
# * modify and redistribute granted by the license, users are provided only
# * with a limited warranty  and the software's author,  the holder of the
# * economic rights,  and the successive licensors  have only  limited
# * liability.
# *
# * In this respect, the user's attention is drawn to the risks associated
# * with loading,  using,  modifying and/or developing or reproducing the
# * software by the user in light of its specific status of free software,
# * that may mean  that it is complicated to manipulate,  and  that  also
# * therefore means  that it is reserved for developers  and  experienced
# * professionals having in-depth computer knowledge. Users are therefore
# * encouraged to load and test the software's suitability as regards their
# * requirements in conditions enabling the security of their systems and/or
# * data to be ensured and,  more generally, to use and operate it in the
# * same conditions as regards security.
# *
# * The fact that you are presently reading this means that you have had
# * knowledge of the CeCILL license and that you accept its terms.
################################################################################
*/

#include "alipile.h"

int compare( const void* a, const void* b) {
  int int_a = * ( (int*) a );
  int int_b = * ( (int*) b );
  
  // an easy expression for comparing
  return (int_a > int_b) - (int_a < int_b);
}

alipile_t *alipile_init(FILE* out) {
  alipile_t* ap = (alipile_t*) calloc(1, sizeof(alipile_t));
  ap->max_cov = 500;
  ap->current_read = -1;
  ap->read_choice = 0;
  ap->current_seq = -1;
  ap->min_cov = 3;
  ap->nb_ali = 0;
  ap->debug = 0;
  ap->nb_changes = 0;
  ap->changes = out;
  ap->nbhapB = 0;
  return ap;
}

void alipile_free(alipile_t *ap) {
  int i;
  for (i = 0; i < ap->nb_ali; i++) free(ap->seqpile[i]);
  //free(ap->seqpile);
  for (i = 0; i < ap->nb_ali; i++) bam_destroy1(ap->pile[i]);
  //free(ap->pile);
  //free(ap->name_seq);
  free(ap->seq);
  free(ap);
}

void clean(alipile_t *ap, int pos) {
  int idread[ap->nb_ali];
  int nbdel = 0, c = 0;
  // retrieve all possibilities
  for( c = 0 ; c < ap->nb_ali ; c++) {
    bam1_t *b = ap->pile[c];
    if(bam_endpos(b) <= pos+10) {
      idread[nbdel] = c;
      nbdel++;
    }
  }
  for( c = nbdel-1 ; c >= 0 ; c-- )
    delete(ap, idread[c]);
}

void add(alipile_t *ap, bam1_t *aln) {
  int c = 0, smallest = -1, value = 0;
  clean(ap, aln->core.pos);

  //need to delete an alignment
  if(ap->nb_ali == ap->max_cov) {
    for( c = 0, value = 0 ; c < ap->nb_ali ; c++) {
      bam1_t *b = ap->pile[c];
      int size = bam_endpos(b) - aln->core.pos;
      if(smallest == -1 || size < value) { smallest = c; value = size; }
    }
    delete(ap, smallest);
  }
  
  //add the new alignment
  ap->pile[ap->nb_ali] = aln;
  ap->seqpile[ap->nb_ali] = get_seqali(aln);
  //printf("->add: %s\tseq=%s\tnbreads=%i\tcurrentpos=%i\n", bam_get_qname(aln), ap->seqpile[ap->nb_ali], ap->nb_ali, aln->core.pos);
  ap->nb_ali++;
  
  // if current_read has been deleted, need to select a new read
  if(ap->current_read == -1 || smallest == ap->current_read)
    change_currentread(ap, aln->core.pos);
}

void change_currentread(alipile_t *ap, int pos) {
  int c = 0, smallest = -1, value = 0;
  for( c = 0 ; c < ap->nb_ali ; c++) {
    bam1_t *b = ap->pile[c];
    int size = bam_endpos(b) - pos;
    if(smallest == -1 || size < value) { smallest = c; value = size; }
  }
  ap->current_read = smallest;
}

void change_currentread2(alipile_t *ap, int pos) {
  int c = 0, smallest = -1;
  float value = 0.0;
  for( c = 0 ; c < ap->nb_ali ; c++) {
    char *nt = get_base(ap, c, pos);
    float ratio = get_ratio_base(ap, nt, pos);
    free(nt);
    if(smallest == -1 || ratio > value) { smallest = c; value = ratio; }
  }
  ap->current_read = smallest;
}

void delete(alipile_t *ap, int c) {
  assert( c < ap->nb_ali );
  bam1_t *remove = ap->pile[c];
  //printf("->remove: %s\tseq=%s\tindex=%i\n", bam_get_qname(remove), get_bamseq(remove), elem);
  bam_destroy1(remove);
  free(ap->seqpile[c]);
  while(c+1 < ap->max_cov) {
    ap->pile[c] = ap->pile[c+1];
    ap->seqpile[c] = ap->seqpile[c+1];
    c++;
  } 
  ap->nb_ali--;
  if(ap->current_read >= ap->nb_ali)
    ap->current_read = 0;
}

int get_lseqali(bam1_t *aln) {
  int k, c;
  uint32_t *cigar = bam_get_cigar(aln);
  //calculate size of ali string
  for (k = 0, c = 0; k < aln->core.n_cigar; ++k) {
    int op, ol;
    op = bam_cigar_op(cigar[k]);
    ol = bam_cigar_oplen(cigar[k]);
    if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || 
       op == BAM_CINS || op == BAM_CDEL) c+=ol;
  }
  return c;
}

char* get_seqali(bam1_t *aln) {
  int i, j, k, c, qlen = 0;
  uint32_t *cigar = bam_get_cigar(aln);
  uint8_t *seq = bam_get_seq(aln);
  int length = get_lseqali(aln);
  char *r = malloc((length+1) * sizeof(char));

  for (k = 0, j = 0, c = 0; k < aln->core.n_cigar; ++k) {
    int op, ol;
    op = bam_cigar_op(cigar[k]);
    ol = bam_cigar_oplen(cigar[k]);
    assert(c <= length);
    if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      for (i = 0; i < ol; ++i, ++j, ++c) r[c] = seq_nt16_str[bam_seqi(seq, j)];
    }
    if(op == BAM_CSOFT_CLIP) {
      for (i = 0; i < ol; ++i, ++j);
    }
    if(op == BAM_CINS) {
      for (i = 0; i < ol; ++i, ++j, ++c) r[c] = tolower(seq_nt16_str[bam_seqi(seq, j)]);
    }
    if(op == BAM_CDEL) {
      for (i = 0; i < ol; ++i, ++c) r[c] = '-';
    }
  }
  r[length]='\0';
  return r;
}

char* get_bamseq(bam1_t *aln) {
  int i = 0;
  const bam1_core_t *c = &aln->core;
  char *cp = malloc(c->l_qseq + 1);
  uint8_t *s = bam_get_seq(aln);
  for (i = 0; i < c->l_qseq; i++)
    cp[i] = seq_nt16_str[bam_seqi(s, i)];
  cp[i] = '\0';
  
  return cp;
}

char* get_base(alipile_t *ap, int read, int pos) {
  char* region = malloc(1000 * sizeof(char));
  region[0]= '\0';
  char* aliseq = ap->seqpile[read];
  bam1_t* aln = ap->pile[read];
  int i = 0, c = 0, len = strlen(aliseq);
  for (i = aln->core.pos , c = 0; c < len ; c++) {
    char nt = aliseq[c];
    if(i == pos) {
      int l = strlen(region);
      if(l == 0 || nt != '-') {
	region[l] = nt;
	region[l+1] = '\0';
      }
    }
    // lowercase indicate insertion in read
    if(!(nt >= 'a' && nt <= 'z')) i++;
  }
  return region;
}

float get_ratio_base(alipile_t *ap, char* str, int pos) {
  int i = 0, nb_id = 0, tot = 0;
  for( i = 0 ; i < ap->nb_ali ; i++ ) {
    char *seq = get_base(ap, i, pos);
    if(strcmp(str, seq) == 0) nb_id++;
    tot++;
    free(seq);
  }
  return ((float)nb_id/(float)tot);
}

void removehapB_read(alipile_t *ap, char* str, int pos, hash_t* readname) {
  int i = 0, j = 0, id = -1, nbdiff = 0, vu = 0;
  if(ap->debug) printf("**** Remove hapB **** hapA=%s read=%s, seqpile=%s aliseq=%s\n", str, bam_get_qname(ap->pile[ap->current_read]), ap->seqpile[ap->current_read], get_seqali(ap->pile[ap->current_read]));
  char *seen[ap->nb_ali];
  int nbocc[ap->nb_ali];
  int idread[ap->nb_ali][ap->nb_ali];
  // retrieve all possibilities
  for( i = 0 ; i < ap->nb_ali ; i++ ) {
    char *seq = get_base(ap, i, pos);
    vu = 0;
    for( j = 0 ; j < nbdiff ; j++) 
      if(strcmp(seen[j], seq) == 0) { vu = 1; idread[j][nbocc[j]]=i; nbocc[j]++; }
    if(!vu) { 
      seen[nbdiff] = seq;
      idread[nbdiff][0]=i; 
      nbocc[nbdiff] = 1; 
      if(strcmp(str, seq) == 0) id = nbdiff;
      nbdiff++;      
    } else { free(seq); }
  }
  // erase reads from hapB
  int nbdelread = 0, delread[ap->nb_ali];
  for( j = 0 ; j < nbdiff ; j++) {
    float ratio = ((float)nbocc[j] / (float)ap->nb_ali);
    if(ap->debug) 
      printf("======> found base=%s\t nbocc=%i\t ratio=%.4f\n", 
	     seen[j], nbocc[j], ratio);
    if( j == id ) continue;
    for( i = nbocc[j]-1 ; i >= 0 ; i-- ) {
      bam1_t *aln = ap->pile[idread[j][i]];
      
      if(ap->debug) 
	printf("======> erase readid= %s\t base=%s\t ratio=%.4f\n", 
	       bam_get_qname(aln), seen[j], ratio);
      
      if(!hash_search(readname, bam_get_qname(aln))) {
	hash_insert(readname, bam_get_qname(aln));
	ap->nbhapB++;
      } else
	hash_delete(readname, bam_get_qname(aln));

      delread[nbdelread++] = idread[j][i];
    }
  }
  qsort(delread, nbdelread, sizeof(int), compare);
  for( i = nbdelread - 1 ; i >= 0 ; i--) delete(ap, delread[i]);

  for( j = 0 ; j < nbdiff ; j++) free(seen[j]);
}

void select_base(alipile_t *ap, alipile_t *allali, int current_pos, 
		 polished_t *s, hash_t* readname) {
  char ref[2] = {toupper(ap->seq[current_pos])};
  char *read = NULL;
  int ret = 0;
  
  // clean the pile
  clean(ap, current_pos);
  clean(allali, current_pos);

  // coverage is too low, keep reference sequence
  if(ap->nb_ali < ap->min_cov) {
    if(ap->debug) 
      printf("==> seq=%s pos=%i base_ref=%s nbali=%i [COV TOO LOW]\n", 
	     ap->name_seq, current_pos, ref, ap->nb_ali);
    add_str(s, ref);
    ret = 1;
  }
  
  if(!ret) {
    read = get_base(ap, ap->current_read, current_pos);
  
    float ratio = get_ratio_base(ap, read, current_pos);
    float ratioall = get_ratio_base(allali, read, current_pos);
    if(ap->debug) 
      printf("==> seq=%s pos=%i base_ref=%s base_read=%s ratio=%f ratioall=%f nbali=%i read=%s\n", 
	     ap->name_seq, current_pos, ref, read, ratio, ratioall, ap->nb_ali,
	     bam_get_qname(ap->pile[ap->current_read]));
    
    // keep read base (homo diff)
    if(!ret && ratioall >= 0.8 && (ratio * ap->nb_ali) > ap->min_cov) {
      if(strcmp(read, "-") != 0) add_str(s, read);
      if(strcmp(ref, read) != 0) {
	fprintf(ap->changes, "%s\t%i\tref=%s\tread=%s\treadname=%s\thomo\tratio1=%.4f\tratio2=%.4f\n", 
		ap->name_seq, current_pos, ref, read, 
		 bam_get_qname(ap->pile[ap->current_read]), ratio, ratioall);

	ap->nb_changes++;
	if(ap->debug) 
	  printf("%s\t%i\tref=%s\tread=%s\treadname=%s\thomo\tratio1=%.4f\tratio2=%.4f\n", 
		ap->name_seq, current_pos, ref, read, 
		 bam_get_qname(ap->pile[ap->current_read]), ratio, ratioall);
      }
      ret = 1;
    }
    
    // keep read base (hetero diff)
    if(!ret && ap->nb_ali > 2*ap->min_cov && ratioall > 0.2 && ratioall < 0.8) {
      // if both ratios are equal (first hetero variation), we select the haplotype with the highest coverage
      float r = (ratio<ratioall) ? ratioall-ratio : ratio-ratioall;
      if(r<0.05) {
	change_currentread2(ap, current_pos);
	free(read);
	read = get_base(ap, ap->current_read, current_pos);
	if(ap->debug) 
	  printf("CHANGE %s\t%i\tref=%s\tread=%s\treadname=%s\thetero\tratio1=%.4f\tratio2=%.4f\n",
		 ap->name_seq, current_pos, ref, read, 
		 bam_get_qname(ap->pile[ap->current_read]), ratio, ratioall);	
      }

      removehapB_read(ap, read, current_pos, readname);
      if(strcmp(read, "-") != 0) add_str(s, read);
      if(strcmp(ref, read) != 0) {
	fprintf(ap->changes, "%s\t%i\tref=%s\tread=%s\treadname=%s\thetero\tratio1=%.4f\tratio2=%.4f\n", 
		ap->name_seq, current_pos, ref, read, 
		bam_get_qname(ap->pile[ap->current_read]), ratio, ratioall);
	ap->nb_changes++;
	if(ap->debug) 
	  printf("%s\t%i\tref=%s\tread=%s\treadname=%s\thetero\tratio1=%.4f\tratio2=%.4f\n",
		 ap->name_seq, current_pos, ref, read, 
		 bam_get_qname(ap->pile[ap->current_read]), ratio, ratioall);
      }
      ret = 1;
    }
    
    // sequencing error, delete the read and select a new read 
    if(!ret && ratioall <= 0.2) {
      delete(ap, ap->current_read);
      ap->current_read = -1;
      change_currentread2(ap, current_pos);
      return select_base(ap, allali, current_pos, s, readname);
    }

    // special case : keep reference base
    if(!ret) add_str(s, ref);
  }

  free(read);
}
