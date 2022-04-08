/// @file polished.c
/// Polished (error-corrected) sequence

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

#include "polished.h"

polished_t *polished_init(void) {
  polished_t* s = (polished_t*) calloc(1, sizeof(polished_t));
  s->chunk = 1000;
  s->size_alloc = s->chunk;
  s->seq = (char *) malloc((s->chunk + 1) * sizeof(char));
  s->seq[s->chunk]='\0';
  s->pos_seq = 0;
  s->buffer_level = 0;
  return s;
}

void polished_free(polished_t *s) {
  if(s == 0) return;
  //free(s->name_seq);
  // on redonne la taille allouee
  //s->seq[strlen(s->seq)] = 'e';
  //s->seq[s->size_alloc] = '\0';
  free(s->seq);
  free(s);
}


void add_base(polished_t *s, char base) {
  char tmp[2] = { base };
  add_str(s, tmp);
  //if(s->seq == NULL) s->seq = (char *) malloc(2 * sizeof(char));
  //else  { 
  //  int oldl = strlen(s->seq);
  //  s->seq = (char *) realloc(s->seq, (oldl + 1) * sizeof(char));
  //}
  //s->seq[strlen(s->seq)-1] = base;
  //s->seq[strlen(s->seq)] = '\0';
}

void add_str(polished_t *s, char *str) {
  assert(str != NULL && s->seq != NULL);
  int len = strlen(str);
  if(s->pos_seq + len >= s->size_alloc) {
    int extend = s->chunk;
    while(len >= extend) { extend += s->chunk; }
    void* tmpseq = (char *) realloc(s->seq, (s->size_alloc + extend) * sizeof(char));  
    assert(tmpseq);
    s->seq = tmpseq;
    s->size_alloc += extend;
  }
  int i = 0;
  for(i = 0 ; i < len ; i++)
    s->seq[s->pos_seq + i ] = str[i];
  s->pos_seq += len;
  s->seq[s->pos_seq]='\0';
}

/*void addBase(polished_t *s, char base) {
  assert(buffer_level >= BUFF_SIZE);
  
  s->buffer[s->buffer_level] = base;
  s->buffer_level++;

  if(buffer_level + 1 == BUFF_SIZE) {
    if(s->seq == NULL) s->seq = (char *) malloc(BUFF_SIZE * sizeof(char));
    else  { 
      int oldl = strlen(s->seq);
      s->seq = (char *) realloc(s->seq, (oldl + BUFF_SIZE) * sizeof(char));
    }
    strcat(s->seq, buffer);
    buffer_level = 0;
  } 
  }*/

void print_seq(polished_t *s, FILE* fo) {
  if(s->seq != NULL) {
    char *cs = s->seq;
    fprintf(fo, ">%s_polished\n", s->name_seq);
    int linesize = 60;
    int len = strlen(s->seq);
    int write = len;
    while(len>0) {
      write = len;
      if(len > linesize) write = linesize;
      char* line = _substr(s->seq, 0, write);
      fprintf(fo, "%s\n", line);
      free(line);
      s->seq += write;
      len -= write;
    }
    s->seq = cs;
  }
}

char* _substr(const char *src, int pos, int len) {
  char *dest = NULL;
  if (len>0) {
    dest = calloc(len+1, sizeof(char));
    if(NULL != dest) strncat(dest, src+pos, len);
  }
  return dest;        
}


void _addchar(char **src, char c) {
  if(*src == NULL) {
    *src = malloc(2 * sizeof(char));
    *src[0] = c;
    *src[1] = '\0';
  } else { 
    int oldl = strlen(*src);
    *src = realloc(*src, (oldl + 1) * sizeof(char));
    //printf("realloc seq= \"%s\"\t char= %c\t len=%i\n", *src, c, strlen(*src));
    *src[strlen(*src)-1] = c;
    *src[strlen(*src)] = '\0';
  }
  //printf("seq= \"%s\"\t char= %c\t len=%i\n", *src, c, strlen(*src));
}
  
