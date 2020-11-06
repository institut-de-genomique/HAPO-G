/// @file hash.c
/// Hash table for string entries

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

#include "hash.h"

const int INIT_SIZE = 268288;
const int EXPAND = 2;
const float MAX_LOAD_FACTOR = 0.7;
const int HASH_FN_OFFSET = 85;


hash_t* hash_init() { 
  return hash_init1(INIT_SIZE); 
}

hash_t* hash_init1(int size) {
  hash_t* h;
  int i = 0;
    
  h = malloc(sizeof(hash_t));
  
  assert(h != 0);
  
  h->msize = size;
  h->size = 0;
  h->htable = malloc(sizeof(elm_t) * h->msize);
  
  assert(h->htable != 0);
  
  for( ; i < h->msize; i++) h->htable[i] = 0;
  
  return h;
}

void hash_free(hash_t* h) {
  int i;
  elm_t *current;
  elm_t *next;
  
  for(i = 0; i < h->msize; i++) {
    for(current = h->htable[i]; current != 0; current = next) {
      next = current->next;
      free(current->key);
      free(current);
    }
  } 
  free(h->htable);
  free(h);
}

unsigned long hash_fn(const char *s) {
  unsigned const char *it = (unsigned const char *) s;
  unsigned long idx = 0;
  
  for( ; *it; it++) idx = HASH_FN_OFFSET * idx + *it;
  
  return idx;
}

void grow(hash_t* h) {
  hash_t* htmp = hash_init1(h->msize * EXPAND);;
  struct hash swap;
  int i = 0;
  elm_t* current;
  
  for( ; i < h->msize; i++)
    for(current = h->htable[i]; current != 0; current = current->next)
      hash_insert(htmp, current->key);
  
  swap = *h;
  *h = *htmp;
  *htmp = swap;
  hash_free(htmp);
}

void hash_insert(hash_t* h, const char *key) {
  elm_t *ins;
  unsigned long hkey;
  
  assert(key);
  
  ins = malloc(sizeof(elm_t));
  
  assert(ins);

  //printf("please delete read= %s\n", key);
  
  ins->key = strdup(key);
  hkey = hash_fn(key) % h->msize;  
  ins->next = h->htable[hkey];
  h->htable[hkey] = ins;
  h->size++;
  
  if(h->size >= (int)(h->msize * MAX_LOAD_FACTOR)) grow(h);
}

int hash_search(hash_t* h, const char *key) {
  unsigned long hkey = hash_fn(key) % h->msize;
  elm_t *se = h->htable[hkey];
  
  for( ; se != 0; se = se->next)
    if(!strcmp(se->key, key)) return 1;

  return 0;
}

void hash_delete(hash_t* h, const char *key) {
  unsigned long hkey = hash_fn(key) % h->msize;
  elm_t **prev = &(h->htable[hkey]);
  elm_t *del;
  
  while(*prev != 0) {
    if(!strcmp((*prev)->key, key)) {
      del = *prev;
      *prev = del->next;
      
      free(del->key);
      free(del);
      h->size--;
      
      return;
    }
    prev = &((*prev)->next);
  }
}
