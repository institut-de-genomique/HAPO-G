/// @file polished.h
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

#ifndef POLISHED_H
#define POLISHED_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "htslib/hts.h"
#include "htslib/sam.h"

//const unsigned int BUFF_SIZE = 1000;

/*! @typedef
 @abstract Structure for alignment pile information.
 @field  max_cov      Maximal coverage of the pile
 @field  min_cov      Minimal coverage required to correct the consensus
 @field  pile         Array of bam alignment
 @field  seqpile      Array of padded alignment
 @field  current_read Index of the selected Read
 @field  current_seq  Index of the reference sequence
 @field  name_seq     Name of the reference sequence
 @field  seq          Reference sequence
 @field  nb_ali       Number of alignment
 */
typedef struct {
  char* name_seq;
  char* seq;
  int size_alloc;
  int pos_seq;
  int chunk;

  char buffer[1000+1];
  int buffer_level;
} polished_t;


// Create a alipile_t structure
/**
   @return An empty alipile_t structure on success, NULL on failure
 */
polished_t *polished_init(void);

// Free a alipile_t structure
/**
   @param b structure to free
   Does nothing if @p b is NULL.
 */
void polished_free(polished_t *s);

// Add nucleotide base to the polished sequence
void add_base(polished_t *s, char base);

// Add string str to the polished sequence
void add_str(polished_t *s, char* str);

// Print polished sequence in FILE
void print_seq(polished_t *s, FILE* fo);

// Substr function
char* _substr(const char *src,int pos,int len);
// Add a char ath the end of src
void _addchar(char **src, char c);

#endif
