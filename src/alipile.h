/// @file alipile.h
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

#ifndef ALIPILE_H
#define ALIPILE_H

#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>
#include <string.h>

#include "polished.h"
#include "hash.h"

#include "htslib/hts.h"
#include "htslib/sam.h"

//const unsigned int MAX_COV = 50;
//const unsigned int MIN_COV = 3;

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
 @field  len_seq      Reference sequence length
 @field  nb_ali       Number of alignment
 @field  nbhapB       Number of removed PE reads from hapB
 @field  changes      FILE output changes
 */
typedef struct {
  int max_cov;
  int min_cov;

  bam1_t* pile[500];
  char* seqpile[500];

  int current_read;
  int read_choice;

  int current_seq;
  char* name_seq;
  char* seq;
  int len_seq;

  int nb_ali;
  int debug;

  int nb_changes;
  int nbhapB;

  FILE* changes;
} alipile_t;


// Create a alipile_t structure
/**
   @return An empty alipile_t structure on success, NULL on failure
 */
alipile_t *alipile_init(FILE* out);

// Free a alipile_t structure
/**
   @param b structure to free
   Does nothing if @p b is NULL.
 */
void alipile_free(alipile_t *ap);

// Clean pile by erasing alignment that do not overlap pos
void clean(alipile_t *ap, int pos);

// Add an alignment to the alipile_t structure
void add(alipile_t *ap, bam1_t *aln);

// Change current read in the alipile_t structure
void change_currentread(alipile_t *ap, int pos);

// Change current read in the alipile_t structure
void change_currentread2(alipile_t *ap, int pos);

// Delete an alignment from the alipile_t structure
void delete(alipile_t *ap, int elem);

// Get padded-alignment length
int get_lseqali(bam1_t *aln);

// Get padded-alignment
char* get_seqali(bam1_t *aln);

// Get sequence from bam alignment
char* get_bamseq(bam1_t *aln);

// Return the base of the read at position pos
char* get_base(alipile_t *ap, int read, int pos);

// return the proportion of ali with base nt at given position
float get_ratio_base(alipile_t *ap, char *nt, int pos);

// Remove alignments from haplotype B
void removehapB_read(alipile_t *ap, char* str, int pos, hash_t* readname);

// Get the corrected nucleotide at position current_pos
void select_base(alipile_t *ap, alipile_t *allali, int current_pos, 
		 polished_t* s, hash_t* readname);


#endif
