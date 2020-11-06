/// @file hash.h
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

#ifndef HASH_H
#define HASH_H

#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>
#include <string.h>

/*! @typedef
 @abstract Structure for a hash entry
 @field  next         Next entry to manage colision
 @field  key          The key of the entry
 */
typedef struct elm {
    struct elm *next;
    char *key;
} elm_t;

/*! @typedef
 @abstract Structure for a hash table
 @field  m_size       Size of the table
 @field  size         The number of entries
 @field  table        The hash table
 */
typedef struct hash {
    int msize;
    int size;
    elm_t **htable;
} hash_t;


/* create and initialize an empty hash table */
hash_t* hash_init(void);
hash_t* hash_init1(int);

/* destroy a hash table */
void hash_free(hash_t*);

/* insert a new key into a hash table */
void hash_insert(hash_t*, const char*);

/* return 1 if the key is present, 0 otherwise */
int hash_search(hash_t*, const char*);

/* delete the record with the given key */
/* if there is no such record, has no effect */
void hash_delete(hash_t*, const char*);

#endif
