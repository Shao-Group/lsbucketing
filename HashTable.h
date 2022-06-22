/*
  A dynamic array based hashtable for k-mers (long unsigned int) with 
  linear probing. Deletion is not implemented.
  By Ke@PSU
  Last modified: 10/16/2021
*/

#ifndef _HASHTABLE_H
#define _HASHTABLE_H 1

#include "util.h"
//#include <stdlib.h>

typedef struct {
    long unsigned* arr;
    size_t size;
    size_t used;
} HashTable;

/*
  Initialize a hashtable with size 64.
*/
void HTableInit(HashTable* table);

/*
  Initialize a hashtable with the given size.
*/
void HTableInitSize(HashTable* table, size_t size);

/*
  Free the entire table.
*/
void HTableFree(HashTable* table);

/*
  Insert a new k-mer into the table. The k-mer itself (long unsigned int % size)
  is used as the key, k-mer+1 is stored as the value (so 0 can be used to
  indicate an empty space). If position is taken, probe the next index in a circular
  manner. 

  If size/used < 2, size of the table will be doubled (and all the entries
  rehashed) before the insertion.
*/
void HTableInsert(HashTable* table, long unsigned enc);

/*
  Search in the table for the given k-mer. (See HTableInsert for hashing and
  probing details.) Return 1 if found, 0 otherwise.
 */
int HTableSearch(HashTable* table, long unsigned enc);

/*
  Dump everything in the table to an ArrayList, if list is NULL, a new
  list will be allocated.
*/
long unsigned* HTableToArray(HashTable* table, long unsigned* list);

#endif // HashTable.h
