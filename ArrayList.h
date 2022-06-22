/*
  A dynamic array based list.
  By Ke@PSU
  Last modified: 10/23/2021
*/

#ifndef _ARRAYLIST_H
#define _ARRAYLIST_H 1

#include "util.h" //use alloc_harder family
//#include <stdlib.h>

typedef struct {
    void** arr;
    size_t size;
    size_t used;
} ArrayList;

/*
  Initialize an array to size 16.
*/
void AListInit(ArrayList* list);

/*
  Initialize an array to the given size.
*/
void AListInitSize(ArrayList* list, size_t size);

/*
  Free the entire array. Data parts are freed by the provided freeData function
*/
void AListFree(ArrayList* list, void (*freeData)(void*));

/*
  Trim the list to the minimum size needed.
*/
void AListTrim(ArrayList* list);

/*
  Insert a new k-mer into the array.
  If already full, the size of the array will be doubled before the insertion.
*/
void AListInsert(ArrayList* list, void* data);

/*
  Clear used count so the array can be used as if a new one.
  Pass in a freeData function if data inside needs to be freed.
*/
void AListClear(ArrayList* list, void (*freeData)(void*));

/*
  Swap two lists.
*/
void AListSwap(ArrayList* la, ArrayList* lb);

#endif // ArrayList.h
