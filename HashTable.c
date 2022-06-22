#include "HashTable.h"

void HTableInit(HashTable* table){
    HTableInitSize(table, 64);
}

void HTableInitSize(HashTable* table, size_t size){
    table->arr = calloc_harder(size, sizeof *table->arr);
    table->size = size;
    table->used = 0;
}

void HTableFree(HashTable* table){
    free(table->arr);
    table->arr = NULL;
    table->used = 0;
    table->size = 0;
}

//return the next available index (including cur_index),
//existance of such a position is guaranteed by the load factor
//if found value along the way, return that position
static inline size_t HTableProbe(HashTable* table, size_t cur_index,
				 long unsigned value){
    while(table->arr[cur_index] != 0 && table->arr[cur_index] != value){
	cur_index += 1;
	if(cur_index >= table->size) cur_index = 0;
    }
    return cur_index;
}

//the hash function
static inline void HTableHash(const size_t size, long unsigned enc, size_t* key, long unsigned* value){
    *key = enc % size;
    *value = enc + 1;
}
static inline void HTableUnhash(long unsigned value, long unsigned* enc){
    *enc = value - 1 ;
}

//add enc to table without checking load factor and incrementing used
//return 1 if added, return 0 if (key, value) is already in table
static inline int HTableAdd(HashTable* table, long unsigned enc){
    size_t key;
    unsigned long value;
    HTableHash(table->size, enc, &key, &value);
    size_t pos = HTableProbe(table, key, value);
    if(table->arr[pos] == 0){
	table->arr[pos] = value;
	return 1;
    }else{
	return 0;
    }
}

void HTableResize(HashTable* table, size_t size){
    long unsigned* old_arr = table -> arr;
    size_t old_size = table -> size;
    table->arr = calloc_harder(size, sizeof *table->arr);
    table->size = size;

    size_t i;
    unsigned long enc;
    for(i=0; i<old_size; i+=1){
	if(old_arr[i]){
	    HTableUnhash(old_arr[i], &enc);
	    HTableAdd(table, enc);
	}
    }

    free(old_arr);
}


void HTableInsert(HashTable* table, long unsigned enc){
    if(table->size < (table->used << 1)){
        HTableResize(table, table->size<<1);
    }
    if(HTableAdd(table, enc)) table->used += 1;
}

int HTableSearch(HashTable* table, long unsigned enc){
    size_t key;
    unsigned long value;
    HTableHash(table->size, enc, &key, &value);
    if(table->arr[HTableProbe(table, key, value)] == 0) return 0;
    else return 1;
}

long unsigned* HTableToArray(HashTable* table, long unsigned* list){
    if(list == NULL){
	list = malloc_harder(sizeof *list *table->used);
    }

    int i, j=0;
    long unsigned cur, enc;
    for(i=0; i<table->size; i+=1){
	cur = table->arr[i];
	if(cur != 0){
	    HTableUnhash(cur, &enc);
	    list[j++] = enc;
	}
    }
    
    return list;
}
