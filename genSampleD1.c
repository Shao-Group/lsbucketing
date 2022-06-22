/*
  Input: k

  Generate a subset S of all length-k sequences satisfying the condition that
  for every length-k sequence s and a position 1<=i<=k, one of the four 
  mutations of s (including s itself) at position i is in S. 

  |S|=4^{k-1}. In fact, the space of all length-k sequences can be partitioned 
  into four equal-size subsets S1 U S2 U S3 U S4 such that each Si satisfies 
  the aforementioned condition. The iterative algorithm is based on 
  this partition.

  The base case is when k=1, the partition is {A}U{C}U{G}U{T}.
  Suppose we already have the partition for k-1, U S'i. Then for k,
  S1=AS1 U CS2 U GS3 U TS4 where AS1 is the set obtained by 
  prepending A to each (k-1)-mer in S1, similar for the rest.
  S2=AS2 U CS3 U GS4 U TS1,
  S3=AS3 U CS4 U GS1 U TS2,
  S4=AS4 U CS1 U GS2 U TS3.

  By: Ke@PSU
  Last edited: 02/08/2022
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ASIZE 4
const char alphabet[ASIZE] = {'A', 'C', 'G', 'T'};

void printSamplesToFile(int k, char* samples, size_t size){
    char filename[50];
    sprintf(filename, "%02d01.sample", k);
    FILE* fout = fopen(filename, "w");
    fprintf(fout, "%zu\n", size);

    size *= k;
    size_t i;
    for(i=0; i<size; i+=k){
	fprintf(fout, "%.*s\n", k, samples+i);
    }

    fclose(fout);
}

void printSamples(int k, char* samples, size_t size){
    size_t i;
    int j;
    for(i=0; i<size; i+=1){
	for(j=0; j<k; j+=1){
	    printf("%c", samples[i*k+j]=='\0'?'_':samples[i*k+j]);
	}
	printf("\n");
    }
}

int main(int argc, char* argv[]){
    if(argc != 2){
	printf("usage: genSampleD1.out k\n");
	return 1;
    }

    int k = atoi(argv[1]);
    if(k <= 1){
	printf("input k should be at least 2\n");
	return 1;
    }

    //NOTE: this should be ASIZE^{k-1} for ASIZE!=4
    size_t NUM_SAMPLES = 1<<((k<<1)-2); 

    char* samples = malloc(sizeof *samples * NUM_SAMPLES*k);

    size_t i;
    for(i=0; i<ASIZE; i+=1){
	samples[i*k+k-1] = alphabet[i];
    }

    int j, cur_idx = k-2, group;
    size_t cur_size = ASIZE * k;
    size_t group_size = k;
    size_t dest, src, ct;
    while(cur_idx > 0){	
	//make another ASIZE-1 copies of the current samples, with rotation
	//i.e., S0,S1,S2,S3 -> S0,S1,S2,S3|S1,S2,S3,S0|S2,S3,S0,S1|S3,S0,S1,S2
	dest = cur_size;
	for(j=1; j<ASIZE; j+=1){
	    for(group = 0; group<ASIZE; group+=1, dest+=group_size){
		src = ((j+group)%ASIZE)*group_size;
		memcpy(samples+dest, samples+src, group_size);
	    }
	}

	//now we have ASIZE groups of rotated samples, each group
	//contains ASIZE subgroups.
	//For each group, prepend the j-th alphabet to the j-th subgroup
	i = 0;
	for(group=0; group<ASIZE; group+=1){
	    for(j=0; j<ASIZE; j+=1){
		for(ct=0; ct<group_size; ct+=k, i+=k){
		    samples[i+cur_idx] = alphabet[j];
		}
	    }
	}

	cur_idx -= 1;
	cur_size <<= 2; //cur_size *= ASIZE;
	group_size <<= 2; //group_size *= ASIZE;
    }

    //now we have partitioned all (k-1)-mers into ASIZE groups
    //prepend the j-th alphabet to the j-th group of samples
    i = 0;
    for(j=0; j<ASIZE; j+=1){
	for(ct=0; ct<group_size; ct+=k, i+=k){
	    samples[i] = alphabet[j];
	}
    }
    
    printSamplesToFile(k, samples, NUM_SAMPLES);

    free(samples);
    return 0;
}
