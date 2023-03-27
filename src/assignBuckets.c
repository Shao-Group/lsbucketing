/*
  Input: k

  Assign each k-mer to a set of buckets (int labels) according to the optimal 
  (1,2)-sensitive bucketing function.
  Each k-mer is assigned to k buckets, each bucket contains |\Sigma| k-mers.

  By: Ke@PSU
  Last edited: 03/25/2023
*/

#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ALPHABETSIZE 4
const char alphabet[ALPHABETSIZE] = {'A', 'C', 'G', 'T'};


static inline void addKMerToBucket(const kmer x, const size_t bucket,
				   const int k, size_t* kmers){
    size_t i = x * k;
    while(kmers[i] > 0){ ++i; }
    kmers[i] = bucket;
}

void printKMerBuckets(FILE* fout, const int k,
		      const size_t* kmers, const size_t num_kmers){
    char buf[k+1];
    buf[k] = '\0';
    kmer i;
    size_t j;
    for(i=0; i<num_kmers; ++i){
	decode(i, k, buf);
	fprintf(fout, "%.*s:", k, buf);
	for(j=i*k; j<(i+1)*k; ++j){
	    fprintf(fout, " %zu", kmers[j]);
	}
	fprintf(fout, "\n");
    }
}

/*
  Assign all the buckets for a given kmer x. Results are stored
  in the buckets array from st_idx to st_idx+k-1.

  See the manuscript for explanation of the algorithm.
*/
void assignBuckets(const kmer x, const int k,
		   size_t* buckets, const size_t st_idx){
    int i;
    size_t num_A[k], val[k], mu[k];
    
    kmer mask = 3lu << ((k-1)<<1);
    size_t p = 1lu << ((k-1)<<1); //ALPHABETSIZE^(k-1)
    size_t cur = x & mask;

    size_t sum_mu;
    
    num_A[0] = 0;
    val[0] = x - cur;
    sum_mu = mu[0] = cur ? p + (cur >> 2)*(k-1) : val[0];
    
    for(i=1; i<k; ++i){
	num_A[i] = num_A[i-1] + (cur ? 0 : 1);
	
	mask >>= 2;
	cur = x & mask;
	p >>= 2;
	
	val[i] = val[i-1] - cur;
	mu[i] = cur ? p + (cur >> 2) * (k-i-1) : val[i];
	sum_mu += mu[i];
    }

    mask = 3lu << ((k-1)<<1);
    size_t j=st_idx, tail = st_idx + k - num_A[k-1] - (cur ? 0 : 1);
    
    for(i=0; i<k; ++i){
	cur = x & mask;
	mask >>= 2;
	p = sum_mu - mu[i] + val[i] - num_A[i]*cur + 1 + num_A[i];
	if(cur){
	    buckets[j] = p;
	    ++j;
	}else{
	    buckets[tail] = p;
	    ++tail;
	}
    }
}


int main(int argc, char* argv[]){
    if(argc != 2){
	printf("usage: assignBuckets.out k\n");
	return 1;
    }

    int k = atoi(argv[1]);

    //NOTE: this should be ALPHABETSIZE^{k} for ALPHABETSIZE!=4
    size_t NUM_KMERS = 1<<(k<<1); 

    //buckets for kmer x are at kmers[x*k .. x*k+k-1]
    size_t* kmers = calloc(NUM_KMERS * k, sizeof *kmers);

    size_t m = 1;
    kmer n, mask, s, t;
    int i;
    
    for(n=0; n<NUM_KMERS; ++n){
	mask = 3lu << ((k-1)<<1);
	for(i=k-1; i>=0; --i){
	    if((n & mask) == 0){
		addKMerToBucket(n, m, k, kmers);
		for(s=1; s<ALPHABETSIZE; ++s){
		    t = n | (s << (i<<1));
		    addKMerToBucket(t, m, k, kmers);
		}
		++m;
	    }
	    mask >>= 2;
	}
    }

    char filename[200];
    sprintf(filename, "buckets-%d.txt", k);
    FILE* fout = fopen(filename, "w");

    printKMerBuckets(fout, k, kmers, NUM_KMERS);

    //test assignBuckets function
    size_t individual[k];
    for(n=0, m=0; n<NUM_KMERS; ++n){
	assignBuckets(n, k, individual, 0);
	for(i=0; i<k; ++i){
	    if(individual[i] != kmers[m]){
		char buf[k+1];
		buf[k] = '\0';
		fprintf(stderr, "Wrong buckets for kmer %.*s, should be %zu, assigned %zu\n", k, decode(n, k, buf), kmers[m], individual[i]);
	    }
	    ++ m;
	}
    }
    
    fclose(fout);
    return 0;
}
