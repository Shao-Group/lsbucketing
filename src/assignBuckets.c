/*
  Input: n

  Assign each n-mer to a set of buckets (int labels) according to the optimal 
  (1,2)-sensitive bucketing function.
  Each n-mer is assigned to n buckets, each bucket contains |\Sigma| n-mers.

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
				   const int n, size_t* nmers){
    size_t i = x * n;
    while(nmers[i] > 0){ ++i; }
    nmers[i] = bucket;
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
void assignBuckets(const kmer x, const int n,
		   size_t* buckets, const size_t st_idx){
    int i;
    size_t num_A[n], val[n], mu[n];
    
    kmer mask = 3lu << ((n-1)<<1);
    size_t p = 1lu << ((n-1)<<1); //ALPHABETSIZE^(n-1)
    size_t cur = x & mask;

    size_t sum_mu;
    
    num_A[0] = 0;
    val[0] = x - cur;
    sum_mu = mu[0] = cur ? p + (cur >> 2)*(n-1) : val[0];
    
    for(i=1; i<n; ++i){
	num_A[i] = num_A[i-1] + (cur ? 0 : 1);
	
	mask >>= 2;
	cur = x & mask;
	p >>= 2;
	
	val[i] = val[i-1] - cur;
	mu[i] = cur ? p + (cur >> 2) * (n-i-1) : val[i];
	sum_mu += mu[i];
    }

    mask = 3lu << ((n-1)<<1);
    size_t j=st_idx, tail = st_idx + n - num_A[n-1] - (cur ? 0 : 1);
    
    for(i=0; i<n; ++i){
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
	printf("usage: assignBuckets.out n\n");
	return 1;
    }

    int n = atoi(argv[1]);

    //NOTE: this should be ALPHABETSIZE^{n} for ALPHABETSIZE!=4
    size_t NUM_KMERS = 1<<(n<<1); 

    //buckets for nmer x are at nmers[x*n .. x*n+n-1]
    size_t* nmers = calloc(NUM_KMERS * n, sizeof *nmers);

    size_t m = 1;
    kmer k, mask, s, t;
    int i;
    
    for(k=0; k<NUM_KMERS; ++k){
	mask = 3lu << ((n-1)<<1);
	for(i=n-1; i>=0; --i){
	    if((k & mask) == 0){
		addKMerToBucket(k, m, n, nmers);
		for(s=1; s<ALPHABETSIZE; ++s){
		    t = k | (s << (i<<1));
		    addKMerToBucket(t, m, n, nmers);
		}
		++m;
	    }
	    mask >>= 2;
	}
    }

    char filename[200];
    sprintf(filename, "buckets-%d.txt", n);
    FILE* fout = fopen(filename, "w");

    printKMerBuckets(fout, n, nmers, NUM_KMERS);

    //test assignBuckets function
    size_t individual[k];
    for(k=0, m=0; k<NUM_KMERS; ++k){
	assignBuckets(k, n, individual, 0);
	for(i=0; i<n; ++i){
	    if(individual[i] != nmers[m]){
		char buf[n+1];
		buf[n] = '\0';
		fprintf(stderr, "Wrong buckets for nmer %.*s, should be %zu, assigned %zu\n", n, decode(k, n, buf), nmers[m], individual[i]);
	    }
	    ++ m;
	}
    }
    
    fclose(fout);
    return 0;
}
