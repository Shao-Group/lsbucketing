#include "util.h"

void printIntArray(const int* x, const int len){
  int i;
  for(i=0; i<len; i+=1){
    printf("%2d ", x[i]);
  }
  printf("\n");
}

int editDist(const kmer s1, const kmer s2, const int k, const int max_id){
    return editDist2(s1, k, s2, k, max_id);
}

int editDist2(const kmer s1, const int k1, const kmer s2, const int k2, const int max_d){
    if(k1 > k2) return editDist2(s2, k2, s1, k1, max_d);
    int diag_index = k2 - k1;
    if(max_d >= 0 && diag_index >= max_d) return diag_index;
    
    int row[k2+1];
    int i, j, diag, cur,tmp;
    for(i=0; i<k2+1; i+=1){
	row[i] = i;
    }

    kmer s1_copy, s2_copy;
    for(i=1, s1_copy=s1; i<k1+1; i+=1, s1_copy>>=2){
	diag_index += 1;
	diag = row[0];
	row[0] = i;
	
	for(j=1, s2_copy=s2; j<k2+1; j+=1, s2_copy>>=2){
	    //substitution
	    cur = diag + ((s1_copy & 3) == (s2_copy & 3) ? 0 : 1);
	    //deletion
	    tmp = row[j] + 1;
	    cur = cur > tmp ? tmp : cur;
	    //insertion
	    tmp = row[j-1] +1;
	    cur = cur > tmp ? tmp : cur;
	    
	    diag = row[j];
	    row[j] = cur;
	}

	if(max_d >= 0 && row[diag_index] >= max_d){
	    break;
	}
    }

    return row[diag_index];
}

int editDist3(const char* s1, const int l1, const char* s2, const int l2, const int max_d){
    if(l1 > l2) return editDist3(s2, l2, s1, l1, max_d);
    int diag_index = l2 - l1;
    if(max_d >= 0 && diag_index >= max_d) return diag_index;
    
    int row[l2+1];
    int i, j, diag, cur,tmp;
    for(i=0; i<l2+1; i+=1){
	row[i] = i;
    }

    for(i=1; i<l1+1; i+=1){
	diag_index += 1;
	diag = row[0];
	row[0] = i;
	
	for(j=1; j<l2+1; j+=1){
	    //substitution
	    cur = diag + (s1[i-1] == s2[j-1] ? 0 : 1);
	    //deletion
	    tmp = row[j] + 1;
	    cur = cur > tmp ? tmp : cur;
	    //insertion
	    tmp = row[j-1] +1;
	    cur = cur > tmp ? tmp : cur;
	    
	    diag = row[j];
	    row[j] = cur;
	}

	if(max_d >= 0 && row[diag_index] >= max_d){
	    break;
	}
    }

    return row[diag_index];
}

kmer encode(const char* str, const int k){
    kmer enc = 0;
    int i, x;
    for(i=0; i<k; i+=1){
	switch(str[i]){
	case 'A': x=0; break;
	case 'C': x=1; break;
	case 'G': x=2; break;
	case 'T': x=3; break;
	}
	enc = (enc << 2)|x;
    }
    return enc;
}

char* decode(const kmer enc, const int k, char* str){
    if(str == NULL){
	str = malloc_harder(sizeof *str *k);
    }
    kmer enc_copy = enc;
    char base[] = {'A', 'C', 'G', 'T'};
    int i;
    for(i=k-1; i>=0; i-=1){
	str[i] = base[enc_copy & 3];
	enc_copy >>= 2;
    }
    return str;
}

kmer* readCentersFromFile(const char* filename, const int k,
			  size_t* numOfCenters){
  FILE* fin = fopen(filename, "r");

  fscanf(fin, "%zu\n", numOfCenters);
  kmer* centers = malloc_harder(sizeof *centers *(*numOfCenters));

  int buffersize = k+2;//for \n and \0
  char kmer_str[buffersize];

  size_t i;
  for(i=0; i<*numOfCenters; i+=1){
      fgets(kmer_str, buffersize, fin);
      centers[i] = encode(kmer_str, k);
  }
  
  fclose(fin);
  return centers;
}

int* readCliquesFromFile(const char* filename, const int k, kmer km1Mask,
			 kmer*** centers, size_t* numOfCenters){
    FILE* fin = fopen(filename, "r");

    fscanf(fin, "%zu\n", numOfCenters);
    int* sizes = malloc_harder(sizeof *sizes *(*numOfCenters));
    kmer** c = malloc_harder(sizeof *c *(*numOfCenters));

    size_t buffersize = k+2;//for \n (or ' ') and \0
    char* kmer_str = malloc_harder(sizeof *kmer_str * buffersize);

    size_t i, j, tmp;
    int* cur;
    for(i=0; i<*numOfCenters; i+=1){
	cur = sizes+i;
	fscanf(fin, "%d ", cur);
	c[i] = malloc_harder(sizeof *(c[i]) *(*cur));
	for(j=0; j<(*cur)-1; j+=1){
	    tmp = getdelim(&kmer_str, &buffersize, ' ', fin);
	    if(tmp == k){// (k-1)-mer + ' '
		c[i][j] = encode(kmer_str, k-1) | km1Mask;
	    }else{// k-mer
		c[i][j] = encode(kmer_str, k);
	    }
	}
	tmp = getline(&kmer_str, &buffersize, fin);
	if(tmp == k){// (k-1)-mer + '\n'
	    c[i][j] = encode(kmer_str, k-1) | km1Mask;
	}else{// k-mer
	    c[i][j] = encode(kmer_str, k);
	}
    }

    free(kmer_str);
    fclose(fin);

    *centers = c;
    return sizes;
}


void readKMerHashFromFile(const char* filename, const int k, int* h){
    FILE* fin = fopen(filename, "r");

    char format[20];
    sprintf(format, "%%%ds %%d\n", k);
    char kmer_str[k+1]; //for \0
    int kmer_hash;
    while(fscanf(fin, format, kmer_str, &kmer_hash) != EOF){
	h[encode(kmer_str, k)] = kmer_hash;
    }

    fclose(fin);
}

void* malloc_harder(size_t size){
    void* ptr;
    while((ptr = malloc(size)) == NULL){
	sleep(30);
    }
    return ptr;
}
void* calloc_harder(size_t num, size_t size){
    void* ptr;
    while((ptr = calloc(num, size)) == NULL){
	sleep(30);
    }
    return ptr;
}
void* realloc_harder(void* ptr, size_t new_size){
    void* new_ptr;
    while((new_ptr = realloc(ptr, new_size)) == NULL){
	sleep(30);
    }
    return new_ptr;
}

kmer randomKMer(int k){
    k <<= 1;

    kmer mask = (1lu<<k)-1;
    int rand_digits = 10;
    int rand_max = 1<<rand_digits;
    int parts = k/rand_digits + 1;
    
    kmer s=0lu, cur;
    int i;
    for(i=0; i<parts; i+=1){
	cur = rand()%rand_max;
	s = (s<<rand_digits)|cur;
    }
    return s&mask;
}

//generate a sequence of n distinct random integers from 0 to max-1
//seq is assumed to have at least n spaces
//user is responsible for seeding
static inline void randSeq(int n, int max, int* seq){
    int i, j, m;
    int cur;
    seq[0] = rand()%max;
    for(i=1; i<n; i+=1){
	max -= 1;
	cur = rand()%max;
	for(j=0; j<i; j+=1){
	    if(cur>=seq[j]) cur+=1;
	    else{//insert at this position to remain sorted
		for(m=i-1; m>=j; m-=1){
		    seq[m+1] = seq[m];
		}
		break;
	    }
	}
	seq[j] = cur;
    }
}

static inline void shuffleIntArray(int n, int* seq){
    int i, j, tmp;
    for(i=n-1; i>=1; i-=1){
	j = rand()%(i+1);
	tmp = seq[i];
	seq[i] = seq[j];
	seq[j] = tmp;
    }
}

//00b-A 01b-C 10b-G 11b-T
static inline int randBase(int avoid){
    int n = 4;
    if(avoid >= 0){
	n -= 1;
    }
    int r = rand()%n;
    if(avoid >= 0 && r>=avoid) r+=1;
    return r;
}

kmer randomEdit(kmer s, int k, int d){
    int done = 0;
    kmer t=s;
    
    int numIndel, numSubst;

    int changed[k];
    int ops[k];
    
    int i, j, body;
    kmer head, tail, mask, new_body;
    
    while(!done){
	if(k > d){
	    for(i=0; i<k; i+=1){
		changed[i] = 0;
	    }
	    
	    numIndel = (rand()%((d>>1)+1))<<1;
	    numSubst = d - numIndel;
	    
	    //numIndel distinct positions
	    randSeq(numIndel, k, ops);
	    //after shuffling, treat ops as indices for
	    //(del, ins, del, ins, ...)
	    shuffleIntArray(numIndel, ops);
	    for(i=0; i<numIndel; i+=2){
		//deletion
		j = ops[i]<<1;
		head = (s>>(j+2))<<j;
		tail = ((1lu<<j)-1) & s;
		s = head|tail;
		//insertion
		j = ops[i+1];
		changed[j] = 1;
		j <<= 1;
		head = (s>>j)<<(j+2);
		tail = ((1lu<<j)-1) & s;
		new_body = (long unsigned) randBase(-1);
		s = head | (new_body<<j) | tail;	    
	    }
	    
	    //substitutions
	    for(i=0; i<numSubst; i+=1){
		do{
		    j = rand()%k;
		}while(changed[j]);
		changed[j] = 1;
	    
		j <<= 1;
		mask = 3lu<<j;
		body = (s & mask)>>j;
		new_body = (long unsigned) randBase(body);
		s = (s & ~mask) | (new_body << j);
	    }
	    
	}else{//k==d, substitute all
	    for(i=0; i<k; i+=1){
		j = i << 1;
		mask = 3lu<<j;
		body = (s & mask)>>j;
		new_body = (long unsigned) randBase(body);
		s = (s & ~mask) | (new_body << j);		
	    }
	}

	if(editDist2(s, k, t, k, -1) == d) done = 1;
	else s = t; //restore and try again
    }
    return s;
}


int isSubsequence(kmer x, int b, kmer s, int k){
    int i=0, j=0;
    int cur = x & 3, tmp;
    while(i<b && j<k){
	tmp = s & 3;
	if(tmp == cur){
	    i += 1;
	    x >>= 2;
	    cur = x & 3;
	}
	j += 1;
	s >>= 2;
    }

    return i==b;
}


int isSubstring(kmer x, int l, kmer s, int k){
    int i;
    kmer mask = (1lu<<(l<<1))-1;
    for(i=0; i<=k-l; i+=1){
	if((s^x)&mask){
	    s >>= 2;
	}else{
	    return 1;
	}
    }
    return 0;
}

int isInSampleD1(kmer x, int k){
    kmer mask = 3lu;
    int cur_partition = x & mask;
    int i = 1;
    int cur_symbol;
    while(i<k){
	i+=1;
	x >>= 2;
	cur_symbol = x & mask;
	if(cur_partition < cur_symbol) cur_partition += 4;
	cur_partition -= cur_symbol;
    }
    return (cur_partition == 0);
}
