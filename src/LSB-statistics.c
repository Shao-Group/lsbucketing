/*
  Input: k r w(hole)|s(ample)

  For d=1, 2, ..., 6, generata N pairs of length-k sequences with 
  edit distance d. A pair (s, t) is said to have a collision if they
  share an r-neighbor [if option s, then the neighbor is required to
  be in the (1,1)-guaranteed sample tested by isInSampleD1()].
  Show the frequencies of collisions.
  
  When r=1, a gap (neither 100% nor 0% collision probability) is expected
  when d=2 for both options whole and sample. So we further distinguish
  two cases:
  i) s differs from t by 2 substituitions (exactly 2 mismatches)
  ii) s differs from t by 1 indels (more than 2 mismatches).

  When r=2, a gap is expected when d=4 for the (1,1)-guaranteed sample set.
  So we further distinguish three cases:
  i) s differs from t by 4 substituitions (exactly 4 mismatches)
  ii) s differs from t by 2 substituitions and 1 indel (more than 4 mismatches,
  exists a mismatch that if changed would result in edit=3)
  iii) s differs from t by 2 indels (more than 4 mismatches, no mismatch can
  be changed to make edit=3).

  By: Ke@PSU
  Last edited: 05/22/2022
*/

#include "util.h"
#include "AVLTree.h"
#include "ArrayList.h"
#include "HashTable.h"
#include <time.h>
#include <string.h>

#define N 100000

//mask for (k-1)-mers
#define luMSB 0x8000000000000000lu

int cmpKMer(const void* a, const void* b){
    kmer s = (kmer) a;
    kmer t = (kmer) b;
    if(s==t) return 0;
    else if(s<t) return -1;
    else return 1;
}

//see beginning comment for details
int getEditType(kmer s, kmer t, int k, int d){
    int i, ct=0;
    kmer mask=3lu;
    if(d==2){
	for(i=0; i<k; i+=1){
	    if((s&mask) != (t&mask)){
		ct += 1;
	    }
	    mask <<= 2;
	}
	if(ct == 2) return 0;
	else return 1;
    }else if(d!=4){
	return -1;
    }
    
    int hasSub=0;
    kmer new_s, tmp;
    for(i=0; i<k; i+=1){
	if((s&mask) != (t&mask)){
	    ct += 1;
	}
	mask <<= 2;
    }
    if(ct == 4) return 0; //4 subs

    mask=3lu;
    //need at least 1 indel
    //check if it's possible to do a sub at certain position i
    //so that the resulting edit distance is 3
    //this sub may happen between i in s and {i-1, i, i+1} in t
    for(i=0; i<k; i+=1){
	if(i>0){//sub with i-1 in t
	    tmp = t & (mask>>2);
	    if((s&mask) != tmp){
		new_s = (s&~mask)|(tmp<<2);
		if(editDist2(new_s, k, t, k, 4) == 3){
		    hasSub = 1;
		    break;
		}
	    }
	}
	//sub with i in t
	tmp = t&mask;
	if((s&mask) != tmp){
	    new_s = (s&~mask)|tmp;
	    if(editDist2(new_s, k, t, k, 4) == 3){
		hasSub = 1;
		break;
	    }
	}
	if(i<k-1){//sub with i+1 in t
	    tmp = t & (mask<<2);
	    if((s&mask) != tmp){
		new_s = (s&~mask)|(tmp>>2);
		if(editDist2(new_s, k, t, k, 4) == 3){
		    hasSub = 1;
		    break;
		}
	    }
	}
	mask <<= 2;
    }
    if(hasSub) return 1; //2 subs 1 indel
    else return 2;//2 indels
}

/*
  Do bfs for r layers from the given kmer cur, for each
  that is isInSampleD1, add to the resulting hs
*/
AVLNode* bfsNeighborsInSampleRadius(kmer cur, int k, int r, int check_sample){
    AVLNode *hs = NULL;
    if(!check_sample || isInSampleD1(cur, k)){
	hs = AVLAdd(hs, (void*)cur, cmpKMer);
    }

    HashTable visited;
    HTableInit(&visited);
    HTableInsert(&visited, cur);

    ArrayList cur_layer, next_layer;
    AListInit(&cur_layer);
    AListInit(&next_layer);
    AListInsert(&cur_layer, (void*)cur);

    
    size_t i, j;
    int depth;
    kmer t, head, body, tail, x, m;

    for(depth=1; depth<=r; depth+=1){
	for(i=0; i<cur_layer.used; i+=1){
	    t = (kmer) cur_layer.arr[i];
	    //(k-1)-mer, no need to ^luMSB as the head will shift MSB out
	    if(t>=luMSB){
		//insertion
		for(j=0; j<k; j+=1){
		    head = (t>>(j<<1))<<((j+1)<<1);
		    tail = ((1lu<<(j<<1))-1) & t;
		    for(m=0; m<4; m+=1){
			body = m<<(j<<1);
			x = head|body|tail;
			//x is a k-mer
			//add to next_layer if not visited
			if(!HTableSearch(&visited, x)){
			    AListInsert(&next_layer, (void*) x);
			    HTableInsert(&visited, x);
			    //add to hs if is in sample
			    if(!check_sample || isInSampleD1(x, k)){
				hs = AVLAdd(hs, (void*)x, cmpKMer);
			    }
			}
		    }
		}
	    }
	    //k-mer
	    else{
		//deletion
		for(j=0; j<k; j+=1){
		    head = (t>>((j+1)<<1))<<(j<<1);
		    tail = ((1lu<<(j<<1))-1) & t;
		    x = head|tail|luMSB;
		    //x is a (k-1)-mer
		    //add to next_layer if not visited
		    if(!HTableSearch(&visited, x)){
			AListInsert(&next_layer, (void*)x);
			HTableInsert(&visited, x);
		    }
		    
		}
		//substitution
		for(j=1; j<=k; j+=1){
		    head = (t>>(j<<1))<<(j<<1);
		    tail = ((1lu<<((j-1)<<1))-1) & t;
		    for(m=0; m<4; m+=1){
			body = m<<((j-1)<<1);
			x = head|body|tail;
			//x is a k-mer
			//add to next_layer if not visited
			if(!HTableSearch(&visited, x)){
			    AListInsert(&next_layer, (void*)x);
			    HTableInsert(&visited, x);
                            //add to hs if is in sample
			    if(!check_sample || isInSampleD1(x, k)){
				hs = AVLAdd(hs, (void*)x, cmpKMer);
			    }
			}
		    }
		}
	    }//end k-mer
	}//end for each in cur_layer

	AListClear(&cur_layer, NULL);
	AListSwap(&cur_layer, &next_layer);
    }//end for depth from 1 to r

    HTableFree(&visited);
    AListFree(&cur_layer, NULL);
    AListFree(&next_layer, NULL);

    return hs;    
}//end bfsNeighborsInSampleRadius

int hasCollision(AVLNode* hs, AVLNode* ht){
    if(hs == NULL) return 0;
    AVLNode* node = AVLSearch(ht, hs->data, cmpKMer);
    if(node)
	return 1;
    if(hasCollision(hs->left, ht)) return 1;
    if(hasCollision(hs->right, ht)) return 1;
    return 0;
}

int main(int argc, char* argv[]){
    if(argc != 4 || (argv[3][0] != 'w' && argv[3][0] != 's')){
	printf("usage: LSB-statistics.out n r w(hole)|s(ample)\n");
	return 1;
    }

    int k = atoi(argv[1]);
    int r = atoi(argv[2]);
    int check_sample = (argv[3][0] == 'w' ? 0 : 1);

    srand(time(0));

    int i, j, d;
    kmer s, t;
    int share_center[2][3] = {{0, 0, 0}, {0, 0, 0}};
    int ct[2][3] = {{0, 0, 0}, {0, 0, 0}};
    int edit_type;
    int col;
    int col_ct;

    AVLNode *hs, *ht;

    printf("edit\t#col\tcol%%\n");
    for(d=1; d<7; d+=1){
	col_ct = 0;
	for(i=0; i<N; i+=1){
	    s = randomKMer(k);
	    t = randomEdit(s, k, d);
	    
	    
	    hs = bfsNeighborsInSampleRadius(s, k, r, check_sample);
	    ht = bfsNeighborsInSampleRadius(t, k, r, check_sample);
	    
	    col = hasCollision(hs, ht);
	    if(col) col_ct += 1;
		
	    if(d==2 || d==4){
		edit_type = getEditType(s, t, k, d);
		ct[d/2-1][edit_type] += 1;
		if(col){
		    share_center[d/2-1][edit_type] += 1;
		}
	    }
	    
	    AVLFreeTree(hs, NULL);
	    AVLFreeTree(ht, NULL);
	}
	printf("%d\t%d\t%.2f%%\n", d, col_ct, col_ct*100.0/N);
    }

    printf("\nedit\tedit_type\t#\t#col\tcol%%\n");

    for(i=0; i<2; i+=1){
	d = (i<<1) + 2;
	for(j=0; j<2+i; j+=1){
	    printf("%d\t%d+%d*2\t\t%d\t%d\t%.2f%%\n",
		   d, d-(j<<1), j, ct[i][j], share_center[i][j],
		   share_center[i][j]*100.0/ct[i][j]);
	}
    }
    
    return 0;
}
