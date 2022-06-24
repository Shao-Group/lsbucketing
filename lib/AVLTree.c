#include "AVLTree.h"

static AVLNode* rotateLeft(AVLNode* root){
	AVLNode* x = root->right;
	root->right = x->left;
	x->left = root;
	int l = AVLNODEHEIGHT(root->left);
	int r = AVLNODEHEIGHT(root->right);
	l = root->height = (l>r?l:r)+1;
	r = AVLNODEHEIGHT(x->right);
	x->height = (l>r?l:r)+1;

	root->size = AVLNODESIZE(root->left)+AVLNODESIZE(root->right)+1;
	x->size = AVLNODESIZE(x->right)+root->size+1;
	return x;
}

static AVLNode* rotateRight(AVLNode* root){
	AVLNode* x = root->left;
	root->left = x->right;
	x->right = root;
	int l = AVLNODEHEIGHT(root->left);
	int r = AVLNODEHEIGHT(root->right);
	r = root->height = (l>r?l:r)+1;
	l = AVLNODEHEIGHT(x->left);
	x->height = (l>r?l:r)+1;

	root->size = AVLNODESIZE(root->left)+AVLNODESIZE(root->right)+1;
	x->size = AVLNODESIZE(x->left)+root->size+1;
	return x;
}

static inline void AVLFreeNode(AVLNode* x, void (*freeData)(void*)){
    if(freeData){
	freeData(x->data);
    }
    free(x);
}

void AVLFreeTree(AVLNode *n, void (*freeData)(void*)){
	if(n==NULL) return;
	AVLFreeTree(n->left, freeData);
	AVLFreeTree(n->right, freeData);
	if(freeData){
	    freeData(n->data);
	}
	free(n);
}

AVLNode* AVLAdd(AVLNode* root, void* d, int (*cmp)(const void*, const void*)){
	AVLNode* x = malloc(sizeof(AVLNode));
	x->left = x->right = NULL;
	x->data = d;
	x->height = 1;
	x->size = 1;
	
	if(root == NULL){
		return x;
	}else{
		AVLNode* path[root->height+1];
		int current=0;
		path[current]=root;
		
		while(path[current]!=NULL){
			if(cmp(path[current]->data, d)>=0){
				path[current+1] = path[current]->left;
			}else{
				path[current+1] = path[current]->right;
			}
			current++;
		}
				
		if(cmp(path[--current]->data, d)>=0) path[current]->left = x;
		else path[current] ->right = x;
		path[current+1] = x;
		
		while(current>=0){
			AVLNode* c = path[current];
			int lh = AVLNODEHEIGHT(c->left);
			int rh = AVLNODEHEIGHT(c->right);
			
			if(lh-rh>1){//L, inserted in left, path[current]->left == path[current+1]
				if(path[current+1]->right==path[current+2]){//LR
					c->left = rotateLeft(path[current+1]);
				}
				c = rotateRight(c);
				if(current>0){
					if(path[current-1]->left==path[current]) path[current-1]->left = c;
					else path[current-1]->right = c;
				}else{
					return c;
				}
				//break;
			}else if(lh-rh<-1){//R
				if(path[current+1]->left==path[current+2]){//RL
					c->right = rotateRight(path[current+1]);
				}
				c = rotateLeft(c);
				if(current>0){
					if(path[current-1]->left==path[current]) path[current-1]->left = c;
					else path[current-1]->right = c;
				}else{
					return c;
				}
				//break;
			}else{
				c->height = (lh>rh?lh:rh)+1;
				c->size = AVLNODESIZE(c->left)+AVLNODESIZE(c->right)+1;
			}
			current--;
		}
		return path[0];
	}
}

AVLNode* AVLSearch(AVLNode* root, const void* d, int (*cmp)(const void*, const void*)){
	while(root!=NULL){
		int result = cmp(root->data, d);
		if(result==0) return root;
		else if(result<0) root = root->right;
		else root=root->left;
	}
	
	return NULL;
}

AVLNode* AVLDeleteNode(AVLNode* root, AVLNode* x, int (*cmp)(const void*, const void*), void (*freeData)(void*)){
	AVLNode* path[root->height];
	int current = 0;
	path[current] = root;
	while(path[current]!=NULL){
		int result = cmp(path[current]->data, x->data);
		if(path[current]==x) break;
		else if(result>=0){
			//path[current]=path[current++]->left;			
			path[current+1] = path[current]->left;
		}else{
			path[current+1] = path[current]->right;
		}
		current++;
	}
	
	if(path[current]==NULL) return root;//no found, unchanged
	if(x->left==NULL){
		if(current==0){
			root = root->right;
			AVLFreeNode(x, freeData);
			return root;
		}else if(path[current-1]->left == x){
			path[current-1]->left = x->right;
		}else{
			path[current-1]->right = x->right;
		}
	}else if(x->right==NULL){
		if(current==0){
			root = root->left;
			AVLFreeNode(x, freeData);
			return root;
		}else if(path[current-1]->left == x){
			path[current-1]->left = x->left;
		}else{
			path[current-1]->right = x->left;
		}
		
	}else{
		int prevx = current-1;
		path[++current]=x->left;
		while(path[current]->right!=NULL){
			path[current+1]=path[current]->right;
			current++;
		}
		
		if(path[current-1]->right==path[current]) path[current-1]->right = path[current]->left;
		else path[current-1]->left = path[current]->left;
		path[current]->left = x->left;
		path[current]->right = x->right;
		if(prevx>=0){
			if(path[prevx]->left==x) path[prevx]->left = path[current];
			else path[prevx]->right = path[current];
		}
		path[prevx+1]=path[current];
	}
	AVLFreeNode(x, freeData);
	current--;

	while(current>=0){
		AVLNode* c = path[current];
		int lh = AVLNODEHEIGHT(c->left);
		int rh = AVLNODEHEIGHT(c->right);
		
		if(lh-rh>1){//R, deleted in right, path[current] has left
			AVLNode* l = c->left;
			if(AVLNODEHEIGHT(l->right)>AVLNODEHEIGHT(l->left)){//LR
				c->left = rotateLeft(l);
			}
			c = rotateRight(c);
			if(current>0){
				if(path[current-1]->left==path[current]) path[current-1]->left = c;
				else path[current-1]->right = c;
			}else{
				return c;
			}
		}else if(lh-rh<-1){//L
			AVLNode* r = c->right;
			if(AVLNODEHEIGHT(r->left)>AVLNODEHEIGHT(r->right)){//LR
				c->right = rotateRight(r);
			}
			c = rotateLeft(c);
			if(current>0){
				if(path[current-1]->left==path[current]) path[current-1]->left = c;
				else path[current-1]->right = c;
			}else{
				return c;
			}
		}else{
			c->height = (lh>rh?lh:rh)+1;
			c->size = AVLNODESIZE(c->left)+AVLNODESIZE(c->right)+1;
		}
		current--;
	}

	return path[0];
}

AVLNode* AVLDeleteValue(AVLNode* root, void* d, int (*cmp)(const void*, const void*), void (*freeData)(void*)){
	AVLNode* x = AVLSearch(root, d, cmp);
	while(x!=NULL){
		root = AVLDeleteNode(root, x, cmp, freeData);
		x = AVLSearch(root, d, cmp);
	}
	return root;
}

void AVLPrint(const AVLNode* root, int depth, int (*getKey)(const void*)){
	if(root==NULL) return;
	int i=0;
	for(;i<depth;i++) printf("-");
	printf("(%d %d %d)\n", getKey(root->data), root->height, root->size);
	AVLPrint(root->left, depth+1, getKey);
	AVLPrint(root->right, depth+1, getKey);
}
