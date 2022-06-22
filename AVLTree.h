/*
  AVL tree for general struct.
  By: Ke@UWM
  Last edited: 12/19/2016
*/

#ifndef _AVLTREE_H
#define _AVLTREE_H 1

#include <stdio.h>
#include <stdlib.h>
//#include <sqlite3.h>
//#include <string.h>
//#include <math.h>
//#include <time.h>

typedef struct AVLNode{
	struct AVLNode* left;
	struct AVLNode* right;
	
	int size;
	int height;
	void* data;
} AVLNode;

#define AVLNODEHEIGHT(node) (node==NULL?0:node->height)
#define AVLNODESIZE(node) (node==NULL?0:node->size)

//provided by user, needs to guarantee no duplicate key
//static int cmp(const void* a, const void* b);

//free the whole tree, data parts are freed by the provided freeData function
void AVLFreeTree(AVLNode *n, void (*freeData)(void*));

//add data d to tree rooted at root
AVLNode* AVLAdd(AVLNode* root, void* d, int (*cmp)(const void*, const void*));

//search in tree rooted at root for data d
AVLNode* AVLSearch(AVLNode* root, const void* d, int (*cmp)(const void*, const void*));


//delete and free a given node in the tree
AVLNode* AVLDeleteNode(AVLNode* root, AVLNode* x, int (*cmp)(const void*, const void*), void (*freeData)(void*));

//delete a given data from the tree
AVLNode* AVLDeleteValue(AVLNode* root, void* d, int (*cmp)(const void*, const void*), void (*freeData)(void*));

//print the tree
void AVLPrint(const AVLNode* root, int depth, int (*getKey)(const void*));

#endif //AVLTree.h
