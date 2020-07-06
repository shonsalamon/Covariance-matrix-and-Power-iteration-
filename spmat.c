/*
 * spmat.c
 *
 *  Created on: 14 May 2020
 *      Author: avikef
 */

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* struct to represent nodes for -list */
typedef struct node {
	int column;
	double value;
    struct node *next;
} node;

/* struct to represent arrays needed to contain data for -array*/
typedef struct arrayStruct {
	double *values;
	int *rowptr, *colind;
} array_struct;

/* struct to represent list of linked list - each for every row of given matrix */
typedef struct linkedListStruct {
	node **node_list;
} linked_List_Struct;

spmat* spmat_allocate_list(int n);
void linked_list_add_row(struct _spmat *A, const double *row, int i);
void free_linked_list(struct _spmat *A);
void linked_list_mult(const struct _spmat *A, const double *v, double *result);

spmat* spmat_allocate_array(int n, int nnz);
void array_add_row(struct _spmat *A, const double *row, int i);
void free_array(struct _spmat *A);
void array_mult(const struct _spmat *A, const double *v, double *result);

void divideVectorByK(double *vector, double K, int N);

/* function to allocate list implementation for spmat */
spmat* spmat_allocate_list(int n){
	spmat *s;
	linked_List_Struct *list_struct;

	s = (spmat*)malloc(sizeof(spmat));
	assert(s != NULL);

	list_struct = (linked_List_Struct*)malloc(sizeof(linked_List_Struct));
	assert(list_struct != NULL);

	s->n = n;
	s->add_row = &linked_list_add_row;
	s->mult = &linked_list_mult;
	s->free = &free_linked_list;
	s->private = (linked_List_Struct*)list_struct;
	assert(s->private != NULL);

	list_struct->node_list = (node**)malloc(sizeof(node*)*(n+1));
	assert(list_struct->node_list != NULL);
	return s;
}

/* function to allocate array implementation for spmat */
spmat* spmat_allocate_array(int n, int nnz){
	spmat *s;
	array_struct *arr_struct;

	arr_struct = (array_struct*)malloc(sizeof(array_struct));
	assert(arr_struct != NULL);

	s = (spmat*)malloc(sizeof(spmat));
	assert(s != NULL);

	s->n = n;
	s->add_row = &array_add_row;
	s->free = &free_array;
	s->mult = &array_mult;

	s->private = (array_struct*)arr_struct;
	assert(s->private != NULL);

	/* values */
	arr_struct->values = (double*)malloc(sizeof(double)*(nnz));
	assert(arr_struct->values != NULL);
	/* colind */
	arr_struct->colind = (int*)malloc(sizeof(int)*(nnz));
	assert(arr_struct->colind != NULL);
	/* rowptr */
	arr_struct->rowptr = (int*)malloc(sizeof(int)*(n + 1));
	assert(arr_struct->rowptr != NULL);

	/* first cell of rowptr always = 0 */
	((int*)arr_struct->rowptr)[0] = 0;

	return s;
}

void linked_list_add_row(struct _spmat *A, const double *row, int i){
	int column, isFirst;
	node *new_node, *current_node;
	linked_List_Struct *list_struct;

	isFirst = 1;
	/* casting A->private into proper list struct */
	list_struct = (linked_List_Struct*)A->private;

	assert(i >= 0);

	/* for each row of list_struct checks if there are nodes.
	   if there isn't -> check next row.
	   else, declaring a new node, copying the proper values into the node,
	   and adding it to the row.
	*/
	for(column = 0; column < A->n; column++){
		if(*row != 0){

			new_node = (node*)malloc(sizeof(node));
			assert(new_node != NULL);

			new_node->column = column;
			new_node->value = *row;
			new_node->next = NULL;

			if(isFirst == 1){
				*(list_struct->node_list) = new_node;
				isFirst = 0;
				current_node = *(list_struct->node_list);
			}
			else{
				current_node->next = new_node;
				current_node = current_node->next;
			}
		}
		row++;
	}
	/* promote list_struct to next row */
	list_struct->node_list++;
}

void array_add_row(struct _spmat *A, const double *row, int i){
	int col, count;
	array_struct *arr_struct;

	assert(i >= 0);

	/* casting A->private into proper arr struct */
	arr_struct = (array_struct*)A->private;
	
	/* needs the count of values\cols added until row i */
	count = *(arr_struct->rowptr);

	/* for each cell in row, in case row != 0,
	   adding proper values to values and colind,
	   and promoting count.
	   after finishing with every row,
	   rowptr[i+1] = count since that is the logic.
	   
	   rowptr[i+1] - rowptr[i] = number of non-zero values in row i.
	 */
	for(col = 0; col < A->n; col++){
		if(*row != 0){

			*(arr_struct->values) = *row;
			*(arr_struct->colind) = col;

			count++;
			(arr_struct->values)++;
			(arr_struct->colind)++;

		}
		row++;
	}
	arr_struct->rowptr += 1;
	*(arr_struct->rowptr) = count;
}

/* goes iteratively through every row of list struct and frees every node */
void free_linked_list(struct _spmat *A){
	node *current, *next_node, **arr;
	int i;
	linked_List_Struct *list_struct;

	/* casting A->private into proper list struct */
	list_struct = (linked_List_Struct*)A->private;

	arr = list_struct->node_list;

	for(i = 0; i < A->n; i++){
		current = *arr;
		if(current != NULL){
			next_node = current->next;
			while(next_node != NULL){
				free(current);
				current = next_node;
				next_node = next_node->next;
			}
			free(current);
		}
		arr++;
	}
	/* frees used structs */
	free(list_struct->node_list);
	free(list_struct);
	free(A);
}

void free_array(struct _spmat *A){
	array_struct *arr_struct;
	/* casting A->private into proper arr struct */
	arr_struct = (array_struct*)A->private;

	/* frees used structs */
	free(arr_struct->values);
	free(arr_struct->colind);
	free(arr_struct->rowptr);
	free(arr_struct);
	free(A);
}

void linked_list_mult(const struct _spmat *A, const double *v, double *result){

	int i;
	node *current, **arr;
	double dot_product, denominator, *tmp_result;
	linked_List_Struct *list_struct;

	/* casting A->private into proper list struct */
	list_struct = (linked_List_Struct*)A->private;
	arr = list_struct->node_list;

	denominator = 0;
	dot_product = 0;
	tmp_result = result;

	current = *arr;

	/* 
	for each row in list struct multiplying node by using it's value\col
	with the proper index of v and adding result to dot_product.
	*/
	for(i = 0; i < A->n; i++){
		while(current != NULL){
			dot_product += (double)((current->value) * v[current->column]);
			current = current->next;
		}
		arr++;
		current = *arr;
		*result = dot_product;
		/* calculating norm of result vector inside mult */
		denominator += (dot_product * dot_product);
		result++;
		dot_product = 0.0;
	}
	/* normalizing vector */
	divideVectorByK(tmp_result, sqrt(denominator), A->n);
}

void array_mult(const struct _spmat *A, const double *v, double *result){
	int i, j, *rowptr, *colind;
	double *values, denominator, dot_product, *tmp_result;
	array_struct *arr_struct;

	denominator = 0;
	tmp_result = result;

	/* casting A->private into proper arr struct */
	arr_struct = (array_struct*)A->private;
	
	/* using pointers to the struct to keep the indexes at its beginning */
	values = (double*)arr_struct->values;
	colind = (int*)arr_struct->colind;
	rowptr = (int*)arr_struct->rowptr;

	/* 
	for each j in [rowptr[i], rowptr[i+1]],
	there are (rowptr[i+1] - rowptr[i]) non-zero values in row i.
	so for each j, we calculate dot product of current line using values\colind arrays contents.
	*/
	for (i = 0 ; i < A->n ; i++){
		dot_product = 0;

		for(j = *rowptr; j < *(rowptr + 1); j++){
			dot_product += (double)((*values) * v[*colind]);
			values++;
			colind++;
		}
		rowptr++;
		*result = dot_product;
		/* calculating norm of result vector inside mult */
		denominator += (dot_product * dot_product);
		result++;
	}
	/* normalizing vector */
	divideVectorByK(tmp_result, sqrt(denominator), A->n);
}
