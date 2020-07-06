/*
 * main.c
 *
 *  Created on: 14 May 2020
 *      Author: avikef
 */

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

double epsilon = 0.00001;
int num = 1;

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

void divideVectorByK(double *vector, double K, int N);
int isLegal(double *vector1, double *vector2, int N);
void randomizeVector(int n, double* vector);

/* divides vector by K */
void divideVectorByK(double *vector, double K, int N){
	int i;

	for(i = 0; i < N; i++){
		*vector = (double)(*vector / K);
		vector++;
	}
}

/* checks for difference using the epsilon constant */
int isLegal(double *vector1, double *vector2, int N){
	int i;
	double result;

	for(i = 0; i < N; i++){
		result = fabs(*vector1 - *vector2);
		if(result > epsilon){
			return 0;
		}
		vector1 = vector1 + 1;
		vector2= vector2 + 1;
	}
	return 1;
}

/* function to randomize a vector in case vector is not given */
void randomizeVector(int n, double* vector){
	int i;

	for(i = 0; i < n; i++){
		*vector = (double)rand();
		vector++;
	}
}

int main(int argc, char* argv[]){
	double *b0, i, j, *input_matrix_row, *result_vector, *tmp, *tmp_values, *tmp_matrix_row;
	FILE *input_matrix_file, *initial_vector_file, *output_vector;
	char *implementation, *matrix_path, *vector_path, *output_file_path;
	int rows, condition, nnz, legal, *tmp_colind, *tmp_rowptr;
	spmat *sparse_matrix;
	/*clock_t start, end;*/
	array_struct *arr_struct;
	node **tmp_nodes_list;
	linked_List_Struct *list_struct;

	/*start = clock();*/

	/* checks for valid size of argc */
	assert(argc == 5 || argc == 4);

	/* check for assertion */
	matrix_path = argv[1];

	/* opening matrix file */
	input_matrix_file = fopen(matrix_path, "rb");
	assert(input_matrix_file != NULL);

	/* read size */
	fread(&rows, sizeof(int), 1, input_matrix_file);
	assert(rows != 0);

	/* places pointer to start of the array */
	fseek(input_matrix_file, sizeof(int), SEEK_CUR);

	/* allocating memory for vector */
	b0 = (double*)malloc(sizeof(double)*rows);
	assert(b0 != NULL);

	/* randomize a vector if one isn't provided */
	if(argc == 4){
		/* this constant is used to differentiate between given\not given vector cases.
		in case vector is given -> argc size is 5.
		else, argc size is.
		*/
		condition = -1;
		srand(time(NULL));

		/* randomize a vector */
		randomizeVector(rows, b0);
		assert(b0 != NULL);
	}
	/* case we are given a vector */
	else{
		condition = 0;
		vector_path = argv[2];

		initial_vector_file = fopen(vector_path, "rb");
		assert(initial_vector_file != NULL);

		fseek(initial_vector_file, sizeof(int)*2, SEEK_SET);
		fread(b0, sizeof(double), rows, initial_vector_file);
		assert(b0 != NULL);

		fclose(initial_vector_file);
	}

	output_file_path = argv[3 + condition];
	output_vector = fopen(output_file_path, "wb");
	assert(output_vector != NULL);

	/* reading implementation choice */
	implementation = argv[4 + condition];
	assert( (strcmp(implementation, "-list") == 0) || (strcmp(implementation, "-array") == 0) );

	input_matrix_row = (double*)malloc(sizeof(double)*(rows));
	assert(input_matrix_row != NULL);

	/* case implementation is array */
	if(strcmp(implementation, "-array") == 0){
		nnz = 0;

		/* counting non zero values */
		for(i = 0; i < rows; i++){
			fread(input_matrix_row, sizeof(double), rows, input_matrix_file);
			assert(input_matrix_row != NULL);
			
			/* saving start index of array to avoid resetting it using -= N */
			tmp_matrix_row = input_matrix_row;

			/* for each row count non zero values */
			for(j = 0; j < rows; j++){
				if(*input_matrix_row != 0){
					nnz++;
				}
				input_matrix_row++;
			}
			input_matrix_row = tmp_matrix_row;
		}

		/* sets pointer on matrix file at the beginning again */
		fseek(input_matrix_file, sizeof(int)*2, SEEK_SET);

		/* allocating spmat using array implementation */
		sparse_matrix = spmat_allocate_array(rows, nnz);

		/* declaring a struct to save needed data to represent a sparse matrix
			using array implementation.
		*/
		arr_struct = (array_struct*)sparse_matrix->private;

		/* saving start indexes for arr_struct fields */
		tmp_colind = arr_struct->colind;
		tmp_values = arr_struct->values;
		tmp_rowptr = arr_struct->rowptr;

		/* building sparse matrix */
		for(i = 0; i < rows; i++){
			fread(input_matrix_row, sizeof(double), rows, input_matrix_file);
			assert(input_matrix_row != NULL);

			sparse_matrix->add_row(sparse_matrix, input_matrix_row, i);
		}

		/* resetting indexes */
		arr_struct->colind = tmp_colind;
		arr_struct->values = tmp_values;
		arr_struct->rowptr = tmp_rowptr;
	}
	/* case implementation is list */
	else{

		/* allocating spmat using list implementation */
		sparse_matrix = spmat_allocate_list(rows);

		/* declaring a struct to save needed data to represent a sparse matrix
			using list implementation.
		*/
		list_struct = (linked_List_Struct*)sparse_matrix->private;
		
		/* saving start indexes for list_struct field */
		tmp_nodes_list = list_struct->node_list;

		/* building sparse matrix */
		for(i = 0; i < rows; i++){
			fread(input_matrix_row, sizeof(double), rows, input_matrix_file);
			assert(input_matrix_row != NULL);

			sparse_matrix->add_row(sparse_matrix, input_matrix_row, i);
		}

		/* resetting indexes */
		list_struct->node_list = tmp_nodes_list;
	}

	/* declaring result vector to output after calculating */
	result_vector = (double*)malloc(rows*sizeof(double));
	assert(result_vector != NULL);

	legal = 0;
	while (legal == 0){
		/* using proper implementation function within spmat module
		to mult and normalize current result vector */
		sparse_matrix->mult(sparse_matrix, b0, result_vector);
		
		/* checking for delta < epsilon */
		legal = isLegal(result_vector, b0, rows);

		tmp = b0;
		b0 = result_vector;
		result_vector = tmp;
	}

	/* writing result to output */
	fwrite(&num, sizeof(int), 1, output_vector);
	fwrite(&rows, sizeof(int), 1, output_vector);
	fwrite(result_vector, sizeof(double), rows, output_vector);

	/* closing files used */
	fclose(output_vector);
	fclose(input_matrix_file);

	/* freeing arrays\spmat */
	free(b0);
	free(input_matrix_row);
	free(result_vector);
	sparse_matrix->free(sparse_matrix);

	/*end = clock();
	printf("%f", (float)(end - start) / CLOCKS_PER_SEC);*/

	return 1;
}
