

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char* argv[]){
	FILE *file, *output_file;
	char *filepath, *output_path;
	int cols, rows, m, n, test_size, total_size, i, j, k, test;
	double *matrix, *cov_row, *ptr_cov_row, *ptr_skip_row, *read_pointer, *current_row;
	double sum_row, avg_row, dot_product;

	assert(argc == 3);

	filepath = argv[1];


	file = fopen(filepath, "rb");

	m = fread(&cols, sizeof(int), 1, file);
	assert (m == 1);
	n = fread(&rows, sizeof(int), 1, file);
	assert (n == 1);

	total_size = rows*cols*sizeof(double);

	matrix = malloc(total_size);

	test_size = fread(matrix, sizeof(double), rows*cols, file);
	assert (test_size == rows*cols);

	output_path = argv[2];

	output_file = fopen(output_path, "wb");

	read_pointer = matrix;
	current_row = matrix;
	sum_row = 0;
/* Normalizing the matrix*/
	for (i = 0; i < rows; i++){
		for (j = 0; j < cols; j++){
			sum_row += *read_pointer;
			read_pointer++;
		}

		avg_row = sum_row / cols;
		for (j = 0; j < cols; j++){
			*current_row = *current_row - avg_row;
			current_row++;
		}
		sum_row = 0;
	}

	read_pointer = matrix;
	current_row = matrix;
	ptr_skip_row = matrix;
	dot_product = 0.0;
	cov_row = malloc(sizeof(double)*rows);
	ptr_cov_row = cov_row;

	fwrite(&rows, sizeof(int), 1, output_file);
	fwrite(&rows, sizeof(int), 1, output_file);
/* calculating the cov matrix.*/
	for (i = 0; i < rows; i++){
		for (j = 0; j < rows; j++){
			for (k = 0; k < cols; k++){
				dot_product = dot_product + (*read_pointer) * (*(ptr_skip_row));
				read_pointer++;
				ptr_skip_row++;
			}
			*ptr_cov_row = dot_product;
			ptr_cov_row++;
			dot_product = 0.0;
			read_pointer = current_row;
		}
		ptr_skip_row = matrix;
		current_row += cols;
		read_pointer = current_row;

		ptr_cov_row -= rows;
		/* writing a row to the matrix.*/
		test = fwrite(ptr_cov_row, sizeof(double), rows, output_file);
		assert (test == rows);
	}

	free(matrix);
	free(cov_row);
	fclose(file);
	fclose(output_file);
	return 1;
}
