#ifndef SMITH_WATERMAN_INT_MATRIX_H
#define SMITH_WATERMAN_INT_MATRIX_H

//zero based matrix data structure containing integers
struct Int_Matrix {
    int **matrix;
    int rows;
    int cols;
};

//returns -1 if the allocation is bad. 0 if ok
int int_matrix_new (struct Int_Matrix *m, int rows, int cols);

//-1 if out of bounds, 0 if in bounds. The result is in ret
int int_matrix_get (struct Int_Matrix *m, int *ret, int row, int col);

//-1 if out of bounds, 0 if in bounds
int int_matrix_set (struct Int_Matrix *m, int set, int row, int col);

void int_matrix_print (struct Int_Matrix *m);

void int_matrix_free (struct Int_Matrix *m);

int int_matrix_rows_size (struct Int_Matrix *m);

int int_matrix_cols_size (struct Int_Matrix *m);


#endif
