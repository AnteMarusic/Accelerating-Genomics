#include "int_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include "utility.h"
#include <limits.h>
#define NEGATIVE_INFINITY INT_MIN  // Minimum value for an int

// Allocate the matrix (s1 in the rows and s2 in the cols)



//returns -1 if the allocation is bad. 0 if ok
int int_matrix_new (struct Int_Matrix *m, int rows, int cols) {
    int i, j;
    if (m == NULL) {
        PRINT("the matrix is NULL");
        return -1;
    }
    m -> matrix = (int**) malloc((rows) * sizeof(int*));
    if (m == NULL) {
    // Memory allocation failed
        PRINT("allocation failed");
        return -1;
    }
    for (i = 0; i < rows; i++) {
        m -> matrix[i] = (int*)malloc((cols) * sizeof(int));
        if (m -> matrix[i] == NULL) {
            // Memory allocation failed
            // Cleanup previously allocated memory
            for (j = 0; j < i; j++) {
                free(m -> matrix[j]);
            }
            free(m);
            PRINT("allocation failed");
            return -1;
        }
    }
    m -> rows = rows;
    m -> cols = cols;
    return 0;
}


//-1 if out of bounds, 0 if in bounds
int int_matrix_get (struct Int_Matrix *m, int *ret, int row, int col) {
    if (row >= m -> rows || col >= m -> cols || row < 0 || col < 0) {
        PRINT("out of bounds get %d, %d\n", row, col);
        return -1;
    }
    *ret = m -> matrix[row][col];
    return 0;
}

//-1 if out of bounds, 0 if in bounds
int int_matrix_set (struct Int_Matrix *m, int set, int row, int col) {
    if (row >= m -> rows || col >= m -> cols || row < 0 || col < 0) {
        PRINT("out of bounds set %d, %d\n", row, col);
        return -1;
    }
    m -> matrix[row][col] = set;
    return 0;
}

void int_matrix_print (struct Int_Matrix *m) {
    int i,j;
    for (i = 0; i < 50; i ++) {
        printf ("_");
    }
    for (i = 0; i < m -> rows; i++) {
        printf("\n");
        for (j = 0; j < m -> cols; j++) {
            (m -> matrix[i][j] == NEGATIVE_INFINITY) ? printf("—∞ ") : printf("%d ", m -> matrix[i][j]);
        }
    }
    printf("\n");
    for (i = 0; i < 50; i ++) {
        printf ("_");
    }
    printf ("\n");
}

void int_matrix_free (struct Int_Matrix *m) {
    int i;
    for (i = 0; i < m -> rows; i++) {
        free(m -> matrix[i]);
    }
    free(m -> matrix);
}

//returns -1 in case the matrix is NULL
int int_matrix_rows_size (struct Int_Matrix *m) {
    if (m != NULL) {
        return m -> rows;
    }
    else {
        return -1;
    }
}

int int_matrix_cols_size (struct Int_Matrix *m) {
    if (m != NULL) {
        return m -> cols;
    }
    else {
        return -1;
    }
}
