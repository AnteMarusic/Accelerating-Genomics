#include "algo.h"
#include "utility.h"
#include "list.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NO_CHOICE 3
#define SX 2
#define ABOVE 1
#define DIAG 0


void matrix_print (int** matrix, int rows, int cols);
void print_res (struct Result res);
struct Result compute_aligment (int **matrix, int **dmatrix, int res_len, int i_max, int j_max, const char *s1, const char *s2);


//todo Not necessary to allocate space for the first row and col in both matrixes
//todo return all possible optimal assignments

struct WholeResult* align(char* s1, char* s2, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
    if (s1 == NULL || s2 == NULL) {
        return NULL;
    }
    int s1_len = strlen(s1);
    int s2_len = strlen(s2);
    PRINT("s1_len: %d, s2_len: %d\n",  s1_len, s2_len);
    int i, j, val, temp, dir, max, k, z, list_len;
    struct Node *head;
    struct Node *cur;
    list_init(&head);
    int **matrix,**dmatrix;

    // Allocate the matrix (s1 in the rows and s2 in the cols)
    matrix = (int**)malloc((s2_len + 1) * sizeof(int*));
    dmatrix = (int**)malloc((s2_len + 1) * sizeof(int*));
    if (matrix == NULL || dmatrix == NULL) {
        // Memory allocation failed
        return NULL;
    }

    for (i = 0; i < s2_len + 1; i++) {
        matrix[i] = (int*)malloc((s1_len + 1) * sizeof(int));
        if (matrix[i] == NULL) {
            // Memory allocation failed
            // Cleanup previously allocated memory
            for (j = 0; j < i; j++) {
                free(matrix[j]);
                free(dmatrix[j]);
            }
            free(matrix);
            free(dmatrix);
            return NULL;
        }
        dmatrix[i] = (int*)malloc((s1_len + 1) * sizeof(int));
        if (dmatrix[i] == NULL) {
            // Memory allocation failed
            // Cleanup previously allocated memory
            free(matrix[i]);
            for (j = 0; j <= i; j++) {
                free(matrix[j]);
                free(dmatrix[j]);
            }
            free(matrix);
            free(dmatrix);
            return NULL;
        }
    }

    //Initialize matrix
    for (i = 0; i < s2_len + 1; i++) {
        for (j = 0; j < s1_len + 1; j++) {
            matrix[i][j] = 0;
            dmatrix[i][j] = NO_CHOICE;
        }
    }
    PRINT("matrix:\n");
    matrix_print(matrix, s2_len + 1, s1_len + 1);
    PRINT("dmatrix:\n");
    matrix_print(dmatrix, s2_len + 1, s1_len + 1);

    //Dynamic Programming
    max = 0;
    for (i = 1; i < s2_len + 1; i++) {
        for (j = 1; j < s1_len + 1; j++) {
            // Calculate scores for neighboring cells and update the current cell
            val = matrix[i][j];
            dir = NO_CHOICE;
            if (s1[j - 1] == s2[i - 1]) {
                temp = matrix[i-1][j-1] + match_score;
                if (temp > val) {
                    val = temp;
                    dir = DIAG;//register in the backtraking
                }

            }
            else {
                temp = matrix[i-1][j-1] + mismatch_score;
                if (temp > val) {
                    val = temp;
                    dir = DIAG; //register in the backtraking
                }
            }
            temp = matrix[i][j-1] + gap_extend_score;
            if (temp > val) {
                val = temp;
                dir = SX;
            }
            temp = matrix[i-1][j] + gap_extend_score;
            if (temp > val) {
                val = temp;
                dir = ABOVE;
            }
            if (val == max) {
                add(&head, create_node(i,j));
            }
            if (val > max) {
                max = val;
                delete_all(&head);
                add(&head, create_node(i,j));
            }
            matrix[i][j] = val;
            dmatrix[i][j] = dir;
        }
    }
    PRINT("matrix:\n");
    matrix_print(matrix, s2_len + 1, s1_len + 1);
    PRINT("dmatrix:\n");
    matrix_print(dmatrix, s2_len + 1, s1_len + 1);

    // Backtracking
    PRINT("max_val: %d\n" , val);
    print_list(head);

    k = s1_len > s2_len ? s1_len : s2_len;

    list_len = get_length(head);
    struct Result *out = (struct Result *)malloc (list_len*sizeof (struct Result));
    for (z = 0; z < list_len; z ++) {
        cur = get(head, z);
        out[z] = compute_aligment (matrix, dmatrix, k, cur->row, cur->col, s1, s2);
        print_res(out[z]);
    }


    // Cleanup allocated memory
    for (i = 0; i < s2_len + 1; i++) {
        free(matrix[i]);
        free(dmatrix[i]);
    }
    free(matrix);
    free(dmatrix);
    //free the list
    delete_all(&head);

    // Output
    struct WholeResult *w_res = (struct WholeResult*) malloc(sizeof (struct WholeResult));
    w_res->res = out;
    w_res->len = list_len;
    return w_res;
}

void matrix_print (int** matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        PRINT("\n");
        for (j = 0; j < cols; j++) {
            PRINT(" %d ", matrix[i][j]);
        }
    }
    PRINT("\n");
}

struct Result compute_aligment (int **matrix, int **dmatrix, int res_len, int i_max, int j_max, const char *s1, const char *s2) {
    int q;
    char *out_1, *out_2;
    int i = i_max, j = j_max, k = res_len;
    struct Result res = {"error", "error"};
    out_1 = (char*) malloc((res_len) * sizeof(char));
    if (out_1 == NULL) {
        /*
        for (q = 0; q < rows + 1; q++) {
            free(matrix[q]);
            free(dmatrix[q]);
        }
        free(matrix);
        free(dmatrix);
         */

        return res;
    }
    out_2 = (char*) malloc((res_len) * sizeof(char));
    if (out_2 == NULL) {
        /*
        for (q = 0; q < res_len + 1; q++) {  //WROOOOONG is not res len but size of s1 add global variables
            free(matrix[q]);
            free(dmatrix[q]);
        }
        free(matrix);
        free(dmatrix);
        free(out_1);
         */

        return res;
    }

    for (q=0; q < res_len; q ++) {
        out_1[q] = '_';
        out_2[q] = '_';
    }

    while (dmatrix[i][j] != NO_CHOICE) {
        if (dmatrix[i][j] == DIAG) {
            out_1[k - 1] = s1[j-1];
            out_2[k - 1] = s2[i-1];
            i--;
            j--;
        }
        else if (dmatrix[i][j] == SX) {
            out_1[k - 1] = s1[j-1];
            out_2[k - 1] = '_';
            j--;
        }
        else {
            out_1[k - 1] = '_';
            out_2[k - 1] = s2[i-1];
            i--;
        }
        k--;
        PRINT("out_1: %s\n", out_1);
        PRINT("out_2: %s\n", out_1);
    }
    res.s1 = out_1;
    res.s2 = out_2;
    return res;
}

void print_res (struct Result res) {
    printf("First string: %s\n", res.s1);
    printf("Second string: %s\n", res.s2);
}

void print_whole_res (struct WholeResult *w_res) {
    for (int i = 0; i < w_res->len; i ++) {
        printf ("alignment n° %d\n", i);
        print_res(w_res->res[i]);
    }
}