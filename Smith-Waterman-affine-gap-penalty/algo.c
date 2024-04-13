#include "algo.h"
#include "utility.h"
#include "int_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

//use the infinity concept. Actually it should add something to infinity, I never subtract
//or am I wrong?

#define NEGATIVE_INFINITY INT_MIN  // Minimum value for an int

int max_of_two_integers(int a, int b) {
    return (a > b) ? a : b;
}

int sum_with_infinity (int a, int b) {
    return (a == NEGATIVE_INFINITY || b == NEGATIVE_INFINITY) ? NEGATIVE_INFINITY : a + b;
}

int max_of_four_integers (int a, int b, int c, int d) {
    return max_of_two_integers (max_of_two_integers (a, b), max_of_two_integers(c, d));
}



//todo Not necessary to allocate space for the first row and col in both matrixes
//todo return all possible optimal assignments

//returns -1 if the strings are NULL or the allocation of the
// matrixes failed, else the result is 0 or positive
int align(char* s1, char* s2, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
    if (s1 == NULL || s2 == NULL) {
        return -1;
    }
    int s1_len = strlen(s1);
    int s2_len = strlen(s2);
    PRINT("s1_len: %d, s2_len: %d\n",  s1_len, s2_len);
    int flag, i, j, max, temp1, temp2, temp3;

    struct Int_Matrix *p, *q, *d;
    p = (struct Int_Matrix *) malloc (sizeof(struct Int_Matrix));
    q = (struct Int_Matrix *) malloc (sizeof(struct Int_Matrix));
    d = (struct Int_Matrix *) malloc (sizeof(struct Int_Matrix));

    // Allocate the matrixes (s2 in the rows and s1 in the cols)
    flag = int_matrix_new (p, s1_len + 1, s2_len + 1);
    if (flag == -1) {    //returns -1 if the allocation is bad. 0 if ok
        return -1;
    }
    flag = int_matrix_new (q, s1_len + 1, s2_len + 1);
    if (flag == -1) {    //returns -1 if the allocation is bad. 0 if ok
        return -1;
    }
    flag = int_matrix_new (d, s1_len + 1, s2_len + 1);
    if (flag == -1) {    //returns -1 if the allocation is bad. 0 if ok
        return -1;
    }

    //Initialize matrixes
    for (i = 0; i < s1_len + 1; i ++) {
        for (j = 0; j < s2_len + 1; j ++) {
            if (int_matrix_set (p, 0, i, j) == -1) return -1;
            if (int_matrix_set (q, 0, i, j) == -1) return -1;
            if (int_matrix_set (d, 0, i, j) == -1) return -1;
        }
    }

    int_matrix_print (p);
    int_matrix_print (q);
    int_matrix_print (d);

    //p has —∞ in the first row and d has all zeros
    for (i = 0; i < s2_len + 1; i++) {
        if (int_matrix_set (p, NEGATIVE_INFINITY, 0, i) == -1) return -1;
    }
    for (i = 0; i < s1_len + 1; i++) {
        if (int_matrix_set (q, NEGATIVE_INFINITY, i, 0) == -1) return -1;
    }

    int_matrix_print (p);
    int_matrix_print (q);
    int_matrix_print (d);

    //Dynamic Programming
    max = 0;
    for (i = 1; i < s1_len + 1; i++) {
        for (j = 1; j < s2_len + 1; j++) {
            // Calculate scores for neighboring cells and update the current cell
            if (int_matrix_get (d, &temp1, i - 1, j) == -1) return -1;
            if (int_matrix_get (p, &temp2, i - 1, j) == -1) return -1;
            int_matrix_set (p, max_of_two_integers(sum_with_infinity(temp1, gap_open_score + gap_extend_score), sum_with_infinity(temp2, gap_extend_score)), i, j);

            if (int_matrix_get (d, &temp1, i, j - 1) == -1) return -1;
            if (int_matrix_get (q, &temp2, i, j - 1) == -1) return -1;
            int_matrix_set (q, max_of_two_integers(sum_with_infinity(temp1, gap_open_score + gap_extend_score), sum_with_infinity(temp2, gap_extend_score)), i, j);

            if (int_matrix_get (p, &temp1, i, j) == -1) return -1;
            if (int_matrix_get (q, &temp2, i, j) == -1) return -1;
            if (int_matrix_get (d, &temp3, i - 1, j - 1) == -1) return -1;
            temp3 = s1[i - 1] == s2[j - 1] ? temp3 + match_score : temp3 + mismatch_score;
            int_matrix_set (d, max_of_four_integers(temp1, temp2, temp3, 0), i, j);

            if (int_matrix_get (d, &temp1, i, j) == -1) return -1;
            max = temp1 > max ? temp1 : max;
        }
    }
    int_matrix_print (p);
    int_matrix_print (q);
    int_matrix_print (d);

    // Cleanup allocated memory
    int_matrix_free (d);
    int_matrix_free (p);
    int_matrix_free (q);

    free(p);
    free(q);
    free(d);

    // Output
    return max;
}
