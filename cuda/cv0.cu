#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Include the string.h header for strcspn function
#include <limits.h>

#define TRUE 1
#define FALSE 0

//#undef DEBUG
#define DEBUG

#ifdef DEBUG
#define PRINT printf
#else
#define PRINT // macros
#endif

#define NEGATIVE_INFINITY INT_MIN  // Minimum value for an int
#define MAX_SIZE 10000
#define SCORE_MATCH 1
#define SCORE_MISMATCH (-1)
#define SCORE_OPEN_GAP (-3)
#define SCORE_EXTEND_GAP (-1)


//zero based matrix data structure containing integers
struct Int_Matrix {
    int **matrix;
    int rows;
    int cols;
};



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


//-1 if out of bounds, 0 if in bounds. The result is in ret
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


//use the infinity concept. Actually it should add something to infinity, I never subtract
//or am I wrong?

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


//todo maybe are necessary controls on the length of the input
//todo what happens if _ are present in the input string?

int main() {
    char s1[MAX_SIZE];
    char s2[MAX_SIZE];

    // Initialize strings with some test data (avoid potential buffer overflows)
    strncpy(s1, "ATGGGCTACTTACGCGGTACTA", sizeof(s1) - 1);  // Safe copy with null terminator
    s1[sizeof(s1) - 1] = '\0';           // Ensure null termination

    strncpy(s2, "ATCGGATTGACTCGTGCGGTACGCAGG", sizeof(s2) - 1);  // Safe copy with null terminator
    s2[sizeof(s2) - 1] = '\0';           // Ensure null termination

    // Remove trailing newline characters (if present)
    s1[strcspn(s1, "\n")] = '\0';
    s2[strcspn(s2, "\n")] = '\0';

    PRINT("First string: %s, Second string: %s\n", s1, s2);
    char *y = strlen(s1) > strlen(s2) ? s1 : s2;
    char *x = strlen(s1) < strlen(s2) ? s1 : s2;

    int result = align(y, x, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP);
    printf("Alignment score: %d\n", result);

    return 0;
}

