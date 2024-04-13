#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Include the string.h header for strcspn function
#include "int_matrix.h"
#include "utility.h"
#include "algo.h"


#define MAX_SIZE 10000
#define SCORE_MATCH 1
#define SCORE_MISMATCH (-1)
#define SCORE_OPEN_GAP (-3)
#define SCORE_EXTEND_GAP (-1)

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

    int result = align(s1, s2, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP);
    printf("Alignment score: %d\n", result);

    return 0;
}
/*
int main() {
    int err;
    int ret;
    struct Int_Matrix *m = (struct Int_Matrix*) malloc(sizeof(struct Int_Matrix));
    err = int_matrix_new(m, 10, 10);
    if (err == -1) {
        return 1;
    }
    for (int i = 0; i < int_matrix_rows_size(m); i ++) {
        for (int j = 0; j < int_matrix_cols_size(m); j ++) {
            int_matrix_set (m, 0, i, j);
        }
    }
    int_matrix_print(m);
    int_matrix_get(m, &ret, 0, -1);
    printf("ret: %d", ret);
    return 0;
}
*/
