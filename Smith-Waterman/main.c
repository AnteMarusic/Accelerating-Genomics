#include <stdio.h>
#include <string.h> // Include the string.h header for strcspn function
#include "algo.h"
#include "utility.h"


#define MAX_SIZE 10000
#define SCORE_MATCH 1
#define SCORE_MISMATCH (-1)
#define SCORE_OPEN_GAP (-2)
#define SCORE_EXTEND_GAP (-2)

//todo maybe are necessary controls on the length of the input
//todo what happens if _ are present in the input string?

int main() {
    char s1[MAX_SIZE];
    char s2[MAX_SIZE];

    printf("____INPUT____\n");
    printf("write two strings to align\n");

    while (fgets(s1, MAX_SIZE, stdin) != NULL && fgets(s2, MAX_SIZE, stdin) != NULL) {
        s1[strcspn(s1, "\n")] = '\0'; // Remove the newline character
        s2[strcspn(s2, "\n")] = '\0'; // Remove the newline character
        PRINT("first string %s, second %s\n", s1, s2);
        struct WholeResult * w_res = align(s1, s2, SCORE_MATCH, SCORE_MISMATCH, SCORE_OPEN_GAP, SCORE_EXTEND_GAP);

        printf("____OUT______\n");
        print_whole_res(w_res);
        printf("____INPUT____\n");
        printf("write two strings to align\n");
    }

    return 0;
}

/*
int main () {
    struct Node * head = NULL;
    add(&head,create_node(0, 1));
    print_list(head);
    PRINT("%d\n", get_length(head));
    PRINT("%d, %d\n", get(head, 0)->row, get(head, 0)->col);
    delete_all(&head);
    print_list(head);
}
 */