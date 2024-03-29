/*  This function performs the aligment of two strings using smith waterman algorithm
    match_score is the only positive int
*/


#ifndef ALGO_H
#define ALGO_H
struct Result {
    char *s1;
    char *s2;
};

struct WholeResult {
    struct Result * res;
    int len;
};

//todo implement the gap open score
struct WholeResult * align(char* s1, char* s2, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score);
void print_whole_res (struct WholeResult *w_res);
#endif
