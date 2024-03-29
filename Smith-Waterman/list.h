#ifndef SMITH_WATERMAN_CLION_LIST_H
#define SMITH_WATERMAN_CLION_LIST_H

struct Node {
    int row;
    int col;
    struct Node* next;
};

void list_init (struct Node ** head);
struct Node* create_node(int row, int col);
void add(struct Node** head, struct Node* node);
void delete_all(struct Node** head);
void print_list(struct Node* head);
int get_length(struct Node* head);
struct Node* get(struct Node* head, int index);

#endif
