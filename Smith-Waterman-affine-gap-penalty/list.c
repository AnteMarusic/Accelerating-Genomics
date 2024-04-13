#include "list.h"
#include <stdio.h>
#include <stdlib.h>
#include "utility.h"

void list_init (struct Node ** head) {
    *head = NULL;
}

struct Node* create_node(int row, int col) {
    struct Node* new_node = (struct Node*)malloc(sizeof(struct Node));
    new_node->row = row;
    new_node->col = col;
    new_node->next = NULL;
    return new_node;
}

void add(struct Node** head, struct Node* node) {
    node->next = *head;
    *head = node;
}

// Delete all elements from the list and free memory
void delete_all(struct Node** head) {
    struct Node* current = *head;
    struct Node* temp;
    while (current != NULL) {
        temp = current;
        current = current->next;
        free(temp);
    }
    *head = NULL;
}

void print_list(struct Node* head) {
    struct Node* current = head;
    PRINT("List: ");
    while (current != NULL) {
        PRINT("(%d, %d) ", current->row, current->col);
        current = current->next;
    }
    PRINT("\n");
}

int get_length(struct Node* head) {
    int count = 0;
    struct Node* current = head;
    while (current != NULL) {
        count++;
        current = current->next;
    }
    return count;
}

struct Node* get(struct Node* head, int index) {
    if (head == NULL || index < 0) {
        // Handle invalid cases (empty list or negative index)
        printf("Error: Invalid index\n");
        return NULL;
    }

    struct Node *current = head;
    int count = 0;
    while (current != NULL && count < index) {
        current = current->next;
        count++;
    }

    if (current == NULL) {
        // Handle case where index is out of bounds (exceeds list length)
        printf("Error: Index out of bounds\n");
        return NULL;
    }
    return current;
}