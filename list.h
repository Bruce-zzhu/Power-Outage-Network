/*
list.h

Visible structs and functions for linked lists.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/

#include <stdbool.h>

/* The linked list. */
struct list;

/* Get a new empty list. */
struct list *newlist(int *item);

/* Add an item to the head of the list. Returns the new list. */
struct list *prependList(struct list *list, int *item);

/* Gets the first item from the list. */
int *peekHead(struct list *list);

/* Takes the first item from the list, updating the list pointer and returns
  the item stored. */
struct list *deleteHead(struct list *list);

/* Free all list items. */
void freeList(struct list *list);

/* chech whether the item is in the list */
bool inList(struct list *list, int *item);

