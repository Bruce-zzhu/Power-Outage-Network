/*
list.c

Implementations for helper functions for linked list construction and
manipulation.

Skeleton written by Grady Fitzpatrick for COMP20007 Assignment 1 2021
*/
#include "list.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



struct list {
  int *item;
  struct list *next;
};

struct list *newlist(int *item){
  struct list *head = (struct list *) malloc(sizeof(struct list));
  assert(head);
  head->item = item;
  head->next = NULL;
  return head;
}

struct list *prependList(struct list *list, int *item){
  struct list *second = (struct list *) malloc(sizeof(struct list));
  assert(second);
  second->item = list->item;
  second->next = list->next;
  list->item = item;
  list->next = second;
  return list;
}

int *peekHead(struct list *list){
  if(! list){
    return NULL;
  }
  return list->item;
}

struct list *deleteHead(struct list *list){
  if(! list){
    return NULL;
  }
  /* Store values we're interested in before freeing list node. */
  if (list->next != NULL) {
      list->item = list->next->item;
      list->next = list->next->next;
      
  }
  else {
      list = NULL;
      free(list);
  }
  
  return list;
}

void freeList(struct list *list){
  struct list *next;
  /* Iterate through list until the end of the list (NULL) is reached. */
  for(next = list; list != NULL; list = next){
    /* Store next pointer before we free list's space. */
    next = list->next;
    free(list);
  }
}


bool inList(struct list *list, int *item) {
  struct list *node = list;
  while (node != NULL) {
    if (*(node->item) == *item) {
      return true;
    }
    node = node->next;
  }
  return false;
}



