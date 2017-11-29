#ifndef QUEUE_H
#define QUEUE_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct QueueElement{
    void * value;
    struct QueueElement * next;
}QueueElement;

typedef struct Queue{
    QueueElement * front;
    QueueElement * back;

    size_t length;
}Queue;

void initQueue(Queue *q);
int enqueue(Queue *q, void * el);
void * dequeue(Queue *q);
bool isEmpty(Queue *q);
bool isFull(Queue *q);
#endif
