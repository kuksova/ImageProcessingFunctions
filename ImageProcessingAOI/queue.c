#include "stdafx.h"
#include <limits.h>
#include "queue.h"

void initQueue(Queue *q)
{
    q->front = NULL;
    q->back = NULL;
    q->length = 0;
}

bool isEmpty(Queue *q)
{
    return q->length == 0;
}

bool isFull(Queue *q)
{
    return false;
}

int enqueue(Queue *q, void * el)
{
    if(!isFull(q)){
        QueueElement * to_add = malloc(sizeof(QueueElement));
        if(to_add == NULL){
            return -1;
        }
    
        to_add -> value = el;
        to_add -> next = NULL;

        if(q->length == 0){
            q->front = to_add;
            q->back = to_add;
        }
        else{
            q->back->next = to_add;
            q->back = to_add;
        }
        (q->length)++;

        return 1;
    }

    return -1;
}

void * dequeue(Queue *q)
{
    void * result;

    if(!isEmpty(q)){
        result= q->front->value;

        QueueElement * removed = q->front;
        q->front = q->front->next;
        free(removed);
        (q->length)--;

        if(q->length == 0){
            q->back = q->front = NULL;
        }
    }

    else{
        result = NULL;
    }

    return result;
}
