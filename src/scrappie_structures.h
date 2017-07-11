#pragma once
#ifndef SCRAPPIE_STRUCTURES_H
#    define SCRAPPIE_STRUCTURES_H

#    include <stddef.h>

typedef struct {
    double start;
    float length;
    float mean;
    float stdv;
    int pos;
    int state;
} event_t;

typedef struct {
    size_t n;
    size_t start;
    size_t end;
    event_t *event;
} event_table;

typedef struct {
    size_t n;
    size_t start;
    size_t end;
    float *raw;
} raw_table;

#endif                          /* SCRAPPIE_DATA_H */
