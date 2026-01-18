#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "config.h"

typedef struct {
    int *parent;
    int *rank;
    int *size;
    int L;
} UnionFind;

UnionFind* create_union_find(int L);
void destroy_union_find(UnionFind* uf);
void reset_union_find(UnionFind* uf);
int find(UnionFind* uf, int x);
void union_sets(UnionFind* uf, int x, int y);

#endif // UNION_FIND_H
