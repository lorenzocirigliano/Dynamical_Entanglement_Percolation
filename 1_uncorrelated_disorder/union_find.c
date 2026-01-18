#include "union_find.h"

UnionFind* create_union_find(int L) {
    UnionFind* uf = malloc(sizeof(UnionFind));
    if (!uf) return NULL;
    
    uf->L = L;
    int size = L * L;
    
    uf->parent = malloc(size * sizeof(int));
    uf->rank = malloc(size * sizeof(int));
    uf->size = malloc(size * sizeof(int));
    
    if (!uf->parent || !uf->rank || !uf->size) {
        destroy_union_find(uf);
        return NULL;
    }
    
    reset_union_find(uf);
    return uf;
}

void destroy_union_find(UnionFind* uf) {
    if (!uf) return;
    
    free(uf->parent);
    free(uf->rank);
    free(uf->size);
    free(uf);
}

void reset_union_find(UnionFind* uf) {
    int size = uf->L * uf->L;
    for (int i = 0; i < size; i++) {
        uf->parent[i] = i;
        uf->rank[i] = 0;
        uf->size[i] = 1;
    }
}

int find(UnionFind* uf, int x) {
    if (uf->parent[x] != x) {
        uf->parent[x] = find(uf, uf->parent[x]); // Path compression
    }
    return uf->parent[x];
}

void union_sets(UnionFind* uf, int x, int y) {
    int root_x = find(uf, x);
    int root_y = find(uf, y);
    
    if (root_x == root_y)
        return;
    
    // Union by rank
    if (uf->rank[root_x] < uf->rank[root_y]) {
        uf->parent[root_x] = root_y;
        uf->size[root_y] += uf->size[root_x];
    } else if (uf->rank[root_x] > uf->rank[root_y]) {
        uf->parent[root_y] = root_x;
        uf->size[root_x] += uf->size[root_y];
    } else {
        uf->parent[root_y] = root_x;
        uf->size[root_x] += uf->size[root_y];
        uf->rank[root_x]++;
    }
}
