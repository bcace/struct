#include "sparse.h"
#include <string.h>
#include <assert.h>
#include <stdio.h>


SparseMatrix::SparseMatrix(int _eqs, double *_data) : eqs(_eqs), data(_data) {
    memset(data, 0, sizeof(double) * eqs * eqs);
}

void SparseMatrix::update_element(int i, int j, double v) {
    fprintf(stderr, "%d %d %d\n", i, j, eqs);
    assert(i < eqs && j < eqs);
    data[i * eqs + j] += v;
}
