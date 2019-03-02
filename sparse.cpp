#include "sparse.h"
#include <string.h>
#include <assert.h>
#include <stdio.h>


SparseMatrix::SparseMatrix(int _eqs, double *_data) : eqs(_eqs), data(_data) {
    memset(data, 0, sizeof(double) * eqs * eqs);
}

void SparseMatrix::update_element(int i, int j, double v) {
    assert(i < eqs && j < eqs);
    data[i * eqs + j] += v;
}

/* TODO: clear up terminology */
void SparseMatrix::solve(double *interface) {

    /* transforming into a triangular matrix */
    for (int i = 0; i < eqs - 1; ++i) {
        for (int j = i + 1; j < eqs; ++j) {
            double elimination_factor = -data[j * eqs + i] / data[i * eqs + i];
            data[j * eqs + i] = 0.0;

            for (int k = i + 1; k < eqs; ++k)
                data[j * eqs + k] += data[i * eqs + k] * elimination_factor;

            interface[j] += interface[i] * elimination_factor;
        }
    }

    /* gauss solve */
    for (int i = eqs - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = eqs - 1; j > i; --j)
            sum += data[i * eqs + j] * interface[j];
        interface[i] = (interface[i] - sum) / data[i * eqs + i];
    }
}

void SparseMatrix::print() {
    for (int i = 0; i < eqs; ++i) {
        for (int j = 0; j < eqs; ++j)
            printf("%g ", data[i * eqs + j]);
        printf("\n");
    }
}
