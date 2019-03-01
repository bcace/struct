#ifndef sparse_h
#define sparse_h


struct SparseMatrix {
    int eqs;
    double *data;

    SparseMatrix(int _eqs, double *_data);

    void update_element(int i, int j, double v);
    void solve(double *interface);
};

#endif
