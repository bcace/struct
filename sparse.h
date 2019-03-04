#ifndef sparse_h
#define sparse_h


struct SparseMatrix {
    int eqs;
    double *data;

    SparseMatrix(int _eqs, double *_data);

    void update_element(int i, int j, double v);
    void solve(double *interface);
    void print();

    // template<int NODES, int DOFS>
    // void insert_element(double Ke[][NODES * DOFS], unsigned node_eqs[][DOFS]) {

    //     for (int col_i = 0; col_i < NODES; ++col_i) {
    //         for (int col_j = 0; col_j < DOFS; ++col_j) {

    //             if (node_eqs[col_i][col_j] > DOFS)
    //                 continue;

    //             int col = col_i * DOFS + col_j;

    //             for (int row_i = 0; row_i < NODES; ++row_i) {
    //                 for (int row_j = 0; row_j < DOFS; ++row_j) {

    //                     if (node_eqs[row_i][row_j] > DOFS)
    //                         continue;

    //                     int row = row_i * DOFS + row_j;

    //                     data[node_eqs[row_i][row_j] * eqs + node_eqs[col_i][col_j]] += Ke[row][col];
    //                 }
    //             }
    //         }
    //     }
    // }
};

#endif
