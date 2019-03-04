#ifndef mitc4_h
#define mitc4_h


struct Mitc4Properties {
    double E;
    double v;
    double h;
    double C[8][8];

    Mitc4Properties(double _E, double _v, double _h);
};

struct Node;
struct SparseMatrix;

void mitc4_add_to_global(SparseMatrix &K, Node &n1, Node &n2, Node &n3, Node &n4, Mitc4Properties &props);

#endif
