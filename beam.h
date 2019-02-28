#ifndef beam_h
#define beam_h


struct BeamProperties {
    double E;       // Young's modulus
    double v;       // Poisson's ratio

    // derived
    double A;       // area
    double G;       // shear modulus
    double Iy, Iz;  // second moments of area (inertia)
    double J;       // polar moment of inertia

    // composites
    double AE;
    double EIy2, EIy3, EIz2, EIz3;
    double GJ;

    void set_rectangular(double _E, double _v, double w, double h);
    void set_circular(double _E, double _v, double r);
    void _update_composites();
};


struct Node;
struct dvec3;
struct SparseMatrix;

void beam_add_to_global(SparseMatrix &K, Node &n1, Node &n2, dvec3 up, BeamProperties &props);

#endif
