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

    void _calculate_composites();
    void set_rectangular(double _E, double _v, double w, double h);
    void set_circular(double _E, double _v, double r);
};

#endif