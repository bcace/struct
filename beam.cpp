#include "beam.h"
#include "node.h"
#include "vec.h"
#include <math.h>


void beam_add_to_global(double *K, Node &n1, Node &n2, dvec3 up, BeamProperties &props) {
    // double (*K)[2] = (double(*)[2])_K;
    // K_spec[3][1] = 4.5;

    double _2a = (n2.pos - n1.pos).length();
    double a = _2a * 0.5;
    double _2a2 = _2a * _2a;
    double _2a3 = _2a2 * _2a;
    double AE_2a = props.AE / _2a;
    double EIy2 = props.EIy2;
    double EIy3 = props.EIy3;
    double EIz2 = props.EIz2;
    double EIz3 = props.EIz3;
    double GJ = props.GJ;

    static double Ke[12][12];
}

void BeamProperties::set_rectangular(double _E, double _v, double w, double h) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = w * h;
    Iy = A * h * h * 0.083333333;
    Iz = A * w * w * 0.083333333;
    J = Iy + Iz;

    /* pre-calculate some frequently used composite values */
    AE = A * E;
    EIy2 = E * Iy * 2.0;
    EIy3 = E * Iy * 3.0;
    EIz2 = E * Iz * 2.0;
    EIz3 = E * Iz * 3.0;
    GJ = G * J;
}

void BeamProperties::set_circular(double _E, double _v, double r) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = r * r * M_PI;
    Iy = Iz = M_PI_4 * r * r * r * r;
    J = Iy + Iz;

    /* pre-calculate some frequently used composite values */
    AE = A * E;
    EIy2 = E * Iy * 2.0;
    EIy3 = E * Iy * 3.0;
    EIz2 = E * Iz * 2.0;
    EIz3 = E * Iz * 3.0;
    GJ = G * J;
}

