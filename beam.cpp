#include "beam.h"
#include <math.h>


void BeamProperties::_calculate_composites() {
    AE = A * E;
    EIy2 = E * Iy * 2.0;
    EIy3 = E * Iy * 3.0;
    EIz2 = E * Iz * 2.0;
    EIz3 = E * Iz * 3.0;
    GJ = G * J;
}

void BeamProperties::set_rectangular(double _E, double _v, double w, double h) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = w * h;
    Iy = A * h * h * 0.083333333;
    Iz = A * w * w * 0.083333333;
    J = Iy + Iz;

    _calculate_composites();
}

void BeamProperties::set_circular(double _E, double _v, double r) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = r * r * M_PI;
    Iy = Iz = M_PI_4 * r * r * r * r;
    J = Iy + Iz;

    _calculate_composites();
}

