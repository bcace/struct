#include "vec.h"
#include <math.h>


dvec3::dvec3() {}

dvec3::dvec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

double dvec3::length() {
    return sqrt(x * x + y * y + z * z);
}

dvec3 dvec3::operator-(const dvec3 &o) {
    return dvec3(x - o.x, y - o.y, z - o.z);
}
