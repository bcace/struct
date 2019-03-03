#include "vec.h"
#include <math.h>
#include <string.h>
#include <stdio.h>


dvec3::dvec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

double dvec3::length() {
    return sqrt(x * x + y * y + z * z);
}

dvec3 dvec3::operator-(dvec3 v) {
    return dvec3(x - v.x, y - v.y, z - v.z);
}

dvec3 dvec3::operator+(dvec3 v) {
    return dvec3(x + v.x, y + v.y, z + v.z);
}

dvec3 dvec3::operator*(double v) {
    return dvec3(x * v, y * v, z * v);
}

dvec3 &dvec3::operator-=(dvec3 v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

double dvec3::dot(dvec3 v) {
    return x * v.x + y * v.y + z * v.z;
}

dvec3 dvec3::cross(dvec3 v) {
    return dvec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

dvec3 &dvec3::normalize() {
    double l = length();
    x /= l;
    y /= l;
    z /= l;
    return *this;
}

double &dvec3::operator[](int i) {
    return v[i];
}

void dvec3::print() {
    fprintf(stderr, "%g, %g, %g\n", x, y, z);
}

/* :dvec2 */

dvec2::dvec2() {}

dvec2::dvec2(double _x, double _y) : x(_x), y(_y) {}

double &dvec2::operator[](int i) {
    return v[i];
}

void dvec2::print() {
    fprintf(stderr, "%g, %g\n", x, y);
}
