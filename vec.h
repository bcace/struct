#ifndef vec_h
#define vec_h


struct dvec3 {
    double x, y, z;

    dvec3();
    dvec3(double _x, double _y, double _z);

    double length();
    dvec3 operator-(const dvec3 &o);
};

#endif
