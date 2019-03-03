#ifndef vec_h
#define vec_h


struct dvec3 {
    union {
        struct {
            double x, y, z;
        };
        double v[3];
    };

    dvec3() {}
    dvec3(double _x, double _y, double _z);

    double length();
    dvec3 operator-(dvec3 v);
    dvec3 operator+(dvec3 v);
    dvec3 operator*(double v);
    dvec3 &operator-=(dvec3 v);
    double dot(dvec3 v);
    dvec3 cross(dvec3 v);
    dvec3 &normalize();
    double &operator[](int i);
    void print();
};


struct dvec2 {
    union {
        struct {
            double x, y;
        };
        double v[2];
    };

    dvec2();
    dvec2(double _x, double _y);
    double &operator[](int i);
    void print();
};

#endif
