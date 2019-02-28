#include "beam.h"
#include "node.h"
#include "vec.h"
#include "mat.h"
#include "sparse.h"
#include <math.h>


void beam_add_to_global(SparseMatrix &K, Node &n1, Node &n2, dvec3 up, BeamProperties &props) {
    dvec3 V21 = n2.pos - n1.pos;
    double _2a = V21.length();
    double a = _2a * 0.5;
    double _2a2 = _2a * _2a;
    double _2a3 = _2a2 * _2a;
    double AE_2a = props.AE / _2a;
    double EIy2 = props.EIy2;
    double EIy3 = props.EIy3;
    double EIz2 = props.EIz2;
    double EIz3 = props.EIz3;
    double GJ = props.GJ;

    static double Ke[12][12]; /* keep static because it initializes to 0 */

    /* fill out upper-right triangle */

    Ke[0][0] = AE_2a;
    Ke[0][6] = -AE_2a;

    Ke[1][1] = EIz3 / _2a3;
    Ke[1][5] = EIz3 / _2a2;
    Ke[1][7] = -EIz3 / _2a3;
    Ke[1][11] = EIz3 / _2a2;

    Ke[2][2] = EIy3 / _2a3;
    Ke[2][4] = -EIy3 / _2a2;
    Ke[2][8] = -EIy3 / _2a3;
    Ke[2][10] = -EIy3 / _2a2;

    Ke[3][3] = GJ / _2a;
    Ke[3][9] = -GJ / _2a;

    Ke[4][4] = EIy2 / a;
    Ke[4][8] = EIy3 / _2a2;
    Ke[4][10] = props.E * props.Iy / a;

    Ke[5][5] = EIz2 / a;
    Ke[5][7] = -EIz3 / _2a2;
    Ke[5][11] = props.E * props.Iz / a;

    Ke[6][6] = AE_2a;

    Ke[7][7] = EIz3 / _2a3;
    Ke[7][11] = -EIz3 / _2a2;

    Ke[8][8] = EIy3 / _2a3;
    Ke[8][10] = EIy3 / _2a2;

    Ke[9][9] = GJ / _2a;

    Ke[10][10] = EIy2 / a;

    Ke[11][11] = EIz2 / a;

    /* copy data to the lower-left triangle */

    for (int row = 0; row < 11; ++row)
        for (int col = row + 1; col < 12; ++col)
            Ke[col][row] = Ke[row][col];

    /* generate local to global coordinate system rotation matrix */

    double t[3][3];

    /* normalized vector along beam */
    t[0][0] = V21.x / _2a;
    t[0][1] = V21.y / _2a;
    t[0][2] = V21.z / _2a;

    // TODO: check why z component is up * 0.5, this matrix should be simply [U, V, W]. Also, who ensures that up is normalized?
    /* up vector */
    t[2][0] = up.x * 0.5;
    t[2][1] = up.y * 0.5;
    t[2][2] = up.z * 0.5;

    /* cross product between the first two */
    t[1][0] = t[2][1] * t[0][2] - t[2][2] * t[0][1];
    t[1][1] = t[2][2] * t[0][0] - t[2][0] * t[0][2];
    t[1][2] = t[2][0] * t[0][1] - t[2][1] * t[0][0];

    /* compose rotation matrix for the entire element */

    static double T[12][12];

    for (int i = 0; i < 12; i += 3)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                T[i + j][i + k] = t[j][k];

    /* rotate element stiffness matrix (Ke = T(transp) * Ke * T) */

    mat_square_transform<12>(Ke, T);

    /* insert local into global stiffness matrix */

    Node nodes[] = { n1, n2 };

    for (int col_i = 0; col_i < 2; ++col_i) {
        for (int col_j = 0; col_j < 6; ++col_j) {

            if (nodes[col_i].eqs.v[col_j] == CONSTR_DOF_EQ)
                continue;

            for (int row_i = 0; row_i < 2; ++row_i) {
                for (int row_j = 0; row_j < 6; ++row_j) {

                    if (nodes[row_i].eqs.v[row_j] == CONSTR_DOF_EQ)
                        continue;

                    int col = col_i * 6 + col_j;
                    int row = row_i * 6 + row_j;
                    K.update_element(nodes[row_i].eqs.v[row_j], nodes[col_i].eqs.v[col_j], Ke[row][col]);
                }
            }
        }
    }
}

void BeamProperties::set_rectangular(double _E, double _v, double w, double h) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = w * h;
    Iy = A * h * h * 0.083333333;
    Iz = A * w * w * 0.083333333;
    J = Iy + Iz;
    _update_composites();
}

void BeamProperties::set_circular(double _E, double _v, double r) {
    E = _E;
    v = _v;
    G = E / (2.0 * (1.0 + v));
    A = r * r * M_PI;
    Iy = Iz = M_PI_4 * r * r * r * r;
    J = Iy + Iz;
    _update_composites();
}

void BeamProperties::_update_composites() {
    AE = A * E;
    EIy2 = E * Iy * 2.0;
    EIy3 = E * Iy * 3.0;
    EIz2 = E * Iz * 2.0;
    EIz3 = E * Iz * 3.0;
    GJ = G * J;
}
