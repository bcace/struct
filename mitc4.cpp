#include "mitc4.h"
#include "node.h"
#include "vec.h"

#define NODE_COUNT 4
#define GAUSS_POINT_COUNT 4


struct Gp {
    double x, y; /* natural CS [-1, 1] */
    double shp_d[3][NODE_COUNT]; /* shape function derivatives at Gauss point */

    Gp(double _x, double _y) : x(_x), y(_y) {
        /* dN / d_ksi */
        shp_d[0][0] = 0.25 * (-1.0 + y);
        shp_d[0][1] = 0.25 * ( 1.0 - y);
        shp_d[0][2] = 0.25 * ( 1.0 + y);
        shp_d[0][3] = 0.25 * (-1.0 - y);

        /* dN / d_eta */
        shp_d[1][0] = 0.25 * (-1.0 + x);
        shp_d[1][1] = 0.25 * (-1.0 - x);
        shp_d[1][2] = 0.25 * ( 1.0 + x);
        shp_d[1][3] = 0.25 * ( 1.0 - x);

        /* dN / d_ksi_eta */
        shp_d[2][0] = 0.5 * (1.0 - x) * (1.0 - y);
        shp_d[2][1] = 0.5 * (1.0 + x) * (1.0 - y);
        shp_d[2][2] = 0.5 * (1.0 + x) * (1.0 + y);
        shp_d[2][3] = 0.5 * (1.0 - x) * (1.0 + y);
    }
} gps[GAUSS_POINT_COUNT] = {
    Gp(-0.577350269189626, -0.577350269189626),
    Gp( 0.577350269189626, -0.577350269189626),
    Gp( 0.577350269189626,  0.577350269189626),
    Gp(-0.577350269189626,  0.577350269189626)
};


void mitc4_add_to_global(SparseMatrix &K, Node &n1, Node &n2, Node &n3, Node &n4) {
    static double Ke[24][24];

    /* basis vectors

        n4------n3
        |       |
        |       |
        n1------n2
    */

    dvec3 U = (n3.pos + n2.pos - n4.pos - n1.pos) * 0.5;
    U.normalize();
    dvec3 V = (n4.pos + n3.pos - n2.pos - n1.pos) * 0.5;
    V -= U * U.dot(V);
    V.normalize();
    dvec3 W = U.cross(V);

    /* node positions in local CS */

    dvec3 pos[] = { n1.pos, n2.pos, n3.pos, n4.pos };

    static dvec2 pos_l[NODE_COUNT];
    for (int i = 0; i < NODE_COUNT; ++i) {
        pos_l[i].x = U.dot(pos[i]);
        pos_l[i].y = V.dot(pos[i]);
    }

    for (int i = 0; i < GAUSS_POINT_COUNT; ++i) {
        Gp &gp = gps[i];

        /* compute jacobian */

        static double J[2][2];
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                J[j][k] = 0.0;
                for (int l = 0; l < NODE_COUNT; ++l)
                    J[j][k] +=  pos_l[l][j] * gp.shp_d[k][l];
                    // J[j][k] +=  pos_l[l][k] * gp.shp_d[j][l];
            }
        }

        /* jacobian determinant */

        double J_det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        /* inverted jacobian */

        static double J_inv[2][2];
        J_inv[0][0] =  J[1][1] / J_det;
        J_inv[0][1] = -J[0][1] / J_det;
        J_inv[1][0] = -J[1][0] / J_det;
        J_inv[1][1] =  J[0][0] / J_det;

        /* shape function derivatives */

        static double shp_d[3][NODE_COUNT];

        for (int j = 0; j < NODE_COUNT; ++j) {
            shp_d[0][j] = gp.shp_d[0][j] * J_inv[0][0] + gp.shp_d[1][j] * J_inv[1][0];
            shp_d[1][j] = gp.shp_d[0][j] * J_inv[0][1] + gp.shp_d[1][j] * J_inv[1][1];
            shp_d[2][j] = gp.shp_d[2][j];
        }


    }
}
