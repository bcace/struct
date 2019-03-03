#include "mitc4.h"
#include "node.h"
#include "vec.h"
#include <math.h>
#include <stdio.h>

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


void _form_B_membrane_matrix(double Bm[][2], double shp_d[][NODE_COUNT], int node) {
    Bm[0][0] = shp_d[0][node];
    Bm[0][1] = 0.0;
    Bm[1][0] = 0.0;
    Bm[1][1] = shp_d[1][node];
    Bm[2][0] = shp_d[1][node];
    Bm[2][1] = shp_d[0][node];
}

void _form_B_positive_bending_matrix(double Bb[][2], double shp_d[][NODE_COUNT], int node) {
    Bb[0][0] = 0.0;
    Bb[0][1] = -shp_d[0][node];
    Bb[1][0] =  shp_d[1][node];
    Bb[1][1] = 0.0;
    Bb[2][0] =  shp_d[0][node];
    Bb[2][1] = -shp_d[1][node];
}

void _form_B_negative_bending_matrix(double Bb[][2], double shp_d[][NODE_COUNT], int node) {
    Bb[0][0] = 0.0;
    Bb[0][1] =  shp_d[0][node];
    Bb[1][0] = -shp_d[1][node];
    Bb[1][1] = 0.0;
    Bb[2][0] = -shp_d[0][node];
    Bb[2][1] =  shp_d[1][node];
}

void _form_B_shear_matrix(double Bs[][3], dvec2 pos_l[], double J_det, int node, dvec2 gps) {
    double dx34 = pos_l[2].x - pos_l[3].x;
    double dy34 = pos_l[2].y - pos_l[3].y;
    double dx21 = pos_l[1].x - pos_l[0].x;
    double dy21 = pos_l[1].y - pos_l[0].y;
    double dx32 = pos_l[2].x - pos_l[1].x;
    double dy32 = pos_l[2].y - pos_l[1].y;
    double dx41 = pos_l[3].x - pos_l[0].x;
    double dy41 = pos_l[3].y - pos_l[0].y;

    static double G[4][12];
    G[0][0]     = -0.5;
    G[0][1]     = -0.25 * dy41;
    G[0][2]     =  0.25 * dx41;
    G[0][9]     =  0.5;
    G[0][10]    = -0.25 * dy41;
    G[0][11]    =  0.25 * dx41;
    G[1][0]     = -0.5;
    G[1][1]     = -0.25 * dy21;
    G[1][2]     =  0.25 * dx21;
    G[1][3]     =  0.5;
    G[1][4]     = -0.25 * dy21;
    G[1][5]     =  0.25 * dx21;
    G[2][3]     = -0.5;
    G[2][4]     = -0.25 * dy32;
    G[2][5]     =  0.25 * dx32;
    G[2][6]     =  0.5;
    G[2][7]     = -0.25 * dy32;
    G[2][8]     =  0.25 * dx32;
    G[3][6]     =  0.5;
    G[3][7]     = -0.25 * dy34;
    G[3][8]     =  0.25 * dx34;
    G[3][9]     = -0.5;
    G[3][10]    = -0.25 * dy34;
    G[3][11]    =  0.25 * dx34;

    double Ax = -pos_l[0].x + pos_l[1].x + pos_l[2].x - pos_l[3].x;
    double Bx =  pos_l[0].x - pos_l[1].x + pos_l[2].x - pos_l[3].x;
    double Cx = -pos_l[0].x - pos_l[1].x + pos_l[2].x + pos_l[3].x;

    double Ay = -pos_l[0].y + pos_l[1].y + pos_l[2].y - pos_l[3].y;
    double By =  pos_l[0].y - pos_l[1].y + pos_l[2].y - pos_l[3].y;
    double Cy = -pos_l[0].y - pos_l[1].y + pos_l[2].y + pos_l[3].y;

    double alph = atan(Ay / Ax);
    double beta = M_PI_2 - atan(Cx / Cy);

    static double Rot[2][2];
    Rot[0][0] =  sin(beta);
    Rot[0][1] = -sin(alph);
    Rot[1][0] = -cos(beta);
    Rot[1][1] =  cos(alph);

    static double Ms[2][4];
    Ms[1][0] = 1.0 - gps.x;
    Ms[0][1] = 1.0 - gps.y;
    Ms[1][2] = 1.0 + gps.x;
    Ms[0][3] = 1.0 + gps.y;

    double r1 = Cx + gps.x * Bx;
    double r2 = Ax + gps.y * Bx;
    double r3 = Cy + gps.x * By;
    double r4 = Ay + gps.y * By;
    r1 = sqrt(r1 * r1 + r3 * r3);
    r2 = sqrt(r2 * r2 + r4 * r4);

    static double Bsv[2][12];
    Bsv.mul_of(Ms, G);

    double J_det_1 = r1 / (8.0 * J_det);
    double J_det_2 = r2 / (8.0 * J_det);
    for (int i = 0; i < 12; ++i) {
        Bsv[0][i] *= J_det_1;
        Bsv[1][i] *= J_det_2;
    }
    static double Bss[2][12];
    Bss.mul_of(Rot, Bsv);

    // TODO: move all the above outside the function since it doesn't use node index.

    // Bs.zero();
    for (int i = 0; i < 3; ++i) {
        Bs[0][i] = Bss[0][node * 3 + i];
        Bs[1][i] = Bss[1][node * 3 + i];
    }
}

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

        /* strain matrices */

        static double Bm[NODE_COUNT][3][2];
        static double Bbp[NODE_COUNT][3][2];
        static double Bbn[NODE_COUNT][3][2];
        static double Bs[NODE_COUNT][2][3];
        static double Bd[NODE_COUNT][6];

        for (int j = 0; j < NODE_COUNT; ++j) {
            _form_B_membrane_matrix(Bm[j], shp_d, j);
            _form_B_positive_bending_matrix(Bbp[j], shp_d, j);
            _form_B_negative_bending_matrix(Bbn[j], shp_d, j);
            _form_B_shear_matrix(Bs[j], pos_l, J_det, j, gp.pos);
            // _form_B_drill(Bd[j], shp_d, j, U, V, W);
        }
    }
}
