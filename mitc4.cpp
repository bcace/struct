#include "mitc4.h"
#include "node.h"
#include "sparse.h"
#include "vec.h"
#include "mat.h"
#include <math.h>
#include <stdio.h>

#define NODE_COUNT 4
#define GAUSS_POINT_COUNT 4
#define SHEAR_CORRECTION_FACTOR 0.833333333


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

void _form_B_shear_matrix(double Bs[][3], dvec2 pos_l[], double J_det, int node, double gp_x, double gp_y) {
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
    Ms[1][0] = 1.0 - gp_x;
    Ms[0][1] = 1.0 - gp_y;
    Ms[1][2] = 1.0 + gp_x;
    Ms[0][3] = 1.0 + gp_y;

    double r1 = Cx + gp_x * Bx;
    double r2 = Ax + gp_y * Bx;
    double r3 = Cy + gp_x * By;
    double r4 = Ay + gp_y * By;
    r1 = sqrt(r1 * r1 + r3 * r3);
    r2 = sqrt(r2 * r2 + r4 * r4);

    static double Bsv[2][12];
    mat_multiply<2, 4, 12>(Ms, G, Bsv);

    double J_det_1 = r1 / (8.0 * J_det);
    double J_det_2 = r2 / (8.0 * J_det);
    for (int i = 0; i < 12; ++i) {
        Bsv[0][i] *= J_det_1;
        Bsv[1][i] *= J_det_2;
    }
    static double Bss[2][12];
    mat_multiply<2, 2, 12>(Rot, Bsv, Bss);

    // TODO: move all the above outside the function since it doesn't use node index.

    // TODO: was: Bs.zero();
    for (int i = 0; i < 3; ++i) {
        Bs[0][i] = Bss[0][node * 3 + i];
        Bs[1][i] = Bss[1][node * 3 + i];
    }
}

void _form_B_drill(double Bd[], double shp_d[][NODE_COUNT], int node, dvec3 &U, dvec3 &V, dvec3 &W) {
    double B1 = -0.5 * shp_d[1][node];
    double B2 =  0.5 * shp_d[0][node];
    double B6 =       -shp_d[2][node];
    Bd[0] = B1 * U.x + B2 * V.x;
    Bd[1] = B1 * U.y + B2 * V.y;
    Bd[2] = B1 * U.z + B2 * V.z;
    Bd[3] = B6 * W.x;
    Bd[4] = B6 * W.y;
    Bd[5] = B6 * W.z;
}

void _form_B_matrix(double B[][6], double Bm[][2], double Bb[][2], double Bs[][3], dvec3 &U, dvec3 &V, dvec3 &W) {

    static double G_membrane_and_bend[2][3];
    G_membrane_and_bend[0][0] = U.x;
    G_membrane_and_bend[0][1] = U.y;
    G_membrane_and_bend[0][2] = U.z;
    G_membrane_and_bend[1][0] = V.x;
    G_membrane_and_bend[1][1] = V.y;
    G_membrane_and_bend[1][2] = V.z;

    static double Bm_shell[3][3];
    mat_multiply<3, 2, 3>(Bm, G_membrane_and_bend, Bm_shell);
    // TODO: was: Bm_shell.mul_of(Bm, G_membrane_and_bend);

    static double Bb_shell[3][3];
    mat_multiply<3, 2, 3>(Bb, G_membrane_and_bend, Bb_shell);
    // TODO: was: Bb_shell.mul_of(Bb, G_membrane_and_bend);

    static double G_shear[3][6];
    // TODO: was: G_shear.zero();
    G_shear[0][0] = W.x;
    G_shear[0][1] = W.y;
    G_shear[0][2] = W.z;
    G_shear[1][3] = U.x;
    G_shear[1][4] = U.y;
    G_shear[1][5] = U.z;
    G_shear[2][3] = V.x;
    G_shear[2][4] = V.y;
    G_shear[2][5] = V.z;

    static double Bs_shell[2][6];
    mat_multiply<2, 3, 6>(Bs, G_shear, Bs_shell);
    // TODO: was: Bs_shell.mul_of(Bs, G_shear);

    // TODO: was: B.zero();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            B[i][j] = Bm_shell[i][j];

    for (int i = 3; i < 6; ++i)
        for (int j = 3; j < 6; ++j)
            B[i][j] = Bb_shell[i - 3][j - 3];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 6; ++j)
            B[i + 6][j] = Bs_shell[i][j];
}

void _insert_element_into_global_stiffness_matrix(SparseMatrix &K, double Ke[][24], Node &n1, Node &n2, Node &n3, Node &n4) {
    Equations eqs[] = { n1.eqs, n2.eqs, n3.eqs, n4.eqs };

    for (int col_i = 0; col_i < 4; ++col_i) {       // nodes
        for (int col_j = 0; col_j < 6; ++col_j) {   // dofs

            if (eqs[col_i].v[col_j] == CONSTR_DOF_EQ)
                continue;

            int col = col_i * 6 + col_j;

            for (int row_i = 0; row_i < 4; ++row_i) {       // nodes
                for (int row_j = 0; row_j < 6; ++row_j) {   // dofs

                    if (eqs[row_i].v[row_j] == CONSTR_DOF_EQ)
                        continue;

                    int row = row_i * 6 + row_j;

                    K.update_element(eqs[row_i].v[row_j], eqs[col_i].v[col_j], Ke[row][col]);
                }
            }
        }
    }
}

void mitc4_add_to_global(SparseMatrix &K, Node &n1, Node &n2, Node &n3, Node &n4, Mitc4Properties &props) {
    static double Ke[24][24];

    /* TODO: this is fishy... */

    double Ktt = props.C[2][2];

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
                    // TODO: was: J[j][k] +=  pos_l[l][k] * gp.shp_d[j][l];
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
            _form_B_shear_matrix(Bs[j], pos_l, J_det, j, gp.x, gp.y);
            _form_B_drill(Bd[j], shp_d, j, U, V, W);
        }

        /* assemble element stiffness matrix */

        for (int j = 0; j < NODE_COUNT; ++j) {
            int j_base = j * 6;

            static double Bj[8][6];
            _form_B_matrix(Bj, Bm[j], Bbn[j], Bs[j], U, V, W);

            for (int k = 0; k < NODE_COUNT; ++k) {
                int k_base = k * 6;

                static double Bk[8][6];
                _form_B_matrix(Bk, Bm[k], Bbp[k], Bs[k], U, V, W);

                static double B[6][6];
                mat_surround(B, Bj, props.C, Bk);

                for (int l = 0; l < 6; ++l)
                    for (int m = 0; m < 6; ++m)
                        Ke[j_base + l][k_base + m] += (B[l][m] + Ktt * Bd[j][l] * Bd[k][m]) * J_det;
            }
        }
    }

    /* insert element stifness matrix into the global one */

    _insert_element_into_global_stiffness_matrix(K, Ke, n1, n2, n3, n4);
}

#include <string.h>

Mitc4Properties::Mitc4Properties(double _E, double _v, double _h) : E(_E), v(_v), h(_h) {
    double M = h * E / (1.0 - v * v);
    double G = h * E / (1.0 + v) * 0.5;
    double D = M * (h * h) / 12.0;
    memset(C, 0, sizeof(double) * 8 * 8); /* TODO: is there a better way to do this? */
    C[0][0] = C[1][1] = M;
    C[2][2] = G;
    C[0][1] = C[1][0] = v * M;
    C[3][3] = C[4][4] = -D;
    C[3][4] = C[4][3] = -v * D;
    C[5][5] = -0.5 * D * (1.0 - v);
    C[6][6] = C[7][7] = G * SHEAR_CORRECTION_FACTOR;
}
