#ifndef mat_h
#define mat_h


/* a = b(transp) * a * b */
template<int S>
void mat_square_transform(double (&a)[S][S], double (&b)[S][S]) {
    static double tmp[S][S];

    /* tmp = b(transp) * a */
    for (int i = 0; i < S; ++i) {
        for (int j = 0; j < S; ++j) {
            double sum = 0.0;
            for (int k = 0; k < S; ++k)
                sum += b[k][i] * a[k][j];
            tmp[i][j] = sum;
        }
    }

    /* a = tmp * b */
    for (int i = 0; i < S; ++i) {
        for (int j = 0; j < S; ++j) {
            double sum = 0.0;
            for (int k = 0; k < S; ++k)
                sum += tmp[i][k] * b[k][j];
            a[i][j] = sum;
        }
    }
}

/* A[Q, Q] = B[P, Q](transp) * C[P, P] * D[P, Q] */
template<int P, int Q>
void mat_surround(double A[][Q], double B[][Q], double C[][P], double D[][Q]) {
    static double tmp[P];

    for (int l = 0; l < Q; ++l) {

        /* tmp = B(transp) * C */
        for (int m = 0; m < P; ++m) {
            tmp[m] = 0.0;
            for (int n = 0; n < P; ++n)
                tmp[m] += B[n][l] * C[n][m];
        }

        /* A = tmp * D */
        for (int o = 0; o < Q; ++o) {
            double sum = 0.0;
            for (int n = 0; n < P; ++n)
                sum += tmp[n] * D[n][o];
            A[l][o] += sum;
        }
    }
}

/* C[P, R] = A[P, Q] * B[Q, R] */
template<int P, int Q, int R>
void mat_multiply(double A[][Q], double B[][R], double C[][R]) {
    for (int i = 0; i < P; ++i) {
        for (int j = 0; j < R; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < Q; ++k)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

#endif
