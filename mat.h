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

#endif
