#include "matrix_generators.h"
#include <cmath>

void generateMatrix(int N, vector<vector<double>> &A, vector<double> &b) {
    A.assign(N, vector<double>(N, 1.0));
    for (int i = 0; i < N; ++i) {
        A[i][i] = 2.0;
    }
    b.assign(N, N + 1.0);
}

double suggestTau(int N, double scale) {
    double lamMax = static_cast<double>(N) + 1.0;
    return scale * 2.0 / lamMax;
}
