
#include "matrix_generators.h"
#include <cmath>
#include <omp.h>

void generateMatrix(int N, vector<vector<double>> &A, vector<double> &b) {
    A.assign(N, vector<double>(N, 1.0));
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < N; ++i) {
        A[i][i] = 2.0;
    }
    b.assign(N, N + 1.0);
}

void generateMatrixModel(int N, vector<vector<double>> &A, vector<double> &b, vector<double> &u) {
    A.assign(N, vector<double>(N, 1.0));
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < N; ++i) {
        A[i][i] = 2.0;
    }
    u.resize(N);
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < N; ++i) {
        u[i] = sin(2 * M_PI * i / N);
    }
    b.assign(N, 0.0);
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j] * u[j];
        }
        b[i] = sum;
    }
}

double suggestTau(int N, double scale) {
    double lamMax = static_cast<double>(N) + 1.0;
    return scale * 2.0 / lamMax;
}
