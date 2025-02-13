#include "matrix_generators.h"
#include <cmath>

// Генерирует матрицу A и вектор b для режима 2
// Здесь A имеет вид: A[i][j] = 2 если i==j, иначе 1,
// а b[i] = N+1 для всех i.
void generateMatrix(int N, vector<vector<double>> &A, vector<double> &b) {
    A.assign(N, vector<double>(N, 1.0));
    for (int i = 0; i < N; ++i) {
        A[i][i] = 2.0;
    }
    b.assign(N, N + 1.0);
}

double suggestTau(int N, double scale) {
    if (N <= 0) return 0.01; // на всякий случай
    double lamMax = static_cast<double>(N) + 1.0;
    return scale * 2.0 / lamMax;
}
