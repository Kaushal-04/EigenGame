#include <iostream>
#include <Eigen/Dense>
#include "EigenGame.hpp"

using namespace Eigen;
using namespace std;

int main() {
    int n = 5;
    MatrixXf matrixA , matrixB;
    matrixA.setRandom(n,n);
    matrixA = matrixA.array().abs(); // Ensure positive values
    matrixB.setRandom(n,n);
    matrixB = matrixB.array().abs(); // Ensure positive values
    solveEigenGame(matrixA, matrixB, n);
    return 0;
}
