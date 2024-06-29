#include <iostream>
#include <Eigen/Dense>
#include "EigenGame.hpp"

using namespace Eigen;
using namespace std;
MatrixXf makeSymmetric(const MatrixXf& A) {
    return 0.5 * (A + A.transpose());
}
MatrixXf eigenValues(const MatrixXf& eigenvectors) {
    SelfAdjointEigenSolver<MatrixXf> solver(eigenvectors.transpose() * eigenvectors);
    return solver.eigenvalues();
}
int main() {
    int n = 4;
    MatrixXf matrixA , matrixB;
    matrixA.setRandom(n,n);
    matrixA = matrixA.array().abs(); // Ensure positive values
    matrixB.setRandom(n,n);
    matrixB = matrixB.array().abs(); // Ensure positive values
    MatrixXf A = makeSymmetric(matrixA);
    MatrixXf B = makeSymmetric(matrixB);
    MatrixXf EigenVector = solveEigenGame(A, B, n);
    MatrixXf EigenValues = eigenValues(EigenVector);
    cout<<"Eigen Vector\n"<<EigenVector<<endl;
    return 0;
}
