#include <iostream>
#include <Eigen/Dense>
#include "EigenGame.hpp"

using namespace Eigen;
using namespace std;

int main() {
    MatrixXf A(3, 3);
    A << 4, 1, 1,
         1, 2, 3,
         1, 3, 10;
    displyMatrix(A);
    return 0;
}
