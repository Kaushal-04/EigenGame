#include<iostream>
#include<iterator>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<random>
#include<cmath>
#include<math.h>
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()
using namespace std;
using namespace Eigen;
#define MAX_SIZE 10000
int n;
void setValueOfn(){
    cout<<"Enter order of matrix:";
    cin>>n;
    if (n > MAX_SIZE) {
        cout << "Error: Order exceeds maximum size. Setting order to maximum size.\n";
        n = MAX_SIZE;
    }
}
MatrixXf makeSymmetric(const MatrixXf& A) {
    return 0.5 * (A + A.transpose());
}
MatrixXf normalizeVector(const MatrixXf& vect) {
    MatrixXf normalizedVector = vect;
    float maxElement = vect.maxCoeff();
    for (int i = 0; i < vect.size(); ++i) {
        normalizedVector(i) /= maxElement;
    }
    return normalizedVector;
}
void insertColumnMatrix(MatrixXf& originalMatrix, const MatrixXf& columnMatrix) {
    if (columnMatrix.rows() != originalMatrix.rows()) {
        cerr << "Dimention of Column & Rows not match" <<endl;
        return;
    }
    originalMatrix.conservativeResize(originalMatrix.rows(), originalMatrix.cols() + columnMatrix.cols());
    originalMatrix.rightCols(columnMatrix.cols()) = columnMatrix;
}
void deleteColumn(MatrixXf& matrix, int columnIndex) {
    if (columnIndex < 0 || columnIndex >= matrix.cols()) {
        cerr << "Error: Column index out of bounds." << endl;
        return;
    }
    for (int i = columnIndex; i < matrix.cols() - 1; ++i) {
        matrix.col(i) = matrix.col(i + 1);
    }
    matrix.conservativeResize(matrix.rows(), matrix.cols() - 1);
}
MatrixXf eigenValues(const MatrixXf& eigenvectors) {
    SelfAdjointEigenSolver<MatrixXf> solver(eigenvectors.transpose() * eigenvectors);
    return solver.eigenvalues();
}
int main(){
    setValueOfn();
    srand(time(0));
    MatrixXf FinalEigenVector(n,1);
    MatrixXf randomMatrixA , randomMatrixB;
    randomMatrixA.setRandom(n,n);
    randomMatrixA = randomMatrixA.array().abs(); // Ensure positive values
    randomMatrixB.setRandom(n,n);
    randomMatrixB = randomMatrixB.array().abs(); // Ensure positive values
    //cout<<"MAtrix A:\n"<<randomMatrixA<<endl;
    //cout<<"Matrix B:\n"<<randomMatrixB<<endl;
    MatrixXf symmetricMatrixA = makeSymmetric(randomMatrixA);
    MatrixXf symmetricMatrixB = makeSymmetric(randomMatrixB);
    cout << "Symmetric Matrix A:\n" << symmetricMatrixA << endl;
    cout << "Symmetric Matrix B:\n" << symmetricMatrixB << endl;

    //Create <float> , dynamic matrix   & Finding EigenValues and Eigen vectors
    MatrixXcf eigenValMatA , eigenValMatB;
    MatrixXcf eigenVectMatA , eigenVectMatB;
    EigenSolver<MatrixXf> eigenValueSolverA(symmetricMatrixA);
    EigenSolver<MatrixXf> eigenValueSolverB(symmetricMatrixB);
    eigenValMatA=eigenValueSolverA.eigenvalues();
    eigenValMatB=eigenValueSolverB.eigenvalues();
    eigenVectMatA=eigenValueSolverA.eigenvectors();
    eigenVectMatB=eigenValueSolverB.eigenvectors();

    MatrixXf EigenVector(n,1);
    insertColumnMatrix(EigenVector , eigenVectMatA.real());
    insertColumnMatrix(EigenVector , eigenVectMatB.real());
    deleteColumn(EigenVector , 0);
    for (int column = 0; column < EigenVector.cols(); ++column) {
        EigenVector.col(column) = normalizeVector(EigenVector.col(column));
    }

    int T=3;  //just for checking otherwise set vale to 10001
    MatrixXf Reward;
    MatrixXf operation(n, n);
    MatrixXf result;
    MatrixXf rewCala , rewCalaR , rewResl, rewCalb , rewCalbR , rewResR;
    float rewa , rewb;
    int colNum;
    MatrixXf vi , vj;
 for (int j = 0; j < T; j++) { //T
     for (int i = 0; i < n; i++) {  //n or k

        colNum = rand()%EigenVector.cols();
        vi = EigenVector.col(colNum);
        colNum = rand()%EigenVector.cols();
        vj = EigenVector.col(colNum);

        float result = (vj.transpose() * symmetricMatrixB * vj)(0, 0);
        if(result < 0) //to fin nan result
            result = result * (-1);
        float sqrtResult = sqrt(result);
        sqrtResult = 1.0 /sqrtResult;
        MatrixXf yj;
        yj = vj / sqrtResult;
        rewCala= vi.transpose() * symmetricMatrixB * vi;
        rewa=rewCala(0,0);
        rewCalaR=symmetricMatrixA * vi;
        rewCalb=vi.transpose() * symmetricMatrixA * vi;
        rewb=rewCalb(0,0);
        rewCalbR=symmetricMatrixB * vi;
        rewResl=rewCalaR * rewa;
        rewResR=rewCalbR * rewb;
        Reward=rewResl - rewResR;
        if(j < i){
            MatrixXf PenARes , PenBRes , diffRes;
            MatrixXf Penalties(n,1) ;
            for(int row=0;  row<n; row++){
                Penalties(row , 0) = 0;
            }
            float PenAScal, PenBScal , PenCScal;
            MatrixXf PenA , PenB , PenC , PenVecA , PenVecB;
            PenA = vi.transpose() * symmetricMatrixA * yj;
            PenAScal=PenA(0,0);
            PenB = vi.transpose() * symmetricMatrixB * vi ;
            PenBScal = PenB(0,0);
            PenC = vi.transpose() * symmetricMatrixB * yj ;
            PenCScal = PenC(0,0);
            PenVecA = symmetricMatrixB * yj;
            PenVecB = symmetricMatrixB * vi;
            PenARes = PenBScal * PenVecA;
            PenBRes = PenCScal * PenVecB;
            diffRes = PenARes - PenBRes;
            Penalties = Penalties + (PenAScal * diffRes);
            MatrixXf delta(n,1);
            delta = Reward - Penalties;
            MatrixXf itadelta(n,1);
            itadelta = -2 * delta;
            MatrixXf wi(n,1);
            wi = vi + itadelta;
            vi = normalizeVector(wi);
            insertColumnMatrix(FinalEigenVector , vi);
         }
     }
 }
    deleteColumn(FinalEigenVector , 0);
    cout<<"Final Eigen Vector:\n"<<FinalEigenVector<<endl;
    MatrixXf EigenValues = eigenValues(FinalEigenVector);
    cout<<"Eigen Values Matrix : \n"<<EigenValues<<endl;
    return 0;
}
