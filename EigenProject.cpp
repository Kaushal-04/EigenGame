#include<iostream>
#include<iterator>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<random>
#include<cmath>
#include<math.h>
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
MatrixXd makeSymmetric(const MatrixXd& A) {
    return 0.5 * (A + A.transpose());
}
MatrixXcd getRandomEigenVector(const MatrixXcd& eigenVectMat) {
    // Get the number of eigen vectors
    int numEigenVec = eigenVectMat.cols();
    random_device rd; 
    mt19937 gen(rd()); 
    uniform_int_distribution<> distrib(0, numEigenVec - 1); // Define the range
    int randomIndex = distrib(gen);
    return eigenVectMat.col(randomIndex);
}
MatrixXd normalizeEigenVector(const MatrixXd& eigenVector) {
    double magnitude = eigenVector.norm(); // Compute the magnitude
    return eigenVector / magnitude; // Normalize by dividing each component by the magnitude
}
void insertColumnMatrix(MatrixXd& originalMatrix, const MatrixXd& columnMatrix) {
    if (columnMatrix.rows() != originalMatrix.rows()) {
        cerr << "Error: Dimensions of the column matrix do not match the number of rows in the original matrix." << std::endl;
        return;
    }
    originalMatrix.conservativeResize(originalMatrix.rows(), originalMatrix.cols() + columnMatrix.cols());
    originalMatrix.rightCols(columnMatrix.cols()) = columnMatrix;
}
void deleteColumn(MatrixXd& matrix, int columnIndex) {
    if (columnIndex < 0 || columnIndex >= matrix.cols()) {
        cerr << "Error: Column index out of bounds." << endl;
        return;
    }
    for (int i = columnIndex; i < matrix.cols() - 1; ++i) {
        matrix.col(i) = matrix.col(i + 1);
    }
    matrix.conservativeResize(matrix.rows(), matrix.cols() - 1);
}
MatrixXd selectRandomColumn(const MatrixXd& matrix) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(0, matrix.cols() - 1);
    int columnIndex = dis(gen);
    return matrix.col(columnIndex);
}
int main(){
    setValueOfn();
    MatrixXd randomMatrixA , randomMatrixB;
    randomMatrixA.setRandom(n,n);
    randomMatrixB.setRandom(n,n);
    MatrixXd symmetricMatrixA = makeSymmetric(randomMatrixA);
    MatrixXd symmetricMatrixB = makeSymmetric(randomMatrixB);
    //cout << "Symmetric Matrix A:\n" << symmetricMatrixA << endl;
    //cout << "Symmetric Matrix B:\n" << symmetricMatrixB << endl;
    //Create <double> , dynamic matrix
    MatrixXcd eigenValMatA , eigenValMatB;
    MatrixXcd eigenVectMatA , eigenVectMatB;
    //Finding EigenValues and Eigen vectors
    EigenSolver<MatrixXd> eigenValueSolverA(symmetricMatrixA);
    EigenSolver<MatrixXd> eigenValueSolverB(symmetricMatrixB);
    eigenValMatA=eigenValueSolverA.eigenvalues();
    eigenValMatB=eigenValueSolverB.eigenvalues();
    eigenVectMatA=eigenValueSolverA.eigenvectors();
    eigenVectMatB=eigenValueSolverB.eigenvectors();
    //cout<<"eigenValMatA:\n"<<eigenValMatA<<endl;
    //cout<<"eigenVectMatA:\n"<<eigenVectMatA<<endl;

    MatrixXd EigenVector(3,1);
    insertColumnMatrix(EigenVector , eigenVectMatA.real());
    insertColumnMatrix(EigenVector , eigenVectMatB.real());
    deleteColumn(EigenVector , 0);
    cout<<"Eigen Vector Matrix:\n"<<EigenVector<<endl;

    //Select a random eigen vector
    MatrixXcd randomEigenVectA = getRandomEigenVector(eigenVectMatA);
    MatrixXcd randomEigenVectB = getRandomEigenVector(eigenVectMatB);
    //Separate real parts of eigen vectors
    MatrixXd eigenVectorA , eigenVectorB;
    eigenVectorA=randomEigenVectA.real();
    eigenVectorB=randomEigenVectB.real();  
    // cout<<eigenVectorA<<endl;
    eigenVectorA=normalizeEigenVector(eigenVectorA);
    eigenVectorB=normalizeEigenVector(eigenVectorB);
    cout<<"Random Eigen Vector from A:\n";
    cout<<eigenVectorA<<endl;
    cout<<eigenVectorA.size()<<endl;
    int T=3;  //just for checking otherwise set vale to 10001
    MatrixXd yj;
    MatrixXd demo;
    MatrixXd Reward;
    MatrixXd operation(n, n);
    MatrixXd result;
    MatrixXd rewCala , rewCalaR , rewResl, rewCalb , rewCalbR , rewResR;
    MatrixXd Penalties(3,1);
    Penalties << 0,0,0;
    MatrixXd PenA , PenB , PenC , PenVecA , PenVecB;
    double PenAScal, PenBScal , PenCScal;
    double rewa , rewb;
    MatrixXd PenARes , PenBRes , diffRes;
// for (int j = 0; j < T; j++) {
//     for (int i = 0; i < n; i++) { 
        MatrixXd vi = selectRandomColumn(EigenVector);
        MatrixXd vj = selectRandomColumn(EigenVector);
        cout<<"Vi\n"<<vi<<endl;
        cout<<"Vj\n"<<vj<<endl;

        demo = eigenVectorB.transpose() * symmetricMatrixB * eigenVectorB;
        demo=demo.inverse();
        yj=eigenVectorB * demo ;
        rewCala=eigenVectorA.transpose() * symmetricMatrixB * eigenVectorA;
        rewa=rewCala(0,0);
        rewCalaR=symmetricMatrixA * eigenVectorA;
        rewCalb=eigenVectorA.transpose() * symmetricMatrixA * eigenVectorA;
        rewb=rewCalb(0,0);
        rewCalbR=symmetricMatrixB * eigenVectorA;
        rewResl=rewCalaR * rewa;
        rewResR=rewCalbR * rewb;
        Reward=rewResl - rewResR;
        //while(j < i){  Think on it
            PenA = eigenVectorA.transpose() * symmetricMatrixA * yj;
            PenAScal=PenA(0,0);
            PenB = eigenVectorA.transpose() * symmetricMatrixB * eigenVectorA ;
            PenBScal = PenB(0,0);
            PenC = eigenVectorA.transpose() * symmetricMatrixB * yj ;
            PenCScal = PenC(0,0);
            PenVecA = symmetricMatrixB * yj;
            PenVecB = symmetricMatrixB * eigenVectorA;
            PenARes = PenBScal * PenVecA;
            PenBRes = PenCScal * PenVecB;
            diffRes = PenARes - PenBRes;
            Penalties = Penalties + (PenAScal * diffRes);
        // }
    return 0;
}