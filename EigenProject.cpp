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
MatrixXd normalizeVector(const MatrixXd& vect) {
    MatrixXd normalizedVector = vect; 
    double maxElement = vect.maxCoeff();
    for (int i = 0; i < vect.size(); ++i) {
        normalizedVector(i) /= maxElement;
    }
    return normalizedVector;
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
MatrixXd eigenValues(const MatrixXd& eigenvectors) {
    SelfAdjointEigenSolver<MatrixXd> solver(eigenvectors.transpose() * eigenvectors);
    return solver.eigenvalues();
}
int main(){
    setValueOfn();
    MatrixXd FinalEigenVector(3,1);
    MatrixXd randomMatrixA , randomMatrixB;
    randomMatrixA.setRandom(n,n);
    randomMatrixA = randomMatrixA.array().abs(); // Ensure positive values
    randomMatrixB.setRandom(n,n);
    randomMatrixB = randomMatrixB.array().abs(); // Ensure positive values
    //cout<<"MAtrix A:\n"<<randomMatrixA<<endl;
    //cout<<"Matrix B:\n"<<randomMatrixB<<endl;
    // MatrixXd symmetricMatrixA = makeSymmetric(randomMatrixA);
    // MatrixXd symmetricMatrixB = makeSymmetric(randomMatrixB);

    MatrixXd symmetricMatrixA(3,3);
    symmetricMatrixA<< 1,2,3,
                        7,8,9,
                        4,5,6;
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
    cout<<"Symmetric Matrix:\n"<<symmetricMatrixA<<endl;
    cout<<"Eigen Value A:\n"<<eigenValMatA.real()<<endl;
    eigenVectMatA=eigenValueSolverA.eigenvectors();
    cout<<"Eigen Vector Matrix A:\n"<<eigenVectMatA.real()<<endl;
    MatrixXd eigenVectTemp = eigenVectMatA.real();
    for(int i=0; i<eigenVectTemp.cols(); i++){
        eigenVectTemp.col(i) = normalizeVector(eigenVectTemp.col(i));
    }
    cout<<"Normalize Eigen Vector:\n"<<eigenVectTemp<<endl; 
    return 0;



    eigenValMatB=eigenValueSolverB.eigenvalues();
    eigenVectMatB=eigenValueSolverB.eigenvectors();
    
    MatrixXd EigenVector(3,1);
    insertColumnMatrix(EigenVector , eigenVectMatA.real());
    insertColumnMatrix(EigenVector , eigenVectMatB.real());
    deleteColumn(EigenVector , 0);
    //cout<<"Eigen Vector Matrix:\n"<<EigenVector<<endl;

    int T=1;  //just for checking otherwise set vale to 10001
    MatrixXd Reward;
    MatrixXd operation(n, n);
    MatrixXd result;
    MatrixXd rewCala , rewCalaR , rewResl, rewCalb , rewCalbR , rewResR;
    double PenAScal, PenBScal , PenCScal;
    double rewa , rewb;
for (int j = 0; j < T; j++) {  // loop to T
    for (int i = 0; i < 2; i++) {  // loop to k or n
        MatrixXd vi = selectRandomColumn(EigenVector);
        MatrixXd vj = selectRandomColumn(EigenVector);
        vi = normalizeVector(vi);
        vj = normalizeVector(vj);
        //cout<<"Vi\n"<<vi<<endl;
        //cout<<"Vj\n"<<vj<<endl;
        double result = (vj.transpose() * symmetricMatrixB * vj)(0, 0);
        //cout<<"Result:"<<result<<endl;
        double sqrtResult = sqrt(result);
        //why nan (Not a number is encountered
        sqrtResult = 1.0 /sqrtResult;
        //cout<<"Square Root : "<<sqrtResult<<endl;
        //double demoInverse = 1.0 / sqrtResult;
        //yj = vj * demoInverse ;
        MatrixXd yj;
        yj = vj * sqrtResult;
        //cout<<"yj\n"<<yj<<endl;
        rewCala= vi.transpose() * symmetricMatrixB * vi;
        rewa=rewCala(0,0);
        rewCalaR=symmetricMatrixA * vi;
        rewCalb=vi.transpose() * symmetricMatrixA * vi;
        rewb=rewCalb(0,0);
        rewCalbR=symmetricMatrixB * vi;
        rewResl=rewCalaR * rewa;
        rewResR=rewCalbR * rewb;
        Reward=rewResl - rewResR;
        //cout<<"Reward:\n"<<Reward<<endl;
        while(j < i){  
            MatrixXd PenARes , PenBRes , diffRes;
            MatrixXd Penalties(3,1);
            Penalties << 0,0,0;
            MatrixXd PenA , PenB , PenC , PenVecA , PenVecB;
            PenA = vi.transpose() * symmetricMatrixA * yj;
            PenAScal=PenA(0,0);
            //cout<<"PenA : "<<PenA<<endl;
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
            //cout<<"Penalties \n"<<Penalties<<endl;
            MatrixXd delta(3,1);
            delta = Reward - Penalties;
            MatrixXd itadelta(3,1);
            itadelta = -2 * delta;
            MatrixXd wi(3,1);
            wi = vi + itadelta;
            //cout<<"wi\n"<<wi<<endl;
            vi = normalizeVector(wi);
            //cout<<"New Vi\n"<<vi <<endl;
            insertColumnMatrix(FinalEigenVector , vi);
            }
        }
    }
    deleteColumn(FinalEigenVector , 0);
    //cout<<"Final Eigen Vector:\n"<<FinalEigenVector<<endl;
    
    MatrixXd EigenValues = eigenValues(FinalEigenVector);
    //cout<<"Eigen Values Matrix : \n"<<EigenValues<<endl;
    return 0;
}