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
    float maxElement = vect.cwiseAbs().maxCoeff();
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
    MatrixXf finalEigVect(n,1);
    MatrixXf randomMatrixA , randomMatrixB;
    randomMatrixA.setRandom(n,n);
    randomMatrixA = randomMatrixA.array().abs(); // Ensure positive values
    randomMatrixB.setRandom(n,n);
    randomMatrixB = randomMatrixB.array().abs(); // Ensure positive values
    MatrixXf A = makeSymmetric(randomMatrixA);
    MatrixXf B = makeSymmetric(randomMatrixB);
          // Eigenvalue and Eigenvector Computation: The GeneralizedEigenSolver class in Eigen is used
          // to solve the generalized eigenvalue problem Ax=λBx. The .compute(A, B) method computes the 
          // eigenvalues and eigenvectors.
    GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(A , B);
    MatrixXf eigenvalues = ges.eigenvalues();
    MatrixXf eigenvectors = ges.eigenvectors();
    for (int i = 0; i < eigenvectors.rows(); ++i) {
        for (int j = 0; j < eigenvectors.cols(); ++j) {
            eigenvectors(i, j) *= -1;
        }
    }
    cout << "Eigenvalues: \n" << eigenvalues << endl;
    cout << "Eigenvectors: \n" << eigenvectors << endl;

    MatrixXf EigenVector(n,1);
    MatrixXf temp;
    for(int i=0; i<eigenvectors.cols(); i++){
        temp = normalizeVector(eigenvectors.col(i));
        insertColumnMatrix(EigenVector , temp);
        deleteColumn(temp , 0);
    }
    deleteColumn(EigenVector , 0);
    cout<<"Normalized Eigen Vector\n"<<EigenVector<<endl;


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

        float result = (vj.transpose() * B * vj)(0, 0);
        if(result < 0) //to fin nan result
            result = result * (-1);
        float sqrtResult = sqrt(result);
        sqrtResult = 1.0 /sqrtResult;
        MatrixXf yj;
        yj = vj / sqrtResult;
        rewCala= vi.transpose() * B * vi;
        rewa = rewCala(0,0);
        rewCalaR = A * vi;
        rewCalb = vi.transpose() * A * vi;
        rewb=rewCalb(0,0);
        rewCalbR = B * vi;
        rewResl=rewCalaR * rewa;
        rewResR=rewCalbR * rewb;
        Reward=rewResl - rewResR;
        int tempj = j;
        MatrixXf PenARes , PenBRes , diffRes;
        MatrixXf Penalties(n,1) ;
        for(int row=0;  row<n; row++){
            Penalties(row , 0) = 0;
        }
        float PenAScal, PenBScal , PenCScal;
        MatrixXf PenA , PenB , PenC , PenVecA , PenVecB;
        while(tempj < i){
            PenA = vi.transpose() * A * yj;
            PenAScal=PenA(0,0);
            PenB = vi.transpose() * B * vi ;
            PenBScal = PenB(0,0);
            PenC = vi.transpose() * B * yj ;
            PenCScal = PenC(0,0);
            PenVecA = B * yj;
            PenVecB = B * vi;
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
            insertColumnMatrix(finalEigVect , vi);
            tempj++;
         }
     }
 }
    deleteColumn(finalEigVect , 0);
    cout<<"Final Eigen Vector:\n"<<finalEigVect<<endl;
    MatrixXf EigenValues = eigenValues(finalEigVect);
    cout<<"Eigen Values Matrix : \n"<<EigenValues<<endl;
    return 0;
}
