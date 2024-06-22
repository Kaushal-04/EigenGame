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
float calculateSigma(const MatrixXf& matrix) {
    EigenSolver<MatrixXf> eigenValueSolver(matrix);
    MatrixXf eigValmatrix = (eigenValueSolver.eigenvalues()).real();
    float maxAbsValue = eigValmatrix.cwiseAbs().maxCoeff();
    return 1 / maxAbsValue;

}
MatrixXf eigenValues(const MatrixXf& eigenvectors) {
    SelfAdjointEigenSolver<MatrixXf> solver(eigenvectors.transpose() * eigenvectors);
    return solver.eigenvalues();
}
int main() {
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
    MatrixXf w = B * EigenVector;
    int T=1; //3 for just check
    int M=2; //check it
    vector<MatrixXf> matA;
    MatrixXf tempMat;
    for(int i=0; i<M; i++){
        tempMat.setRandom(n,n);
        tempMat = tempMat.array().abs(); // Ensure positive values
        matA.push_back(tempMat);
    }
    vector<MatrixXf> matB;
    for(int i=0; i<M; i++){
        tempMat.setRandom(n,n);
        tempMat = tempMat.array().abs(); // Ensure positive values
        matB.push_back(tempMat);
    }
    MatrixXf Atm = MatrixXf::Zero(n,n);
    MatrixXf Btm = MatrixXf::Zero(n,n);
    MatrixXf rewCala , rewCalaR , rewResl, rewCalb , rewCalbR , rewResR;
    float rewa , rewb;
    MatrixXf Reward;
    int colNum;
    MatrixXf vi , vj;
    for(int j=0; j<T; j++){
        for(int i=0; i<n; i++){
            MatrixXf Deltaim , DeltaBvim;
            MatrixXf DelAvg(n,1) ;
                for(int row=0;  row<n; row++){
                    DelAvg(row , 0) = 0;
                }
            MatrixXf DelBvAvg(n,1) ;
                for(int row=0;  row<n; row++){
                    DelBvAvg(row , 0) = 0;
                }
            for(int m=0; m<M; m++){
                Atm += matA[m];
                Atm /= m+1;
                Btm += matB[m];
                Btm /= m+1;
                colNum = rand()%EigenVector.cols();
                vi = EigenVector.col(colNum);
                colNum = rand()%EigenVector.cols();
                vj = EigenVector.col(colNum);

                float result = (vj.transpose() * B * vj)(0, 0);
                if(result < 0) //to fin nan result
                    result = result * (-1);
                float sigma = calculateSigma(B);
                float maxi = max(result , sigma);
                float sqrtResult = sqrt(maxi);
                sqrtResult = 1.0 /sqrtResult;
                MatrixXf yj;
                yj = vj * sqrtResult;
                MatrixXf Byj = (B * vj )* sqrtResult;
                //Calculate Reward
                rewCala= vi.transpose() * Btm * vi;
                rewa=rewCala(0,0);
                rewCalaR = Atm * vi;
                rewCalb = vi.transpose() * Atm * vi;
                rewb = rewCalb(0,0);
                rewCalbR = Btm * vi;
                rewResl = rewCalaR * rewa;
                rewResR = rewCalbR * rewb;
                Reward = rewResl - rewResR;
                //Calculate Penalties
                int tempj = j;
                MatrixXf Ba , Bb;
                float BbL;
                MatrixXf diffRes;
                MatrixXf Penalties(n,1) ;
                for(int row=0;  row<n; row++){
                    Penalties(row , 0) = 0;
                }
                float PenA , PenB;
                while(tempj < i){
                    PenA = (vi.transpose() * Atm * yj)(0,0);
                    PenB = (vi.transpose() * Btm * vi)(0,0) ;
                    Ba = PenB * Byj ;
                    BbL = (vi.transpose() * Byj)(0,0);
                    Bb = (Btm * vi ) * BbL ;
                    diffRes = Ba - Bb;
                    Penalties = Penalties + (PenA * diffRes);
                    tempj++;
                }
                Deltaim = Reward - Penalties;
                DeltaBvim = (Btm * vi);
                int randomCol = rand()%w.cols();
                MatrixXf tempi = w.col(randomCol);
                DeltaBvim -= tempi;
                //Average of Deltaim
                DelAvg += Deltaim;
                DelBvAvg += DeltaBvim;
            } 
            DelAvg = DelAvg / M;      
            MatrixXf newVi = vi + (-2 * DelAvg);
            vi = normalizeVector(newVi);
            insertColumnMatrix(finalEigVect , vi);
            DelBvAvg = DelBvAvg / M;
            w.col(i) = w.col(i) + (-2 * DelBvAvg);
        }
    }
    deleteColumn(finalEigVect , 0);
    cout<<"Final Eigen Vector\n"<<finalEigVect<<endl;
    MatrixXf EigenValues = eigenValues(finalEigVect);
    cout<<"Eigen Values Matrix : \n"<<EigenValues<<endl;
    return 0;
}
