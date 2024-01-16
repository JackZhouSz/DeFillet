//
// Created by 13900K on 2024/1/16.
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>


int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> C;
    V.resize(4, 3);
    V.row(0) = Eigen::Vector3d(0.0,1.0, 0.0);
    V.row(1) = Eigen::Vector3d(-1.0,0.0, 0.0);
    V.row(2) = Eigen::Vector3d(1.0,0.0, 0.0);
    V.row(3) = Eigen::Vector3d(0.0,-1.0, 0.0);

    F.resize(2, 3);
    F.row(0) = Eigen::Vector3i(0,1,2);
    F.row(1) = Eigen::Vector3i(1,3,2);

    igl::cotmatrix(V, F, C);
    for(int i = 0; i < C.rows(); i++) {
        for(int j = 0; j < C.cols(); j++) {
            std::cout << C.coeff(i, j) << ' ';
        }
        std::cout <<std::endl;
    }
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(-C);

    if (ldlt.info() == Eigen::Success)
    {
        Eigen::VectorXd eigenvalues = ldlt.vectorD();

        // 打印特征值
        std::cout << "Eigenvalues of the sparse matrix:\n" << eigenvalues << "\n";
    }
    else {
        std::cerr << "Decomposition failed! Matrix may not be positive definite.\n";
    }
    return 0;
}

