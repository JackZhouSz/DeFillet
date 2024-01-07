//
// Created by 13900K on 2024/1/7.
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>

int main() {
    Eigen::SparseMatrix<double> A(3, 3);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.emplace_back(Eigen::Triplet<double>(0, 0, 1));
    triplets.emplace_back(Eigen::Triplet<double>(1, 1, 2));
    triplets.emplace_back(Eigen::Triplet<double>(2, 2, 3));
    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
//    Eigen::SparseMatrix<double> Q = Eigen::KroneckerProductSparse<Eigen::SparseMatrix<double>, Eigen::MatrixXd>(A, I);
    Eigen::SparseMatrix<double> Q = Eigen::KroneckerProductSparse<Eigen::MatrixXd, Eigen::SparseMatrix<double>>(I, A);
    for(int i = 0; i < Q.rows(); i++) {
        for(int j = 0; j < Q.cols(); j++) {
            std::cout << Q.coeff(i, j) << ' ';
        }
        std::cout << std::endl;
    }
    return 0;
}
