//
// Created by xiaowuga on 1/1/2024.
//
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

int main() {
    // 定义稀疏矩阵 A 和向量 b
    Eigen::SparseMatrix<double> A(3, 3); // 3x3 稀疏矩阵
    Eigen::VectorXd b(3);
    b << 1, 0, 1;// 长度为3的向量
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.emplace_back(Eigen::Triplet<double>(0, 0 , 1));
    triplets.emplace_back(Eigen::Triplet<double>(0, 1 , -1));
    triplets.emplace_back(Eigen::Triplet<double>(0, 2 ,  1));

    triplets.emplace_back(Eigen::Triplet<double>(1, 0 ,  2));
    triplets.emplace_back(Eigen::Triplet<double>(1, 1 ,  -1));
    triplets.emplace_back(Eigen::Triplet<double>(1, 2 ,  -1));

    triplets.emplace_back(Eigen::Triplet<double>(2, 0 ,  1));
    triplets.emplace_back(Eigen::Triplet<double>(2, 1 ,  1));
    triplets.emplace_back(Eigen::Triplet<double>(2, 2 ,  -1));

    A.setFromTriplets(triplets.begin(), triplets.end());
    // 填充稀疏矩阵 A 和向量 b 的数据，这里只是示例，需要提供具体的数据

    // 设置初始解向量
    Eigen::VectorXd initialGuess(3); // 长度为3的向量作为初始解
    initialGuess << 1,0,1;        // 设置初始解向量的初始值，可以根据需要修改

    // 定义共轭梯度法求解器
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    // 设置求解器的参数
    solver.setMaxIterations(1000); // 设置最大迭代次数
    solver.setTolerance(1e-6);    // 设置迭代收敛容许误差
    // 求解线性方程组 Ax = b，使用设置了初始值的解向量 initialGuess
    Eigen::VectorXd x = solver.compute(A).solveWithGuess(b, initialGuess);
    std::cout << solver.error() << std::endl;
    if (solver.info() != Eigen::Success) {
        // 求解失败
        std::cout << "Eigen CG solver failed to converge." << std::endl;

        return 1;
    }

    // 打印解向量 x
    std::cout << "Solution vector x:\n" << x << std::endl;




    return 0;
}

//#include <iostream>
//#include <Eigen/Dense>

//int main() {
//    // 定义矩阵 A 和向量 b
//    Eigen::Matrix3d A; // 3x3 矩阵
//    Eigen::Vector3d b; // 长度为3的向量
//
//    // 填充矩阵 A 和向量 b 的数据，这里只是示例，需要提供具体的数据
//    A << 1, -1, 1,
//            2, -1, -1,
//            1, 1, -1;
//
//    b << 1, 0, 1;
//
//    // 解线性方程 Ax = b
//    Eigen::Vector3d x = A.lu().solve(b);
//
//    // 打印解向量 x
//    std::cout << "Solution vector x:\n" << x << std::endl;
//
//    return 0;
//}

