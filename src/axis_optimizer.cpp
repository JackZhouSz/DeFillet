//
// Created by xiaowuga on 2025/8/25.
//

#include <axis_optimizer.h>


AxisOptimizer::AxisOptimizer(const std::vector<vec3>& vectors, const std::vector<vec4> points, float weight, int max_iters, float learning_rate) {
    int num_vector = vectors.size();
    vectors_.resize(num_vector, 3);
    for(int i = 0; i < num_vector; i++) {
        vectors_.row(i) << vectors[i].x, vectors[i].y, vectors[i].z;
    }
    int num_points = points.size();
    Eigen::MatrixXd X(num_points, 3);
    for(int i = 0; i < num_points; i++) {
        X.row(i) << points[i].x,points[i].y, points[i].z;
    }
    Eigen::RowVector3d mu = X.colwise().mean();
    points_ = X.rowwise() - mu;

    principal_component_ = computePCA(points_);
    orthogonal_matrix_ = gramSchmidt(vectors_);

}


VectorXd AxisOptimizer::computePCA(const MatrixXd& points) {
    // 计算协方差矩阵
    MatrixXd cov_matrix = points.transpose() * points / points.rows();

    // 特征值分解
    SelfAdjointEigenSolver<MatrixXd> solver(cov_matrix);

    // 返回主成分（特征向量）
    return solver.eigenvectors().col(cov_matrix.cols() - 1);
}


MatrixXd AxisOptimizer::gramSchmidt(const MatrixXd& vectors) {
    int k = vectors.rows();
    int n = vectors.cols();

    MatrixXd orthogonal_matrix(k, n);

    for (int i = 0; i < k; ++i) {
        VectorXd vi = vectors.row(i).transpose();

        for (int j = 0; j < i; ++j) {
            double dot = orthogonal_matrix.row(j).transpose().dot(vi);
            vi -= dot * orthogonal_matrix.row(j).transpose();
        }

        orthogonal_matrix.row(i) = vi.transpose();
    }

    return orthogonal_matrix;
}

VectorXd AxisOptimizer::optimize() {
    // 使用PCA的主成分作为初始向量
    VectorXd w = principal_component_;
    w.normalize();  // 确保初始向量长度为1

    for (int iter = 0; iter < max_iters_; ++iter) {
        // 计算目标函数的梯度
        VectorXd grad(w.size());

        // 计算正交性梯度（简单近似）
        grad.setZero();
        for (int i = 0; i < orthogonal_matrix_.rows(); ++i) {
            grad -= weight_ * (orthogonal_matrix_.row(i).transpose() * std::abs(orthogonal_matrix_.row(i).transpose().dot(w)));
        }

        // 计算对齐性梯度
        grad += (1 - weight_) * (principal_component_ / (w.norm() * principal_component_.norm()) - w.transpose() * principal_component_ / w.norm());

        // 更新w
        w -= learning_rate_ * grad;
        w.normalize();  // 确保每次更新后w的长度保持为1
    }

    return w;
}

vec3 AxisOptimizer::getOptimalAxis() {
    VectorXd res = optimize();
    return vec3(res[0], res[1], res[2]);
}