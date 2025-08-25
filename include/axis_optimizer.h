//
// Created by xiaowuga on 2025/8/25.
//

#ifndef AXISOPTIMIZER_H
#define AXISOPTIMIZER_H

#include <Eigen/Dense>
#include <easy3d/core/types.h>

using namespace easy3d;
using namespace Eigen;
using namespace std;

class AxisOptimizer {
public:
    AxisOptimizer(const std::vector<vec3>& vectors, const std::vector<vec4> points, float weight = 0, int max_iters = 0, float learning_rate = 0.01);

    VectorXd computePCA(const MatrixXd& points);

    MatrixXd gramSchmidt(const MatrixXd& vectors);

    double objectiveFunction(const VectorXd& w);

    VectorXd optimize();

    vec3 getOptimalAxis();
private:
    MatrixXd vectors_;              // 给定的向量集合
    MatrixXd points_;               // 给定的点集合
    VectorXd principal_component_;  // PCA主成分
    MatrixXd orthogonal_matrix_;    // 正交化后的向量集合
    float weight_;                 // 正交性和对齐性的加权系数
    int max_iters_;                 // 最大迭代次数
    float learning_rate_;          // 学习率
};



#endif //AXISOPTIMIZER_H
