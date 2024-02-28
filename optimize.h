//
// Created by xiaowuga on 2/26/2024.
//

#ifndef DEFILLET_OPTIMIZE_H
#define DEFILLET_OPTIMIZE_H

#include "common.h"
#include <easy3d/core/surface_mesh.h>

namespace DEFILLET {
    class Optimize {
    public:
        Optimize(const std::vector<Eigen::Vector3d>& points,
                 const std::vector<std::vector<size_t>>& faces,
                 const std::vector<size_t>& point_ancestors,
                 const std::vector<size_t>& face_ancestors,
                 const std::vector<Eigen::Vector3d>& point_tar_normals,
                 const std::vector<Eigen::Vector3d>& face_tar_normals,
                 const std::vector<size_t>& fixed_points,
                 std::string type = "edge-based",
                 double beta = 10);

        bool solve();
        void remesh();
        void get_points(std::vector<Eigen::Vector3d>& points);
        void get_faces(std::vector<std::vector<size_t>> & faces);
    private:
        void edge_init(const std::vector<size_t>& face_ancestors,
                       const std::vector<Eigen::Vector3d>& face_tar_normals,
                       const std::vector<size_t>& fixed_points);

        void point_init(const std::vector<Eigen::Vector3d>& normals,
                        const std::vector<size_t>& fixed_points);

        void face_init(const std::vector<Eigen::Vector3d>& normals,
                       const std::vector<size_t>& fixed_points);

    private:
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
        double beta_;
        easy3d::SurfaceMesh* mesh_;
        Eigen::VectorXd d_;

    };
}

#endif //DEFILLET_OPTIMIZE_H
