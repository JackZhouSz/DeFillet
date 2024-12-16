//
// Created by xiaowuga on 2024/12/15.
//

#ifndef DEFILLET_V2_H
#define DEFILLET_V2_H

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

class DeFilletv2 {
public:
    DeFilletv2(easy3d::SurfaceMesh* mesh
             , double angle_thr = 40
             , double beta = 1.0
             , double gamma = 1.0
             , int num_opt_iter = 5);
    void run();
    void initialize();
    void optimize();

    double dihedral_angle(const easy3d::vec3& n1
                                 , const easy3d::vec3& n2
                                 , bool rad = false);
public:
    easy3d::SurfaceMesh* mesh_;
    easy3d::SurfaceMesh* fillet_mesh_;
    easy3d::SurfaceMesh* non_fillet_mesh_;
    easy3d::SurfaceMesh* focus_area_;
    std::vector<easy3d::vec3> tar_nor; // target normals
    std::vector<easy3d::vec3> s_; // sites;

    double angle_thr_;
    double beta_;
    double gamma_;
    int num_opt_iter_;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
};



#endif //DEFILLET_V2_H
