//
// Created by 13900K on 2024/3/18.
//

#ifndef DEFILLET_DEFILLET_H
#define DEFILLET_DEFILLET_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

class DeFillet {
public:
    DeFillet() {
        angle_ = 45;
        mesh_ = nullptr;
        fillet_mesh_ = nullptr;
        beta_ = 1.0;
        gamma_ = 1.0;
        num_opt_iter_ = 5;
    }

    void run_geodesic();

    void refine_target_normal();

    bool run_defillet();

    void init_opt();

    bool opt();

    void extract_fillet_region();

    // geter and seter
    void set_mesh(easy3d::SurfaceMesh* mesh) {
        mesh_ = mesh;
        box = mesh_->bounding_box();
    }
    void set_fillet_mesh(easy3d::SurfaceMesh* fillet_mesh) {
        fillet_mesh_ = fillet_mesh;
    }
    void set_angle(double angle) {angle_ = angle;}
    void set_beta(double beta) { beta_ = beta;}
    void set_gamma(double gamma) {gamma_ = gamma;}
    void set_num_opt_iter(double num_opt_iter) { num_opt_iter_ = num_opt_iter;}

    easy3d::SurfaceMesh* get_fillet_mesh() { return fillet_mesh_;}

    double get_geodesic_time() { return geodesic_time_;}
    double get_defillet_init_time() { return defillet_init_time_;}
    double get_defillet_iter_time() { return defillet_iter_time_;}
    double get_target_normals_refine_time() { return target_normals_refine_time_;}
private:
    easy3d::Box3 box;
    easy3d::SurfaceMesh* mesh_;
    easy3d::SurfaceMesh* fillet_mesh_;
    std::vector<int> sources_;
    std::vector<easy3d::vec3> sources_normals_;
    Eigen::VectorXd d_;

    double angle_;
    double beta_;
    double gamma_;
    int num_opt_iter_;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;

    double geodesic_time_;
    double boundary_refine_time_;
    double target_normals_refine_time_;
    double defillet_init_time_;
    double defillet_iter_time_;

};


#endif //DEFILLET_DEFILLET_H
