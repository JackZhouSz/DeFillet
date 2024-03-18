//
// Created by 13900K on 2024/3/18.
//

#ifndef DEFILLET_DEFILLET_H
#define DEFILLET_DEFILLET_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

class DeFillet {
public:
    DeFillet() {
        angle_ = 45;

    }

    void run_geodesic();

    void refine_target_normal();

    void run_defillet();

    void optimize();

    void extract_interest_region();

    // geter and seter
    void set_mesh(easy3d::SurfaceMesh* mesh) {
        mesh_ = mesh;
        box = mesh_->bounding_box();
    }
    void set_angle(double angle) {angle_ = angle;}
    void set_beta(double beta) { beta_ = beta;}
    void set_num_opt_iter(double num_opt_iter) { num_opt_iter_ = num_opt_iter;}

    double get_geodesic_time() { return geodesic_time_;}
    double get_defillet_time() { return defillet_time_;}
    double get_target_normals_refine_time() { return target_normals_refine_time_;}
private:
    easy3d::Box3 box;
    easy3d::SurfaceMesh* mesh_;
    easy3d::SurfaceMesh* fillet_mesh_;
    double angle_;

    double beta_;
    int num_opt_iter_;

    double geodesic_time_;
    double boundary_refine_time_;
    double target_normals_refine_time_;
    double defillet_time_;

};


#endif //DEFILLET_DEFILLET_H
