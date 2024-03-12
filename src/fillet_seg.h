//
// Created by xiaowuga on 3/11/2024.
//

#ifndef DEFILLET_FILLET_SEG_H
#define DEFILLET_FILLET_SEG_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

class FilletSeg {
public:
    FilletSeg(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
        eps_ = 0.03; s_ = 10;
        scoring_time_ = 0.0; gcp_time_ = 0.0;
        boundary_refine_time_ = 0.0; target_normals_refine_time_ = 0.0;
        radius_ = 0.1;
        box = mesh_->bounding_box();
    }

    void seg();

    // main steps
    void face_scoring();
    void run_gcp();
    void run_geodesic();
    void refine_fillet_boundary();
    void refine_target_normal();

    // tools
    double cal_vertex_score(int vid);
    void boardcast_vertex_domain(int vid, double  score);


    // geter and seter
    void set_eps(double eps) {eps_ = eps;}
    void set_s(double s) {s_ = s;}
    void set_radius(double radius) { radius_ = radius;}
    void set_min_score(double min_score) { min_score_ = min_score;}

    double get_runtime() {
        return scoring_time_ + gcp_time_ + geodesic_time_ + boundary_refine_time_ + target_normals_refine_time_;
    }
    double get_scoring_time() { return scoring_time_;}
    double get_gcp_time() { return gcp_time_;}
    double get_boundary_time() { return boundary_refine_time_;}
    double get_target_normals_refine_time() { return target_normals_refine_time_;}

    std::vector<double> get_vertices_scores() { return vertices_scores_; }

private:
    easy3d::SurfaceMesh* mesh_;
    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vertices_;
    std::vector<std::vector<int>> site_of_vertices_;

    std::vector<double> vertices_scores_;

    easy3d::Box3 box;
    double eps_;
    double s_;
    double radius_;
    double min_score_;

    double scoring_time_;
    double gcp_time_;
    double geodesic_time_;
    double boundary_refine_time_;
    double target_normals_refine_time_;
};


#endif //DEFILLET_FILLET_SEG_H
