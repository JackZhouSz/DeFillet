//
// Created by xiaowuga on 3/11/2024.
//

#ifndef DEFILLET_FILLET_SEG_H
#define DEFILLET_FILLET_SEG_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

class FilletSeg {
public:
    FilletSeg(){
        eps_ = 0.03; s_ = 10;
        scoring_time_ = 0.0; gcp_time_ = 0.0; sor_time_ = 0.0;
        radius_ = 0.05;
        alpha_ = 0.5;
        angle_ = 40;
        min_score_ = 0.5;
        mesh_ = nullptr;
        std_ratio_ = 0.3;
        nb_neighbors_ = 30;
        num_sor_iter_ = 3;
        w_convex_ = 0.08;
        w_concave_ = 1.0;
        w1_ = 0.3;
        w2_ = 0.4;
    }

    void reset() {
        eps_ = 0.03; s_ = 10;
        scoring_time_ = 0.0; gcp_time_ = 0.0; sor_time_ = 0.0;
        radius_ = 0.1;
        alpha_ = 0.5;
        angle_ = 40;
        min_score_ = 0.5;
        box.clear();
        std_ratio_ = 1.0;
        nb_neighbors_ = 30;
        num_sor_iter_ = 3;
        mesh_ = nullptr;
        w_convex_ = 0.08;
        w_concave_ = 1.0;
        w1_ = 0.3;
        w2_ = 0.4;
    }

    void seg();

    // main steps
    void run_sor();
    void run_scoring();
    void run_gcp();

    // tools
    double cal_vertex_score(int vid);
    void boardcast_vertex_domain(int vid, double  score);


    // geter and seter
    void set_mesh(easy3d::SurfaceMesh* mesh) {
        mesh_ = mesh;
        box = mesh_->bounding_box();
    }
    void set_eps(double eps) {eps_ = eps;}
    void set_s(double s) {s_ = s;}
    void set_radius(double radius) { radius_ = radius;}
    void set_min_score(double min_score) { min_score_ = min_score;}
    void set_alpha(double alpha) {alpha_ = alpha;}
    void set_angle(double angle) {angle_ = angle;}
    void set_num_sor_iter(int num_sor_iter) { num_sor_iter_ = num_sor_iter; }
    void set_std_ratio(double std_ratio) { std_ratio_ = std_ratio; }
    void set_nb_neighbors(double nb_neighbors) { nb_neighbors_ = nb_neighbors; }
    void set_w_convex(double w_convex) {w_convex_ = w_convex;}
    void set_w_concave(double w_concave) {w_concave_ = w_concave;}
    void set_w1(double w1) {w1_ = w1;}
    void set_w2(double w2) {w2_ = w2;}

    double get_runtime() {
        return scoring_time_ + gcp_time_;
    }
    double get_sor_time() const { return sor_time_;}
    double get_scoring_time() const { return scoring_time_;}
    double get_gcp_time() const { return gcp_time_;}


    const std::vector<float>& get_vertices_scores() const { return vertices_scores_; }
    const std::vector<easy3d::vec3>& get_sites() const {return sites_;}
    const std::vector<easy3d::vec3>& get_vertices() const {return vertices_;}
    const easy3d::Box3& get_box() const { return box;}
    const std::vector<int> get_vertices_index_inside_box() const {
        std::vector<int> res;
        int nb_vertices = vertices_.size();
        for(int i = 0; i < nb_vertices; i++) {
            if(box.contains(vertices_[i])) {
                res.emplace_back(i);
            }
        }
        return res;
    }
    const std::vector<int> get_sor_labels() const { return sor_labels_; };

private:
    easy3d::SurfaceMesh* mesh_;
    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vertices_;
    std::vector<std::vector<int>> site_of_vertices_;

    std::vector<float> vertices_scores_;
    std::vector<int> sor_labels_;

    easy3d::Box3 box;
    int num_sor_iter_;
    double std_ratio_;
    int nb_neighbors_;
    double eps_;
    double s_;
    double radius_;
    double min_score_;
    double alpha_;
    double angle_;
    double w_convex_;
    double w_concave_;
    double w1_, w2_;
    double sor_time_;
    double scoring_time_;
    double gcp_time_;

};


#endif //DEFILLET_FILLET_SEG_H
