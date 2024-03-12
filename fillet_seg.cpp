//
// Created by xiaowuga on 3/11/2024.
//

#include "fillet_seg.h"
#include "MeshVoronoi3d.h"

void FilletSeg::seg() {
    site_scoring();
    run_gcp();
    run_geodesic();
    refine_fillet_boundary();
}

void FilletSeg::site_scoring() {
    MeshVoronoi3d mv3d(mesh_);
    sites_ = mv3d.get_sites();
    vertices_ = mv3d.get_vertices();
    site_of_vertices_ = mv3d.get_vertices2sites();
    int nb_sites = sites_.size();
    int nb_vertices = vertices_.size();
    auto score = mesh_->face_property<double>("f:score", 0.0);
    auto count = mesh_->face_property<int>("f:count", 0);
    vertex_scores_.resize(nb_vertices, 0);
    double len = radius_ * box.diagonal_length();

    for(int i = 0; i < nb_vertices; i++) {
        double R = (vertices_[i] - sites_[site_of_vertices_[i][0]]).norm();
        if(R < len && box.contains(vertices_[i])) {
            double score = cal_vertex_score(i);
            if(score > min_score_) {
                vertex_scores_[i] = score;
                boardcast_vertex_domain(i, score);
            }
        }
    }

}

void FilletSeg::run_gcp() {

}

void FilletSeg::run_geodesic() {

}

void FilletSeg::refine_fillet_boundary() {

}

void FilletSeg::refine_target_normal() {

}

double FilletSeg::cal_vertex_score(int vid) {
    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;

    std::map<int,int> root;

    for(auto f : site_of_vertices_[vid]) {
        easy3d::SurfaceMesh::Face face_handle(f);
        que.push(face_handle);
        root[f] = f;
    }

    std::function<int(int)> get_root = [&](int k)-> int {
        return root[k] == k ? k : (root[k] = get_root(root[k]));
    };

    double R = (sites_[site_of_vertices_[vid][0]] - vertices_[vid]).norm();
    double tol = R * eps_;
    double err = 0.0;

}

void FilletSeg::boardcast_vertex_domain(int vid, double score) {

}



