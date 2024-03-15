//
// Created by xiaowuga on 3/11/2024.
//

#include "fillet_seg.h"
#include "MeshVoronoi3d.h"
#include <map>

void FilletSeg::seg() {
    run_scoring();
    run_gcp();
    run_geodesic();
    refine_fillet_boundary();
}

void FilletSeg::run_scoring() {
    MeshVoronoi3d mv3d(mesh_);
    sites_ = mv3d.get_sites();
    vertices_ = mv3d.get_vertices();
    site_of_vertices_ = mv3d.get_vertices2sites();
    int nb_sites = sites_.size();
    int nb_vertices = vertices_.size();
    auto scores = mesh_->face_property<double>("f:scores");
    auto counts = mesh_->face_property<int>("f:counts");
    for(auto f : mesh_->faces()) {
        scores[f] = 0.0;
        counts[f] = 0;
    }
    vertices_scores_.resize(nb_vertices, 0);
    double len = radius_ * box.diagonal_length();
    for(int i = 0; i < nb_vertices; i++) {
        double R = (vertices_[i] - sites_[site_of_vertices_[i][0]]).norm();
        if(R < len && box.contains(vertices_[i])) {
            double score = cal_vertex_score(i);
            vertices_scores_[i] = 0.0;
            if(score > min_score_) {
                vertices_scores_[i] = score;
                boardcast_vertex_domain(i, score);
            }
        } else {
            vertices_scores_[i] = 0.0;
        }
    }

    for(auto f : mesh_->faces()) {
        if(counts[f] > 0) {
            scores[f] /= counts[f];
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

    for(auto f : site_of_vertices_[vid]) {
        easy3d::SurfaceMesh::Face face_handle(f);
        que.push(face_handle);
    }


    double R = (sites_[site_of_vertices_[vid][0]] - vertices_[vid]).norm();
    double tol = R * eps_;
    double err = 0.0;
    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        if (vis.find(cur) != vis.end()) {
            continue;
        }
        vis.insert(cur);
        err += fabs((sites_[cur.idx()] - vertices_[vid]).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid()) {
                double len = (sites_[opp_f.idx()] - vertices_[vid]).norm();
                if (len < R + tol) {
                    if (vis.find(opp_f) == vis.end()) {
                        que.push(opp_f);
                    }
                }

            }
        }
    }

    std::set<easy3d::SurfaceMesh::Face> vis2;

    easy3d::SurfaceMesh::Face face_handle(site_of_vertices_[vid][0]);

    que.push(face_handle);
    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        if (vis2.find(cur) != vis2.end()) {
            continue;
        }
        vis2.insert(cur);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && vis.find(opp_f) != vis.end()) {
                que.push(opp_f);
            }
        }
    }
    if(vis2.size() != vis.size()) {
        return 0.0;
    }

    return exp(-s_ * err / vis.size());
}


void FilletSeg::boardcast_vertex_domain(int vid, double score) {
    auto scores = mesh_->get_face_property<double>("f:scores");
    auto counts = mesh_->get_face_property<int>("f:counts");

    std::set<easy3d::SurfaceMesh::Face> vis;
    int st = site_of_vertices_[vid][0];
    double axis_len = (vertices_[vid] - sites_[st]).norm();
    easy3d::SurfaceMesh::Face sf(st);
    std::queue<easy3d::SurfaceMesh::Face>que;
    que.push(sf);

    while(!que.empty()) {
        auto cur = que.front(); que.pop();
        if(vis.find(cur) == vis.end()) {
            vis.insert(cur);
        }
        else {
            continue;
        }
        for(auto h : mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if(opp_f.is_valid() && vis.find(opp_f) == vis.end()) {
                double len = (vertices_[vid] - sites_[opp_f.idx()]).norm();
                if(fabs(axis_len - len) / axis_len < eps_ && vis.find(opp_f) == vis.end()) {
                    que.push(opp_f);
                }
            }
        }
    }

    for (auto f: vis) {
        scores[f] += score;
        counts[f] += 1;
    }
}



