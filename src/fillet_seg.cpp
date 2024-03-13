//
// Created by xiaowuga on 3/11/2024.
//

#include "fillet_seg.h"
#include "MeshVoronoi3d.h"

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
    auto scores = mesh_->face_property<double>("f:scores", 0.0);
    auto counts = mesh_->face_property<int>("f:counts", 0);
    vertices_scores_.resize(nb_vertices, 0);
    double len = radius_ * box.diagonal_length();

    for(int i = 0; i < nb_vertices; i++) {
        double R = (vertices_[i] - sites_[site_of_vertices_[i][0]]).norm();
        if(R < len && box.contains(vertices_[i])) {
            double score = cal_vertex_score(i);
            if(score > min_score_) {
                vertices_scores_[i] = score;
                boardcast_vertex_domain(i, score);
            }
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
    while(!que.empty()) {
        auto cur = que.front(); que.pop();
        if(vis.find(cur) != vis.end()) {
            continue;
        }
        vis.insert(cur);

        for(auto h : mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if(opp_f.is_valid()) {
                double len = (sites_[opp_f.idx()] - vertices_[vid]).norm();
                if(len < R + tol) {
                    if(vis.find(opp_f) == vis.end()) {
                        err += fabs(len - R);
                        que.push(opp_f);
                        int r = get_root(cur.idx());
                        root[opp_f.idx()] = r;
                    } else {
                        int r1 = get_root(cur.idx());
                        int r2 = get_root(opp_f.idx());
                        root[r1] = r2;
                    }
                }

            }
        }

        int nb_component = 0;
        for(auto item : root) {
            root[item.first] = get_root(item.first);
            if(item.first == item.second) {
                nb_component++;
                if(nb_component > 1) break;
            }
        }

        if(nb_component > 1) {
            return 0.0;
        }

        return exp(-s_ * err / vis.size());
    }
}

void FilletSeg::boardcast_vertex_domain(int vid, double score) {
    auto scores = mesh_->face_property<double>("f:scores");
    auto counts = mesh_->face_property<int>("f:counts");

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
                if(abs(axis_len - len) / axis_len < eps_ && vis.find(opp_f) == vis.end()) {
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



