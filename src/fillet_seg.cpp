//
// Created by xiaowuga on 3/11/2024.
//

#include "fillet_seg.h"
#include "MeshVoronoi3d.h"
#include "knn.h"
#include "gcp.h"

#include <map>
#include <omp.h>

void FilletSeg::seg() {
    run_scoring();
    run_gcp();
}

void FilletSeg::run_sor() {
    int nb_vertices = vertices_.size();
    sor_labels_.resize(nb_vertices);
    if(nb_neighbors_ < 1 || std_ratio_ <= 0) {
        std::cout << "Illegal input parameters, the number of neighbors and "
                     "standard deviation ratio must be positive." << std::endl;
    }
    auto start_time = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for(int i = 0; i < nb_vertices; i++) {
        if(box.contains(vertices_[i])) {
            sor_labels_[i] = 1;
        } else {
            sor_labels_[i] = 0;
        }
    }

    for(int i = 0; i < num_sor_iter_; i++) {
        std::vector<KNN::Point> knn_points;
        std::vector<int> indices;
#pragma omp parallel for
        for(int j = 0; j < nb_vertices; j++) {
            if(sor_labels_[j]) {
                #pragma omp critical
                {
                    knn_points.emplace_back(KNN::Point(vertices_[j].x,
                                                       vertices_[j].y, vertices_[j].z));
                    indices.emplace_back(j);
                }
            }
        }
        KNN::KdSearch kds(knn_points);
        int num = knn_points.size();
        size_t valid_distances = 0;
        std::vector<double> avg_distances(num);
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(knn_points[j], nb_neighbors_, tmp_indices, dist);
            double mean = -1.0;

            if(dist.size() > 0u) {
                valid_distances++;
                std::for_each(dist.begin(), dist.end(),
                              [](double &d) { d = std::sqrt(d); });
                mean = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();
            }
            avg_distances[j] = mean;
        }
        if(valid_distances == 0) {
            continue;
        }
        double cloud_mean = std::accumulate(
                avg_distances.begin(), avg_distances.end(), 0.0,
                [](double const &x, double const &y) { return y > 0 ? x + y : x; });

        cloud_mean /= valid_distances;
        double sq_sum = std::inner_product(
                avg_distances.begin(), avg_distances.end(), avg_distances.begin(),
                0.0, [](double const &x, double const &y) { return x + y; },
                [cloud_mean](double const &x, double const &y) {
                    return x > 0 ? (x - cloud_mean) * (y - cloud_mean) : 0;
                });
        double std_dev = std::sqrt(sq_sum / (valid_distances - 1));
        double distance_threshold = cloud_mean + std_ratio_ * std_dev;
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            if(avg_distances[j] > 0 && avg_distances[j] < distance_threshold) {
                sor_labels_[indices[j]] = 1;
            } else {
                sor_labels_[indices[j]] = 0;
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    sor_time_ = 1.0 * duration.count() / 1000000.0;
}

void FilletSeg::run_scoring() {
    omp_set_num_threads(10);

    MeshVoronoi3d mv3d(mesh_);
    sites_ = mv3d.get_sites();
    vertices_ = mv3d.get_vertices();
    site_of_vertices_ = mv3d.get_vertices2sites();
    int nb_sites = sites_.size();
    int nb_vertices = vertices_.size();
    auto scores = mesh_->face_property<float>("f:scores");
    auto counts = mesh_->face_property<int>("f:counts");
    run_sor();
    for(auto f : mesh_->faces()) {
        scores[f] = 0.0;
        counts[f] = 0;
    }
    vertices_scores_.resize(nb_vertices, 0);
    double len = radius_ * box.diagonal_length();
    auto start_time = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for (int i = 0; i < nb_vertices; i++) {
        double R = (vertices_[i] - sites_[site_of_vertices_[i][0]]).norm();
        if (R < len && sor_labels_[i]) {
            double score = cal_vertex_score(i);
            vertices_scores_[i] = 0.0;
            if (score > min_score_) {
                vertices_scores_[i] = score;
                boardcast_vertex_domain(i, score);
            }
        } else {
            vertices_scores_[i] = 0.0;
        }
    }
    int nb_faces = mesh_->n_faces();
#pragma omp parallel for
    for(int i = 0; i < nb_faces; i++) {
        easy3d::SurfaceMesh::Face f(i);
        if (counts[f] > 0) {
            if(counts[f] < 10) {
                scores[f] = 0.0;
            }
            else {
                scores[f] /= 1.0 * counts[f];
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    scoring_time_ = 1.0 * duration.count() / 1000000.0;
}


void FilletSeg::run_gcp() {
    auto start_time = std::chrono::high_resolution_clock::now();
    auto scores = mesh_->face_property<float>("f:scores");
    auto normals = mesh_->face_property<easy3d::vec3>("f:normal");
    auto gcp = mesh_->face_property<int>("f:gcp_labels");

    for(auto f : mesh_->faces()) {
        normals[f] = mesh_->compute_face_normal(f);
    }

    int nb_face = mesh_->n_faces();
    GCoptimizationGeneralGraph gc(nb_face, 2);

    std::vector<double> data_cost(2 * nb_face);
    for(auto f : mesh_->faces()) {
        data_cost[f.idx()] = scores[f];
        data_cost[f.idx()] = std::max(data_cost[f.idx()], w1_);
        data_cost[f.idx() + nb_face] = (1.0 - scores[f]);
        data_cost[f.idx() + nb_face] = std::max(data_cost[f.idx() + nb_face], w1_);
    }
    std::vector<std::pair<int,int>> edges;
    std::vector<double> edge_weights;
    for(auto e : mesh_->edges()) {
        auto f0 = mesh_->face(e, 0);
        auto f1 = mesh_->face(e, 1);
        if(f0.is_valid() && f1.is_valid()) {
            double dot_val = easy3d::dot(normals[f0], normals[f1]);
            dot_val = std::clamp(dot_val, -0.99999, 0.99999);
            double w = acos(dot_val) / M_PI;
            if(w < 0) {
                w = std::fabs(w) * w_convex_ * alpha_;
            } else {
                w = std::fabs(w) * w_concave_ * alpha_;
            }
            w = std::max(w, w2_);
            int id1 = f0.idx(), id2 = f1.idx();
            if(id1 > id2)
                std::swap(id1, id2);
            gc.setNeighbors(id1, id2);
            edges.emplace_back(std::make_pair(id1, id2));
            edge_weights.emplace_back(w);
        }
    }

    GCP::DataCost data_item(data_cost, nb_face, 2);
    GCP::SmoothCost smooth_item(edges, edge_weights);
    gc.setDataCostFunctor(&data_item);
    gc.setSmoothCostFunctor(&smooth_item);
    std::cout << "Before optimization energy is " << gc.compute_energy() << std::endl;
    gc.expansion(4);
    std::cout << "After optimization energy is " << gc.compute_energy() << std::endl;

#pragma omp parallel for
    for(int i = 0; i < nb_face; i++) {
        easy3d::SurfaceMesh::Face f(i);
        if(gc.whatLabel(i) == 0) {
            gcp[f] = 0;
        } else {
            gcp[f] = 1;
        }
    }

    auto sources = mesh_->vertex_property<int>("v:sources");

    for(auto e : mesh_->edges()) {
        auto f0 = mesh_->face(e, 0);
        auto f1 = mesh_->face(e, 1);
        auto v0 = mesh_->vertex(e, 0);
        auto v1 = mesh_->vertex(e, 1);
        if(f0.is_valid() && f1.is_valid()) {
            double dot_val = easy3d::dot(normals[f0], normals[f1]);
            dot_val = std::clamp(dot_val, -0.99999, 0.99999);
            double di_angle = acos(dot_val) / M_PI;
            if(di_angle > 0.5) {
                sources[v0] = 1; sources[v1] = 1;
            }
            else {
                sources[v0] = 0; sources[v1] = 0;
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    gcp_time_ = 1.0 * duration.count() / 1000000.0;
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
    auto scores = mesh_->get_face_property<float>("f:scores");
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



