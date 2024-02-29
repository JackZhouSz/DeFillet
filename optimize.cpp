//
// Created by xiaowuga on 2/26/2024.
//

#include "optimize.h"

#include <igl/cat.h>

#include "visualization.h"

namespace DEFILLET {
    Optimize::Optimize(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       const std::vector<size_t>& point_ancestors,
                       const std::vector<size_t>& face_ancestors,
                       const std::vector<Eigen::Vector3d>& point_tar_normals,
                       const std::vector<Eigen::Vector3d>& face_tar_normals,
                       const std::vector<size_t>& fixed_points,
                       std::string type,
                       double beta) : beta_(beta) {
        mesh_ = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh_->add_vertex(easy3d::vec3(points[i].x(),
                                          points[i].y(), points[i].z()));
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh_->add_triangle(v1, v2, v3);
        }
        int nb_edges = mesh_->n_edges();
        auto tar_nomrals = mesh_->add_face_property<easy3d::vec3>("f:tar_normals");
        auto src_points = mesh_->add_face_property<int>("f:src_idx");
        for(auto f : mesh_->faces()) {
            int idx = f.idx();
            tar_nomrals[f] = easy3d::vec3(face_tar_normals[idx].x(),
                                          face_tar_normals[idx].y(), face_tar_normals[idx].z());
            src_points[f] = face_ancestors[idx];
        }

        refine_tar_normals();

        std::vector<Eigen::Triplet<double>> triplets;
        int nb_fixed_points = fixed_points.size();
        d_ = Eigen::VectorXd(nb_fixed_points * 3);
        E.resize(nb_fixed_points * 3, nb_points * 3);
//        Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
        for(int i = 0; i < nb_fixed_points; i++) {
            int id = fixed_points[i];
            auto v = mesh_->position(easy3d::SurfaceMesh::Vertex(id));
            triplets.emplace_back(Eigen::Triplet<double>(3 * i, id , 1.0));
            d_[ 3 * i] = v.x;
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1, id + nb_points, 1.0));
            d_[3 * i + 1] = v.y;
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  id + nb_points * 2, 1.0));
            d_[3 * i + 2] = v.z;
        }
        E.setFromTriplets(triplets.begin(), triplets.end());
        nb_fixed_points = fixed_points.size();
        update_solver();
    }

    void Optimize::update_solver() {
        int nb_points = mesh_->n_vertices();
        int nb_edges = mesh_->n_edges();
        int nb_faces = mesh_->n_faces();
        auto tar_nomrals = mesh_->get_face_property<easy3d::vec3>("f:tar_normals");
        auto src_points = mesh_->get_face_property<int>("f:src_idx");
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::SparseMatrix<double> FNC(nb_edges * 2, nb_points * 3);
        for(auto e : mesh_->edges()) {
            int id = e.idx();
            int v0 = mesh_->vertex(e, 0).idx();
            int v1 = mesh_->vertex(e, 1).idx();
            auto f1 = mesh_->face(e, 0);
            double nx = tar_nomrals[f1].x, ny = tar_nomrals[f1].y, nz = tar_nomrals[f1].z;
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v0, nx));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v0 + nb_points, ny));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v0 + 2 * nb_points, nz));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v1, -nx));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v1 + nb_points, -ny));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id, v1 + 2 * nb_points, -nz));

            auto f2 = mesh_->face(e, 1);
            nx = tar_nomrals[f2].x, ny = tar_nomrals[f2].y, nz = tar_nomrals[f2].z;
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v0, nx));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v0 + nb_points, ny));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v0 + 2 * nb_points, nz));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v1, -nx));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v1 + nb_points, -ny));
            triplets.emplace_back(Eigen::Triplet<double>(2 * id + 1, v1 + 2 * nb_points, -nz));
        }
        FNC.setFromTriplets(triplets.begin(), triplets.end());
        triplets.clear();
        Eigen::SparseMatrix<double> FCC(nb_faces, nb_points * 3);
        for(auto f : mesh_->faces()) {
            std::vector<easy3d::vec3> pos;
            std::vector<int> pidx;
            double nx = tar_nomrals[f].x, ny = tar_nomrals[f].y, nz = tar_nomrals[f].z;
            for(auto v : mesh_->vertices(f)) {
                pos.emplace_back(mesh_->position(v));
                pidx.emplace_back(v.idx());
            }

            int num = pidx.size();
            easy3d::vec3 center = easy3d::vec3(0,0,0);
            for(int i = 0; i < num; i++) {
                center += pos[i];
            }


            for(int i = 0; i < num; i++) {
                int id = pidx[i];
                triplets.emplace_back(Eigen::Triplet<double>(f.idx(), id, 2.0 * nx / num));
                triplets.emplace_back(Eigen::Triplet<double>(f.idx(), id + nb_points, 2.0 * ny / num));
                triplets.emplace_back(Eigen::Triplet<double>(f.idx(), id + 2 * nb_points, 2 * nz / num));
            }
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), src_points[f], -nx * 2));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), src_points[f] + nb_points, -ny * 2));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), src_points[f] + 2 * nb_points, -nz * 2));

        }
        FCC.setFromTriplets(triplets.begin(), triplets.end());


//        igl::cat(1, FNC, FCC, A);
        Eigen::SparseMatrix<double> A;
        igl::cat(1, FNC, FCC, A);
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> ATA = AT * A;
        Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
        I = beta_ * I;
        Eigen::SparseMatrix<double> Q = I + ATA;
        Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
        Eigen::SparseMatrix<double> ET = E.transpose();
        Eigen::SparseMatrix<double> tempMat1;
        Eigen::SparseMatrix<double> tempMat2;
        Eigen::SparseMatrix<double> S;
        igl::cat(1, Q, E, tempMat1);
        igl::cat(1, ET, zero, tempMat2);
        igl::cat(0, tempMat1, tempMat2, S);

        solver_.compute(S);
        if(solver_.info()!= Eigen::Success) {
            std::cout << "decomposition failed" << std::endl;
        }
    }

    void Optimize::recursive_refine_normals(easy3d::SurfaceMesh::Face f) {
        auto tar_nomrals = mesh_->get_face_property<easy3d::vec3>("f:tar_normals");
        auto src_points = mesh_->get_face_property<int>("f:src_idx");
        auto vis = mesh_->get_face_property<bool>("f:vis");

        vis[f] = true;

        for(auto h : mesh_->halfedges(f)) {
            auto nxt_h = mesh_->next(h);
            auto f1 = mesh_->face(mesh_->opposite(h));
            auto f2 = mesh_->face(mesh_->opposite(nxt_h));
            if(f1.is_valid() && f2.is_valid()) {
                auto v1 = mesh_->position(mesh_->target(h)) - mesh_->position(easy3d::SurfaceMesh::Vertex(src_points[f]));
                auto v2 = mesh_->position(mesh_->source(h)) - mesh_->position(easy3d::SurfaceMesh::Vertex(src_points[f1]));
                auto v3 = mesh_->position(mesh_->target(nxt_h)) - mesh_->position(easy3d::SurfaceMesh::Vertex(src_points[f2]));
                if(easy3d::dot(v1, v2) < 0 && easy3d::dot(v1, v3) < 0 && easy3d::dot(v2, v3) > 0) {
                    tar_nomrals[f] = (tar_nomrals[f1] + tar_nomrals[f2]).normalize();
                    src_points[f] = src_points[f1];
                }
            }
        }
    }

    void Optimize::refine_tar_normals() {
        auto tar_nomrals = mesh_->get_face_property<easy3d::vec3>("f:tar_normals");
        auto src_points = mesh_->get_face_property<int>("f:src_idx");
        auto vis = mesh_->face_property<bool>("f:vis", false);
        bool state = false;
        int num = 0;
        int ct = 1;
//        do {
//            num = 0;
//            for (auto cur_f: mesh_->faces()) {
//                if (vis[cur_f] == state) {
//                    std::queue<easy3d::SurfaceMesh::Face> que;
//                    que.push(cur_f);
//                    vis[cur_f] = (!state);
//                    while (!que.empty()) {
//                        auto f = que.front();
//                        que.pop();
//                        if (vis[f] != state) continue;
//                        vis[f] = (!state);
//                        for (auto h: mesh_->halfedges(f)) {
//                            auto nxt_h = mesh_->next(h);
//                            auto f1 = mesh_->face(mesh_->opposite(h));
//                            auto f2 = mesh_->face(mesh_->opposite(nxt_h));
//                            if (f1.is_valid() && f2.is_valid()) {
//                                auto v1 = tar_nomrals[f];
//                                auto v2 = tar_nomrals[f1];
//                                auto v3 = tar_nomrals[f2];
//                                if (easy3d::dot(v1, v2) < 0.9 && easy3d::dot(v1, v3) < 0.9 && easy3d::dot(v2, v3) > 0.9) {
////                            tar_nomrals[f] = ((tar_nomrals[f1] + tar_nomrals[f2]) / 2).normalize();
//                                    tar_nomrals[f] = tar_nomrals[f1];
//                                    src_points[f] = src_points[f1];
//                                    num++;
//                                }
//                            }
//                            if (f1.is_valid()) {
//                                que.push(f1);
//                            }
//                        }
//                    }
//                }
//            }
//            ct++;
//            std::cout << "num=" << num <<std::endl;
//            state = (!state);
//        } while(num != 0);
////
//        while(true) {
//            num = 0;
//            for (auto f: mesh_->faces()) {
//                for (auto h: mesh_->halfedges(f)) {
//                    auto nxt_h = mesh_->next(h);
//                    auto f1 = mesh_->face(mesh_->opposite(h));
//                    auto f2 = mesh_->face(mesh_->opposite(nxt_h));
//                    if (f1.is_valid() && f2.is_valid()) {
//                        auto v1 = tar_nomrals[f];
//                        auto v2 = tar_nomrals[f1];
//                        auto v3 = tar_nomrals[f2];
//                        if (easy3d::dot(v1, v2) < 0.9 && easy3d::dot(v1, v3) < 0.9 && easy3d::dot(v2, v3) > 0.9) {
////                            tar_nomrals[f] = ((tar_nomrals[f1] + tar_nomrals[f2]) / 2).normalize();
//                            tar_nomrals[f] = tar_nomrals[f1];
//                            src_points[f] = src_points[f1];
//                            num++;
//                        }
//                    }
//                }
//            }
//            std::cout << "asd_num=" << num << std::endl;
//            if(num == 0) break;
//        }

        auto vis_1 = mesh_->vertex_property<bool>("f:vis1", false);
        state = true;
        do {
            num = 0;
            for (auto cur_v: mesh_->vertices()) {
                if(vis_1[cur_v] == state) continue;
                std::queue<easy3d::SurfaceMesh::Vertex>que;
                que.push(cur_v);
                while(!que.empty()) {
                    auto v = que.front(); que.pop();
                    if(vis_1[v] == state) continue;
                    vis_1[v] = state;
                    auto st_h = mesh_->out_halfedge(v);
                    auto it = st_h;
                    do {
                        auto cur_f = mesh_->face(it);
                        auto prev_f = mesh_->face(mesh_->prev_around_source(it));
                        auto nxt_f = mesh_->face(mesh_->next_around_source(it));
                        if (cur_f.is_valid() && prev_f.is_valid() && nxt_f.is_valid()) {
                            auto v1 = tar_nomrals[cur_f];
                            auto v2 = tar_nomrals[prev_f];
                            auto v3 = tar_nomrals[nxt_f];
                            if (easy3d::dot(v1, v2) < 0.9 && easy3d::dot(v1, v3) < 0.9 && easy3d::dot(v2, v3) > 0.9) {
//                            tar_nomrals[f] = ((tar_nomrals[f1] + tar_nomrals[f2]) / 2).normalize();
                                tar_nomrals[cur_f] = tar_nomrals[prev_f];
                                src_points[cur_f] = src_points[prev_f];
                                num++;
                            }
                        }
                        auto tar_v = mesh_->target(it);
                        if(vis_1[tar_v] != state) {
                            que.push(tar_v);
                        }
                        it = mesh_->next_around_source(it);
                    } while (it != st_h);
                }
            }
            state = (!state);
            std::cout << "num=" << num <<std::endl;
        } while(num != 0);

        for(int i = 0; i < 20; i++) {
            for (auto cur_v: mesh_->vertices()) {
                if(vis_1[cur_v] == state) continue;
                std::queue<easy3d::SurfaceMesh::Vertex>que;
                que.push(cur_v);
                while(!que.empty()) {
                    auto v = que.front(); que.pop();
                    if(vis_1[v] == state) continue;
                    vis_1[v] = state;
                    auto st_h = mesh_->out_halfedge(v);
                    auto it = st_h;
                    do {
                        auto cur_f = mesh_->face(it);
                        auto prev_f = mesh_->face(mesh_->prev_around_source(it));
                        auto nxt_f = mesh_->face(mesh_->next_around_source(it));
                        if (cur_f.is_valid() && prev_f.is_valid() && nxt_f.is_valid()) {
                            auto v1 = tar_nomrals[cur_f];
                            auto v2 = tar_nomrals[prev_f];
                            auto v3 = tar_nomrals[nxt_f];
                            if (easy3d::dot(v1, v2) > 0.9 && easy3d::dot(v1, v3) > 0.9 && easy3d::dot(v2, v3) > 0.9) {
//                            tar_nomrals[f] = ((tar_nomrals[f1] + tar_nomrals[f2]) / 2).normalize();
                                tar_nomrals[cur_f] = (( tar_nomrals[cur_f] + tar_nomrals[prev_f] + tar_nomrals[nxt_f]) / 3).normalize();
//                                src_points[cur_f] = src_points[prev_f];
                                num++;
                            }
                        }
                        auto tar_v = mesh_->target(it);
                        if(vis_1[tar_v] != state) {
                            que.push(tar_v);
                        }
                        it = mesh_->next_around_source(it);
                    } while (it != st_h);
                }
            }
            state = (!state);
        }
        


    }

    void Optimize::edge_init(const std::vector<size_t>& face_ancestors,
                             const std::vector<Eigen::Vector3d>& face_tar_normals,
                             const std::vector<size_t>& fixed_points) {
        int nb_points = mesh_->n_vertices();

    }

    void Optimize::point_init(const std::vector<Eigen::Vector3d>& normals,
                             const std::vector<size_t>& fixed_points) {

    }

    void Optimize::face_init(const std::vector<Eigen::Vector3d>& normals,
                              const std::vector<size_t>& fixed_points) {

    }

    bool Optimize::solve() {

        auto& points = mesh_->points();
        int nb_points = points.size();
        Eigen::VectorXd p(nb_points * 3);
        for(int i = 0; i < nb_points; i++) {
            p[i] = points[i].x;
            p[i + nb_points] = points[i].y;
            p[i + 2 * nb_points] = points[i].z;
        }
        Eigen::VectorXd b(nb_points * 3 + d_.size());
        b << beta_ * p, d_;
        Eigen::VectorXd x = solver_.solve(b);
        if(solver_.info()!= Eigen::Success) {
            // solving failed
            std::cout << "solving failed" << std::endl;
            return false;
        }

        for(auto v : mesh_->vertices()) {
            int id = v.idx();
            mesh_->position(v) = easy3d::vec3(x[id],
                                              x[id + nb_points], x[id + 2 * nb_points]);
        }

    }

    void Optimize::get_points(std::vector<Eigen::Vector3d>& points) {
        int nb_points = mesh_->n_vertices();
        points.resize(nb_points);

        for(auto v : mesh_->vertices()) {
            auto p = mesh_->position(v);
            points[v.idx()] = Eigen::Vector3d(p.x, p.y, p.z);
        }
    }

    void Optimize::get_faces(std::vector<std::vector<size_t>> &faces) {
        for(auto f : mesh_->faces()) {
            std::vector<size_t> tmp;
            for(auto v : mesh_->vertices(f)) {
                tmp.emplace_back(v.idx());
            }
            faces.emplace_back(tmp);
        }
        auto tar = mesh_->get_face_property<easy3d::vec3>("f:tar_normals");
        easy3d::Viewer viewer("vector_field");
        viewer.add_model(mesh_);
        const Box3 &box = mesh_->bounding_box();
        float length = norm(box.max_point() - box.min_point()) * 0.01f;
        std::vector<vec3> tmp;
        for (auto f : mesh_->faces()) {
            easy3d::vec3 center = easy3d::vec3(0,0,0);
            int ct =0;
            for(auto v : mesh_->vertices(f)) {
                center += mesh_->position(v); ct++;
            }
            center /= ct;
            const easy3d::vec3 s = center;
            easy3d::vec3 v = tar[f];
            const vec3 t = center + v * length;
            tmp.push_back(s);
            tmp.push_back(t);
        }
        auto drawable = mesh_->renderer()->add_lines_drawable("normals");
        drawable->update_vertex_buffer(tmp);
        drawable->set_uniform_coloring(vec4(0.0f, 1.0f, 0.0f, 1.0f));
        drawable->set_line_width(3.0f);
        viewer.run();
    }

    void Optimize::remesh() {
        bool is_update = false;
        auto tar_nomrals = mesh_->get_face_property<easy3d::vec3>("f:tar_normals");
        auto src_points = mesh_->get_face_property<int>("f:src_idx");


        for(auto f : mesh_->faces()) {
            std::vector<easy3d::SurfaceMesh::Halfedge> hal;
            std::vector<double> len;
            double p = 0;
            for(auto h : mesh_->halfedges(f)) {
                hal.emplace_back(h);
                len.emplace_back(mesh_->edge_length(mesh_->edge(h)));
                p += len.back();
            }
            p /= 2;
            double s = std::sqrt(p * (p - len[0]) * (p - len[1]) * (p - len[2]));
            if(s < 1e-2) {
                is_update = true;
                easy3d::vec3 center = easy3d::vec3(0, 0, 0);
                int ct = 0;
                if(len[0] > len[1] && len[0] > len[2]) {
                    auto ver = mesh_->target(mesh_->next(hal[0]));
                    for(auto v : mesh_->vertices(ver)) {
                        center += mesh_->position(v);ct++;
                    }
                    mesh_->position(ver) = center / ct;
                    auto ff = mesh_->face(mesh_->opposite(hal[1]));
                    tar_nomrals[f] = tar_nomrals[ff];
                    src_points[f] = src_points[ff];
                }
                else if(len[1] > len[0] && len[1] > len[2]) {
                    auto ver = mesh_->target(mesh_->next(hal[1]));
                    for(auto v : mesh_->vertices(ver)) {
                        center += mesh_->position(v);ct++;
                    }
                    mesh_->position(ver) = center / ct;
                    auto ff = mesh_->face(mesh_->opposite(hal[2]));
                    tar_nomrals[f] = tar_nomrals[ff];
                    src_points[f] = src_points[ff];
                } else {
                    auto ver = mesh_->target(mesh_->next(hal[0]));
                    for(auto v : mesh_->vertices(ver)) {
                        center += mesh_->position(v);ct++;
                    }
                    mesh_->position(ver) = center / ct;
                    auto ff = mesh_->face(mesh_->opposite(hal[0]));
                    tar_nomrals[f] = tar_nomrals[ff];
                    src_points[f] = src_points[ff];
                }
            }
        }
        if(is_update) {
            std::cout << "SAD" <<std::endl;
            update_solver();
        }
    }
}