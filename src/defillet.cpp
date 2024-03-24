//
// Created by 13900K on 2024/3/18.
//

#include "defillet.h"
#include "surafce_mesh_segmenter.h"

#include <Xin_Wang.h>

#include <igl/cat.h>

#include <easy3d/util/stop_watch.h>


void DeFillet::run_geodesic() {
    easy3d::StopWatch sw; sw.start();
    extract_fillet_region();
    sources_.clear();
    sources_normals_.clear();
    auto original_point_index = fillet_mesh_->vertex_property<int>("v:original_index");
    auto fillet_sources = fillet_mesh_->vertex_property<int>("v:sources");
    std::map<int, int> mp;
    for(auto v : fillet_mesh_->vertices()) {
        easy3d::SurfaceMesh::Vertex vv(original_point_index[v]);
        if(fillet_sources[v] == 1) {
            sources_.emplace_back(v.idx());
            easy3d::SurfaceMesh::Vertex origin_v(original_point_index[v]);
            mp[v.idx()] = sources_normals_.size();
            sources_normals_.emplace_back(mesh_->compute_vertex_normal(origin_v));

        }
    }
    std::cout << sources_.size() << std::endl;
    int nb_points = fillet_mesh_->n_vertices();
    int nb_faces = fillet_mesh_->n_faces();
    std::vector<CPoint3D> xin_points;
    for(auto v : fillet_mesh_->vertices()) {
        auto& p = fillet_mesh_->position(v);
        xin_points.emplace_back(CPoint3D(p.x, p.y, p.z));
    }

    std::vector<CBaseModel::CFace> xin_faces;

    for(auto f : fillet_mesh_->faces()) {
        int num = 0;
        easy3d::vec3 center = easy3d::vec3(0,0,0);
        std::vector<size_t> indices;
        for(auto v : fillet_mesh_->vertices(f)) {
            center += fillet_mesh_->position(v);
            indices.emplace_back(v.idx());
            num++;
        }
        center /= num;
        xin_points.emplace_back(CPoint3D(center.x, center.y, center.z));

        for(int i = 0; i < num; i++) {
            xin_faces.emplace_back(CBaseModel::CFace(indices[i],
                                                     indices[(i + 1) % num], nb_points + f.idx()));
        }
    }
    CRichModel xin_mesh(xin_points, xin_faces);
    std::set<int> xin_sources(sources_.begin(), sources_.end());
    CXin_Wang alg(xin_mesh, xin_sources);
    alg.Execute();

    auto face_sources = fillet_mesh_->face_property<int>("f:sources");
    auto face_tar_normals = fillet_mesh_->face_property<easy3d::vec3>("f:tar_normals");
    auto point_geo_dis = fillet_mesh_->vertex_property<float>("v:geo_dis");
    for(int i = 0; i < nb_points; i++) {
        easy3d::SurfaceMesh::Vertex v(i);
        point_geo_dis[v] = alg.GetDistanceField()[i];
    }

    for(int i = 0; i < nb_faces; i++) {
        easy3d::SurfaceMesh::Face f(i );
        face_sources[f] = alg.GetAncestor(i + nb_points);
        face_tar_normals[f] = sources_normals_[mp[face_sources[f]]];
    }
    geodesic_time_ = sw.elapsed_seconds(3);
    refine_target_normal();
}

void DeFillet::refine_target_normal() {
    easy3d::StopWatch sw; sw.start();
    auto face_sources = fillet_mesh_->face_property<int>("f:sources");
    auto face_tar_normals = fillet_mesh_->face_property<easy3d::vec3>("f:tar_normals");
    int nb_points = fillet_mesh_->n_vertices();
    std::vector<bool>vis(nb_points, false);
    bool s = true;
    int num = 0;
    double thr = cos(angle_ * M_PI / 180.0);
    do {
        num = 0;
        for (auto cur_v: fillet_mesh_->vertices()) {

            if(vis[cur_v.idx()] == s) continue;
            std::queue<easy3d::SurfaceMesh::Vertex>que;
            que.push(cur_v);
            while(!que.empty()) {
                auto v = que.front(); que.pop();
                if(vis[v.idx()] == s) continue;
                vis[v.idx()] = s;
                auto st_h = fillet_mesh_->out_halfedge(v);
                auto it = st_h;
                do {
                    auto cur_f = fillet_mesh_->face(it);
                    auto prev_f = fillet_mesh_->face(fillet_mesh_->prev_around_source(it));
                    auto nxt_f = fillet_mesh_->face(fillet_mesh_->next_around_source(it));
                    if (cur_f.is_valid() && prev_f.is_valid() && nxt_f.is_valid()) {
                        auto v1 = face_tar_normals[cur_f];
                        auto v2 = face_tar_normals[prev_f];
                        auto v3 = face_tar_normals[nxt_f];
                        if (easy3d::dot(v1, v2) < thr && easy3d::dot(v1, v3) < thr && easy3d::dot(v2, v3) > thr) {
                            face_tar_normals[cur_f] = face_tar_normals[prev_f];
                            face_sources[cur_f] = face_sources[prev_f];
                            num++;
                        }
                    }
                    auto tar_v = fillet_mesh_->target(it);
                    if(vis[tar_v.idx()] != s) {
                        que.push(tar_v);
                    }
                    it = fillet_mesh_->next_around_source(it);
                } while (it != st_h);
            }
        }
        s = (!s);
        std::cout << "num=" << num <<std::endl;
    } while(num != 0);
    target_normals_refine_time_ = sw.elapsed_seconds(3);
}

bool DeFillet::run_defillet() {

    init_opt();

    easy3d::StopWatch sw;
    for(int i = 0; i < num_opt_iter_; i++) {
        std::cout << "iter " << i + 1 << " is processing." << std::endl;
        if(!opt()) {
            std::cout << "opt error." << std::endl;
            return false;
        }
    }
    defillet_iter_time_ = sw.elapsed_seconds(3);
    return true;
}

void DeFillet::init_opt(){
    easy3d::StopWatch sw; sw.start();
    int nb_points = fillet_mesh_->n_vertices();
    int nb_faces = fillet_mesh_->n_faces();
    auto nx = fillet_mesh_->face_property<float>("f:tar_normals_x");
    auto ny = fillet_mesh_->face_property<float>("f:tar_normals_y");
    auto nz = fillet_mesh_->face_property<float>("f:tar_normals_z");
    auto face_sources = fillet_mesh_->face_property<int>("f:sources");
    auto original_point_index = fillet_mesh_->vertex_property<int>("v:original_index");

    std::vector<Eigen::Triplet<double>> triplets;

    int row = 0;
    for(auto e : fillet_mesh_->edges()) {
        int v0 = fillet_mesh_->vertex(e, 0).idx();
        int v1 = fillet_mesh_->vertex(e, 1).idx();
        auto f0 = fillet_mesh_->face(e, 0);
        auto f1 = fillet_mesh_->face(e, 1);
        if(f0.is_valid() && f1.is_valid()) {
            triplets.emplace_back(Eigen::Triplet<double>(row, v0, nx[f0]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, ny[f0]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, nz[f0]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, -nx[f0]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, -ny[f0]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, -nz[f0]));

            row++;

            triplets.emplace_back(Eigen::Triplet<double>(row, v0, nx[f1]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, ny[f1]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, nz[f1]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, -nx[f1]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, -ny[f1]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, -nz[f1]));

            row++;
        }
    }

    Eigen::SparseMatrix<double> E1(row, nb_points * 3);
    E1.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    for(auto f : fillet_mesh_->faces()) {
        int num = 0;
        for(auto v : fillet_mesh_->vertices(f)) {
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx(), beta_ * nx[f]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + nb_points, beta_ * ny[f]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + 2 * nb_points, beta_ * nz[f]));
            num++;
        }

        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), face_sources[f], -nx[f] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), face_sources[f] + nb_points, -ny[f] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), face_sources[f] + 2 * nb_points, -nz[f] * beta_ * num));
    }
    Eigen::SparseMatrix<double> E2(nb_faces, nb_points * 3);
    E2.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<easy3d::SurfaceMesh::Vertex> fixed_points;
    auto fillet_sources = fillet_mesh_->vertex_property<int>("v:sources");
    row = 0;
    d_.resize(nb_points * 3);
    triplets.clear();
    for(auto v : fillet_mesh_->vertices()) {
        easy3d::SurfaceMesh::Vertex vv(original_point_index[v]);
        if(fillet_sources[v] == 1) {
            fixed_points.emplace_back(v);
            auto pos = fillet_mesh_->position(v);
            triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() , 1.0));
            d_[row++] = pos.x;
            triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() + nb_points, 1.0));
            d_[row++] = pos.y;
            triplets.emplace_back(Eigen::Triplet<double>(row,  v.idx() + nb_points * 2, 1.0));
            d_[row++] = pos.z;
        }
    }
    d_.conservativeResize(row);

    Eigen::SparseMatrix<double> FIX(row, nb_points * 3);
    FIX.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> A;
    igl::cat(1, E1, E2, A);
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> ATA = AT * A;
    Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
    I = gamma_ * I;

    Eigen::SparseMatrix<double> Q = I + ATA;
    Eigen::SparseMatrix<double> zero(row, row);
    Eigen::SparseMatrix<double> FIXT = FIX.transpose();
    Eigen::SparseMatrix<double> tempMat1;
    Eigen::SparseMatrix<double> tempMat2;
    Eigen::SparseMatrix<double> E;
    igl::cat(1, Q, FIX, tempMat1);
    igl::cat(1, FIXT, zero, tempMat2);
    igl::cat(2, tempMat1, tempMat2, E);

    solver_.compute(E);
    if(solver_.info()!= Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
    }
    defillet_init_time_ = sw.elapsed_seconds(3);
}

bool DeFillet::opt() {
    auto& points = fillet_mesh_->points();
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

    for(auto v : fillet_mesh_->vertices()) {
        int id = v.idx();
        fillet_mesh_->position(v) = easy3d::vec3(x[id],
                                          x[id + nb_points], x[id + 2 * nb_points]);
    }
    return true;
}

void DeFillet::extract_fillet_region() {
    auto gcp = mesh_->face_property<int>("f:gcp_labels");
    easy3d::SurfaceMeshSegmenter seg(mesh_);
    if(fillet_mesh_)
        delete fillet_mesh_;

    fillet_mesh_ = seg.segment<int>(gcp, 1);

    auto mesh_sources = mesh_->vertex_property<int>("v:sources");
    auto fillet_sources = fillet_mesh_->vertex_property<int>("v:sources");
    auto original_point_index = fillet_mesh_->vertex_property<int>("v:original_index");

    for(auto v : fillet_mesh_->vertices()) {
        easy3d::SurfaceMesh::Vertex vv(original_point_index[v]);
        fillet_sources[v] = mesh_sources[vv];
    }
}