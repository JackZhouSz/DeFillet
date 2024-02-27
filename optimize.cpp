//
// Created by xiaowuga on 2/26/2024.
//

#include "optimize.h"

#include <igl/cat.h>

namespace DEFILLET {
    Optimize::Optimize(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       const std::vector<Eigen::Vector3d>& normals,
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

        if(type == "edge-based")
            edge_init(normals, fixed_points);
        else if(type == "point-based")
            point_init(normals, fixed_points);
        else if(type == "face-based")
            face_init(normals, fixed_points);
    }

    void Optimize::edge_init(const std::vector<Eigen::Vector3d>& normals,
                             const std::vector<size_t>& fixed_points) {
        int nb_points = mesh_->n_vertices();
        int nb_edges = mesh_->n_edges();
        int nb_faces = mesh_->n_faces();

        auto tar_nomrals = mesh_->add_face_property<easy3d::vec3>("f:tar_normals");
        for(auto f : mesh_->faces()) {
            int idx = f.idx();
            tar_nomrals[f] = easy3d::vec3(normals[idx].x(),
                                          normals[idx].y(), normals[idx].z());
        }
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
        Eigen::SparseMatrix<double> FDC(nb_faces * 3, nb_points * 3);
        for(auto f : mesh_->faces()) {
            std::vector<easy3d::vec3> pos;
            std::vector<int> idx;
            for(auto v : mesh_->vertices(f)) {
                pos.emplace_back(mesh_->position(v));
                idx.emplace_back(v.idx());
            }
            easy3d::vec3 tmp = pos[0] - pos[1];
            double e1 = (pos[0] - pos[1]).;
            double e2 = (pos[0] - pos[2]).norm();
            double e3 = (pos[1] - pos[2]).norm();
            double s1 = e1 / e2;
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[0], 1.0 - s1));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[0] + nb_points, 1.0 - s1));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[0] + 2 * nb_points, 1.0 - s1));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[1], -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[1] + nb_points, -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[1] + 2 * nb_points, -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[2], s1));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[2] + nb_points, s1));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx(), idx[2] + 2 * nb_points, s1));

            double s2 = e1 / e3;
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[0], 1.0 - s2));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[0] + nb_points, 1.0 - s2));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[0] + 2 * nb_points, 1.0 - s2));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[2], -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[2] + nb_points, -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[2] + 2 * nb_points, -1.0));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[1], s2));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[1] + nb_points, s2));
            triplets.emplace_back(Eigen::Triplet<double>(3 * f.idx() + 1, idx[1] + 2 * nb_points, s2));
        }

        int nb_fixed_points = fixed_points.size();
        d_ = Eigen::VectorXd(nb_fixed_points * 3);
        triplets.clear();
        Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
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

        Eigen::SparseMatrix<double> A = FNC;
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
}