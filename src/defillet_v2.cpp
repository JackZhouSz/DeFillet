//
// Created by xiaowuga on 2024/12/15.
//

#include "util.h"
#include "defillet_v2.h"
#include "surafce_mesh_segmenter.h"



#include <Xin_Wang.h>

#include <easy3d/util/stop_watch.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/algo/surface_mesh_geometry.h>

#include <igl/cat.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD>

DeFilletv2::DeFilletv2(easy3d::SurfaceMesh* mesh
                       , double angle_thr
                       , double beta
                       , double gamma
                       , int num_opt_iter)
                        : mesh_(mesh)
                        , angle_thr_(angle_thr)
                        , beta_(beta)
                        , gamma_(gamma)
                        , num_opt_iter_(num_opt_iter) {

}

void DeFilletv2::run() {
    initialize();
    optimize();
}


void DeFilletv2::initialize() {
    std::cout << "Start defillet initialization..." <<std::endl;
    easy3d::StopWatch sw; sw.start();
    auto fillet_labels = mesh_->face_property<int>("f:fillet_labels");

    auto labels = mesh_->face_property<int>("f:labels");
    for(auto f : mesh_->faces()) {
        labels[f] = 1;
    }
    for(auto e : mesh_->edges()) {
        auto f0 = mesh_->face(e, 0);
        auto f1 = mesh_->face(e, 1);
        if(fillet_labels[f0] + fillet_labels[f1] < 2) {
            labels[f0] = 0; labels[f1] = 0;
            // std::cout << "ASD" <<std::endl;
        }
        // if(fillet_labels[f0] == 0)
        // if(fillet_labels[f0] + fillet_labels[f1] < 2) {
        //     labels[f0] = 0; labels[f1] = 0;
        // }
    }
    easy3d::SurfaceMeshSegmenter sms(mesh_);
    focus_area_ =  sms.segment<int>(labels, 0);
    std::cout << focus_area_->n_vertices() << std::endl;
    // easy3d::SurfaceMesh* non_focus_area_ = sms.segment<int>(labels, 1);
    // easy3d::SurfaceMeshIO::save("../out/non_fillet.ply", non_focus_area_);
    easy3d::SurfaceMeshIO::save("../out/focus_area.ply", focus_area_);
    auto ofi = focus_area_->face_property<int>("f:original_index"); //original face index
    std::map<int, easy3d::vec3> fixed;
    for(auto v : focus_area_->vertices()) {
        int ct = 0, flag = 0;
        for(auto h : focus_area_->halfedges(v)) {
            auto oh = focus_area_->opposite(h);
            auto f0 = focus_area_->face(h);
            auto f1 = focus_area_->face(oh);
            if(f0.is_valid()) {
                int id = ofi[f0];
                auto ff = easy3d::SurfaceMesh::Face(id);
                if(mesh_->is_border(ff)) {
                    flag =1; break;
                }
            }
            if(f1.is_valid() && mesh_->is_border(f1) == true) {
                int id = ofi[f1];
                auto ff = easy3d::SurfaceMesh::Face(id);
                if(mesh_->is_border(ff)) {
                    flag =1; break;
                }
            }
            if(!f0.is_valid() || !f1.is_valid()) {
                ct++;
            }
            else if(f0.is_valid() && f1.is_valid()) {
                easy3d::vec3 n0 = focus_area_->compute_face_normal(f0);
                easy3d::vec3 n1 = focus_area_->compute_face_normal(f1);
                double angle = dihedral_angle(n0, n1);
                if(angle > angle_thr_) {
                    flag = 1;
                }
            }
        }
        if(ct == 2 && flag == 0) {

            fixed[v.idx()] = focus_area_->compute_vertex_normal(v);
        }
    }
    int nb_points = focus_area_->n_vertices();
    int nb_faces = focus_area_->n_faces();

    std::vector<CPoint3D> xin_points;
    for(auto v : focus_area_->vertices()) {
        auto& p = focus_area_->position(v);
        xin_points.emplace_back(CPoint3D(p.x, p.y, p.z));
    }
    std::vector<CBaseModel::CFace> xin_faces;

    for(auto f : focus_area_->faces()) {
        int num = 0;
        easy3d::vec3 center = easy3d::vec3(0,0,0);
        std::vector<size_t> indices;
        for(auto v : focus_area_->vertices(f)) {
            center += focus_area_->position(v);
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
    std::set<int> special;
    std::set<int> xin_sources;

    for(auto e : focus_area_->edges()) {
        auto f0 = focus_area_->face(e, 0);
        auto f1 = focus_area_->face(e, 1);
        if(f0.is_valid() && f1.is_valid()) {
            easy3d::vec3 n0 = focus_area_->compute_face_normal(f0);
            easy3d::vec3 n1 = focus_area_->compute_face_normal(f1);
            double angle = dihedral_angle(n0, n1);
            if(angle > angle_thr_) {
                auto f = easy3d::SurfaceMesh::Face(ofi[f0]);
                if(fillet_labels[f] == 1) {
                    f = f0;
                }
                else {
                    f = f1;
                }
                for(auto v : focus_area_->vertices(f)) {
                    if(fixed.find(v.idx()) != fixed.end()) {
                        special.insert(v.idx());
                    }
                }
            }
        }
    }
    for(auto id : fixed) {
        if(special.find(id.first) == special.end()) {
            xin_sources.insert(id.first);
        }
    }
    auto fr = focus_area_->face_property<int>("f:face_root"); // face root
    auto ftn = focus_area_->face_property<easy3d::vec3>("f:tar_normals"); // face target normal
    auto pgd = focus_area_->vertex_property<float>("v:geo_dis");
    /*std::vector<easy3d::vec3> ff;
    for(auto id : xin_sources) {
        ff.emplace_back(focus_area_->points()[id]);
    }

    save_point_set(ff, "../out/ff.ply");*/
    CXin_Wang alg(xin_mesh, xin_sources);
    alg.Execute();
    std::cout << "AS11D" <<std::endl;
    for(int i = 0; i < nb_points; i++) {
        easy3d::SurfaceMesh::Vertex v(i);
        pgd[v] = alg.GetDistanceField()[i];
    }
    std::vector<float> geo;
    for(int i = 0; i < nb_faces; i++) {
        easy3d::SurfaceMesh::Face f(i );
        fr[f] = alg.GetAncestor(i + nb_points);
        ftn[f] = fixed[fr[f]];

    }
    std::vector<easy3d::vec3> sp;
    for(auto id : special) {
        auto v = easy3d::SurfaceMesh::Vertex(id);
        pgd[v] = 0.0;
        sp.emplace_back(focus_area_->position(v));
        easy3d::vec3 n = focus_area_->compute_vertex_normal(v);
        for(auto f : focus_area_->faces(v)) {
            fr[f] = v.idx();
            ftn[f] = n;
        }
    }
    save_point_set(sp, "../out/sp.ply");
    std::vector<easy3d::vec3> s;
    for(auto id : xin_sources) {
        auto v = easy3d::SurfaceMesh::Vertex(id);
        if(special.find(id) == special.end()) {
            s.emplace_back(focus_area_->position(v));
        }
    }
    save_point_set(s, "../out/s.ply");
    std::vector<bool>vis(nb_points, false);
    bool state = true;
    int num;
    do {
        num = 0;
        for(auto v : focus_area_->vertices()) {
            if(vis[v.idx()] == state) continue;
            std::queue<easy3d::SurfaceMesh::Vertex>que;
            que.push(v);
            while(!que.empty()) {
                auto cur = que.front(); que.pop();
                auto st_h = focus_area_->out_halfedge(cur);
                auto it = st_h;
                do {
                    auto cur_f = focus_area_->face(it);
                    auto prev_f = focus_area_->face(focus_area_->prev_around_source(it));
                    auto nxt_f = focus_area_->face(focus_area_->next_around_source(it));
                    if (cur_f.is_valid() && prev_f.is_valid() && nxt_f.is_valid()) {
                        auto n1 = ftn[cur_f];
                        auto n2 = ftn[prev_f];
                        auto n3 = ftn[nxt_f];
                        if (dihedral_angle(n1, n2) > angle_thr_ && dihedral_angle(n1, n3) > angle_thr_ && dihedral_angle(n2, n3) < angle_thr_) {
                            ftn[cur_f] = ftn[prev_f];
                            fr[cur_f] = fr[prev_f];
                            num++;
                        }
                    }
                    auto tar_v = focus_area_->target(it);
                    if(vis[tar_v.idx()] != state) {
                        que.push(tar_v);
                        vis[tar_v.idx()] = state;
                    }
                    it = focus_area_->next_around_source(it);
                } while (it != st_h);
            }
        }
        std::cout << "number of faca normals flip = " << num <<std::endl;
    }while(num != 0);
    init_time_ = sw.elapsed_seconds(5);
    std::cout << "Defillet initialize sussessfully! time="<< init_time_ <<std::endl;
    // easy3d::SurfaceMeshIO::save("../out/focus_area.obj", focus_area_);
    save_mesh_field_with_mtl1(focus_area_, pgd.vector(), "../out/focus.obj", "../out/focus.mtl");
    save_mesh_normals(focus_area_, ftn.vector(), "../out/ftn.obj");

}


void DeFilletv2::optimize() {
    e_e.clear(), e_f.clear();
    easy3d::StopWatch sw; sw.start();
    int nb_points = focus_area_->n_vertices();
    int nb_faces = focus_area_->n_faces();
    std::vector<float> nx(nb_faces);
    std::vector<float> ny(nb_faces);
    std::vector<float> nz(nb_faces);
    auto fr = focus_area_->face_property<int>("f:face_root");
    auto ftn = focus_area_->face_property<easy3d::vec3>("f:tar_normals");
#pragma omp parallel for
    for(int i = 0; i < nb_faces; i++) {
        auto f = easy3d::SurfaceMesh::Face(i);
        nx[i] = ftn[f].x;
        ny[i] = ftn[f].y;
        nz[i] = ftn[f].z;
    }
    double beta1 = 1.0;
    std::vector<Eigen::Triplet<double>> triplets;
    int row = 0;
    for(auto e : focus_area_->edges()) {
        int v0 = focus_area_->vertex(e, 0).idx();
        int v1 = focus_area_->vertex(e, 1).idx();
        auto f0 = focus_area_->face(e, 0);
        auto f1 = focus_area_->face(e, 1);

        if(f0.is_valid() && f1.is_valid()) {
            triplets.emplace_back(Eigen::Triplet<double>(row, v0, beta1 * nx[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, beta1 *ny[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, beta1 *nz[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, beta1 *(-nx[f0.idx()])));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, beta1 *(-ny[f0.idx()])));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, beta1 *(-nz[f0.idx()])));

            row++;

            triplets.emplace_back(Eigen::Triplet<double>(row, v0, beta1 * nx[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, beta1 * ny[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, beta1 * nz[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, beta1 * (-nx[f1.idx()])));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, beta1 * (-ny[f1.idx()])));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, beta1 * (-nz[f1.idx()])));

            row++;
        }
    }
    Eigen::SparseMatrix<double> E(row, nb_points * 3); // Edge Directional Constraint Energy
    E.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    std::set<int> fixed;
    for(auto f : focus_area_->faces()) {
        int num = 0;
        for(auto v : focus_area_->vertices(f)) {
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx(), beta_ * nx[f.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + nb_points, beta_ * ny[f.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + 2 * nb_points, beta_ * nz[f.idx()]));
            num++;
        }

        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f], -nx[f.idx()] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f] + nb_points, -ny[f.idx()] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f] + 2 * nb_points, -nz[f.idx()] * beta_ * num));
        fixed.insert(fr[f]);
    }
    Eigen::SparseMatrix<double> F(nb_faces, nb_points * 3); // Face Normal Constraint Energy
    F.setFromTriplets(triplets.begin(), triplets.end());

    row = 0;
    d_.resize(nb_points * 3);
    triplets.clear();
    std::vector<easy3d::vec3> ff;
    for(auto id : fixed) {
        easy3d::SurfaceMesh::Vertex v(id);
        // ff.emplace_back(focus_area_->position(v));
        auto pos = focus_area_->position(v);
        triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() , 1.0));
        d_[row++] = pos.x;
        triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() + nb_points, 1.0));
        d_[row++] = pos.y;
        triplets.emplace_back(Eigen::Triplet<double>(row,  v.idx() + nb_points * 2, 1.0));
        d_[row++] = pos.z;
    }
    // save_point_set(ff, "../out/asd.ply");
    d_.conservativeResize(row);

    Eigen::SparseMatrix<double> D(row, nb_points * 3);
    D.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> A;
    igl::cat(1, E, F, A);
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> ATA = AT * A;
    Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
    I = gamma_ * I;

    Eigen::SparseMatrix<double> Q = I + ATA;
    Eigen::SparseMatrix<double> zero(row, row);
    Eigen::SparseMatrix<double> DT = D.transpose();
    Eigen::SparseMatrix<double> tempMat1;
    Eigen::SparseMatrix<double> tempMat2;
    Eigen::SparseMatrix<double> M;
    igl::cat(1, Q, D, tempMat1);
    igl::cat(1, DT, zero, tempMat2);
    igl::cat(2, tempMat1, tempMat2, M);
    solver_.compute(M);
    if(solver_.info()!= Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
        return;
    }


    for(int iter = 0; iter < num_opt_iter_; iter++) {
        computeEe_and_E_f();
        std::cout << "iter " << iter + 1 << " is processing." << std::endl;
        auto& points = focus_area_->points();

        Eigen::VectorXd p(nb_points * 3);
        for(int i = 0; i < nb_points; i++) {
            p[i] = points[i].x;
            p[i + nb_points] = points[i].y;
            p[i + 2 * nb_points] = points[i].z;
        }

        Eigen::VectorXd b(nb_points * 3 + d_.size());
        b << gamma_ * p, d_;
        Eigen::VectorXd x = solver_.solve(b);
        if(solver_.info()!= Eigen::Success) {
            // solving failed
            std::cout << "solving failed" << std::endl;
            return;
        }
        double ec = 0;
        for(auto v : focus_area_->vertices()) {

            int id = v.idx();
            easy3d::vec3 new_p = easy3d::vec3(x[id],
                                              x[id + nb_points], x[id + 2 * nb_points]);
            ec += (new_p - focus_area_->position(v)).length2();

            focus_area_->position(v) = new_p;
        }
        e_c.emplace_back(ec);

    }

    opt_time_ = sw.elapsed_seconds(5);
    std::cout << "Defillet optimize sussessfully! time="<< opt_time_ <<std::endl;
    easy3d::SurfaceMeshIO::save("../out/result.ply", focus_area_);

    auto opi = focus_area_->vertex_property<int>("v:original_index");
    for(auto v : focus_area_->vertices()) {
        int id = opi[v];
        easy3d::SurfaceMesh::Vertex vv(id);
        mesh_->position(vv) = focus_area_->position(v);
    }
    auto gcp = mesh_->face_property<int>("f:fillet_labels");
    std::vector<float> label;
    for(auto f : mesh_->faces()) {
        label.emplace_back(gcp[f]);
    }
    save_mesh_field(mesh_, label, "../out/a.ply");
    // easy3d::SurfaceMeshIO::save("../out/a.ply", mesh_);
    // save_mesh_field_with_mtl(mesh_, label, "../out/defillet.obj", "../out/defillet.mtl");

    for(int i = 0; i < num_opt_iter_; i++) {
        std::cout << e_e[i] << ", ";
    }
    std::cout << std::endl;
    for(int i = 0; i < num_opt_iter_; i++) {
        std::cout << e_f[i] << ", ";
    }
    std::cout << std::endl;
    for(int i = 0; i < num_opt_iter_; i++) {
        std::cout << e_c[i] << ", ";
    }
    std::cout << std::endl;

}

Eigen::VectorXd solve_ldlt(const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > &ldlt,
    const Eigen::VectorXd& right,
    double pinvtoler)
{
    Eigen::VectorXd X_0 = ldlt.permutationP() * right;
    Eigen::VectorXd X_1 = ldlt.matrixL().solve(X_0);
    Eigen::VectorXd X_2(ldlt.vectorD().size());
    X_2.setZero();
    for (int i = 0; i < ldlt.vectorD().size(); ++i)
        if (abs(ldlt.vectorD()(i)) > pinvtoler)
            X_2[i] = X_1[i] / ldlt.vectorD()(i);
    Eigen::VectorXd X_3 = ldlt.matrixU().solve(X_2);
    return ldlt.permutationPinv() * X_3;
}

void DeFilletv2::optimize2() {
    easy3d::StopWatch sw; sw.start();
    int nb_points = focus_area_->n_vertices();
    int nb_faces = focus_area_->n_faces();
    std::vector<float> nx(nb_faces);
    std::vector<float> ny(nb_faces);
    std::vector<float> nz(nb_faces);
    auto fr = focus_area_->face_property<int>("f:face_root");
    auto ftn = focus_area_->face_property<easy3d::vec3>("f:tar_normals");
#pragma omp parallel for
    for(int i = 0; i < nb_faces; i++) {
        auto f = easy3d::SurfaceMesh::Face(i);
        nx[i] = ftn[f].x;
        ny[i] = ftn[f].y;
        nz[i] = ftn[f].z;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    int row = 0;
    for(auto e : focus_area_->edges()) {
        int v0 = focus_area_->vertex(e, 0).idx();
        int v1 = focus_area_->vertex(e, 1).idx();
        auto f0 = focus_area_->face(e, 0);
        auto f1 = focus_area_->face(e, 1);

        if(f0.is_valid() && f1.is_valid()) {
            triplets.emplace_back(Eigen::Triplet<double>(row, v0, nx[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, ny[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, nz[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, -nx[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, -ny[f0.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, -nz[f0.idx()]));

            row++;

            triplets.emplace_back(Eigen::Triplet<double>(row, v0, nx[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + nb_points, ny[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v0 + 2 * nb_points, nz[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1, -nx[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + nb_points, -ny[f1.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(row, v1 + 2 * nb_points, -nz[f1.idx()]));

            row++;
        }
    }
    Eigen::SparseMatrix<double> E(row, nb_points * 3); // Edge Directional Constraint Energy
    E.setFromTriplets(triplets.begin(), triplets.end());

    triplets.clear();
    std::set<int> fixed;
    for(auto f : focus_area_->faces()) {
        int num = 0;
        for(auto v : focus_area_->vertices(f)) {
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx(), beta_ * nx[f.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + nb_points, beta_ * ny[f.idx()]));
            triplets.emplace_back(Eigen::Triplet<double>(f.idx(), v.idx() + 2 * nb_points, beta_ * nz[f.idx()]));
            num++;
        }

        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f], -nx[f.idx()] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f] + nb_points, -ny[f.idx()] * beta_ * num));
        triplets.emplace_back(Eigen::Triplet<double>(f.idx(), fr[f] + 2 * nb_points, -nz[f.idx()] * beta_ * num));
        fixed.insert(fr[f]);
    }
    Eigen::SparseMatrix<double> F(nb_faces, nb_points * 3); // Face Normal Constraint Energy
    F.setFromTriplets(triplets.begin(), triplets.end());

    row = 0;
    d_.resize(nb_points * 3);
    triplets.clear();
    std::vector<easy3d::vec3> ff;
    for(auto id : fixed) {
        easy3d::SurfaceMesh::Vertex v(id);
        // ff.emplace_back(focus_area_->position(v));
        auto pos = focus_area_->position(v);
        triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() , 1.0));
        d_[row++] = pos.x;
        triplets.emplace_back(Eigen::Triplet<double>(row, v.idx() + nb_points, 1.0));
        d_[row++] = pos.y;
        triplets.emplace_back(Eigen::Triplet<double>(row,  v.idx() + nb_points * 2, 1.0));
        d_[row++] = pos.z;
    }
    // save_point_set(ff, "../out/asd.ply");
    d_.conservativeResize(row);

    Eigen::SparseMatrix<double> D(row, nb_points * 3);
    D.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> A;
    igl::cat(1, E, F, A);
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> ATA = AT * A;
    Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
    I = gamma_ * I;

    Eigen::SparseMatrix<double> Q = ATA;
    Eigen::SparseMatrix<double> zero(row, row);
    Eigen::SparseMatrix<double> DT = D.transpose();
    Eigen::SparseMatrix<double> tempMat1;
    Eigen::SparseMatrix<double> tempMat2;
    Eigen::SparseMatrix<double> M;
    igl::cat(1, Q, D, tempMat1);
    igl::cat(1, DT, zero, tempMat2);
    igl::cat(2, tempMat1, tempMat2, M);
    // bool is_symmetric = Q.isApprox(Q.transpose());
    // if (!is_symmetric) {
    //     std::cerr << "Matrix is not symmetric!" << std::endl;
    //     // return -1;
    // }
    // solver_.compute(M);
    // if(solver_.info()!= Eigen::Success) {
    //     std::cout << "decomposition failed" << std::endl;
    //     return;
    // }



    auto& points = focus_area_->points();

    Eigen::VectorXd p(nb_points * 3);
    for(int i = 0; i < nb_points; i++) {
        p[i] = points[i].x;
        p[i + nb_points] = points[i].y;
        p[i + 2 * nb_points] = points[i].z;
    }
    std::cout << "ASD" <<std::endl;

    Eigen::VectorXd b(nb_points * 3 + d_.size());
    b << 1.0 * p, d_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    Eigen::VectorXd x = solve_ldlt(solver, b, 1e-6);
     // = solver.solve(b);
    std::cout << "ASD" <<std::endl;

    for(auto v : focus_area_->vertices()) {
        int id = v.idx();
        focus_area_->position(v) = easy3d::vec3(x[id],
                                          x[id + nb_points], x[id + 2 * nb_points]);
    }
    // for(int iter = 0; iter < num_opt_iter_; iter++) {
    //     std::cout << "iter " << iter + 1 << " is processing." << std::endl;
    //     auto& points = focus_area_->points();
    //
    //     Eigen::VectorXd p(nb_points * 3);
    //     for(int i = 0; i < nb_points; i++) {
    //         p[i] = points[i].x;
    //         p[i + nb_points] = points[i].y;
    //         p[i + 2 * nb_points] = points[i].z;
    //     }
    //
    //     Eigen::VectorXd b(nb_points * 3 + d_.size());
    //     b << gamma_ * p, d_;
    //     Eigen::VectorXd x = solver_.solve(b);
    //     if(solver_.info()!= Eigen::Success) {
    //         // solving failed
    //         std::cout << "solving failed" << std::endl;
    //         return;
    //     }
    //
    //     for(auto v : focus_area_->vertices()) {
    //         int id = v.idx();
    //         focus_area_->position(v) = easy3d::vec3(x[id],
    //                                           x[id + nb_points], x[id + 2 * nb_points]);
    //     }
    // }

    opt_time_ = sw.elapsed_seconds(5);
    std::cout << "Defillet optimize sussessfully! time="<< opt_time_ <<std::endl;
    easy3d::SurfaceMeshIO::save("../out/result.ply", focus_area_);

}



void DeFilletv2::computeEe_and_E_f() {

    auto ftn = focus_area_->face_property<easy3d::vec3>("f:tar_normals"); // face target normal
    double ee = 0, ef = 0;
    for(auto e : focus_area_->edges()) {
        auto f1 = focus_area_->face(e, 0);
        auto f2 = focus_area_->face(e, 1);
        auto v1 = focus_area_->position(focus_area_->vertex(e, 0));
        auto v2 = focus_area_->position(focus_area_->vertex(e, 1));
        easy3d::vec3 tmp = (v1 - v2);
        if(f1.is_valid() && f2.is_valid()) {
            ee += fabs(easy3d::dot( ftn[f1], tmp)) + fabs(easy3d::dot( ftn[f2], tmp));
            // easy3d::vec3 p1 = cross(ftn[f1], ftn[f2]).normalize();
            //
            // ee += cross(tmp, p1).length();
        }
    }
    e_e.emplace_back(ee);

    for(auto f : focus_area_->faces()) {
        auto n = focus_area_->compute_face_normal(f);
        // ef += 1.0 -  fabs(easy3d::dot(ftn[f], n));
        ef += easy3d::cross(ftn[f], n).length();
    }
    e_f.emplace_back(ef);
}



