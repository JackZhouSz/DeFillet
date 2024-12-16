//
// Created by xiaowuga on 2024/12/15.
//

#include "io.h"
#include "defillet_v2.h"
#include "surafce_mesh_segmenter.h"

#include <easy3d/util/stop_watch.h>

#include <Xin_Wang.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/algo/surface_mesh_geometry.h>

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
}


void DeFilletv2::initialize() {
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
        }
    }
    easy3d::SurfaceMeshSegmenter sms(mesh_);
    focus_area_ =  sms.segment<int>(labels, 0);
    easy3d::SurfaceMeshIO::save("../out/focus_area.ply", focus_area_);
    auto normal = focus_area_->face_property<easy3d::vec3>("f:normals");
    auto ofi = focus_area_->face_property<int>("f:original_index"); //original face index

    std::map<int, easy3d::vec3> fixed;
    for(auto v : focus_area_->vertices()) {
        int ct = 0, flag = 0;
        for(auto h : focus_area_->halfedges(v)) {
            auto oh = focus_area_->opposite(h);
            auto f0 = focus_area_->face(h);
            auto f1 = focus_area_->face(oh);
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

    std::vector<easy3d::vec3> ff;
    for(auto id : xin_sources) {
        ff.emplace_back(focus_area_->points()[id]);
    }
    save_point_set(ff, "../out/ff.ply");
}

void DeFilletv2::optimize() {

}


double DeFilletv2::dihedral_angle(const easy3d::vec3& n1
                             , const easy3d::vec3& n2
                             , bool rad) {

    // double radians = abs(acos(dot(n1, n2)));
    double radians = easy3d::geom::angle(n1, n2);
    if(rad) {
        return radians;
    } else {
        double degrees = radians * 180.0 / M_PI;
        return degrees;
    }
                             }


