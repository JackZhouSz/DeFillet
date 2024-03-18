//
// Created by 13900K on 2024/3/18.
//

#include "defillet.h"
#include "surafce_mesh_segmenter.h"

#include <Xin_Wang.h>

void DeFillet::run_geodesic() {
    extract_fillet_region();
    sources_.clear();
    for(auto v : fillet_mesh_->vertices()) {
        if(mesh_->is_border(v)) {
            sources_.insert(v.idx());
        }
    }
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
    CXin_Wang alg(xin_mesh, sources_);
    alg.Execute();

    auto normals = fillet_mesh_->face_property<easy3d::vec3>("f:normal");

    for(auto f : fillet_mesh_->faces()) {
        normals[f] = fillet_mesh_->compute_face_normal(f);
    }

    auto face_sources = fillet_mesh_->face_property<int>("f:sources");
    auto face_tar_normals = fillet_mesh_->face_property<easy3d::vec3>("f:tar_normals");
    for(int i = nb_points; i < nb_points + nb_faces; i++) {
        easy3d::SurfaceMesh::Face f(i - nb_points);
        face_sources[f] = alg.GetAncestor(i);
        mesh_->compute_vertex_normal()
        face_tar_normals[f] = normals[easy3d::SurfaceMesh::Face(face_sources[f])]
    }
}

void DeFillet::refine_target_normal() {

}

void DeFillet::run_defillet() {

}

void DeFillet::optimize(){

}

void DeFillet::extract_fillet_region() {
    auto gcp = mesh_->face_property<int>("f:gcp_labels");
    easy3d::SurfaceMeshSegmenter seg(mesh_);
    if(fillet_mesh_)
        delete fillet_mesh_;
    fillet_mesh_ = seg.segment(gcp, 1);

}