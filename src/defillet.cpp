//
// Created by 13900K on 2024/3/18.
//

#include "defillet.h"


void DeFillet::run_geodesic() {

}

void DeFillet::refine_target_normal() {

}

void DeFillet::run_defillet() {

}

void DeFillet::optimize(){

}

void DeFillet::extract_interest_region() {
    std::set<easy3d::SurfaceMesh::Face> faces_set;
    std::set<easy3d::SurfaceMesh::Vertex> vertex_set;
    auto gcp = mesh_->face_property<int>("f:gcp_labels");
    for(auto f : mesh_->faces()) {
        if(gcp[f]) {
            faces_set.insert(f);
            for(auto v : mesh_->vertices(f)) {
                vertex_set.insert(v);
            }
        }
    }
}