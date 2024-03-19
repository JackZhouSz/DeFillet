//
// Created by 13900K on 2024/3/18.
//

#include "defillet.h"
#include "surafce_mesh_segmenter.h"

#include <Xin_Wang.h>

void DeFillet::run_geodesic() {
    extract_fillet_region();
    sources_.clear();
    sources_normals_.clear();
    auto original_point_index = fillet_mesh_->vertex_property<int>("v:original_index");
    std::map<int, int> mp;
    for(auto v : fillet_mesh_->vertices()) {
        if(fillet_mesh_->is_border(v)) {
            sources_.emplace_back(v.idx());
            easy3d::SurfaceMesh::Vertex origin_v(original_point_index[v]);
            mp[v.idx()] = sources_normals_.size();
            sources_normals_.emplace_back(mesh_->compute_vertex_normal(origin_v));

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
//    std::cout << "ASD" <<std::endl;
//    refine_target_normal();
}

void DeFillet::refine_target_normal() {
    auto face_sources = fillet_mesh_->face_property<int>("f:sources");
    auto face_tar_normals = fillet_mesh_->face_property<easy3d::vec3>("f:tar_normals");
    int nb_points = fillet_mesh_->n_vertices();
    std::vector<bool>vis(nb_points);
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
                auto st_h = mesh_->out_halfedge(v);
                auto it = st_h;
                do {
                    auto cur_f = mesh_->face(it);
                    auto prev_f = mesh_->face(mesh_->prev_around_source(it));
                    auto nxt_f = mesh_->face(mesh_->next_around_source(it));
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
                    auto tar_v = mesh_->target(it);
                    if(vis[tar_v.idx()] != s) {
                        que.push(tar_v);
                    }
                    it = mesh_->next_around_source(it);
                } while (it != st_h);
            }
        }
        s = (!s);
        std::cout << "num=" << num <<std::endl;
    } while(num != 0);
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

    fillet_mesh_ = seg.segment<int>(gcp, 1);
}