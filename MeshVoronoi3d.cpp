//
// Created by xiaowuga on 3/6/2024.
//

#include "MeshVoronoi3d.h"

MeshVoronoi3d::MeshVoronoi3d(easy3d::SurfaceMesh* mesh) : mesh_(mesh){

    int nb_sites = mesh_->n_faces();
    sites_.resize(nb_sites);
    auto face_normals = mesh_->face_property<easy3d::vec3>("f:normal");
    for(auto f : mesh_->faces()) {
        easy3d::vec3 center = easy3d::vec3(0,0,0);
        int ct = 0;
        for(auto v : mesh_->vertices(f)) {
            center += mesh_->position(v); ct++;
        }
        face_normals[f] = mesh_->compute_face_normal(f);
        sites_[f.idx()] = center / ct;
    }

    easy3d::Box3 box = mesh_->bounding_box();


    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T_.insert(CGAL_Point(sites_[i].x,
                                               sites_[i].y, sites_[i].z ));
        mp_sites[v] = i;
    }

    double diag = box.diagonal_length() * 10;
    double xmin = box.min_coord(0) - diag;
    double xmax = box.max_coord(0) + diag;
    double ymin = box.min_coord(1) - diag;
    double ymax = box.max_coord(1) + diag;
    double zmin = box.min_coord(2) - diag;
    double zmax = box.max_coord(2) + diag;
    T_.insert(CGAL_Point(xmin, ymin, zmax));
    T_.insert(CGAL_Point(xmax, ymin, zmax));
    T_.insert(CGAL_Point(xmin, ymax, zmax));
    T_.insert(CGAL_Point(xmax, ymax, zmax));
    T_.insert(CGAL_Point(xmin, ymin, zmin));
    T_.insert(CGAL_Point(xmax, ymin, zmin));
    T_.insert(CGAL_Point(xmin, ymax, zmin));
    T_.insert(CGAL_Point(xmax, ymax, zmin));

    vertices_.clear(); regions_.clear(); poles_.clear();
    final_pole_.clear(); vertices2sites_.clear();
    regions_.resize(nb_sites);
    poles_.resize(nb_sites);
    final_pole_.resize(nb_sites);

    for(auto cell = T_.finite_cells_begin(); cell != T_.finite_cells_end(); cell++) {
        auto v = T_.dual(cell);
        int vid = vertices_.size();
        vertices_.emplace_back(easy3d::vec3(v.x(), v.y(), v.z()));
        std::vector<int> tmp;
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            regions_[id].emplace_back(vid);
            tmp.emplace_back(id);
        }
        vertices2sites_.emplace_back(tmp);
    }


    for(int i = 0; i < nb_sites; i++) {
        int num = regions_[i].size();
        double maxx = 0;
        int pole = -1;
        easy3d::SurfaceMesh::Face f(i);
        for(int j = 0; j < num; j++) {
            easy3d::vec3 dir = vertices_[regions_[i][j]] - sites_[i];
            double len = dir.norm();
            if(box.contains(vertices_[regions_[i][j]]) &&easy3d::dot(face_normals[f], dir) / len < -0.5 && len > maxx) {
                maxx = len;
                pole = regions_[i][j];
            }
        }
        final_pole_[i] = pole;
    }





}