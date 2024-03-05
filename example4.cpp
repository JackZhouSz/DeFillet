//
// Created by 13900K on 2024/3/5.
//
#include "common.h"
#include "voronoi.h"
#include "io.h"
#include "utils.h"
#include "visualization.h"
#include "defillet.h"

int main() {
    std::string file_path = "../data/bottle_fillet_remeshed.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    auto face_normal = mesh->face_property<easy3d::vec3>("f:normal");
    for(auto f : mesh->faces()) {
        face_normal[f] = mesh->compute_face_normal(f);
    }
    int nb_points = mesh->n_vertices();
    int nb_faces = mesh->n_faces();
    std::vector<CGAL_Point> sites(nb_points + nb_faces);
    for(auto v : mesh->vertices()) {
        easy3d::vec3 pos = mesh->position(v);
        sites[v.idx()] = CGAL_Point(pos.x, pos.y, pos.z);
    }
    for(auto f : mesh->faces()) {
        easy3d::vec3 pos = easy3d::vec3(0, 0, 0);
        int ct = 0;
        for(auto v : mesh->vertices(f)) {
            pos += mesh->position(v); ct++;
        }
        pos /= ct;
        sites[f.idx() + nb_points] = CGAL_Point(pos.x, pos.y, pos.z);
    }

    Voronoi vor(sites);
    vor.cal();
    const std::vector<CGAL_Point>& vertices = vor.get_vertices();
    const std::vector<std::vector<int>>&  vor.get_regions();
//    RENDER::points_scalar_field_visualization(pp, field);

    return 0;
}
