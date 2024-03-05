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
    std::string file_path = "../data/20440_27177360_1_part_0.ply";
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
    const std::vector<CGAL_Point>& cgal_finite_vertices = vor.get_finite_vertices();
    const std::vector<std::vector<int>>& cgal_finite_regions = vor.get_finite_regions();
    const std::vector<std::vector<int>>& finite_cell_pole = vor.get_finite_cell_pole();

    int nb_sites = nb_points + nb_faces;
    std::vector<double> field(nb_faces);
    for(int i = nb_points; i < nb_sites; i++) {
        easy3d::SurfaceMesh::Face f(i - nb_points);
        easy3d::vec3 v1 = face_normal[f];
        double minn = 1.0;
        for(int j = 0; j < finite_cell_pole[i].size(); j++) {
            int id = finite_cell_pole[i][j];
            easy3d::vec3 s = easy3d::vec3( cgal_finite_vertices[id].x(),cgal_finite_vertices[id].y(),
                                           cgal_finite_vertices[id].z());
            easy3d::vec3 t = easy3d::vec3(sites[i].x(), sites[i].y(), sites[i].z());
            easy3d::vec3 v2 = (t - s).normalize();
            double val = 1.0 - std::fabs(easy3d::dot(v1, v2));
            if(minn > val) {
                minn = val;
            }
        }
        field[i - nb_points] = minn;
    }

    std::vector<easy3d::vec3> pp;
    for(int i = 0; i < nb_faces; i++) {
        int id = i + nb_points;
        easy3d::vec3 pos = easy3d::vec3(sites[id].x(),sites[id].y(), sites[id].z());
        pp.emplace_back(pos);
    }
    RENDER::points_scalar_field_visualization(pp, field);

    return 0;
}
