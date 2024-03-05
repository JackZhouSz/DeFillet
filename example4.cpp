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
    std::vector<CGAL_Point> cgal_sites(nb_points);
    for(auto v : mesh->vertices()) {
        easy3d::vec3 pos = mesh->position(v);
        cgal_sites[v.idx()] = CGAL_Point(pos.x, pos.y, pos.z);
    }
    std::vector<easy3d::vec3> sites;
    UTILS::cgal_points_to_easy3d_points(cgal_sites, sites);
    easy3d::Box3 box = mesh->bounding_box();
    double diag = box.diagonal_length();

    Voronoi vor(cgal_sites);
    vor.cal();
    const std::vector<CGAL_Point>& cgal_vertices = vor.get_finite_vertices();

    std::vector<easy3d::vec3> vertices;
    UTILS::cgal_points_to_easy3d_points(cgal_vertices, vertices);

    const std::vector<std::vector<int>>& regions = vor.get_finite_regions();


    std::vector<std::vector<int>> poles(nb_points);
    for(int i = 0; i < nb_points; i++) {
        int nb_region = regions[i].size();
        int flag = 2;
        easy3d::vec3 vv;
        for(int j = 0; j < nb_region; j++) {

            if(!box.contains(vertices[regions[i][j]])) {
                vv = (vertices[regions[i][j]] - sites[i]);
                flag = 1;
                break;
            }
        }
        if(flag == 1) {
            double maxx = 0;
            int id = -1;
            for (int j = 0; j < nb_region; j++) {
                easy3d::vec3 tmp = vertices[regions[i][j]] - sites[i];
                if(easy3d::dot(tmp, vv) < 0) {
                    double len = tmp.norm();
                    if(maxx > len) {
                        maxx = len;
                        id = regions[i][j];
                    }
                }
            }
            poles[i].emplace_back(id);
        }
        else if(flag == 2) {
            double maxx = 0;
            int id1 = -1;
            easy3d::vec3 v1;
            for (int j = 0; j < nb_region; j++) {
                easy3d::vec3 tmp = vertices[regions[i][j]] - sites[i];
                double len = tmp.norm();
                if(len > maxx) {
                    maxx = len;
                    id1 = regions[i][j];
                    v1 = tmp.normalize();
                }
            }
            maxx = 0;
            int id2 = -1;
            for (int j = 0; j < nb_region; j++) {
                easy3d::vec3 tmp = vertices[regions[i][j]] - sites[i];
                double len = tmp.norm();
                if(len > maxx && easy3d::dot(tmp, v1) < 0) {
                    maxx = len;
                    id2 = regions[i][j];
                }
            }
            poles[i].emplace_back(id1);
            poles[i].emplace_back(id2);
        }
    }

    easy3d::PointCloud* cloud = new PointCloud;
    for(int i = 0; i < nb_points; i++) {
        for(int j = 0; j < poles[i].size(); j++) {
            cloud->add_vertex(vertices[poles[i][j]]);
        }
    }
    easy3d::Viewer viewer("ASD");
    viewer.add_model(cloud);
    viewer.run();
    return 0;
}
