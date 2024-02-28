//
// Created by 小乌嘎 on 2024/2/10.
//
#include "common.h"
#include "voronoi.h"
#include "io.h"
#include "utils.h"
#include "visualization.h"
#include "defillet.h"

int main() {
    std::string file_path = "../data/bottle_fillet_remeshed.ply";
    std::vector<std::array<double, 3>> array_points;
    std::vector<std::vector<size_t>> array_faces;
    IO::read_ply_mesh(file_path, array_points, array_faces);

    std::vector<CGAL_Point> cgal_points;
    UTILS::array_points_to_cgal_points(array_points, cgal_points);
    Voronoi vor(cgal_points);

    vor.cal();
    const std::vector<CGAL_Point>& cgal_finite_vertices = vor.get_finite_vertices();
    const std::vector<std::vector<int>>& cgal_finite_regions = vor.get_finite_regions();

    std::vector<easy3d::vec3> easy3d_points;
//    UTILS::cgal_points_to_easy3d_points(cgal_finite_vertices, easy3d_points);

    const std::vector<std::vector<int>>& finite_cell_pole = vor.get_finite_cell_pole();

    std::vector<Eigen::Vector3d> eigen_points;
    UTILS::array_points_to_eigen_points(array_points, eigen_points);
    std::vector<Eigen::Vector3d> eigen_finite_vertices;
    UTILS::cgal_points_to_eigen_points(cgal_finite_vertices, eigen_finite_vertices);
    std::vector<double> density_field;
    std::vector<double> fillet_field;
    std::vector<size_t> fillet_poles;
    DEFILLET::compute_fillet_field(eigen_points, array_faces,
                                   eigen_finite_vertices, finite_cell_pole, density_field,
                                   fillet_field, fillet_poles, 1e-2, 1.0);

    UTILS::array_points_to_easy3d_points(array_points, easy3d_points);

//    int nb_points = eigen_points.size();
//    std::vector<double> pole_fields(eigen_points.size());
//    for(int i = 0; i < nb_points; i++) {
//        pole_fields[i] = (eigen_finite_vertices[fillet_poles[i]] - eigen_points[i]).norm();
//    }
//    RENDER::points_scalar_field_visualization(easy3d_points, pole_fields);
//    return 0;

//    RENDER::points_scalar_field_visualization(easy3d_points, fillet_field);
    std::vector<double> labels;
    DEFILLET::run_graph_cut(eigen_points, array_faces, eigen_finite_vertices,
                            fillet_poles, fillet_field,labels, 0.5, 0.24);

//    RENDER::points_scalar_field_visualization(easy3d_points, labels);
    std::vector<Eigen::Vector3d> fillet_points;
    std::vector<std::vector<size_t>> fillet_faces;
    std::vector<int> point_map;
    std::vector<std::vector<size_t>> non_fillet_faces;
    std::vector<size_t> fillet_bounding;
    std::vector<Eigen::Vector3d> fillet_bounding_normals;
    DEFILLET::extract_interest_region(eigen_points, array_faces, labels,
                                      fillet_points, fillet_faces, point_map, non_fillet_faces,
                                      fillet_bounding, fillet_bounding_normals);

    std::vector<easy3d::vec3> easy3d_fillet_points;

    UTILS::eigen_points_to_easy3d_points(fillet_points, easy3d_fillet_points);
    std::vector<size_t> point_ancestors;
    std::vector<size_t> face_ancestors;
    std::vector<double> distance;
    std::vector<Eigen::Vector3d> point_tar_normals;
    std::vector<Eigen::Vector3d> face_tar_normals;
    DEFILLET::run_geodesic(fillet_points, fillet_faces, fillet_bounding,
                           fillet_bounding_normals,
                           point_ancestors, face_ancestors, distance, point_tar_normals, face_tar_normals);



    std::vector<Eigen::Vector3d> new_fillet_points;
    std::vector<std::vector<size_t>> new_fillet_faces;
    if(DEFILLET::iterative_optimize(fillet_points, fillet_faces, point_ancestors, face_ancestors,  point_tar_normals,
                                    face_tar_normals, fillet_bounding, new_fillet_points,new_fillet_faces, "edge-based", 1.0, 20)) {
        for(int i = 0; i < new_fillet_points.size(); i++) {
            int id = point_map[i];
            eigen_points[id] = new_fillet_points[i];
        }
        for(int i = 0; i < new_fillet_faces.size(); i++) {
            std::vector<size_t> tmp;
            int num = new_fillet_faces[i].size();
            for(int j = 0; j < num; j++) {
                tmp.emplace_back(point_map[new_fillet_faces[i][j]]);
            }
            non_fillet_faces.emplace_back(tmp);
        }
        UTILS::eigen_points_to_easy3d_points(eigen_points, easy3d_points);
        RENDER::mesh_visualization(easy3d_points, non_fillet_faces);
    }

//    UTILS::eigen_points_to_easy3d_points(eigen_points, easy3d_points);
//    RENDER::mesh_visualization(easy3d_points, array_faces);
    return 0;
}