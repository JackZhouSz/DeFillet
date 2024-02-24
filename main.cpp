//
// Created by xiaowuga on 12/9/2023.
//

#include <iostream>
#include <ostream>
#include <fstream>
#include "voronoi.h"
#include "common.h"
#include "defillet1.h"
#include "visualization.h"
#include "utils.h"
#include <chrono>
#include <Eigen/IterativeLinearSolvers>

using namespace std;

double threshold_otsu(std::vector<double>& data, int nb_bin) {
    std::vector<double> dh = data;
    double_vector1D_data_normalize(dh);
    double thresh = -1;

    int nb_data = data.size();
    std::vector<double> bins(nb_bin + 1, 0);
    double u = 0;

    for(int i = 0; i < nb_data; i++) {
        int id = int(dh[i] * nb_bin);
        bins[id] += 1;
    }
    double p = 1.0 / nb_data;
    for(int i = 0 ; i < nb_bin + 1; i++) {
        bins[i] *= p;
        u += i * bins[i];
    }
    std::cout << "ASD" <<std::endl;
    double wk = 0, uk = 0, gmax = std::numeric_limits<double>::min();
//    std::cout << "ASD" <<std::endl;
    std::cout << nb_bin <<std::endl;
    for(int i = 0; i < nb_bin + 1; i++) {
//        std::cout << i <<std::endl;
        wk += bins[i];
        uk += i * bins[i];
        double tmp = uk - u * wk;
        double sigma = tmp * tmp / (wk * (1.0 - wk) + 1e-6);
        if(sigma >= gmax) {
            thresh = i;
            gmax = sigma;
        }
    }

    std::cout << thresh <<std::endl;
    return thresh / nb_bin;
}

int main() {
    std::string file_path = "../data/bottle_fillet_remeshed.ply";
    std::vector<std::array<double, 3>> v_pos;
    std::vector<std::vector<size_t>> f_ind;
    read_ply_mesh(file_path,v_pos, f_ind);

    std::vector<Point> sites;
    my_points_convert_to_cgal_points(v_pos, sites);
    Voronoi vor(sites);
    std::vector<Point> vertices;
    std::vector<std::vector<int>> sites_regions;
    vor.cal(vertices, sites_regions);

    std::vector<Point> selected_vertices;
    std::vector<std::vector<int>>  selected_regions;

    select_vertices_and_regions_inside_sites_bounding_box(sites, vertices, sites_regions,
                                                          selected_vertices, selected_regions);

    points_visualization(selected_vertices);
    return 0;
    std::vector<double> vertices_density_field;
    DEFILLET::computes_vertices_density_field(selected_vertices,1.0, vertices_density_field);
//    double_vector1D_data_normalize(vertices_density_field);
//    points_field_visualization(selected_vertices, vertices_density_field);
//    return 0;
    std::vector<double> site_density_field;
    std::vector<double> site_pole_radius_field;
    DEFILLET::compute_sites_pole_radius_field(sites, selected_vertices,
                                              selected_regions, vertices_density_field,
                                              site_density_field, site_pole_radius_field);

    int nb_site = sites.size();
    std::vector<double> density_over_radius(nb_site);
    for(int i = 0; i < nb_site; i++) {
        density_over_radius[i] = site_density_field[i] / site_pole_radius_field[i];
    }
    double_vector1D_data_normalize(density_over_radius);
    double thr = threshold_otsu(density_over_radius, 255);

//    std::cout << thrs << std::endl;
//    return 0;
//

//    std::ofstream outputFile("../output.txt");
//
//    if (outputFile.is_open()) {
//        // 将向量元素写入文件
//        for (const auto &element : density_over_radius) {
//            outputFile << element << "\n";
//        }
//
//        // 关闭文件流
//        outputFile.close();
//        std::cout << "Vector elements written to output.txt." << std::endl;
//    } else {
//        std::cerr << "Unable to open output.txt for writing." << std::endl;
//    }
//    return 0;

    std::vector<std::vector<std::pair<int,double>>> mesh_graph;
    extract_vertex_graph_from_mesh(sites, f_ind, mesh_graph);
    std::vector<int> island;
    do {

        DEFILLET::find_island(mesh_graph, density_over_radius, island, thr, 150);
        std::cout << "island size = " << island.size() << std::endl;
        DEFILLET::corrosion_island(sites, island, 1.0, thr, density_over_radius);
    } while(island.size() > 0);
    
//    points_field_visualization(sites, density_over_radius);
//    return 0;

    std::vector<std::pair<int,int>> edge;
    std::vector<double> edge_weight;
    extract_edge_from_mesh(sites, f_ind, edge, edge_weight);
    std::vector<int> labels;
    DEFILLET::run_graph_cut(density_over_radius, edge, edge_weight, labels, 0.8,thr);
    std::vector<double> ss(labels.size());
    for(int i = 0; i < labels.size(); i++) {
        ss[i] = labels[i];
    }
//    points_field_visualization(sites,ss);
//    return 0;



    std::vector<bool> is_fillet;
    DEFILLET::classify_defillet_triangle_via_vertices_voting(f_ind, density_over_radius, thr, is_fillet);
    binary_mesh_segmentation_visualization(sites, f_ind, is_fillet);
//    return 0;
    std::vector<Point> fillet_points;
    std::vector<std::vector<size_t>> fillet_faces;
    std::vector<int> point_map_table;
    std::vector<int> faces_map_table;
    std::vector<bool> face_labels;
    DEFILLET::extract_interest_region(sites, f_ind, is_fillet, fillet_points,
                                      fillet_faces, point_map_table, faces_map_table,face_labels);

//    std::cout << fillet_faces.size() <<std::endl;

//    binary_mesh_segmentation_visualization(fillet_points, fillet_faces, face_labels);
    std::vector<Point> dense_points;
    std::vector<std::vector<size_t>> dense_faces;
    std::vector<int> centroid;
    std::vector<bool> centroid_label;
    densify_mesh_via_centroid_insertion(fillet_points, fillet_faces,
                                        face_labels, dense_points, dense_faces,
                                        centroid, centroid_label);
    std::vector<int> ancestor;
    std::vector<double> distance;
    DEFILLET::run_geodesic(dense_points, dense_faces, centroid,
                           centroid_label, ancestor, distance);

//    std::set<int> tmp;
//
//    for(size_t i = 0; i < ancestor.size(); i++) {
//        int xx = ancestor[i] - fillet_points.size();
//        tmp.insert(xx);
//    }
//    std::vector<bool> bin(fillet_faces.size(), false);
//
//    for(auto i : tmp) {
//        bin[i] = true;
//    }
//    binary_mesh_segmentation_visualization(fillet_points, fillet_faces, bin);
//    return 0;
    std::vector<Vector_3> target_normals;
    std::vector<int> fixed_points;
    std::vector<std::pair<int,std::pair<int,Vector_3>>> edge_vector;

    DEFILLET::geodesic_post_processing(fillet_points, fillet_faces, ancestor, face_labels,target_normals,
                                       fixed_points, edge_vector, 10.0);

    double_vector1D_data_normalize(distance);
    mesh_face_normals_vector_field(fillet_points, fillet_faces, target_normals);
//    return 0;
    std::vector<Point> new_fillet_points;
    auto start_time = std::chrono::high_resolution_clock::now();
    if(DEFILLET::iterative_optimize(fillet_points, fillet_faces, target_normals, new_fillet_points, fixed_points, 2, 10)) {

        std::vector<std::array<double, 3>> my_points;
        cgal_points_convert_to_my_points(new_fillet_points, my_points);
//        write_ply_points("../data/all3.ply", my_points);
        for(int i = 0; i < new_fillet_points.size(); i++) {
            int id = point_map_table[i];
            sites[id] = new_fillet_points[i];
        }

    } else {
        return 0;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    mesh_visualization(sites, f_ind);
//    if(DEFILLET::optimize_with_fixed_solve(fillet_points, fillet_faces, target_normals, new_fillet_points, fixed_points, edge_vector)) {
//        std::vector<std::array<double, 3>> my_points;
//        cgal_points_convert_to_my_points(new_fillet_points, my_points);
//        write_ply_points("../data/all.ply", my_points);
////        mesh_visualization(new_fillet_points, fillet_faces);
//    }


    return 0;
}