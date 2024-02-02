//
// Created by xiaowuga on 12/9/2023.
//

#include <iostream>
#include <ostream>
#include <fstream>
#include "voronoi.h"
#include "common.h"
#include "defillet.h"
#include "visualization.h"
#include <chrono>
#include <Eigen/IterativeLinearSolvers>

#include <nlohmann/json.hpp>





using json = nlohmann::json;


using namespace std;


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


    std::vector<double> vertices_density_field;
    DEFILLET::computes_vertices_density_field(selected_vertices,1.0, vertices_density_field);

    std::vector<double> site_density_field;
    std::vector<double> site_pole_radius_field;
    std::vector<int> site_to_vertives;
    DEFILLET::compute_sites_pole_radius_field1(sites, selected_vertices,
                                              selected_regions, vertices_density_field,
                                              site_density_field, site_pole_radius_field,
                                              site_to_vertives);

    int nb_site = sites.size();
    std::vector<double> density_over_radius(nb_site);
    for(int i = 0; i < nb_site; i++) {
        density_over_radius[i] = site_density_field[i] / site_pole_radius_field[i];
    }
    double_vector1D_data_normalize(density_over_radius);
//    points_field_visualization(sites, density_over_radius);
//    return 0;

    // 打开文件流

    std::vector<std::vector<std::pair<int,double>>> mesh_graph;
    extract_vertex_graph_from_mesh(sites, f_ind, mesh_graph);
    double thr = 0.25;
    std::vector<int> island;
    std::vector<int> asd;
    do {

        DEFILLET::find_island(mesh_graph, density_over_radius, island, thr, 150);
        std::cout << "island size = " << island.size() << std::endl;
        asd.insert(asd.end(), island.begin(), island.end());
        DEFILLET::corrosion_island(sites, island, 1.0, thr, density_over_radius);
    } while(island.size() > 0);
    points_field_visualization(sites, density_over_radius);
    return 0;
    std::cout << asd.size() <<std::endl;
    std::vector<easy3d::vec3> pp;
    for(int i = 0; i < sites.size(); i++) {
        pp.emplace_back(sites[i].x(), sites[i].y(), sites[i].z());
    }
    for(int i = 0; i < selected_vertices.size(); i++) {
        pp.emplace_back(selected_vertices[i].x(), selected_vertices[i].y(), selected_vertices[i].z());
    }
    easy3d::Viewer viewer("scalar_field");
    easy3d::PointsDrawable* pd = new easy3d::PointsDrawable("pd");
    pd->update_vertex_buffer(pp);
    viewer.add_drawable(pd);

    easy3d::LinesDrawable* ld = new easy3d::LinesDrawable("ld");
    std::vector<easy3d::vec3> b;
    for(int i = 0; i < asd.size(); i++) {
        int id = asd[i];
        b.emplace_back(pp[id]);
        int id2 = site_to_vertives[id];
        b.emplace_back(pp[id2 + sites.size()]);

    }
    ld->update_vertex_buffer(b);
    viewer.add_drawable(ld);
    viewer.run();
    return 0;
}