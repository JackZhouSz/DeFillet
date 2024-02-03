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


    std::vector<double> vertices_density_field;
    DEFILLET::computes_vertices_density_field(selected_vertices,1.0, vertices_density_field);
    double_vector1D_data_normalize(vertices_density_field);
//    points_field_visualization(selected_vertices, vertices_density_field);
//    return 0;

    std::vector<double> site_density_field;
    std::vector<double> site_pole_radius_field;
    std::vector<int> site_to_vertives;
    std::vector<std::vector<int>> das;
    DEFILLET::compute_sites_pole_radius_field1(sites, selected_vertices,
                                              selected_regions, vertices_density_field,
                                              site_density_field, site_pole_radius_field,
                                              site_to_vertives,das);

//    int nb_site = sites.size();
//    std::vector<double> density_over_radius(nb_site);
//    for(int i = 0; i < nb_site; i++) {
//        density_over_radius[i] = site_density_field[i] / site_pole_radius_field[i];
//    }
    double_vector1D_data_normalize(vertices_density_field);
    std::vector<double>sss(selected_vertices.size(), 0);
    for(int i = 0; i < das.size(); i++) {
        double len = 0;
        for(int j = 0; j < das[i].size(); j++) {
            len += std::sqrt(CGAL::squared_distance(sites[das[i][j]], selected_vertices[i]));
        }
        if(das[i].size() > 0) {
            len /= das[i].size();
        }
        sss[i] = len;

    }
    double_vector1D_data_normalize(sss);
    points_field_visualization(selected_vertices, sss);
    return 0;


    return 0;
}