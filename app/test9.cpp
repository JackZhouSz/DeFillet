//
// Created by xiaowuga on 2024/12/2.
//
#include "fillet_seg_v7.h"

#include <omp.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <nlohmann/json.hpp>
#include <CLI/CLI.hpp>
using json = nlohmann::json;

int main(int argc, char **argv) {
    easy3d::logging::initialize(false, false,
             false,
            false);
    std::string config_path;
    CLI::App app{"DEFILLET Command Line"};
    app.add_option("-c,--config", config_path, "Configure file")->required();
    CLI11_PARSE(app, argc, argv);
    std::ifstream file(config_path);
    json j;
    file >> j;
    std::string path = j["path"];
    float eps = j["eps"];
    int num_samples = j["num_samples"];
    float radius_thr = j["radius_thr"];
    float angle_thr = j["angle_thr"];
    bool local_voronoi = j["local_voronoi"];
    bool sor_filter = j["sor_filter"];
    bool radius_filter = j["radius_filter"];
    int num_sor_iter = j["num_sor_iter"];
    int num_neighbors = j["num_neighbors"];
    float sor_std_radio = j["sor_std_radio"];
    float h = j["h"];
    int num_drift_iter = j["num_drift_iter"];
    float mean_shift_eps = j["mean_shift_eps"];
    float dis_thr = j["dis_thr"];
    int min_support = j["min_support"];
    int max_blending_angle = j["max_blending_angle"];
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(path);
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    printf("Number of threads: %d\n", num_threads);
    std::cout << path << std::endl;
    FilletSegV7 fillet_seg_v7(mesh, eps, num_samples
              , radius_thr, angle_thr, local_voronoi, sor_filter, radius_filter
              , num_sor_iter, num_neighbors, sor_std_radio, h);
    fillet_seg_v7.run();
    return 0;
}