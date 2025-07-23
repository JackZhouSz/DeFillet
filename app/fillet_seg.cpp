//
// Created by xiaowuga on 2024/12/2.
//
#include "fillet_seg_v8.h"
#include "defillet_v2.h"
#include "util.h"
#include <omp.h>
#include <voro_viewer.h>

#include <nlohmann/json.hpp>
#include <CLI/CLI.hpp>

#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/util/file_system.h>
using json = nlohmann::json;

void save_fillet_seg_result(FilletSegV8& fillet_seg_v8, std::string out_dir, bool flag = true);
void save_defillet_result(DeFilletv2& de_filletv2, std::string out_dir, bool flag = true);

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
    std::string out_dir = j["out_dir"];
    float eps = j["eps"];
    int num_samples = j["num_samples"];
    float radius_thr = j["radius_thr"];
    float angle_thr = j["angle_thr"];
    bool local_voronoi = j["local_voronoi"];
    bool radius_filter = j["radius_filter"];
    int num_neighbors = j["num_neighbors"];
    float h = j["h"];
    float beta = j["beta"];
    float gamma = j["gamma"];
    int num_smmoth_iter = j["num_smmoth_iter"];
    float lamdba = j["lamdba"];
    int num_opt_iter = j["num_opt_iter"];
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(path);
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    printf("Number of threads: %d\n", num_threads);
    FilletSegV8 fillet_seg_v8(mesh, eps, num_samples
              , radius_thr, angle_thr, local_voronoi, radius_filter
              , num_neighbors,  h,  num_smmoth_iter, lamdba);
    fillet_seg_v8.run();
    // return 0;
    save_fillet_seg_result(fillet_seg_v8, out_dir, true);
    easy3d::VoroViewer viewer("ASD");
    auto mavv = fillet_seg_v8.mavv_;
    auto mavvns = fillet_seg_v8.mavvns_;
    auto mavvcs = fillet_seg_v8.mavvcs_;
    viewer.init(mavv, mesh, mavvns, mavvcs);
    viewer.run();
    // save_fillet_seg_result()
    DeFilletv2 de_filletv2(mesh, angle_thr, beta, gamma, num_opt_iter);
    de_filletv2.run();
    return 0;
}


void save_fillet_seg_result(FilletSegV8& fillet_seg_v8, std::string out_dir, bool flag) {
    easy3d::SurfaceMesh* mesh = fillet_seg_v8.mesh_;
    std::string base_name = easy3d::file_system::base_name(mesh->name());
    if(flag) {
        auto vv = fillet_seg_v8.vv_;
        auto vvr = fillet_seg_v8.vvr_;
        save_point_field(vv, vvr,out_dir + base_name + "_vvr.ply");
        auto vv_density = fillet_seg_v8.vv_density_;
        save_point_field(vv, vv_density, out_dir + base_name + "_vvd.ply");
        auto mavv = fillet_seg_v8.mavv_;
        auto mavvr = fillet_seg_v8.mavvr_;
        save_point_field(mavv, mavvr,out_dir + base_name + "_mavvr.ply");
        auto ss = fillet_seg_v8.ss_;
        save_mesh_field(mesh, ss, out_dir + base_name + "_ss.ply");
        save_mesh_field_with_mtl(mesh, ss, out_dir + base_name + "_ss.obj",out_dir + base_name + "_ss.mtl");
    }
    int num = mesh->n_faces();
    std::vector<float> fillet(mesh->n_faces(), 0.0);
    auto labels = fillet_seg_v8.fillet_labels_;
    for(int i = 0; i < num; i++) {
        fillet[i] = labels[i];
    }
    save_mesh_field_with_mtl(mesh, fillet, out_dir + base_name + "_fillet.obj", out_dir + base_name + "_fillet.mtl");


}


void save_defillet_result(DeFilletv2& de_filletv2, std::string out_dir, bool flag) {

}