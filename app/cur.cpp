//
// Created by xiaowuga on 2025/1/16.
//
#include <easy3d/algo/surface_mesh_curvature.h>

#include <easy3d/fileio/surface_mesh_io.h>

#include "util.h"

int main() {
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load("D:\\Sig25_DeFillet\\cur\\high_res_1r.ply");
    easy3d::SurfaceMeshCurvature smc(mesh);
    smc.analyze();
    auto cur_max = mesh->get_vertex_property<float>("v:curv-max");
    std::vector<float> ff = cur_max.vector();
    for(int i = 0; i < ff.size(); i++) {
        ff[i] = abs(ff[i]);
    }
    save_mesh_field_with_mtl1(mesh,ff, "../out/cur.obj", "../out/cur.mtl");

    return 0;
}