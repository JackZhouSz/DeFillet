//
// Created by xiaowuga on 2025/8/14.
//


#include <utils.h>
#include <decompose.h>

#include <easy3d/fileio/surface_mesh_io.h>


int main() {
    std::string case_path = "../data/091_wheel_assembly_part.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);
    DeFillet::Decompose decompose(mesh, 60);

    decompose.apply();

    std::vector<easy3d::SurfaceMesh*> components = decompose.all_components();

    DeFillet::save_components(mesh, components, "../out3/components.ply");

}
