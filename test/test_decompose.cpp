//
// Created by xiaowuga on 2025/8/14.
//


#include <decompose.h>
#include <easy3d/fileio/surface_mesh_io.h>


int main() {
    std::string case_path = "../data/091_wheel_assembly_part.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);
    DeFillet::Decompose decompose(mesh, 60);

    decompose.apply();

    std::vector<easy3d::SurfaceMesh*> components = decompose.all_components();


    for(size_t i = 0; i < components.size(); i++) {
        easy3d::SurfaceMeshIO::save("../out3/" + std::to_string(i) + ".ply", components[i]);
    }

}
