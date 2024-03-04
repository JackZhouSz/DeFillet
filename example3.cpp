//
// Created by 13900K on 2024/3/4.
//

#include <iostream>
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>

using namespace std;

int main() {
    std::string file_path = "../data/20440_27177360_1.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    auto face_normal = mesh->add_face_property<easy3d::vec3>("f:normal");
    for(auto f : mesh->faces()) {
        face_normal[f] = mesh->compute_face_normal(f);
    }

//    
    return 0;
}