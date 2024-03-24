//
// Created by 13900K on 2024/3/21.
//


#include <CLI/CLI.hpp>
#include <easy3d/core/surface_mesh.h>
#include <surafce_mesh_segmenter.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/util/file_system.h>
#include <queue>
#include <cmath>


int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string part_mesh;
        std::string output_dir;
    } args;
//    CLI::App app{"DEFILLET Command Line"};
//    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
//    app.add_option("-p,--part_mesh", args.part_mesh, "Part Mesh")->required();
//    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
//    CLI11_PARSE(app, argc, argv);
    args.input_mesh = "D:\\code\\defillet_case\\20501_ae5143b3_0\\compose.ply";
    args.part_mesh = "D:\\code\\defillet_case\\20501_ae5143b3_0\\2024-03-24_11-31-30_20501_ae5143b3_0_part_17\\defillet.ply";
    args.output_dir = "D:\\code\\defillet_case\\20501_ae5143b3_0\\";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    easy3d::SurfaceMesh* part_mesh = easy3d::SurfaceMeshIO::load(args.part_mesh);

    auto original_point_index = part_mesh->vertex_property<int>("v:original_index");
    auto original_face_index = part_mesh->face_property<int>("f:original_index");
    auto fillet_gcp = part_mesh->face_property<int>("f:gcp_labels");
    auto mesh_gcp = mesh->face_property<int>("f:gcp_labels");
    auto coloring_mesh = mesh->face_property<easy3d::vec3>("f:color");
    auto coloring_part_mesh = part_mesh->face_property<easy3d::vec3>("f:color");
    for(auto f : part_mesh->faces()) {
        easy3d::SurfaceMesh::Face ff(original_face_index[f]);
        mesh_gcp[ff] = fillet_gcp[f];
    }
    for(auto v : part_mesh->vertices()) {
        easy3d::SurfaceMesh::Vertex vv(original_point_index[v]);
        mesh->position(vv) = part_mesh->position(v);
    }
    for(auto f : mesh->faces()) {
        if(mesh_gcp[f] == 1) {
            coloring_mesh[f] = easy3d::vec3(1.0, 0,0);
        }
        else {
            coloring_mesh[f] = easy3d::vec3(0, 0,1.0);
        }
    }
    std::string out_path = args.output_dir + "compose.ply";
    easy3d::io::save_ply(out_path, mesh);
    return 0;
}