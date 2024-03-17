#include <CLI/CLI.hpp>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/util/file_system.h>
#include <fillet_seg.h>
#include <igl/jet.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string output_dir;
        double alpha;
    } args;
    CLI::App app{"DEFILLET Command Line"};

    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    app.add_option("-a, --alpha", args.alpha, "Alpha")->default_val(0.5);
    CLI11_PARSE(app, argc, argv);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    for(auto p : mesh->face_properties()) {
        std::cout << p << std::endl;
    }
    return 0;
}