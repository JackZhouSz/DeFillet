//
// Created by xiaowuga on 3/19/2024.
//

#include <defillet.h>


#include <CLI/CLI.hpp>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/util/file_system.h>
#include <igl/jet.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string fillet_mesh;
        std::string output_dir;
        double beta;
        double gamma;
        int num_opt_iter;
    } args;
    CLI::App app{"DEFILLET Command Line"};

    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-f,--fillet_mesh", args.fillet_mesh, "Fillet Mesh")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    app.add_option("--beta", args.beta, "Beta")->default_val(0.1);
    app.add_option("--gamma", args.gamma, "Gamma")->default_val(1.0);
    app.add_option("--num_opt_iter", args.num_opt_iter, "num_opt_iter")->default_val(10);
    CLI11_PARSE(app, argc, argv);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    easy3d::SurfaceMesh* fillet_mesh = easy3d::SurfaceMeshIO::load(args.fillet_mesh);
    DeFillet defillet;
    defillet.set_mesh(mesh);
    defillet.set_fillet_mesh(fillet_mesh);
    defillet.set_beta(args.beta);
    defillet.set_gamma(args.gamma);
    defillet.set_num_opt_iter(args.num_opt_iter);
    defillet.run_defillet();
    std::string out_defillet_path = args.output_dir +  "defillet.ply";
    easy3d::io::save_ply(out_defillet_path, fillet_mesh, true);
    return 0;
}