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
        double w_convex;
        double w_concave;
        double w1;
        double w2;
    } args;
    CLI::App app{"DEFILLET Command Line"};

    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    app.add_option("-a, --alpha", args.alpha, "Alpha")->default_val(0.5);
    app.add_option("--w_convex", args.w_convex, "w_convex")->default_val(0.08);
    app.add_option("--w_concave", args.w_concave, "w_convex")->default_val(1.0);
    app.add_option("--w1", args.w1, "w1")->default_val(0.3);
    app.add_option("--w2", args.w2, "w2")->default_val(0.4);
    CLI11_PARSE(app, argc, argv);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    bool scores = false;
    bool normal = false;
    for(auto p : mesh->face_properties()) {
        if(p == "f:normal") {
            normal = true;
        }
        if(p == "f:scores") {
            scores = true;
        }
    }
    if(!scores) {
        std::cout << "Please compute face score first." << std::endl;
        return 0;
    }
    if(!normal) {
        std::cout << "Please compute face normal first." << std::endl;
        return 0;
    }
    FilletSeg fillet_seg;
    fillet_seg.set_mesh(mesh);
    fillet_seg.set_alpha(args.alpha);
    fillet_seg.set_w_convex(args.w_convex);
    fillet_seg.set_w_concave(args.w_concave);
    fillet_seg.set_w1(args.w1);
    fillet_seg.set_w2(args.w2);
    fillet_seg.run_gcp();
    auto gcp = mesh->face_property<int>("f:gcp_labels");
    auto coloring_mesh = mesh->face_property<easy3d::vec3>("f:color", easy3d::vec3(0, 0, 0));
    for(auto f : mesh->faces()) {
        if(gcp[f] == 1)
            coloring_mesh[f] = easy3d::vec3(1.0, 0.0, 0.0);
        else
            coloring_mesh[f] = easy3d::vec3(0.0, 0.0, 1.0);
    }
    std::string base_name = easy3d::file_system::base_name(args.input_mesh);
    std::string out_gcp_path = args.output_dir + "gcp.ply";
    easy3d::io::save_ply(out_gcp_path, mesh, true);

    std::string scoring_info_path = args.output_dir +  "gcp_info.json";
    json info;
    info["input_mesh"] = args.input_mesh;
    info["output_dir"] = args.output_dir;
    info["alpha"] = args.alpha;
    info["w_convex"] = args.w_convex;
    info["w_concave"] = args.w_concave;
    info["w1"] = args.w1;
    info["w2"] = args.w2;
    info["gcp_time"] = fillet_seg.get_gcp_time();
    std::ofstream file(scoring_info_path.c_str());
    file << std::setw(4) << info << std::endl;
    return 0;
}