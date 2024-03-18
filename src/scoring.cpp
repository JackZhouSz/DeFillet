//
// Created by 13900K on 2024/3/15.
//

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
        double eps;
        double radius;
        double s;
        double min_score;
        double std_ratio;
        int num_sor_iter;
        int nb_neighbors;
    } args;

    CLI::App app{"DEFILLET Command Line"};

    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    app.add_option("-e,--eps", args.eps, "Eps")->default_val(0.03);
    app.add_option("-r,--radius", args.radius, "Radius")->default_val(0.06);
    app.add_option("-s,--s", args.s, "S")->default_val(10.0);
    app.add_option("--min_score", args.min_score, "min_score")->default_val(0.5);
    app.add_option("--nb_neighbors", args.nb_neighbors, "nb_neighbors")->default_val(30);
    app.add_option("--num_sor_iter", args.num_sor_iter, "num_sor_iter")->default_val(3);
    app.add_option("--std_ratio", args.std_ratio, "std_ratio")->default_val(0.3);
    CLI11_PARSE(app, argc, argv);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);

    FilletSeg fillet_seg;
    fillet_seg.set_mesh(mesh);
    fillet_seg.set_eps(args.eps);
    fillet_seg.set_radius(args.radius);
    fillet_seg.set_s(args.s);
    fillet_seg.set_min_score(args.min_score);
    fillet_seg.run_scoring();
    std::string base_name = easy3d::file_system::base_name(args.input_mesh);
    std::string out_mesh_path = args.output_dir + "mesh_scoring.ply";
    auto scores = mesh->face_property<float>("f:scores");
    int nb_faces = mesh->n_faces();
    Eigen::VectorXd Z(nb_faces);
    for (int i = 0; i < nb_faces; i++) {
        Z[i] = scores.vector()[i];
    }
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto coloring_mesh = mesh->face_property<easy3d::vec3>("f:color", easy3d::vec3(0, 0, 0));
    for(auto f : mesh->faces()) {
        coloring_mesh[f] = easy3d::vec3(Ct(f.idx(), 0),
                                        Ct(f.idx(), 1), Ct(f.idx(), 2));
    }
    easy3d::io::save_ply(out_mesh_path, mesh, true);
    easy3d::PointCloud* sites = new easy3d::PointCloud;
    auto& tmp_sites = fillet_seg.get_sites();
    for(auto v : tmp_sites) {
        sites->add_vertex(v);
    }
    auto coloring_sites= sites->vertex_property<easy3d::vec3>("v:color", easy3d::vec3(0, 0, 0));
    for(auto v : sites->vertices()) {
        coloring_sites[v] = easy3d::vec3(Ct(v.idx(), 0),
                                         Ct(v.idx(), 1), Ct(v.idx(), 2));
    }
    std::string out_sites_path = args.output_dir + "sites_scoring.ply";
    easy3d::io::save_ply(out_sites_path, sites, true);
    easy3d::PointCloud* vertices = new easy3d::PointCloud;
    auto& tmp_vertices = fillet_seg.get_vertices();
    auto& tmp_vertices_scores = fillet_seg.get_vertices_scores();
    auto& tmp_index = fillet_seg.get_vertices_index_inside_box();

    int num = tmp_index.size();
    Z.resize(num);
    for(int i = 0; i < num; i++) {
        vertices->add_vertex(tmp_vertices[tmp_index[i]]);
        Z[i] = tmp_vertices_scores[tmp_index[i]];
    }
    igl::jet(Z, true, Ct);
    auto coloring_vertices = vertices->vertex_property<easy3d::vec3>("v:color", easy3d::vec3(0, 0, 0));
    for(auto v : vertices->vertices()) {
        coloring_vertices[v] = easy3d::vec3(Ct(v.idx(), 0),
                                            Ct(v.idx(), 1), Ct(v.idx(), 2));
    }
    std::string out_vertices_path = args.output_dir + "vertices_scoring.ply";
    easy3d::io::save_ply(out_vertices_path, vertices, true);

    easy3d::PointCloud* sor_vertices = new easy3d::PointCloud;
    const std::vector<int> sor_labels = fillet_seg.get_sor_labels();
    for(int i = 0; i < sor_labels.size(); i++) {
        if(sor_labels[i]) {
//            std::cout <<"ASD" <<std::endl;
            sor_vertices->add_vertex(tmp_vertices[i]);
        }
    }
    std::string out_vertices_sor_path = args.output_dir +  "vertices_sor.ply";
    easy3d::io::save_ply(out_vertices_sor_path, sor_vertices, true);

    std::string scoring_info_path = args.output_dir + "scoring_info.json";
    json info;
    info["input_mesh"] = args.input_mesh;
    info["output_dir"] = args.output_dir;
    info["eps"] = args.eps;
    info["radius"] = args.radius;
    info["s"] = args.s;
    info["min_score"] = args.min_score;
    info["std_ratio"] = args.std_ratio;
    info["nb_neighbors"] = args.nb_neighbors;
    info["num_sor_iter"] = args.num_sor_iter;
    info["sor_time"] = fillet_seg.get_sor_time();
    info["scoring_time"] = fillet_seg.get_scoring_time();
    std::ofstream file(scoring_info_path.c_str());
    file << std::setw(4) << info << std::endl;
    return 0;
}