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

int seg_via_sharp_feature(easy3d::SurfaceMesh* mesh, double angle) {
    auto segs = mesh->add_face_property<int>("f:seg_label", -1);
    auto normals = mesh->face_property<easy3d::vec3>("f:normal");
    for(auto f : mesh->faces()) {
        segs[f] = -1;
        normals[f] = mesh->compute_face_normal(f);
    }
    int seg_id = 0;
    double maxx = 0;
    for(auto f : mesh->faces()) {
        if(segs[f] == -1) {
            std::queue<easy3d::SurfaceMesh::Face> que;
            que.push(f);
            while(!que.empty()) {
                auto cur = que.front(); que.pop();
                if(segs[cur] != -1)
                    continue;
                segs[cur] = seg_id;
                for(auto h : mesh->halfedges(cur)) {
                    auto nxt_f = mesh->face(mesh->opposite(h));
                    if(nxt_f.is_valid() && segs[nxt_f] == -1) {
                        double dot_val = easy3d::dot(normals[cur], normals[nxt_f]);
                        dot_val = std::clamp(dot_val, -0.99999,0.99999);
                        double di_angle = std::acos(dot_val) * 180 / M_PI;
                        if(di_angle < angle) {
                            que.push(nxt_f);
                        }
                    }
                }
            }
            seg_id++;
        }
    }

    return seg_id;
}

int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string output_dir;
    } args;

//    CLI::App app{"DEFILLET Command Line"};
//
//    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
//    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
//    CLI11_PARSE(app, argc, argv);
    args.input_mesh = "E:\\defillet_case\\21767_53eaa558_12\\21767_53eaa558_12.ply";
    args.output_dir = "E:\\defillet_case\\21767_53eaa558_12\\21767_53eaa558_12_parts\\";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    int num = seg_via_sharp_feature(mesh, 30);
    auto segs = mesh->get_face_property<int>("f:seg_label");
    std::string base_name = easy3d::file_system::base_name(args.input_mesh);
    easy3d::SurfaceMeshSegmenter segmenter(mesh);
    std::cout << num << std::endl;
    for(int i = 0; i < num; i++) {
        easy3d::SurfaceMesh* part_mesh = segmenter.segment(segs, i);
        std::string out_path = args.output_dir + base_name + "_part_" + std::to_string(i) + ".ply";
        easy3d::io::save_ply(out_path, part_mesh, false);
    }
    return 0;
}
