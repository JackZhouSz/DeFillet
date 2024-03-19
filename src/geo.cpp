//
// Created by 13900K on 2024/3/19.
//

#include "defillet.h"


#include <CLI/CLI.hpp>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/util/file_system.h>
#include <igl/jet.h>

#include <nlohmann/json.hpp>


#include <easy3d/viewer/viewer.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/texture_manager.h>

using json = nlohmann::json;

int main(int argc, char **argv) {
    struct {
        std::string input_mesh;
        std::string output_dir;
    } args;
    CLI::App app{"DEFILLET Command Line"};

    app.add_option("-i,--input_mesh", args.input_mesh, "Input Mesh")->required();
    app.add_option("-o,--output_dir", args.output_dir, "Output Directory")->required();
    CLI11_PARSE(app, argc, argv);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(args.input_mesh);
    bool gcp_labels = false;
    for(auto p : mesh->face_properties()) {
        if(p == "f:gcp_labels") {
            gcp_labels = true;
        }
    }
    if(!gcp_labels) {
        std::cout << "Please gcp_label first." << std::endl;
        return 0;
    }
    DeFillet defillet;
    defillet.set_mesh(mesh);
    defillet.run_geodesic();
    easy3d::SurfaceMesh* fillet_mesh = defillet.get_fillet_mesh();
    std::string out_geo_path = args.output_dir +  "geo.ply";
    easy3d::io::save_ply(out_geo_path, fillet_mesh, false);
    easy3d::Viewer viewer("asdsad");
    viewer.add_model(fillet_mesh);
    auto geo_dis = fillet_mesh->get_vertex_property<float>("v:geo_dis");

    auto drawable = fillet_mesh->renderer()->get_triangles_drawable("faces");
    drawable->set_scalar_coloring(easy3d::State::VERTEX, "v:geo_dis", nullptr, 0.0f, 0.0f);

    const std::string texture_file = "D:\\code\\defillet\\lib-easy3d\\resources\\colormaps\\rainbow.png";
    easy3d::Texture *texture = easy3d::TextureManager::request(texture_file);
    if (!texture) {
        return 0;
    }
    std::cout << "ASD" <<std::endl;
    // Use the texture
    drawable->set_texture(texture);
    auto face_tar_normals = fillet_mesh->face_property<easy3d::vec3>("f:tar_normals");
    const easy3d::Box3 &box = fillet_mesh->bounding_box();
    float length = norm(box.max_point() - box.min_point()) * 0.02f;
    std::vector<easy3d::vec3> tmp;
    for(auto f : fillet_mesh->faces()) {
        int num = 0;
        easy3d::vec3 center(0,0,0);
        for(auto v : fillet_mesh->vertices(f)) {
            center += fillet_mesh->position(v);
            num++;
        }
        center /= num;
        tmp.emplace_back(center);
        tmp.emplace_back(center + face_tar_normals[f] * length);
    }
    easy3d::LinesDrawable* line = new easy3d::LinesDrawable;
    line->update_vertex_buffer(tmp);
    viewer.add_drawable(line);
    viewer.update();
    viewer.run();
    


    return 0;
}