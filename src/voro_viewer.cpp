//
// Created by xiaowuga on 2024/11/22.
//

#include "voro_viewer.h"

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/gui/picker_point_cloud.h>
#include <easy3d/gui/picker_surface_mesh.h>
#include <GLFW/glfw3.h>
#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#define EASY3D_MOD_ALT GLFW_MOD_ALT

namespace easy3d {
    void VoroViewer::init(std::vector<easy3d::vec3>& vor_points
                     ,easy3d::SurfaceMesh* mesh
                     ,std::vector<std::vector<int>>& vvns
                     ,std::vector<std::vector<int>>& vvcs) {
        mesh_ = mesh;
        cloud_ = new easy3d::PointCloud;
        int num = vor_points.size();
        for(int i = 0; i < num; i++) {
            cloud_->add_vertex(vor_points[i]);
        }
        vvns_ = vvns;
        vvcs_ = vvcs;
        cloud_->set_name("cloud");
        mesh_->set_name("mesh");
        add_model(mesh_);
        add_model(cloud_);
        cloud_->renderer()->get_points_drawable("vertices")->set_impostor_type(easy3d::PointsDrawable::SPHERE);
        mesh_->renderer()->get_lines_drawable("edges")->set_visible(true);
    }

    bool VoroViewer::mouse_press_event(int x, int y, int button, int modifiers) {
        if (modifiers == EASY3D_MOD_CONTROL) {
            // auto cloud = dynamic_cast<easy3d::PointCloud *>(current_model());
            if (cloud_) {
                easy3d::PointCloudPicker picker(camera());
                // picker.pick_vertice(cloud, polygon_, button == GLFW_MOUSE_BUTTON_RIGHT);
                auto v = picker.pick_vertex(cloud_, x,y);
                mark_selection1(v);

            }
        }
        else if(modifiers == EASY3D_MOD_ALT) {
            if(mesh_) {
                easy3d::SurfaceMeshPicker picker(camera());
                auto face = picker.pick_face(mesh_, x, y);
            }
        }

        return Viewer::mouse_press_event(x, y, button, modifiers);
    }
    bool VoroViewer::key_press_event(int key, int modifiers) {
        if (key == GLFW_KEY_SPACE) {
            bool s = mesh_->renderer()->is_visible();
            mesh_->renderer()->set_visible(!s);
            update();
        }
        else
            return Viewer::key_press_event(key, modifiers);
    }

    void VoroViewer::mark_selection1(easy3d::PointCloud::Vertex v){
        if(v.idx() < 0 ) {
            return;
        }
        auto drawable1 = cloud_->renderer()->get_points_drawable("vertices");
        // auto select = cloud->vertex_property<bool>("v:select");
        auto colors1 = cloud_->vertex_property<easy3d::vec3>("v:color");
        for(auto vv : cloud_->vertices()) {
            colors1[vv] = drawable1->color();
        }
        if(v.idx() < cloud_->n_vertices()) {
            // select[v] = true;
            colors1[v] = easy3d::vec3(1,0,0);
        }
        drawable1->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::VERTEX, "v:color");
        drawable1->update();

        auto drawable2 = mesh_->renderer()->get_triangles_drawable("faces");
        // auto triangle_range = mesh_->get_face_property<std::pair<int, int> >("f:triangle_range");
        auto colors2 = mesh_->face_property<easy3d::vec3>("f:color");
        // for(auto f : mesh_->faces()) {
        //     colors2[f] = drawable1->color();
        // }

        for(auto f : mesh_->faces()) {
            colors2[f] = drawable2->color();
        }
        int num = vvcs_[v.idx()].size();
        std::cout << num << std::endl;
        for(int i = 0; i < num; i++) {
            easy3d::SurfaceMesh::Face face(vvcs_[v.idx()][i]);

            colors2[face] = easy3d::vec3(0,1,0);
        }
        num = vvns_[v.idx()].size();
        std::cout << num << std::endl;
        for(int i = 0; i < num; i++) {
            easy3d::SurfaceMesh::Face face(vvns_[v.idx()][i]);

            // const auto& range = triangle_range[face];
            // drawable2->set_highlight_range(range);
            std::cout << "picked face " << face << std::endl;
            colors2[face] = easy3d::vec3(1,0,0);
        }
        // drawable2->set_highlight(true);
        drawable2->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::FACE, "f:color");
        drawable2->update();

    }


}