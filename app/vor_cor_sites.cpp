//
// Created by xiaowuga on 2024/10/30.
//
#include <fillet_seg_v2.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/gui/picker_point_cloud.h>
#include <easy3d/renderer/shape.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/util/logging.h>

#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>

#include <GLFW/glfw3.h>
#ifdef __APPLE__
#define EASY3D_MOD_CONTROL GLFW_MOD_SUPER
#else
#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#endif
#define  USE_LASSO
namespace easy3d {
    class PointCloud;
}

class VorCorSites : public easy3d::Viewer
{
public:
    VorCorSites(const std::string& title = "PointSelection") : Viewer(title){

    }
    void init(std::vector<easy3d::vec3>& vor_points
              , easy3d::SurfaceMesh* mesh
              , std::vector<std::vector<int>>& vor_cor_sites);
    void mark_selection(easy3d::PointCloud::Vertex v);

private:
    /// Mouse button press event handler
    bool mouse_press_event(int x, int y, int button, int modifiers) override;


    // std::string usage() const override;

private:
    easy3d::PointCloud* cloud_;
    easy3d::SurfaceMesh* mesh_;
    std::vector<std::vector<int>> vor_cor_sites_;
    std::vector<std::vector<int>>
};

void VorCorSites::init(std::vector<easy3d::vec3>& vor_points
                 ,easy3d::SurfaceMesh* mesh
                 ,std::vector<std::vector<int>>& vor_cor_sites) {
    mesh_ = mesh;
    cloud_ = new easy3d::PointCloud;
    int num = vor_points.size();
    for(int i = 0; i < num; i++) {
        cloud_->add_vertex(vor_points[i]);
    }
    vor_cor_sites_ = vor_cor_sites;
    cloud_->set_name("cloud");
    mesh_->set_name("mesh");
    add_model(cloud_);
    add_model(mesh_);
    mesh_->renderer()->get_lines_drawable("edges")->set_visible(true);
}

bool VorCorSites::mouse_press_event(int x, int y, int button, int modifiers) {
    if (modifiers == EASY3D_MOD_CONTROL) {
        // auto cloud = dynamic_cast<easy3d::PointCloud *>(current_model());
        if (cloud_) {
            easy3d::PointCloudPicker picker(camera());
            // picker.pick_vertice(cloud, polygon_, button == GLFW_MOUSE_BUTTON_RIGHT);
            auto v = picker.pick_vertex(cloud_, x,y);
            mark_selection(v);

        } else
            return Viewer::mouse_press_event(x, y, button, modifiers);
    }
}

void VorCorSites::mark_selection(easy3d::PointCloud::Vertex v){
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
    int num = vor_cor_sites_[v.idx()].size();
    std::cout << num << std::endl;
    for(auto f : mesh_->faces()) {
        colors2[f] = drawable2->color();
    }
    for(int i = 0; i < num; i++) {
        easy3d::SurfaceMesh::Face face(vor_cor_sites_[v.idx()][i]);

        // const auto& range = triangle_range[face];
        // drawable2->set_highlight_range(range);
        std::cout << "picked face " << face << std::endl;
        colors2[face] = easy3d::vec3(1,0,0);
    }
    // drawable2->set_highlight(true);
    drawable2->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::FACE, "f:color");
    drawable2->update();

}


int main () {
    std::string mesh_path = "D:\\code\\sig25_defillet\\fillet.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(mesh_path);
    FilletSegV2 seg(mesh);
    seg.probability_of_sites();
    auto vv = seg.get_vor_vertices();
    auto vor_cor_sites = seg.vor_cor_sites_;
    // for(int i = 0; i < vor_cor_sites.size(); i++) {
    //     if(vor_cor_sites[i].size() != 0) {
    //         std::cout << vor_cor_sites[i].size() << std::endl;
    //     }
    // }
    VorCorSites viewer("VorCorSites");
    viewer.init(vor, mesh, vor_cor_sites);
    viewer.run();
    return 0;
}
