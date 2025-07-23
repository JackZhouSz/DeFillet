#include <iostream>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/renderer.h>
#include <GLFW/glfw3.h>
class MyViewer : public easy3d::Viewer {
public:
    MyViewer(easy3d::SurfaceMesh* mesh, const std::string property_name, const std::string& title = "MyViewer") : Viewer(title), mesh_(mesh){
        pid = mesh_->get_face_property<int>(property_name);
        type = mesh_->face_property<int>("f:type");
        for(auto f : mesh_->faces()) {
            type[f] = 0;
        }
        std::set<int> s(pid.vector().begin(), pid.vector().end());
        num = s.size();
        id = -1;
        add_model(mesh_);
    }

private:
    /// Mouse button press event handler
    // bool mouse_press_event(int x, int y, int button, int modifiers) override;
    bool key_press_event(int key, int modifiers) override;
private:
    easy3d::SurfaceMesh* mesh_;
    easy3d::SurfaceMesh::FaceProperty<int> pid;
    easy3d::SurfaceMesh::FaceProperty<int> type;
    int num;
    int id;
};
bool MyViewer::key_press_event(int key, int modifiers) {
    if (key == GLFW_KEY_DOWN) {
        id++;
        std::cout << id << ' ' << num << std::endl;
        if(id == num) {
            exit();
        }
        auto drawable = mesh_->renderer()->get_triangles_drawable("faces");
        auto colors = mesh_->face_property<easy3d::vec3>("f:color");
        for(auto f : mesh_->faces()) {
            colors[f] = drawable->color();
        }
        for(auto f : mesh_->faces()) {
            if(pid[f] == id) {
                colors[f] = easy3d::vec3(1.0,0,0);
            }
        }
        drawable->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::FACE, "f:color");
        drawable->update();
    }
    else if(key == GLFW_KEY_0) {
        for(auto f : mesh_->faces()) {
            if(pid[f] == id) {
                type[f] = 0;
            }
        }
        std::cout << "plane" << ' ' << id << " done." <<std::endl;
    }
    else if(key == GLFW_KEY_1) {
        for(auto f : mesh_->faces()) {
            if(pid[f] == id) {
                type[f] = 1;
            }
        }
        std::cout << "cyinder" << ' ' << id << " done." <<std::endl;
    }
    else if(key == GLFW_KEY_2) {
        for(auto f : mesh_->faces()) {
            if(pid[f] == id) {
                type[f] = 2;
            }
        }
        std::cout << "sphere" << ' ' << id << " done." <<std::endl;
    }
    else if(key == GLFW_KEY_3) {
        for(auto f : mesh_->faces()) {
            if(pid[f] == id) {
                type[f] = 3;
            }
        }
        std::cout << "qudric" << ' ' << id << " done." <<std::endl;
    }
    else
        return Viewer::key_press_event(key, modifiers);
}


bool load_clr(easy3d::SurfaceMesh* model,
              std::string property_name,
              std::string path) {
    auto pid = model->face_property<int>(property_name);
    std::vector<int> pfid;
    std::vector<easy3d::vec3> colors;
    std::ifstream in(path);
    if (in.good()) {
        int nv, nf, np, m;
        in >> nv >> nf >> np >> m;
        double box;
        for (int i = 0; i < 7; i++)
            in >> box;
        for (int i = 0; i < np; i++) {
            int id, proxy, pnf;
            in >> id >> proxy >> pnf;
            pfid.push_back(proxy);
            for (int j = 0; j < pnf; j++) {
                int fid;
                in >> fid;
                pid[easy3d::SurfaceMesh::Face(fid)] = id;
            }
            easy3d::vec3 c;
            in >> c[0] >> c[1] >> c[2];
            colors.push_back(c);
        }
    }
    bool success = in.good();
    in.close();
    if(success) {
        return true;
    }
    else {
        model->remove_face_property(property_name);
        return false;
    }
}



int main() {
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load("D:\\Sig25_DeFillet\\compareFI\\rgb\\0_sed.obj");
    load_clr(mesh, "f:pid", "D:\\Sig25_DeFillet\\compareFI\\0_qse.clr");
    MyViewer viewer( mesh, "f:pid");
    viewer.run();
    auto type = mesh->face_property<int>("f:type");
    auto color = mesh->face_property<easy3d::vec3>("f:color");
    for(auto f : mesh->faces()) {
        // std::cout << type[f] <<std::endl;
        if(type[f] == 0) {
            color[f] = easy3d::vec3(0.0,0,1.0);
        }
        else if(type[f] == 1) {
            color[f] = easy3d::vec3(1.0,1.0,0);
        }
        else if(type[f] == 2) {
            color[f] = easy3d::vec3(0.0,1.0,0);
        }
        else {
            color[f] = easy3d::vec3(1.0,0.0,0);
        }
    }
    easy3d::SurfaceMeshIO::save("../out/4_qse.ply", mesh);
    return 0;
}