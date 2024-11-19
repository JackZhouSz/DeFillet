//
// Created by xiaowuga on 2024/11/6.
//

#include <dkm_parallel.hpp>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/point_cloud_io.h>

#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/renderer/shape.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/gui/picker_point_cloud.h>
#include <easy3d/gui/picker_surface_mesh.h>
#include <fillet_seg_v4.h>
#include <igl/jet.h>
#include <GLFW/glfw3.h>

#include <numeric>
#include <knn.h>

#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#define EASY3D_MOD_ALT GLFW_MOD_ALT

void save_point_cloud(std::vector<easy3d::vec3>& points, std::string out_path) {
    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(auto p : points) {
        cloud->add_vertex(p);
    }
    easy3d::PointCloudIO::save(out_path, cloud);
}

void save_point_cloud_with_scalar_field(std::vector<easy3d::vec3>& points
                                       ,std::vector<float>& field
                                       ,std::string out_path) {
    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(auto p : points) {
        cloud->add_vertex(p);
    }
    int num = cloud->n_vertices();
    Eigen::VectorXd Z(num);
    for (int i = 0; i < num; i++) {
        Z[i] = field[i];
    }
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto coloring = cloud->vertex_property<easy3d::vec3>("v:color", easy3d::vec3(0, 0, 0));
    for(auto v : cloud->vertices()) {
        coloring[v] = easy3d::vec3(Ct(v.idx(), 0),
                                        Ct(v.idx(), 1), Ct(v.idx(), 2));
    }
    easy3d::PointCloudIO::save(out_path, cloud);
}

void tetrahedron_circumcenter(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius) {
    int num = points.size();
    if(num != 4) {
        std::cerr << "points of tetrahedron do not equal 4" << std::endl;
        return;
    }
    Eigen::MatrixXd d_matrix(4, 4); d_matrix.setOnes();
    // Multiplied by (x²+y²+z²)
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = points[i].x;
        d_matrix(i, 1) = points[i].y;
        d_matrix(i, 2) = points[i].z;
    }
    float M11 = d_matrix.determinant();
    // Multiplied by x
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].y;
        d_matrix(i, 2) = points[i].z;
    }
    float M12 = d_matrix.determinant();

    // std::cout << M11 << std::endl;
    // std::cout << "bb" << std::endl;

    // Multiplied by y
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].z;
    }
    float M13 = d_matrix.determinant();

    // Multiplied by z
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].y;
    }
    float M14 = d_matrix.determinant();

    // Multiplied by 1
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].y;
        d_matrix(i, 3) = points[i].z;
    }

    float M15 = d_matrix.determinant();

    center = easy3d::vec3(0.5 * (M12 / M11), -0.5 * (M13 / M11), 0.5 * (M14 / M11));
    radius = std::sqrt(easy3d::dot(center, center) - (M15 / M11));

}


namespace easy3d {
    class PointCloud;
}

class ASDViewer : public easy3d::Viewer
{
public:
    ASDViewer(const std::string& title = "PointSelection") : Viewer(title){

    }
    void init(std::vector<easy3d::vec3>& vor_points
              , easy3d::SurfaceMesh* mesh
              , std::vector<std::vector<int>>& ccs);
    void mark_selection1(easy3d::PointCloud::Vertex v);
    // void mark_selection2(easy3d::SurfaceMesh::Face f);

private:
    /// Mouse button press event handler
    bool mouse_press_event(int x, int y, int button, int modifiers) override;


    // std::string usage() const override;

private:
    easy3d::PointCloud* cloud_;
    easy3d::SurfaceMesh* mesh_;
    std::vector<std::vector<int>> ccs_;
    KNN::KdSearch* kdt;
    std::vector<KNN::Point> knn_points;
};

void fit_sphere(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius, int max_iter = 100) {
    int num = points.size();
    std::vector<int> numbers(num);
    std::iota(numbers.begin(), numbers.end(), 0);  // 填充范围内的整数
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::vector<easy3d::vec3> samples(4);
    float loss = 1e9;
    for(int i = 0; i < max_iter; i++) {
        std::shuffle(numbers.begin(), numbers.end(),std::default_random_engine(seed));
        for(int j = 0; j < 4; j++) {
            samples[j] = points[numbers[j]];
        }
        easy3d::vec3 center_candidate;
        float radius_candidate = -1;

        tetrahedron_circumcenter(samples, center_candidate, radius_candidate);
        float err = 0;
        for(int j = 0; j < num; j++) {
            err += fabs((points[j] - center_candidate).norm() - radius_candidate);
        }
        err /= num;
        if(err < loss) {
            loss = err;
            center = center_candidate;
            radius = radius_candidate;
        }
    }
}

void ASDViewer::init(std::vector<easy3d::vec3>& vor_points
                 ,easy3d::SurfaceMesh* mesh
                 ,std::vector<std::vector<int>>& ccs) {
    mesh_ = mesh;
    cloud_ = new easy3d::PointCloud;
    int num = vor_points.size();
    for(int i = 0; i < num; i++) {
        cloud_->add_vertex(vor_points[i]);
    }

    cloud_->set_name("cloud");
    mesh_->set_name("mesh");
    ccs_ = ccs;
    add_model(mesh_);
    add_model(cloud_);
    cloud_->renderer()->get_points_drawable("vertices")->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    mesh_->renderer()->get_lines_drawable("edges")->set_visible(true);


    for(int i = 0; i < vor_points.size(); i++) {
        knn_points.emplace_back(KNN::Point(vor_points[i].x,
                                                       vor_points[i].y, vor_points[i].z));
    }
    kdt = new KNN::KdSearch(knn_points);
}

bool ASDViewer::mouse_press_event(int x, int y, int button, int modifiers) {
    if (modifiers == EASY3D_MOD_CONTROL) {
        // auto cloud = dynamic_cast<easy3d::PointCloud *>(current_model());
        if (cloud_) {
            easy3d::PointCloudPicker picker(camera());
            // picker.pick_vertice(cloud, polygon_, button == GLFW_MOUSE_BUTTON_RIGHT);
            auto v = picker.pick_vertex(cloud_, x, y);
            mark_selection1(v);

        }
    }
    else if(modifiers == EASY3D_MOD_ALT) {
        if(mesh_) {
            easy3d::SurfaceMeshPicker picker(camera());
            auto face = picker.pick_face(mesh_, x, y);
            // mark_selection2(face);
        }
    }

    return Viewer::mouse_press_event(x, y, button, modifiers);
}

void ASDViewer::mark_selection1(easy3d::PointCloud::Vertex v){
    if(v.idx() < 0 ) {
        return;
    }
    std::vector<easy3d::vec3> points;

    std::vector<easy3d::vec3> pp;
    float radius;
    easy3d::vec3 center = cloud_->position(v);
    for(int iter = 0; iter < 10; iter++) {
        KNN::Point knn_p(center.x, center.y, center.z);
        std::vector<size_t> tmp_indices;
        std::vector<double> dist;
        kdt->kth_search(knn_p, 10, tmp_indices, dist);

        for(int i = 0; i < tmp_indices.size(); i++) {
            for(int j = 0; j < ccs_[tmp_indices[i]].size(); j++) {
                auto f = easy3d::SurfaceMesh::Face(ccs_[tmp_indices[i]][j]);
                easy3d::vec3 cc(0,0,0);
                int cct = 0;
                for(auto vv : mesh_->vertices(f)) {
                    cc += mesh_->position(vv);
                    cct++;
                }
                cc /= cct;
                points.emplace_back(cc);
            }
        }

        // std::cout << "fitting.." << std::endl;
        fit_sphere(points, center, radius);
        pp.emplace_back(center);
        // std::cout << center << ' ' << radius << std::endl;
    }
    auto drawable3 = new easy3d::PointsDrawable("DD");
    drawable3->update_vertex_buffer({center});
    drawable3->set_point_size(10);
    drawable3->set_color(easy3d::vec4(1,1,0,1));
    drawable3->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    add_drawable(drawable3);

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

    for(auto f : mesh_->faces()) {
        colors2[f] = drawable2->color();
    }
    for(int j = 0; j < ccs_[v.idx()].size(); j++) {
        easy3d::SurfaceMesh::Face f(ccs_[v.idx()][j]);
            colors2[f] = easy3d::vec3(0,1,0);

    }
    drawable2->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::FACE, "f:color");
    drawable2->update();

}

int main() {
    // easy3d::PointCloudIO::save("../dkm_cloud_clean_kmeans.ply", center);
    std::string path = "D:\\code\\sig25_defillet\\20440_27177360_1.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(path);

    FilletSegV4 fillet_seg_v4(mesh);
    fillet_seg_v4.seg();
    // auto cc = fillet_seg_v4.cc_;
    // auto ccr = fillet_seg_v4.cr_;
    // auto ccvvn = fillet_seg_v4.ccvvn_;
    // auto ccsn = fillet_seg_v4.ccsn_;
    // auto ccs = fillet_seg_v4.ccs_;
    auto vv = fillet_seg_v4.vv_;
    auto vvcs = fillet_seg_v4.vvcs_;
    // save_point_cloud(cc, "../../sig25_defillet/cc.ply");
    ASDViewer viewer("ASD");
    viewer.init(vv, mesh,vvcs);
    viewer.run();
    // save_point_cloud_with_scalar_field(cc, ccr,"../../sig25_defillet/22181_526fb582_0_ccr.ply");
    // save_point_cloud_with_scalar_field(cc, ccvvn,"../../sig25_defillet/22181_526fb582_0_ccvvn.ply");
    // save_point_cloud_with_scalar_field(cc, ccsn,"../../sig25_defillet/22181_526fb582_0_ccsn.ply");

    return 0;
}