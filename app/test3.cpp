//
// Created by xiaowuga on 2024/10/29.
//


#include <fillet_seg_v3.h>
#include <numeric>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>


#include <easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/gui/picker_point_cloud.h>
#include <easy3d/gui/picker_surface_mesh.h>
#include <easy3d/renderer/shape.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/util/logging.h>
#include <knn.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>

#include <GLFW/glfw3.h>

#include <igl/jet.h>

#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#define EASY3D_MOD_ALT GLFW_MOD_ALT

void sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels) {
    int nb_points = points.size();
    int nb_neighbors = 10, num_sor_iter = 3;
    float std_ratio = 0.3;
    for(int i = 0; i < num_sor_iter; i++) {
        std::vector<KNN::Point> knn_points;
        std::vector<int> indices;
#pragma omp parallel for
        for(int j = 0; j < nb_points; j++) {
            if(labels[j]) {
                #pragma omp critical
                {
                    knn_points.emplace_back(KNN::Point(points[j].x,
                                                       points[j].y, points[j].z));
                    indices.emplace_back(j);
                }
            }
        }
        KNN::KdSearch kds(knn_points);
        int num = knn_points.size();
        size_t valid_distances = 0;
        std::vector<double> avg_distances(num);
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(knn_points[j], nb_neighbors, tmp_indices, dist);
            double mean = -1.0;

            if(dist.size() > 0u) {
                valid_distances++;
                std::for_each(dist.begin(), dist.end(),
                              [](double &d) { d = std::sqrt(d); });
                mean = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();
            }
            avg_distances[j] = mean;
        }
        if(valid_distances == 0) {
            continue;
        }
        double cloud_mean = std::accumulate(
        avg_distances.begin(), avg_distances.end(), 0.0,
        [](double const &x, double const &y) { return y > 0 ? x + y : x; });

        cloud_mean /= valid_distances;
        double sq_sum = std::inner_product(
        avg_distances.begin(), avg_distances.end(), avg_distances.begin(),
        0.0, [](double const &x, double const &y) { return x + y; },
        [cloud_mean](double const &x, double const &y) {
            return x > 0 ? (x - cloud_mean) * (y - cloud_mean) : 0;
        });
        double std_dev = std::sqrt(sq_sum / (valid_distances - 1));
        double distance_threshold = cloud_mean + std_ratio * std_dev;
        // std::cout << distance_threshold << std::endl;
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            // std::cout << avg_distances[j] << std::endl;
            if(avg_distances[j] > 0 && avg_distances[j] < distance_threshold) {
                labels[indices[j]] = true;
            } else {
                labels[indices[j]] = false;
            }
        }
    }

}

void find_cluster_center(const std::vector<easy3d::vec3>& data, easy3d::vec3& centroid) {
    int num = data.size();
    if(num < 10) {
        std::cout << "point too few" <<std::endl;
        return;
    }
    std::vector<bool> labels(num, true);
    sor(data, labels);
    centroid = easy3d::vec3(0,0,0);
    std::vector<easy3d::vec3> clean_data;
    for(int i = 0; i < num; i++) {
        if(labels[i]) {
            centroid += data[i];
            clean_data.emplace_back(data[i]);
        }
    }
    num = clean_data.size();
    centroid /= num;

    std::vector<float> dist(num);
    for(int i = 0; i < num; i++) {
        dist[i] = (clean_data[i] - centroid).length2();
    }

    std::vector<float> tmp = dist;
    std::nth_element(tmp.begin(), tmp.end(), tmp.begin() + num / 2);
    double med = tmp[num / 2];
    for(int i = 0; i < num; i++) {
        tmp[i] = fabs(dist[i] - med);
    }
    std::nth_element(tmp.begin(), tmp.end(), tmp.begin() + num / 2);
    double mad = tmp[num / 2];
    easy3d::vec3 new_centroid(0,0,0);
    int ct = 0;
    std::vector<easy3d::vec3> dd;
    for(int i = 0; i < num; i++) {
        if(dist[i] - med < 5.2 * mad) {
            new_centroid += clean_data[i];
            dd.emplace_back(clean_data[i]);
            ct++;
        }
    }
    std::cout << ct << ' ' << new_centroid << std::endl;
    new_centroid /= ct;
    centroid = new_centroid;
    float xx = 0;
    for(size_t i = 0; i < dd.size(); i++) {
        xx += (dd[i] - new_centroid).length2();
    }
    xx /= dd.size();
    // std::cout << xx << std::endl;

}

void sqdist_omp(Eigen::MatrixXd& a,Eigen::MatrixXd& b,
                   Eigen::MatrixXd& res) {
    int n = b.cols(), m = a.cols();
    res.resize(n, m);
#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res(i,j) = (b.col(i) - a.col(j)).squaredNorm();
        }
    }
}

bool fuzzy_clustering(Eigen::MatrixXd& P, Eigen::MatrixXd& C, int nb_cluster, bool robust = true) {
    int nb_points = P.cols();
    int dim = P.rows();
    if(nb_cluster < 1 || nb_cluster > nb_points) {
        std::cerr << "Incorrect setting for nb_cluster" << std::endl;
        return false;
    }

    Eigen::VectorXd onesDPx1 = Eigen::VectorXd::Ones(dim); // DP: dim of points
    Eigen::RowVectorXd ones1xNP = Eigen::RowVectorXd::Ones(nb_points); // NP: number of points
    Eigen::VectorXd onesNCx1 = Eigen::VectorXd::Ones(nb_cluster); // NC: number of cluster

    int maxIter = 100; // maximum number of iterations
    float minDiff = 1e-6; // minimum difference

    Eigen::MatrixXd U = Eigen::MatrixXd::Random(nb_cluster, nb_points);
    U = U.array() / (onesNCx1 * U.colwise().sum()).array();
    Eigen::MatrixXd U2 = U.array().square();
    Eigen::VectorXd rowSumU2 = U2.rowwise().sum(); // nb_cluster x 1
    C = (P * U2.transpose()).array() / (onesDPx1 * rowSumU2.transpose()).array();
    Eigen::MatrixXd dist2; // nb_cluster x nb_points
    sqdist_omp(P, C, dist2);

    // return true;
    double obj = (U2.array() * dist2.array()).sum();
    for(int iter = 0; iter < maxIter; iter++) {
        float obj_old = obj;
        Eigen::MatrixXd tmp = dist2.array().inverse();
        U = tmp.array() / (onesNCx1 * tmp.colwise().sum()).array();
        U2 = U.array().square();
        rowSumU2 = U2.rowwise().sum();
        C = (P * U2.transpose()).array() / (onesDPx1 * rowSumU2.transpose()).array();
        sqdist_omp(P, C, dist2);
        Eigen::MatrixXd U2D2 = U2.array() * dist2.array();
        obj = U2D2.sum();
        if(fabs(obj - obj_old) < minDiff) {
            break;
        }

    }
    return true;
}

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
              , std::vector<std::vector<int>>& vvns
              , std::vector<std::vector<int>>& vvcs
              , std::vector<std::vector<int>>& scvv);
    void mark_selection1(easy3d::PointCloud::Vertex v);
    void mark_selection2(easy3d::SurfaceMesh::Face f);

private:
    /// Mouse button press event handler
    bool mouse_press_event(int x, int y, int button, int modifiers) override;


    // std::string usage() const override;

private:
    easy3d::PointCloud* cloud_;
    easy3d::SurfaceMesh* mesh_;
    std::vector<std::vector<int>> vvns_;
    std::vector<std::vector<int>> vvcs_;
    std::vector<std::vector<int>> scvv_;
};

void VorCorSites::init(std::vector<easy3d::vec3>& vor_points
                 ,easy3d::SurfaceMesh* mesh
                 ,std::vector<std::vector<int>>& vvns
                 ,std::vector<std::vector<int>>& vvcs
                 ,std::vector<std::vector<int>>& scvv) {
    mesh_ = mesh;
    cloud_ = new easy3d::PointCloud;
    int num = vor_points.size();
    for(int i = 0; i < num; i++) {
        cloud_->add_vertex(vor_points[i]);
    }
    vvns_ = vvns;
    vvcs_ = vvcs;
    scvv_ = scvv;
    cloud_->set_name("cloud");
    mesh_->set_name("mesh");

    add_model(mesh_);
    add_model(cloud_);
    cloud_->renderer()->get_points_drawable("vertices")->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    mesh_->renderer()->get_lines_drawable("edges")->set_visible(true);
}

bool VorCorSites::mouse_press_event(int x, int y, int button, int modifiers) {
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
            mark_selection2(face);
        }
    }

    return Viewer::mouse_press_event(x, y, button, modifiers);
}

void VorCorSites::mark_selection1(easy3d::PointCloud::Vertex v){
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


void VorCorSites::mark_selection2(easy3d::SurfaceMesh::Face f){
    if(f.idx() < 0 ) {
        return;
    }
    auto drawable1 = cloud_->renderer()->get_points_drawable("vertices");
    // auto select = cloud->vertex_property<bool>("v:select");
    auto colors1 = cloud_->vertex_property<easy3d::vec3>("v:color");
    for(auto vv : cloud_->vertices()) {
        colors1[vv] = drawable1->color();
    }
    easy3d::PointCloud* dsa = new easy3d::PointCloud;
    for(size_t i = 0; i < scvv_[f.idx()].size(); i++) {
        easy3d::PointCloud::Vertex v(scvv_[f.idx()][i]);
        colors1[v] = easy3d::vec3(0,1,0);
        dsa->add_vertex(cloud_->position(v));
    }
    easy3d::vec3 s(0,0,0);
    int ct = 0;
    for(auto v : mesh_->vertices(f)) {
        s += mesh_->position(v);
        ct++;
    }
    s /= ct;


    // Eigen::MatrixXd data, centroids;
    // Eigen::VectorXi labels;
    // data.resize(3,dsa->n_vertices());
    // for(auto v : dsa->vertices()) {
    //     auto p = dsa->position(v);
    //     float r =  (p - s).norm();
    //     data.col(v.idx()) = Eigen::Vector3d(p.x,p.y,p.z);
    // }
    //
    // if(dsa->n_vertices() > 10) {
    //     fuzzy_clustering(data, centroids, 1);
    //     easy3d::PointsDrawable* dd = new easy3d::PointsDrawable("ASD");
    //     std::vector<easy3d::vec3> c;
    //     for(int i = 0; i < centroids.cols(); i++) {
    //         c.emplace_back(easy3d::vec3(centroids(0, i), centroids(1, i), centroids(2, i)));
    //         std::cout << (c.back() - s).norm() << std::endl;
    //         // c.emplace_back(easy3d::vec3(centroids(i, 0), centroids(i, 1), centroids(i, 2)));
    //     }
    //     dd->update_vertex_buffer({c});
    //     dd->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    //     dd->set_point_size(15);
    //     dd->set_color(easy3d::vec4(0,0,1,1));
    //     add_drawable(dd);
    //     easy3d::PointCloudIO::save("../scvv.ply",dsa);
    //
    // }
    easy3d::vec3 centroid;
    if(dsa->n_vertices() > 10) {
        find_cluster_center(dsa->points(), centroid);
        easy3d::PointsDrawable* dd = new easy3d::PointsDrawable("ASD");
        dd->update_vertex_buffer({centroid});
        dd->set_impostor_type(easy3d::PointsDrawable::SPHERE);
        dd->set_point_size(15);
        dd->set_color(easy3d::vec4(0,0,1,1));
        add_drawable(dd);
        std::cout << (centroid - s).norm() << std::endl;
        easy3d::PointCloudIO::save("../scvv.ply",dsa);

    }
    drawable1->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::VERTEX, "v:color");
    drawable1->update();
    auto drawable2 = mesh_->renderer()->get_triangles_drawable("faces");
    // auto triangle_range = mesh_->get_face_property<std::pair<int, int> >("f:triangle_range");
    auto colors2 = mesh_->face_property<easy3d::vec3>("f:color");
    // for(auto f : mesh_->faces()) {
    //     colors2[f] = drawable1->color();
    // }

    for(auto ff : mesh_->faces()) {
        colors2[ff] = drawable2->color();
    }
    colors2[f] = easy3d::vec3(1,0,0);
    // drawable2->set_highlight(true);
    drawable2->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::FACE, "f:color");
    drawable2->update();

}

int main() {
    std::string mesh_path = "D:\\code\\sig25_defillet\\asd.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(mesh_path);
    FilletSegV3 fillet_seg_v3(mesh);
    fillet_seg_v3.solve();
    auto vv = fillet_seg_v3.vv_;
    auto vvns = fillet_seg_v3.vvns_;
    auto vvcs = fillet_seg_v3.vvcs_;
    auto scvv = fillet_seg_v3.scvv_;
    auto vvs = fillet_seg_v3.vvs_;
    auto sr = fillet_seg_v3.sr_;
    int num = mesh->n_faces();
    Eigen::VectorXd Z(num);
    std::vector<int> face_indices;
    int ct = 0;
    for (int i = 0; i < num; i++) {
        if(sr[i] > 0) {
            Z[ct++] = sr[i];
            face_indices.emplace_back(i);
        }
    }
    Z.conservativeResize(ct);
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto mesh_color = mesh->face_property<easy3d::vec3>("f:color");
    ct = 0;
    for(auto f : mesh->faces()) {
        if(sr[f.idx()] > 0) {
            mesh_color[f] = easy3d::vec3(Ct(ct, 0),
                                            Ct(ct, 1), Ct(ct ,2));
            ct++;
        }else {
            mesh_color[f] = easy3d::vec3(0,0,0);
        }
    }
    easy3d::SurfaceMeshIO::save("../sr.ply", mesh);

    VorCorSites viewer("VorCorSites");
    viewer.init(vv, mesh, vvns, vvcs, scvv);
    viewer.run();
    // easy3d::PointCloud* cloud = new easy3d::PointCloud;
    // int idx = 0;
    // for(auto p : vv) {
    //     // if(score[idx] < 80)
    //     cloud->add_vertex(p);
    //     // idx++;
    // }
    // int num = cloud->n_vertices();
    // Eigen::VectorXd Z(cloud->n_vertices());
    // for (int i = 0; i < num; i++) {
    //     Z[i] = vvs[i];
    // }
    //
    // Eigen::MatrixXd Ct;
    // igl::jet(Z, true, Ct);
    // auto coloring = cloud->vertex_property<easy3d::vec3>("v:color", easy3d::vec3(0, 0, 0));
    // for(auto f : cloud->vertices()) {
    //     coloring[f] = easy3d::vec3(Ct(f.idx(), 0),
    //                                     Ct(f.idx(), 1), Ct(f.idx(), 2));
    // }
    // easy3d::PointCloudIO::save("../wuga1.ply", cloud);
    return 0;
}