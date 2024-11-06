//
// Created by xiaowuga on 2024/10/29.
//
#include<easy3d/viewer/viewer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/util/resource.h>

#include <Eigen/Dense>

int ct = 0;

void solve(easy3d::Viewer* viewer, std::string file_path) {
    easy3d::PointCloud* cloud = easy3d::PointCloudIO::load(file_path);
    int num = cloud->n_vertices();
    std::vector<easy3d::vec3> pts = cloud->points();
    Eigen::MatrixXd S(num - 1, 3);
    for(int i = 0; i < num - 1; i++) {
        S.row(i) = Eigen::Vector3d(pts[i].x, pts[i].y, pts[i].z) - Eigen::Vector3d(pts.back().x, pts.back().y, pts.back().z);
    }
    Eigen::MatrixXd M = S.transpose() * S;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigen solver failed!" << std::endl;
    }
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

    easy3d::vec3 aa(0,0,1);
    easy3d::vec3 bb(0,0,-1);
    std::cout << easy3d::cross(aa, bb).norm() << std::endl;
    // easy3d::vec3 aixs2 = easy3d::vec3(eigenvectors.col(2)[0], eigenvectors.col(2)[1], eigenvectors.col(2)[2]).normalize();

    easy3d::vec3 axis1 = easy3d::vec3(eigenvectors.col(0)[0], eigenvectors.col(0)[1], eigenvectors.col(0)[2]).normalize();
    easy3d::vec3 axis2 = easy3d::vec3(eigenvectors.col(1)[0], eigenvectors.col(1)[1], eigenvectors.col(1)[2]).normalize();
    easy3d::vec3 axis3 = easy3d::vec3(eigenvectors.col(2)[0], eigenvectors.col(2)[1], eigenvectors.col(2)[2]).normalize();
    //
    // std::cout << pts.back() << std::endl;
    // std::cout << axis1 << std::endl;
    // std::cout << axis2 << std::endl;
    // std::cout << axis3 << std::endl;

    const std::vector<easy3d::vec3> points = {
        pts.back(),
        pts.back() + axis1,
        pts.back(),
        pts.back() + axis2 * 0.5,
        pts.back(),
        pts.back() + axis3 * 0.1,
    };
    auto lines = new easy3d::LinesDrawable(std::to_string(ct++));
    lines->update_vertex_buffer(points);
    viewer->add_drawable(lines);
    viewer->add_model(cloud);
}
int main() {
    // easy3d::resource::initialize();
    easy3d::Viewer* viewer = new easy3d::Viewer("ASD");
    // solve(viewer, "D:\\code\\defillet\\data2\\19.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\1.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\2.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\3.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\4.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\5.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\6.ply");
    // solve(viewer, "D:\\code\\defillet\\data2\\7.ply");
    for(int i = 0; i < 100; i++) {
        solve(viewer, "D:\\code\\defillet\\data2\\" + std::to_string(i) + ".ply");

        // solve(viewer, "D:\\code\\defillet\\data2\\19.ply");
        // solve(viewer, "D:\\code\\defillet\\data2\\20.ply");
        // solve(viewer, "D:\\code\\defillet\\data2\\20.ply");
    }
    // solve(viewer, "D:\\code\\defillet\\data2\\18.ply");
    viewer->run();

    return 0;
}
