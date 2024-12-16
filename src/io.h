//
// Created by xiaowuga on 2024/12/1.
//

#ifndef IO_H
#define IO_H
#include <igl/jet.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>

static void save_point_set(std::vector<easy3d::vec3>& points,const std::string& path) {
    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(auto p : points) {
        cloud->add_vertex(p);
    }
    easy3d::PointCloudIO::save(path, cloud);
}
static void save_point_field(std::vector<easy3d::vec3>& points
                   , std::vector<float>& field, const std::string& path)  {
    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(auto p : points) {
        cloud->add_vertex(p);
    }
    int num = field.size();
    Eigen::VectorXd Z(num);

    for (int i = 0; i < num; i++) {
        Z[i] = field[i];
    }
    // Z.conservativeResize(ct);
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto color = cloud->vertex_property<easy3d::vec3>("v:color");
    for(auto v : cloud->vertices()) {
        color[v] = easy3d::vec3(Ct(v.idx(), 0),Ct(v.idx(), 1), Ct(v.idx() ,2));
    }
    easy3d::PointCloudIO::save(path, cloud);
}

static void save_mesh_field(easy3d::SurfaceMesh* mesh
                   , std::vector<float>& field, const std::string& path)  {

    auto color = mesh->face_property<easy3d::vec3>("f:color");
    int num = field.size();
    Eigen::VectorXd Z(num);
    for (int i = 0; i < num; i++) {
        Z[i] = field[i];
    }
    // Z.conservativeResize(ct);
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    for(auto f : mesh->faces()) {
        color[f] = easy3d::vec3(Ct(f.idx(), 0),Ct(f.idx(), 1), Ct(f.idx() ,2));
    }

    easy3d::SurfaceMeshIO::save( path, mesh);
}

#endif //IO_H
