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
#include <easy3d/algo/surface_mesh_geometry.h>
#include <kernel.h>
#include <easy3d/core/random.h>

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

    // std::vector<float> field2 = field;
    // float maxx = *std::max_element(field.begin(), field.end());
    // float minn = *std::min_element(field.begin(), field.end());
    // for(int i = 0; i < num; i++) {
    //     float val = (field2[i] - minn) / (maxx - minn);
    //     val = field[i];
    //     // val = gaussian_kernel(val, 0.1);
    //     // if(val < 1e-3) {
    //     //     val = 1e-3;
    //     // }
    //     // if(val > 1.0 - 1e-3) {
    //     //     val = 0.99;
    //     // }
    //     // if(val < 0.9) {
    //     //     val *= 0.5;
    //     // }
    //     field2[i] = val;
    //
    // }
    Eigen::VectorXd Z(num);

    for (int i = 0; i < num; i++) {
        Z[i] = field[i];
    }
    // Z.conservativeResize(ct);
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto color = cloud->vertex_property<easy3d::vec3>("v:color");
    for(auto v : cloud->vertices()) {

        // std::cout << Ct(v.idx(), 0) <<std::endl;
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
        if(field[f.idx()] > 0.5) {
            color[f] = easy3d::vec3(0.0,0.0, 1.0);
        }
        else {

            color[f] = easy3d::vec3(1.0,0.0, 0.0);
        }
        // color[f] = easy3d::vec3(Ct(f.idx(), 0),Ct(f.idx(), 1), Ct(f.idx() ,2));
    }

    easy3d::SurfaceMeshIO::save(path, mesh);
}

static void save_mesh_field_with_mtl(easy3d::SurfaceMesh* mesh
                   , std::vector<float>& field, const std::string& mesh_path, const std::string& mtl_path) {
    std::vector<float> field2 = field;
    float maxx = *std::max_element(field.begin(), field.end());
    float minn = *std::min_element(field.begin(), field.end());
    int num = field.size();
    for(int i = 0; i < num; i++) {
        // std::cout <<
        // if(field2[i] > 0.999)
        //     field2[i] = 0;
        float val = (field2[i] - minn) / (maxx - minn);
        // float val = (field2[i] - minn) / (maxx - minn);
        //  {
        //     val = 0.5;
        // }
        // else
        // float val;
        // val = field[i];
       // float val = gaussian_kernel(field2[i], 0.3);
        if(val < 1e-3) {
            val = 1e-3;
        }
        if(val > 1.0 - 1e-3) {
            val = 0.999;
        }
        field2[i] = val;
    }
    std::ofstream mtlFile(mtl_path);
    mtlFile << "newmtl texture_material\n";
    mtlFile << "map_Kd turbo.png\n";
    mtlFile.close();

    std::ofstream objFile(mesh_path);
    objFile << "mtllib " << mtl_path << "\n";
    objFile << "usemtl texture_material\n";
    int ct = 0;
    for(auto f : mesh->faces()) {
        std::vector<int> tmp;
        // if(field2[f.idx()] < 0.9) {
            for(auto v : mesh->vertices(f)) {
                auto p = mesh->position(v);
                objFile << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
                objFile << "vt " << field2[f.idx()] << " " << 1e-3 << "\n";
                tmp.emplace_back(++ct);
            }
            objFile << "f";
            for(auto id : tmp) {
                objFile << " "<< id << "/" << id;
            }
            objFile << "\n";
        // }
    }
    objFile.close();
}

static void save_mesh_field_with_mtl1(easy3d::SurfaceMesh* mesh
                   , std::vector<float>& field, const std::string& mesh_path, const std::string& mtl_path) {
    std::vector<float> field2 = field;
    float maxx = *std::max_element(field.begin(), field.end());
    float minn = *std::min_element(field.begin(), field.end());
    int num = field.size();
    for(int i = 0; i < num; i++) {
        float val = (field2[i] - minn) / (maxx - minn);
        if(val < 1e-3) {
            val = 1e-3;
        }
        if(val > 1.0 - 1e-3) {
            val = 0.99;
        }
        field2[i] = val;
    }
    std::ofstream mtlFile(mtl_path);
    mtlFile << "newmtl texture_material\n";
    mtlFile << "map_Kd turbo.png\n";
    mtlFile.close();

    std::ofstream objFile(mesh_path);
    objFile << "mtllib " << mtl_path << "\n";
    objFile << "usemtl texture_material\n";
    int ct = 0;
    for(auto v : mesh->vertices()) {
        easy3d::vec3 p = mesh->position(v);
        objFile << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
        objFile << "vt " << field2[v.idx()] << " " << 1e-3 << "\n";
        // tmp.emplace_back(++ct);
    }
    for(auto f : mesh->faces()) {
        objFile << "f";
        for(auto v : mesh->vertices(f)) {

            objFile << " " << v.idx() + 1 << "/" << v.idx() + 1;
            // tmp.emplace_back(++ct);
        }
        objFile << "\n";
    }
    objFile.close();
}


static void save_mesh_normals(easy3d::SurfaceMesh* mesh
                            , std::vector<easy3d::vec3> normals
                            , const std::string& path)  {
    std::ofstream outFile(path);
    double len = mesh->bounding_box().diagonal_length() * 0.02;
    int ct = 1;
    for(auto f : mesh->faces()) {
        auto p1 =  easy3d::geom::centroid(mesh, f);
        auto p2 =  p1 + normals[f.idx()] * len;
        outFile << "v " << p1.x << " " << p1.y << " " << p1.z << std::endl;
        outFile << "v " << p2.x << " " << p2.y << " " << p2.z << std::endl;
        outFile << "l " << ct << " " << ct + 1 << std::endl;
        ct += 2;
    }
    outFile.close();
}

static void save_lines(std::vector<easy3d::vec3> lines
                            , const std::string& path)  {
    std::ofstream outFile(path);
    int ct = 1;
    int num = lines.size();
    for(int i = 0; i < num; i += 2) {
        outFile << "v " << lines[i].x << " " << lines[i].y << " " << lines[i].z << std::endl;
        outFile << "v " << lines[i + 1].x << " " << lines[i + 1].y << " " << lines[i + 1].z << std::endl;
        outFile << "l " << ct << " " << ct + 1 << std::endl;
        // break;
        ct += 2;
    }
    outFile.close();
}


static double dihedral_angle(const easy3d::vec3& n1
                                 , const easy3d::vec3& n2
                                 , bool rad = false) {

    // double radians = abs(acos(dot(n1, n2)));
    double radians = easy3d::geom::angle(n1, n2);
    if(rad) {
        return radians;
    } else {
        double degrees = radians * 180.0 / M_PI;
        return degrees;
    }
}


static void save_components(const easy3d::SurfaceMesh* mesh,
                            const std::vector<easy3d::SurfaceMesh*>components,
                            const std::string path) {

    easy3d::SurfaceMesh* out_mesh = new easy3d::SurfaceMesh(*mesh);

    auto color_prop = out_mesh->face_property<easy3d::vec3>("f:color");

    for(size_t i = 0; i < components.size(); i++) {

        auto color = easy3d::random_color();

        auto component = components[i];
        auto original_face_index = component->face_property<int>("f:original_index");

        for(auto face : component->faces()) {
            int index = original_face_index[face];
            auto f = easy3d::SurfaceMesh::Face(index);
            color_prop[f] = color;
        }
    }

    easy3d::SurfaceMeshIO::save(path, out_mesh);
}

#endif //IO_H
