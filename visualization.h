//
// Created by xiaowuga on 1/1/2024.
//

#ifndef DEFILLET_VISUALIZATION_H
#define DEFILLET_VISUALIZATION_H

#include "common.h"
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/texture_manager.h>
#include <easy3d/fileio/surface_mesh_io.h>
//#include <easy3d/util/resource.h>
//#include <easy3d/util/initializer.h>
#include <igl/jet.h>

#include <easy3d/fileio/surface_mesh_io.h>
using namespace easy3d;

namespace RENDER {
    inline void points_visualization(const std::vector<easy3d::vec3> &points) {
        easy3d::PointCloud *cloud = new easy3d::PointCloud;
        int num = points.size();
        for (int i = 0; i < num; i++) {
            cloud->add_vertex(points[i]);
        }

        Viewer viewer("points");
        viewer.add_model(cloud);
        viewer.run();
    }

    inline void points_scalar_field_visualization(const std::vector<easy3d::vec3> &points,
                                                  std::vector<double> &scalar_field) {
        using namespace easy3d;
        easy3d::PointCloud *cloud = new easy3d::PointCloud;
        int num = points.size();
        for (int i = 0; i < num; i++) {
            cloud->add_vertex(points[i]);
        }
        Eigen::VectorXd Z(num);
        for (int i = 0; i < num; i++) {
            Z[i] = scalar_field[i];
        }
        Eigen::MatrixXd Ct;
        igl::jet(Z, true, Ct);
        Viewer viewer("points_scalar_field");
        viewer.add_model(cloud);
        const std::string color_name = "v:color-segments";
        auto coloring = cloud->vertex_property<vec3>(color_name, vec3(0, 0, 0));
        auto drawable = cloud->renderer()->get_points_drawable("vertices");
        for (auto v: cloud->vertices()) {
            coloring[v] = easy3d::vec3(Ct(v.idx(), 0), Ct(v.idx(), 1), Ct(v.idx(), 2));
        }

        drawable->set_impostor_type(PointsDrawable::SPHERE);
        drawable->set_property_coloring(State::VERTEX, color_name);
        drawable->update();
        viewer.update();
        viewer.run();
    }

    inline void  mesh_scalar_field_visualization(const std::vector<easy3d::vec3>& points,
                                                 const std::vector<std::vector<size_t>>& faces,
                                                 std::vector<double>& scalar_field) {
        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(points[i]);
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh->add_triangle(v1, v2, v3);
        }

        Eigen::VectorXd Z(nb_points);
        for (int i = 0; i < nb_points; i++) {
            Z[i] = scalar_field[i];
        }
        Eigen::MatrixXd Ct;
        igl::jet(Z, true, Ct);
        Viewer viewer("mesh_points_scalar_field");
        viewer.add_model(mesh);
        const std::string color_name = "v:color-segments";
        auto coloring = mesh->vertex_property<vec3>(color_name, vec3(0, 0, 0));
        auto drawable_faces = mesh->renderer()->get_triangles_drawable("faces");
        for (auto v: mesh->vertices()) {
            coloring[v] = easy3d::vec3(Ct(v.idx(), 0), Ct(v.idx(), 1), Ct(v.idx(), 2));
        }
        drawable_faces->set_property_coloring(State::VERTEX, color_name);
        drawable_faces->update();

        auto drawable_points = mesh->renderer()->get_points_drawable("vertices");
        drawable_points->set_impostor_type(PointsDrawable::SPHERE);
        drawable_points->set_property_coloring(State::VERTEX, color_name);
        viewer.update();
        viewer.run();
    }

    inline void edge_visualization(std::vector<easy3d::vec3> fillet_points,
                                   std::vector<std::vector<size_t>> fillet_edges) {

        Viewer viewer("edges");
        easy3d::PointsDrawable* pointsDrawable = new easy3d::PointsDrawable("pointsDrawable");
        pointsDrawable->update_vertex_buffer(fillet_points);
        easy3d::LinesDrawable* linesDrawable = new easy3d::LinesDrawable("linesDrawable");
        linesDrawable->update_vertex_buffer(fillet_points);
        std::vector<unsigned int> indices;
        int nb_fillet_edges = fillet_edges.size();
        for(int i = 0; i < nb_fillet_edges; i++) {
            indices.emplace_back(fillet_edges[i][0]);
            indices.emplace_back(fillet_edges[i][1]);
        }
        linesDrawable->update_element_buffer(indices);
        viewer.add_drawable(pointsDrawable);
        viewer.add_drawable(linesDrawable);
        viewer.update();
        viewer.run();
    }

    inline void mesh_visualization(const std::vector<easy3d::vec3>& points,
                                   const std::vector<std::vector<size_t>>& faces) {
        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(points[i]);
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh->add_triangle(v1, v2, v3);
        }
//        mesh->position(easy3d::SurfaceMesh::Vertex(0)) = easy3d::vec3(0,0,0);

        Viewer viewer("meshes");
        viewer.add_model(mesh);
        viewer.run();
    }

    inline void mesh_point_vector_field(const std::vector<easy3d::vec3>& points,
                                           const std::vector<std::vector<std::size_t>>& faces,
                                           const std::vector<easy3d::vec3>& vector_field) {

        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(points[i]);
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh->add_triangle(v1, v2, v3);
        }

        const Box3 &box = mesh->bounding_box();
        float length = norm(box.max_point() - box.min_point()) * 0.01f;
        std::vector<vec3> tmp;
        for (auto ver : mesh->vertices()) {

            const easy3d::vec3 s = mesh->position(ver);
            easy3d::vec3 v = vector_field[ver.idx()];
            const vec3 t = s + v * length;
            tmp.push_back(s);
            tmp.push_back(t);
        }
        Viewer viewer("vector_field");
        viewer.add_model(mesh);
        // Create a drawable for rendering the normal vectors.
        auto drawable = mesh->renderer()->add_lines_drawable("normals");
        // Upload the data to the GPU.
        drawable->update_vertex_buffer(tmp);

        // We will draw the normal vectors in a uniform green color
        drawable->set_uniform_coloring(vec4(0.0f, 1.0f, 0.0f, 1.0f));

        // Set the line width
        drawable->set_line_width(3.0f);

        // Also show the standard "edges"
        mesh->renderer()->get_lines_drawable("edges")->set_visible(true);

        // Run the viewer
        viewer.run();

    }


    inline void binary_mesh_segmentation_visualization(const std::vector<easy3d::vec3>& points,
                                                       const std::vector<std::vector<std::size_t>>& faces,
                                                       const std::vector<double>& label) {
        using namespace easy3d;
        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(points[i]);
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh->add_triangle(v1, v2, v3);
        }

        Viewer viewer("binary_segmentation");
        viewer.add_model(mesh);
        auto drawable = mesh->renderer()->get_triangles_drawable("faces");
        std::string color_name = "f:color";
        auto coloring = mesh->face_property<vec3>(color_name, vec3(0, 0, 0));
        for(auto f : mesh->faces()) {
            int flag = int(label[f.idx()] + 0.1);
            if(flag) {
                coloring[f] = easy3d::vec3(1.0, 0.0, 0.0);
            } else {
                coloring[f] = easy3d::vec3(0.0, 0.0, 1.0);
            }
        }

        const std::string save_file_name = "../bin_seg.ply";
        easy3d::SurfaceMeshIO::save(save_file_name, mesh);
        drawable->set_property_coloring(State::FACE, color_name);
        drawable->update();
        viewer.update();
        viewer.run();
    }

}
//
//
//inline void points_connection_visualization(const std::vector<Point>& points, std::vector<double>& scalar_field,
//                                            std::vector<unsigned int>& indices) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::PointCloud* cloud = new easy3d::PointCloud;
//    int num = points.size();
//    for(int i = 0; i < num; i++) {
//        cloud->add_vertex(easy3d::vec3(points[i].x()
//                , points[i].y(), points[i].z()));
//    }
//
//    Eigen::VectorXd Z(num);
//    for(int i = 0; i < num; i++) {
//        Z[i] = scalar_field[i];
//    }
//    Eigen::MatrixXd Ct;
//    igl::jet(Z, true, Ct);
//    Viewer viewer("scalar_field");
//    viewer.add_model(cloud);
//    const std::string color_name = "v:color-segments";
//    auto coloring = cloud->vertex_property<vec3>(color_name, vec3(0, 0, 0));
//    cv::ColormapTypes selectedColormap = cv::COLORMAP_JET;
//    auto drawable = cloud->renderer()->get_points_drawable("vertices");
//    for(auto v : cloud->vertices()) {
////        cv::Vec3b color = colorImage.at<cv::Vec3b>(v.idx(), 0);
////        cv::Vec3b color = cvCt()
////        OpenCV colorspace is  BGR
////        coloring[v] = easy3d::vec3(color[2]/ 255.0, color[1] / 255.0, color[0] / 255.0);
//        coloring[v] = easy3d::vec3(Ct(v.idx(), 0), Ct(v.idx(), 1), Ct(v.idx(), 2));
//    }
//
//    drawable->set_impostor_type(PointsDrawable::SPHERE);
//    drawable->set_property_coloring(State::VERTEX, color_name);
//    drawable->update();
//    viewer.update();
//    viewer.run();
//}
//
//inline void points_field_range_visualization(const std::vector<Point>& points, std::vector<double>& scalar_field,
//                                      double lower_bound, double upper_bound) {
//    using namespace easy3d;
//
//    easy3d::initialize();
//    easy3d::PointCloud* cloud = new easy3d::PointCloud;
//    int num = points.size();
//
//    for(int i = 0; i < num; i++) {
//        if(lower_bound <= scalar_field[i]
//           && scalar_field[i] <= upper_bound) {
//            cloud->add_vertex(easy3d::vec3(points[i].x(), points[i].y(), points[i].z()));
//        }
//    }
//    Viewer viewer("scalar_field_in_range");
//    viewer.add_model(cloud);
//    viewer.run();
////    delete cloud;
//}
//
//inline void points_visualization(const std::vector<Point>& points) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::PointCloud* cloud = new easy3d::PointCloud;
//    int num = points.size();
//    for(int i = 0; i < num; i++) {
//        cloud->add_vertex(easy3d::vec3(points[i].x()
//                , points[i].y(), points[i].z()));
//    }
//
//    Viewer viewer("scalar_field");
//    viewer.add_model(cloud);
//    viewer.run();
//}
//
//inline void mesh_visualization(const std::vector<Point>& points,
//                        const std::vector<std::vector<std::size_t>>& faces) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
//    int nb_points = points.size();
//    for(int i = 0; i < nb_points; i++) {
//        mesh->add_vertex(easy3d::vec3(points[i].x(),
//                                      points[i].y(), points[i].z()));
//    }
//
//    int nb_faces = faces.size();
//    for(int i = 0; i < nb_faces; i++) {
//        auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
//        auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
//        auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
//        mesh->add_triangle(v1, v2, v3);
//    }
//    Viewer viewer("binary_segmentation");
//    viewer.add_model(mesh);
//
//    viewer.run();
//
//}
//
//inline void binary_mesh_segmentation_visualization(const std::vector<Point>& vertices,
//                                            const std::vector<std::vector<size_t>>& faces,
//                                            const std::vector<bool>& binary) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
//    int nb_vertices = vertices.size();
//    for(int i = 0; i < nb_vertices; i++) {
//        mesh->add_vertex(easy3d::vec3(vertices[i].x(),
//                                      vertices[i].y(), vertices[i].z()));
//    }
//
//    int nb_faces = faces.size();
//    for(int i = 0; i < nb_faces; i++) {
//        auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
//        auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
//        auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
//        mesh->add_triangle(v1, v2, v3);
//    }
//
//    Viewer viewer("binary_segmentation");
//    viewer.add_model(mesh);
//    auto drawable = mesh->renderer()->get_triangles_drawable("faces");
//    std::string color_name = "face_color";
//    auto coloring = mesh->add_face_property<vec3>(color_name, vec3(0, 0, 0));
//    for(auto f : mesh->faces()) {
//        if(binary[f.idx()]) {
//            coloring[f] = easy3d::vec3(1.0, 0.0, 0.0);
//        } else {
//            coloring[f] = easy3d::vec3(0.0, 0.0, 1.0);
//        }
//    }
//
//    drawable->set_property_coloring(State::FACE, color_name);
//    drawable->update();
//    viewer.update();
//    viewer.run();
//}
//
//
//inline void mesh_field_visualization(const std::vector<Point>& points,
//                              const std::vector<std::vector<std::size_t>>& faces,
//                              const std::vector<double>& scalar_field) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
//    int nb_points = points.size();
//    for(int i = 0; i < nb_points; i++) {
//        mesh->add_vertex(easy3d::vec3(points[i].x(),
//                                      points[i].y(), points[i].z()));
//    }
//
//    int nb_faces = faces.size();
//    for(int i = 0; i < nb_faces; i++) {
//        auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
//        auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
//        auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
//        mesh->add_triangle(v1, v2, v3);
//    }
//    Viewer viewer("scalar_field");
//    viewer.add_model(mesh);
//
//    auto drawable = mesh->renderer()->get_triangles_drawable("faces");
//    auto elevation = mesh->add_vertex_property<float>("v:elevation");
//    for (auto v : mesh->vertices())
//        elevation[v] = scalar_field[v.idx()];
//
//
//    drawable->set_scalar_coloring(State::VERTEX, "v:elevation", nullptr, 0.0f, 0.0f);
//
//    const std::string texture_file = resource::directory() + "/colormaps/rainbow.png";
//    Texture *texture = TextureManager::request(texture_file);
//    if (!texture) {
//        return;
//    }
//
//    // Use the texture
//    drawable->set_texture(texture);
//
//    // Run the viewer
//    viewer.run();
//}
//
//
//inline void mesh_face_normals_vector_field(const std::vector<Point>& points,
//                                           const std::vector<std::vector<std::size_t>>& faces,
//                                           const std::vector<Vector_3>& vector_field) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
//    int nb_points = points.size();
//    for(int i = 0; i < nb_points; i++) {
//        mesh->add_vertex(easy3d::vec3(points[i].x(),
//                                      points[i].y(), points[i].z()));
//    }
//
//    int nb_faces = faces.size();
//    for(int i = 0; i < nb_faces; i++) {
//        auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
//        auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
//        auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
//        mesh->add_triangle(v1, v2, v3);
//    }
//
//    const Box3 &box = mesh->bounding_box();
//    float length = norm(box.max_point() - box.min_point()) * 0.05f;
//    std::vector<vec3> tmp;
//    int ct = 0;
//    for (auto f : mesh->faces()) {
////        if(ct++ % 10 != 0) {
////            continue;
////        }
//        vec3 center(0, 0, 0); // face center
//        int count = 0;
//        for (auto v : mesh->vertices(f)) {
//            center += mesh->position(v);
//            ++count;
//        }
//
//        const vec3 s = center / count;
//        vec3 v(vector_field[f.idx()].x(), vector_field[f.idx()].y(), vector_field[f.idx()].z());
//        const vec3 t = s + v * length;
//        tmp.push_back(s);
//        tmp.push_back(t);
//    }
//    Viewer viewer("vector_field");
//    viewer.add_model(mesh);
//    // Create a drawable for rendering the normal vectors.
//    auto drawable = mesh->renderer()->add_lines_drawable("normals");
//    // Upload the data to the GPU.
//    drawable->update_vertex_buffer(tmp);
//
//    // We will draw the normal vectors in a uniform green color
//    drawable->set_uniform_coloring(vec4(0.0f, 1.0f, 0.0f, 1.0f));
//
//    // Set the line width
//    drawable->set_line_width(3.0f);
//
//    // Also show the standard "edges"
//    mesh->renderer()->get_lines_drawable("edges")->set_visible(true);
//
//    // Run the viewer
//    viewer.run();
//
//}
//
//
//inline void points_normals_vector_field(const std::vector<Point>& points,
//                                           const std::vector<Vector_3>& normals) {
//    using namespace easy3d;
//    easy3d::initialize();
//    easy3d::PointCloud* cloud = new easy3d::PointCloud;
//    int nb_points = points.size();
//    for(int i = 0; i < nb_points; i++) {
//        cloud->add_vertex(easy3d::vec3(points[i].x(),
//                                      points[i].y(), points[i].z()));
//    }
//
//
//    const Box3 &box = cloud->bounding_box();
//    float length = norm(box.max_point() - box.min_point()) * 0.02f;
//    std::vector<vec3> tmp;
//    int ct = 0;
//    for (auto v : cloud->vertices()) {
//
//        const vec3 s = cloud->position(v);
//        vec3 d(normals[v.idx()].x(), normals[v.idx()].y(), normals[v.idx()].z());
//        const vec3 t = s + d * length;
//        tmp.push_back(s);
//        tmp.push_back(t);
//    }
//    Viewer viewer("normals_field");
//    viewer.add_model(cloud);
//    // Create a drawable for rendering the normal vectors.
//    auto drawable = new easy3d::LinesDrawable;
//    // Upload the data to the GPU.
//    drawable->update_vertex_buffer(tmp);
//
//    // We will draw the normal vectors in a uniform green color
//    drawable->set_uniform_coloring(vec4(0.0f, 1.0f, 0.0f, 1.0f));
//
//    // Set the line width
//    drawable->set_line_width(3.0f);
//    viewer.add_drawable(drawable);
//    // Also show the standard "edges"
////    mesh->renderer()->get_lines_drawable("edges")->set_visible(true);
//
//    // Run the viewer
//    viewer.run();
//
//}
#endif //DEFILLET_VISUALIZATION_H
