//
// Created by xiaowuga on 2025/8/15.
//

#ifndef UTILS_H
#define UTILS_H

#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/core/random.h>




namespace DeFillet {


    /**
    * @brief Compute the (unsigned) angle between two 3D vectors in degrees.
    *
    * Uses atan2(||n1 × n2||, n1 · n2), which is generally more numerically stable
    * than acos of the normalized dot product.
    *
    * @param n1 First vector (e.g., a face normal).
    * @param n2 Second vector (e.g., a face normal).
    * @return Angle in degrees, in the range [0, 180].
    */
    static float angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2) {
        const double dot = easy3d::dot(n1, n2);
        const double cross_norm = easy3d::cross(n1, n2).norm();

        return std::atan2(cross_norm, dot) * 180.0 / M_PI;
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

    static void save_fillet_regions(const easy3d::SurfaceMesh* mesh,
                                    const easy3d::SurfaceMesh::FaceProperty<int>& fillet_label,
                                    const std::string path) {
        easy3d::SurfaceMesh* out_mesh = new easy3d::SurfaceMesh(*mesh);

        auto color_prop = out_mesh->face_property<easy3d::vec3>("f:color");
        easy3d::vec3 fillet_color = easy3d::vec3(1.0, 0,0);
        easy3d::vec3 non_fillet_color = easy3d::vec3(0, 0,1.0);

        for(auto face : out_mesh->faces()) {
            int label = fillet_label[face];
            color_prop[face] = label != 0 ? fillet_color : non_fillet_color;
        }
        easy3d::SurfaceMeshIO::save(path, out_mesh);

    }


}

#endif //UTILS_H
