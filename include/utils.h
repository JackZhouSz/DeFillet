//
// Created by xiaowuga on 2025/8/15.
//

#ifndef UTILS_H
#define UTILS_H

#include <easy3d/core/point_cloud.h>
#include <easy3d/core/surface_mesh.h>




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
    float angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2);



    void save_components(const easy3d::SurfaceMesh* mesh,
                                const std::vector<easy3d::SurfaceMesh*>components,
                                const std::string path);


    easy3d::SurfaceMesh* split_component(const easy3d::SurfaceMesh* mesh,
                         easy3d::SurfaceMesh::FaceProperty<int>& component_labels,
                         int label);


    void save_fillet_regions(const easy3d::SurfaceMesh* mesh,
                             const easy3d::SurfaceMesh::FaceProperty<int>& fillet_label,
                             const std::string path);


}

#endif //UTILS_H
