//
// Created by xiaowuga on 2024/12/15.
//

#ifndef DEFILLET_V2_H
#define DEFILLET_V2_H

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

class DeFilletv2 {
public:
    easy3d::SurfaceMesh* mesh_;
    easy3d::SurfaceMesh* fillet_mesh_;
    easy3d::SurfaceMesh* non_fillet_mesh_;
};



#endif //DEFILLET_V2_H
