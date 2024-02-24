//
// Created by xiaowuga on 12/17/2023.
//

//#include "io.h"
//#include "data_convert.h"
//#include "utils.h"
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Cell_handle  Cell_handle;
typedef K::Point_3          CGAL_Point;
typedef K::Vector_3         CGAL_Vector_3;
typedef K::Aff_transformation_3 Transformation;
typedef CGAL::Surface_mesh<CGAL_Point> Surface_mesh;
typedef Surface_mesh::Vertex_index Vertex_index;
typedef Surface_mesh::Face_index Face_index;
typedef Surface_mesh::Edge_index Edge_index;
typedef Surface_mesh::Halfedge_index Halfedge_index;