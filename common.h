//
// Created by xiaowuga on 12/17/2023.
//
#include "happly.h"
#include "io.h"
#include "data_convert.h"
#include "utils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef K::Point_3          Point;
typedef K::Vector_3         Vector_3;
typedef K::Aff_transformation_3 Transformation;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef Surface_mesh::Vertex_index Vertex_index;
typedef Surface_mesh::Face_index Face_index;
typedef Surface_mesh::Edge_index Edge_index;
typedef Surface_mesh::Halfedge_index Halfedge_index;