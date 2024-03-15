//
// Created by xiaowuga on 3/6/2024.
//

#ifndef DEFILLET_MESHVORONOI3D_H
#define DEFILLET_MESHVORONOI3D_H

#include "easy3d/core/surface_mesh.h"

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


class MeshVoronoi3d {
public:
    MeshVoronoi3d(easy3d::SurfaceMesh* mesh);

    const std::vector<easy3d::vec3>& get_sites() const { return sites_;}
    const std::vector<easy3d::vec3>& get_vertices() const { return vertices_;}
    const std::vector<std::vector<int>>& get_regions() const { return regions_;}
    const std::vector<std::vector<int>>& get_poles() const { return poles_;}
    const std::vector<int>& get_final_pole() const {return final_pole_;}
    const std::vector<std::vector<int>>& get_vertices2sites() const {return vertices2sites_;}
private:
    easy3d::SurfaceMesh* mesh_;
    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vertices_;
    std::vector<std::vector<int>> regions_;
    std::vector<std::vector<int>> poles_;
    std::vector<int> final_pole_;
    std::vector<std::vector<int>> vertices2sites_;

    Triangulation T_;
};


#endif //DEFILLET_MESHVORONOI3D_H
