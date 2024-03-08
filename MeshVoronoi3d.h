//
// Created by xiaowuga on 3/6/2024.
//

#ifndef DEFILLET_MESHVORONOI3D_H
#define DEFILLET_MESHVORONOI3D_H

#include "common.h"
#include <easy3d/core/surface_mesh.h>

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
