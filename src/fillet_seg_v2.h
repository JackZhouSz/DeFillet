//
// Created by xiaowuga on 2024/10/9.
//

#ifndef FILLET_SEG_V2_H
#define FILLET_SEG_V2_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

class FilletSegV2 {
public:
    FilletSegV2() {
        // fillet segmentation para
        eps_ = 0.03; sigma_ = 10; min_score_ = 0.5;
        // sor para
        std_ratio_ = 1.0; nb_neighbors_ = 30; num_sor_iter_ = 3;

    }
    void FilletSegV2::seg(){}

    void voronoi3d(const std::vector<easy3d::vec3>& sites
                 , std::vector<easy3d::vec3>& vor_vertices
                 , std::vector<std::vector<int>>& sites_nearby_vertices
                 , std::vector<std::vector<int>>& vertices_nearby_sites);
    void sor(const std::vector<easy3d::vec3>& points
           , std::vector<bool>& labels);
private:
    easy3d::SurfaceMesh* mesh_;

    easy3d::Box3 box;
    int num_sor_iter_;
    double std_ratio_;
    int nb_neighbors_;
    double eps_;
    double sigma_;
    double radius_;
    double min_score_;

    double scoring_time_;
};



#endif //FILLET_SEG_V2_H
