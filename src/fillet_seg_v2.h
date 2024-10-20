//
// Created by xiaowuga on 2024/10/9.
//

#ifndef FILLET_SEG_V2_H
#define FILLET_SEG_V2_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"

class FilletSegV2 {
public:
    FilletSegV2(easy3d::SurfaceMesh* mesh);
    const std::vector<easy3d::vec3>& get_vor_vertices() {return vor_vertices_;}
    const std::vector<float>& get_vor_scores() {return vor_vertices_score_;}
    void voronoi3d(const std::vector<easy3d::vec3>& sites
                 , std::vector<easy3d::vec3>& vor_vertices
                 , std::vector<std::vector<int>>& sites_nearby_vertices
                 , std::vector<std::vector<int>>& vertices_nearby_sites);
    void sor(const std::vector<easy3d::vec3>& points
           , std::vector<bool>& labels);

    void seg();
    void crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
                        , std::vector<easy3d::SurfaceMesh::Face>& patch);
    void probability_of_sites();
    float probability_of_vertex(easy3d::vec3& vertex
                , std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);
    void density(std::vector<easy3d::vec3>& points, std::vector<float>& avg_distances);
    void decompose();
    void deduplication(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels);

private:
    easy3d::Box3 box_;
    easy3d::SurfaceMesh::FaceProperty<float> face_area;
    easy3d::SurfaceMesh::FaceProperty<easy3d::vec3> face_normals;
    easy3d::SurfaceMesh::FaceProperty<float> face_score;
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> mp_;
    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vor_vertices_;
    std::vector<std::vector<int>> site_cconfilt_;
    std::vector<float> vor_vertices_score_;
    easy3d::SurfaceMesh* mesh_;



    // sor para
    int num_sor_iter_;
    double std_ratio_;
    int nb_neighbors_;


    double eps_;
    double sigma_;
    double radius_;
    double angle_thr;
    double min_score_;

    double scoring_time_;
};



#endif //FILLET_SEG_V2_H
