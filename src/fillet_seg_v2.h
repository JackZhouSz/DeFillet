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
    const std::vector<double>& get_vor_scores() {return vor_vertices_score_;}
    void voronoi3d(const std::vector<easy3d::vec3>& sites
                 , std::vector<easy3d::vec3>& vor_vertices
                 , std::vector<std::vector<int>>& sites_nearby_vertices
                 , std::vector<std::vector<int>>& vertices_nearby_sites);
    void sor(const std::vector<easy3d::vec3>& points
           , std::vector<bool>& labels);

    void seg();
    void crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
                        , std::vector<easy3d::SurfaceMesh::Face>& patch);

    void eigen_value(bool sor_flag = true);
    void probability_of_sites();
    void probability_of_sites2();
    void sites_probablity(bool sor_flag = true);

    double transition_angle(easy3d::vec3& vertex,
        std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);

    float probability_of_vertex(easy3d::vec3& vertex
                , std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);

    bool vor_vertex_aixs(easy3d::vec3& vertex, std::vector<int>& vertex_nearby_sites_indices, easy3d::vec3& axis);

    float FilletSegV2::probability_of_vertex_eigen_value(easy3d::vec3& vertex
                , std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);


    float FilletSegV2::probability_of_vertex_area(easy3d::vec3& vertex
            , std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);

    float FilletSegV2::probability_of_vertex_normal(easy3d::vec3& vertex
        , std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites);

    void density(std::vector<easy3d::vec3>& points, std::vector<double>& avg_distances);
    void decompose();
    void deduplication(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels);

public:
    easy3d::Box3 box_;
    easy3d::SurfaceMesh::FaceProperty<float> face_area;
    easy3d::SurfaceMesh::FaceProperty<easy3d::vec3> face_normals;
    easy3d::SurfaceMesh::FaceProperty<float> face_score;
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> mp_;


    std::vector<easy3d::vec3> total_vor_vertices_;
    std::vector<std::vector<int>> total_vor_nearby_sites_;

    std::vector<double> vor_radius_;
    std::vector<double> vor_transition_angle_;


    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vor_vertices_;
    std::vector<std::vector<int>> site_cconfilt_;
    std::vector<double> vor_vertices_score_;
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
    int dsa;


    std::vector<std::vector<int>> vor_cor_sites_;
};



#endif //FILLET_SEG_V2_H
