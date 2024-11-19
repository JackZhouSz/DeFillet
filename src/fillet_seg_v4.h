//
// Created by xiaowuga on 2024/11/6.
//

#ifndef FILLET_SEG_V4_H
#define FILLET_SEG_V4_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"
#include <Eigen/Dense>

class FilletSegV4 {
public:
    FilletSegV4(easy3d::SurfaceMesh* mesh);

    void crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
         , std::vector<easy3d::SurfaceMesh::Face>& patch);

    void generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv
                                       , std::vector<float>& vvr
                                       , std::vector<std::vector<int>>& vvns
                                       , std::vector<std::vector<int>>& vvcs
                                       , bool sor_flag = true
                                       , bool radius_thr = true);



    void voronoi_vertices_drift(std::vector<easy3d::vec3>& vv
                              , std::vector<float>& vvr
                              , std::vector<std::vector<int>>& vvcs);

    void seg();
    void clustering(std::vector<easy3d::vec3>& vv
              , std::vector<float>& vvr
              , std::vector<std::vector<int>>& vvcs
              , std::vector<easy3d::vec3>& cc
              , std::vector<float>& cr
              , std::vector<std::vector<int>>& ccs
              , std::vector<std::vector<int>>& ccvv
              , std::vector<float>& ccvvn
              , std::vector<float>& ccsn
              , int num_cluster);
private:
    void sor(const std::vector<easy3d::vec3>& points
        , std::vector<bool>& labels
        , int nb_neighbors, int num_sor_iter, float std_ratio);
    void voronoi3d(const std::vector<easy3d::vec3>& s
                 , std::vector<easy3d::vec3>& vv
                 , std::vector<float>& vvr
                 , std::vector<std::vector<int>>& snv
                 , std::vector<std::vector<int>>& vns);
    double dihedral_angle(easy3d::SurfaceMesh* mesh
                      , easy3d::SurfaceMesh::Face f1
                      , easy3d::SurfaceMesh::Face f2
                      , bool rad = false);

    bool search_voronoi_vertex_corresponding_sites(easy3d::vec3& vv
                                                 , float r
                                                 , std::vector<int>& vvns
                                                 , std::vector<int>& vvcs);

    void fit_sphere(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius, int max_iter = 100);

    void tetrahedron_circumcenter(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius);
public:
    easy3d::Box3 box_;

    std::vector<double> sites_area_;
    std::vector<easy3d::vec3> sites_normals_;


    std::vector<easy3d::vec3> s_;   // sites
    std::vector<easy3d::vec3> vv_; //voronoi vertices;
    std::vector<float> vvr_; // voronoi vertices radius
    std::vector<std::vector<int>> vvns_; // voronoi vertices neighboring sites
    std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites
    std::vector<easy3d::vec3> cc_;  // cluster centroids
    std::vector<float> cr_; // centroids radius
    std::vector<std::vector<int>> ccs_; // centroids corresponding sites
    std::vector<std::vector<int>> ccvv_; // centroids corresponding voronoi vertices
    std::vector<float> ccvvn_; // centroids corresponding voronoi vertices num
    std::vector<float> ccsn_; // centroids corresponding sites num
    std::vector<std::vector<int>> scvv_; // sites corresponding voronoi vertices
    std::vector<float> vvs_; // voronoi vertices score
    std::vector<float> sr_; // sites radius

    easy3d::SurfaceMesh* mesh_;

    double eps_;
    double sigma_;
    double radius_;
    double angle_thr_;

    int nb_neighbors_;
    float std_radio_;
    int num_sor_iter_;

    int num_cluster_;

};


#endif //FILLET_SEG_V4_H
