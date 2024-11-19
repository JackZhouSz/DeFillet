//
// Created by xiaowuga on 2024/11/14.
//

#ifndef FILLET_SEG_V5_H
#define FILLET_SEG_V5_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"
#include <Eigen/Dense>

class FilletSegV5 {
public:
    FilletSegV5(easy3d::SurfaceMesh* mesh);
    void run();

    // core functions
    void generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv
                                       , std::vector<float>& vvr
                                       , std::vector<std::vector<int>>& vvns
                                       , bool sor_flag = true
                                       , bool radius_thr = true);

    void voronoi_density_drift(std::vector<easy3d::vec3>& vv
                             , std::vector<float>& vvr
                             , std::vector<std::vector<int>>& vvns
                             , std::vector<std::vector<int>>& scvv);



    void farthest_point_sampling(int nb_samples, std::vector<int>& indices, float& max_gap);




    // utils
    void sor(const std::vector<easy3d::vec3>& points
           , std::vector<bool>& labels
           , int nb_neighbors, int num_sor_iter, float std_ratio);
    void crop_local_patch(int f, float max_gap
                        , std::vector<int>& indices);
    void voronoi3d(const std::vector<easy3d::vec3>& s
                 , std::vector<easy3d::vec3>& v
                 , std::vector<float>& vvr
                 , std::vector<std::vector<int>>& snv
                 , std::vector<std::vector<int>>& vns);
    double dihedral_angle(easy3d::SurfaceMesh* mesh
                  , easy3d::SurfaceMesh::Face f1
                  , easy3d::SurfaceMesh::Face f2
                  , bool rad = false);
    int calculate_samples_num(float radius);

    void save_point_set(std::vector<easy3d::vec3>& points, std::string& path);
    //test
public:
    easy3d::Box3 box_;
    easy3d::SurfaceMesh* mesh_;

    std::vector<double> sa_; //site area
    std::vector<easy3d::vec3> s_;   // sites
    std::vector<easy3d::vec3> vv_; //voronoi vertices;
    std::vector<float> vvr_; // voronoi vertices radius
    std::vector<std::vector<int>> vvns_; // voronoi vertices neighboring sites
    std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites
    std::vector<std::vector<int>> scvv_; // sites corresponding voronoi vertices

    double eps_;
    double sigma_;
    float radius_;
    double angle_thr_;

    int nb_neighbors_;
    float std_radio_;
    int num_sor_iter_;
};



#endif //FILLET_SEG_V5_H
