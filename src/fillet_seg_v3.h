//
// Created by xiaowuga on 2024/10/31.
//

#ifndef FILLET_SEG_V3_H
#define FILLET_SEG_V3_H

#include "easy3d/core/surface_mesh.h"
#include "easy3d/core/point_cloud.h"
#include <Eigen/Dense>

class FilletSegV3 {

private:
    void sor();

public:
    FilletSegV3(easy3d::SurfaceMesh* mesh);
    void voronoi3d(const std::vector<easy3d::vec3>& sites
             , std::vector<easy3d::vec3>& vor_vertices
             , std::vector<float>& vor_vertices_radius
             , std::vector<std::vector<int>>& sites_nearby_vertices
             , std::vector<std::vector<int>>& vertices_nearby_sites);
    void crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
             , std::vector<easy3d::SurfaceMesh::Face>& patch);

    void solve();

    void generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv,std::vector<float>& vvr, std::vector<std::vector<int>>& vvns, bool sor_flag = true);
    void generate_global_voronoi_vertices(std::vector<easy3d::vec3>& vv, std::vector<std::vector<int>>& vvns);

    void search_voronoi_vertices_corresponding_sites(std::vector<easy3d::vec3>& vv, std::vector<std::vector<int>>& vvns, std::vector<std::vector<int>>& vvcs);

    void search_sites_corresponding_voronoi_vertices(std::vector<std::vector<int>>& vvcs, std::vector<std::vector<int>>& scvv);

    void compute_voronoi_vertices_score(std::vector<easy3d::vec3>& vv,std::vector<float>& vvr, std::vector<float>& vvs);

    void esimate_sites_radius(std::vector<easy3d::vec3>& sites, std::vector<easy3d::vec3>& vv
                            ,std::vector<float>& vvr, std::vector<std::vector<int>>& scvv
                            ,std::vector<float>& sr);

    void sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels);

    void kmeans(const Eigen::MatrixXd& data, int num_clusters, Eigen::MatrixXd& centroids, Eigen::VectorXi& labels);

    bool find_cluster_center(const std::vector<easy3d::vec3>& data, easy3d::vec3& centroid);

    bool find_centroid_via_kmeans(std::vector<easy3d::vec3>& points, std::vector<float>& radius, easy3d::vec3& centroid);

public:
    easy3d::Box3 box_;

    std::vector<double> sites_area_;
    std::vector<easy3d::vec3> sites_normals_;


    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vv_; //voronoi vertices;
    std::vector<float> vvr_; // voronoi vertices radius
    // std::vector<std::vector<int>> vvns_; // voronoi vertices nearby sites
    std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites
    std::vector<std::vector<int>> scvv_; // sites corresponding voronoi vertices
    std::vector<float> vvs_; // voronoi vertices score
    std::vector<float> sr_; // sites radius
    std::vector<std::pair<float, float>> sar_; // span_and_radius;

    easy3d::SurfaceMesh* mesh_;

    double eps_;
    double sigma_;
    double radius_;
    double angle_thr;

    int nb_neighbors_;
};

#endif //FILLET_SEG_V3_H
