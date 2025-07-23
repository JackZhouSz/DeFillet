//
// Created by xiaowuga on 2024/12/9.
//

#ifndef FILLET_SEG_V8_H
#define FILLET_SEG_V8_H

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>

class FilletSegV8 {
public:
    FilletSegV8(easy3d::SurfaceMesh* mesh
               , float eps = 0.02
               , int num_samples = 500
               , float radius_thr = 0.07
               , float angle_thr = 40
               , bool local_voronoi = true
               , bool radius_filter = true
               , int num_neighbors = 30
               , float h = 0.1
               , int num_smooth_iter = 10
               , float lamdba = 0.5);

    void run();

    void generate_voronoi_vertices(std::vector<easy3d::vec3>& vv
                             , std::vector<float>& vvr
                             , std::vector<std::vector<int>>& vvns
                             , std::vector<std::vector<int>>& vvcs
                             , std::vector<std::vector<int>>& scvv
                             , std::vector<float>& vv_density);

    void generate_voronoi_vertices_via_global_sites(std::vector<easy3d::vec3>& vv
                                                  , std::vector<float>& vvr
                                                  , std::vector<std::vector<int>>& vvns);

    void generate_voronoi_vertices_via_local_patch(std::vector<easy3d::vec3>& vv
                                                 , std::vector<float>& vvr
                                                 , std::vector<std::vector<int>>& vvns);


    void marching_axis_transform(std::vector<easy3d::vec3>& mavv
                               , std::vector<std::vector<int>>& mavvns
                               , std::vector<std::vector<int>>& mavvcs
                               , std::vector<float>& mavvr
                               , std::vector<int>& scma);

    void run_graph_cut(std::vector<float>& score, std::vector<int>& fillet_labels);


    void post_processing();
    //utils
    void crop_local_patch(int sid
                        , float max_gap
                        , std::vector<int>& indices);
    bool search_vvcs(easy3d::vec3& vv,  float r, std::vector<int>& vvns, std::vector<int>& vvcs);


    void farthest_point_sampling(int num_samples
                               , std::vector<easy3d::vec3>& s
                               , std::vector<std::vector<int>>& sg
                               , std::vector<int>& indices);

    void sor(const std::vector<easy3d::vec3>& points
           , std::vector<bool>& labels
           , int nb_neighbors, int num_sor_iter, float std_ratio);

    void voronoi3d(const std::vector<easy3d::vec3>& s
                 , std::vector<easy3d::vec3>& v
                 , std::vector<float>& vvr
                 , std::vector<std::vector<int>>& snv
                 , std::vector<std::vector<int>>& vns);

    // double dihedral_angle(const easy3d::vec3& n1
    //                     , const easy3d::vec3& n2
    //                     , bool rad = false);

public:
    easy3d::SurfaceMesh* mesh_;
    std::vector<easy3d::vec3> s_; // sites
    std::vector<easy3d::vec3> sn_; // sites normals
    std::vector<std::vector<int>> sg_; // sites graph
    std::vector<float> ss_; // site score
    std::vector<int> fillet_labels_; // 0 for fillet, 1 for non-fillet
    std::vector<std::pair<int,int>> sge_; // sites graph edges
    std::vector<float> sgel_; // sites graph edges length
    std::vector<easy3d::vec3> vv_; //voronoi vertices;
    std::vector<float> vvr_; // voronoi vertices radius
    std::vector<float> vv_density_; //voronoi vertices density
    std::vector<std::vector<int>> vvns_; // voronoi vertices neighboring sites
    std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites
    std::vector<std::vector<int>> scvv_; // sites corresponding voronoi vertices
    std::vector<float> vvs_; // voronoi vertices score
    std::vector<easy3d::vec3> mavv_; //marching axis voronoi vertices
    std::vector<float> mavvr_; // marching axis voronoi vertices radius
    std::vector<std::vector<int>> mavvns_; // marching axis voronoi vertices neighboring sites
    std::vector<std::vector<int>> mavvcs_; // marching axis voronoi vertices corresponding sites
    std::vector<int> scma_; // site corresponding marching axis

    float eps_;    // tolerance to noise
    float radius_thr_;  // radius threshold
    float angle_thr_; // angle threshold

    int num_samples_; // samples number of FPS

    // sor para
    int num_neighbors_; // k-neighbors parameter of kd-tree
    float sor_std_radio_;   // standard deviation
    int num_sor_iter_;  // iteration number of sor


    float h_;   // kernel bandwidth
    int num_smooth_iter_; // iteration number of smmoth
    float lamdba_; // graph cut lamdba

    // time
    float gvv_time_; //generate voronoi vertices timme
    float mat_time_; //voronoi dennity drift time
    float gcp_time_; //graph cut time
    float seg_tot_time_; //segmentation total time

    bool local_voronoi_; // generate Voronoi vertices using global or local strategy (true for local, and false for global)
    bool radius_filter_; // use radius threshhold or not
    bool sor_filter_; // use sor filter or not


};



#endif //FILLET_SEG_V8_H
