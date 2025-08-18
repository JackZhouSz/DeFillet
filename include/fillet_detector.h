//
// Created by xiaowuga on 2025/8/16.
//

#ifndef DEFILLET_FILET_DETECTOR_H
#define DEFILLET_FILET_DETECTOR_H

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>

using namespace easy3d;
using namespace std;


namespace DeFillet {

    class FilletDetectorParameters {

    public:

        float epsilon, radius_thr, angle_thr, sigma, lamdba;

        int num_patches, num_neighbors, num_smooth_iter;
    };

    class VoronoiVertices {
    public:
        vec3 pos;

        float radius;

        std::vector<SurfaceMesh::Face> neigh_sites; // neighboring sites

        std::vector<SurfaceMesh::Face> corr_sites; //corresponding sites

        bool flag;
    };


    class FilletDetector {

    public:

        FilletDetector(SurfaceMesh *mesh, FilletDetectorParameters& parameters);

        void apply();

        void generate_voronoi_vertices();

        void filter_voronoi_vertices();

        void density_driven_voronoi_drift();
        //
        // void compute_voronoi_vertices_density_field();
        //
        // void rolling_ball_trajectory_transform();
        //
        // void graph_cut();

    public: // setor or getor functions

        PointCloud* voronoi_vertices() const;

    private: //utils functions

        float farthest_point_sampling(int num_samples,  std::vector<SurfaceMesh::Face>& indices);

        void crop_local_patch(SurfaceMesh::Face cface, float max_gap, std::vector<SurfaceMesh::Face>& patch);

        int check_osculation_condition(VoronoiVertices& vor, float thickness) const;

        int count_boundaries_components(std::vector<SurfaceMesh::Edge>& boundaries) const;

        std::vector<vec3> face_centroids(std::vector<SurfaceMesh::Face>& faces) const;


    public:

        // std::vector<easy3d::vec3> vv_; //voronoi vertices
        // std::vector<float> vvr_; // voronoi vertices radius
        //
        // std::vector<std::vector<int>> vvns_; // voronoi vertices neighboring sites
        // std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites

        SurfaceMesh* mesh_;

        Box3 box_;

        FilletDetectorParameters parameters_;

    private:

        SurfaceMesh::FaceProperty<vec3> face_normals_;

        SurfaceMesh::EdgeProperty<float> dihedral_angle_;

        SurfaceMesh::EdgeProperty<float> centroid_distance_;

        std::vector<VoronoiVertices> voronoi_vertices_;


    };
}


#endif //DEFILLET_FILET_DETECTOR_H
