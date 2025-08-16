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
        easy3d::vec3 pos;

        float radius;

        std::vector<int> neigh_sites; // neighboring sites

        std::vector<int> corr_sites;

    };


    class FilletDetector {

    public:

        FilletDetector(SurfaceMesh *mesh, FilletDetectorParameters& parameters);

        void apply();


    private:

        float farthest_point_sampling(int num_samples,  std::vector<SurfaceMesh::Face>& indices);

        void crop_local_patch(SurfaceMesh::Face cface, float max_gap, std::vector<SurfaceMesh::Face>& patch);

        std::vector<vec3> face_centroids(std::vector<SurfaceMesh::Face>& faces) const;


    public:

        std::vector<easy3d::vec3> vv_; //voronoi vertices
        std::vector<float> vvr_; // voronoi vertices radius

        std::vector<std::vector<int>> vvns_; // voronoi vertices neighboring sites
        std::vector<std::vector<int>> vvcs_; // voronoi vertices corresponding sites



    private:

        easy3d::SurfaceMesh* mesh_;

        FilletDetectorParameters parameters_;



    };
}


#endif //DEFILLET_FILET_DETECTOR_H
