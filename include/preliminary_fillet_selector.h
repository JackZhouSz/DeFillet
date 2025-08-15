//
// Created by xiaowuga on 2025/8/15.
//

#ifndef PRELIMINARY_FILLET_SELECTOR_H
#define PRELIMINARY_FILLET_SELECTOR_H


#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/point_cloud.h>

using namespace easy3d;
using namespace std;

namespace DeFillet {
    class PreliminaryFilletSelector {

    public:

        /**
         * @brief Construct a PreliminaryFilletSelector and precompute geometric properties.
         *
         * @param mesh         Pointer to the mesh on which the selector operates.
         * @param num_patches  Number of patches considered during selection.
         * @param radius_thr   Threshold for radius filtering.
         * @param angle_thr    Threshold for angle filtering.
         */
        PreliminaryFilletSelector(SurfaceMesh* mesh, int num_patches, float radius_thr, float angle_thr);

        ~PreliminaryFilletSelector();

        void apply();

        PointCloud* patch_centroid() const;

        SurfaceMesh::FaceProperty<int>  fillet_labels() const {return fillet_label_;}

    private:

        float farthest_point_sampling(int num_samples,  std::vector<SurfaceMesh::Face>& indices);

        void crop_local_patch(SurfaceMesh::Face cface, float max_gap, std::vector<SurfaceMesh::Face>& patch);

        std::vector<vec3> face_centroids(std::vector<SurfaceMesh::Face>& faces) const;

    public:

        int num_samples_, num_patches_;

        float radius_thr_, angle_thr_;

        Box3 box_;

        SurfaceMesh* mesh_;

        std::vector<SurfaceMesh::Face> patch_centroid_;

        // Face property storing the normal vector for each face in the input mesh.
        SurfaceMesh::FaceProperty<vec3> face_normals_;

        // Edge property storing the dihedral angle (in degrees) between the two faces in the input mesh.
        SurfaceMesh::EdgeProperty<float> dihedral_angle_;


        SurfaceMesh::EdgeProperty<float> centroid_distance_;

        SurfaceMesh::FaceProperty<int> fillet_label_;
    };
}


#endif //PRELIMINARY_FILLET_SELECTOR_H
