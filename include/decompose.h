//
// Created by xiaowuga on 2025/8/14.
//

#ifndef DECOMPOSE_H
#define DECOMPOSE_H

#include <easy3d/core/surface_mesh.h>

using namespace easy3d;
using namespace std;
namespace DeFillet {

/**
 * @brief Dihedral-angle based mesh decomposer: splits a surface mesh into
 *        face-connected components by cutting along "sharp" edges.
 */
class Decompose {
public:
    /**
     * @brief Constructor.
     * @param mesh      Input mesh (non-owning pointer). It must remain valid and
     *                  read-only for the lifetime of this object.
     *                  A near-manifold triangle mesh is recommended. For polygon meshes,
     *                  dihedral angles are computed from adjacent face normals.
     *                  Degenerate/zero-area faces may yield unstable normals and end up
     *                  as tiny singleton components.
     * @param angle_thr Dihedral angle threshold in **degrees** ([0, 180]).
     */
    Decompose(SurfaceMesh* mesh, float angle_thr);

   /**
    * @brief Destructor.
    *
    * Removes the auxiliary properties created during the mesh decomposition
    * process (e.g., per-edge dihedral angles) to free memory and prevent
    * leaving stale data in mesh_.
    */
    ~Decompose();


    /**
     * @brief Perform the decomposition.
     *
     * Behavior:
     * - Uses current angle_thr_ and mesh_ to compute face-connected components.
     * - Boundary and non-manifold edges are treated as cuts.
     * - Results are cached for export via all_components().
     *
     * Note: call apply() again after set_angle_thr() to recompute with the new threshold.
     */
     void apply();

    /**
     * @brief Export all components as independent meshes.
     * @return A vector of heap-allocated SurfaceMesh* (one per component).
     *         **Ownership** transfers to the caller (who must `delete` them).
     *         - Vertices are reindexed per component.
     *         - Common attributes (positions/normals/colors/UVs) may be copied/interpolated
     *           where applicable (implementation-dependent).
     *         - Custom attributes that cannot be safely transferred may be omitted unless
     *           explicitly handled in the implementation.
     *
     */
     std::vector<SurfaceMesh*> all_components() {return components_;}



    /**
     * @brief Set a new dihedral angle threshold in degrees.
     * @param angle_thr New threshold in [0, 180]. This does not recompute immediately.
     *                  Call apply() after changing the threshold to refresh the result.
     */
    void set_angle_thr(float angle_thr) {angle_thr_ = angle_thr;}

private:

    SurfaceMesh* split_component(int label);

    std::vector<SurfaceMesh*> split_alone_egde(SurfaceMesh* component,
                          SurfaceMesh::Edge e, int component_id,
                          int new_id0, int new_id1);

    void clear_all_components();


     /**
     * @brief Compute the (unsigned) angle between two 3D vectors in degrees.
     *
     * Uses atan2(||n1 × n2||, n1 · n2), which is generally more numerically stable
     * than acos of the normalized dot product.
     *
     * @param n1 First vector (e.g., a face normal).
     * @param n2 Second vector (e.g., a face normal).
     * @return Angle in degrees, in the range [0, 180].
     */
     double angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2);

private:
    SurfaceMesh* mesh_; //< Input mesh (non-owning; lifetime managed by the caller).
    float angle_thr_;   //< Dihedral angle threshold in degrees; controls where to cut.

    std::vector<SurfaceMesh*> components_;

    SurfaceMesh::FaceProperty<int> component_labels_;
    SurfaceMesh::FaceProperty<vec3> face_normals_;

    SurfaceMesh::EdgeProperty<float> dihedral_angle_;

};
}



#endif //DECOMPOSE_H
