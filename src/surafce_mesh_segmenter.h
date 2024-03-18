//
// Created by 13900K on 2024/3/18.
//

#ifndef DEFILLET_SURAFCE_MESH_SEGMENTER_H
#define DEFILLET_SURAFCE_MESH_SEGMENTER_H

#include <easy3d/core/surface_mesh.h>

namespace easy3d {
    class SurfaceMeshSegmenter {
    public:
        SurfaceMeshSegmenter(SurfaceMesh* mesh);
        ~SurfaceMeshSegmenter();

        template<typename FT>
        SurfaceMesh* segment(const SurfaceMesh::FaceProperty<FT>& segments, FT label);
    private:
        easy3d::SurfaceMesh* mesh_;
    };
}


#endif //DEFILLET_SURAFCE_MESH_SEGMENTER_H
