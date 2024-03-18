//
// Created by 13900K on 2024/3/18.
//

#include "surafce_mesh_segmenter.h"

namespace easy3d {
    SurafceMeshSegmenter::SurafceMeshSegmenter(SurfaceMesh *mesh) : mesh_(mesh) {
    }
    SurafceMeshSegmenter::~SurafceMeshSegmenter() {

    }

    template<typename FT>
    SurfaceMesh* SurafceMeshSegmenter::segment(const SurfaceMesh::FaceProperty<FT>& segments, FT label) {
        for(auto f : mesh_-)
    }
}