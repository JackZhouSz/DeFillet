//
// Created by 13900K on 2024/3/12.
//



#include "visualization.h"
#include "src/fillet_seg.h"



int main() {
    std::string file_path = "../data/21681_cbbaafd2_3.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    FilletSeg filletSeg(mesh);
    filletSeg.face_scoring();

    easy3d::Viewer viewer("ASD");
}