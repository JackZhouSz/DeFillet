//
// Created by xiaowuga on 2025/8/15.
//

#include <preliminary_fillet_selector.h>
#include <utils.h>

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>

int main() {

    std::string case_path = "../data/22162_439af02a_0.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);

    DeFillet::PreliminaryFilletSelector pfs(mesh, 500, 0.05, 30);

    pfs.apply();

    easy3d::PointCloud* samples = pfs.patch_centroid();

    easy3d::PointCloudIO::save("../out3/samples.ply" , samples);

    auto labels = pfs.fillet_labels();

    DeFillet::save_fillet_regions(mesh, labels, "../out3/fillet_seg.ply");

    return 0;
}