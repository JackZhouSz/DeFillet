//
// Created by xiaowuga on 2025/8/15.
//

#include <fillet_detector.h>
#include <utils.h>

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>

int main() {

    std::string case_path = "../data/20440_27177360_1_remesh.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);

    DeFillet::FilletDetectorParameters parameters;
    parameters.epsilon = 0.01;
    parameters.radius_thr = 0.05;
    parameters.angle_thr = 40;
    parameters.sigma = 1.0;
    parameters.lamdba = 0.5;
    parameters.num_patches = 800;
    parameters.num_neighbors = 30;
    parameters.num_smooth_iter = 10;

    DeFillet::FilletDetector detector(mesh, parameters);

    detector.generate_voronoi_vertices();
    detector.filter_voronoi_vertices();

    easy3d::PointCloud* vv = detector.voronoi_vertices();
    easy3d::PointCloudIO::save("../out3/samples.ply" , vv);
//    DeFillet::PreliminaryFilletSelector pfs(mesh, 500, 0.05, 30);
//
//    pfs.apply();
//
//    easy3d::PointCloud* samples = pfs.patch_centroid();
//
//    easy3d::PointCloudIO::save("../out3/samples.ply" , samples);
//
//    auto labels = pfs.fillet_labels();
//
//    DeFillet::save_fillet_regions(mesh, labels, "../out3/fillet_seg.ply");

    return 0;
}