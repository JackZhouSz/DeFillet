//
// Created by xiaowuga on 2025/8/15.
//

#include <fillet_detector.h>
#include <utils.h>

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>

int main() {

    std::string case_path = "../data/22162_439af02a_0.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);

    DeFillet::FilletDetectorParameters parameters;
    parameters.epsilon = 0.02;
    parameters.radius_thr = 0.04;
    parameters.angle_thr = 40;
    parameters.sigma = 1.0;
    parameters.lamdba = 0.5;
    parameters.num_patches = 800;
    parameters.num_neighbors = 30;
    parameters.num_smooth_iter = 10;

    DeFillet::FilletDetector detector(mesh, parameters);

    detector.generate_voronoi_vertices();
    detector.compute_voronoi_vertices_density_field();
    detector.filter_voronoi_vertices();
    detector.rolling_ball_trajectory_transform();

    detector.compute_fillet_radius_rate_field();

    // for(int i = 0; i < 1; i++) {
        detector.rate_field_smoothing();
    // }
    // detector.density_driven_voronoi_drift();

    easy3d::PointCloud* vv = detector.voronoi_vertices();
    easy3d::PointCloudIO::save("../out3/samples.ply" , vv);

    std::vector<float> field = detector.radius_rate_field();

    DeFillet::save_rate_field(mesh, field, "../out3/field.ply");

    return 0;
}