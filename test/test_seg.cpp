//
// Created by xiaowuga on 2025/8/15.
//

#include <fillet_detector.h>
#include <utils.h>

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>

#include <omp.h>

int main() {

    std::string case_path = "../data/22181_526fb582_0_remesh.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(case_path);

    easy3d::vec3 mean(0,0,0);
    for(auto vertex : mesh->vertices()) {
        mean += mesh->position(vertex);
    }

    mean /= mesh->n_vertices();

    for(auto vertex : mesh->vertices()) {
        mesh->position(vertex) = mesh->position(vertex) - mean;
    }
    float maxx = 0;
    for(auto vertex : mesh->vertices()) {
        float len = (mesh->position(vertex)).norm();
        maxx = max(len, maxx);
    }

    for(auto vertex : mesh->vertices()) {
        mesh->position(vertex) = mesh->position(vertex) / maxx;
    }

    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);

    DeFillet::FilletDetectorParameters parameters;
    parameters.epsilon = 0.02;
    parameters.radius_thr = 0.05;
    parameters.angle_thr = 30;
    parameters.sigma = 1.0;
    parameters.lamdba = 10000;
    parameters.num_patches = 2000;
    parameters.num_neighbors = 300;
    parameters.num_smooth_iter = 10;

    parameters.num_sor_iter = 3;
    parameters.num_sor_neighbors = 10;
    parameters.num_sor_std_ratio = 1.0;

    DeFillet::FilletDetector detector(mesh, parameters);


    detector.generate_voronoi_vertices();
    detector.filter_voronoi_vertices();


    detector.compute_voronoi_vertices_density_field();

    detector.rolling_ball_trajectory_transform();

    detector.compute_fillet_radius_field();

    detector.compute_fillet_radius_rate_field();

    detector.rate_field_smoothing();

    detector.graph_cut();


    easy3d::PointCloud* vv = detector.voronoi_vertices();
    easy3d::PointCloudIO::save("../out3/voronoi_vertices.ply" , vv);

    easy3d::PointCloud* center = detector.rolling_ball_centers();
    easy3d::PointCloudIO::save("../out3/rolling_ball_centers.ply" , center);

    std::vector<float> radius_field = detector.radius_field();
    DeFillet::save_field(mesh, radius_field, "../out3/radius_field.ply");

     std::vector<float> rate_field = detector.radius_rate_field();

     DeFillet::save_field(mesh, rate_field, "../out3/rate_field.ply");

    std::vector<int> fillet_labels = detector.fillet_labels();

    DeFillet::save_fillet_segmentation(mesh, fillet_labels, "../out3/segmentation.ply");

    return 0;
}