//
// Created by xiaowuga on 2024/11/14.
//
#include <fillet_seg_v5.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/point_cloud_io.h>

#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>

#include <mean_shift.h>
#include <voro_viewer.h>

int main() {
    easy3d::logging::initialize(false, true,
                 true,
                true);
    // int num_threads = omp_get_num_threads();
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    printf("Number of threads: %d\n", num_threads);
    std::string path = "D:\\code\\sig25_defillet\\22162_439af02a_0.ply";
    // std::string cloud_path = "D:\\code\\sig25_defillet\\res3.ply";
    // easy3d::PointCloud* cloud = easy3d::PointCloudIO::load(cloud_path);
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(path);
    FilletSegV5 fillet_seg_v5(mesh);

    fillet_seg_v5.run();
    auto vv = fillet_seg_v5.vv_;
    auto vvns = fillet_seg_v5.vvns_;
    auto vvcs = fillet_seg_v5.vvcs_;
    auto scvv = fillet_seg_v5.scvv_;
    easy3d::VoroViewer viewer("VoroViewer");
    viewer.init(vv, mesh, vvns, vvcs, scvv);
    viewer.run();
    // MeanShift ms(cloud->points());
    // std::vector<std::pair<easy3d::vec3,int>> res= ms.run(cloud->points(), 0.05);
    // easy3d::PointCloud* out = new easy3d::PointCloud;
    // for(auto i : res) {
    //     out->add_vertex(i.first);
    // }
    // easy3d::PointCloudIO::save("../as1d.ply", out);
    // std::cout << "done." <<std::endl;
    return 0;
}
