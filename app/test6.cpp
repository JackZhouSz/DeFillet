//
// Created by xiaowuga on 2024/11/6.
//

#include <dkm_parallel.hpp>
#include <easy3d/core/point_cloud.h>
#include <easy3d/fileio/point_cloud_io.h>


int main() {
    std::string path = "D:\\code\\defillet\\cloud_clean.ply";
    easy3d::PointCloud* cloud = easy3d::PointCloudIO::load(path);
    std::vector<std::array<float, 3>> data;
    for(auto v : cloud->vertices()) {
        easy3d::vec3 p = cloud->position(v);
        data.emplace_back(std::array<float, 3>{p.x, p.y, p.z});
    }
    auto cluster_data = dkm::kmeans_lloyd_parallel(data, 1000);
    std::vector<easy3d::vec3> means;

    for (const auto& mean : std::get<0>(cluster_data)) {
        means.emplace_back(easy3d::vec3(mean[0],mean[1],mean[2]));
    }

    easy3d::PointCloud* center = new easy3d::PointCloud;
    for(auto p : means) {
        center->add_vertex(p);
    }
    easy3d::PointCloudIO::save("../dkm_cloud_clean_kmeans.ply", center);
    return 0;
}