//
// Created by xiaowuga on 2024/11/22.
//
#include <iostream>
#include <mean_shift.h>
#include <easy3d/fileio/point_cloud_io.h>
int main() {
    easy3d::PointCloud* cloud = easy3d::PointCloudIO::load("../scvv.ply");
    MeanShift ms(cloud->points());
    std::vector<std::pair<easy3d::vec3,int>> res = ms.run(cloud->points(), 0.1);
    easy3d::PointCloud* cloud1 = new easy3d::PointCloud;
    for(auto it : res) {
        cloud1->add_vertex(it.first);
        std::cout << it.second << std::endl;
    }
    easy3d::PointCloudIO::save("../asd.ply", cloud1);
    return 0;
}