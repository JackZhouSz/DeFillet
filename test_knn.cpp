


#include "knn.h"



int main()
{
    // 创建数据集
//    PointCloud cloud;
//    cloud.points.emplace_back(1.0, 2.0, 3.0);
//    cloud.points.emplace_back(4.0, 5.0, 6.0);
//    cloud.points.emplace_back(7.0, 8.0, 9.0);
    std::vector<KNN::Point> points;
    points.emplace_back(KNN::Point(1.0, 2.0,3.0));
    points.emplace_back(KNN::Point(4.0, 5.0,6.0));
    points.emplace_back(KNN::Point(7.0, 8.0,9.0));

    KNN::KdSearch kds(points);
    KNN::Point query(2.0, 4.0, 6.0);
    std::vector<int> neighbors;
    std::vector<double> squared_distances;
    kds.radius_search(query, 10 * 10, neighbors, squared_distances);
    for(int i = 0; i < neighbors.size(); i++) {
        std::cout <<neighbors[i] << ' ';
    }
    std::cout << std::endl;
    for(int i = 0; i < squared_distances.size(); i++) {
        std::cout << squared_distances[i] << ' ';
    }
    std::cout << std::endl;
    return 0;
}
