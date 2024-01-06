//
// Created by xiaowuga on 12/17/2023.
//

#ifndef DEFILLET_DATA_CONVERT_H
#define DEFILLET_DATA_CONVERT_H

#include "common.h"
#include "knn.h"

inline void my_points_convert_to_cgal_points(const std::vector<std::array<double, 3>>& my_points, std::vector<Point>& cgal_points) {
    cgal_points.clear();
    int nb_points = my_points.size();
    cgal_points.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        cgal_points[i] = Point(my_points[i][0],my_points[i][1], my_points[i][2]);
    }
}

inline void cgal_points_convert_to_my_points(const std::vector<Point>& cgal_points, std::vector<std::array<double, 3>>& my_points) {
    my_points.clear();
    int nb_points = cgal_points.size();
    my_points.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        my_points[i] = {cgal_points[i].x(), cgal_points[i].y(), cgal_points[i].z()};
    }
}

inline void cgal_points_convert_to_knn_points(const std::vector<Point>& cgal_points, std::vector<KNN::Point>& knn_points) {
    knn_points.clear();
    int nb_points = cgal_points.size();
    knn_points.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        knn_points[i] = KNN::Point(cgal_points[i].x(), cgal_points[i].y(), cgal_points[i].z());
    }
}



#endif //DEFILLET_DATA_CONVERT_H
