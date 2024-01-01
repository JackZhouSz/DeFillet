//
// Created by xiaowuga on 12/16/2023.
//

#ifndef DEFILLET_KNN_H
#define DEFILLET_KNN_H

#include "nanoflann.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

namespace KNN {
    struct Point {
        double p[3];
        Point(double x = 0.0, double y = 0.0, double z = 0.0) {
            p[0] = x;p[1] = y;p[2] = z;
        }
    };

    class KdSearch {
    public:
        KdSearch(std::vector<Point>& points);
        ~KdSearch();
        int radius_search(const Point& p, double squared_radius, std::vector<int> &neighbors,
                std::vector<double> &squared_distances) const;
    protected:
        std::vector<Point>* points_;
        void* tree_;
    };
}


#endif //DEFILLET_KNN_H
