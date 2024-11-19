//
// Created by xiaowuga on 2024/11/7.
//

#ifndef KNN4D_H
#define KNN4D_H

#include "nanoflann.hpp"
#include <vector>

namespace KNN {
    struct Point4D {
        double p[4];
        Point4D(double x = 0.0, double y = 0.0, double z = 0.0, double w = 0.0) {
            p[0] = x;p[1] = y;p[2] = z; p[3] = w;
        }
    };

    class KdSearch4D {
    public:
        KdSearch4D(std::vector<Point4D>& points);
        ~KdSearch4D();
        // int radius_search(const Point4D& p, double squared_radius, std::vector<size_t> &neighbors,
        //         std::vector<double> &squared_distances) const;
        void kth_search(const Point4D& p, int k, std::vector<size_t> &neighbors,
                        std::vector<double> &squared_distances) const;
    protected:
        std::vector<Point4D>* points_;
        void* tree_;
    };
}



#endif //KNN4D_H
