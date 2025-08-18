//
// Created by xiaowuga on 2024/11/19.
//

#ifndef MEAN_SHIFT_H
#define MEAN_SHIFT_H


#include<easy3d/core/point_cloud.h>

namespace DeFillet {
    class MeanShift {
    public:
        MeanShift() {
            set_kernel(NULL);
        }

        void shift_point(easy3d::vec3 &point, const std::vector<easy3d::vec3> &data, double kernel_bandwidth,
                         easy3d::vec3 &shift_point);

        std::pair<easy3d::vec3, int>
        run(const easy3d::vec3 &points, const std::vector<easy3d::vec3> &data, double kernel_bandwidth,
            double eps = 1e-4);

        std::vector<std::pair<easy3d::vec3, int>>
        run(const std::vector<easy3d::vec3> &points, const std::vector<easy3d::vec3> &data, double kernel_bandwidth,
            double eps = 1e-4);

    public:
        double (*kernel_func)(double, double);

        void set_kernel(double (*_kernel_func)(double, double));

    private:
        // std::vector<easy3d::vec3> data_;
    };

}


#endif //MEAN_SHIFT_H
