//
// Created by xiaowuga on 2024/11/19.
//

#include "mean_shift.h"


double gaussian_kernel(double distance, double kernel_bandwidth){
    double temp =  exp(-1.0/2.0 * (distance*distance) / (kernel_bandwidth*kernel_bandwidth));
    return temp;
}

// void MeanShift::set_kernel( double (*_kernel_func)(double,double) ) {
//     if(!_kernel_func){
//         kernel_func = gaussian_kernel;
//     } else {
//         kernel_func = _kernel_func;
//     }
// }

void MeanShift::shift_point(easy3d::vec3& point
                            ,double kernel_bandwidth
                            ,easy3d::vec3& shift_point) {
    shift_point = easy3d::vec3(0,0,0);
    double tw = 0;  // total weight
    int num = data_.size();
#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        easy3d::vec3 tmp = data_[i];
        double d = (point - tmp).norm(); //distance
        double w = gaussian_kernel(d, kernel_bandwidth); // weight
#pragma omp critical
        {
            shift_point += w * tmp;
            tw += w;
        }
    }
    double tw_inv = 1.0 / tw;
    shift_point *= tw_inv;
}


std::vector<std::pair<easy3d::vec3,int>> MeanShift::run(const std::vector<easy3d::vec3>& points
                                                       ,double kernel_bandwidth
                                                       ,double eps) {
    const double eps_sqr = eps * eps;
    int num = points.size();
    // std::cout << num <<std::endl;
    std::vector<bool> stop_moving(num, false);
    std::vector<easy3d::vec3> shift_points = points;
    double max_shift_distance;
    int max_iter_num = 100, iter = 0;
    do {
        max_shift_distance = 0;
        for(int i = 0; i < num; i++) {
            if(!stop_moving[i]) {
                easy3d::vec3 tmp;
                shift_point(shift_points[i], kernel_bandwidth,  tmp);
                double shift_distance_sqr = (tmp - shift_points[i]).length2();
                max_shift_distance = std::max(max_shift_distance, shift_distance_sqr);
                if(shift_distance_sqr < eps_sqr) {
                    stop_moving[i] = true;
                }
                shift_points[i] = tmp;
            }
        }
        // printf("max_shift_distance: %f\n", sqrt(max_shift_distance));
        iter++;
    } while(max_shift_distance > eps_sqr && iter < max_iter_num);

    std::vector<std::pair<easy3d::vec3,int>> res;
    for(int i = 0; i < num; i++) {
        double minn = std::numeric_limits<double>::max();
        int idx = 0;
        for(size_t j = 0; j < data_.size(); j++) {
            double dis = (data_[j] - shift_points[i]).norm();
            if(minn > dis) {
                minn = dis; idx = j;
            }
        }
        res.emplace_back(std::make_pair(shift_points[i], idx));
    }
    return res;
}

std::pair<easy3d::vec3,int> MeanShift::run(const easy3d::vec3& point
                                          ,double kernel_bandwidth
                                          ,double eps) {
    std::vector<easy3d::vec3> points = {point};
    std::vector<std::pair<easy3d::vec3, int>> res = run(points, kernel_bandwidth, eps);
    return res[0];
}