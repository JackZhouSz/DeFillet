//
// Created by xiaowuga on 2024/11/25.
//

#ifndef KERNEL_H
#define KERNEL_H



static double gaussian_kernel(double distance, double kernel_bandwidth){
    double temp =  exp(-0.5 * (distance*distance) / (kernel_bandwidth*kernel_bandwidth));
    return temp;
}

#endif //KERNEL_H
