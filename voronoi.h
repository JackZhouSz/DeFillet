//
// Created by xiaowuga on 12/9/2023.
//

#ifndef DEFILLET_VORONOI_H
#define DEFILLET_VORONOI_H

#include "common.h"



class Voronoi {
public:
    Voronoi(std::vector<Point>& sites) : sites_(sites) {}
    void cal(std::vector<Point>& vertices, std::vector<std::vector<int>>& sites_regions);
private:
    std::vector<Point> sites_;
    Triangulation T_;
};


#endif //DEFILLET_VORONOI_H
