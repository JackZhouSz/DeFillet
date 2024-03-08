//
// Created by xiaowuga on 12/9/2023.
//

#ifndef DEFILLET_VORONOI_H
#define DEFILLET_VORONOI_H

#include "common.h"

class Voronoi {
public:
    Voronoi(std::vector<CGAL_Point>& sites, double scale = 0.0001) : sites_(sites), scale_(scale) {}
    void cal();
    void cal_v1();
    void set_scale(double scale) {scale_ = scale;};
    const std::vector<CGAL_Point>& get_sites() const { return sites_;}

    const std::vector<CGAL_Point>& get_vertices() const { return vertices_;}

    const std::vector<std::vector<int>>& get_regions() const {return regions_;}

    const std::vector<CGAL_Point>& get_finite_vertices() const { return finite_vertices_;}

    const std::vector<std::vector<int>>& get_finite_regions() const { return finite_regions_;}

    const std::vector<std::vector<int>>& get_finite_cell_pole() const {return finite_cell_pole_;}
    const std::vector<std::vector<int>>& get_cell_pole() const {return cell_pole_;}
    std::vector<double> cal_radius();
private:
    void cal_finite_vertices_and_regions();
    void cal_finite_cell_pole();
    void cal_cell_pole();
private:
    std::vector<CGAL_Point> sites_;
    std::vector<CGAL_Point> vertices_;
    std::vector<std::vector<int>> regions_;
    std::vector<std::vector<int>> vertices_graph_;

    std::vector<CGAL_Point> finite_vertices_;
    std::vector<std::vector<int>> finite_regions_;
    std::vector<std::vector<int>> finite_cell_pole_;
    std::vector<std::vector<int>> finite_vertices_graph_;

    std::vector<std::vector<int>> cell_pole_;

    Triangulation T_;

    double scale_;
};


#endif //DEFILLET_VORONOI_H
