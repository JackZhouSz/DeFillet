//
// Created by xiaowuga on 3/11/2024.
//

#ifndef DEFILLET_FILLET_SEG_H
#define DEFILLET_FILLET_SEG_H

#include <easy3d/core/surface_mesh.h>

class FilletSeg {
public:
    FilletSeg(easy3d::SurfaceMesh* mesh) : mesh_(mesh),eps_(0.03), s_(10.0), runtime_(0.0){}
    void seg();
    void site_scoring();
    void run_gc();
    void boardcast_vertex_domain(int vid, double  score);
    double cal_vertex_score(int vid);
    void set_eps(double eps) {eps_ = eps;}
    void set_s(double s) {s_ = s;}
    double get_runtime() {}
private:
    easy3d::SurfaceMesh* mesh_;
    std::vector<easy3d::vec3> sites_;
    std::vector<easy3d::vec3> vertices_;
    std::vector<std::vector<int>> site_of_vertices_;
    double eps_;
    double s_;
    double runtime_;
};


#endif //DEFILLET_FILLET_SEG_H
