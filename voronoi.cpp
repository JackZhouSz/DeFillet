//
// Created by xiaowuga on 12/9/2023.
//

#include "voronoi.h"
#include <algorithm>


void Voronoi::cal() {
    int nb_sites = sites_.size();
    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T_.insert(sites_[i]);
        mp_sites[v] = i;
    }

    vertices_.clear(); regions_.clear();
    regions_.resize(nb_sites);

    for(auto cell = T_.finite_cells_begin(); cell != T_.finite_cells_end(); cell++) {
        auto v = T_.dual(cell);
        int vid = vertices_.size();
        vertices_.emplace_back(v);
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            regions_[id].emplace_back(vid);
        }
    }

    cal_finite_vertices_and_regions();
    cal_finite_cell_pole();
}

void Voronoi::cal_v1() {
    int nb_sites = sites_.size();
    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T_.insert(sites_[i]);
        mp_sites[v] = i;
    }

    vertices_.clear(); regions_.clear();
    regions_.resize(nb_sites);

    for(auto cell = T_.finite_cells_begin(); cell != T_.finite_cells_end(); cell++) {
        auto v = T_.dual(cell);
        int vid = vertices_.size();
        vertices_.emplace_back(v);
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            regions_[id].emplace_back(vid);
        }
    }

    cal_cell_pole();
}

void Voronoi::cal_finite_vertices_and_regions() {
    finite_vertices_.clear();
    finite_regions_.clear();
    CGAL_Point bboxMin, bboxMax;
    auto bbox = CGAL::bounding_box(sites_.begin(), sites_.end());

    bboxMin = bbox.min(); bboxMax = bbox.max();

    double diagonal = std::sqrt(CGAL::squared_distance(bboxMin, bboxMax));
    double eps = scale_ * diagonal;


    int nb_vertices = vertices_.size();
    std::vector<int> is_inside(nb_vertices);
    for(int i = 0; i < nb_vertices; i++) {
        const CGAL_Point& p = vertices_[i];
        if(p.x() >= bboxMin.x() - eps && p.x() <= bboxMax.x() + eps
           && p.y() >= bboxMin.y() - eps && p.y() <= bboxMax.y() + eps
           && p.z() >= bboxMin.z() - eps && p.z() <= bboxMax.z() + eps) {
            is_inside[i] = finite_vertices_.size();
            finite_vertices_.emplace_back(p);
        }
        else {
            is_inside[i] = -1;
        }
    }

    int nb_sites = sites_.size();
    for(int i = 0; i < nb_sites; i++) {
        const std::vector<int>& region = regions_[i];
        std::vector<int> tmp;
        int nb_region = region.size();
        for(int j = 0; j < nb_region; j++) {
            int id = is_inside[region[j]];
            if(id != -1) {
                tmp.emplace_back(id);
            }
        }
        finite_regions_.emplace_back(tmp);
    }

//    int nb_finite_vertices = finite_vertices_.size();
//    finite_vertices_graph_.clear();
//    finite_vertices_graph_.resize(nb_finite_vertices);
//    for(int i = 0; i < nb_vertices; i++) {
//        int uid = is_inside[i];
//        if(uid != -1) {
//            int nb_nei = vertices_graph_[i].size();
//            for(int j = 0; j < nb_nei; j++) {
//                int vid = is_inside[vertices_graph_[i][j]];
//                if(vid != -1) {
//                    finite_vertices_graph_[uid].emplace_back(vid);
//                }
//            }
//        }
//    }
}

void Voronoi::cal_finite_cell_pole() {
    int nb_sites = sites_.size();
    finite_cell_pole_.clear();
    finite_cell_pole_.resize(nb_sites);
    for(int i = 0; i < nb_sites; i++) {
        const CGAL_Point& s = sites_[i];
        const std::vector<int>& region = finite_regions_[i];
        int nb_region = region.size();
        int pole_id1 = -1;
        double max_value1 = std::numeric_limits<double>::min();
        for(int j = 0; j < nb_region; j++) {
            double len = std::sqrt(CGAL::squared_distance(finite_vertices_[region[j]], s));
            if(max_value1 < len) {
                max_value1 = len;
                pole_id1 = region[j];
            }
        }

        int pole_id2 = -1;
        double max_value2 = std::numeric_limits<double>::min();
        for(int j = 0; j < nb_region; j++) {
            double len = std::sqrt(CGAL::squared_distance(finite_vertices_[region[j]], s));
            auto v1 = finite_vertices_[pole_id1] - s;
            auto v2 =  finite_vertices_[region[j]] - s;
            double dot_value = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();

            if(dot_value < 0 && max_value2 < len) {
                max_value2 = len;
                pole_id2 = region[j];
            }
        }
        if(pole_id1 != -1) {
            finite_cell_pole_[i].emplace_back(pole_id1);
        }
        if(pole_id2 != -1) {
            finite_cell_pole_[i].emplace_back(pole_id2);
        }
    }
}


void Voronoi::cal_cell_pole() {
    int nb_sites = sites_.size();
    cell_pole_.clear();
    cell_pole_.resize(nb_sites);
    for(int i = 0; i < nb_sites; i++) {
        const CGAL_Point& s = sites_[i];
        const std::vector<int>& region = regions_[i];
        int nb_region = region.size();
        int pole_id1 = -1;
        double max_value1 = std::numeric_limits<double>::min();
        for(int j = 0; j < nb_region; j++) {
            double len = std::sqrt(CGAL::squared_distance(vertices_[region[j]], s));
            if(max_value1 < len) {
                max_value1 = len;
                pole_id1 = region[j];
            }
        }

        int pole_id2 = -1;
        double max_value2 = std::numeric_limits<double>::min();
        for(int j = 0; j < nb_region; j++) {
            double len = std::sqrt(CGAL::squared_distance(vertices_[region[j]], s));
            auto v1 = vertices_[pole_id1] - s;
            auto v2 =  vertices_[region[j]] - s;
            double dot_value = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
            if(dot_value < 0 && max_value2 < len) {
                max_value2 = len;
                pole_id2 = region[j];
            }
        }
        if(pole_id1 != -1) {
            cell_pole_[i].emplace_back(pole_id1);
        }
        if(pole_id2 != -1) {
            cell_pole_[i].emplace_back(pole_id2);
        }
    }
}

std::vector<double> Voronoi::cal_radius() {
    int num = finite_vertices_graph_.size();
    std::vector<double> res(num, 0);
    double x = 100;
    for(int i = 0; i < num; i++) {
        res[i] = 0;
        double maxx = std::numeric_limits<double>::min();
        double minn = std::numeric_limits<double>::max();
        for(size_t j = 0; j < finite_vertices_graph_[i].size(); j++) {
            int id = finite_vertices_graph_[i][j];
            double len = std::sqrt(CGAL::squared_distance(finite_vertices_[i],
                                                          finite_vertices_[id]));

            res[i] += len;
            maxx = std::max(maxx, len);
            minn = std::min(minn, len);
        }
        res[i] -= (maxx + minn);
        if(finite_vertices_graph_[i].size() > 2) {
            res[i] /= (finite_vertices_graph_[i].size() - 2);
        } else {
            res[i] = minn;
        }
//        res[i] /= finite_vertices_graph_[i].size();
    }
    return res;
}

