//
// Created by xiaowuga on 12/9/2023.
//

#include "voronoi.h"
#include <algorithm>


void Voronoi::cal(std::vector<Point> &vertices, std::vector<std::vector<int>> &sites_regions) {
    int nb_sites = sites_.size();
    std::map<Vertex_handle, int> mp;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T_.insert(sites_[i]);
        mp[v] = i;
    }

    vertices.clear(); sites_regions.clear();
    sites_regions.resize(nb_sites);

    for(auto cell = T_.finite_cells_begin(); cell != T_.finite_cells_end(); cell++) {
        auto v = T_.dual(cell);
        int vid = vertices.size();
        vertices.emplace_back(v);
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp[vh];
            sites_regions[id].emplace_back(vid);
        }
    }

}

