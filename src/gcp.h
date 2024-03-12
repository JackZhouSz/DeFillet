//
// Created by 小乌嘎 on 2024/2/14.
//

#ifndef DEFILLET_GCP_H
#define DEFILLET_GCP_H

#include <iostream>
#include <vector>
#include <map>
#include "GCoptimization.h"

namespace GCP {
    class DataCost : public GCoptimization::DataCostFunctor {
    public:
        std::vector<double> data_;
        int nb_nodes_, label_num_;
        DataCost(std::vector<double>& data, int nb_node, int label_num) : data_(data),
                                                                          nb_nodes_(nb_node), label_num_(label_num) {}

        virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s,
                                                       GCoptimization::LabelID l) override {
            return data_[l * nb_nodes_ + s];
        }
    };

    class SmoothCost : public GCoptimization::SmoothCostFunctor{
    public:
        std::map<int, std::map<int, double>> data_;

        SmoothCost(const std::vector<std::pair<int,int>> & edges,
                   const std::vector<double>& edges_weight) {
            int nb_data = edges.size();
            for(int i = 0; i < nb_data; i++) {
                int x = edges[i].first, y = edges[i].second;
                if(x > y) {
                    std::swap(x, y);
                }
                data_[x][y] = edges_weight[i];
            }
        }

        virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s1, GCoptimization::SiteID s2,
                                                       GCoptimization::LabelID l1, GCoptimization::LabelID l2) override {
            if(l1 == l2) {
                return 0;
            }
            else {
                int id1 = s1, id2 = s2;
                if(id1 > id2) {
                    std::swap(id1, id2);
                }
                double res = data_[id1][id2];
                return res;
            }
        }

    };
}


#endif //DEFILLET_GCP_H
