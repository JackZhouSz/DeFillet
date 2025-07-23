//
// Created by xiaowuga on 2025/3/17.
//

#include <iostream>
#include <vector>
#include <map>
#include <easy3d/core/surface_mesh.h>

#include "GCoptimization.h"
#include <easy3d/fileio/surface_mesh_io.h>
#include "util.h"

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


int main() {
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load("C:\\Users\\xiaowuga\\Desktop\\error_values.ply");
    for(auto fp : mesh->face_properties()) {
        std::cout << fp << std::endl;
    }
    double lamdba = 1000000;
    auto value = mesh->get_face_property<float>("f:value");
    const std::vector<float> b = value.vector();
    std::vector<float> a = b;
    float maxx = *std::max_element(a.begin(), a.end());
    float minn = *std::min_element(a.begin(), a.end());
    std::cout <<maxx << ' ' << minn <<std::endl;
    for(auto f : mesh->faces()) {
        value[f] = (value[f] - minn) / (maxx - minn);
    }
    int nb_face = mesh->n_faces();
    GCoptimizationGeneralGraph gc(nb_face, 2);

    std::vector<double> data_cost(2 * nb_face);
    for(int i = 0; i < nb_face; i++) {
        easy3d::SurfaceMesh::Face f(i);
        data_cost[i] = value[f] / nb_face / 2;
        data_cost[i + nb_face] = (1.0 - value[f]) / nb_face / 2;
    }
    std::cout << "ASD" << std::endl;
    int num_edges = mesh->n_edges();
    std::vector<std::pair<int,int>> edges;
    std::vector<double> edge_weights;
    double sum = 0.0;

    for(auto e : mesh->edges()) {
        int x = mesh->vertex(e, 0).idx();
        int y = mesh->vertex(e, 1).idx();
        double ew = lamdba * mesh->edge_length(e);
        edges.emplace_back(std::make_pair(x,y));
        edge_weights.emplace_back(ew);
        sum += ew;
    }

    for(int i = 0; i < num_edges; i++) {
        edge_weights[i] /= sum;
    }
    GCP::DataCost data_item(data_cost, nb_face, 2);
    GCP::SmoothCost smooth_item(edges, edge_weights);
    gc.setDataCostFunctor(&data_item);
    gc.setSmoothCostFunctor(&smooth_item);
    std::cout << "Before optimization energy is " << gc.compute_energy() << std::endl;
    gc.expansion(10);
    std::cout << "After optimization energy is " << gc.compute_energy() << std::endl;

    auto gcp = mesh->face_property<int>("f:fillet_labels");

    for(int i = 0; i < nb_face; i++) {
        easy3d::SurfaceMesh::Face f(i);
        if(gc.whatLabel(i) == 0) {
            gcp[f] = 0;
        } else {
            gcp[f] = 1;
        }
    }
    auto color = mesh->face_property<easy3d::vec3>("f:color");
    for(auto f : mesh->faces()) {
        if(gcp[f] == 0) {
            color[f] = easy3d::vec3(0,0,1.0);
        }
        else {
            color[f] = easy3d::vec3(1.0,0,0);
        }
    }
    easy3d::SurfaceMeshIO::save("../error_values.ply", mesh);
    save_mesh_field(mesh, value.vector(),"../aa.ply");
    return 0;
}
