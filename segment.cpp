//
// Created by xiaowuga on 3/7/2024.
//

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <GCoptimization.h>
#include <map>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/viewer/viewer.h>
class DataCost : public GCoptimization::DataCostFunctor {
public:
    std::vector<std::vector<double>> data_;
    int nb_nodes_, label_num_;
    DataCost(std::vector<std::vector<double>>& data, int nb_node, int label_num) : data_(data),
                                                                      nb_nodes_(nb_node), label_num_(label_num) {}

    virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s,
                                                   GCoptimization::LabelID l) override {
        return data_[s][l];
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


int main() {
    std::ifstream file("C:\\Users\\xiaowuga\\Documents\\WeChat Files\\wxid_s2unj7vcku6p21\\FileStorage\\File\\2024-03\\soft_label\\01J07922_upper.txt"); // 打开名为data.txt的文件
    int num_class = 17;
    std::vector<std::vector<double>> data;
    std::vector<double> dd;
    if (file.is_open()) {
        double num;
        while (file >> num) { // 逐个数字读取文件中的数据
            dd.emplace_back(num);
        }
        int nb_points = dd.size() / num_class;
        data.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            for(int j = 0; j < num_class; j++) {
                data[i].emplace_back(dd[i * num_class + j]);
            }
        }
        file.close(); // 关闭文件
    } else {
        std::cout << "Unable to open file" << std::endl;
    }

    int nb_points = data.size();
    double a = 0;
    for(int i = 0; i < nb_points; i++) {
        double sum = 0;
        for(int j = 0; j < num_class; j++) {
            sum += data[i][j];
        }
        double b = 0.0;
        for(int j = 0; j < num_class; j++) {
            data[i][j] =  1.0 - data[i][j] / sum;
            if(data[i][j] > b) {
                b = data[i][j];
            }
        }
        a = std::min(b, a);
    }
    std::cout << a <<std::endl;

    std::string file_path = "D:\\code\\GraphCutExample\\test_for_GC\\01J07922_upper.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    std::map<int, std::map<int, double>> mp;
    double maxx = 0;
    double minn = std::numeric_limits<double>::max();
    int nb_faces = mesh->n_faces();
    std::vector<std::vector<double>>data_cost(nb_faces);

    for(auto f : mesh->faces()) {
        std::vector<double> tmp(num_class, 0);
        int ct = 0;
        for(auto v : mesh->vertices(f)) {
            for(int i = 0; i < num_class; i++) {
                tmp[i] += data[v.idx()][i];
            }
            ct++;
        }
        double sum = 0;
        for(int i = 0; i < num_class; i++) {
            tmp[i] /= ct;
            sum += tmp[i];
        }
        for(int i = 0; i < num_class; i++) {
            tmp[i] = 1.0 - tmp[i] / sum;
        }
        data_cost[f.idx()] = tmp;
    }

    for(auto e : mesh->edges()) {
        double len = mesh->edge_length(e);
        maxx = std::max(len, maxx);
        minn = std::min(len, minn);
    }
    DataCost data_item(data_cost, nb_faces, num_class);




    double alpha = 0.5;

    GCoptimizationGeneralGraph gcp(nb_faces, num_class);
    std::vector<std::pair<int,int>> edges;
    std::vector<double> edge_weights;
    for(auto e : mesh->edges()) {
        auto f1 = mesh->face(e, 0);
        auto f2 = mesh->face(e, 1);
        if(f1.is_valid() && f2.is_valid()) {
            int id1 = f1.idx();
            int id2 = f2.idx();
            if (id1 > id2) {
                std::swap(id1, id2);
            }
            gcp.setNeighbors(id1, id2);
            double len = mesh->edge_length(e);
            edges.emplace_back(std::make_pair(id1, id2));
            double val = alpha * (len - minn) / (maxx - minn);
            edge_weights.emplace_back(val);
        }
    }
//    std::cout << "ASD" <<std::endl;*/
    SmoothCost smooth_item(edges, edge_weights);

    gcp.setDataCostFunctor(&data_item);
    gcp.setSmoothCostFunctor(&smooth_item);

    std::cout << "Before optimization energy is " << gcp.compute_energy() << std::endl;
    gcp.expansion(10);
    std::cout << "After optimization energy is " << gcp.compute_energy() << std::endl;
    std::vector<int> labels;
    labels.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        labels[i] = gcp.whatLabel(i);
        std::cout << labels[i] <<std::endl;
    }

    std::vector<easy3d::vec3> color_tables(num_class);
    for(int i = 0; i < num_class; i++) {
        color_tables[i] = easy3d::random_color();
    }
    easy3d::Viewer viewer("vector_field");
    viewer.add_model(mesh);
    auto drawable = mesh->renderer()->get_triangles_drawable("faces");
    std::string color_name = "f:color";
    auto coloring = mesh->face_property<easy3d::vec3>(color_name, easy3d::vec3(0, 0, 0));
    for(auto f : mesh->faces()) {
        std::vector<int> hist(num_class, 0);
        for(auto v : mesh->vertices(f)) {
            hist[labels[v.idx()]]++;
        }
        int maxx = 0;
        int id = 0;
        for(int i = 0; i < num_class; i++) {
            if(hist[i] > maxx) {
                maxx = hist[i];
                id = i;
            }
        }
//        std::cout << id << std::endl;
        coloring[f] = color_tables[id];
    }
    drawable->set_property_coloring(easy3d::State::FACE, color_name);
    drawable->update();
    viewer.update();
    viewer.run();
    return 0;
}