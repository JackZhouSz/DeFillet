//
// Created by xiaowuga on 12/17/2023.
//

#ifndef DEFILLET_DEFILLET_H
#define DEFILLET_DEFILLET_H

#include "common.h"
#include <GCoptimization.h>


namespace DEFILLET {

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

    void computes_vertices_density_field(const std::vector<Point>& vertices, double radius,
                                         std::vector<double>& density_field);

    void compute_sites_pole_radius_field(const std::vector<Point>& sites,
                                         const std::vector<Point>& vertices,
                                         const std::vector<std::vector<int>>& site_regions,
                                         const std::vector<double>& density_field,
                                         std::vector<double>& site_density_field,
                                         std::vector<double>& site_pole_radius_field);


    void run_graph_cut(const std::vector<double>& node_weight,
                       const std::vector<std::pair<int,int>>& edges,
                       const std::vector<double>& edge_weight,
                       std::vector<int>& labels,
                       double alpha);

    void find_island(const std::vector<std::vector<std::pair<int,double>>>& mesh_graph,
                          std::vector<double>& scalar_field,
                          std::vector<int>& island,
                          double scalar_thr,
                          double island_thr);

    void corrosion_island(const std::vector<Point>& points,
                          const std::vector<int> island,
                          double radius,
                          double scalar_thr,
                          std::vector<double>& scalar_field);

    void classify_defillet_triangle_via_vertices_voting(const std::vector<std::vector<size_t>>& f_ind,
                                                        const std::vector<double>& vertices_field,
                                                        double thr,
                                                        std::vector<bool>& is_fillet);

    void extract_interest_region(const std::vector<Point>& points,
                               const std::vector<std::vector<size_t>>& faces,
                               const std::vector<bool>& is_fillet,
                               std::vector<Point>& fillet_points,
                               std::vector<std::vector<size_t>>& fillet_faces,
                               std::vector<int>& point_map_table,
                               std::vector<int>& faces_map_table,
                               std::vector<bool>& face_labels);

    void run_geodesic(const std::vector<Point>& points,
                      const std::vector<std::vector<size_t>>& faces,
                      const std::vector<int>& centroid,
                      const std::vector<bool>& centroid_labels,
                      std::vector<int>& ancestor,
                      std::vector<double>& distance);

    void geodesic_post_processing(const std::vector<Point>& points,
                                  const std::vector<std::vector<size_t>>& faces,
                                  const std::vector<int>& ancestor,
                                  const std::vector<bool>& face_labels,
                                  std::vector<Vector_3>& target_normals,
                                  std::vector<int>& fixed_points,
                                  std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector,
                                  double angle_thr);


    bool optimize_with_fixed_solve(const std::vector<Point>& points,
                  const std::vector<std::vector<size_t>>& faces,
                  std::vector<Vector_3>& normals,
                  std::vector<Point>& new_points,
                  std::vector<int>& fixed_points,
                  std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector);

    bool optimize_sparseLU(const std::vector<Point>& points,
                           const std::vector<std::vector<size_t>>& faces,
                           std::vector<Vector_3>& normals,
                           std::vector<Point>& new_points,
                           std::vector<int>& fixed_points,
                           std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector);

    bool optimize_spareQR(const std::vector<Point>& points,
                          const std::vector<std::vector<size_t>>& faces,
                          std::vector<Vector_3>& normals,
                          std::vector<Point>& new_points,
                          std::vector<int>& fixed_points,
                          double beta);

    bool iterative_optimize(const std::vector<Point>& points,
                            const std::vector<std::vector<size_t>>& faces,
                            std::vector<Vector_3>& normals,
                            std::vector<Point>& new_points,
                            std::vector<int>& fixed_points,
                            double beta,
                            int nb_optimize);
}


#endif //DEFILLET_DEFILLET_H
