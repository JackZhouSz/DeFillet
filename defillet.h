//
// Created by 小乌嘎 on 2024/2/10.
//

#ifndef DEFILLET_DEFILLET_H
#define DEFILLET_DEFILLET_H
#include <iostream>
#include "utils.h"

namespace DEFILLET {
    void compute_fillet_field(const std::vector<Eigen::Vector3d>& points,
                              const std::vector<std::vector<size_t>>& faces,
                              const std::vector<Eigen::Vector3d>& vor_vertices,
                              const std::vector<std::vector<int>>& vor_site_poles,
                              std::vector<double>& density_field,
                              std::vector<double>& fillet_field,
                              std::vector<size_t>& fillet_poles,
                              double eps = 1e-2,
                              double radius = -1.0);
    void compute_fillet_field_v1(const std::vector<Eigen::Vector3d>& points,
                              const std::vector<std::vector<size_t>>& faces,
                              const std::vector<Eigen::Vector3d>& vor_vertices,
                              const std::vector<std::vector<int>>& vor_site_poles,
                              std::vector<double>& density_field,
                              std::vector<double>& fillet_field,
                              std::vector<size_t>& fillet_poles,
                              double eps = 1e-2);

    void run_graph_cut(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       const std::vector<Eigen::Vector3d>& vor_vertices,
                       const std::vector<size_t>& fillet_poles,
                       std::vector<double>& fillet_field,
                       std::vector<double>& labels,
                       double alpha = 0.5,
                       double thr = -1.0);

    void run_graph_cut_v1(std::vector<Eigen::Vector3d>& points,
                       std::vector<std::vector<size_t>>& faces,
                       std::vector<double>& fillet_field,
                       std::vector<double>& labels,
                       double alpha_1 = 0.5,
                       double alpha_2 = 0.5);

    void extract_interest_region(const std::vector<Eigen::Vector3d>& points,
                                 const std::vector<std::vector<size_t>>& faces,
                                 const std::vector<double>& label,
                                 std::vector<Eigen::Vector3d>& fillet_points,
                                 std::vector<std::vector<size_t>>& fillet_faces,
                                 std::vector<int>& point_map,
                                 std::vector<size_t>& fillet_bounding,
                                 std::vector<Eigen::Vector3d>& fillet_bounding_normals);



    void run_geodesic(const std::vector<Eigen::Vector3d>& points,
                      const std::vector<std::vector<size_t>>& faces,
                      const std::vector<size_t>& sources,
                      const std::vector<Eigen::Vector3d>& sources_normals,
                      std::vector<size_t>& ancestor,
                      std::vector<double>& distance,
                      std::vector<Eigen::Vector3d>& point_tar_normals,
                      std::vector<Eigen::Vector3d>& face_tar_normals);

    bool optimize_qem(const std::vector<Eigen::Vector3d>& points,
                  const std::vector<std::pair<int,int>>& edges,
                  std::vector<Eigen::Vector3d>& normals,
                  std::vector<Eigen::Vector3d>& cur_normals,
                  std::vector<Eigen::Vector3d>& new_points,
                  std::vector<size_t>& fixed_points,
                  double beta);

    bool optimize_centroid(const std::vector<Eigen::Vector3d>& points,
                  const std::vector<std::pair<int,int>>& edges,
                  std::vector<Eigen::Vector3d>& normals,
                  std::vector<Eigen::Vector3d>& cur_normals,
                  std::vector<Eigen::Vector3d>& new_points,
                  std::vector<size_t>& fixed_points,
                  double beta);

    bool iterative_optimize(const std::vector<Eigen::Vector3d>& points,
                            const std::vector<std::vector<size_t>>& faces,
                            std::vector<Eigen::Vector3d>& normals,
                            std::vector<Eigen::Vector3d>& new_points,
                            std::vector<size_t>& fixed_points,
                            double beta,
                            int num_iterations = 10);
}


#endif //DEFILLET_DEFILLET_H
