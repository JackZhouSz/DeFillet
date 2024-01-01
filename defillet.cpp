//
// Created by xiaowuga on 12/17/2023.
//

#include "defillet.h"
#include <limits>
#include <functional>
#include <queue>
#include <Xin_Wang.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


#include <algorithm>
namespace DEFILLET {

    void computes_vertices_density_field(const std::vector<Point>& vertices,
                                         double radius,
                                         std::vector<double>& density_field) {
        int nb_vertices = vertices.size();
        std::vector<KNN::Point> knn_vertices;
        cgal_points_convert_to_knn_points(vertices, knn_vertices);
        KNN::KdSearch kds(knn_vertices);
        density_field.resize(nb_vertices);
        for(int i = 0; i < nb_vertices; i++) {
            std::vector<int> neighbors;
            std::vector<double> squared_distances;
            int num = kds.radius_search(knn_vertices[i], radius, neighbors, squared_distances);
            density_field[i] = num;
        }

//        double_vector1D_data_normalize(density_field);
    }

    void compute_sites_pole_radius_field(const std::vector<Point>& sites,
                                         const std::vector<Point>& vertices,
                                         const std::vector<std::vector<int>>& site_regions,
                                         const std::vector<double>& vertices_density_field,
                                         std::vector<double>& site_density_field,
                                         std::vector<double>& site_pole_radius_field) {
        int nb_sites = sites.size();
        int nb_vertices = vertices.size();

        site_density_field.resize(nb_sites);
        site_pole_radius_field.resize(nb_sites);

        for(int i = 0; i < nb_sites; i++) {
            const Point& s = sites[i];
            const std::vector<int>& region = site_regions[i];
            int nb_region = region.size();
            int pole_id1 = 0;
            double max_value1 = std::numeric_limits<double>::min();
            for(int j = 0; j < nb_region; j++) {
                double len = std::sqrt(CGAL::squared_distance(vertices[region[j]], s));
                if(max_value1 < len) {
                    max_value1 = len;
                    pole_id1 = region[j];
                }
            }

            int pole_id2 = -1;
            double max_value2 = std::numeric_limits<double>::min();
            for(int j = 0; j < nb_region; j++) {
                double len = std::sqrt(CGAL::squared_distance(vertices[region[j]], s));
                auto v1 = vertices[pole_id1] - s;
                auto v2 =  vertices[region[j]] - s;
                double dot_value = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();

                if(dot_value < 0 && max_value2 < len) {
                    max_value2 = len;
                    pole_id2 = region[j];
                }
            }
            double max_value = 0;
            int id = pole_id1;
            if(pole_id2 == -1) {
                id = pole_id1;
                max_value = max_value1;
            }
            else {
                if(vertices_density_field[pole_id1] > vertices_density_field[pole_id2]) {
                    id = pole_id1;
                    max_value = max_value1;
                }
                else {
                    id = pole_id2;
                    max_value = max_value2;
                }
            }

            site_pole_radius_field[i] = max_value;
            site_density_field[i] = vertices_density_field[id];
        }

    }

    void run_graph_cut(const std::vector<double>& node_weight,
                       const std::vector<std::pair<int,int>>& edges,
                       const std::vector<double>& edge_weight,
                       std::vector<int>& labels,
                       double alpha) {

        int nb_node = node_weight.size();
        int nb_edges = edges.size();
        int nb_labels = 2;

        std::vector<double> data_cost(2 * nb_node);
        for(int i = 0; i < nb_node; i++) {
            data_cost[i] = node_weight[i] / (nb_node);
            data_cost[i + nb_node] = (1.0 - node_weight[i]) / (nb_node);
        }

        std::vector<double> final_edge_weight(nb_edges);
        double sum_edge_weight = 0;
        for(int i = 0; i < nb_edges; i++) {
            sum_edge_weight += edge_weight[i];
        }

        for(int i = 0; i < nb_edges; i++) {
            final_edge_weight[i] = alpha * edge_weight[i] / sum_edge_weight;
        }

        DataCost data_item(data_cost, nb_node, nb_labels);
        SmoothCost smooth_item(edges, final_edge_weight);

        GCoptimizationGeneralGraph gcp(nb_node, 2);
        gcp.setVerbosity(0);
        for(int i = 0; i < nb_edges; i++) {
            int x = edges[i].first;
            int y = edges[i].second;
            gcp.setNeighbors(x, y);
        }

        gcp.setDataCostFunctor(&data_item);
        gcp.setSmoothCostFunctor(&smooth_item);

        std::cout << "Before optimization energy is " << gcp.compute_energy() << std::endl;
        gcp.expansion(100);
        std::cout << "After optimization energy is " << gcp.compute_energy() << std::endl;

        labels.resize(nb_node);
        for(int i = 0; i < nb_node; i++) {
            if(gcp.whatLabel(i) == 0) {
                labels[i] = 0;
            } else {
                labels[i] = 1;
            }
        }
    }


    void find_island(const std::vector<std::vector<std::pair<int,double>>>& mesh_graph,
                          std::vector<double>& scalar_field,
                          std::vector<int>& island,
                          double scalar_thr,
                          double island_thr) {
        int nb_node = mesh_graph.size();
        std::vector<bool> vis(nb_node ,false);
        std::vector<int> cur_island;
        island.clear();

        std::function<void(int)> dfs1 = [&](int k) {
            vis[k] = true;
            cur_island.emplace_back(k);
            int num = mesh_graph[k].size();
            for(int i = 0; i < num; i++){
                int cur_node = mesh_graph[k][i].first;
                if(!vis[cur_node]) {
                    if(scalar_field[cur_node] < scalar_thr) {

                        dfs1(cur_node);
                    }
                }
            }
        };

        for(int i = 0; i < nb_node; i++) {
            if(!vis[i] && scalar_field[i] < scalar_thr) {
                cur_island.clear();
                dfs1(i);
                if(cur_island.size() < island_thr) {
                    island.insert(island.end(), cur_island.begin(), cur_island.end());
                }
            }
        }

        std::function<void(int)> dfs2 = [&](int k) {
            vis[k] = true;
            cur_island.emplace_back(k);
            int num = mesh_graph[k].size();
            for(int i = 0; i < num; i++){
                int cur_node = mesh_graph[k][i].first;
                if(!vis[cur_node]) {
                    if(scalar_field[cur_node] > scalar_thr) {

                        dfs2(cur_node);
                    }
                }
            }
        };


        vis.resize(nb_node, false);

        for(int i = 0; i < nb_node; i++) {
            if(!vis[i] && scalar_field[i] > scalar_thr) {
                cur_island.clear();
                dfs2(i);
                if(cur_island.size() < island_thr) {
                    island.insert(island.end(), cur_island.begin(), cur_island.end());
                }
            }
        }


    }

    void corrosion_island(const std::vector<Point>& points,
                          const std::vector<int> island,
                          double radius,
                          double scalar_thr,
                          std::vector<double>& scalar_field) {
        std::vector<KNN::Point> knn_points;
        cgal_points_convert_to_knn_points(points, knn_points);
        KNN::KdSearch kds(knn_points);
        int nb_island = island.size();
        for(int i = 0; i < nb_island; i++) {
            int id = island[i];
            std::vector<int> neighbors;
            std::vector<double> squared_distances;
            int num = kds.radius_search(knn_points[id], radius, neighbors, squared_distances);
            double val = scalar_field[id] > scalar_thr ? 1.0 : 0;
            double max_value = 0.0, min_value = 1.0;
            for(int i = 0; i < num; i++) {
                max_value = std::max(scalar_field[neighbors[i]], max_value);
                min_value = std::min(scalar_field[neighbors[i]], min_value);
            }
            scalar_field[id] = scalar_field[id] > scalar_thr ? min_value : max_value;
        }
    }

    void classify_defillet_triangle_via_vertices_voting(const std::vector<std::vector<size_t>>& f_ind,
                                                        const std::vector<double>& vertices_field,
                                                        double thr,
                                                        std::vector<bool>& is_fillet) {
        int nb_faces = f_ind.size();

        is_fillet.resize(nb_faces);
        for(int i = 0; i < nb_faces; i++) {
            int ct  = 0, num = f_ind[i].size();
            for(int j = 0 ; j < num; j++) {
                int id = f_ind[i][j];
                if(vertices_field[id] > thr) {
                    ct++;
                }
                else {
                    ct--;
                }
            }
            if(ct > 0) {
                is_fillet[i] = true;
            } else {
                is_fillet[i] = false;
            }
        }

    }

    void extract_interest_region(const std::vector<Point>& points,
                               const std::vector<std::vector<size_t>>& faces,
                               const std::vector<bool>& is_fillet,
                               std::vector<Point>& fillet_points,
                               std::vector<std::vector<size_t>>& fillet_faces,
                               std::vector<int>& point_map_table,
                               std::vector<int>& faces_map_table,
                               std::vector<bool>& face_labels) {

        Surface_mesh mesh;
        int nb_points = points.size(), nb_faces = faces.size();

        for(int i = 0; i < nb_points; i++) {
            mesh.add_vertex(points[i]);
        }

        for(int i = 0; i < nb_faces; i++) {
            mesh.add_face(Vertex_index(faces[i][0])
                    ,Vertex_index(faces[i][1]), Vertex_index(faces[i][2]));
        }

        std::set<int> interest_face;

        for(auto f : mesh.faces()) {
            if(is_fillet[f.idx()]) {
                for(auto af : mesh.faces_around_face(mesh.halfedge(f))) {
                    interest_face.insert(af.idx());
                }
                interest_face.insert(f.idx());
            }
        }

        std::set<int> interest_vertex;

        for(auto fi : interest_face) {
            auto f = Face_index (fi);
            for(auto fv : mesh.vertices_around_face(mesh.halfedge(f))) {
                interest_vertex.insert(fv.idx());
            }
        }

        fillet_points.clear();
        point_map_table.clear();
        std::map<int,int> point_map;
        for(auto vi : interest_vertex) {
            point_map[vi] = point_map_table.size();
            point_map_table.emplace_back(vi);
            fillet_points.emplace_back(points[vi]);
        }

        fillet_faces.clear();
        faces_map_table.clear();
        face_labels.clear();
        for(auto fi : interest_face) {
            std::vector<size_t> tmp;
            int num = faces[fi].size();
            for(int i = 0; i < num; i++) {
                tmp.emplace_back(point_map[faces[fi][i]]);
            }
            fillet_faces.emplace_back(tmp);
            faces_map_table.emplace_back(fi);
            face_labels.emplace_back(is_fillet[fi]);
        }


    }

    void run_geodesic(const std::vector<Point>& points,
                      const std::vector<std::vector<size_t>>& faces,
                      const std::vector<int>& centroid,
                      const std::vector<bool>& centroid_labels,
                      std::vector<int>& ancestor,
                      std::vector<double>& distance) {
        int nb_points = points.size();
        std::vector<CPoint3D> xin_points;
        for(int i = 0; i < nb_points; i++) {
            xin_points.emplace_back(CPoint3D(points[i].x(),
                                             points[i].y(), points[i].z()));
        }

        int nb_faces = faces.size();
        std::vector<CBaseModel::CFace> xin_faces;

        for(int i = 0; i < nb_faces; i++) {
            xin_faces.emplace_back(CBaseModel::CFace(faces[i][0],
                                                     faces[i][1], faces[i][2]));
        }

        CRichModel mesh(xin_points, xin_faces);
        std::set<int> sources;

        int nb_centroid = centroid.size();
        for(int i = 0; i < nb_centroid; i++) {
            if(!centroid_labels[i]) {
                sources.insert(centroid[i]);
            }
        }
        CXin_Wang alg(mesh, sources);
        alg.Execute();

//        std::cout << nb_points << ' ' << mesh.GetNumOfVerts() << std::endl;
        ancestor.resize(nb_points);
        distance.resize(nb_points);

        for(int i = 0; i < nb_points; i++) {
            ancestor[i] = alg.GetAncestor(i);
            distance[i] = alg.GetDistanceField()[i];
        }
    }


    void geodesic_post_processing(const std::vector<Point>& points,
                                  const std::vector<std::vector<size_t>>& faces,
                                  const std::vector<int>& ancestor,
                                  std::vector<Vector_3>& target_normals) {
        std::vector<int> new_ancestor = ancestor;
        int nb_points = points.size();
        int nb_faces = faces.size();

        for(int i = 0; i < nb_faces; i++) {
            new_ancestor[i] -= nb_points;
        }


        Surface_mesh mesh;
        for(int i = 0; i < nb_points; i++) {
            mesh.add_vertex(Point(points[i][0], points[i][1], points[i][2]));
        }

        for(int i = 0; i < nb_faces; i++) {
            mesh.add_face(Vertex_index(faces[i][0])
                    ,Vertex_index(faces[i][1]), Vertex_index(faces[i][2]));
        }
        typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
        auto fnormals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
        CGAL::Polygon_mesh_processing::compute_face_normals(mesh,  fnormals);

        target_normals.resize(nb_faces);
        for(int i = 0; i < nb_faces; i++) {
            int id = new_ancestor[i];
            target_normals[i] = fnormals[face_descriptor(id)];
        }
    }

    void optimize(const std::vector<Point>& points,
                  const std::vector<std::vector<size_t>>& faces,
                  std::vector<Vector_3>& normals,
                  std::vector<Point>& new_points) {

        int nb_faces = faces.size();
        int nb_points = points.size();

        Eigen::SparseMatrix<double> A(nb_faces * 2, nb_points * 3);
        std::cout << "ASD" << std::endl;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(nb_points * 3);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < nb_faces; i++) {
            double n0 = normals[i].x(), n1 = normals[i].y(), n2 = normals[i].z();
            int v0 = faces[i][0], v1 = faces[i][1], v2 = faces[i][2];
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0, n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 1, n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 2, n2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1, -n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 1, -n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 2, -n2));

            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0, n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 1, n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 2, n2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2, -n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 1, -n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 2, -n2));
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
        std::cout << "ASD" << std::endl;

        Eigen::VectorXd init(3 * nb_points);
        for(int i = 0; i < nb_points; i++) {
            init[3 * i] = points[i].x();
            init[3 * i + 1] = points[i].y();
            init[3 * i + 2] = points[i].z();
        }

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
        solver.setMaxIterations(10); // 设置最大迭代次数
        solver.setTolerance(1e-3);    // 设置迭代收敛容许误差
        std::cout << "solove" <<std::endl;
        Eigen::VectorXd x = solver.compute(A).solveWithGuess(b, init);
        std::cout << "solove" <<std::endl;
        if (solver.info() != Eigen::Success) {
            // 求解失败
            std::cout << "Eigen CG solver failed to converge." << std::endl;
            return;
        } else {
            std::cout << "Eigen CG solver success to converge." << std::endl;
        }
        new_points.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            new_points[i] = Point(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
        }
    }
}