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
#include <Eigen/SVD>
#include <igl/cat.h>
#include <igl/min_quad_with_fixed.h>
#include <Eigen/SparseLU>
#include <algorithm>
#include <easy3d/core/surface_mesh.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <unsupported/Eigen/KroneckerProduct>


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
                                  const std::vector<bool>& face_labels,
                                  std::vector<Vector_3>& target_normals,
                                  std::vector<int>& fixed_points,
                                  std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector,
                                  double angle_thr) {

        int nb_points = points.size();
        int nb_faces = faces.size();
        std::vector<int> new_ancestor(ancestor.begin() + points.size(), ancestor.end());

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

        edge_vector.clear();
        std::set<int> fixed_faces_set;
        for (const Edge_index& edge : mesh.edges()) {

            // 获取边的两个半边
            Halfedge_index he1 = mesh.halfedge(edge, 0);
            Halfedge_index he2 = mesh.halfedge(edge, 1);
//            std::cout << "ASD1" << std::endl;
            Face_index f1 = mesh.face(he1);
            Face_index f2 = mesh.face(he2);
            if(f1.idx() >= nb_faces || f2.idx() >= nb_faces) {
//                if(f1.idx() < nb_faces) {
//                    fixed_faces_set.insert(f1.idx());
//                }
//                if(f2.idx() < nb_faces) {
//                    fixed_faces_set.insert(f2.idx());
//                }
                continue;
            }
            Eigen::Vector3d tnf1 = Eigen::Vector3d(target_normals[f1.idx()].x(),
                                                  target_normals[f1.idx()].y(), target_normals[f1.idx()].z()).normalized();
            Eigen::Vector3d tnf2 = Eigen::Vector3d(target_normals[f2.idx()].x(),
                                                  target_normals[f2.idx()].y(), target_normals[f2.idx()].z()).normalized();

            double dot_product = tnf1.dot(tnf2);

            // 计算夹角（弧度）

            double angle_radians = std::acos(dot_product);
//            std::cout << angle_radians << std::endl;
            // 将弧度转换为角度
            double angle_degrees = angle_radians * 180.0 / M_PI;

            if(angle_degrees < angle_thr) {
                Eigen::Vector3d onf1 = Eigen::Vector3d(fnormals[f1].x(),
                                                       fnormals[f1].y(), fnormals[f1].z()).normalized();

                int id1 = mesh.source(he1).idx();
                int id2 = mesh.target(he1).idx();

                // 计算旋转四元数
                Eigen::Quaterniond quaternion;
                quaternion.setFromTwoVectors(onf1, tnf1);

                // 将四元数转换为旋转矩阵
                Eigen::Matrix3d R = quaternion.toRotationMatrix();
                Vector_3 cgal_v = points[id1] - points[id2];
                Eigen::Vector3d v = Eigen::Vector3d(cgal_v.x(), cgal_v.y(), cgal_v.z()).normalized();
                v = R * v;

                edge_vector.emplace_back(std::make_pair(id1, std::make_pair(id2, Vector_3(v.x(), v.y(), v.z()))));
            }


        }

        fixed_points.clear();
        std::set<int> fixed_points_set;
        for(int i = 0; i < nb_faces; i++) {
            if(!face_labels[i]) {
                int num = faces[i].size();
                for (int j = 0; j < num; j++) {
                    fixed_points_set.insert(faces[i][j]);
                }
            }
        }
        fixed_points.insert(fixed_points.end(), fixed_points_set.begin(), fixed_points_set.end());
    }

    bool optimize_with_fixed_solve(const std::vector<Point>& points,
                  const std::vector<std::vector<size_t>>& faces,
                  std::vector<Vector_3>& normals,
                  std::vector<Point>& new_points,
                  std::vector<int>& fixed_points,
                  std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector) {
        int nb_faces = faces.size();
        int nb_points = points.size();
        int nb_fixed_points = fixed_points.size();
        int nb_edges = edge_vector.size();
        std::cout << "run optimize..." << std::endl;
        // normals constraint
        std::cout << "add normals constraint" << std::endl;
        Eigen::SparseMatrix<double> NC(nb_faces * 2, nb_points * 3);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(nb_faces * 2 + nb_edges);
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
        NC.setFromTriplets(triplets.begin(), triplets.end());

        // edges vector constraint
        std::cout << "add edges vector constraint" << std::endl;

        Eigen::SparseMatrix<double> EVC(nb_edges, nb_points * 3);
        triplets.clear();
        for(int i = 0; i < nb_edges; i++) {
            int id1 = edge_vector[i].first;
            int id2 = edge_vector[i].second.first;
            const Vector_3& v = edge_vector[i].second.second;
            double c1 = v.x(), c2 = v.y(), c3 = v.z();
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1, (-c2 + c3)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2, (c2 - c3)));

            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 1, (c1 - c3)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 1, (-c1 + c3)));

            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 2, (-c1 + c2)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 2, (c1 - c2)));
        }
        EVC.setFromTriplets(triplets.begin(), triplets.end());


        // fixed points constraint
        std::cout << "add fixed points constraint" << std::endl;
        Eigen::VectorXd d = Eigen::VectorXd(nb_fixed_points * 3);
        std::cout << nb_fixed_points << std::endl;
        triplets.clear();
        Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
        std::cout << nb_points << std::endl;
        for(int i = 0; i < nb_fixed_points; i++) {
            int id = fixed_points[i];
            triplets.emplace_back(Eigen::Triplet<double>(3 * i,  3 * id, 1.0)); d[ 3 * i] = points[id].x();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1,  3 * id + 1, 1.0)); d[3 * i + 1] = points[id].y();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  3 * id + 2, 1.0)); d[3 * i + 2] = points[id].z();
        }
        E.setFromTriplets(triplets.begin(), triplets.end());

        std::cout << "run solver..." << std::endl;


        Eigen::SparseMatrix<double> A;
        igl::cat(1, NC, EVC, A);
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> Q = AT * A;
        Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
        Eigen::SparseMatrix<double> ET = E.transpose();
        Eigen::SparseMatrix<double> tempMat1;
        Eigen::SparseMatrix<double> tempMat2;
        Eigen::SparseMatrix<double> L;
        igl::cat(1, Q, E, tempMat1);
        igl::cat(1, ET, zero, tempMat2);
        igl::cat(0, tempMat1, tempMat2, L);
        Eigen::VectorXd b(nb_points * 3 + d.size());
        b << Eigen::VectorXd::Zero(nb_points * 3), d;

        std::cout << L.rows() << ' ' << L.cols() << std::endl;
        igl::min_quad_with_fixed_data<double> solver_data;
        igl::min_quad_with_fixed_precompute(L, Eigen::VectorXi(), Eigen::SparseMatrix<double>(), true, solver_data);
        Eigen::VectorXd x(b.size(), 0);
        for(int i = 0; i < nb_points; i++) {
            x[3 * i] = points[i].x();
            x[3 * i + 1] = points[i].y();
            x[3 * i + 2] = points[i].z();
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        Eigen::MatrixXd Y(0, 3), Beq(0, 3);
        igl::min_quad_with_fixed_solve(solver_data, b, Y, Beq, x);
        new_points.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            new_points[i] = Point(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        std::cout << "done." << std::endl;
    }

    bool optimize_sparseLU(const std::vector<Point>& points,
                           const std::vector<std::vector<size_t>>& faces,
                           std::vector<Vector_3>& normals,
                           std::vector<Point>& new_points,
                           std::vector<int>& fixed_points,
                           std::vector<std::pair<int,std::pair<int,Vector_3>>>& edge_vector) {
        int nb_faces = faces.size();
        int nb_points = points.size();
        int nb_fixed_points = fixed_points.size();
        int nb_edges = edge_vector.size();
        std::cout << "run optimize..." << std::endl;
        // normals constraint
        std::cout << "add normals constraint" << std::endl;
        Eigen::SparseMatrix<double> NC(nb_faces * 2, nb_points * 3);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(nb_faces * 2 + nb_edges);
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
        NC.setFromTriplets(triplets.begin(), triplets.end());

        // edges vector constraint
        std::cout << "add edges vector constraint" << std::endl;

        Eigen::SparseMatrix<double> EVC(nb_edges, nb_points * 3);
        triplets.clear();
        for(int i = 0; i < nb_edges; i++) {
            int id1 = edge_vector[i].first;
            int id2 = edge_vector[i].second.first;
            const Vector_3& v = edge_vector[i].second.second;
            double c1 = v.x(), c2 = v.y(), c3 = v.z();
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1, (-c2 + c3)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2, (c2 - c3)));

            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 1, (c1 - c3)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 1, (-c1 + c3)));

            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 2, (-c1 + c2)));
            triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 2, (c1 - c2)));
        }
        EVC.setFromTriplets(triplets.begin(), triplets.end());


        // fixed points constraint
        std::cout << "add fixed points constraint" << std::endl;
        Eigen::VectorXd d = Eigen::VectorXd(nb_fixed_points * 3);
        std::cout << nb_fixed_points << std::endl;
        triplets.clear();
        Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
        std::cout << nb_points << std::endl;
        for(int i = 0; i < nb_fixed_points; i++) {
            int id = fixed_points[i];
            triplets.emplace_back(Eigen::Triplet<double>(3 * i,  3 * id, 1.0)); d[ 3 * i] = points[id].x();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1,  3 * id + 1, 1.0)); d[3 * i + 1] = points[id].y();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  3 * id + 2, 1.0)); d[3 * i + 2] = points[id].z();
        }
        E.setFromTriplets(triplets.begin(), triplets.end());

        std::cout << "run solver..." << std::endl;


        Eigen::SparseMatrix<double> A;
        igl::cat(1, NC, EVC, A);
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> Q = AT * A;
        Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
        Eigen::SparseMatrix<double> ET = E.transpose();
        Eigen::SparseMatrix<double> tempMat1;
        Eigen::SparseMatrix<double> tempMat2;
        Eigen::SparseMatrix<double> L;
        igl::cat(1, Q, E, tempMat1);
        igl::cat(1, ET, zero, tempMat2);
        igl::cat(0, tempMat1, tempMat2, L);
        Eigen::VectorXd b(nb_points * 3 + d.size());
        b << Eigen::VectorXd::Zero(nb_points * 3), d;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        std::cout << L.rows() << ' ' << L.cols() << std::endl;
        solver.compute(L);
        if(solver.info()!= Eigen::Success) {
            std::cout << "decomposition failed" << std::endl;
            return false;
        }
        Eigen::VectorXd x = solver.solve(b);
        if(solver.info()!= Eigen::Success) {
            // solving failed
            std::cout << "solving failed" << std::endl;
            return false;
        }
        new_points.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            new_points[i] = Point(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        std::cout << "done." << std::endl;
    }

    bool optimize_spareQR(const std::vector<Point>& points,
                          const std::vector<std::vector<size_t>>& faces,
                          std::vector<Vector_3>& normals,
                          std::vector<Point>& new_points,
                          std::vector<int>& fixed_points,
                          double beta) {
        int nb_faces = faces.size();
        int nb_points = points.size();
        int nb_fixed_points = fixed_points.size();
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        V.resize(nb_points, 3);
        for(int i = 0; i < nb_points; i++) {
            V(i, 0) = points[i].x();V(i, 1) = points[i].y();V(i, 2) = points[i].z();
        }
        F.resize(nb_faces, 3);
        for(int i = 0; i < nb_faces; i++) {
            F(i, 0) = faces[i][0]; F(i, 1) = faces[i][1]; F(i, 2) = faces[i][2];
        }

        Eigen::SparseMatrix<double> C,M,Minv,L,LC;
        igl::cotmatrix(V,F, C);
        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
        igl::invert_diag(M,Minv);
        L = Minv * C;

        LC = Eigen::KroneckerProductSparse<Eigen::MatrixXd, Eigen::SparseMatrix<double>>(Eigen::MatrixXd::Identity(3, 3), L);


        Eigen::SparseMatrix<double> NC(nb_faces * 2, nb_points * 3);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(nb_faces * 2);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < nb_faces; i++) {
            double n0 = normals[i].x(), n1 = normals[i].y(), n2 = normals[i].z();
            int v0 = faces[i][0], v1 = faces[i][1], v2 = faces[i][2];
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0, n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points, n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points * 2, n2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1, -n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points, -n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points * 2, -n2));

            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0, n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points, n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points * 2, n2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v2, -n0));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v2 + nb_points, -n1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v2 + nb_points * 2, -n2));
        }
        NC.setFromTriplets(triplets.begin(), triplets.end());



        // fixed points constraint
        std::cout << "add fixed points constraint" << std::endl;
        Eigen::VectorXd d = Eigen::VectorXd(nb_fixed_points * 3);
        triplets.clear();
        Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
        for(int i = 0; i < nb_fixed_points; i++) {
            int id = fixed_points[i];
            triplets.emplace_back(Eigen::Triplet<double>(3 * i,  id, 1.0)); d[ 3 * i] = points[id].x();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1,  id + nb_points, 1.0)); d[3 * i + 1] = points[id].y();
            triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  id + nb_points * 2, 1.0)); d[3 * i + 2] = points[id].z();
        }
        E.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SparseMatrix<double> LTL = beta * LC.transpose() * LC;


        Eigen::SparseMatrix<double> A = NC;
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> ATA = AT * A;
        Eigen::SparseMatrix<double> Q = LTL + ATA;
        Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
        Eigen::SparseMatrix<double> ET = E.transpose();
        Eigen::SparseMatrix<double> tempMat1;
        Eigen::SparseMatrix<double> tempMat2;
        Eigen::SparseMatrix<double> S;
        igl::cat(1, Q, E, tempMat1);
        igl::cat(1, ET, zero, tempMat2);
        igl::cat(0, tempMat1, tempMat2, S);
        Eigen::VectorXd ab = Eigen::VectorXd::Zero(nb_points * 3);
        Eigen::VectorXd p(nb_points * 3);
        for(int i = 0; i < nb_points; i++) {
            p[i] = points[i].x();
            p[i + nb_points] = points[i].y();
            p[i + 2 * nb_points] = points[i].z();
        }
        Eigen::VectorXd lb = LTL  * p;
        Eigen::VectorXd b(nb_points * 3 + d.size());
        b << ab + lb, d;
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver; // 创建SparseQR求解器
        std::cout << S.rows() << ' ' << S.cols() << std::endl;
        solver.compute(S);
        if(solver.info()!= Eigen::Success) {
            std::cout << "decomposition failed" << std::endl;
            return false;
        }
        Eigen::VectorXd x = solver.solve(b);
        if(solver.info()!= Eigen::Success) {
            // solving failed
            std::cout << "solving failed" << std::endl;
            return false;
        }
        new_points.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            new_points[i] = Point(x[i], x[i + nb_points], x[i + 2 * nb_points]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        std::cout << "done." << std::endl;
        return true;
    }

    bool iterative_optimize(const std::vector<Point>& points,
                            const std::vector<std::vector<size_t>>& faces,
                            std::vector<Vector_3>& normals,
                            std::vector<Point>& new_points,
                            std::vector<int>& fixed_points,
                            double beta,
                            int nb_optimize) {
        std::vector<Point> old_points = points;
        double zz = beta;
        for(int i = 0; i < nb_optimize; i++) {

            optimize_spareQR(old_points, faces, normals, new_points, fixed_points, zz);
            old_points = new_points;
            zz /= 10;
//            new_points.clear();
        }
    }
}