//
// Created by 小乌嘎 on 2024/2/10.
//

#include "defillet.h"
#include "gcp.h"
#include "optimize.h"
#include <numeric>
#include <easy3d/core/surface_mesh.h>
#include <Xin_Wang.h>

//#include <Eigen/Sparse>
#include <Eigen/SparseLU>
//#include <unsupported/Eigen/KroneckerProduct>

#include <igl/cat.h>
#include "common.h"

#include <queue>
#include "voronoi.h"

namespace DEFILLET {

    void compute_fillet_field_v1(const std::vector<Eigen::Vector3d>& points,
                                 const std::vector<std::vector<size_t>>& faces,
                                 std::vector<double>& fillet_field,
                                 int nei) {

        Eigen::AlignedBox3d boundingBox;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            boundingBox.extend(points[i]);
        }
        double diagonal = boundingBox.diagonal().norm();

        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(easy3d::vec3(points[i].x(), points[i].y(), points[i].z()));
        }
        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            easy3d::SurfaceMesh::Vertex v0(faces[i][0]);
            easy3d::SurfaceMesh::Vertex v1(faces[i][1]);
            easy3d::SurfaceMesh::Vertex v2(faces[i][2]);
            mesh->add_triangle(v0, v1, v2);
        }
        fillet_field.resize(nb_points);
        for(auto cur_v : mesh->vertices()) {
            std::set<easy3d::SurfaceMesh::Vertex> k_near;
            std::set<easy3d::SurfaceMesh::Vertex> vis;
            std::priority_queue<std::pair<double, easy3d::SurfaceMesh::Vertex>> que;
            que.push(make_pair(0.0, cur_v));
            while((!que.empty()) && k_near.size() < nei) {
                auto v = que.top().second; que.pop();
                k_near.insert(v);
                for(auto vat : mesh->vertices(v)) {
                    if(vis.find(vat) == vis.end()) {
                        vis.insert(vat);
                        double len = (mesh->position(vat) - mesh->position(v)).norm();
                        que.push(std::make_pair(-len, vat));
                    }
                }

            }

            std::vector<CGAL_Point> sites;
            std::vector<int> sites_id;
            for(auto v : k_near) {
                sites.emplace_back(CGAL_Point(points[v.idx()].x(),
                                            points[v.idx()].y(),points[v.idx()].z()));
                sites_id.emplace_back(v.idx());
            }
            int nb_sites = sites.size();
            Voronoi vor(sites);
            vor.cal_v1();
            const std::vector<CGAL_Point>& vor_vertices = vor.get_vertices();
            std::vector<Eigen::Vector3d> eigen_vertices;
            UTILS::cgal_points_to_eigen_points(vor_vertices, eigen_vertices);
            const std::vector<std::vector<int>>& cell_pole = vor.get_cell_pole();
            int nb_vertices = vor_vertices.size();
            std::vector<double>tmp;
            for(int i = 0; i < nb_sites; i++) {
                int num = cell_pole[i].size();
                double minn = std::numeric_limits<double>::max();
                bool flag = false;
                for(int j = 0; j < num; j++) {
                    int id = cell_pole[i][j];
                    if(boundingBox.contains(eigen_vertices[id])) {
                        double len = (eigen_vertices[id] - points[sites_id[i]]).norm();
                        if(minn > len) {
                             minn = len; flag = true;
                        }
                    }
                }
                if(flag) {
                    tmp.emplace_back(minn);
                }
            }

            int k = tmp.size();
            if(k > 0) {
                double avg = 0;
                for (int i = 0; i < k; i++) {
                    avg += tmp[i];
                }
                avg /= k;
                double val = 0;
                for (int i = 0; i < k; i++) {
                    val += (tmp[i] - avg) * (tmp[i] - avg);
                }
                val /= k;
//                std::cout <<cur_v.idx() << ' '<<  k << ' ' << val <<' ' << nb_sites << std::endl;
                fillet_field[cur_v.idx()] = k * (1.0 - exp(-val)) / nb_sites;
            } else {
                fillet_field[cur_v.idx()] = 0.001;

            }
//            std::cout << cur_v.idx() << ' ' << fillet_field[cur_v.idx()] <<std::endl;
        }
    }

    void compute_fillet_field(const std::vector<Eigen::Vector3d>& points,
                              const std::vector<std::vector<size_t>>& faces,
                              const std::vector<Eigen::Vector3d>& vor_vertices,
                              const std::vector<std::vector<int>>& vor_site_poles,
                              std::vector<double>& density_field,
                              std::vector<double>& fillet_field,
                              std::vector<size_t>& fillet_poles,
                              double eps,
                              double radius) {

        Eigen::AlignedBox3d boundingBox;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            boundingBox.extend(points[i]);
        }
        double diagonal = boundingBox.diagonal().norm();
        if(radius < 0.0)
            radius = diagonal * eps;
        std::cout << "eps = " << eps << ", " << "diagonal = " << diagonal << ", " << "radius = " <<  radius << std::endl;
        int nb_vertices = vor_vertices.size();
        std::vector<KNN::Point> knn_vertices;
        UTILS::eigen_points_to_knn_points(vor_vertices, knn_vertices);
        KNN::KdSearch kds(knn_vertices);
        density_field.resize(nb_vertices);
        for(int i = 0; i < nb_vertices; i++) {
            std::vector<size_t> neighbors;
            std::vector<double> squared_distances;
            int num = kds.radius_search(knn_vertices[i], radius, neighbors, squared_distances);
            density_field[i] = num;
        }
        UTILS::vector_normalize(density_field);
//        double max_len = std::numeric_limits<double>::min();
//        double min_len = std::numeric_limits<double>::max();
//        for(int i = 0; i < nb_points; i++) {
//            for(size_t j = 0; j < vor_site_poles[i].size(); j++) {
//                int id = vor_site_poles[i][j];
//                double len = (vor_vertices[id] - points[i]).norm();
//                max_len = std::max(max_len, len);
//                min_len = std::min(min_len, len);
//            }
//        }
//        std::cout << max_len << ' ' << min_len << std::endl;
        fillet_field.resize(nb_points);
        fillet_poles.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            double maxx = 0;
            size_t max_id = vor_site_poles[i][0];
            for(size_t j = 0 ; j < vor_site_poles[i].size(); j++) {
                int id = vor_site_poles[i][j];
                double len = (vor_vertices[id] - points[i]).norm();
                double val = density_field[id];
                if(maxx < val) {
                    maxx = val;
                    max_id = id;
                }
            }
            fillet_field[i] = maxx;
            fillet_field[i] = 1.0;
            fillet_poles[i] = max_id;
        }

        double minn = *std::min_element(fillet_field.begin(),fillet_field.end());

        UTILS::vector_normalize(fillet_field);
    }

    void run_graph_cut(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       const std::vector<Eigen::Vector3d>& vor_vertices,
                       const std::vector<size_t>& fillet_poles,
                       std::vector<double>& fillet_field,
                       std::vector<double>& labels,
                       double alpha,
                       double thr) {
        UTILS::vector_normalize(fillet_field);
        if(thr < 0.0) {
            thr = UTILS::threshold_otsu(fillet_field, 255);
            std::cout << "thr = " << thr << std::endl;
        }
//        thr = 0.24;
        int nb_points = points.size();
        std::vector<double> data_cost(2 * nb_points);

        std::vector<std::vector<std::pair<int,double>>> graph;
        UTILS::extract_vertex_graph_from_mesh(points, faces, graph);
        for(int i = 0; i < nb_points; i++) {
            std::set<size_t> neighbors;
            UTILS::extract_k_ring_neighboring(graph, i, 1, neighbors);
            int num = neighbors.size();
            int ct = 0;
            for(auto j : neighbors) {
                if(fillet_field[j] > thr)
                    ct++;
            }
            data_cost[i] = 1.0 * ct / (2 * nb_points);
            data_cost[i + nb_points] = 1.0 * (num - ct) / (2 * nb_points);
        }

        std::vector<std::pair<int,int>> edges;
        std::vector<double> edge_weights;
        UTILS::extract_edge_from_mesh(points, faces, edges, edge_weights);
        int nb_edges = edges.size();
        double edges_sum = std::accumulate(edge_weights.begin(), edge_weights.end(), 0.0f);
//        std::cout << edges_sum << std::endl;
        std::vector<double> smooth_cost(nb_edges);

        for(int i = 0; i < nb_edges; i++) {
            double val = alpha * edge_weights[i] / edges_sum;
            smooth_cost[i] = val;
        }

        GCP::DataCost data_item(data_cost, nb_points, 2);
        GCP::SmoothCost smooth_item(edges, smooth_cost);

        GCoptimizationGeneralGraph gcp(nb_points, 2);
        gcp.setVerbosity(0);
        for(int i = 0; i < nb_edges; i++) {
            int x = edges[i].first;
            int y = edges[i].second;
            gcp.setNeighbors(x, y);
        }

        gcp.setDataCostFunctor(&data_item);
        gcp.setSmoothCostFunctor(&smooth_item);

        std::cout << "Before optimization energy is " << gcp.compute_energy() << std::endl;
        gcp.expansion(10);
        std::cout << "After optimization energy is " << gcp.compute_energy() << std::endl;
        std::cout << gcp.giveDataEnergy() << ' ' << gcp.giveSmoothEnergy() << std::endl;
        labels.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            if(gcp.whatLabel(i) == 0) {
                labels[i] = 0.0;
            } else {
                labels[i] = 1.0;
            }
        }
    }

    void run_graph_cut_v1(std::vector<Eigen::Vector3d>& points,
                       std::vector<std::vector<size_t>>& faces,
                       std::vector<double>& fillet_field,
                       std::vector<double>& labels,
                       double alpha) {
        std::vector<std::pair<int,int>> edges;
        std::vector<double> edge_weights1;
        UTILS::extract_edge_from_mesh(points, faces, edges, edge_weights1);
        UTILS::vector_normalize(fillet_field);
        UTILS::vector_normalize(edge_weights1);
        double thr = UTILS::threshold_otsu(fillet_field, 100);
//        thr *= 0.9;
        std::cout << "thr = " << thr << std::endl;
        int nb_points = points.size();
        std::vector<double> data_cost(2 * nb_points);
        for(int i = 0; i < nb_points; i++) {
            if(fillet_field[i] > thr) {
                data_cost[i] = 1.0;
                data_cost[i + nb_points] = 0.0;
//                data_cost[i] = 0.5 + 0.5 * (fillet_field[i] - thr) / (1 - thr);
//                data_cost[i + nb_points] = 1.0 - data_cost[i];
            } else {
                data_cost[i] = 0;
                data_cost[i + nb_points] = 1.0;
//                data_cost[i] = 0.5 * fillet_field[i] / thr;
//                data_cost[i + nb_points] = 1.0 - data_cost[i];
            }
        }
        int nb_edges = edges.size();
        std::vector<double> edge_weights2(nb_edges);
        for(int i = 0; i < nb_edges; i++) {
            double len = (points[edges[i].first] - points[edges[i].second]).norm();
            edge_weights2[i] = len * alpha;
        }
        UTILS::vector_normalize(edge_weights2);
        for(int i = 0; i < nb_edges; i++) {
            edge_weights1[i] *= alpha;
        }
        std::vector<double> smooth_cost(nb_edges);
        for(int i = 0; i < nb_edges; i++) {
            smooth_cost[i] = edge_weights1[i];
//            smooth_cost[i] = (edge_weights1[i] + edge_weights2[i]) / 2;
        }
        GCP::DataCost data_item(data_cost, nb_points, 2);
        GCP::SmoothCost smooth_item(edges, smooth_cost);

        GCoptimizationGeneralGraph gcp(nb_points, 2);
        gcp.setVerbosity(0);
        for(int i = 0; i < nb_edges; i++) {
            int x = edges[i].first;
            int y = edges[i].second;
            gcp.setNeighbors(x, y);
        }

        gcp.setDataCostFunctor(&data_item);
        gcp.setSmoothCostFunctor(&smooth_item);

        std::cout << "Before optimization energy is " << gcp.compute_energy() << std::endl;
        gcp.expansion(10);
        std::cout << "After optimization energy is " << gcp.compute_energy() << std::endl;
        std::cout << gcp.giveDataEnergy() << ' ' << gcp.giveSmoothEnergy() << std::endl;
        labels.resize(nb_points);
        for(int i = 0; i < nb_points; i++) {
            if(gcp.whatLabel(i) == 0) {
                labels[i] = 0.1;
            } else {
                labels[i] = 1.0;
            }
        }
    }

    void extract_interest_region(const std::vector<Eigen::Vector3d>& points,
                                 const std::vector<std::vector<size_t>>& faces,
                                 const std::vector<double>& label,
                                 std::vector<Eigen::Vector3d>& fillet_points,
                                 std::vector<std::vector<size_t>>& fillet_faces,
                                 std::vector<int>& point_map,
                                 std::vector<std::vector<size_t>>& non_fillet_faces,
                                 std::vector<size_t>& fillet_bounding,
                                 std::vector<Eigen::Vector3d>& fillet_bounding_normals) {
        easy3d::SurfaceMesh* mesh = new easy3d::SurfaceMesh;
        int nb_points = points.size();
        for(int i = 0; i < nb_points; i++) {
            mesh->add_vertex(easy3d::vec3(points[i].x(),
                                          points[i].y(), points[i].z()));
        }

        int nb_faces = faces.size();
        for(int i = 0; i < nb_faces; i++) {
            auto v1 = easy3d::SurfaceMesh::Vertex(faces[i][0]);
            auto v2 = easy3d::SurfaceMesh::Vertex(faces[i][1]);
            auto v3 = easy3d::SurfaceMesh::Vertex(faces[i][2]);
            mesh->add_triangle(v1, v2, v3);
        }

        std::set<easy3d::SurfaceMesh::Face> faces_set;
        std::set<easy3d::SurfaceMesh::Vertex> vertex_set;
        std::set<easy3d::SurfaceMesh::Vertex> bounding_set;
        for(auto f : mesh->faces()) {
            int flag = 0;
            for(auto v : mesh->vertices(f)) {
                flag |= int(label[v.idx()] + 0.1);
            }
            if(flag) {
                faces_set.insert(f);
                for(auto v : mesh->vertices(f)) {
                    vertex_set.insert(v);
                    flag = int(label[v.idx()] + 0.1);
                    if(flag == 0) {
                        bounding_set.insert(v);
                    }
                }
            } else {
                non_fillet_faces.emplace_back(faces[f.idx()]);
            }
        }

        point_map.resize(vertex_set.size());
        std::map<int, size_t> mp;
        size_t id = 0;
        for(auto v : vertex_set) {
            point_map[id] = v.idx();
            mp[v.idx()] = id++;
            easy3d::vec3 p = mesh->position(v);
            fillet_points.emplace_back(Eigen::Vector3d(p.x, p.y, p.z));
        }

        for(auto f : faces_set) {
            std::vector<size_t> tmp;
            for(auto v : mesh->vertices(f)) {
                tmp.emplace_back(mp[v.idx()]);
            }
            fillet_faces.emplace_back(tmp);
        }

        for(auto v : bounding_set) {
            int id = mp[v.idx()];
            easy3d::vec3 n = mesh->compute_vertex_normal(v);
            fillet_bounding.emplace_back(id);
            fillet_bounding_normals.emplace_back(Eigen::Vector3d(n.x, n.y, n.z));
        }
    }

    void run_geodesic(const std::vector<Eigen::Vector3d>& points,
                      const std::vector<std::vector<size_t>>& faces,
                      const std::vector<size_t>& sources,
                      const std::vector<Eigen::Vector3d>& sources_normals,
                      std::vector<size_t>& point_ancestors,
                      std::vector<size_t>& face_ancestors,
                      std::vector<double>& distance,
                      std::vector<Eigen::Vector3d>& point_tar_normals,
                      std::vector<Eigen::Vector3d>& face_tar_normals) {

        int nb_points = points.size();
        int nb_faces = faces.size();

        std::vector<CPoint3D> xin_points;
        for(int i = 0; i < nb_points; i++) {
            xin_points.emplace_back(CPoint3D(points[i].x(),
                                             points[i].y(), points[i].z()));
        }
        std::vector<CBaseModel::CFace> xin_faces;
        for(int i = 0; i < nb_faces; i++) {
            int num = faces[i].size();
            Eigen::Vector3d center; center.setZero();
            for(int j = 0; j < num; j++) {
                center += points[faces[i][j]];
            }
            center /= num;
            xin_points.emplace_back(CPoint3D(center.x(),
                                             center.y(), center.z()));
            for(int j = 0; j < num; j++) {
                xin_faces.emplace_back(CBaseModel::CFace(faces[i][j],
                                                         faces[i][(j + 1) % num], nb_points + i));
            }
        }

        CRichModel mesh(xin_points, xin_faces);
        std::set<int> xin_sources(sources.begin(), sources.end());
        CXin_Wang alg(mesh, xin_sources);
        alg.Execute();
        point_ancestors.resize(nb_points);
        face_ancestors.resize(nb_faces);
        distance.resize(nb_points);
        point_tar_normals.resize(nb_points);
        face_tar_normals.resize(nb_faces);

        int nb_sources = sources.size();
        std::map<int,int> mp;
        for(int i = 0; i < nb_sources; i++) {
            mp[sources[i]] = i;
        }
        for(int i = 0; i < nb_points; i++) {
            point_ancestors[i] = alg.GetAncestor(i);
            point_tar_normals[i] = sources_normals[mp[point_ancestors[i]]];
            distance[i] = alg.GetDistanceField()[i];
        }

        for(int i = 0; i < nb_faces; i++) {
            face_ancestors[i] = alg.GetAncestor(i + nb_points);
            face_tar_normals[i] = sources_normals[mp[face_ancestors[i]]];
        }
    }

    bool optimize_qem(const std::vector<Eigen::Vector3d>& points,
                  const std::vector<std::pair<int,int>>& edges,
                  std::vector<Eigen::Vector3d>& normals,
                  std::vector<Eigen::Vector3d>& cur_normals,
                  std::vector<Eigen::Vector3d>& new_points,
                  std::vector<size_t>& fixed_points,
                  double beta) {
        int nb_points = points.size(), nb_edges = edges.size();
        Eigen::SparseMatrix<double> NC(nb_edges * 2, nb_points * 3);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < nb_edges; i++) {

            int v0 = edges[i].first, v1 = edges[i].second;
            double nx1 = normals[v0].x(), ny1 = normals[v0].y(), nz1 = normals[v0].z();
            double nx2 = cur_normals[v0].x(), ny2 = cur_normals[v0].y(), nz2 = cur_normals[v0].z();
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0, nx1 + nx2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points, ny1 + ny2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points * 2, nz1 + nz2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1, -(nx1 + nx2)));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points, -(ny1 + ny2)));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points * 2, -(nz1 + nz2)));

            nx1 = normals[v1].x(), ny1 = normals[v1].y(), nz1 = normals[v1].z();
            nx2 = cur_normals[v1].x(), ny2 = cur_normals[v1].y(), nz2 = cur_normals[v1].z();
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0, nx1 + nx2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points, ny1 + ny2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points * 2, nz1 + nz2));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1, -(nx1 + nx2)));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1 + nb_points, -(ny1 + ny2)));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1 + nb_points * 2, -(nz1 + nz2)));
        }
        NC.setFromTriplets(triplets.begin(), triplets.end());

        int nb_fixed_points = fixed_points.size();
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

        Eigen::SparseMatrix<double> A = NC;
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> ATA = AT * A;
        Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
        I = beta * I;
        Eigen::SparseMatrix<double> Q = I + ATA;
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

        Eigen::VectorXd b(nb_points * 3 + d.size());
        b << beta * p, d;

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

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
            new_points[i] = Eigen::Vector3d(x[i], x[i + nb_points], x[i + 2 * nb_points]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        return true;
    }

    bool optimize_centroid(const std::vector<Eigen::Vector3d>& points,
                      const std::vector<std::pair<int,int>>& edges,
                      std::vector<Eigen::Vector3d>& normals,
                      std::vector<Eigen::Vector3d>& cur_normals,
                      std::vector<Eigen::Vector3d>& new_points,
                      std::vector<size_t>& fixed_points,
                      double beta) {
        int nb_points = points.size(), nb_edges = edges.size();
        Eigen::SparseMatrix<double> NC(nb_edges * 2, nb_points * 3);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < nb_edges; i++) {
            int v0 = edges[i].first, v1 = edges[i].second;
            double nx1 = normals[v0].x(), ny1 = normals[v0].y(), nz1 = normals[v0].z();
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0, nx1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points, ny1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v0 + nb_points * 2, nz1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1, -nx1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points, -ny1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i, v1 + nb_points * 2, -nz1));

            nx1 = normals[v1].x(), ny1 = normals[v1].y(), nz1 = normals[v1].z();
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0, nx1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points, ny1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v0 + nb_points * 2, nz1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1, -nx1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1 + nb_points, -ny1));
            triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, v1 + nb_points * 2, -nz1));
        }
        NC.setFromTriplets(triplets.begin(), triplets.end());

        std::vector<std::vector<int>> graph(nb_points);
        for(int i = 0; i < nb_edges; i++) {
            int x = edges[i].first, y = edges[i].second;
            graph[x].emplace_back(y);
            graph[y].emplace_back(x);
        }
        triplets.clear();
        Eigen::SparseMatrix<double> CC(nb_points, nb_points * 3);
        for(int i = 0; i < nb_points; i++) {
            int num = graph[i].size();
            for(int j = 0; j < num; j++) {
                int id = graph[i][j];
                triplets.emplace_back(Eigen::Triplet<double>(i, id, 1.0));
                triplets.emplace_back(Eigen::Triplet<double>(i, id + nb_points, 1.0));
                triplets.emplace_back(Eigen::Triplet<double>(i, id + 2 * nb_points, 1.0));
            }
            triplets.emplace_back(Eigen::Triplet<double>(i, i, -num));
            triplets.emplace_back(Eigen::Triplet<double>(i, i + nb_points, -num));
            triplets.emplace_back(Eigen::Triplet<double>(i, i + 2 * nb_points, -num));
        }
        CC.setFromTriplets(triplets.begin(), triplets.end());

        int nb_fixed_points = fixed_points.size();
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

        Eigen::SparseMatrix<double> A;
        igl::cat(1, NC, CC, A);
        Eigen::SparseMatrix<double> AT = A.transpose();
        Eigen::SparseMatrix<double> ATA = AT * A;
        Eigen::SparseMatrix<double> I(ATA.rows(), ATA.cols()); I.setIdentity();
        I = beta * I;
        Eigen::SparseMatrix<double> Q = I + ATA;
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
//        p = AT * p;
        Eigen::VectorXd b(nb_points * 3 + d.size());
        b << beta * p, d;

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

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
            new_points[i] = Eigen::Vector3d(x[i], x[i + nb_points], x[i + 2 * nb_points]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
        }
        return true;
    }

    bool iterative_optimize(const std::vector<Eigen::Vector3d>& points,
                            const std::vector<std::vector<size_t>>& faces,
                            const std::vector<size_t>& point_ancestors,
                            const std::vector<size_t>& face_ancestors,
                            const std::vector<Eigen::Vector3d>& point_tar_normals,
                            const std::vector<Eigen::Vector3d>& face_tar_normals,
                            const std::vector<size_t>& fixed_points,
                            std::vector<Eigen::Vector3d>& new_points,
                            std::vector<std::vector<size_t>>& new_faces,
                            std::string type,
                            double beta,
                            int num_iterations) {

        Optimize opt(points, faces, point_ancestors, face_ancestors, point_tar_normals, face_tar_normals, fixed_points, type, beta);

        for(int i = 0; i < num_iterations; i++) {

            std::cout << i << std::endl;
            if(!opt.solve()) {
                return false;
            }

        }
//        opt.remesh();
//        opt.refine_tar_normals();
//        opt.solve();
        opt.get_points(new_points);
        opt.get_faces(new_faces);
        return true;
    }
}