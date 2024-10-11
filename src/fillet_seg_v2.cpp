//
// Created by xiaowuga on 2024/10/9.
//

#include "fillet_seg_v2.h"
#include "knn.h"

#include <omp.h>
#include <numeric>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>

void FilletSegV2::voronoi3d(const std::vector<easy3d::vec3>& sites
                 , std::vector<easy3d::vec3>& vor_vertices
                 , std::vector<std::vector<int>>& sites_nearby_vertices
                 , std::vector<std::vector<int>>& vertices_nearby_sites) {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
    typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
    typedef Triangulation::Vertex_handle  Vertex_handle;
    typedef Triangulation::Cell_handle  Cell_handle;
    typedef K::Point_3          CGAL_Point;
    typedef K::Vector_3         CGAL_Vector_3;
    typedef K::Aff_transformation_3 Transformation;
    typedef CGAL::Surface_mesh<CGAL_Point> Surface_mesh;
    typedef Surface_mesh::Vertex_index Vertex_index;
    typedef Surface_mesh::Face_index Face_index;
    typedef Surface_mesh::Edge_index Edge_index;
    typedef Surface_mesh::Halfedge_index Halfedge_index;

    int nb_sites = sites.size();
    easy3d::Box3 box;

    for(auto p : sites) {box.grow(p);}
    Triangulation T;

    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T.insert(CGAL_Point(sites[i].x, sites[i].y, sites[i].z ));
        mp_sites[v] = i;
    }

    double diag = box.diagonal_length() * 10;
    double xmin = box.min_coord(0) - diag;
    double xmax = box.max_coord(0) + diag;
    double ymin = box.min_coord(1) - diag;
    double ymax = box.max_coord(1) + diag;
    double zmin = box.min_coord(2) - diag;
    double zmax = box.max_coord(2) + diag;
    T.insert(CGAL_Point(xmin, ymin, zmax));
    T.insert(CGAL_Point(xmax, ymin, zmax));
    T.insert(CGAL_Point(xmin, ymax, zmax));
    T.insert(CGAL_Point(xmax, ymax, zmax));
    T.insert(CGAL_Point(xmin, ymin, zmin));
    T.insert(CGAL_Point(xmax, ymin, zmin));
    T.insert(CGAL_Point(xmin, ymax, zmin));
    T.insert(CGAL_Point(xmax, ymax, zmin));

    vor_vertices.clear(); sites_nearby_vertices.clear(); vertices_nearby_sites.clear();
    sites_nearby_vertices.resize(nb_sites);
    for(auto cell = T.finite_cells_begin(); cell != T.finite_cells_end(); cell++) {
        auto v = T.dual(cell);
        int vid = vor_vertices.size();
        vor_vertices.emplace_back(easy3d::vec3(v.x(), v.y(), v.z()));
        std::vector<int> tmp;
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            sites_nearby_vertices[id].emplace_back(vid);
            tmp.emplace_back(id);
        }
        vertices_nearby_sites.emplace_back(tmp);
    }

}
void FilletSegV2::sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels) {
    int nb_points = points.size();
    if(nb_neighbors_ < 1 || std_ratio_ <= 0) {
        std::cout << "Illegal input parameters, the number of neighbors and "
                     "standard deviation ratio must be positive." << std::endl;
    }
    for(int i = 0; i < num_sor_iter_; i++) {
        std::vector<KNN::Point> knn_points;
        std::vector<int> indices;
#pragma omp parallel for
        for(int j = 0; j < nb_points; j++) {
            if(labels[j]) {
                #pragma omp critical
                {
                    knn_points.emplace_back(KNN::Point(points[j].x,
                                                       points[j].y, points[j].z));
                    indices.emplace_back(j);
                }
            }
        }
        KNN::KdSearch kds(knn_points);
        int num = knn_points.size();
        size_t valid_distances = 0;
        std::vector<double> avg_distances(num);
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(knn_points[j], nb_neighbors_, tmp_indices, dist);
            double mean = -1.0;

            if(dist.size() > 0u) {
                valid_distances++;
                std::for_each(dist.begin(), dist.end(),
                              [](double &d) { d = std::sqrt(d); });
                mean = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();
            }
            avg_distances[j] = mean;
        }
        if(valid_distances == 0) {
            continue;
        }
        double cloud_mean = std::accumulate(
        avg_distances.begin(), avg_distances.end(), 0.0,
        [](double const &x, double const &y) { return y > 0 ? x + y : x; });

        cloud_mean /= valid_distances;
        double sq_sum = std::inner_product(
        avg_distances.begin(), avg_distances.end(), avg_distances.begin(),
        0.0, [](double const &x, double const &y) { return x + y; },
        [cloud_mean](double const &x, double const &y) {
            return x > 0 ? (x - cloud_mean) * (y - cloud_mean) : 0;
        });
        double std_dev = std::sqrt(sq_sum / (valid_distances - 1));
        double distance_threshold = cloud_mean + std_ratio_ * std_dev;

#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            if(avg_distances[j] > 0 && avg_distances[j] < distance_threshold) {
                labels[indices[j]] = true;
            } else {
                labels[indices[j]] = false;
            }
        }
    }

}
