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

#include <easy3d/algo/surface_mesh_geometry.h>
#include <easy3d/core/random.h>


#include <Eigen/Dense>

double dihedral_angle(easy3d::SurfaceMesh* mesh
                      , easy3d::SurfaceMesh::Face f1
                      , easy3d::SurfaceMesh::Face f2
                      , bool rad = false) {
    easy3d::vec3 n1 = mesh->compute_face_normal(f1).normalize();
    easy3d::vec3 n2 = mesh->compute_face_normal(f2).normalize();

    double radians = abs(acos(dot(n1, n2)));

    if(rad) {
        return radians;
    } else {
        double degrees = radians * 180.0 / M_PI;
        return degrees;
    }
}

bool point_inside_tetrahedron(easy3d::vec3 p, std::vector<easy3d::vec3> tetra) {
    Eigen::Vector3d P(p.x, p.y, p.z);
    Eigen::Vector3d A(tetra[0].x, tetra[0].y, tetra[0].z);
    Eigen::Vector3d B(tetra[1].x, tetra[1].y, tetra[1].z);
    Eigen::Vector3d C(tetra[2].x, tetra[2].y, tetra[2].z);
    Eigen::Vector3d D(tetra[3].x, tetra[3].y, tetra[3].z);
    auto signedVolume = [](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4) -> double {
        Eigen::Matrix3d M;
        M.col(0) = p2 - p1;
        M.col(1) = p3 - p1;
        M.col(2) = p4 - p1;
        return M.determinant() / 6.0; // 带符号的体积
    };
    double volumeABCD = signedVolume(A, B, C, D);
    if (std::abs(volumeABCD) < 1e-6) {
        // std::cerr << "Degenerate tetrahedron." << std::endl;
        return false;
    }

    // 计算点P与四面体每个面形成的四个带符号的体积
    double volumePBCD = signedVolume(P, B, C, D);
    double volumePACD = signedVolume(P, A, C, D);
    double volumePABD = signedVolume(P, A, B, D);
    double volumePABC = signedVolume(P, A, B, C);

    // 检查所有体积的符号是否一致
    bool sameSign = (volumeABCD > 0 && volumePBCD > 0 && volumePACD > 0 && volumePABD > 0 && volumePABC > 0) ||
                    (volumeABCD < 0 && volumePBCD < 0 && volumePACD < 0 && volumePABD < 0 && volumePABC < 0);

    return sameSign;
}

FilletSegV2::FilletSegV2(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
    // fillet segmentation para
    eps_ = 0.005; sigma_ = 10; min_score_ = 0.5; angle_thr = 30; radius_ = 0.03;
    // sor para
    std_ratio_ = 1.0; nb_neighbors_ = 30; num_sor_iter_ = 3;
    face_normals = mesh_->face_property<easy3d::vec3>("f:normal");
    face_area = mesh_->face_property<float>("f:area");
    for(auto f : mesh_->faces()) {
        face_area[f] = easy3d::geom::triangle_area(mesh_, f);
        face_normals[f] = mesh_->compute_face_normal(f);
        int ct = 0;
        easy3d::vec3 p(0,0,0);
        for(auto v : mesh_->vertices(f)) {
            p += mesh_->position(v); ct++;
        }
        if(ct)
            sites_.emplace_back(p / ct);
    }
    face_score = mesh_->face_property<float>("f:scores");
    box_ = mesh_->bounding_box();
}


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


void FilletSegV2::density(std::vector<easy3d::vec3>& points, std::vector<float>& avg_distances) {
    int num = points.size();
    std::vector<KNN::Point> knn_points(num);
    std::vector<int> indices;
    for(int i = 0; i < num; i++) {
        knn_points.emplace_back(KNN::Point(points[i].x,
                                                       points[i].y, points[i].z));
        indices.emplace_back(i);
    }
    KNN::KdSearch kds(knn_points);
    avg_distances.resize(num);
    // std::vector<double> avg_distances(num);
#pragma omp parallel for
    for(int j = 0; j < num; j++) {
        std::vector<size_t> tmp_indices;
        std::vector<float> dist;
        kds.kth_search(knn_points[j], nb_neighbors_, tmp_indices, dist);
        float mean = -1.0;

        if(dist.size() > 0u) {
            std::for_each(dist.begin(), dist.end(),
                          [](float &d) { d = std::sqrt(d); });
            mean = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();
        }
        avg_distances[j] = mean;
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
            std::vector<float> dist;
            kds.kth_search(knn_points[j], nb_neighbors_, tmp_indices, dist);
            double mean = -1.0;

            if(dist.size() > 0u) {
                valid_distances++;
                std::for_each(dist.begin(), dist.end(),
                              [](float &d) { d = std::sqrt(d); });
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


float FilletSegV2::probability_of_vertex(easy3d::vec3& vertex, std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites) {


    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    for(auto f : vertex_nearby_sites) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites[0].idx()] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;


    // easy3d::vec3 ab = sites_[vertex_nearby_sites[0].idx()] - sites_[vertex_nearby_sites[1].idx()];
    // easy3d::vec3 ac = sites_[vertex_nearby_sites[0].idx()] - sites_[vertex_nearby_sites[2].idx()];
    // easy3d::vec3 ad = sites_[vertex_nearby_sites[0].idx()] - sites_[vertex_nearby_sites[3].idx()];
    // return 1.0 / 6.0 * fabs(easy3d::dot(ab, easy3d::cross(ac, ad))) / (4.0 / 3.0 * M_PI * R * R * R);

    double tmp = radius_ * box_.diagonal_length();
    // if(R > tmp) {
    //     return 0;
    // }
    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        // if (vis.find(cur) != vis.end()) {
        //     continue;
        // }
        err += fabs((sites_[cur.idx()] - vertex).norm() - R);
        // err += face_area[cur];
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid()) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (len < R + tol && dihedral_angle(mesh_, cur, opp_f) < angle_thr
                    && vis.find(opp_f) == vis.end()) {
                        que.push(opp_f);
                        vis.insert(opp_f);
                }
            }
        }
    }

    std::set<easy3d::SurfaceMesh::Face> vis2;
    que.push(vertex_nearby_sites[0]);
    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        if (vis2.find(cur) != vis2.end()) {
            continue;
        }
        vis2.insert(cur);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && vis.find(opp_f) != vis.end() && vis2.find(opp_f) == vis2.end()) {
                que.push(opp_f);
                vis2.insert(cur);
            }
        }
    }
    // // return vis2.size();
    if(vis2.size() != vis.size()) {
        return 0.0;
    }
    // std::cout << vis.size() << std::endl;
    return exp(-sigma_ * err / vis.size());
    // return err / (4 * M_PI * R * R);
}


void FilletSegV2::crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
                        , std::vector<easy3d::SurfaceMesh::Face>& patch) {
    radius = box_.diagonal_length() * radius;
    double area_thr = 2 * M_PI * radius * radius;
    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face> que;
    que.push(f); vis.insert(f);
    double area_sum = 0;
    while(area_sum < area_thr && !que.empty()) {
        auto cur_f = que.front(); que.pop();
        area_sum += face_area[cur_f];
        for(auto h : mesh_->halfedges(cur_f)) {
            auto nxt_f = mesh_->face(mesh_->opposite(h));
            if(nxt_f.is_valid() && vis.find(nxt_f) == vis.end()
               && dihedral_angle(mesh_, cur_f, nxt_f) < angle_thr) {
                vis.insert(nxt_f);
                que.push(nxt_f);
               }
        }
    }

    patch.insert(patch.end(), vis.begin(), vis.end());
}

void FilletSegV2::probability_of_sites() {
    // std::vector<std::vector<int>> sites_nearby_vertices;
    // std::vector<std::vector<int>> vertices_nearby_sites;
    std::vector<easy3d::vec3> tmp;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 50 == 0) {
            std::vector<easy3d::SurfaceMesh::Face> patch;
            crop_local_patch(f, radius_, patch);
            std::vector<easy3d::vec3> sites;
            for(size_t i = 0; i < patch.size(); i++) {
                sites.emplace_back(sites_[patch[i].idx()]);
            }
            std::vector<std::vector<int>> sites_nearby_vertices;
            std::vector<std::vector<int>> vertices_nearby_sites;
            std::vector<easy3d::vec3> vor_vertices;
            voronoi3d(sites,  vor_vertices,sites_nearby_vertices, vertices_nearby_sites);
            for(size_t i = 0; i < vor_vertices.size(); i++) {
                if(box_.contains(vor_vertices[i])) {
                    tmp.emplace_back(vor_vertices[i]);
                }
            }
        }
    }

    std::vector<bool> label(tmp.size(), true);
    deduplication(tmp, label);
    std::vector<easy3d::vec3> tmp2;
    for(size_t i = 0; i < tmp.size(); i++) {
        if(label[i]) {
            tmp2.emplace_back(tmp[i]);
        }
    }

    std::vector<bool> label2(tmp2.size(), true);
    sor(tmp2, label2);
    for(size_t i = 0; i < tmp2.size(); i++) {
        if(label2[i]) {
            vor_vertices_.emplace_back(tmp2[i]);
        }
    }
    std::cout << "ASD" <<std::endl;
    // std::cout << "ASD" <<std::endl;
    density(vor_vertices_, vor_vertices_score_);
}

void FilletSegV2::deduplication(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels) {
    int nb_points = points.size();
    std::vector<KNN::Point> knn_points;
    std::vector<int> indices;

    for(int j = 0; j < nb_points; j++) {
        knn_points.emplace_back(KNN::Point(points[j].x,
                                           points[j].y, points[j].z));
        indices.emplace_back(j);
    }
    KNN::KdSearch kds(knn_points);
    for(int i = 0; i < nb_points; i++) {
        if(labels[i]) {
            std::vector<size_t> neighbors;
            std::vector<float> squared_distances;
            int num = kds.radius_search(knn_points[i], 1e-3, neighbors, squared_distances);
            for(int j = 0; j < num; j++) {
                // double d = sqrt(squared_distances[j]);
                if(neighbors[j] != i && squared_distances[j] < 1e-6) {
                    labels[neighbors[j]] = false;
                }
            }
        }
    }
}

void FilletSegV2::decompose() {
    std::vector<std::vector<int>> sites_nearby_vertices;
    std::vector<std::vector<int>> vertices_nearby_sites;
    std::vector<easy3d::vec3> vor_vertices;
    voronoi3d(sites_,  vor_vertices, sites_nearby_vertices, vertices_nearby_sites);

    std::vector<bool> labels(vor_vertices.size(), true);
    for(size_t i = 0; i < vor_vertices.size(); i++) {
        if(!box_.contains(vor_vertices[i])) {
            labels[i] = false;
        }
    }
    sor(vor_vertices, labels);
    double thr = box_.diagonal_length();
    site_cconfilt_.resize(sites_.size());
    for(size_t i = 0; i < vor_vertices.size(); i++) {
        double R = (sites_[vertices_nearby_sites[i][0]] - vor_vertices[i]).norm();
        if(labels[i] && R < thr) {
            for(int j = 0; j < 4; j++) {
                easy3d::vec3 n1 = face_normals[easy3d::SurfaceMesh::Face(vertices_nearby_sites[i][j])];
                for(int k = j + 1; k < 4; k++) {
                    easy3d::vec3 n2 = face_normals[easy3d::SurfaceMesh::Face(vertices_nearby_sites[i][k])];
                    double angle = acos(easy3d::dot(n1, n2)) * 180.0 / M_PI;
                    if(angle > 90) {
                        site_cconfilt_[vertices_nearby_sites[i][j]].emplace_back(vertices_nearby_sites[i][k]);
                        site_cconfilt_[vertices_nearby_sites[i][k]].emplace_back(vertices_nearby_sites[i][j]);
                    }
                }
            }

        }
    }
    auto face_label = mesh_->face_property<int>("f:labels");
    for(auto f : mesh_->faces()) {
        face_label[f] = -1;
    }
    int ct = 0;
    for(auto f : mesh_->faces()) {
        if(face_label[f] == -1) {
            std::queue<easy3d::SurfaceMesh::Face> que;
            que.push(f); face_label[f] = ct;
            std::set<int> confilct;
            confilct.insert(site_cconfilt_[f.idx()].begin(), site_cconfilt_[f.idx()].end());
            while(!que.empty()) {
                auto cur_f = que.front(); que.pop();
                for(auto h : mesh_->halfedges(cur_f)) {
                    auto opp_f = mesh_->face(mesh_->opposite(h));
                    if(opp_f.is_valid() && face_label[opp_f] == -1 && confilct.find(opp_f.idx()) == confilct.end()
                        && dihedral_angle(mesh_, cur_f, opp_f) < 30) {
                        face_label[opp_f] = ct;
                        confilct.insert(site_cconfilt_[opp_f.idx()].begin(), site_cconfilt_[opp_f.idx()].end());
                        que.push(opp_f);
                    }
                }
            }
            ct++;
        }
    }
    auto face_color = mesh_->face_property<easy3d::vec3>("f:color");
    std::vector<easy3d::vec3> color(ct);
    for(int i = 0; i < ct; i++) {
        color[i] = easy3d::random_color();
    }

    for(auto f : mesh_->faces()) {
        face_color[f] = color[face_label[f]];
    }
}
