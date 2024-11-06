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
#include <easy3d/fileio/point_cloud_io.h>


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

void randperm(Eigen::VectorXi& dex, int a, int b) {
    int n = b - a;
    dex.resize(n);
    std::iota(dex.data(), dex.data() + n, a);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    std::shuffle(dex.data(),dex.data()+n, rng);
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

bool is_inside_triangle(easy3d::vec3 p, easy3d::vec3 a, easy3d::vec3 b, easy3d::vec3 c) {
    easy3d::vec3 n = easy3d::cross((b - a), (c - a)).normalize();
    easy3d::vec3 aix1 = (b - a).normalize();
    easy3d::vec3 aix2 = easy3d::cross(n, aix1).normalize();
    easy3d::vec2 a2d(0,0);
    easy3d::vec2 b2d( dot((b - a), aix1), dot((b - a), aix2));
    easy3d::vec2 c2d( dot((c - a), aix1), dot((c - a), aix2));
    easy3d::vec2 p2d( dot((p - a), aix1), dot((p - a), aix2));
    double x = easy3d::dot(a2d - p2d, b2d - p2d);
    double y = easy3d::dot(b2d - p2d, c2d - p2d);
    double z = easy3d::dot(c2d - p2d, a2d - p2d);
    if((x > 0 && y > 0 && z > 0) || (x < 0 && y < 0 && z < 0)) {
        return true;
    } else {
        return false;
    }
}
bool check(easy3d::vec3 p, easy3d::vec3 a, easy3d::vec3 b, easy3d::vec3 c, easy3d::vec3 d) {
    if(!is_inside_triangle(d, a,b,c)) {
        return is_inside_triangle(p, a, b, c) || is_inside_triangle(p, a,d,c);
    } else {
        return is_inside_triangle(p, a, b, c) || is_inside_triangle(p, a,d, b);
    }
}

FilletSegV2::FilletSegV2(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
    // fillet segmentation para
    eps_ = 0.01; sigma_ = 10; min_score_ = 0.5; angle_thr = 40; radius_ = 0.03;
    // sor para
    std_ratio_ = 0.5; nb_neighbors_ = 30; num_sor_iter_ = 4;
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

    dsa = 0;
}


void elkan_kmeans(Eigen::MatrixXd& data, int max_iter,
            Eigen::VectorXi& idx, Eigen::MatrixXd& center, int& m) {
    size_t num = data.rows(), dim = data.cols();
    Eigen::VectorXi dex;

    randperm(dex, 0, num);
    // for(int i = 0; i < dex.size(); i++) {
    //     std::cout << dex[i] << std::endl;
    // }
    // std::cout <<std::endl;
    center = data(dex.segment(0, m), Eigen::all);
    // return;
    Eigen::VectorXd tmp(m);
    idx.resize(num);
    idx.setConstant(-1);
    Eigen::VectorXi clust_size(m); clust_size.setZero();
    Eigen::MatrixXd center_sum(center.rows(), center.cols());
    center_sum.setZero();

    for(int iter = 0; iter < max_iter; iter++) {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < m; i++) {
            tmp[i] = std::numeric_limits<double>::max();
            for (int j = 0; j < m; j++) {
                if (i == j) continue;
                double dis = (center.row(i) - center.row(j)).squaredNorm();
                if (dis < tmp[i]) {
                    tmp[i] = dis;
                }
            }
        }

//        clust_size.setZero();
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < num; i++) {
            if(idx[i] == -1) {
                double minn = std::numeric_limits<double>::max();
                int id = 0;
                for(int j = 0; j < m; j++) {
                    double dis = (data.row(i) - center.row(j)).squaredNorm();
                    if(dis < minn) {
                        minn = dis; id = j;
                        if(2 * dis < tmp[j]) break;
                    }
                }
                idx[i] = id;
                clust_size[id]++;
                center_sum.row(id) += data.row(i);
            }
            else {
                double dis = (data.row(i) - center.row(idx[i])).squaredNorm();
                if(2 * dis <= tmp[idx[i]]) {
                    continue;
                } else {
                    double minn = std::numeric_limits<double>::max();
                    int id = 0;
                    for(int j = 0; j < m; j++) {
                        dis = (data.row(i) - center.row(j)).squaredNorm();
                        if(dis < minn) {
                            minn = dis; id = j;
                            if(2 * dis <= tmp[j]) break;
                        }
                    }
                    clust_size[idx[i]]--;
                    center_sum.row(idx[i]) -= data.row(i);

                    idx[i] = id;
                    clust_size[idx[i]]++;
                    center_sum.row(idx[i]) += data.row(i);
                }
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < m; i++) {
            // center.row(i) = center_sum.row(i) / clust_size[i];
            if(clust_size[i] == 0) {
                std::cout << "ASD"<<std::endl;
            }
            center.row(i) = data.row(idx[i]);
//            std::cout << center_sum.row(i) << ' ' << clust_size[i] << std::endl;
        }
    }
}

void sqdist_omp(Eigen::MatrixXd& a,Eigen::MatrixXd& b,
                   Eigen::MatrixXd& res) {
    int n = a.rows(), m = b.rows();
    //    res.resize(n, m);
#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res(i,j) = (a.row(i) - b.row(j)).squaredNorm();
        }
    }
}

void kmeans(Eigen::MatrixXd& data, int max_iter,
                Eigen::VectorXi& idx, Eigen::MatrixXd& center, int& m) {
    size_t num = data.rows(), dim = data.cols();
    Eigen::VectorXi dex;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    // igl::randperm(num, dex, rng);
    randperm(dex, 0, num);

    center = data(dex.segment(0, m), Eigen::all);
    Eigen::MatrixXd tmp(center.rows(), data.cols());


    for(int i = 0; i < max_iter; i++) {
        Eigen::VectorXd ct = Eigen::VectorXd::Zero(m);
        tmp.resize(center.rows(), data.rows());
        sqdist_omp(center, data, tmp);

        std::vector<double> min_values(tmp.cols(), std::numeric_limits<double>::max());
        idx.resize(tmp.cols());

#pragma omp parallel for
        for(int j = 0; j < tmp.cols(); j++) {
            for(int k = 0; k < tmp.rows(); k++) {
                if(min_values[j] > tmp(k, j)) {
                    min_values[j] = tmp(k, j);
                    idx[j] = k;
                }
            }
        }

        center.setZero();
#pragma omp parallel for
        for(int j = 0; j < num ; j++) {
#pragma omp critical
            {
                center.row(idx[j]) += data.row(j);
                ct[idx[j]]++;
            }
        }

#pragma omp parallel for
        for(int j = 0; j < m; j++) {
            center.row(j) = data.row(idx[j]);
        }
    }
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


void FilletSegV2::density(std::vector<easy3d::vec3>& points, std::vector<double>& avg_distances) {
    int nb_points = points.size();
    std::vector<KNN::Point> knn_points;
    std::vector<int> indices;
    avg_distances.resize(nb_points, 0);
#pragma omp parallel for
    for(int i = 0; i < nb_points; i++) {
#pragma omp critical
        {
            knn_points.emplace_back(KNN::Point(points[i].x,
                                               points[i].y, points[i].z));
            indices.emplace_back(i);
        }
    }
    KNN::KdSearch kds(knn_points);
    std::vector<double> avg_dist(nb_points);
    for(int i = 0; i < nb_points; i++) {
        std::vector<size_t> tmp_indices;
        std::vector<double> dist;
        kds.kth_search(knn_points[i], nb_neighbors_, tmp_indices, dist);
        double mean = -1.0;

        if(dist.size() > 0u) {
            std::for_each(dist.begin(), dist.end(),
                          [](double &d) { d = std::sqrt(d); });
            mean = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();
        }
        avg_distances[i] = mean;
        // std::cout << avg_dist[i] << std::endl;
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
        // std::cout << distance_threshold << std::endl;
#pragma omp parallel for
        for(int j = 0; j < num; j++) {
            // std::cout << avg_distances[j] << std::endl;
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

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        vor_cor_sites_.emplace_back(std::vector<int>());
        return 2.0;
    }

    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        err += abs((sites_[cur.idx()] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f) == vis.end() && len < R + tol) {
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
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f.idx());
        vor_cor_sites_.emplace_back(indices);
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        // easy3d::PointCloud* cloud = new easy3d::PointCloud;
        // for(size_t i = 0; i < indices.size(); i++) {
        //     cloud->add_vertex(sites_[indices[i]]);
        // }
        // cloud->add_vertex(vertex);
        // easy3d::PointCloudIO::save("../data/" + std::to_string(dsa++) + ".ply", cloud);
        return 2.0;
    }


    Eigen::MatrixXd S(num, 3);
    for(int i = 0; i < num; i++) {
        S.row(i) = Eigen::Vector3d(sites_[i].x - vertex.x, sites_[i].y - vertex.y, sites_[i].z - vertex.z);
    }
    Eigen::MatrixXd M = S.transpose() * S;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigen solver failed!" << std::endl;
        return 1e9;
    }

    if(eigensolver.eigenvalues()[0] > 1e-6) {
        return exp(-0.1 *  eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[0]) + (1 -  exp(- 0.1 * (eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[1] - 1.0)));
    } else {
        return 2.0;
    }
    // return eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[0];
}


float FilletSegV2::probability_of_vertex_eigen_value(easy3d::vec3& vertex, std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites) {

    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    for(auto f : vertex_nearby_sites) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites[0].idx()] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        return 1.0;
    }

    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        err += abs((sites_[cur.idx()] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f) == vis.end() && len < R + tol) {
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
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f.idx());
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        return 1.0;
    }


    Eigen::MatrixXd S(num, 3);
    for(int i = 0; i < num; i++) {
        S.row(i) = Eigen::Vector3d(sites_[i].x - vertex.x, sites_[i].y - vertex.y, sites_[i].z - vertex.z);
    }
    Eigen::MatrixXd M = S.transpose() * S;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigen solver failed!" << std::endl;
        return 1.0;
    }

    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(size_t i = 0; i < indices.size(); i++) {
        cloud->add_vertex(sites_[indices[i]]);
    }
    cloud->add_vertex(vertex);
    if((eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[0]) > 10.0) {
        // easy3d::PointCloudIO::save("../data2/" + std::to_string(dsa++) + ".ply", cloud);;
    std::cout << eigensolver.eigenvalues()[0] << ' ' << eigensolver.eigenvalues()[1] << ' ' << eigensolver.eigenvalues()[2] << std::endl;
    }
    if(eigensolver.eigenvalues()[0] > 1e-6) {
        return exp(-0.1 *  eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[0]);
    } else {
        return 1.0;
    }
    // return eigensolver.eigenvalues()[2] / eigensolver.eigenvalues()[0];
}

float FilletSegV2::probability_of_vertex_area(easy3d::vec3& vertex, std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites) {


    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    for(auto f : vertex_nearby_sites) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites[0].idx()] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        return 0.0;
    }

    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        err += abs((sites_[cur.idx()] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f) == vis.end() && len < R + tol) {
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
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f.idx());
        err += face_area[f];
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        // easy3d::PointCloud* cloud = new easy3d::PointCloud;
        // for(size_t i = 0; i < indices.size(); i++) {
        //     cloud->add_vertex(sites_[indices[i]]);
        // }
        // cloud->add_vertex(vertex);
        // easy3d::PointCloudIO::save("../data/" + std::to_string(dsa++) + ".ply", cloud);
        return 0.0;
    }
    return err / (4 * R * R);
}

float FilletSegV2::probability_of_vertex_normal(easy3d::vec3& vertex, std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites) {


    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    for(auto f : vertex_nearby_sites) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites[0].idx()] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        return 1.0;
    }

    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        err += abs((sites_[cur.idx()] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f) == vis.end() && len < R + tol) {
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
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f.idx());
        easy3d::vec3 asd = (sites_[f.idx()] - vertex).normalize();
        err += easy3d::cross(asd, face_normals[f]).length2();
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        // easy3d::PointCloud* cloud = new easy3d::PointCloud;
        // for(size_t i = 0; i < indices.size(); i++) {
        //     cloud->add_vertex(sites_[indices[i]]);
        // }
        // cloud->add_vertex(vertex);
        // easy3d::PointCloudIO::save("../data/" + std::to_string(dsa++) + ".ply", cloud);
        return 1.0;
    }
    return err / vis.size();
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

void FilletSegV2::eigen_value(bool sor_flag) {
    std::vector<easy3d::vec3> tmp;
    int idx = 0;
    int ct = 0;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 100 == 0) {
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
                    ct++;
                    std::vector<int> indices(4);
                    for(int j = 0; j < 4; j++) {
                        indices[j] = patch[vertices_nearby_sites[i][j]].idx();
                    }
                    std::sort(indices.begin(), indices.end());
                    bool flag = true;
                    if(mp_.find(indices[0]) != mp_.end()) {
                        auto& sub_mp1 = mp_[indices[0]];
                        if(sub_mp1.find(indices[1]) != sub_mp1.end()) {
                            auto& sub_mp2 = sub_mp1[indices[1]];
                            if(sub_mp2.find(indices[2]) != sub_mp2.end()) {
                                auto& sub_mp3 = sub_mp2[indices[2]];
                                if(sub_mp3.find(indices[3]) != sub_mp3.end()) {
                                    flag = false;
                                }
                            }
                        }
                    }
                    if(flag) {
                        total_vor_vertices_.emplace_back(vor_vertices[i]);
                        total_vor_nearby_sites_.emplace_back(indices);
                        mp_[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                    }
                }
            }
        }
    }
    std::cout << idx << ' ' << ct << std::endl;
    std::vector<bool> labels(total_vor_vertices_.size(), true);
    if(sor_flag)
        sor(total_vor_vertices_, labels);
    for(size_t i = 0; i < labels.size(); i++) {
        if(labels[i]) {
            std::vector<easy3d::SurfaceMesh::Face> ss = {
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][0]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][1]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][2]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][3])
            };

            double score = probability_of_vertex(total_vor_vertices_[i], ss);
            vor_vertices_.emplace_back(total_vor_vertices_[i]);
            vor_vertices_score_.emplace_back(score);
        }
    }
    auto maxElement = std::max_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    auto minElement = std::min_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    std::cout << *maxElement << ' ' << *minElement << std::endl;
}

void FilletSegV2::probability_of_sites() {
    // std::vector<std::vector<int>> sites_nearby_vertices;
    // std::vector<std::vector<int>> vertices_nearby_sites;
    std::vector<easy3d::vec3> tmp;
    int idx = 0;
    int ct = 0;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 100 == 0) {
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
                    ct++;
                    std::vector<int> indices(4);
                    for(int j = 0; j < 4; j++) {
                        indices[j] = patch[vertices_nearby_sites[i][j]].idx();
                    }
                    std::sort(indices.begin(), indices.end());
                    bool flag = true;
                    if(mp_.find(indices[0]) != mp_.end()) {
                        auto& sub_mp1 = mp_[indices[0]];
                        if(sub_mp1.find(indices[1]) != sub_mp1.end()) {
                            auto& sub_mp2 = sub_mp1[indices[1]];
                            if(sub_mp2.find(indices[2]) != sub_mp2.end()) {
                                auto& sub_mp3 = sub_mp2[indices[2]];
                                if(sub_mp3.find(indices[3]) != sub_mp3.end()) {
                                    flag = false;
                                }
                            }
                        }
                    }
                    if(flag) {
                        total_vor_vertices_.emplace_back(vor_vertices[i]);
                        total_vor_nearby_sites_.emplace_back(indices);
                        mp_[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                    }
                }
            }
        }
    }
    std::cout << idx << ' ' << ct << std::endl;
    std::vector<bool> labels(total_vor_vertices_.size(), true);
    // sor(total_vor_vertices_, labels);
    for(size_t i = 0; i < labels.size(); i++) {
        if(labels[i]) {
            std::vector<easy3d::SurfaceMesh::Face> ss = {
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][0]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][1]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][2]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][3])
            };

            double score = probability_of_vertex(total_vor_vertices_[i], ss);
            // if(score < 5) {
                vor_vertices_.emplace_back(total_vor_vertices_[i]);
                vor_vertices_score_.emplace_back(score);
            // }
        }

    }
    // auto maxElement = std::max_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    // auto minElement = std::min_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    // std::cout << *maxElement << ' ' << *minElement << std::endl;
   //  Eigen::MatrixXd data(labels.size(), 3);
   //  Eigen::VectorXi indices;
   // Eigen::MatrixXd center;
   //  int id = 0;
   //  for(size_t i = 0; i < labels.size(); i++) {
   //      if(labels[i]) {
   //          data.row(id++) = Eigen::Vector3d(total_vor_vertices_[i].x, total_vor_vertices_[i].y, total_vor_vertices_[i].z);
   //      }
   //  }
   //  data.conservativeResize(id, 3);
   //  int m = 1000;
   //  elkan_kmeans(data, 100, indices, center, m);
   //  // kmeans(data, 100, indices, center, m);
   //  for(size_t i = 0; i < center.rows(); i++) {
   //      easy3d::vec3 p(center(i,0), center(i,1), center(i,2));
   //      std::cout << p.x << ' ' << p.y << ' ' << p.z << std::endl;
   //      vor_vertices_.emplace_back(p);
   //      vor_vertices_score_.emplace_back(0.0);
   //  }
}

void FilletSegV2::probability_of_sites2() {
    // std::vector<std::vector<int>> sites_nearby_vertices;
    // std::vector<std::vector<int>> vertices_nearby_sites;
    std::vector<easy3d::vec3> tmp;
    int idx = 0;
    int ct = 0;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 100 == 0) {
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
                    ct++;
                    std::vector<int> indices(4);
                    for(int j = 0; j < 4; j++) {
                        indices[j] = patch[vertices_nearby_sites[i][j]].idx();
                    }
                    std::sort(indices.begin(), indices.end());
                    bool flag = true;
                    if(mp_.find(indices[0]) != mp_.end()) {
                        auto& sub_mp1 = mp_[indices[0]];
                        if(sub_mp1.find(indices[1]) != sub_mp1.end()) {
                            auto& sub_mp2 = sub_mp1[indices[1]];
                            if(sub_mp2.find(indices[2]) != sub_mp2.end()) {
                                auto& sub_mp3 = sub_mp2[indices[2]];
                                if(sub_mp3.find(indices[3]) != sub_mp3.end()) {
                                    flag = false;
                                }
                            }
                        }
                    }
                    if(flag) {
                        total_vor_vertices_.emplace_back(vor_vertices[i]);
                        total_vor_nearby_sites_.emplace_back(indices);
                        mp_[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                    }
                }
            }
        }
    }
    std::cout << idx << ' ' << ct << std::endl;
    std::vector<bool> labels(total_vor_vertices_.size(), true);
    // sor(total_vor_vertices_, labels);
    for(size_t i = 0; i < labels.size(); i++) {
        if(labels[i]) {
            std::vector<easy3d::SurfaceMesh::Face> ss = {
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][0]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][1]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][2]),
                easy3d::SurfaceMesh::Face(total_vor_nearby_sites_[i][3])
            };

            double angle = transition_angle(total_vor_vertices_[i], ss);
            if(angle > 0) {
                vor_vertices_.emplace_back(total_vor_vertices_[i]);
                vor_transition_angle_.emplace_back(angle);
                double r = (sites_[total_vor_nearby_sites_[i][0]] - total_vor_vertices_[i]).norm();
                vor_radius_.emplace_back(r);
            }
        }
    }

    std::vector<KNN::Point> knn_points;
#pragma omp parallel for
    for(int i = 0; i < vor_vertices_.size(); i++) {
#pragma omp critical
            {
                knn_points.emplace_back(KNN::Point(vor_vertices_[i].x,
                                                   vor_vertices_[i].y, vor_vertices_[i].z));
            }
    }

    KNN::KdSearch kds(knn_points);
    int k = 30;
    for(int i = 0; i < vor_vertices_.size(); i++) {
        std::vector<size_t> tmp_indices;
        std::vector<double> dist;
        kds.kth_search(knn_points[i], k, tmp_indices, dist);
        double tot_r = 0;
        double tot_a = 0;
        double tot_d = 0;
        for(int j = 0; j < tmp_indices.size(); j++) {
            double r_dis = std::abs(vor_radius_[tmp_indices[j]] - vor_radius_[i]);
            double a_dis = std::abs(vor_transition_angle_[tmp_indices[j]] - vor_transition_angle_[i]);
            double d_dis = (vor_vertices_[tmp_indices[j]] - vor_vertices_[i]).norm();
            tot_r += r_dis;
            tot_a += a_dis;
            tot_d += d_dis;
        }
        tot_r = tot_r / k;
        tot_a = tot_a / k;
        tot_d = tot_d / k;
        // if(sum / 10 > 1.0) {
        //     vor_vertices_score_.emplace_back(1.0);
        // } else
        vor_vertices_score_.emplace_back( tot_a);
    }
    auto max_ele = std::max_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    auto min_ele = std::min_element(vor_vertices_score_.begin(), vor_vertices_score_.end());
    // for(size_t i = 0; i < vor_vertices_score_.size(); i++) {
    //
    // }
    std::cout << *max_ele << ' ' << *min_ele << std::endl;
    // std::vector<double> fff = vor_vertices_score_;
    // std::sort(fff.begin(), fff.end());


}

bool FilletSegV2::vor_vertex_aixs(easy3d::vec3& vertex, std::vector<int>& vertex_nearby_sites_indices, easy3d::vec3& axis) {
    std::set<int> vis;
    std::queue<int>que;
    for(auto f : vertex_nearby_sites_indices) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites_indices[0]] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        return 1.0;
    }

    while(!que.empty()) {
        auto cur_id = que.front();
        easy3d::SurfaceMesh::Face cur(cur_id);
        que.pop();
        err += abs((sites_[cur_id] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f.idx()) == vis.end() && len < R + tol) {
                        que.push(opp_f.idx());
                        vis.insert(opp_f.idx());
                }
            }
        }
    }

    std::set<int> vis2;
    que.push(vertex_nearby_sites_indices[0]);
    while(!que.empty()) {
        auto cur_id = que.front();
        easy3d::SurfaceMesh::Face cur(cur_id);
        que.pop();
        if (vis2.find(cur.idx()) != vis2.end()) {
            continue;
        }
        vis2.insert(cur.idx());
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && vis.find(opp_f.idx()) != vis.end() && vis2.find(opp_f.idx()) == vis2.end()) {
                que.push(opp_f.idx());
                vis2.insert(cur.idx());
            }
        }
    }
    // // return vis2.size();
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f);
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        return false;
    }


    Eigen::MatrixXd S(num, 3);
    for(int i = 0; i < num; i++) {
        S.row(i) = Eigen::Vector3d(sites_[indices[i]].x - vertex.x, sites_[indices[i]].y - vertex.y, sites_[indices[i]].z - vertex.z);
    }
    // Eigen::MatrixXd S(num, 3);
    Eigen::MatrixXd M = S.transpose() * S;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigen solver failed!" << std::endl;
        return false;
    }
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
    axis = easy3d::vec3(eigenvectors.col(0)[0], eigenvectors.col(0)[1], eigenvectors.col(0)[2]).normalize();
    return true;
}

void FilletSegV2::sites_probablity(bool sor_flag) {
    int idx = 0, ct = 0;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 100 == 0) {
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
                    ct++;
                    std::vector<int> indices(4);
                    for(int j = 0; j < 4; j++) {
                        indices[j] = patch[vertices_nearby_sites[i][j]].idx();
                    }
                    std::sort(indices.begin(), indices.end());
                    bool flag = true;
                    if(mp_.find(indices[0]) != mp_.end()) {
                        auto& sub_mp1 = mp_[indices[0]];
                        if(sub_mp1.find(indices[1]) != sub_mp1.end()) {
                            auto& sub_mp2 = sub_mp1[indices[1]];
                            if(sub_mp2.find(indices[2]) != sub_mp2.end()) {
                                auto& sub_mp3 = sub_mp2[indices[2]];
                                if(sub_mp3.find(indices[3]) != sub_mp3.end()) {
                                    flag = false;
                                }
                            }
                        }
                    }
                    if(flag) {
                        total_vor_vertices_.emplace_back(vor_vertices[i]);
                        total_vor_nearby_sites_.emplace_back(indices);
                        mp_[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                    }
                }
            }
        }
    }
    std::cout << idx << ' ' << ct << std::endl;
    std::vector<bool> labels(total_vor_vertices_.size(), true);

    if(sor_flag)
        sor(total_vor_vertices_, labels);

    std::vector<KNN::Point> knn_points;
#pragma omp parallel for
    for(int i = 0; i < vor_vertices_.size(); i++) {
#pragma omp critical
        {
            // if(labels[i])
                knn_points.emplace_back(KNN::Point(total_vor_vertices_[i].x,
                                               total_vor_vertices_[i].y, total_vor_vertices_[i].z));
        }
    }

    KNN::KdSearch kds(knn_points);
    int k = 100;
    vor_vertices_ = total_vor_vertices_;
    vor_vertices_score_.resize(vor_vertices_.size());
    for(size_t i = 0; i < labels.size(); i++) {
        if(labels[i]) {
            std::vector<int> sites_indices = {total_vor_nearby_sites_[i][0],  total_vor_nearby_sites_[i][1],
            total_vor_nearby_sites_[i][2], total_vor_nearby_sites_[i][3]};
            easy3d::vec3 aixs1;
            // std::cout << "ASD" << std::endl;
            if(!vor_vertex_aixs(total_vor_vertices_[i], sites_indices, aixs1)) {
                vor_vertices_score_[i] = 1.0;
                continue;
            }

            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(knn_points[i], k, tmp_indices, dist);
            Eigen::MatrixXd S(k, 3);
            int cc = 0;
            for(int j = 0; j < k; j++) {
                // if(tmp_indices[j] != i)
                S.row(cc++) = Eigen::Vector3d(total_vor_vertices_[tmp_indices[j]].x - total_vor_vertices_[i].x
                                           ,total_vor_vertices_[tmp_indices[j]].y - total_vor_vertices_[i].y
                                           , total_vor_vertices_[tmp_indices[j]].z - total_vor_vertices_[i].z);
            }
            // std::cout << "ASD" << std::endl;
            Eigen::MatrixXd M = S.transpose() * S;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(M);
            if (eigensolver.info() != Eigen::Success) {
                std::cerr << "Eigen solver failed!" << std::endl;
                vor_vertices_score_[i] = 1.0;
            }
            Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
            Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
            easy3d::vec3 aixs2 = easy3d::vec3(eigenvectors.col(2)[0], eigenvectors.col(2)[1], eigenvectors.col(2)[2]).normalize();
            double vv = easy3d::cross(aixs1, aixs2).norm();
            vor_vertices_score_[i] = vv;
            // std::cout << vv <<std::endl;
        }else {
            vor_vertices_score_[i] = 1.0;
        }
        // std::cout << i << std::endl;
    }
}


double FilletSegV2::transition_angle(easy3d::vec3& vertex,
        std::vector<easy3d::SurfaceMesh::Face>& vertex_nearby_sites) {
        std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    for(auto f : vertex_nearby_sites) {
        vis.insert(f);
        que.push(f);
    }

    double R = (sites_[vertex_nearby_sites[0].idx()] - vertex).norm();
    double tol = R * eps_;
    double err = 0.0;

    double thr = box_.diagonal_length() * 0.07;
    if(R > thr) {
        return -1.0;
    }

    while(!que.empty()) {
        auto cur = que.front();
        que.pop();
        err += abs((sites_[cur.idx()] - vertex).norm() - R);
        for (auto h: mesh_->halfedges(cur)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur, opp_f) < 30) {
                double len = (sites_[opp_f.idx()] - vertex).norm();
                if (vis.find(opp_f) == vis.end() && len < R + tol) {
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
    std::vector<int> indices;
    for(auto f : vis) {
        indices.emplace_back(f.idx());
    }
    int num = indices.size();
    if(vis2.size() != vis.size()) {
        return -1.0;
    }
    double maxx = 0;
    for(int i = 0;  i < num; i++) {
        easy3d::vec3 n1 = (sites_[indices[i]] - vertex).normalize();
        // easy3d::vec3 n1 = face_normals[easy3d::SurfaceMesh::Face(indices[i])];
        for(int j = i + 1; j < num; j++) {
            easy3d::vec3 n2 = (sites_[indices[j]] - vertex).normalize();
            // easy3d::vec3 n2 = face_normals[easy3d::SurfaceMesh::Face(indices[j])];
            double radians = abs(acos(dot(n1, n2)));
            double degrees = radians * 180.0 / M_PI;
            if(degrees > maxx) {
                maxx = degrees;
            }
        }
    }
    // if(maxx < 30) {
    //     return -1.0;
    // }
    return maxx;
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
            std::vector<double> squared_distances;
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
