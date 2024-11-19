//
// Created by xiaowuga on 2024/11/14.
//

#include "fillet_seg_v5.h"

#include <knn.h>
#include <mean_shift.h>
#include <random>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>

#include <easy3d/algo/surface_mesh_geometry.h>
#include <easy3d/fileio/point_cloud_io.h>

FilletSegV5::FilletSegV5(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
    // fillet segmentation para
    eps_ = 0.03; sigma_ = 10; angle_thr_ = 40;

    // sor segmentation parameters
    nb_neighbors_ = 30;
    num_sor_iter_ = 3;
    std_radio_ = 0.3;

    sa_.resize(mesh_->n_faces());
    for(auto f : mesh_->faces()) {
        sa_[f.idx()] = easy3d::geom::triangle_area(mesh_, f);
        int ct = 0;
        easy3d::vec3 p(0,0,0);
        for(auto v : mesh_->vertices(f)) {
            p += mesh_->position(v); ct++;
        }
        if(ct)
            s_.emplace_back(p / ct);
    }
    box_ = mesh_->bounding_box();
    radius_ = box_.diagonal_length() * 0.08;

}

void FilletSegV5::run() {
    std::cout << "Start generate voronoi vertices..." <<std::endl;
    generate_local_voronoi_vertices(vv_, vvr_, vvns_, false);
    std::cout << "Generate voronoi vertices sussessfully!" <<std::endl;
    // std::string vv_path = "../vv.ply";
    voronoi_density_drift(vv_, vvr_, vvns_, scvv_);
    // save_point_set(vv_, vv_path);


}

void FilletSegV5::generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv, std::vector<float>& vvr
                                                , std::vector<std::vector<int>>& vvns
                                                ,bool sor_flag, bool radius_thr) {
    int idx = 0;
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> mp_;
    std::vector<easy3d::vec3> tmp_vv;
    std::vector<float> tmp_vvr;
    std::vector<std::vector<int>> tmp_vvns;

    int num_samples = 5000;
    std::vector<int> ind;
    float max_gap = 0.5 * M_PI * radius_;
    farthest_point_sampling(num_samples, ind, max_gap);
    for(int k = 0; k < num_samples; k++) {
        std::vector<int> patch;
        // std::cout << k <<' ' << num_samples<< std::endl;
        crop_local_patch(ind[k], max_gap, patch);

        std::vector<easy3d::vec3> ls; // local sites;
        for(size_t i = 0; i < patch.size(); i++) {
            ls.emplace_back(s_[patch[i]]);
        }

        std::vector<std::vector<int>> lsnv;   // local sites neighboring vertices;
        std::vector<float> lvvr; // vor_vertices_radius
        std::vector<std::vector<int>> lvns ; // local vertices neighboring sites;
        std::vector<easy3d::vec3> lvv; // local voronoi vertices;
        voronoi3d(ls, lvv, lvvr,lsnv, lvns);

        for(size_t i = 0; i < lvv.size(); i++) {
            std::vector<int> indices(4);
            for(int j = 0; j < 4; j++) {
                indices[j] = patch[lvns[i][j]];
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
            if(flag && box_.contains(lvv[i])) {
                tmp_vv.emplace_back(lvv[i]);
                tmp_vvr.emplace_back(lvvr[i]);
                tmp_vvns.emplace_back(indices);
                mp_[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
            }
        }
    }

    std::vector<bool> labels(tmp_vv.size(), true);
    if(radius_thr) {
        for(size_t i = 0; i < tmp_vvr.size(); i++) {
            if(tmp_vvr[i] > radius_) {
                labels[i] = false;
            }
        }
    }
    if(sor_flag) {
        sor(tmp_vv, labels, nb_neighbors_, num_sor_iter_, std_radio_);
    }

    for(size_t i = 0; i < labels.size(); i++) {
        std::vector<int> tmp_indices;
        if(labels[i]) {
            vv.emplace_back(tmp_vv[i]);
            vvr.emplace_back(tmp_vvr[i]);
            vvns.emplace_back(tmp_vvns[i]);
        }
    }
}

void FilletSegV5::voronoi_density_drift(std::vector<easy3d::vec3>& vv
                                      , std::vector<float>& vvr
                                      , std::vector<std::vector<int>>& vvns
                                      , std::vector<std::vector<int>>& scvv) {


    int max_threads, num_threads;

    // 获取当前计算环境中可用的最大线程数
    max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    // printf("Max threads: %d\n", max_threads);
    std::string path = "../vv.ply";
    save_point_set(vv_, path);
    int num = s_.size();
    scvv.resize(num);
    num = vv.size();
    for(int i = 0; i < num; i++) {
        for(int j = 0; j < 4; j++) {
            int idx = vvns[i][j];
            scvv[idx].emplace_back(i);
        }
    }
    for(int iter = 0; iter < 10; iter++) {
#pragma omp parallel for
        for(int i = 0; i < num; i++) {
            std::set<int> unique;
            for(int j = 0; j < 4; j++) {
                int idx = vvns[i][j];
                unique.insert(scvv[idx].begin(), scvv[idx].end());
            }
            std::vector<easy3d::vec3> points;
            std::vector<int> indices;
            for(auto id : unique) {
                points.emplace_back(vv[id]);
                indices.emplace_back(id);
            }
            MeanShift ms(points);
            std::pair<easy3d::vec3, int> res = ms.run(vv[i], 1.0);
            int idx = res.second;
            idx = indices[idx];
            vv[i] = vv[idx];
            vvr[i] = vvr[idx];
        }
        path = "../vv_" + std::to_string(iter) +".ply";
        save_point_set(vv_, path);
    }
}


void FilletSegV5::farthest_point_sampling(int nb_samples, std::vector<int>& indices, float& max_gap) {
    int num = mesh_->n_faces();
    std::vector<float> dis(num, std::numeric_limits<float>::max());
    indices.clear();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> distribution(0, num - 1); // 生成1到100之间的随机整数

    for(int i = 0; i < nb_samples; i++) {
        std::priority_queue<std::pair<float, easy3d::SurfaceMesh::Face>> que;
        easy3d::SurfaceMesh::Face f;
        if(i == 0) {
            f = easy3d::SurfaceMesh::Face(distribution(gen));
        } else {
            float maxx = 0.0;
            int idx = 0;
            for(int j = 0; j < num; j++) {
                if(dis[j] > maxx) {
                    maxx = dis[j]; idx = j;
                }
            }
            f = easy3d::SurfaceMesh::Face(idx);
        }
        indices.emplace_back(f.idx());
        dis[f.idx()] = 0.0;
        que.push(std::make_pair(0.0, f));

        while(!que.empty()) {
            auto cur_dis = que.top().first;
            auto cur_f = que.top().second;
            que.pop();
            if(dis[cur_f.idx()] < -cur_dis) {
                continue;
            }
            dis[cur_f.idx()] = -cur_dis;
            for(auto h : mesh_->halfedges(cur_f)) {
                auto nxt_f = mesh_->face(mesh_->opposite(h));
                if(nxt_f.is_valid()) {
                    float val = (s_[cur_f.idx()] - s_[nxt_f.idx()]).norm();
                    if(dis[nxt_f.idx()] > -cur_dis + val) {
                        que.push(std::make_pair(cur_dis - val, nxt_f));
                    }
                }
            }
        }
    }

    max_gap = 0.0;
    for(int j = 0; j < num; j++) {
        if(dis[j] > max_gap) {
            max_gap = dis[j];
        }
    }
}

void FilletSegV5::voronoi3d(const std::vector<easy3d::vec3>& s
                         , std::vector<easy3d::vec3>& vv
                         , std::vector<float>& vvr
                         , std::vector<std::vector<int>>& snv
                         , std::vector<std::vector<int>>& vns) {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
    typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
    typedef Triangulation::Vertex_handle  Vertex_handle;
    typedef K::Point_3          CGAL_Point;

    int nb_sites = s.size();
    easy3d::Box3 box;

    for(auto p : s) {box.grow(p);}
    Triangulation T;

    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T.insert(CGAL_Point(s[i].x, s[i].y, s[i].z ));
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

    vv.clear(); snv.clear(); vns.clear();
    snv.resize(nb_sites);
    for(auto cell = T.finite_cells_begin(); cell != T.finite_cells_end(); cell++) {
        auto v = T.dual(cell);
        int vid = vv.size();
        vv.emplace_back(easy3d::vec3(v.x(), v.y(), v.z()));
        std::vector<int> tmp;
        float radius = 0;
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            snv[id].emplace_back(vid);
            radius += (s[id] - vv.back()).norm();
            tmp.emplace_back(id);
        }
        radius /= 4;
        vvr.emplace_back(radius);
        vns.emplace_back(tmp);
    }

}


double FilletSegV5::dihedral_angle(easy3d::SurfaceMesh* mesh
                      , easy3d::SurfaceMesh::Face f1
                      , easy3d::SurfaceMesh::Face f2
                      , bool rad) {
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


void FilletSegV5::crop_local_patch(int f, float max_gap
                                 , std::vector<int>& indices) {

    easy3d::SurfaceMesh::Face ff(f);
    std::set<int> vis;
    std::priority_queue<std::pair<float, easy3d::SurfaceMesh::Face>> que;
    que.push(std::make_pair(0.0, ff));
    while(!que.empty()) {
        auto cur_dis = que.top().first;
        auto cur_f = que.top().second;
        que.pop();
        vis.insert(cur_f.idx());
        for(auto h : mesh_->halfedges(cur_f)) {
            auto nxt_f = mesh_->face(mesh_->opposite(h));
            if(nxt_f.is_valid()) {
                float val = (s_[cur_f.idx()] - s_[nxt_f.idx()]).norm();
                if(max_gap > -cur_dis + val && vis.find(nxt_f.idx()) == vis.end()
                    && dihedral_angle(mesh_, cur_f, nxt_f) < angle_thr_) {
                    que.push(std::make_pair(cur_dis - val, nxt_f));
                }
            }
        }
    }
    indices.insert(indices.end(), vis.begin(), vis.end());
}

void FilletSegV5::sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels
       , int nb_neighbors, int num_sor_iter, float std_ratio) {
    int nb_points = points.size();
    // int nb_neighbors = 30, num_sor_iter = 3;
    // float std_ratio = 0.3;
    for(int i = 0; i < num_sor_iter; i++) {
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
            kds.kth_search(knn_points[j], nb_neighbors, tmp_indices, dist);
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
        double distance_threshold = cloud_mean + std_ratio * std_dev;
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


int FilletSegV5::calculate_samples_num(float radius) {
    float area = easy3d::geom::surface_area(mesh_);
    return int(std::ceil(area / (M_PI * radius * radius)));
}

void FilletSegV5::save_point_set(std::vector<easy3d::vec3>& points, std::string& path) {
    easy3d::PointCloud* cloud = new easy3d::PointCloud;
    for(auto p : points) {
        cloud->add_vertex(p);
    }
    easy3d::PointCloudIO::save(path, cloud);
}