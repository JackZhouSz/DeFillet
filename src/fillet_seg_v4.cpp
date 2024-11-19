//
// Created by xiaowuga on 2024/11/6.
//

#include "fillet_seg_v4.h"

#include <knn.h>
#include <knn4d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>

#include <easy3d/algo/surface_mesh_geometry.h>
#include <dkm_parallel.hpp>

FilletSegV4::FilletSegV4(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
    // fillet segmentation para
    eps_ = 0.03; sigma_ = 10; angle_thr_ = 40; radius_ = 0.08;
    nb_neighbors_ = 30;
    num_sor_iter_ = 2;
    std_radio_ = 1.0;
    // sites_normals_ = mesh_->face_property<easy3d::vec3>("f:normal");
    // sites_area_ = mesh_->face_property<float>("f:area");
    sites_area_.resize(mesh_->n_faces());
    sites_normals_.resize(mesh_->n_faces());
    for(auto f : mesh_->faces()) {
        sites_area_[f.idx()] = easy3d::geom::triangle_area(mesh_, f);
        sites_normals_[f.idx()] = mesh_->compute_face_normal(f);
        int ct = 0;
        easy3d::vec3 p(0,0,0);
        for(auto v : mesh_->vertices(f)) {
            p += mesh_->position(v); ct++;
        }
        if(ct)
            s_.emplace_back(p / ct);
    }
    box_ = mesh_->bounding_box();

    num_cluster_ = 1000;
}

void FilletSegV4::seg() {
    // std::vector<easy3d::vec3> vv;
    // std::vector<float> vvr;
    // std::vector<std::vector<int>> vvns;
    std::cout << "Start generate voronoi vertices..." <<std::endl;
    generate_local_voronoi_vertices(vv_, vvr_,vvns_, vvcs_,true);
    std::cout << "Generate voronoi vertices sussessfully!" <<std::endl;
    voronoi_vertices_drift(vv_, vvr_, vvcs_);
    // std::cout << "Start clustering..." <<std::endl;
    // clustering(vv_, vvr_, vvcs_, cc_, cr_, ccs_, ccvv_, ccvvn_, ccsn_, num_cluster_);
    // std::cout << "Clustering sussessfully!" <<std::endl;


}
void FilletSegV4::generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv, std::vector<float>& vvr
                                                , std::vector<std::vector<int>>& vvns
                                                , std::vector<std::vector<int>>& vvcs
                                                ,bool sor_flag, bool radius_thr) {
    int idx = 0;
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> mp_;
    std::vector<easy3d::vec3> tmp_vv;
    std::vector<float> tmp_vvr;
    std::vector<std::vector<int>> tmp_vvns;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 500 != 0) continue;
        std::vector<easy3d::SurfaceMesh::Face> patch;
        crop_local_patch(f, radius_, patch);

        std::vector<easy3d::vec3> ls; // local sites;
        for(size_t i = 0; i < patch.size(); i++) {
            ls.emplace_back(s_[patch[i].idx()]);
        }

        std::vector<std::vector<int>> lsnv;   // local sites neighboring vertices;
        std::vector<float> lvvr; // vor_vertices_radius
        std::vector<std::vector<int>> lvns ; // local vertices neighboring sites;
        std::vector<easy3d::vec3> lvv; // local voronoi vertices;
        voronoi3d(ls, lvv, lvvr,lsnv, lvns);

        for(size_t i = 0; i < lvv.size(); i++) {
            std::vector<int> indices(4);
            for(int j = 0; j < 4; j++) {
                indices[j] = patch[lvns[i][j]].idx();
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
        float thr = box_.diagonal_length() * radius_;
        for(size_t i = 0; i < tmp_vvr.size(); i++) {
            if(tmp_vvr[i] > thr) {
                labels[i] = false;
            }
        }
    }
    if(sor_flag) {
        sor(tmp_vv, labels, nb_neighbors_, num_sor_iter_, std_radio_);
    }

    for(size_t i = 0; i < labels.size(); i++) {
        std::vector<int> tmp_indices;
        if(labels[i] && search_voronoi_vertex_corresponding_sites(tmp_vv[i], tmp_vvr[i], tmp_vvns[i],tmp_indices)) {
            vv.emplace_back(tmp_vv[i]);
            vvr.emplace_back(tmp_vvr[i]);
            vvcs.emplace_back(tmp_indices);
            vvns.emplace_back(tmp_vvns[i]);
        }
    }
}

void FilletSegV4::voronoi_vertices_drift(std::vector<easy3d::vec3>& vv
                                       , std::vector<float>& vvr
                                       , std::vector<std::vector<int>>& vvcs) {
    int num = vv.size();
    std::cout << num << std::endl;
    double area_thr = box_.diagonal_length() * radius_;
    for(int i = 0; i < num; i++) {
        if(i % 10 == 0) {
            std::cout << i << std::endl;
        }
        double area_sum = 0;
        std::set<int> vis;
        std::queue<easy3d::SurfaceMesh::Face> que;
        for(size_t j = 0; j < vvcs[i].size(); j++) {
            vis.insert(vvcs[i][j]);
            easy3d::SurfaceMesh::Face f(vvcs[i][j]);
            que.push(f);
            area_sum += sites_area_[f.idx()];

        }
        while(area_sum < area_thr && !que.empty()) {
            auto cur_f = que.front(); que.pop();
            area_sum += sites_area_[cur_f.idx()];
            for(auto h : mesh_->halfedges(cur_f)) {
                auto nxt_f = mesh_->face(mesh_->opposite(h));
                if(nxt_f.is_valid() && vis.find(nxt_f.idx()) == vis.end()
                   && dihedral_angle(mesh_, cur_f, nxt_f) < angle_thr_) {
                    vis.insert(nxt_f.idx());
                    que.push(nxt_f);
                }
            }
        }

        std::vector<easy3d::vec3> ls;
        for(auto id : vis) {
            ls.emplace_back(s_[id]);
        }
        std::vector<std::vector<int>> lsnv;   // local sites neighboring vertices;
        std::vector<float> lvvr; // vor_vertices_radius
        std::vector<std::vector<int>> lvns ; // local vertices neighboring sites;
        std::vector<easy3d::vec3> lvv; // local voronoi vertices;
        voronoi3d(ls, lvv, lvvr,lsnv, lvns);

        int ct = 0;
        std::vector<KNN::Point> knn_points;
        for(int j = 0; j < lvv.size(); j++) {
            if(box_.contains(lvv[j])) {
                knn_points.emplace_back(KNN::Point(lvv[j].x, lvv[j].y, lvv[j].z));
                ct++;
            }
        }
        if(ct < 20) {
            continue;
        }
        KNN::KdSearch kds(knn_points);
        double minn = 1e8;
        int kk = 5;
        for(int j = 0; j < ct; j++) {
            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(knn_points[j], kk, tmp_indices, dist);
            double avg_dis = std::accumulate(dist.begin(), dist.end(), 0.0) / kk;
            if(avg_dis < minn) {
                vv[i] = easy3d::vec3(knn_points[j].p[0], knn_points[j].p[1], knn_points[j].p[2]);
                minn = avg_dis;
            }
        }
    }
}

void FilletSegV4::clustering(std::vector<easy3d::vec3>& vv
                           , std::vector<float>& vvr
                           , std::vector<std::vector<int>>& vvcs
                           , std::vector<easy3d::vec3>& c
                           , std::vector<float>& cr
                           , std::vector<std::vector<int>>& ccs
                           , std::vector<std::vector<int>>& ccvv
                           , std::vector<float>& ccvvn
                           , std::vector<float>& ccsn
                           , int num_cluster) {
    int num = vv.size();
    std::vector<std::array<float, 4>> data;
    for(int i = 0; i < num; i++) {
        easy3d::vec3& p = vv[i];
        data.emplace_back(std::array<float, 4>{p.x, p.y, p.z, vvr[i]});
    }
    auto cluster_data = dkm::kmeans_lloyd_parallel(data, num_cluster);
    for (const auto& mean : std::get<0>(cluster_data)) {
        c.emplace_back(easy3d::vec3(mean[0],mean[1],mean[2]));
        cr.emplace_back(mean[3]);
    }
    ccs.resize(num_cluster);
    ccvv.resize(num_cluster);
    ccvvn.resize(num_cluster);
    ccsn.resize(num_cluster);
    int idx = 0;
    for (const auto& label : std::get<1>(cluster_data)) {
        ccvv[label].emplace_back(idx++);
    }
    for(size_t i = 0; i < ccvv.size(); i++) {
        std::set<int> tmp_s;
        ccvvn[i] = ccvv[i].size();
        for(size_t j = 0; j < ccvv[i].size(); j++) {
            int vvid = ccvv[i][j];
            for(size_t k = 0; k < vvcs[vvid].size(); k++) {
                tmp_s.insert(vvcs[vvid][k]);
            }
        }
        ccs[i] = std::vector<int>(tmp_s.begin(), tmp_s.end());
        ccsn[i] = ccs[i].size();
    }
}


void FilletSegV4::crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
                        , std::vector<easy3d::SurfaceMesh::Face>& patch) {
    radius = box_.diagonal_length() * radius;
    double area_thr = 2 * M_PI * radius * radius;
    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face> que;
    que.push(f); vis.insert(f);
    double area_sum = 0;
    while(area_sum < area_thr && !que.empty()) {
        auto cur_f = que.front(); que.pop();
        area_sum += sites_area_[cur_f.idx()];
        for(auto h : mesh_->halfedges(cur_f)) {
            auto nxt_f = mesh_->face(mesh_->opposite(h));
            if(nxt_f.is_valid() && vis.find(nxt_f) == vis.end()
               && dihedral_angle(mesh_, cur_f, nxt_f) < angle_thr_) {
                vis.insert(nxt_f);
                que.push(nxt_f);
               }
        }
    }
    patch.insert(patch.end(), vis.begin(), vis.end());
}


void FilletSegV4::voronoi3d(const std::vector<easy3d::vec3>& s
                         , std::vector<easy3d::vec3>& vv
                         , std::vector<float>& vvr
                         , std::vector<std::vector<int>>& snv
                         , std::vector<std::vector<int>>& vns) {

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


void FilletSegV4::sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels
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


double FilletSegV4::dihedral_angle(easy3d::SurfaceMesh* mesh
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


bool FilletSegV4::search_voronoi_vertex_corresponding_sites(easy3d::vec3& vv, float r,std::vector<int>& vvns
                                                  , std::vector<int>& vvcs) {
    std::queue<easy3d::SurfaceMesh::Face> que;
    std::set<int> vis;
    for(size_t i = 0; i < vvns.size(); i++) {
        vis.insert(vvns[i]);
        que.push(easy3d::SurfaceMesh::Face(vvns[i]));
    }
    double tol = r * eps_;
    while(!que.empty()) {
        auto cur_f = que.front();
        que.pop();
        for (auto h: mesh_->halfedges(cur_f)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur_f, opp_f) < angle_thr_) {
                double len = (s_[opp_f.idx()] - vv).norm();
                if (vis.find(opp_f.idx()) == vis.end() && len < r + tol && len > r - tol) {
                    que.push(opp_f);
                    vis.insert(opp_f.idx());
                }
            }
        }
    }
    std::set<int> vis2;
    vis2.insert(vvns[0]);
    que.push(easy3d::SurfaceMesh::Face(vvns[0]));
    while(!que.empty()) {
        auto cur_f = que.front();
        que.pop();
        for (auto h: mesh_->halfedges(cur_f)) {
            auto opp_f = mesh_->face(mesh_->opposite(h));
            if (opp_f.is_valid() && dihedral_angle(mesh_, cur_f, opp_f) < angle_thr_) {
                double len = (s_[opp_f.idx()] - vv).norm();
                if (vis2.find(opp_f.idx()) == vis2.end() && len < r + tol && len > r - tol) {
                    que.push(opp_f);
                    vis2.insert(opp_f.idx());
                }
            }
        }
    }
    if(vis2.size() != vis.size()) {
        return false;
    }
    vvcs.insert(vvcs.end(), vis.begin(), vis.end());
    return true;
}


void FilletSegV4::fit_sphere(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius, int max_iter) {
    int num = points.size();
    std::vector<int> numbers(num);
    std::iota(numbers.begin(), numbers.end(), 0);  // 填充范围内的整数
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::vector<easy3d::vec3> samples(4);
    float loss = 1e9;
    for(int i = 0; i < max_iter; i++) {
        std::shuffle(numbers.begin(), numbers.end(),std::default_random_engine(seed));
        for(int j = 0; j < 4; j++) {
            samples[j] = points[numbers[j]];
        }
        easy3d::vec3 center_candidate;
        float radius_candidate = -1;

        tetrahedron_circumcenter(samples, center_candidate, radius_candidate);
        float err = 0;
        for(int j = 0; j < num; j++) {
            err += fabs((points[j] - center_candidate).norm() - radius_candidate);
        }
        err /= num;
        if(err < loss) {
            loss = err;
            center = center_candidate;
            radius = radius_candidate;
        }
    }
}


void FilletSegV4::tetrahedron_circumcenter(std::vector<easy3d::vec3>& points, easy3d::vec3& center, float& radius) {
    int num = points.size();
    if(num != 4) {
        std::cerr << "points of tetrahedron do not equal 4" << std::endl;
        return;
    }
    Eigen::MatrixXd d_matrix(4, 4); d_matrix.setOnes();
    // Multiplied by (x²+y²+z²)
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = points[i].x;
        d_matrix(i, 1) = points[i].y;
        d_matrix(i, 2) = points[i].z;
    }
    float M11 = d_matrix.determinant();
    // Multiplied by x
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].y;
        d_matrix(i, 2) = points[i].z;
    }
    float M12 = d_matrix.determinant();

    // std::cout << M11 << std::endl;
    // std::cout << "bb" << std::endl;

    // Multiplied by y
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].z;
    }
    float M13 = d_matrix.determinant();

    // Multiplied by z
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].y;
    }
    float M14 = d_matrix.determinant();

    // Multiplied by 1
    for(int i  = 0; i < 4; i++) {
        d_matrix(i, 0) = easy3d::dot(points[i], points[i]);
        d_matrix(i, 1) = points[i].x;
        d_matrix(i, 2) = points[i].y;
        d_matrix(i, 3) = points[i].z;
    }

    float M15 = d_matrix.determinant();

    center = easy3d::vec3(0.5 * (M12 / M11), -0.5 * (M13 / M11), 0.5 * (M14 / M11));
    radius = std::sqrt(easy3d::dot(center, center) - (M15 / M11));

}