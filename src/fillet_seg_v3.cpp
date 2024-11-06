//
// Created by xiaowuga on 2024/10/31.
//

#include "fillet_seg_v3.h"

#include <knn.h>
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


void FilletSegV3::sor(const std::vector<easy3d::vec3>& points, std::vector<bool>& labels) {
    int nb_points = points.size();
    int nb_neighbors = 30, num_sor_iter = 3;
    float std_ratio = 0.3;
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

bool FilletSegV3::find_cluster_center(const std::vector<easy3d::vec3>& data, easy3d::vec3& centroid) {
    int num = data.size();
    if(num < 20) {
        // std::cout << "point too few" <<std::endl;
        return false;
    }
    std::vector<bool> labels(num, true);
    sor(data, labels);
    centroid = easy3d::vec3(0,0,0);
    std::vector<easy3d::vec3> clean_data;
    for(int i = 0; i < num; i++) {
        if(labels[i]) {
            centroid += data[i];
            clean_data.emplace_back(data[i]);
        }
    }

    num = clean_data.size();
    if(num < 20) {
        return false;
    }
    centroid /= num;

    std::vector<float> dist(num);
    for(int i = 0; i < num; i++) {
        dist[i] = (clean_data[i] - centroid).length2();
    }

    std::vector<float> tmp = dist;
    std::nth_element(tmp.begin(), tmp.end(), tmp.begin() + num / 2);
    double med = tmp[num / 2];
    for(int i = 0; i < num; i++) {
        tmp[i] = fabs(dist[i] - med);
    }
    std::nth_element(tmp.begin(), tmp.end(), tmp.begin() + num / 2);
    double mad = tmp[num / 2];
    easy3d::vec3 new_centroid(0,0,0);
    int ct = 0;
    // std::vector<easy3d::vec3> dd;
    for(int i = 0; i < num; i++) {
        // if(dist[i] - med < 1.4 * mad) {
            new_centroid += clean_data[i];
            // dd.emplace_back(clean_data[i]);
            ct++;
        // }
    }
    // std::cout << ct << ' ' << new_centroid << std::endl;
    if(ct > 1) {
        new_centroid /= ct;
        centroid = new_centroid;
        return true;
    }
    else {
        return false;
    }

    // std::cout << xx << std::endl;

}

FilletSegV3::FilletSegV3(easy3d::SurfaceMesh* mesh) : mesh_(mesh){
    // fillet segmentation para
    eps_ = 0.005; sigma_ = 10; angle_thr = 60; radius_ = 0.08;
    nb_neighbors_ = 50;
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
            sites_.emplace_back(p / ct);
    }
    box_ = mesh_->bounding_box();
}

void FilletSegV3::voronoi3d(const std::vector<easy3d::vec3>& sites
                 , std::vector<easy3d::vec3>& vor_vertices
                 , std::vector<float>& vor_vertices_radius
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
        float radius = 0;
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            sites_nearby_vertices[id].emplace_back(vid);
            radius += (sites[id] - vor_vertices.back()).norm();
            tmp.emplace_back(id);
        }
        radius /= 4;
        vor_vertices_radius.emplace_back(radius);
        vertices_nearby_sites.emplace_back(tmp);
    }

}

void FilletSegV3::crop_local_patch(easy3d::SurfaceMesh::Face f, double radius
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
               && dihedral_angle(mesh_, cur_f, nxt_f) < angle_thr) {
                vis.insert(nxt_f);
                que.push(nxt_f);
            }
        }
    }
    patch.insert(patch.end(), vis.begin(), vis.end());
}

void FilletSegV3::generate_local_voronoi_vertices(std::vector<easy3d::vec3>& vv, std::vector<float>& vvr, std::vector<std::vector<int>>& vvns, bool sor_flag) {
    // gvv: global voronoi vertices
    // gvvns: global voronoi vertices nearby sites
    int idx = 0;
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> mp_;
    std::vector<easy3d::vec3> tmp_vv;
    std::vector<float> tmp_vvr;
    std::vector<std::vector<int>> tmp_vvns;
    for(auto f : mesh_->faces()) {
        if(f.idx() % 150 != 0) continue;
        std::vector<easy3d::SurfaceMesh::Face> patch;
        crop_local_patch(f, radius_, patch);

        std::vector<easy3d::vec3> ls; // local sites;
        for(size_t i = 0; i < patch.size(); i++) {
            ls.emplace_back(sites_[patch[i].idx()]);
        }

        std::vector<std::vector<int>> lsnv;   // local sites_nearby_vertices;
        std::vector<float> lvvr; // vor_vertices_radius
        std::vector<std::vector<int>> lvns ; // local vertices_nearby_sites;
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

    if(sor_flag) {
        std::vector<bool> labels(tmp_vv.size(), true);
        sor(tmp_vv, labels);
        for(size_t i = 0; i < labels.size(); i++) {
            if(labels[i]) {
                vv.emplace_back(tmp_vv[i]);
                vvr.emplace_back(tmp_vvr[i]);
                vvns.emplace_back(tmp_vvns[i]);
            }
        }
    }else {
        vv = tmp_vv;
        vvr = tmp_vvr;
        vvns = tmp_vvns;
    }
}

void FilletSegV3::generate_global_voronoi_vertices(std::vector<easy3d::vec3>& vv, std::vector<std::vector<int>>& vvns) {

}

void FilletSegV3::search_voronoi_vertices_corresponding_sites(std::vector<easy3d::vec3>& vv, std::vector<std::vector<int>>& vvns
                                               , std::vector<std::vector<int>>& vvcs) {
    int num = vv.size();
    for(int i = 0; i < num; i++) {
        std::set<int> vis;
        std::queue<easy3d::SurfaceMesh::Face> que;
        for(size_t j = 0; j < vvns[i].size(); j++) {
            vis.insert(vvns[i][j]);
            que.push(easy3d::SurfaceMesh::Face(vvns[i][j]));
        }
        double R = (sites_[vvns[i][0]] - vv[i]).norm();
        double tol = R * eps_;
        if(R > mesh_->bounding_box().diagonal_length() * 0.5) {
            vvcs.emplace_back(std::vector<int>());
            continue;;
        }
        while(!que.empty()) {
            auto cur_f = que.front();
            que.pop();
            for (auto h: mesh_->halfedges(cur_f)) {
                auto opp_f = mesh_->face(mesh_->opposite(h));
                if (opp_f.is_valid() && dihedral_angle(mesh_, cur_f, opp_f) < 30) {
                    double len = (sites_[opp_f.idx()] - vv[i]).norm();
                    if (vis.find(opp_f.idx()) == vis.end() && len < R + tol && len > R - tol) {
                        que.push(opp_f);
                        vis.insert(opp_f.idx());
                    }
                }
            }
        }

        vvcs.emplace_back(std::vector<int>());
        vvcs.back().insert(vvcs.back().end(), vis.begin(), vis.end());
    }
}


void FilletSegV3::search_sites_corresponding_voronoi_vertices(std::vector<std::vector<int>>& vvcs, std::vector<std::vector<int>>& scvv) {
    int num = sites_.size();
    scvv.resize(num);
    for(size_t i = 0; i < vvcs.size(); i++) {
        for(size_t j = 0;  j < vvcs[i].size(); j++) {
            int idx = vvcs[i][j];
            scvv[idx].emplace_back(i);
        }
    }
}


void FilletSegV3::compute_voronoi_vertices_score(std::vector<easy3d::vec3>& vv, std::vector<float>& vvr, std::vector<float>& vvs) {
    int num = vv.size();
    std::vector<KNN::Point> knn_points;
    for(int i = 0; i < num; i++) {
        knn_points.emplace_back(KNN::Point(vv[i].x, vv[i].y, vv[i].z));
    }
    KNN::KdSearch kdt(knn_points);
    vvs.resize(vv.size(), 0);
    std::vector<float> count(vv.size(), 0);
    for(int i = 0; i < num; i++) {
        std::vector<size_t> tmp_indices;
        std::vector<double> dist;
        kdt.kth_search(knn_points[i], nb_neighbors_, tmp_indices, dist);
        for(int j = 0; j < tmp_indices.size(); j++) {
            float v_dis = (vv[i] - vv[tmp_indices[j]]).norm();
            float r_dis = fabs(vvr[i] - vvr[tmp_indices[j]]);
            if(v_dis  <= r_dis) {
                count[i] += 1.0;
            }
        }
    }
    // vvs = count;
    for(int i = 0; i < num; i++) {
        // std::cout <<i << ' ' <<  count[i] << std::endl;
        std::vector<size_t> tmp_indices;
        std::vector<double> dist;
        kdt.kth_search(knn_points[i], nb_neighbors_, tmp_indices, dist);
        std::vector<std::pair<float,int>> tmp;
        for(int j = 0; j < tmp_indices.size(); j++) {
            tmp.emplace_back(std::make_pair(count[tmp_indices[j]],tmp_indices[j]));
        }
        std::sort(tmp.begin(), tmp.end());
        for(size_t j = 0; j < tmp.size(); j++) {
            if(tmp[j].second == i) {
                if(tmp[j].first > 1) {
                    vvs[i] = 1.0 * (j + 1) / nb_neighbors_;
                }
                else vvs[i] = 0.0;
                break;
            }
        }
    }

}

void FilletSegV3::esimate_sites_radius(std::vector<easy3d::vec3>& sites, std::vector<easy3d::vec3>& vv
                        ,std::vector<float>& vvr, std::vector<std::vector<int>>& scvv
                        ,std::vector<float>& sr) {
    int num = sites.size();
    sr.resize(num);
    for(int i = 0; i < num; i++) {
        // std::cout << "start: " <<i << std::endl;
        std::vector<easy3d::vec3> scvv_i;
        std::vector<float> raidus_i;
        for(size_t j = 0; j < scvv[i].size(); j++) {
            scvv_i.emplace_back(vv[scvv[i][j]]);
            raidus_i.emplace_back((vv[scvv[i][j]] - sites[i]).norm());
        }
        easy3d::vec3 center;
        // if(find_cluster_center(scvv_i, center)) {
        //     sr[i] = (sites[i] - center).norm();
        // }
        // else {
        //     // sr[i] = mesh_->bounding_box().radius();
        //     sr[i] = -1.0;
        // }
        if(find_centroid_via_kmeans(scvv_i, raidus_i,center)) {
            sr[i] = (sites[i] - center).norm();
        }else {
            sr[i] = -1.0;
            if(scvv_i[i].size() == 1) {
                sr[i] = 0.0;
            }
        }
        // std::cout << "done: " <<i << std::endl;
    }
}



void FilletSegV3::solve() {
    generate_local_voronoi_vertices(vv_, vvr_,vvns_, true);
    search_voronoi_vertices_corresponding_sites(vv_, vvns_, vvcs_);
    search_sites_corresponding_voronoi_vertices(vvcs_, scvv_);
    // compute_voronoi_vertices_score(vv_, vvr_, vvs_);
    esimate_sites_radius(sites_, vv_, vvr_, scvv_, sr_);
    // count_span_and_radius_of_voronoi_vertices(vv_,vvns_,vvcs_, sar_);
    // std::ofstream outfile("../aa.txt");
    // for(size_t i = 0; i < sar_.size(); i++) {
    //     outfile << sar_[i].first << ' ' << sar_[i].second << std::endl;
    // }
    // outfile.close();
}



void FilletSegV3::kmeans(const Eigen::MatrixXd& data, int num_clusters, Eigen::MatrixXd& centroids, Eigen::VectorXi& labels) {
    // 随机初始化聚类中心
    int num = data.rows();
    Eigen::VectorXd distances = Eigen::VectorXd::Constant(num, std::numeric_limits<float>::max());
    centroids.resize(num_clusters, data.cols());
    Eigen::VectorXd cur = data.colwise().sum() / num;
    for(int i = 0; i < num_clusters; i++) {
        float maxx = 0;
        int idx = 0;
        for(int j = 0; j < num; j++) {
            double dis = (data.row(j) - cur).norm();
            if(distances[j] > dis) {
                distances[j] = dis;
            }
            if(maxx < distances[j]) {
                maxx = distances[j];
                idx = j;
            }
        }
        cur = data.row(idx);
        centroids.row(i) = cur;
    }
    labels.resize(num);
    // return;
    distances.resize(num_clusters);
    const int max_iterations = 10;
    for (int iter = 0; iter < max_iterations; ++iter) {
        // std::cout << iter << std::endl;
        // 分配每个样本到最近的聚类中心
        for (int i = 0; i < data.rows(); ++i) {
            float minn = 1e9;
            for (int k = 0; k < num_clusters; ++k) {
                distances(k) = (data.row(i) - centroids.row(k)).squaredNorm();
                if(minn > distances[k]) {
                    labels[i] = k;
                    minn = distances[k];
                }
            }
            // labels(i) = distances.minCoeff();
        }

        // 更新聚类中心
        Eigen::MatrixXd new_centroids = Eigen::MatrixXd::Zero(num_clusters, data.cols());
        Eigen::VectorXi cluster_sizes = Eigen::VectorXi::Zero(num_clusters);

        for (int i = 0; i < data.rows(); ++i) {
            new_centroids.row(labels(i)) += data.row(i);
            cluster_sizes(labels(i))++;
        }

        for (int k = 0; k < num_clusters; ++k) {
            if (cluster_sizes(k) > 0) {
                centroids.row(k) = new_centroids.row(k) / cluster_sizes(k);
            } else {
                centroids.row(k) = data.row(rand() % num);
            }
        }
    }
}

bool FilletSegV3::find_centroid_via_kmeans(std::vector<easy3d::vec3>& points, std::vector<float>& radius, easy3d::vec3& centroid) {
    int nb_points = points.size();
    if(nb_points < 2) {
        return false;
    }
    Eigen::MatrixXd data(nb_points, 4);
    Eigen::MatrixXd centroids(2, 4); centroids.setZero();
    for(int i = 0; i < nb_points; i++) {
        data.row(i) = Eigen::RowVector4d(points[i].x, points[i].y, points[i].z, radius[i]);
        centroids.row(0) += data.row(i);
    }
    centroids /= nb_points;
    float maxx = 0;
    for(int i = 0; i < nb_points; i++) {
        float dis = (data.row(i) - centroids.row(0)).norm();
        if(dis > maxx) {
            centroids.row(1) = data.row(i);
            maxx = dis;
        }
    }
    Eigen::VectorXi labels;

    kmeans(data, 2, centroids, labels);
    int ct0 = 0, ct1 = 0;
    for(int i = 0; i < labels.size(); i++) {
        if(labels[i] == 0) {
            ct0++;
        }
        else {
            ct1++;
        }
    }

    if(ct0 > ct1) {
        centroid = easy3d::vec3(centroids(0,0), centroids(0,1), centroids(0,2));
    }else {
        centroid = easy3d::vec3(centroids(1,0), centroids(1,1), centroids(1,2));
    }
}