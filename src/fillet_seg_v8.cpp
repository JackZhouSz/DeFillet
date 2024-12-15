//
// Created by xiaowuga on 2024/12/9.
//

#include "fillet_seg_v8.h"

#include "kernel.h"
#include "knn4d.h"
#include "knn.h"
#include "io.h"
#include "gcp.h"

#include <random>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>

#include <easy3d/algo/surface_mesh_geometry.h>
#include <easy3d/util/stop_watch.h>

#include <iostream>


FilletSegV8::FilletSegV8(easy3d::SurfaceMesh* mesh
                       , float eps
                       , int num_samples
                       , float radius_thr
                       , float angle_thr
                       , bool local_voronoi
                       , bool sor_filter
                       , bool radius_filter
                       , int num_sor_iter
                       , int num_neighbors
                       , float sor_std_radio
                       , float h
                       , int num_smooth_iter
                       , float lamdba)
                       :   mesh_(mesh), eps_(eps), num_samples_(num_samples)
                         , angle_thr_(angle_thr), local_voronoi_(local_voronoi)
                         , radius_filter_(radius_filter), h_(h),  sor_filter_(sor_filter)
                         , num_sor_iter_(num_sor_iter), num_neighbors_(num_neighbors), num_smooth_iter_(num_smooth_iter)
                         , sor_std_radio_(sor_std_radio), lamdba_(lamdba){
    int num = mesh->n_faces();
    s_.resize(num);
    sn_.resize(num);
    sg_.resize(num);

    for(auto f : mesh->faces()) {
        s_[f.idx()] = easy3d::geom::centroid(mesh, f);
        sn_[f.idx()] = mesh->compute_face_normal(f);
    }

    for(auto f : mesh->faces()) {
        for(auto he : mesh->halfedges(f)) {
            auto of = mesh->face(mesh->opposite(he));
            if(of.is_valid()) {
                sg_[f.idx()].emplace_back(of.idx());
            }
        }
    }

    for(auto e : mesh->edges()) {
        auto f0 = mesh_->face(e, 0);
        auto f1 = mesh_->face(e, 1);
        if(f0.is_valid() && f1.is_valid()) {
            int x = f0.idx(), y = f1.idx();
            if(x > y) {
                std::swap(x, y);
            }
            sge_.emplace_back(std::make_pair(x, y));
            sgel_.emplace_back(mesh_->edge_length(e));
        }
    }

    easy3d::Box3 box = easy3d::geom::bounding_box<easy3d::Box3>(s_);
    radius_thr_ = box.diagonal_length() * radius_thr;
    eps_ = box.diagonal_length() * eps;
}

void FilletSegV8::run() {
    generate_voronoi_vertices(vv_, vvr_, vvns_, vvcs_, scvv_, vv_density_);
    if(vv_.size() != 0) {
        marching_axis_transform(mavv_, mavvns_, mavvcs_, mavvr_, scma_);
        run_graph_cut(ss_, fillet_labels_);
        std::cout << "Fillet segmentation total time="<< gvv_time_ + mat_time_ + gcp_time_ <<std::endl;
    }
}


void FilletSegV8::generate_voronoi_vertices(std::vector<easy3d::vec3>& vv
                                          , std::vector<float>& vvr
                                          , std::vector<std::vector<int>>& vvns
                                          , std::vector<std::vector<int>>& vvcs
                                          , std::vector<std::vector<int>>& scvv
                                          , std::vector<float>& vv_density) {
    // std::cout << radius_thr << std::endl;
    std::cout << "Start generate voronoi vertices..." <<std::endl;
    easy3d::StopWatch sw;
    std::vector<easy3d::vec3> tmp_vv;
    std::vector<float> tmp_vvr;
    std::vector<std::vector<int>> tmp_vvns;
    if(local_voronoi_) {
        generate_voronoi_vertices_via_local_patch(tmp_vv,tmp_vvr,tmp_vvns);
    }
    else {
        generate_voronoi_vertices_via_global_sites(tmp_vv,tmp_vvr,tmp_vvns);
    }

    int num = tmp_vv.size();
    std::vector<bool> labels(tmp_vv.size(), true);
    if(radius_filter_) {
#pragma omp parallel for
        for(int i = 0; i < num; i++) {
            if(tmp_vvr[i] > radius_thr_) {
                labels[i] = false;
            }
        }
    }
    if(sor_filter_)
        sor(tmp_vv, labels,num_neighbors_, num_sor_iter_, sor_std_radio_);

    omp_lock_t lock;
    omp_init_lock(&lock);
#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        std::vector<int> indices;
        if(labels[i] &&  search_vvcs(tmp_vv[i], tmp_vvr[i], tmp_vvns[i], indices)) {
                omp_set_lock(&lock);
                vv.emplace_back(tmp_vv[i]);
                vvr.emplace_back(tmp_vvr[i]);
                vvns.emplace_back(tmp_vvns[i]);
                vvcs.emplace_back(indices);
                omp_unset_lock(&lock);
        }
    }
    num = vv.size();
    scvv.resize(s_.size());

#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        for(size_t j = 0; j < vvns[i].size(); j++) {
            int id = vvns[i][j];
            omp_set_lock(&lock);
            scvv[id].emplace_back(i);
            omp_unset_lock(&lock);
        }
    }
    omp_destroy_lock(&lock);


    int vv_num = vv_.size();
    std::vector<KNN::Point4D> data;
    std::vector<easy3d::vec4> vv4d;
    for(int i = 0; i < vv_num; i++) {
        data.emplace_back(KNN::Point4D(vv[i].x,vv[i].y, vv[i].z, vvr[i]));
        vv4d.emplace_back(easy3d::vec4(vv[i].x,vv[i].y, vv[i].z, vvr[i]));
    }
    KNN::KdSearch4D kds(data);
    vv_density.resize(vv_num);

#pragma omp parallel for
    for(int i = 0; i < vv_num; i++) {
        std::vector<size_t> indices;
        std::vector<double> dist;
        kds.kth_search(data[i], num_neighbors_, indices, dist);
        vv_density[i] = 0.0;
        double w = 0;
        for(size_t j = 0; j < indices.size(); j++) {
            int idx = indices[j];
            double dis = fabs((vv4d[idx] - vv4d[i]).norm());
            double val = gaussian_kernel(dis, 1.0);
            vv_density[i] += val * dis;
            w += val;
        }
        vv_density[i] /= w;
    }
    gvv_time_ = sw.elapsed_seconds(5);
    if(vv.size() != 0) {
        std::cout << "Generate voronoi vertices sussessfully! time="<<gvv_time_ <<std::endl;
        // save_point_field(vv, vvr, "../out/vv.ply");
    }
    else {
        std::cout << "Generate voronoi vertices failly! time="<<gvv_time_ <<std::endl;
        std::cout << "Increasing the value of radius_thr may solve the problem." << std::endl;
    }

}


void FilletSegV8::generate_voronoi_vertices_via_global_sites(std::vector<easy3d::vec3>& vv
                                                           , std::vector<float>& vvr
                                                           , std::vector<std::vector<int>>& vvns) {
    easy3d::Box3 box = easy3d::geom::bounding_box<easy3d::Box3>(s_);
    std::vector<std::vector<int>> tmp_snvv;   // sites neighboring vertices;
    std::vector<float> tmp_vvr; // vor_vertices_radius
    std::vector<std::vector<int>> tmp_vvns ; // local vertices neighboring sites;
    std::vector<easy3d::vec3> tmp_vv; // local voronoi vertices;
    voronoi3d(s_, tmp_vv, tmp_vvr,tmp_snvv, tmp_vvns);
    int num = tmp_vv.size();
    for(int i = 0; i < num; i++) {
        if(box.contains(tmp_vv[i])) {
            vv.emplace_back(tmp_vv[i]);
            vvr.emplace_back(tmp_vvr[i]);
            vvns.emplace_back(tmp_vvns[i]);
        }
    }
}

void FilletSegV8::generate_voronoi_vertices_via_local_patch(std::vector<easy3d::vec3>& vv
                                                          , std::vector<float>& vvr
                                                          , std::vector<std::vector<int>>& vvns) {

    int idx = 0;
    easy3d::Box3 box = easy3d::geom::bounding_box<easy3d::Box3>(s_);
    std::map<int, std::map<int, std::map<int,std::map<int,int>>>> unique;

    std::vector<int> ssi; // sampling sites indices
    farthest_point_sampling(num_samples_, s_, sg_, ssi);

    float max_gap = 0.5 * M_PI * radius_thr_;
    omp_lock_t lock;
    omp_init_lock(&lock);
#pragma omp parallel for
    for(int k = 0; k < num_samples_; k++) {
        std::vector<int> patch;

        crop_local_patch(ssi[k], max_gap, patch);

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
            if(unique.find(indices[0]) != unique.end()) {
                auto& sub_mp1 = unique[indices[0]];
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
            if(flag && box.contains(lvv[i])) {
                omp_set_lock(&lock);
                vv.emplace_back(lvv[i]);
                vvr.emplace_back(lvvr[i]);
                vvns.emplace_back(indices);
                unique[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                omp_unset_lock(&lock);
            }
        }
    }

}

void FilletSegV8::marching_axis_transform(std::vector<easy3d::vec3>& mavv
                                        , std::vector<std::vector<int>>& mavvns
                                        , std::vector<std::vector<int>>& mavvcs
                                        , std::vector<float>& mavvr
                                        , std::vector<int>& scma) {
    std::cout << "Start marching axis transform..." <<std::endl;
    easy3d::StopWatch sw;
     int vv_num = vv_.size();

    int s_num = s_.size();
    std::vector<float> smd(s_num, std::numeric_limits<float>::max()); // sites minimum density
    std::vector<int> smavv(s_num, -1); // sites marching axis voronoi vertices
    for(int i = 0; i < vv_num; i++) {
        for(size_t j = 0; j < vvcs_[i].size(); j++) {
            int idx = vvcs_[i][j];
            if(smd[idx] > vv_density_[i]) {
                smd[idx] = vv_density_[i];
                smavv[idx] = i;
            }
        }
    }


    std::set<int> unique;
    for(int i = 0; i < s_num; i++) {
        if(smavv[i] != -1) {
            unique.insert(smavv[i]);
        }
    }
    std::map<int, int> mp;

    for(auto i : unique) {
        mp[i] = mavv.size();
        mavv.emplace_back(vv_[i]);
        mavvns.emplace_back(vvns_[i]);
        mavvcs.emplace_back(vvcs_[i]);
        mavvr.emplace_back(vvr_[i]);
    }

    scma.resize(s_num);
    std::fill(scma.begin(), scma.end(), -1);
    for(int i = 0; i < s_num; i++) {
        if(smavv[i] != -1) {
            scma[i] = mp[smavv[i]];
        }
    }
    mat_time_ = sw.elapsed_seconds(5);
    std::cout << "Marching axis transform sussessfully! time="<< mat_time_ <<std::endl;

}

void FilletSegV8::run_graph_cut(std::vector<float>& score, std::vector<int>& fillet_labels) {
    std::cout << "Start perfrom graph cut..." <<std::endl;
    easy3d::StopWatch sw;
    int s_num = s_.size();
    score.resize(s_num);
#pragma omp parallel for
    for(int i = 0; i < s_num; i++) {
        if(scma_[i] == -1) {
            score[i] = 1.0;
            continue;
        }
        int ct = 0;
        int mai = scma_[i];
        easy3d::vec3 d1 = (s_[i] - mavv_[mai]).normalize();
        float val = 0.0;
        for(size_t j = 0; j < sg_[i].size(); j++) {
            int idx = sg_[i][j];
            int maj = scma_[idx];
            easy3d::vec3 d2 = (s_[idx] - mavv_[maj]).normalize();
            if(maj != -1
            && dihedral_angle(sn_[i], sn_[idx]) < angle_thr_
            && easy3d::dot(d1, d2) > 0) {
                double vr;
                if(mai != maj)
                    vr = fabs(mavvr_[mai] - mavvr_[maj]);
                else {
                    vr = fabs((s_[i] - mavv_[mai]).norm() - (s_[idx] - mavv_[mai]).norm());
                }
                double vd = (s_[i] - s_[idx]).norm();

                val += vr / vd;
                ct++;
            }
        }
        if(ct > 0)
            score[i] = val / ct;
        else
            score[i] = 1.0;
        score[i] = min(score[i], 1.0f);
    }

    for(int iter = 0; iter < num_smooth_iter_; iter++) {
        std::vector<float> tmp = score;
        for(int i = 0; i < s_num; i++) {
            if(scma_[i] == -1) continue;
            double val = tmp[i];
            // val = 0.0;
            int ct = 0;
            for(size_t j = 0; j < sg_[i].size(); j++) {
                int idx = sg_[i][j];
                if(dihedral_angle(sn_[i], sn_[idx]) < angle_thr_) {
                    val += tmp[idx]; ct++;
                }
            }
            score[i] = val / (ct + 1);
        }
    }
    // save_mesh_field(mesh_, score, "../out/score.ply");


    auto gcp = mesh_->face_property<int>("f:gcp_labels");


    int nb_face = s_num;
    GCoptimizationGeneralGraph gc(nb_face, 2);

    std::vector<double> data_cost(2 * nb_face);
    for(int i = 0; i < s_num; i++) {
        data_cost[i] = score[i];
        data_cost[i + nb_face] = (1.0 - score[i]);
    }
    int num_edges = sge_.size();
    std::vector<double> edge_weights;
    double sum = 0.0;
    for(int i = 0; i < num_edges; i++) {
        int x = sge_[i].first, y = sge_[i].second;
        gc.setNeighbors(x, y);
        sum += sgel_[i];
        if(dihedral_angle(sn_[x], sn_[y]) < angle_thr_)
            edge_weights.emplace_back(lamdba_ * sgel_[i]);
        else
            edge_weights.emplace_back(0);
    }
    // for(int i = 0; i < num_edges; i++) {
    //     edge_weights[i] /= sum;
    // }
    GCP::DataCost data_item(data_cost, nb_face, 2);
    GCP::SmoothCost smooth_item(sge_, edge_weights);
    gc.setDataCostFunctor(&data_item);
    gc.setSmoothCostFunctor(&smooth_item);
    std::cout << "Before optimization energy is " << gc.compute_energy() << std::endl;
    gc.expansion(10);
    std::cout << "After optimization energy is " << gc.compute_energy() << std::endl;
#pragma omp parallel for
    for(int i = 0; i < nb_face; i++) {
        easy3d::SurfaceMesh::Face f(i);
        if(gc.whatLabel(i) == 0) {
            gcp[f] = 0;
        } else {
            gcp[f] = 1;
        }
    }
    fillet_labels = gcp.vector();
    gcp_time_ = sw.elapsed_seconds(5);
    std::cout << "Graph cut perfrom sussessfully! time="<< gcp_time_ <<std::endl;
}




void FilletSegV8::crop_local_patch(int sid
                                 , float max_gap
                                 , std::vector<int>& indices) {
    std::set<int> vis;
    std::priority_queue<std::pair<float, int>> que;
    que.push(std::make_pair(0.0, sid));
    while(!que.empty()) {
        auto cur_dis = que.top().first;
        auto cur_sid = que.top().second;
        que.pop();
        int num = sg_[cur_sid].size();
        for(int i = 0; i < num; i++) {
            int nxt_sid = sg_[cur_sid][i];
            float val = (s_[cur_sid] - s_[nxt_sid]).norm();
            if(max_gap > -cur_dis + val && vis.find(nxt_sid) == vis.end()
                && dihedral_angle( sn_[cur_sid], sn_[nxt_sid]) < angle_thr_) {
                que.push(std::make_pair(cur_dis - val, nxt_sid));
                vis.insert(nxt_sid);
                }
        }
    }
    indices.insert(indices.end(), vis.begin(), vis.end());
}


bool FilletSegV8::search_vvcs(easy3d::vec3& vv,  float r
                         , std::vector<int>& vvns, std::vector<int>& vvcs) {
    std::set<int> vis(vvns.begin(), vvns.end());
    std::queue<int> que;
    que.push(vvns[0]);vis.erase(vvns[0]);
    while(!que.empty()) {
        int u = que.front(); que.pop();
        auto f = easy3d::SurfaceMesh::Face(u);
        for(auto v : mesh_->vertices(f)) {
            for(auto t : mesh_->faces(v)) {
                if(vis.find(t.idx()) != vis.end()) {
                    que.push(t.idx()); vis.erase(t.idx());
                }
            }
        }
    }
    if(vis.empty()) {
        return false;
    }

    while(!que.empty()) que.pop();
    vis.clear();
    vis.insert(vvns[0]);
    que.push(vvns[0]);
    // double tol = r * eps_;
    double tol = eps_;
    double err = 0;
    while(!que.empty()) {
        auto cur_sid = que.front();
        que.pop();
        for(int i = 0; i < sg_[cur_sid].size(); i++) {
            int nxt_sid = sg_[cur_sid][i];
            if(dihedral_angle(sn_[cur_sid], sn_[nxt_sid]) < angle_thr_) {
                double len = (s_[nxt_sid] - vv).norm();
                if(vis.find(nxt_sid) == vis.end() && len < r + tol) {
                    que.push(nxt_sid);
                    vis.insert(nxt_sid);
                }
            }
        }
    }
    for(int i = 0; i < 4; i++) {
        if(vis.find(vvns[i]) == vis.end()) {
            return false;
        }
    }
    // if(vis.size() <= 20) {
    //     return false;
    // }
    for(auto i : vis) {
        double len = fabs((s_[i] - vv).norm() - r);
        if(len > tol) {
            return false;
        }
    }
    vvcs = std::vector<int>(vis.begin(), vis.end());
    return true;
}

void FilletSegV8::farthest_point_sampling(int num_samples
                                        , std::vector<easy3d::vec3>& s
                                        , std::vector<std::vector<int>>& sg
                                        , std::vector<int>& indices) {
    int num = s.size();
    std::vector<float> dis(num, std::numeric_limits<float>::max());
    indices.clear();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> distribution(0, num - 1); // 生成1到100之间的随机整数

    std::cout << "num_samples=" << num_samples << std::endl;
    for(int i = 0; i < num_samples; i++) {
        // std::cout << i << std::endl;
        std::priority_queue<std::pair<float, int>> que;
        int sid;
        if(i == 0) {
            sid = distribution(gen);
        } else {
            float maxx = 0.0;
            int idx = 0;
            for(int j = 0; j < num; j++) {
                if(dis[j] > maxx) {
                    maxx = dis[j]; idx = j;
                }
            }
            sid = idx;
        }
        indices.emplace_back(sid);
        dis[sid] = 0.0;
        que.push(std::make_pair(0.0, sid));
        std::set<int> vis;
        while(!que.empty()) {
            auto cur_dis = que.top().first;
            auto cur_sid = que.top().second;
            que.pop();
            if(dis[cur_sid] < -cur_dis) {
                continue;
            }

            dis[cur_sid] = -cur_dis;

            for(size_t j = 0; j < sg[cur_sid].size(); j++) {
                int nxt_sid = sg[cur_sid][j];
                if(vis.find(nxt_sid) == vis.end()) {
                    float val = (s[cur_sid] - s[nxt_sid]).norm();
                    if(dis[nxt_sid] > -cur_dis + val) {
                        que.push(std::make_pair(cur_dis - val, nxt_sid));
                        vis.insert(nxt_sid);
                    }
                }
            }
        }
    }
}



void FilletSegV8::voronoi3d(const std::vector<easy3d::vec3>& s
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

void FilletSegV8::sor(const std::vector<easy3d::vec3>& points
                    , std::vector<bool>& labels
                    , int nb_neighbors
                    , int num_sor_iter
                    , float std_ratio) {
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


double FilletSegV8::dihedral_angle(const easy3d::vec3& n1
                                 , const easy3d::vec3& n2
                                 , bool rad) {

    // double radians = abs(acos(dot(n1, n2)));
    double radians = easy3d::geom::angle(n1, n2);
    if(rad) {
        return radians;
    } else {
        double degrees = radians * 180.0 / M_PI;
        return degrees;
    }
}


// std::vector<std::vector<int>> FilletSegV8::tarjan(const std::vector<std::vector<int>> &G) {
//     std::vector<int>dfn
// }
