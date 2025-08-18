//
// Created by xiaowuga on 2025/8/16.
//

#include <random>

#include <fillet_detector.h>
#include <voronoi3d.h>
#include <utils.h>

#include <easy3d/util/stop_watch.h>
#include <easy3d/algo/surface_mesh_geometry.h>

namespace DeFillet {

    FilletDetector::FilletDetector(SurfaceMesh* mesh, FilletDetectorParameters& parameters) {

        mesh_ = new SurfaceMesh(*mesh);

        box_ = mesh_->bounding_box();

        parameters_ = parameters;

        // Precompute and cache: dihedral angles for all edges.
        // The property key "e:dihedral-angle" stores the dihedral angle per edge.
        dihedral_angle_  =  mesh_->edge_property<float>("e:dihedral-angle");

        // Precompute and cache: centroid distance between the two adjacent triangles of each edge.
        // The property key "e:centroid_distance" stores the distance between face centroids.
        centroid_distance_ = mesh_->edge_property<float>("e:centroid_distance");

        // Precompute and cache: face normals for all faces.
        face_normals_ = mesh_->face_property<easy3d::vec3>("f:normal");

        for(auto face : mesh_->faces()) {
            face_normals_[face] = mesh_->compute_face_normal(face);
        }

        // Iterate through all edges to compute dihedral angles.
        for(auto edge : mesh_->edges()){
            dihedral_angle_[edge] = -1;     // Default/invalid marker.
            centroid_distance_[edge] = -1;  // Default/invalid marker.

            // Get the two faces adjacent to this edge.
            // For a boundary edge, one of these may be invalid.
            auto f0 = mesh_->face(edge, 0);
            auto f1 = mesh_->face(edge, 1);


            // Only compute if both adjacent faces are valid (non-boundary edge).
            if(f0.is_valid() && f1.is_valid()) {

                vec3 n0 = face_normals_[f0];
                vec3 n1 = face_normals_[f1];

                // Compute the angle between the two face normals
                // and store it as the dihedral angle for this edge.
                dihedral_angle_[edge] = angle_between(n0, n1);

                vec3 c0 = easy3d::geom::centroid(mesh_, f0);
                vec3 c1 = easy3d::geom::centroid(mesh_, f1);

                // Compute the centroid distance between the two face normals
                // and store it as centroid distance for this edge
                centroid_distance_[edge] = (c0 - c1).norm();
            }
        }
    }

    void FilletDetector::generate_voronoi_vertices() {

        std::cout << "Start generate Voronoi vertices..." <<std::endl;

        easy3d::StopWatch sw;

        std::map<SurfaceMesh::Face, std::map<SurfaceMesh::Face, std::map<SurfaceMesh::Face,std::map<SurfaceMesh::Face,int>>>> unique;

        int num_patches = parameters_.num_patches;
        std::vector<SurfaceMesh::Face> patch_centroids;

        float max_gap = farthest_point_sampling(num_patches, patch_centroids);

        max_gap = parameters_.radius_thr * box_.diagonal_length();
        int idx = 0;

        omp_lock_t lock;
        omp_init_lock(&lock);

#pragma omp parallel for
        for(int i = 0; i < num_patches; i++) {

            std::vector<SurfaceMesh::Face> patch;

            crop_local_patch(patch_centroids[i], 2 * max_gap, patch);

            std::vector<vec3> ls = face_centroids(patch); // local sites;

            std::vector<std::vector<int>> lsnv;   // local sites neighboring vertices;
            std::vector<float> lvvr; // vor_vertices_radius
            std::vector<std::vector<int>> lvns ; // local vertices neighboring sites;
            std::vector<easy3d::vec3> lvv; // local voronoi vertices;
            voronoi3d(ls, box_, lvv, lvvr,lsnv, lvns);

            for(size_t j = 0; j < lvv.size(); j++) {

                std::vector<SurfaceMesh::Face> indices(4);
                for(int k = 0; k < 4; k++) {
                    indices[k] = patch[lvns[j][k]];
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
                if(flag && box_.contains(lvv[j])) {
                    omp_set_lock(&lock);

                    VoronoiVertices vv;
                    vv.pos = lvv[j];
                    vv.radius = lvvr[j];
                    vv.neigh_sites = indices;
                    voronoi_vertices_.emplace_back(vv);

                    unique[indices[0]][indices[1]][indices[2]][indices[3]] = idx++;
                    omp_unset_lock(&lock);
                }
            }
        }


        std::cout << "Generate voronoi vertices successfully! Voronoi vertices number = " <<voronoi_vertices_.size() <<std::endl;
        std::cout << "The consumption time of Voronoi vertices generatioin: " << sw.elapsed_seconds(5) << std::endl;

    }


    void FilletDetector::filter_voronoi_vertices() {

        std::cout << "Start filter Voronoi_vertices..." <<std::endl;

        easy3d::StopWatch sw;

        std::vector<VoronoiVertices> tmp;

        int num_radius = 0;

        //Rule_1: radius threshold
        float radius_thr = parameters_.radius_thr * box_.diagonal_length();

        for(size_t i = 0; i < voronoi_vertices_.size(); i++) {

            float radius = voronoi_vertices_[i].radius;

            if(radius < radius_thr) {
                tmp.emplace_back(voronoi_vertices_[i]);
            }
            else
                ++num_radius;

        }

        voronoi_vertices_ = tmp;

        tmp.clear();

        //Rule_2 & Rule_3: osculation condition----connected, non-tangented, genus = 1
        float eps = parameters_.epsilon;
        int num_connected = 0, num_tangented = 0, num_genus = 0;
        for(size_t i = 0; i < voronoi_vertices_.size(); i++) {

            float thickness = voronoi_vertices_[i].radius * eps;
            int state_code = check_osculation_condition(voronoi_vertices_[i], thickness);

            if(state_code == 0) {
                tmp.emplace_back(voronoi_vertices_[i]);
            }
            else if(state_code == 1) {
                ++num_connected;
            }
            else if(state_code == 2) {
                ++num_tangented;
            }
            else if(state_code == 3) {
                ++num_genus;
            }
        }

        voronoi_vertices_ = tmp;

        std::cout << "Filtered out by radius threshold: " << num_radius <<std::endl;
        std::cout << "Filtered out by connected condition: " << num_connected <<std::endl;
        std::cout << "Filtered out by tangented condition: " << num_tangented <<std::endl;
        std::cout << "Filtered out by genus condition: " << num_genus <<std::endl;
        std::cout << "Remaining Voronoi vertices: " << tmp.size() <<std::endl;
        std::cout << "The consumption time of filtering Voronoi vertices: " << sw.elapsed_seconds(5) << std::endl;

    }


    void FilletDetector::density_driven_voronoi_drift() {


        std::vector<easy3d::vec3>& vv
                , std::vector<float>& vvr
                , std::vector<std::vector<int>>&vvns
                , std::vector<std::vector<int>>&vvcs
                , std::vector<std::vector<int>>& scvv
                , std::vector<float>& vvs

        std::cout << "Start voronoi density drift..." <<std::endl;
        easy3d::StopWatch ss;
        std::vector<easy3d::vec3> ovv = vv;// original voronoi vertices
        std::vector<float> ovvr = vvr; // original voronoi vertiecs radius
        std::vector<std::vector<int>> ovvns = vvns; // original voronoi vertices neighboring sites
        std::vector<std::vector<int>> ovvcs = vvcs; // voronoi vertices corresponding sites
        std::vector<float> ovvs = vvs; // original voronoi vertiecs scores

        for(int iter = 1; iter <= 10; iter++) {
            std::cout << "iteration " << iter << " is running..."<< std::endl;
            easy3d::StopWatch sw;

            int num = vv.size();
            std::vector<easy3d::vec3> tmp_vv = vv;
            std::vector<float> tmp_vvr = vvr;
            omp_lock_t lock;
            omp_init_lock(&lock);

#pragma omp parallel for
            for(int i = 0; i < num; i++) {
                std::set<int> unique;
                for(size_t j = 0; j < vvns[i].size(); j++) {
                    int idx = vvns[i][j];
                    unique.insert(scvv[idx].begin(), scvv[idx].end());
                }
                std::vector<easy3d::vec3> points;
                std::vector<int> indices;
                for(auto id : unique) {
                    points.emplace_back(tmp_vv[id]);
                    indices.emplace_back(id);
                }
                MeanShift ms;

                std::pair<easy3d::vec3, int> res = ms.run(vv[i],points, h_, mean_shift_eps_);
                int idx = res.second;
                idx = indices[idx];
                vv[i] = res.first;
                vvr[i] = tmp_vvr[idx];
            }

            float t = sw.elapsed_seconds(5);
            std::cout << "done. time=" << t <<std::endl;

            std::string path = "../out/vv_" + std::to_string(iter) +".ply";
            save_point_set(vv_, path);
        }

        std::cout  << "reprojection..." << std::endl;
        std::vector<KNN::Point> data;
        int num = ovv.size();
        for(int i = 0; i < num; i++) {
            data.emplace_back(KNN::Point(ovv[i].x, ovv[i].y, ovv[i].z));
        }
        KNN::KdSearch kds(data);
#pragma omp parallel for
        for(int i = 0; i < num; i++) {
            std::vector<size_t> tmp_indices;
            std::vector<double> dist;
            kds.kth_search(KNN::Point(vv[i].x, vv[i].y, vv[i].z), 1, tmp_indices, dist);
            vv[i] = ovv[tmp_indices[0]];
            vvr[i] = ovvr[tmp_indices[0]];
            vvns[i] = ovvns[tmp_indices[0]];
            vvs[i] = ovvs[tmp_indices[0]];
        }

        // std::cout << "DAS" <<std::endl;
        num = scvv.size();
        scvv.clear(); scvv.resize(num);
        omp_lock_t lock;
        omp_init_lock(&lock);
        num = vv.size();
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
        vdd_time_ = ss.elapsed_seconds(5);
        std::cout << "Voronoi density drift sussessfully! time="<<vdd_time_<<std::endl;
        std::string path = "../out/vv_rep.ply";
    }


    PointCloud* FilletDetector::voronoi_vertices() const {

        PointCloud* cloud = new PointCloud;

        int num = voronoi_vertices_.size();

        for(int i = 0; i < num; i++) {
            cloud->add_vertex(voronoi_vertices_[i].pos);
        }

        return cloud;
    }


    float FilletDetector::farthest_point_sampling(int num_samples,
                                                  std::vector<SurfaceMesh::Face>& indices) {

        int num = mesh_->n_faces();

        indices.clear();
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int> distribution(0, num - 1);

        std::vector<float> dis(num, std::numeric_limits<float>::max());

        for(int i = 0; i < num_samples; i++) {

            std::priority_queue<std::pair<float, SurfaceMesh::Face>> que;
            SurfaceMesh::Face sface;
            if(i == 0) {
                sface = SurfaceMesh::Face(distribution(gen));
            } else {
                float maxx = 0.0;
                int idx = 0;
                for(int j = 0; j < num; j++) {
                    if(dis[j] > maxx) {
                        maxx = dis[j]; idx = j;
                    }
                }
                sface = SurfaceMesh::Face(idx);
            }

            indices.emplace_back(sface);

            dis[sface.idx()] = 0.0;
            que.push(std::make_pair(0.0, sface));

            std::set<SurfaceMesh::Face> vis;
            while(!que.empty()) {

                auto cur_dis = que.top().first;
                auto cur_face = que.top().second;
                que.pop();

                int cur_fid = cur_face.idx();

                if(dis[cur_fid] < -cur_dis) {
                    continue;
                }

                dis[cur_fid] = -cur_dis;


                for(auto halfdege : mesh_->halfedges(cur_face)) {
                    auto opp_face = mesh_->face(mesh_->opposite(halfdege));
                    auto edge = mesh_->edge(halfdege);

                    if(vis.find(opp_face) == vis.end()) {
                        float val = centroid_distance_[edge];
                        if(dis[opp_face.idx()] > -cur_dis + val) {
                            que.push(std::make_pair(cur_dis - val, opp_face));
                            vis.insert(opp_face);
                        }
                    }
                }

            }
        }

        float max_gap = 0.0;
        for(int i = 0; i < num; i++) {
            if(dis[i] > max_gap) {
                max_gap = dis[i];
            }
        }

        return max_gap;

    }

    void FilletDetector::crop_local_patch(SurfaceMesh::Face cface,
                                          float max_gap, std::vector<SurfaceMesh::Face>& patch) {

        std::set<SurfaceMesh::Face> vis;
        std::priority_queue<std::pair<float, SurfaceMesh::Face>> que;
        que.push(std::make_pair(0.0, cface));

        float angle_thr = parameters_.angle_thr;

        while(!que.empty()) {

            auto cur_dis = que.top().first;
            auto cur_face = que.top().second;
            que.pop();

            for(auto halfdege : mesh_->halfedges(cur_face)) {
                auto opp_face = mesh_->face(mesh_->opposite(halfdege));

                if(vis.find(opp_face) == vis.end()) {
                    auto edge = mesh_->edge(halfdege);
                    float val = centroid_distance_[edge];
                    if(max_gap > -cur_dis + val && dihedral_angle_[edge] < angle_thr) {
                        que.push(std::make_pair(cur_dis - val, opp_face));
                        vis.insert(opp_face);
                    }
                }
            }

        }

        patch.insert(patch.end(), vis.begin(), vis.end());
    }

    int FilletDetector::check_osculation_condition(VoronoiVertices& vor, float thickness) const {

        std::vector<SurfaceMesh::Face>& neigh_sites = vor.neigh_sites;
        std::vector<SurfaceMesh::Face>& corr_sites = vor.corr_sites;
        vec3 pos = vor.pos;
        float radius = vor.radius;
        float thr = thickness * 0.5;

        std::set<SurfaceMesh::Face> vis;
        std::queue<SurfaceMesh::Face> que;

        que.push(neigh_sites[0]);

        while(!que.empty()) {

            SurfaceMesh::Face cur_face = que.front(); que.pop();

            if(vis.find(cur_face) != vis.end())
                continue;

            vis.insert(cur_face);

            for(auto halfdege : mesh_->halfedges(cur_face)) {
                auto opp_face = mesh_->face(mesh_->opposite(halfdege));

                if(opp_face.is_valid() && vis.find(opp_face) == vis.end()) {
                    vec3 opp_pos = easy3d::geom::centroid(mesh_, opp_face);
                    float err = fabs((pos - opp_pos).norm() - radius);
                    if(err < thr) {
                        que.push(opp_face);
                    }
                }
            }
        }

        // Check if corr_stites are connected
        for(size_t i = 0; i < neigh_sites.size(); i++) {
            if(vis.find(neigh_sites[i]) == vis.end()) {
                return 1;
            }
        }

        // Check if corr_stites are tangented to surface
        std::vector<SurfaceMesh::Edge> boundaries;

        for(size_t i = 0; i < neigh_sites.size(); i++) {
            for(auto halfedge : mesh_->halfedges(neigh_sites[i])) {
                SurfaceMesh::Edge edge  = mesh_->edge(halfedge);
                boundaries.emplace_back(edge);
            }
        }

        if(count_boundaries_components(boundaries) == 1) {
            return 2;
        }


        boundaries.clear();
        //Check if corr_stites's genus = 1
        std::vector<SurfaceMesh::Face> tmp(vis.begin(), vis.end());
        for(int i = 0; i < tmp.size(); i++) {
            SurfaceMesh::Face cur = tmp[i];
            for(auto halfedge : mesh_->halfedges(cur)) {

                auto opp_face = mesh_->face(mesh_->opposite(halfedge));

                if(opp_face.is_valid() && vis.find(opp_face) == vis.end()) {

                    SurfaceMesh::Edge edge  = mesh_->edge(halfedge);
                    boundaries.emplace_back(edge);

                }
            }
        }

        if(count_boundaries_components(boundaries) > 1) {
            return 3;
        }


        corr_sites = tmp;

        return 0;

    }


    int FilletDetector::count_boundaries_components(std::vector<SurfaceMesh::Edge>& boundaries) const {

        std::map<SurfaceMesh::Vertex, int> mp;
        int id = 0;

        for(size_t i = 0; i < boundaries.size(); i++) {
            SurfaceMesh::Edge edge  = boundaries[i];

            SurfaceMesh::Vertex v0 = mesh_->vertex(edge, 0);
            SurfaceMesh::Vertex v1 = mesh_->vertex(edge, 1);

            if(mp.find(v0) == mp.end()) {
                mp[v0] = id++;
            }

            if(mp.find(v1) == mp.end()) {
                mp[v1] = id++;
            }
        }

        std::vector<int> bcj(id);

        function<int(int)> find = [&](int x){ return bcj[x]==x ? x : bcj[x]=find(bcj[x]); };

        std::iota(bcj.begin(), bcj.end(), 0);

        for(size_t i = 0; i < boundaries.size(); i++) {
            SurfaceMesh::Edge edge  = boundaries[i];

            SurfaceMesh::Vertex v0 = mesh_->vertex(edge, 0);
            SurfaceMesh::Vertex v1 = mesh_->vertex(edge, 1);

            int id0 = mp[v0], id1 = mp[v1];

            bcj[find(id0)] = find(id1);
        }

        int num_components = 0;
        for(size_t i = 0; i < bcj.size(); i++) {
            if(bcj[i] == i)
                 num_components++;
        }

        return num_components;
    }


    std::vector<vec3> FilletDetector::face_centroids(std::vector<SurfaceMesh::Face>& faces) const{

        std::vector<vec3> centroids;

        for(size_t i = 0; i < faces.size(); i++) {
            vec3 p = easy3d::geom::centroid(mesh_, faces[i]);
            centroids.emplace_back(p);
        }

        return centroids;
    }

}
