//
// Created by xiaowuga on 2025/8/16.
//

#include <random>

#include <fillet_detector.h>
#include <voronoi3d.h>
#include <utils.h>
#include <knn4d.h>

#include <easy3d/util/stop_watch.h>
#include <easy3d/algo/surface_mesh_geometry.h>
#include <easy3d/fileio/point_cloud_io.h>


#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/hessian_energy.h>
#include <igl/curved_hessian_energy.h>
#include <igl/avg_edge_length.h>


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>





using namespace Eigen;
typedef Eigen::SparseMatrix<double> SparseMat;

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

        for(auto face : mesh_->faces()) {

            Sites site;

            site.pos = easy3d::geom::centroid(mesh_, face);
            site.face = face;
            site.radius = box_.diagonal_length();
            site.center = easy3d::vec3(0,0,0);
            site.rate = -1.0;
            site.flag = false;
            sites_.emplace_back(site);
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

                vec3 c0 = sites_[f0.idx()].pos;
                vec3 c1 = sites_[f1.idx()].pos;

                // Compute the centroid distance between the two face normals
                // and store it as centroid distance for this edge
                centroid_distance_[edge] = (c0 - c1).norm();
            }
        }

    }

    void FilletDetector::generate_voronoi_vertices() {

        std::cout << "Start generating Voronoi vertices..." <<std::endl;

        easy3d::StopWatch sw;

        voronoi_vertices_.clear();

        std::map<SurfaceMesh::Face, std::map<SurfaceMesh::Face, std::map<SurfaceMesh::Face,std::map<SurfaceMesh::Face,int>>>> unique;

        int num_patches = parameters_.num_patches;
        std::vector<SurfaceMesh::Face> patch_centroids;

        float max_gap = farthest_point_sampling(num_patches, patch_centroids);


        max_gap = 0.5 * M_PI * parameters_.radius_thr * box_.diagonal_length();
        int idx = 0;

        omp_lock_t lock;
        omp_init_lock(&lock);

#pragma omp parallel for
        for(int i = 0; i < num_patches; i++) {

            std::vector<SurfaceMesh::Face> patch;

            crop_local_patch(patch_centroids[i], max_gap, patch);



            std::vector<vec3> ls = face_centroids(patch); // local sites;

            if(ls.size() < 10)
                continue;

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


        std::cout << "Successfully generated Voronoi vertices! Voronoi vertices number = " <<voronoi_vertices_.size() <<std::endl;
        std::cout << "The consumption time of Voronoi vertices generatioin: " << sw.elapsed_seconds(5) << std::endl;

    }

    void FilletDetector::compute_voronoi_vertices_density_field() {

        std::cout << "Start computing Voronoi vertices density field..." <<std::endl;

        easy3d::StopWatch sw;

        int num = voronoi_vertices_.size();

        std::vector<KNN::Point4D> knn_data;
        std::vector<easy3d::vec4> data;

        for(int i = 0; i < num; i++) {
            auto& pos = voronoi_vertices_[i].pos;
            float radius = voronoi_vertices_[i].radius;
            knn_data.emplace_back(KNN::Point4D(pos.x,pos.y, pos.z, radius));
            data.emplace_back(easy3d::vec4(pos.x,pos.y, pos.z, radius));
        }

        KNN::KdSearch4D kds(knn_data);
        int num_neighbors = parameters_.num_neighbors;
        int sigma = parameters_.sigma;

        float max_value = 0, min_value = box_.diagonal_length();

#pragma omp parallel for
        for(int i = 0; i < num; i++) {
            std::vector<size_t> indices;
            std::vector<double> dist;
            kds.kth_search(knn_data[i], num_neighbors, indices, dist);
            float density = 0;
            double w = 0;
            std::vector<vec4> points;
            for(size_t j = 0; j < indices.size(); j++) {
                int idx = indices[j];
                double dis = fabs((data[idx] - data[i]).norm());
                double val = gaussian_kernel(dis, sigma);
                density += val * dis;
                w += val;
                points.emplace_back(data[idx]);
            }
            voronoi_vertices_[i].density =  density / w;
            max_value = max(max_value, voronoi_vertices_[i].density);
            min_value = min(min_value, voronoi_vertices_[i].density);

            voronoi_vertices_[i].axis = axis_direction(points);


        }

        std::cout << "max_value: " << max_value << std::endl;
        std::cout << "min_value: " << min_value << std::endl;

        std::cout << "The consumption time of computing Voronoi vertices density field: " << sw.elapsed_seconds(5) << std::endl;
    }

    void FilletDetector::filter_voronoi_vertices() {

        std::cout << "Start filtering Voronoi_vertices..." <<std::endl;

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
            else if(state_code == 2 || state_code == 4) {
                ++num_tangented;
            }
            else if(state_code == 3) {
                ++num_genus;
            }
        }

        std::cout << "Filtered out by radius threshold: " << num_radius <<std::endl;
        std::cout << "Filtered out by connected condition: " << num_connected <<std::endl;
        std::cout << "Filtered out by tangented condition: " << num_tangented <<std::endl;
        std::cout << "Filtered out by genus condition: " << num_genus <<std::endl;

        // voronoi_vertices_ = tmp;
        std::vector<vec4> points;
        std::vector<bool>labels;

        for(int i = 0; i < tmp.size(); i++) {
            points.emplace_back(vec4(tmp[i].pos.x,tmp[i].pos.y, tmp[i].pos.z, tmp[i].radius));
            labels.emplace_back(true);
        }
        sor(points, labels, 30,3, 0.3);
        voronoi_vertices_.clear();
        for(int i = 0; i < tmp.size(); i++) {
            if(labels[i])
                voronoi_vertices_.emplace_back(tmp[i]);
        }

        std::cout << "Remaining Voronoi vertices: " << voronoi_vertices_.size() <<std::endl;
        std::cout << "The consumption time of filtering Voronoi vertices: " << sw.elapsed_seconds(5) << std::endl;

    }


    void FilletDetector::rolling_ball_trajectory_transform() {

        std::cout << "Start rolling-ball trajectory transform..." <<std::endl;

        easy3d::StopWatch sw;


        for(size_t i = 0; i < voronoi_vertices_.size(); i++) {

            auto& corr_sites = voronoi_vertices_[i].corr_sites;

            for(int j = 0; j < corr_sites.size(); j++) {
                auto idx = corr_sites[j].idx();
                sites_[idx].corr_vv.emplace_back(i);
            }
        }

        std::vector<Sites> sites = sites_;

        for(size_t i = 0; i < sites_.size(); i++) {

            auto& corr_vv = sites_[i].corr_vv;

            if(!corr_vv.empty()) {
                sites[i].flag = true;
                float min_value = 1e9;
                for(int j = 0; j < corr_vv.size(); j++) {
                    int idx = corr_vv[j];

                    if(min_value > voronoi_vertices_[idx].density) {
                        min_value = voronoi_vertices_[idx].density;
                        sites[i].center = voronoi_vertices_[idx].pos;
                        sites[i].vvid = idx;
                        sites[i].radius = (sites_[i].pos - voronoi_vertices_[idx].pos).norm();
                    }
                }
            }
        }
        sites_ = sites;
        float angle_thr = parameters_.angle_thr;
        for(size_t i = 0; i < sites_.size(); i++) {
            auto& corr_vv = sites_[i].corr_vv;
            vec3 p0 = sites_[i].pos;
            if(corr_vv.empty()) {
                std::priority_queue<pair<int,SurfaceMesh::Face>>que;
                que.push(make_pair(0, SurfaceMesh::Face(i)));
                while(!que.empty()) {
                    int step = que.top().first;
                    SurfaceMesh::Face cur = que.top().second; que.pop();

                    if(-step > 5) continue;

                    if(!sites_[cur.idx()].corr_vv.empty()) {
                        vec3 p1 = sites_[cur.idx()].center;
                        vec3 d = voronoi_vertices_[sites_[cur.idx()].vvid].axis;
                        sites[i].center = p1 + dot((p0 - p1), d) * d;
                        sites[i].radius = (p0 - sites[i].center).norm();
                        sites[i].flag = true;
                        break;
                    }


                    for(auto halfedge : mesh_->halfedges(cur)) {
                        auto opp_face = mesh_->face(mesh_->opposite(halfedge));
                        auto edge = mesh_->edge(halfedge);

                        if(opp_face.is_valid() && dihedral_angle_[edge] < angle_thr) {
                            que.push(make_pair(step-1, opp_face));
                        }
                    }

                 }

            }
        }

        sites_ = sites;
        std::cout << "The consumption time of rolling-ball trajectorytransform: " << sw.elapsed_seconds(5) << std::endl;


        PointCloud* cloud = new PointCloud;
        for(int i = 0; i < sites_.size(); i++) {
            if(sites_[i].flag) {
                cloud->add_vertex(sites_[i].center);
                // std::cout << sites_[i].center.x << ' ' << sites_[i].center.y << ' ' << sites_[i].center.y << std::endl;
            }
        }
        easy3d::PointCloudIO::save("../out3/center.ply", cloud);


    }

    void FilletDetector::compute_fillet_radius_rate_field() {
        std::cout << "Start computing fillet radius rate field..." <<std::endl;

        easy3d::StopWatch sw;
        float angle_thr = parameters_.angle_thr;
        for(auto face : mesh_->faces()) {
            if(sites_[face.idx()].flag == false) {
                sites_[face.idx()].rate = 1.0;
                continue;
            }
            // if(sites_[face.idx()].corr_vv.empty()) {
            //     sites_[face.idx()].rate = 1.0;
            //     continue;
            // }

            float rate = 0;
            int count = 0;
            for(auto halfedge : mesh_->halfedges(face)) {

                auto opp_face = mesh_->face(mesh_->opposite(halfedge));
                auto edge = mesh_->edge(halfedge);

                if(opp_face.is_valid() && dihedral_angle_[edge] < angle_thr) {

                    float dis = centroid_distance_[edge];
                    float radius_diff = fabs(sites_[face.idx()].radius - sites_[opp_face.idx()].radius);
                    rate += radius_diff / dis;
                    ++count;
                }
            }

            if(count) {
                sites_[face.idx()].rate = min(rate / count, 1.0f);
                // sites_[face.idx()].rate = rate / count;
            }
            else {
                sites_[face.idx()].rate = 1.0;
            }
        }

        std::cout << "The consumption time of computing fillet radius rate field: " << sw.elapsed_seconds(5) << std::endl;
    }

    void FilletDetector::rate_field_smoothing() {

        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd noisy;
        int num_vertices = mesh_->n_vertices();
        V.resize(num_vertices, 3);
        noisy.resize(num_vertices);

        for(auto vertex : mesh_->vertices()) {
            int vid = vertex.idx();
            auto& pos = mesh_->position(vertex);
            V(vid , 0) = pos.x;
            V(vid , 1) = pos.y;
            V(vid , 2) = pos.x;

            noisy[vid] = 0.0;
            int ct = 0;
            for(auto face : mesh_->faces(vertex)) {
                int idx = face.idx();
                noisy[vid] += sites_[idx].rate;
                ++ct;
            }
            if(ct)
                noisy[vid] /= ct;
        }

        int num_faces = mesh_->n_faces();
        F.resize(num_faces, 3);
        for(auto face : mesh_->faces()) {

            int fid = face.idx();

            std::vector<int> indices;
            for(auto vertex : mesh_->vertices(face)) {
                indices.emplace_back(vertex.idx());
            }
            F(fid , 0) = indices[0];
            F(fid , 1) = indices[1];
            F(fid , 2) = indices[2];
        }

        SparseMat L, M;
        igl::cotmatrix(V, F, L);
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        Eigen::SimplicialLDLT<SparseMat> solver(M);
        SparseMat MinvL = solver.solve(L);
        SparseMat QL = L.transpose()*MinvL;
        SparseMat QH;
        igl::hessian_energy(V, F, QH);
        SparseMat QcH;
        igl::curved_hessian_energy(V, F, QcH);

        const double al = 0.5;
        SparseMat A = al*QL + (1.-al)*M;
        SparseMat I(A.rows(), A.cols()); I.setIdentity();
        Eigen::SparseLU<SparseMat> lapSolver(A + I);
        Eigen::VectorXd denoisy = lapSolver.solve(al*M*noisy);


        for(auto face : mesh_->faces()) {
            float val = 0;
            int ct = 0;
            for(auto vertex : mesh_->vertices(face)) {
                int idx = vertex.idx();
                val += denoisy[idx];
                ++ct;
            }
            if(ct)
                val /= ct;
            sites_[face.idx()].rate = val;
        }

    }


    PointCloud* FilletDetector::voronoi_vertices() const {

        PointCloud* cloud = new PointCloud;

        int num = voronoi_vertices_.size();

        for(int i = 0; i < num; i++) {
            cloud->add_vertex(voronoi_vertices_[i].pos);
        }

        auto nor = cloud->add_vertex_property<vec3>("v:normal");

        for(auto  v : cloud->vertices()) {
            nor[v] = voronoi_vertices_[v.idx()].axis;
        }

        return cloud;
    }

    std::vector<float> FilletDetector::radius_rate_field() const {

        std::vector<float> field;
        for(int i = 0; i < sites_.size(); i++) {
            field.emplace_back(sites_[i].rate);
        }

        return field;
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
//                sface = SurfaceMesh::Face(distribution(gen));
                for(auto face : mesh_->faces()) {
                    if(labels_[face]) {
                        sface = face; break;
                    }
                }
            } else {
                float maxx = 0.0;
                int idx = 0;
                for(int j = 0; j < num; j++) {
                    SurfaceMesh::Face f(j);
                    if(labels_[f] && dis[j] > maxx) {
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


                for(auto halfedge : mesh_->halfedges(cur_face)) {
                    auto opp_face = mesh_->face(mesh_->opposite(halfedge));
                    auto edge = mesh_->edge(halfedge);

                    if(labels_[opp_face] && vis.find(opp_face) == vis.end()) {
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
            SurfaceMesh::Face f(i);
            if(labels_[f] && dis[i] > max_gap) {
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

            for(auto halfedge : mesh_->halfedges(cur_face)) {
                auto opp_face = mesh_->face(mesh_->opposite(halfedge));

                if(vis.find(opp_face) == vis.end()) {
                    auto edge = mesh_->edge(halfedge);
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
                    vec3 opp_pos = sites_[opp_face.idx()].pos;
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

        if(count_boundaries_components(boundaries) < 3 ) {
            return 2;
        }


        boundaries.clear();
        //Check if corr_stites's genus = 1
        std::vector<SurfaceMesh::Face> tmp(vis.begin(), vis.end());
        for(int i = 0; i < tmp.size(); i++) {
            SurfaceMesh::Face cur = tmp[i];
            for(auto halfedge : mesh_->halfedges(cur)) {

                auto opp_face = mesh_->face(mesh_->opposite(halfedge));

                if(labels_[opp_face]  && opp_face.is_valid() && vis.find(opp_face) == vis.end()) {

                    SurfaceMesh::Edge edge  = mesh_->edge(halfedge);
                    boundaries.emplace_back(edge);

                }
            }
        }

        if(count_boundaries_components(boundaries) > 1) {
            return 3;
        }


//         float maxx = 0;
//         for(int i = 0; i < neigh_sites.size(); i++) {
//             for(int j = i + 1; j < neigh_sites.size(); j++) {
//                 vec3 n1 = face_normals_[neigh_sites[i]];
//                 vec3 n2 = face_normals_[neigh_sites[j]];
//                 maxx = max(maxx, angle_between(n1, n2));
//             }
//         }
//
//         if(maxx < 50) {
// //            return 4;
//         }

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
            vec3 p = sites_[faces[i].idx()].pos;
            centroids.emplace_back(p);
        }

        return centroids;
    }

    SurfaceMesh::FaceProperty<int> FilletDetector::labels() const {
        return labels_;
    }

}
