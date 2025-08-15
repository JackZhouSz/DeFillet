//
// Created by xiaowuga on 2025/8/15.
//

#include <preliminary_fillet_selector.h>
#include <voronoi3d.h>
#include <utils.h>
#include <voronoi3d.h>

#include <easy3d/algo/surface_mesh_geometry.h>


namespace DeFillet {

    PreliminaryFilletSelector::PreliminaryFilletSelector(SurfaceMesh* mesh,int num_patches,
                                                         float radius_thr, float angle_thr) : mesh_(mesh), num_patches_(num_patches), angle_thr_(angle_thr) {

        // Compute and store the axis-aligned bounding box (AABB) of the mesh.
        box_ = mesh_->bounding_box();

        radius_thr_ = box_.diagonal_length() * radius_thr;

        // Precompute and cache: dihedral angles for all edges.
        // The property key "e:dihedral-angle" stores the dihedral angle per edge.
        dihedral_angle_  =  mesh_->edge_property<float>("e:dihedral-angle");

        // Precompute and cache: centroid distance between the two adjacent triangles of each edge.
        // The property key "e:centroid_distance" stores the distance between face centroids.
        centroid_distance_ = mesh_->edge_property<float>("e:centroid_distance");

        // The property key "f:fillet_label" stores store a label for each face.
        // Label semantics for each face:  0 = non-fillet, 1 = fillet
        fillet_label_ = mesh_->face_property<int>("f:fillet_label");

        // Precompute and cache: face normals for all faces.
        face_normals_ = mesh_->face_property<easy3d::vec3>("f:normal");

        for(auto face : mesh_->faces()) {
            face_normals_[face] = mesh_->compute_face_normal(face);
            fillet_label_[face] = 0;
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

    PreliminaryFilletSelector::~PreliminaryFilletSelector() {

        if(dihedral_angle_)
            mesh_->remove_edge_property<float>(dihedral_angle_);


        if(face_normals_)
            mesh_->remove_face_property(face_normals_);

        if(centroid_distance_)
            mesh_->remove_edge_property<float>(centroid_distance_);



    }

    void PreliminaryFilletSelector::apply() {

        float max_gap = farthest_point_sampling(num_patches_, patch_centroid_);

        omp_lock_t lock;
        omp_init_lock(&lock);
#pragma omp parallel for
        for(int i = 0; i < num_patches_; i++) {

            std::vector<SurfaceMesh::Face> patch;
            crop_local_patch(patch_centroid_[i], max_gap * 2, patch);

            auto ls = face_centroids(patch);

            std::vector<std::vector<int>> lsnv;   // local sites neighboring vertices;
            std::vector<float> lvvr; // vor_vertices_radius
            std::vector<std::vector<int>> lvns ; // local vertices neighboring sites;
            std::vector<easy3d::vec3> lvv; // local voronoi vertices;

            voronoi3d(ls, box_, lvv, lvvr,lsnv, lvns);

            for(size_t j = 0; j < lvv.size(); j++) {
                if(lvvr[j] < radius_thr_)
                    for(int k = 0; k < 4; k++) {
                        auto face = patch[lvns[j][k]];
                        omp_set_lock(&lock);
                        fillet_label_[face] = 1;
                        omp_unset_lock(&lock);
                    }
            }

            // std::cout << "ASD" << std::endl;

        }


    }

    PointCloud* PreliminaryFilletSelector::patch_centroid() const{

        PointCloud* cloud = new PointCloud;

        for(size_t i = 0; i < patch_centroid_.size(); i++) {
            auto centroid = easy3d::geom::centroid(mesh_, patch_centroid_[i]);
            cloud->add_vertex(centroid);
        }

        return cloud;
    }


    float PreliminaryFilletSelector::farthest_point_sampling(int num_samples,
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


    void PreliminaryFilletSelector::crop_local_patch(SurfaceMesh::Face cface,
                                                     float max_gap, std::vector<SurfaceMesh::Face>& patch) {

        std::set<SurfaceMesh::Face> vis;
        std::priority_queue<std::pair<float, SurfaceMesh::Face>> que;
        que.push(std::make_pair(0.0, cface));

        while(!que.empty()) {

            auto cur_dis = que.top().first;
            auto cur_face = que.top().second;
            que.pop();

            for(auto halfdege : mesh_->halfedges(cur_face)) {
                auto opp_face = mesh_->face(mesh_->opposite(halfdege));

                if(vis.find(opp_face) == vis.end()) {
                    auto edge = mesh_->edge(halfdege);
                    float val = centroid_distance_[edge];
                    if(max_gap > -cur_dis + val && dihedral_angle_[edge] < angle_thr_) {
                        que.push(std::make_pair(cur_dis - val, opp_face));
                        vis.insert(opp_face);
                    }
                }
            }

        }

        patch.insert(patch.end(), vis.begin(), vis.end());
    }

    std::vector<vec3> PreliminaryFilletSelector::face_centroids(std::vector<SurfaceMesh::Face>& faces) const{

        std::vector<vec3> centroids;

        for(size_t i = 0; i < faces.size(); i++) {
            vec3 p = easy3d::geom::centroid(mesh_, faces[i]);
            centroids.emplace_back(p);
        }

        return centroids;
    }



}
