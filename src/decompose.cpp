//
// Created by xiaowuga on 2025/8/14.
//

#include "decompose.h"


namespace DeFillet {

    /**
     * @brief Compute the (unsigned) angle between two 3D vectors in degrees.
     *
     * Uses atan2(||n1 × n2||, n1 · n2), which is generally more numerically stable
     * than acos of the normalized dot product.
     *
     * @param n1 First vector (e.g., a face normal).
     * @param n2 Second vector (e.g., a face normal).
     * @return Angle in degrees, in the range [0, 180].
     */
    double angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2) {
        const double dot = easy3d::dot(n1, n2);
        const double cross_norm = easy3d::cross(n1, n2).norm();
        return std::atan2(cross_norm, dot) * 180.0 / M_PI;
    }

    Decompose::Decompose(SurfaceMesh* mesh, float angle_thr) : mesh_(mesh), angle_thr_(angle_thr) {

        // Precompute (and cache) dihedral angles for all edges.
        dihedral_angle_  =  mesh_->edge_property<float>("e:dihedral-angle");

        mesh_->update_face_normals();
        face_normals_ = mesh_->face_property<easy3d::vec3>("f:nromal");

        for(auto edge : mesh_->edges()){
            dihedral_angle_[edge] = -1;
            auto f0 = mesh_->face(edge, 0);
            auto f1 = mesh_->face(edge, 1);
            if(f0.is_valid() && f1.is_valid()) {
                easy3d::vec3 n0 = face_normals_[f0];
                easy3d::vec3 n1 = face_normals_[f1];
                dihedral_angle_[edge] = angle_between(n0, n1);
            }
        }

        component_labels_ = mesh_->face_property<int>("f:component_label");

        for(auto face : mesh_->faces()) {
            component_labels_[face] = -1;
        }
    }

    Decompose::~Decompose() {

        // auto dihedral_angle  = mesh_->get_edge_property<float>("e:dihedral-angle");
        if(dihedral_angle_)
            mesh_->remove_edge_property<float>(dihedral_angle_);

        // auto face_normals = mesh_->get_face_property<easy3d::vec3>("f:nromal");
        if(face_normals_)
            mesh_->remove_face_property(face_normals_);

        // auto component_labels = mesh_->get_face_property<int>("f:component_label");
        if(component_labels_)
            mesh_->remove_face_property(component_labels_);

        clear_all_components();

    }


    void Decompose::apply() {
        // auto component_labels = mesh_->add_face_property<int>("f:component_label", -1);
        // auto dihedral_angle  = mesh_->get_edge_property<float>("e:dihedral-angle");

        int num_component = 0;

        // decompose via dihedral_angle
        for(auto face : mesh_->faces()) {
            if(component_labels_[face] == -1) {
                std::cout << "ASD" << std::endl;
                std::queue<easy3d::SurfaceMesh::Face> que;
                que.push(face);
                while(!que.empty()) {
                    auto cur = que.front(); que.pop();
                    if(component_labels_[cur] != -1)
                        continue;
                    component_labels_[cur] = num_component;
                    for(auto halfdege : mesh_->halfedges(cur)) {
                        auto opp_face = mesh_->face(mesh_->opposite(halfdege));
                        auto edge = mesh_->edge(halfdege);
                        if(opp_face.is_valid() && component_labels_[opp_face] == -1
                            && dihedral_angle_[edge] < angle_thr_) {
                                que.push(opp_face);
                        }
                    }
                }
                num_component++;
            }
        }

        std::cout << num_component << std::endl;
        //clear old components
        clear_all_components();

        std::queue<std::pair<SurfaceMesh*,int>> que;

        for(int i = 0; i < num_component; i++) {
            auto component = split_component(i);
            que.push( std::make_pair(component, i));
        }

        while(!que.empty()) {

            auto component = que.front().first;
            int component_id = que.front().second;
            que.pop();

            bool flag = true;
            for(auto edge : component->edges()) {

                auto f0 = component->face(edge, 0);
                auto f1 = component->face(edge, 1);

                if(f0.is_valid() && f1.is_valid()) {

                    vec3 n0 = component->compute_face_normal(f0);
                    vec3 n1 = component->compute_face_normal(f1);

                    float angle = angle_between(n0, n1);

                    if(angle > angle_thr_) {

                        int new_id0 = num_component++;
                        int new_id1 = num_component++;

                        auto res = split_alone_egde(component, edge, component_id, new_id0, new_id1);

                        que.push(std::make_pair(res[0], new_id0));
                        que.push(std::make_pair(res[1], new_id1));

                        flag = false;
                        break;
                    }
                }
            }

            if(flag) {
                components_.emplace_back(component);
            }
            else {
                delete component;
            }
        }

    }

    SurfaceMesh* Decompose::split_component(int label) {
        // Collect all faces with the given label, and all vertices used by them
        std::set<easy3d::SurfaceMesh::Face> faces_set;
        std::set<easy3d::SurfaceMesh::Vertex> points_set;


        // Select faces with matching label and collect their vertices
        for(auto f : mesh_->faces()) {
            if(component_labels_[f] == label) {
                faces_set.insert(f);
                for(auto v : mesh_->vertices(f)) {
                    points_set.insert(v);
                }
            }
        }


        // Map original vertex index -> new vertex index, and inverse map
        std::map<int, size_t> mp;
        int nb_points = points_set.size();
        std::vector<int> point_map(nb_points); // new vertex id -> original vertex id

        size_t id = 0;
        SurfaceMesh* component = new SurfaceMesh;

        // Insert vertices into the new mesh
        for(auto v : points_set) {
            point_map[id] = v.idx(); // new -> original
            mp[v.idx()] = id;      // original -> new
            ++id;

            // Copy vertex position
            easy3d::vec3 p = mesh_->position(v);
            component->add_vertex(p);
        }


        // Map original face index -> new face index, and insert faces
        int nb_faces = faces_set.size();
        std::vector<int> face_map(nb_faces); // new face id -> original face id
        id = 0;
        for(auto f : faces_set) {
            face_map[id] = f.idx(); // new -> original
            ++id;

            // Remap each vertex in the face to the new vertex index
            std::vector<easy3d::SurfaceMesh::Vertex> tmp;
            for(auto v : mesh_->vertices(f)) {
                tmp.emplace_back(easy3d::SurfaceMesh::Vertex(mp[v.idx()]));
            }

            // Add the face to the new mesh (assumes triangles)
            component->add_triangle(tmp[0], tmp[1], tmp[2]);
        }

        // Create properties in the new mesh to store original indices
        auto original_point_index = component->vertex_property<int>("v:original_index");
        auto original_face_index = component->face_property<int>("f:original_index");

        // Fill vertex original indices
        for(auto v : component->vertices()) {
            original_point_index[v] = point_map[v.idx()];
        }

        // Fill face original indices
        for(auto f : component->faces()) {
            original_face_index[f] = face_map[f.idx()];
        }

        return component;
    }


    std::vector<SurfaceMesh*> Decompose::split_alone_egde(SurfaceMesh* component, SurfaceMesh::Edge e,
                                                          int component_id, int new_id0, int new_id1) {

        auto original_face_index = component->face_property<int>("f:original_index");


        SurfaceMesh::Face f0 = component->face(e, 0);
        f0 = SurfaceMesh::Face(original_face_index[f0]);

        SurfaceMesh::Face f1 = component->face(e, 1);
        f1 = SurfaceMesh::Face(original_face_index[f1]);

        std::queue<pair<SurfaceMesh::Face, int>> que;
        que.push(make_pair(f0, new_id0));
        que.push(make_pair(f1, new_id1));

        while(!que.empty()) {

            auto cur = que.front().first;
            int id = que.front().second;
            que.pop();

            if(component_labels_[cur] == id)
                continue;

            component_labels_[cur] = id;

            for(auto halfdege : mesh_->halfedges(cur)) {

                auto opp_face = mesh_->face(mesh_->opposite(halfdege));
                auto edge = mesh_->edge(halfdege);

                if(opp_face.is_valid() && component_labels_[opp_face] == component_id

                    && dihedral_angle_[edge] < angle_thr_) {
                    que.push(make_pair(opp_face, id));

                }

            }

        }

        easy3d::SurfaceMesh* res0 = split_component(new_id0);
        easy3d::SurfaceMesh* res1 = split_component(new_id1);

        return vector<SurfaceMesh*>{res0, res1};

    }


    void Decompose::clear_all_components() {

        for(auto item : components_) {
            delete item;
        }

        components_.clear();
    }

    double Decompose::angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2) {

        const double dot = easy3d::dot(n1, n2);
        const double cross_norm = easy3d::cross(n1, n2).norm();

        return std::atan2(cross_norm, dot) * 180.0 / M_PI;
    }



}

