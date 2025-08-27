//
// Created by xiaowuga on 2025/8/17.
//

#include <numeric>

#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/core/random.h>


#include <utils.h>
#include <knn4d.h>

#include <Eigen/Dense>

#include <igl/jet.h>

using namespace easy3d;
using namespace std;

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
    float angle_between(const easy3d::vec3& n1, const easy3d::vec3& n2) {
        const double dot = easy3d::dot(n1, n2);
        const double cross_norm = easy3d::cross(n1, n2).norm();

        return atan2(cross_norm, dot) * 180.0 / M_PI;
    }

    double gaussian_kernel(double distance, double kernel_bandwidth){
        double temp =  exp(-0.5 * (distance*distance) / (kernel_bandwidth*kernel_bandwidth));
        return temp;
    }

    vec3 axis_direction(std::vector<easy3d::vec4> points) {

        int num = points.size();
        Eigen::MatrixXd X(num, 3);

        for(int i = 0; i < num; i++) {
            X.row(i) << points[i].x, points[i].y, points[i].z;
        }

        Eigen::RowVector3d mu = X.colwise().mean();
        Eigen::MatrixXd Y = X.rowwise() - mu;

        Eigen::Matrix3d C = (Y.transpose() * Y) / double(num);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);

        Eigen::Matrix3d evecs = es.eigenvectors();

        float x = evecs.col(2).x();
        float y = evecs.col(2).y();
        float z = evecs.col(2).z();

        vec3 res(x, y, z);
        res.normalize();

        return res;

    }


    easy3d::vec3 project_to_line(easy3d::vec3 pos, easy3d::vec3 point, easy3d::vec3 dir) {
        easy3d::vec3 vdir = pos - point;
        float t = easy3d::dot(vdir, dir) / dir.norm();
        return point + t * dir;
    }

    easy3d::SurfaceMesh* split_component(const easy3d::SurfaceMesh* mesh,
                         easy3d::SurfaceMesh::FaceProperty<int>& component_labels,
                         int label) {
        // Collect all faces with the given label, and all vertices used by them
        std::set<easy3d::SurfaceMesh::Face> faces_set;
        std::set<easy3d::SurfaceMesh::Vertex> points_set;


        // Select faces with matching label and collect their vertices
        for(auto f : mesh->faces()) {
            if(component_labels[f] == label) {
                faces_set.insert(f);
                for(auto v : mesh->vertices(f)) {
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
            easy3d::vec3 p = mesh->position(v);
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
            for(auto v : mesh->vertices(f)) {
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


    void sor(const std::vector<easy3d::vec4>& points
                    , std::vector<bool>& labels
                    , int nb_neighbors
                    , int num_sor_iter
                    , float std_ratio) {
    int nb_points = points.size();
    // int nb_neighbors = 30, num_sor_iter = 3;
    // float std_ratio = 0.3;
    // labels.resize(nb_points, fal);
    for(int i = 0; i < num_sor_iter; i++) {
        std::vector<KNN::Point4D> knn_points;
        std::vector<int> indices;
#pragma omp parallel for
        for(int j = 0; j < nb_points; j++) {
            if(labels[j]) {
                #pragma omp critical
                {
                    knn_points.emplace_back(KNN::Point4D(points[j].x,
                                                       points[j].y, points[j].z, points[j].w));
                    indices.emplace_back(j);
                }
            }
        }
        KNN::KdSearch4D kds(knn_points);
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


    void save_components(const easy3d::SurfaceMesh* mesh,
                            const std::vector<easy3d::SurfaceMesh*>components,
                            const std::string path) {

        SurfaceMesh* out_mesh = new easy3d::SurfaceMesh(*mesh);

        auto color_prop = out_mesh->face_property<easy3d::vec3>("f:color");

        for(size_t i = 0; i < components.size(); i++) {

            auto color = easy3d::random_color();

            auto component = components[i];
            auto original_face_index = component->face_property<int>("f:original_index");

            for(auto face : component->faces()) {
                int index = original_face_index[face];
                auto f = easy3d::SurfaceMesh::Face(index);
                color_prop[f] = color;
            }
        }

        SurfaceMeshIO::save(path, out_mesh);
    }

    void save_field(const easy3d::SurfaceMesh* mesh,
                     const std::vector<float>& field,
                     const std::string path) {
        SurfaceMesh* out_mesh = new easy3d::SurfaceMesh(*mesh);

        auto color_prop = out_mesh->face_property<easy3d::vec3>("f:color");

        int num = field.size();

        Eigen::VectorXd Z(num);
        float max_value = 0, min_value = 1e9;
        for (int i = 0; i < num; i++) {
            Z[i] = field[i];
            max_value = max(field[i], max_value);
            min_value = min(field[i], min_value);
        }
        std::cout << "max_value: " << max_value << std::endl;
        std::cout << "min_value: " << min_value << std::endl;
        // Z.conservativeResize(ct);
        Eigen::MatrixXd Ct;
        igl::jet(Z, true, Ct);
        for(auto f : out_mesh->faces()) {
            color_prop[f] = easy3d::vec3(Ct(f.idx(), 0),Ct(f.idx(), 1), Ct(f.idx() ,2));
        }

        easy3d::SurfaceMeshIO::save(path, out_mesh);
    }

    void save_fillet_segmentation(const SurfaceMesh* mesh,
                             const std::vector<int>& fillet_label,
                             const string path) {
        SurfaceMesh* out_mesh = new easy3d::SurfaceMesh(*mesh);

        auto color_prop = out_mesh->face_property<easy3d::vec3>("f:color");
        vec3 fillet_color = easy3d::vec3(1.0, 0,0);
        vec3 non_fillet_color = easy3d::vec3(0, 0,1.0);

        for(auto face : out_mesh->faces()) {
            int label = fillet_label[face.idx()];
            color_prop[face] = label != 0 ? fillet_color : non_fillet_color;
        }
        SurfaceMeshIO::save(path, out_mesh);

    }

}

