//
// Created by xiaowuga on 2025/1/11.
//

#include<easy3d/core/surface_mesh.h>
#include<easy3d/algo/surface_mesh_sampler.h>
#include <easy3d/kdtree/kdtree_search_nanoflann.h>
#include <easy3d/core/point_cloud.h>
#include <numeric>


#include <knn.h>
#include <easy3d/fileio/surface_mesh_io.h>

struct Eval {
    double cd;
    double nc;
    double f1;
    double ecd;
    double ef1;
};
void threshold_percentage(const std::vector<double>& data, std::vector<double>& is_thr, const std::vector<double>& thr) {
    int num_thr = thr.size();
    is_thr.clear();
    int num = data.size();
    for(int i = 0; i < num_thr; i++) {
        int ct = 0;
        for(int j = 0; j < num; j++) {
            if(data[j] <= thr[i]) {
                ct++;
            }
        }
        if(ct == 0) {
            is_thr.emplace_back(0.0);
        }
        else {
            is_thr.emplace_back(1.0 * ct / num);
        }
    }
}

void distance_p2p(const std::vector<easy3d::vec3>& src_p, const std::vector<easy3d::vec3>& tar_p,
                  std::vector<double>& pdist) {
    // easy3d::KdTreeSearch_NanoFLANN kdt(tar);
    int num = src_p.size();
    pdist.clear();
    pdist.resize(num);

    std::vector<KNN::Point> knn_points;
    for(size_t i = 0; i < tar_p.size(); i++) {
        knn_points.emplace_back(KNN::Point(tar_p[i].x, tar_p[i].y, tar_p[i].z));
    }
    KNN::KdSearch kdt(knn_points);
#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        std::vector<size_t> indices;
        std::vector<double> dist;
        KNN::Point p(src_p[i].x, src_p[i].y, src_p[i].z);
        kdt.kth_search(p, 1, indices, dist);
        int idx = indices[0];
        pdist[i] = (src_p[i] - tar_p[idx]).norm();
    }
}

void distance_p2p(const easy3d::PointCloud* src, const easy3d::PointCloud* tar,
                  std::vector<double>& pdist, std::vector<double>& ndist) {
    // easy3d::KdTreeSearch_NanoFLANN kdt(tar);
    const std::vector<easy3d::vec3>& src_p = src->points();
    const std::vector<easy3d::vec3>& src_n = src->get_vertex_property<easy3d::vec3>("v:normal").vector();
    const std::vector<easy3d::vec3>& tar_p = tar->points();
    const std::vector<easy3d::vec3>& tar_n = tar->get_vertex_property<easy3d::vec3>("v:normal").vector();
    int num = src_p.size();
    pdist.resize(num);
    ndist.resize(num);

    std::vector<KNN::Point> knn_points;
    for(size_t i = 0; i < tar_p.size(); i++) {
        knn_points.emplace_back(KNN::Point(tar_p[i].x, tar_p[i].y, tar_p[i].z));
    }
    KNN::KdSearch kdt(knn_points);
#pragma omp parallel for
    for(int i = 0; i < num; i++) {
        std::vector<size_t> indices;
        std::vector<double> dist;
        KNN::Point p(src_p[i].x, src_p[i].y, src_p[i].z);
        kdt.kth_search(p, 1, indices, dist);
        int idx = indices[0];
        pdist[i] = (src_p[i] - tar_p[idx]).norm();
        easy3d::vec3 n1 = normalize(src_n[i]);
        easy3d::vec3 n2 = normalize(tar_n[idx]);
        ndist[i] = fabs(dot(n1, n2));
    }
}

void edge_sampling(const easy3d::PointCloud* cloud, std::vector<easy3d::vec3>& edge_samples, double radius = 0.004, double dot_thr = 0.2) {
    // easy3d::KdTreeSearch_NanoFLANN kdt(cloud);
    const std::vector<easy3d::vec3>& points = cloud->points();
    const std::vector<easy3d::vec3>& normals = cloud->get_vertex_property<easy3d::vec3>("v:normal").vector();
    int num = cloud->n_vertices();
    std::vector<KNN::Point> knn_points;
    for(size_t i = 0; i < num; i++) {
        knn_points.emplace_back(KNN::Point(points[i].x, points[i].y, points[i].z));
    }
    KNN::KdSearch kdt(knn_points);
    for(int i = 0; i < num; i++) {
        std::vector<size_t> indices;
        std::vector<double> dist;
        KNN::Point p(points[i].x, points[i].y, points[i].z);
        kdt.radius_search(p, radius, indices, dist);
        float minn = 1.0;
        for(size_t j = 0; j < indices.size(); j++) {
            int idx = indices[j];
            minn = std::min(minn, easy3d::dot(normals[idx], normals[i]));
        }
        if(minn < dot_thr) {
            edge_samples.emplace_back(points[i]);
        }
    }
}

void eval_mesh(const easy3d::SurfaceMesh* src, const easy3d::SurfaceMesh* tar, Eval& result, double ef1_thr = 0.004,std::vector<double> cd_thr = std::vector<double>{0.005}) {
    easy3d::SurfaceMeshSampler sampler;
    const easy3d::PointCloud* src_p = sampler.apply(src, 50000);
    const easy3d::PointCloud* tar_p = sampler.apply(tar, 50000);
    // std::cout << "ASD"<< std::endl;
    std::vector<double> completeness, completeness_normals;
    distance_p2p(tar_p, src_p, completeness, completeness_normals);
    std::vector<double> recall;
    threshold_percentage(completeness, recall, cd_thr);
    double mean_completeness = std::accumulate(completeness.begin(), completeness.end(), 0.0) / completeness.size();
    double mean_completeness_normals = std::accumulate(completeness_normals.begin(), completeness_normals.end(), 0.0) / completeness_normals.size();

    std::vector<double> accuracy, accuracy_normals;
    distance_p2p(src_p, tar_p, accuracy, accuracy_normals);
    std::vector<double> precision;
    threshold_percentage(accuracy, precision, cd_thr);
    double mean_accuracy = std::accumulate(accuracy.begin(), accuracy.end(), 0.0) / accuracy.size();
    double mean_accuracy_normals = std::accumulate(accuracy_normals.begin(), accuracy_normals.end(), 0.0) / accuracy_normals.size();
    result.cd = 0.5 * (mean_completeness + mean_accuracy);
    result.cd = mean_completeness;
    result.nc = 0.5 * (mean_accuracy_normals + mean_completeness_normals);
    result.nc = mean_accuracy_normals;
    std::vector<double> F;
    for(size_t i = 0; i < cd_thr.size(); i++) {
        double val = 2 * precision[i] * recall[i] / ( precision[i] + recall[i] + 1e-12);
        F.emplace_back(val);
    }
    result.f1 = F[0];


    std::vector<easy3d::vec3> src_ep, tar_ep;
    edge_sampling(src_p, src_ep);
    edge_sampling(tar_p, tar_ep);

    if(src_ep.size() == 0) {
        for(int i = 0; i < 486; i++) {
            src_ep.emplace_back(easy3d::vec3(0,0,0));
        }
    }

    if(tar_ep.size() == 0) {
        result.ecd = 0.0;
        result.ef1 = 1.0;
        return;
    }
    std::vector<double> dist;
    distance_p2p(tar_ep, src_ep, dist);
    std::vector<double> tmp;
    threshold_percentage(dist, tmp,{ef1_thr});
    double r = tmp[0];
    double gt2pred_mean_ecd =  std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();

    distance_p2p(src_ep, tar_ep, dist);
    threshold_percentage(dist,  tmp, {ef1_thr});
    double p = tmp[0];
    double pred2gt_mean_ecd = std::accumulate(dist.begin(), dist.end(), 0.0) / dist.size();

    result.ecd = 0.5 * (gt2pred_mean_ecd + pred2gt_mean_ecd);
    result.ecd = gt2pred_mean_ecd;
    if(r + p > 1e-12) {
        result.ef1 = 2 * r * p / (r + p);
    }
    else
        result.ef1 = 0.0;

}
void normalize_mesh(easy3d::SurfaceMesh* mesh) {
    easy3d::Box3 box = mesh->bounding_box();
    double val = box.max_range();
    easy3d::vec3 center = box.center();
    for(auto v : mesh->vertices()) {
        mesh->position(v) = (mesh->position(v) - center) / val;
    }
}
int main() {
    std::string dir_path = "D:\\Sig25_DeFillet\\result\\ours\\bianbanjing\\";
    std::vector<std::string> names = {"gt", "a"};
    std::vector<easy3d::SurfaceMesh*> meshes;
    for(int i = 0; i < names.size(); i++) {
        easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(dir_path + names[i] + ".ply");
        normalize_mesh(mesh);
        meshes.emplace_back(mesh);

    }
    for(int i = 1; i < names.size(); i++) {
        Eval result;
        eval_mesh(meshes[i], meshes[0], result);
        std::cout << names[i] <<" "<<result.cd << ' ' << result.nc << ' ' << result.f1  << ' ' << result.ecd << ' ' << result.ef1 << std::endl;
    }


    return 0;
}
