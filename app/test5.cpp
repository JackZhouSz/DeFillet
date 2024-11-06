//
// Created by xiaowuga on 2024/11/3.
//
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/renderer/drawable_points.h>
#include <Eigen/Dense>

void kmeans(const Eigen::MatrixXd& data, int num_clusters, Eigen::MatrixXd& centroids, Eigen::VectorXi& labels) {
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
        // std::cout << idx << ' ' <<cur.x() << ' ' << cur.y() << ' ' << cur.z() <<std::endl;
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


void sqdist_omp(Eigen::MatrixXd& a,Eigen::MatrixXd& b,
                   Eigen::MatrixXd& res) {
    int n = b.cols(), m = a.cols();
    res.resize(n, m);
#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res(i,j) = (b.col(i) - a.col(j)).squaredNorm();
        }
    }
}

bool fuzzy_clustering(Eigen::MatrixXd& P, Eigen::MatrixXd& C, int nb_cluster, bool robust = true) {
    int nb_points = P.cols();
    int dim = P.rows();
    if(nb_cluster < 1 || nb_cluster > nb_points) {
        std::cerr << "Incorrect setting for nb_cluster" << std::endl;
        return false;
    }

    Eigen::VectorXd onesDPx1 = Eigen::VectorXd::Ones(dim); // DP: dim of points
    Eigen::RowVectorXd ones1xNP = Eigen::RowVectorXd::Ones(nb_points); // NP: number of points
    Eigen::VectorXd onesNCx1 = Eigen::VectorXd::Ones(nb_cluster); // NC: number of cluster

    int maxIter = 100; // maximum number of iterations
    float minDiff = 1e-6; // minimum difference

    Eigen::MatrixXd U = Eigen::MatrixXd::Random(nb_cluster, nb_points);
    U = U.array() / (onesNCx1 * U.colwise().sum()).array();
    Eigen::MatrixXd U2 = U.array().square();
    Eigen::VectorXd rowSumU2 = U2.rowwise().sum(); // nb_cluster x 1
    C = (P * U2.transpose()).array() / (onesDPx1 * rowSumU2.transpose()).array();
    Eigen::MatrixXd dist2; // nb_cluster x nb_points
    sqdist_omp(P, C, dist2);

    // return true;
    double obj = (U2.array() * dist2.array()).sum();
    for(int iter = 0; iter < maxIter; iter++) {
        float obj_old = obj;
        Eigen::MatrixXd tmp = dist2.array().inverse();
        U = tmp.array() / (onesNCx1 * tmp.colwise().sum()).array();
        U2 = U.array().square();
        rowSumU2 = U2.rowwise().sum();
        C = (P * U2.transpose()).array() / (onesDPx1 * rowSumU2.transpose()).array();
        sqdist_omp(P, C, dist2);
        Eigen::MatrixXd U2D2 = U2.array() * dist2.array();
        obj = U2D2.sum();
        if(fabs(obj - obj_old) < minDiff) {
            break;
        }

    }
    return true;
}


void find_cluster_center(const std::vector<easy3d::vec3>& data, easy3d::vec3& centroid, std::vector<easy3d::vec3>& inter_data) {
    int num = data.size();
    if(num < 10) {
        std::cout << "point too few" <<std::endl;
        return;
    }
    inter_data.clear();
    std::vector<bool> labels(num, true);
    // sor(data, labels);
    // centroid = easy3d::vec3(0,0,0);
    std::vector<easy3d::vec3> clean_data = data;
    // for(int i = 0; i < num; i++) {
    //     if(labels[i]) {
    //         centroid += data[i];
    //         clean_data.emplace_back(data[i]);
    //     }
    // }
    // num = clean_data.size();
    // centroid /= num;

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
    for(int i = 0; i < num; i++) {
        if(dist[i] - med < 5.2 * mad) {
            new_centroid += clean_data[i];
            inter_data.emplace_back(clean_data[i]);
            ct++;
        }
    }
    // std::cout << ct << ' ' << new_centroid << std::endl;
    new_centroid /= ct;
    centroid = new_centroid;
}

int main() {

    easy3d::PointCloud* cloud = easy3d::PointCloudIO::load("D:\\code\\defillet\\scvv.ply");
    easy3d::vec3 centroid(0,0,0);
    for(auto v : cloud->vertices()) {
        centroid += cloud->position(v);
    }
    centroid /= cloud->n_vertices();
    std::vector<easy3d::vec3> data = cloud->points();
    for(int i = 0; i < 10; i++) {

        std::vector<easy3d::vec3> inter_data;
        find_cluster_center(data, centroid, inter_data);
        // std::cout << centroid << std::endl;
        easy3d::Viewer viewer("SAD");
        std::cout << centroid << std::endl;
        std::vector<easy3d::vec3> p1 = cloud->points();
        easy3d::PointsDrawable* pd1 = new easy3d::PointsDrawable("c1");
        pd1->update_vertex_buffer(p1); pd1->set_color(easy3d::vec4(1.0,0,0,1.0)); pd1->set_point_size(5);
        pd1->set_impostor_type(easy3d::PointsDrawable::SPHERE);

        std::vector<easy3d::vec3> p2 = {centroid};
        easy3d::PointsDrawable* pd2 = new easy3d::PointsDrawable("c2");
        pd2->update_vertex_buffer(p2); pd2->set_color(easy3d::vec4(0.0,1,0,1.0)); pd2->set_point_size(20);
        pd2->set_impostor_type(easy3d::PointsDrawable::SPHERE);

        std::vector<easy3d::vec3> p3 = inter_data;
        easy3d::PointsDrawable* pd3 = new easy3d::PointsDrawable("c3");
        pd3->update_vertex_buffer(p3); pd3->set_color(easy3d::vec4(0.0,0,1,1.0)); pd3->set_point_size(10);
        pd3->set_impostor_type(easy3d::PointsDrawable::SPHERE);

        for(int i = 0; i < data.size(); i++) {

        }
        viewer.add_drawable(pd1);
        viewer.add_drawable(pd2);
        viewer.add_drawable(pd3);

        viewer.run();
        data = inter_data;
    }
    // Eigen::MatrixXd data, centroids;
    // Eigen::VectorXi labels;
    // data.resize(3,cloud->n_vertices());
    // for(auto v : cloud->vertices()) {
    //     auto p = cloud->position(v);
    //     data.col(v.idx()) = Eigen::Vector3d(p.x,p.y,p.z);
    // }
    // fuzzy_clustering(data, centroids, 2);
    // // kmeans(data, 2, centroids, labels);
    // // std::cout << centroids.rows() << ' ' << centroids.cols() << std::endl;
    // // std::cout <<centroids(0,0) <<' ' << centroids(1,0) << ' ' << centroids(2, 0) << std::endl;
    // // std::cout <<centroids(0,1) <<' ' << centroids(1,1) << ' ' << centroids(2, 1) << std::endl;
    // std::vector<easy3d::vec3> c;
    // for(int i = 0; i < centroids.cols(); i++) {
    //     c.emplace_back(easy3d::vec3(centroids(0, i), centroids(1, i), centroids(2, i)));
    //     // c.emplace_back(easy3d::vec3(centroids(i, 0), centroids(i, 1), centroids(i, 2)));
    // }
    // easy3d::Viewer viewer("SAD");
    // easy3d::PointsDrawable* cen = new easy3d::PointsDrawable("c");
    // cen->update_vertex_buffer(c);cen->set_color(easy3d::vec4(0.0,0,1,1.0)); cen->set_point_size(20);
    // cen->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    // viewer.add_model(cloud);
    // viewer.add_drawable(cen);
    // viewer.run();
    // for(int i = 0; i < 10; i++) {
    //     kmeans(data, 2, centroids, labels);
    //     easy3d::Viewer viewer("SAD");
    //     std::vector<easy3d::vec3> c1;
    //     std::vector<easy3d::vec3> c2;
    //     std::vector<easy3d::vec3> c = {easy3d::vec3(centroids(0, 0), centroids(0, 1), centroids(0, 2)),
    //                                    easy3d::vec3(centroids(1, 0), centroids(1, 1), centroids(1, 2)), };
    //     easy3d::PointsDrawable* pd1 = new easy3d::PointsDrawable("c1");
    //     easy3d::PointsDrawable* pd2 = new easy3d::PointsDrawable("c2");
    //     easy3d::PointsDrawable* cen = new easy3d::PointsDrawable("c");
    //     easy3d::PointsDrawable* points = new easy3d::PointsDrawable("pp");
    //     for(int j = 0; j < labels.size(); j++) {
    //         if(labels[j] == 0) {
    //             c1.emplace_back(easy3d::vec3(data(j, 0), data(j, 1), data(j, 2)));
    //         }
    //         else {
    //             c2.emplace_back(easy3d::vec3(data(j, 0), data(j, 1), data(j, 2)));
    //         }
    //     }
    //
    //     pd1->update_vertex_buffer(c1); pd1->set_color(easy3d::vec4(1.0,0,0,1.0)); pd1->set_point_size(10);
    //     pd1->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    //     pd2->update_vertex_buffer(c2); pd2->set_color(easy3d::vec4(0.0,1,0,1.0));pd2->set_point_size(10);
    //     pd2->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    //     cen->update_vertex_buffer(c);cen->set_color(easy3d::vec4(0.0,0,1,1.0)); cen->set_point_size(20);
    //     cen->set_impostor_type(easy3d::PointsDrawable::SPHERE);
    //     viewer.add_drawable(pd1);
    //     viewer.add_drawable(pd2);
    //     viewer.add_drawable(cen);
    //     points->update_vertex_buffer(cloud->points());
    //     points->set_point_size(5);
    //     points->set_color(easy3d::vec4(0.0,1,1,1.0));
    //     viewer.add_drawable(points);
    //     viewer.run();
    //     std::vector<easy3d::vec3> dd;
    //     if(c1.size() > c2.size()) {
    //         dd = c1;
    //     } else {
    //         dd = c2;
    //     }
    //     std::cout << c1.size() << ' ' << c2.size() << ' ' << c1.size() + c2.size() << ' ' << cloud->n_vertices() << ' ' << 1.0 * (c1.size() + c2.size()) / cloud->n_vertices() << std::endl;
    //     data.resize(dd.size(), 3);
    //     for(int j = 0; j < dd.size(); j++) {
    //         data.row(j) = Eigen::Vector3d(dd[j].x, dd[j].y,dd[j].z);
    //     }
    // }
    return 0;
}