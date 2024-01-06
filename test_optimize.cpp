//
// Created by 13900K on 2024/1/6.
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/cat.h>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <cmath>
#include "visualization.h"
#include "defillet.h"

bool optimize_sparseLU(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       std::vector<Eigen::Vector3d>& normals,
                       std::vector<Eigen::Vector3d>& new_points,
                       std::vector<int>& fixed_points,
                       std::vector<std::pair<int,std::pair<int,Eigen::Vector3d>>>& edge_vector) {
    int nb_faces = faces.size();
    int nb_points = points.size();
    int nb_fixed_points = fixed_points.size();
    int nb_edges = edge_vector.size();
    std::cout << "run optimize..." << std::endl;
    // normals constraint
    std::cout << "add normals constraint" << std::endl;
    Eigen::SparseMatrix<double> NC(nb_faces * 2, nb_points * 3);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(nb_faces * 2 + nb_edges);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < nb_faces; i++) {
        double n0 = normals[i].x(), n1 = normals[i].y(), n2 = normals[i].z();
        int v0 = faces[i][0], v1 = faces[i][1], v2 = faces[i][2];
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0, n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 1, n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 2, n2));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1, -n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 1, -n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 2, -n2));

        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0, n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 1, n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 2, n2));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2, -n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 1, -n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 2, -n2));
    }
    NC.setFromTriplets(triplets.begin(), triplets.end());

    // edges vector constraint
    std::cout << "add edges vector constraint" << std::endl;

    Eigen::SparseMatrix<double> EVC(nb_edges, nb_points * 3);
    triplets.clear();
    for(int i = 0; i < nb_edges; i++) {
        int id1 = edge_vector[i].first;
        int id2 = edge_vector[i].second.first;
        const Eigen::Vector3d& v = edge_vector[i].second.second;
        double c1 = v.x(), c2 = v.y(), c3 = v.z();
        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1, (-c2 + c3)));
        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2, (c2 - c3)));

        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 1, (c1 - c3)));
        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 1, (-c1 + c3)));

        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id1 + 2, (-c1 + c2)));
        triplets.emplace_back(Eigen::Triplet<double>(i,  3 * id2 + 2, (c1 - c2)));
    }
    EVC.setFromTriplets(triplets.begin(), triplets.end());


    // fixed points constraint
    std::cout << "add fixed points constraint" << std::endl;
    Eigen::VectorXd d = Eigen::VectorXd(nb_fixed_points * 3);
    std::cout << nb_fixed_points << std::endl;
    triplets.clear();
    Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
    std::cout << nb_points << std::endl;
    for(int i = 0; i < nb_fixed_points; i++) {
        int id = fixed_points[i];
        triplets.emplace_back(Eigen::Triplet<double>(3 * i,  3 * id, 1.0)); d[ 3 * i] = points[id].x();
        triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1,  3 * id + 1, 1.0)); d[3 * i + 1] = points[id].y();
        triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  3 * id + 2, 1.0)); d[3 * i + 2] = points[id].z();
    }
    E.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "run solver..." << std::endl;


    Eigen::SparseMatrix<double> A;
    igl::cat(1, NC, EVC, A);
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> Q = AT * A;
    Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
    Eigen::SparseMatrix<double> ET = E.transpose();
    Eigen::SparseMatrix<double> tempMat1;
    Eigen::SparseMatrix<double> tempMat2;
    Eigen::SparseMatrix<double> L;
    igl::cat(1, Q, E, tempMat1);
    igl::cat(1, ET, zero, tempMat2);
    igl::cat(0, tempMat1, tempMat2, L);
    Eigen::VectorXd b(nb_points * 3 + d.size());
    b << Eigen::VectorXd::Zero(nb_points * 3), d;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver; // 创建SparseQR求解器
    std::cout << L.rows() << ' ' << L.cols() << std::endl;
    solver.compute(L);
    if(solver.info()!= Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
        return false;
    }
    Eigen::VectorXd x = solver.solve(b);
    if(solver.info()!= Eigen::Success) {
        // solving failed
        std::cout << "solving failed" << std::endl;
        return false;
    }
    new_points.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        new_points[i] = Eigen::Vector3d (x[3 * i], x[3 * i + 1], x[3 * i + 2]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
    }
    std::cout << "done." << std::endl;
    return true;
}

bool optimize_sparseLU2(const std::vector<Eigen::Vector3d>& points,
                       const std::vector<std::vector<size_t>>& faces,
                       std::vector<Eigen::Vector3d>& normals,
                       std::vector<Eigen::Vector3d>& new_points,
                       std::vector<int>& fixed_points) {
    int nb_faces = faces.size();
    int nb_points = points.size();
    int nb_fixed_points = fixed_points.size();
    std::cout << "run optimize..." << std::endl;
    // normals constraint
    std::cout << "add normals constraint" << std::endl;
    Eigen::SparseMatrix<double> NC(nb_faces * 2, nb_points * 3);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(nb_faces * 2);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < nb_faces; i++) {
        double n0 = normals[i].x(), n1 = normals[i].y(), n2 = normals[i].z();
        int v0 = faces[i][0], v1 = faces[i][1], v2 = faces[i][2];
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0, n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 1, n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v0 + 2, n2));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1, -n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 1, -n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i, 3 * v1 + 2, -n2));

        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0, n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 1, n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v0 + 2, n2));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2, -n0));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 1, -n1));
        triplets.emplace_back(Eigen::Triplet<double>(2 * i + 1, 3 * v2 + 2, -n2));
    }
    NC.setFromTriplets(triplets.begin(), triplets.end());



    // fixed points constraint
    std::cout << "add fixed points constraint" << std::endl;
    Eigen::VectorXd d = Eigen::VectorXd(nb_fixed_points * 3);
    std::cout << nb_fixed_points << std::endl;
    triplets.clear();
    Eigen::SparseMatrix<double> E(nb_fixed_points * 3, nb_points * 3);
    std::cout << nb_points << std::endl;
    for(int i = 0; i < nb_fixed_points; i++) {
        int id = fixed_points[i];
        triplets.emplace_back(Eigen::Triplet<double>(3 * i,  3 * id, 1.0)); d[ 3 * i] = points[id].x();
        triplets.emplace_back(Eigen::Triplet<double>(3 * i + 1,  3 * id + 1, 1.0)); d[3 * i + 1] = points[id].y();
        triplets.emplace_back(Eigen::Triplet<double>(3 * i + 2,  3 * id + 2, 1.0)); d[3 * i + 2] = points[id].z();
    }
    E.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "run solver..." << std::endl;


    Eigen::SparseMatrix<double> A = NC;
    Eigen::SparseMatrix<double> AT = A.transpose();
    Eigen::SparseMatrix<double> Q = AT * A;
    Eigen::SparseMatrix<double> zero(nb_fixed_points * 3, nb_fixed_points * 3);
    Eigen::SparseMatrix<double> ET = E.transpose();
    Eigen::SparseMatrix<double> tempMat1;
    Eigen::SparseMatrix<double> tempMat2;
    Eigen::SparseMatrix<double> L;
    igl::cat(1, Q, E, tempMat1);
    igl::cat(1, ET, zero, tempMat2);
    igl::cat(0, tempMat1, tempMat2, L);
    Eigen::VectorXd b(nb_points * 3 + d.size());
    b << Eigen::VectorXd::Zero(nb_points * 3), d;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver; // 创建SparseQR求解器
    std::cout << L.rows() << ' ' << L.cols() << std::endl;
    for (int i = 0; i < L.rows(); ++i) {
        for (int j = 0; j < L.cols(); ++j) {
            double value = L.coeff(i, j);
            double value1 = L.coeff(j, i);
            if(value != value1) {
                std::cout << "ASD" << std::endl;
            }
//            std::cout  << value << " ";
        }
//        std::cout << std::endl;
    }
    solver.compute(L);
    if(solver.info()!= Eigen::Success) {
        std::cout << "decomposition failed" << std::endl;
        return false;
    }
    Eigen::VectorXd x = solver.solve(b);
    if(solver.info()!= Eigen::Success) {
        // solving failed
        std::cout << "solving failed" << std::endl;
        return false;
    }
    new_points.resize(nb_points);
    for(int i = 0; i < nb_points; i++) {
        new_points[i] = Eigen::Vector3d (x[3 * i], x[3 * i + 1], x[3 * i + 2]);
//            std::cout << x[3 * i] << ' ' << x[3 * i + 1] << ' ' << x[3 * i + 2] << std::endl;
    }
    std::cout << "done." << std::endl;
    return true;
}

void test_data1(std::vector<Eigen::Vector3d>& points,
                std::vector<std::vector<size_t>>& faces,
                std::vector<Eigen::Vector3d>& normals,
                std::vector<int>& fixed_points,
                std::vector<std::pair<int,std::pair<int,Eigen::Vector3d>>>& edge_vector) {
    Eigen::Vector3d O( 0, 0.5, 0);
    Eigen::Vector3d A(-std::sqrt(3)/ 2, 1.0 / 2.0, 0);
    Eigen::Vector3d B(-std::sqrt(3)/ 2, -1.0 / 2.0, 0);
    Eigen::Vector3d C(std::sqrt(3)/ 2, 1.0 / 2.0, 0);
    Eigen::Vector3d D(std::sqrt(3)/ 2, -1.0 / 2.0, 0);
    Eigen::Vector3d E(0, 0, 1.0);
    Eigen::Vector3d on1 = (A - E).cross(B - A).normalized();
    Eigen::Vector3d on2 = (E - C).cross(D - E).normalized();
    Eigen::Vector3d tn1 = Eigen::Vector3d(0,0,1.0);
    Eigen::Vector3d tn2 = Eigen::Vector3d(0,0,1.0);

    Eigen::Vector3d evAE = (A - O).normalized();
    Eigen::Vector3d evAB = (A - B).normalized();
    Eigen::Vector3d evBE = (B - O).normalized();

    Eigen::Vector3d evCE = (C - O).normalized();
    Eigen::Vector3d evCD = (C - D).normalized();
    Eigen::Vector3d evDE = (D - O).normalized();

    points.clear();
    points.emplace_back(A);
    points.emplace_back(B);
    points.emplace_back(C);
    points.emplace_back(D);
    points.emplace_back(E);
    faces = {{0, 1, 4}, {2,4,3}};
    normals.clear();
    normals.emplace_back(tn1);
    normals.emplace_back(tn2);
    fixed_points = {0, 1, 2, 3};
    edge_vector;

    edge_vector.emplace_back(std::make_pair(0, std::make_pair(4, evAE)));
    edge_vector.emplace_back(std::make_pair(0, std::make_pair(1, evAB)));
    edge_vector.emplace_back(std::make_pair(1, std::make_pair(4, evBE)));

    edge_vector.emplace_back(std::make_pair(2, std::make_pair(4, evCE)));
    edge_vector.emplace_back(std::make_pair(2, std::make_pair(3, evCD)));
    edge_vector.emplace_back(std::make_pair(3, std::make_pair(4, evDE)));
}

void test_data2(std::vector<Eigen::Vector3d>& points,
                std::vector<std::vector<size_t>>& faces,
                std::vector<Eigen::Vector3d>& normals,
                std::vector<int>& fixed_points,
                std::vector<std::pair<int,std::pair<int,Eigen::Vector3d>>>& edge_vector) {
    Eigen::Vector3d TB( 0, 0.5, 0);
    Eigen::Vector3d TC( 0, -0.5, 0);
    Eigen::Vector3d A(-std::sqrt(3)/ 2, 0, 0);
    Eigen::Vector3d B(0, -1.0 / 2.0, 1.0);
    Eigen::Vector3d C(0, 1.0 / 2.0, 1.0);
    Eigen::Vector3d D(std::sqrt(3)/ 2, 0, 0);

//    Eigen::Vector3d on1 = (A - E).cross(B - A).normalized();
//    Eigen::Vector3d on2 = (E - C).cross(D - E).normalized();
    Eigen::Vector3d tn1 = Eigen::Vector3d(0,0,1.0);
    Eigen::Vector3d tn2 = Eigen::Vector3d(0,0,1.0);

    Eigen::Vector3d evAB = (A - TB).normalized();
    Eigen::Vector3d evAC = (TC - A).normalized();

    Eigen::Vector3d evBD = (TB - D).normalized();
    Eigen::Vector3d evCD = (D - TC).normalized();

    points.clear();
    points.emplace_back(A);
    points.emplace_back(B);
    points.emplace_back(C);
    points.emplace_back(D);
    faces = {{0, 2, 1}, {2,3, 1}};
    normals.clear();
    normals.emplace_back(tn1);
    normals.emplace_back(tn2);
    fixed_points = {0, 3};

    edge_vector.emplace_back(std::make_pair(0, std::make_pair(1, evAB)));
    edge_vector.emplace_back(std::make_pair(0, std::make_pair(2, evAC)));

    edge_vector.emplace_back(std::make_pair(1, std::make_pair(3, evBD)));
    edge_vector.emplace_back(std::make_pair(2, std::make_pair(3, evCD)));
}

void test_data3(int num,
                std::vector<Eigen::Vector3d>& points,
                std::vector<std::vector<size_t>>& faces,
                std::vector<Eigen::Vector3d>& normals,
                std::vector<int>& fixed_points,
                std::vector<std::pair<int,std::pair<int,Eigen::Vector3d>>>& edge_vector) {
    double r = 0.5;
    points.clear();
    Eigen::Vector3d v1(r, -0.5, 0), offset(0, 0, -0.5);
    Eigen::Vector3d v2 = v1 + offset;
    points.emplace_back(v1); points.emplace_back(v2);
    for(int i = 0; i <= num; i++) {
        double angle = static_cast<double>(i) * M_PI / (2.0 * num);
        double x = r * std::cos(angle); // 计算 x 坐标
        double y = r * std::sin(angle); // 计算 y 坐标
        Eigen::Vector3d u(x, y, 0);
        Eigen::Vector3d v = u + offset;
//        std::cout << u.x() << ' ' << u.y() << ' ' << u.z() << std::endl;
//        std::cout << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
        points.emplace_back(u); points.emplace_back(v);
    }
    Eigen::Vector3d v3(-0.5, r, 0);
    Eigen::Vector3d v4 = v3 + offset;
    points.emplace_back(v3); points.emplace_back(v4);
    faces.clear();
    int nb_points = points.size();
    std::cout << nb_points << std::endl;
    for(int i = 0; i < nb_points - 2; i += 2) {
        std::vector<size_t> tmp1 = {(size_t)i, (size_t)i + 1, (size_t)i + 3};
        faces.emplace_back(tmp1);
        std::vector<size_t> tmp2 = {(size_t)i, (size_t)i + 3, (size_t)i + 2};
        faces.emplace_back(tmp2);
    }
    std::vector<Point> pp;
    for(int i = 0; i < nb_points; i++) {
        pp.emplace_back(Point(points[i].x(), points[i].y(), points[i].z()));
    }
//    points_visualization(pp);
//    mesh_visualization(pp, faces);

    int nb_faces = faces.size();
    std::vector<Point> dense_points;
    std::vector<std::vector<size_t>> dense_faces;
    std::vector<int> centroid;
    std::vector<bool> centroid_label;
    std::vector<bool> face_labels(faces.size(), true);
    face_labels[0] = false; face_labels[1] = false;
    face_labels[nb_faces - 1] = false; face_labels[nb_faces - 2] = false;
    densify_mesh_via_centroid_insertion(pp, faces,
                                        face_labels, dense_points, dense_faces,
                                        centroid, centroid_label);

    std::vector<int> ancestor;
    std::vector<double> distance;
    DEFILLET::run_geodesic(dense_points, dense_faces, centroid,
                           centroid_label, ancestor, distance);
//    double_vector1D_data_normalize(distance);
//    mesh_field_visualization(pp, faces, distance);
//    mesh_visualization(dense_points, dense_faces);
    std::vector<Vector_3> target_normals;

    std::vector<std::pair<int,std::pair<int,Vector_3>>> asd;

    DEFILLET::geodesic_post_processing(pp, faces, ancestor, target_normals,
                                       fixed_points, asd, 10.0);
    fixed_points = {0, 1, 2, 3, nb_points - 1, nb_points - 2, nb_points - 3, nb_points - 4};
    int evf = asd.size();
    for(int i = 0; i < nb_faces; i++) {
        normals.emplace_back(Eigen::Vector3d(target_normals[i].x(), target_normals[i].y(),target_normals[i].z()));
    }
    for(int i = 0; i < evf; i++) {
        int x = asd[i].first;
        int y = asd[i].second.first;
        Vector_3 p = asd[i].second.second;
        edge_vector.emplace_back(std::make_pair(x, std::make_pair(y,
                                                                  Eigen::Vector3d(p.x(), p.y(), p.z()))));
    }
    mesh_face_normals_vector_field(pp, faces, target_normals);

}

int main() {
    std::vector<Eigen::Vector3d> points;
    std::vector<std::vector<size_t>> faces;
    std::vector<Eigen::Vector3d> normals;
    std::vector<int> fixed_points;
    std::vector<std::pair<int,std::pair<int, Eigen::Vector3d>>> edge_vector;
//    test_data1( points,faces,normals,fixed_points,edge_vector);
//    test_data2( points,faces,normals,fixed_points,edge_vector);
    test_data3(2, points,faces,normals,fixed_points,edge_vector);
//    return 0;
    std::vector<Eigen::Vector3d> new_points;
//    optimize_sparseLU(points, faces, normals, new_points, fixed_points, edge_vector);
    optimize_sparseLU2(points, faces, normals, new_points, fixed_points);

    std::cout << std::endl;
    for(int i = 0; i < new_points.size(); i++) {
        std::cout << new_points[i].x() << ' '<< new_points[i].y() << ' ' << new_points[i].z() << std::endl;
    }
    std::vector<Point> pp;
    for(int i = 0; i < points.size(); i++) {
        pp.emplace_back(Point(new_points[i].x(), new_points[i].y(), new_points[i].z()));
    }
    mesh_visualization(pp, faces);
    return 0;
}