//
// Created by xiaowuga on 12/17/2023.
//

#ifndef DEFILLET_UTILS_H
#define DEFILLET_UTILS_H

#include "common.h"
#include <limits>
#include <CGAL/squared_distance_3.h>



#include <opencv2/opencv.hpp>

inline void select_points_inside_bounding_box(std::vector<Point>& sites, const std::vector<Point>& points, std::vector<Point>& inside_points) {
    inside_points.clear();
    int nb_points = points.size();
    auto bbox = CGAL::bounding_box(sites.begin(), sites.end());
    for(int i = 0; i < nb_points; i++) {

        if(points[i].x() >= bbox.xmin() && points[i].x() <= bbox.xmax()
           && points[i].y() >= bbox.ymin() && points[i].y() <= bbox.ymax()
           && points[i].z() >= bbox.zmin() && points[i].z() <= bbox.zmax()) {
            inside_points.emplace_back(points[i]);
        }
    }
}

inline void select_vertices_and_regions_inside_sites_bounding_box(std::vector<Point>& sites,
                                                                  const std::vector<Point>& vertices,
                                                                  const std::vector<std::vector<int>>& regions,
                                                                  std::vector<Point>& selected_vertices,
                                                                  std::vector<std::vector<int>> & selected_regions) {
    selected_vertices.clear();
    selected_regions.clear();

    Point bboxMin, bboxMax;
    auto bbox = CGAL::bounding_box(sites.begin(), sites.end());

    bboxMin = bbox.min(); bboxMax = bbox.max();

    double diagonal = std::sqrt(CGAL::squared_distance(bboxMin, bboxMax));
    double scale = 0.01;
    double eps = scale * diagonal;

    int nb_vertices = vertices.size();
    std::vector<int> is_inside(nb_vertices);
    for(int i = 0; i < nb_vertices; i++) {
        const Point& p = vertices[i];
        if(p.x() >= bboxMin.x() - eps && p.x() <= bboxMax.x() + eps
           && p.y() >= bboxMin.y() - eps && p.y() <= bboxMax.y() + eps
           && p.z() >= bboxMin.z() - eps && p.z() <= bboxMax.z() + eps) {
            is_inside[i] = selected_vertices.size();
            selected_vertices.emplace_back(p);

        }
        else {
            is_inside[i] = -1;
        }
    }

    int nb_sites = sites.size();
    for(int i = 0; i < nb_sites; i++) {
        const std::vector<int>& region = regions[i];
        std::vector<int> tmp;
        int nb_region = region.size();
        for(int j = 0; j < nb_region; j++) {
            int id = is_inside[region[j]];
            if(id != -1) {
                tmp.emplace_back(id);
            }
        }
        selected_regions.emplace_back(tmp);
    }
}


inline void double_vector1D_data_normalize(std::vector<double>& data) {
    double min_value = std::numeric_limits<double>::max();
    double max_value = std::numeric_limits<double>::min();

    int nb_data = data.size();
    for(int i = 0; i < nb_data; i++) {
        min_value = std::min(min_value, data[i]);
        max_value = std::max(max_value, data[i]);
    }

    for(int i = 0; i < nb_data; i++) {
        data[i] = (data[i] - min_value) / (max_value - min_value);
    }
}


inline cv::Vec3b getColorFromColormap(double value, cv::ColormapTypes colormap) {
    cv::Mat colorMap(1, 1, CV_8U, cv::Scalar(value));
    cv::Mat colorImage;

    // 应用颜色映射
    cv::applyColorMap(colorMap, colorImage, colormap);

    // 获取生成的颜色
    cv::Vec3b color = colorImage.at<cv::Vec3b>(0, 0);

    return color;
}


inline cv::Mat getColorImageFromValues(const std::vector<double>& values, cv::ColormapTypes colormap) {
    cv::Mat grayImage(values.size(), 1, CV_8U);

    // 将数值填充到灰度图中
    for (int i = 0; i < values.size(); ++i) {
        grayImage.at<uchar>(i, 0) = static_cast<uchar>(values[i] * 255); // 将标量值转换为 0-255 的范围
    }

    // 应用颜色映射
    cv::Mat colorImage;
    cv::applyColorMap(grayImage, colorImage, colormap);

    return colorImage;
}

inline void extract_edge_from_mesh(std::vector<Point> v_pos,
                                    std::vector<std::vector<size_t>>& f_ind,
                                    std::vector<std::pair<int,int>>& edge,
                                    std::vector<double>& edge_weight) {
    int nb_v = v_pos.size();

    std::map<int, std::set<int>> mp;

    for(size_t i = 0; i < f_ind.size(); i++) {
        int num = f_ind[i].size();
        for(size_t j = 0; j < num; j++) {
            int s = f_ind[i][j], t = f_ind[i][(j + 1) % num];
            if(s > t) {
                 std::swap(s, t);
            }
            mp[s].insert(t);
        }
    }

    edge.clear();
    edge_weight.resize(nb_v);

    for(auto& item : mp) {
        int x = item.first;
        for(auto& y : item.second) {
            double len = std::sqrt(CGAL::squared_distance(v_pos[x], v_pos[y]));
            edge.emplace_back(std::make_pair(x, y));
            edge_weight.emplace_back(len);
        }
    }

}

inline void extract_vertex_graph_from_mesh(std::vector<Point>& v_pos,
                                           std::vector<std::vector<size_t>>& f_ind,
                                           std::vector<std::vector<std::pair<int,double>>>& graph) {
    int nb_v = v_pos.size();

    std::map<int, std::set<int>> mp;

    for(size_t i = 0; i < f_ind.size(); i++) {
        int num = f_ind[i].size();
        for(size_t j = 0; j < num; j++) {
            int s = f_ind[i][j], t = f_ind[i][(j + 1) % num];
            mp[s].insert(t);
            mp[t].insert(s);
        }
    }

    graph.clear();
    graph.resize(nb_v);
    for(auto& item : mp) {
        int x = item.first;
        for(auto& y : item.second) {
            double len = std::sqrt(CGAL::squared_distance(v_pos[x], v_pos[y]));
            graph[x].emplace_back(std::make_pair(y, len));
        }
    }
}

inline void densify_mesh_via_centroid_insertion(const std::vector<Point>& points,
                                                const std::vector<std::vector<size_t>>& faces,
                                                const std::vector<bool>& faces_labels,
                                                std::vector<Point>& dense_points,
                                                std::vector<std::vector<size_t>>& dense_faces,
                                                std::vector<int>& centroid,
                                                std::vector<bool>& centroid_label) {
    dense_points = points;
    int nb_faces = faces.size();
    for(int i = 0; i < nb_faces; i++) {
        double sum_x = 0, sum_y = 0, sum_z = 0;
        int num = faces[i].size();
        for(int j = 0; j < num; j++) {
            sum_x += points[faces[i][j]].x();
            sum_y += points[faces[i][j]].y();
            sum_z += points[faces[i][j]].z();
        }
        size_t cen_id = dense_points.size();
        centroid.emplace_back(cen_id);
        dense_points.emplace_back(Point(sum_x / num, sum_y / num, sum_z / num));
        centroid_label.emplace_back(faces_labels[i]);
        for(int j = 0; j < num; j++) {
            std::vector<size_t> tmp = {faces[i][j], faces[i][(j + 1) % num], cen_id};
            dense_faces.emplace_back(tmp);
        }
    }

}
#endif //DEFILLET_UTILS_H
