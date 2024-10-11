//
// Created by xiaowuga on 2024/10/9.
//

#include <fillet_seg.h>
#include <fillet_seg_v2.h>
#include <surafce_mesh_segmenter.h>

#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/algo/surface_mesh_geometry.h>

#include "surafce_mesh_segmenter.h"
#include "fillet_seg.h"

double dihedral_angle(easy3d::SurfaceMesh* mesh
                      , easy3d::SurfaceMesh::Face f1
                      , easy3d::SurfaceMesh::Face f2
                      , bool rad = false) {
    easy3d::vec3 n1 = mesh->compute_face_normal(f1).normalize();
    easy3d::vec3 n2 = mesh->compute_face_normal(f2).normalize();

    double radians = abs(acos(dot(n1, n2)));

    if(rad) {
        return radians;
    } else {
        double degrees = radians * 180.0 / M_PI;
        return degrees;
    }
}

std::vector<easy3d::vec3> load_patch(easy3d::SurfaceMesh* mesh, easy3d::SurfaceMesh::Face f, double radius) {
    double area_thr = 2 * M_PI * radius * radius;
    // std::set<easy3d::SurfaceMesh::Face> vis;
    std::vector<easy3d::vec3> sites;
    double area_sum = 0;

    auto face_area = mesh->face_property<double>("f:area");
    auto vis = mesh->face_property<int>("f:label");
    for(auto ff : mesh->faces()) {
        vis[ff] = 0;
    }
    std::queue<easy3d::SurfaceMesh::Face> que;
    que.push(f); vis[f] = 1;

    while(area_sum < area_thr && !que.empty()) {
        auto cur_f = que.front(); que.pop();
        area_sum += face_area[cur_f];
        for(auto h : mesh->halfedges(cur_f)) {
            auto nxt_f = mesh->face(mesh->opposite(h));
            if(nxt_f.is_valid() && vis[nxt_f] ==  0
               && dihedral_angle(mesh, cur_f, nxt_f) < 15.0) {
                vis[nxt_f] = 1;
                que.push(nxt_f);
            }
        }
    }

    easy3d::SurfaceMeshSegmenter sms(mesh);
    easy3d::SurfaceMesh* patch = sms.segment(vis, 1);
    FilletSeg fillet_seg;
    fillet_seg.set_mesh(patch);
    fillet_seg.run_scoring();
    easy3d::PointCloud* sor_vertices = new easy3d::PointCloud;
    auto& tmp_vertices = fillet_seg.get_vertices();
    const std::vector<int> sor_labels = fillet_seg.get_sor_labels();
    std::vector<easy3d::vec3> res;
    for(int i = 0; i < sor_labels.size(); i++) {
        if(sor_labels[i]) {
            res. emplace_back(tmp_vertices[i]);
        }
    }
    return res;
}


int main(int argc, char **argv) {

    std::string mesh_path = "D:\\code\\defillet_case\\21164_de08f61e_1\\21164_de08f61e_1.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(mesh_path);

    auto face_area = mesh->add_face_property<double>("f:area");
    for(auto f : mesh->faces()) {
        face_area[f] = easy3d::geom::triangle_area(mesh, f);
    }
    double radius = 0.04;
    easy3d::Box3 box = mesh->bounding_box();
    radius = box.diagonal_length() * radius;
    easy3d::PointCloud* pp = new easy3d::PointCloud;
    for(auto f : mesh->faces()) {
        if(f.idx() % 100 == 0) {
            std::cout << f.idx() << std::endl;
            std::vector<easy3d::vec3> asd = load_patch(mesh, f, radius);
            for(auto v : asd) {
                pp->add_vertex(v);
            }
        }
    }

    easy3d::PointCloudIO::save("../out2.ply", pp);
    return 0;
}