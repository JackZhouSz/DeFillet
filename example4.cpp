//
// Created by 13900K on 2024/3/5.
//
#include "common.h"
#include "MeshVoronoi3d.h"
#include "io.h"
#include "utils.h"
#include "visualization.h"
#include "defillet.h"

#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>

#include <easy3d/algo/surface_mesh_curvature.h>

int maxx = 0;
bool is_connected(const easy3d::SurfaceMesh* mesh, const std::vector<int>& faces) {
    std::set<easy3d::SurfaceMesh::Face> face_set;
    std::set<easy3d::SurfaceMesh::Face> vis;
    for(auto f : faces) {
        easy3d::SurfaceMesh::Face face_handle(f);
        face_set.insert(face_handle);
    }
    easy3d::SurfaceMesh::Face st(faces[0]);
    std::queue<easy3d::SurfaceMesh::Face>que;
    que.push(st);
    while(!que.empty()) {
        auto cur = que.front(); que.pop();
        if(vis.find(cur) == vis.end()) {
            vis.insert(cur);
        }
        else {
            continue;
        }
        for(auto h : mesh->halfedges(cur)) {
            auto opp_f = mesh->face(mesh->opposite(h));
            if(opp_f.is_valid() && face_set.find(opp_f) != face_set.end()) {
                que.push(opp_f);
            }
        }
    }

    return vis.size() == faces.size();

}

void domain_region(easy3d::SurfaceMesh* mesh, easy3d::vec3 center,
                   const std::vector<easy3d::vec3>& sites,  int st) {
    std::set<easy3d::SurfaceMesh::Face> vis;
    auto fillet = mesh->face_property<double>("f:fillet");
    double axis_len = (center - sites[st]).norm();
    easy3d::SurfaceMesh::Face sf(st);
    std::queue<easy3d::SurfaceMesh::Face>que;
    que.push(sf);
    while(!que.empty()) {
        auto cur = que.front(); que.pop();
        if(vis.find(cur) == vis.end()) {
            vis.insert(cur);
//            if(vis.size() > 40) {
//                break;
//            }
        }
        else {
            continue;
        }
        for(auto h : mesh->halfedges(cur)) {
            auto opp_f = mesh->face(mesh->opposite(h));
            if(opp_f.is_valid() && vis.find(opp_f) == vis.end()) {
                double len = (center - sites[opp_f.idx()]).norm();
                if(abs(axis_len - len) / axis_len < 0.05) {
                    que.push(opp_f);
                }
            }
        }
    }
    maxx = std::max(maxx , (int)vis.size());
    if(vis.size() > 50 || 1) {
        for (auto v: vis) {
            fillet[v] += 1.0;
        }
    }
}

int main() {
    std::string file_path = "../data/20431_27e7f151_2.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    MeshVoronoi3d mv3d(mesh);

    const std::vector<easy3d::vec3>& vertices = mv3d.get_vertices();
    const std::vector<easy3d::vec3>& sites = mv3d.get_sites();
    const std::vector<std::vector<int>>& regions = mv3d.get_regions();
    const std::vector<int>& final_pole = mv3d.get_final_pole();
    const std::vector<std::vector<int>> vertices2sites = mv3d.get_vertices2sites();
    int nb_sites = sites.size();
    int nb_vertices = vertices.size();
    auto fillet = mesh->add_face_property<double>("f:fillet", 0);
    easy3d::Viewer viewer("ASD");
    viewer.add_model(mesh);
//    mesh->renderer()->get_triangles_drawable("faces")->set_visible(false);
//    mesh->renderer()->get_lines_drawable("edges")->set_visible(true);
    easy3d::PointsDrawable* pd_vertices = new easy3d::PointsDrawable;
    easy3d::PointsDrawable* pd_sites = new easy3d::PointsDrawable;
    easy3d::LinesDrawable* lDrawable = new easy3d::LinesDrawable;
    std::vector<easy3d::vec3> a;
    std::vector<easy3d::vec3> c;
    std::vector<easy3d::vec3> b;
    Box3 box = mesh->bounding_box();
    int ct = 0;
    for(int i = 0; i < nb_vertices; i++) {
        int x = i;
        double len = (vertices[x] - sites[vertices2sites[x][0]]).norm();
        if(len < 1.5 && box.contains(vertices[x]) &&is_connected(mesh, vertices2sites[x])) {
            domain_region(mesh,  vertices[x], sites, vertices2sites[x][0]);
            a.emplace_back(vertices[x]);
            for(int j = 0; j < vertices2sites[x].size(); j++) {
                int id = vertices2sites[x][j];
                c.emplace_back(sites[id]);
                b.emplace_back(vertices[id]);
                b.emplace_back(sites[id]);
            }
            ct++;
            if(ct >10) break;
        }
    }

    Eigen::VectorXd Z(nb_sites);
    for (int i = 0; i < nb_sites; i++) {
        Z[i] = fillet.vector()[i];
    }
    Eigen::MatrixXd Ct;
    igl::jet(Z, true, Ct);
    auto drawable = mesh->renderer()->get_triangles_drawable("faces");
    std::string color_name = "f:color";
    auto coloring = mesh->face_property<vec3>(color_name, vec3(0, 0, 0));
    for(auto f : mesh->faces()) {
        coloring[f] = easy3d::vec3(Ct(f.idx(), 0),Ct(f.idx(), 1), Ct(f.idx(), 2));
    }
    drawable->set_property_coloring(State::FACE, color_name);
    drawable->update();
//    pd_vertices->update_vertex_buffer(a);
//    pd_vertices->set_impostor_type(easy3d::PointsDrawable::SPHERE);
//    pd_vertices->set_point_size(10);
//    pd_vertices->set_uniform_coloring(easy3d::vec4(0.5,0.5,0.5,1.0));
//    pd_sites->update_vertex_buffer(c);
//    pd_sites->set_impostor_type(easy3d::PointsDrawable::SPHERE);
//    pd_sites->set_point_size(10);
//    lDrawable->update_vertex_buffer(b);
//
//    viewer.add_drawable(pd_vertices);
//    viewer.add_drawable(pd_sites);
//    viewer.add_drawable(lDrawable);
    std::cout << maxx <<std::endl;
    viewer.run();
    return 0;
}
