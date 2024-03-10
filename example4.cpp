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


double compute_score(const easy3d::SurfaceMesh* mesh,  const std::vector<easy3d::vec3>& sites,
                   easy3d::vec3 center,  const std::vector<int>& faces) {
    std::set<easy3d::SurfaceMesh::Face> vis;
    std::queue<easy3d::SurfaceMesh::Face>que;
    std::map<int,int> root;
    for(auto f : faces) {
        easy3d::SurfaceMesh::Face face_handle(f);
        que.push(face_handle);
        root[f] = f;
    }
    std::function<int(int)> get_root = [&](int k)-> int {
        return root[k] == k ? k : (root[k] = get_root(root[k]));
    };
    double R = (sites[faces[0]] - center).norm();
    double eps = R * 0.02;
    double err = 0;
    while(!que.empty()) {
        auto cur = que.front(); que.pop();
        if(vis.find(cur) != vis.end()) {
            continue;
        }
        vis.insert(cur);
        for(auto h : mesh->halfedges(cur)) {
            auto opp_f = mesh->face(mesh->opposite(h));
            if(opp_f.is_valid()) {
                double len = (sites[opp_f.idx()] - center).norm();
                if(len < R + eps) {
                    if(vis.find(opp_f) == vis.end()) {
                        err += fabs(len - R);
                        que.push(opp_f);
                        int r = get_root(cur.idx());
                        root[opp_f.idx()] = r;
                    } else {
                        int r1 = get_root(cur.idx());
                        int r2 = get_root(opp_f.idx());
                        root[r1] = r2;
                    }
                }

            }
        }
    }

    int ct = 0;
    for(auto item : root) {
        root[item.first] = get_root(item.first);
        if(item.first == item.second) {
            ct++;
            if(ct > 1) break;
        }
    }
    if(ct > 1) {
        return 0.0;
    }
    double x = 1.0 * err / vis.size();
    return exp(-10* x);
}




void domain_region(easy3d::SurfaceMesh* mesh, easy3d::vec3 center,
                   const std::vector<easy3d::vec3>& sites,  int st, double score) {
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
                if(abs(axis_len - len) / axis_len < 0.02 && vis.find(opp_f) == vis.end()) {
                    que.push(opp_f);
                }
            }
        }
    }
    maxx = std::max(maxx , (int)vis.size());
    if(1) {
        for (auto v: vis) {
            fillet[v] += score;
        }
    }
}

int main() {
    std::string file_path = "../data/21681_cbbaafd2_3.ply";
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
    double minn = 1e9;
    Box3 box = mesh->bounding_box();

    for(int i = 0; i < nb_vertices; i++) {
        double len = (vertices[i] - sites[vertices2sites[i][0]]).norm();
        if(len < 3.0 && box.contains(vertices[i])) {
            double score = compute_score(mesh, sites, vertices[i], vertices2sites[i]);
            minn = min(score, minn);
            if(score > 1e-3) {
                domain_region(mesh, vertices[i], sites, vertices2sites[i][0], score);
            }
        } else {

        }
    }

    std::cout << "minn=" <<minn <<std::endl;
//    mesh->renderer()->get_triangles_drawable("faces")->set_visible(false);
//    mesh->renderer()->get_lines_drawable("edges")->set_visible(true);
    easy3d::PointsDrawable* pd_vertices = new easy3d::PointsDrawable;
    easy3d::PointsDrawable* pd_sites = new easy3d::PointsDrawable;
    easy3d::LinesDrawable* lDrawable = new easy3d::LinesDrawable;
    std::vector<easy3d::vec3> a;
    std::vector<easy3d::vec3> c;
    std::vector<easy3d::vec3> b;

    int ct = 0;

//    for(int i = 0; i < nb_vertices; i++) {
//        int x = i;
//        double len = (vertices[x] - sites[vertices2sites[x][0]]).norm();
//        if(box.contains(vertices[x]) &&is_connected(mesh, vertices2sites[x])) {
//            domain_region(mesh,  vertices[x], sites, vertices2sites[x][0]);
//            a.emplace_back(vertices[x]);
//            for(int j = 0; j < vertices2sites[x].size(); j++) {
//                int id = vertices2sites[x][j];
//                c.emplace_back(sites[id]);
//                b.emplace_back(vertices[id]);
//                b.emplace_back(sites[id]);
//            }
//            ct++;
////            if(ct >1000) break;
//        }
//    }

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
