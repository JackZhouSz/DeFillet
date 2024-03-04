//
// Created by 13900K on 2024/3/4.
//

#include <iostream>
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <queue>
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/renderer/drawable_triangles.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/texture_manager.h>


using namespace std;

int main() {
    std::string file_path = "../data/20440_27177360_1.ply";
    easy3d::SurfaceMesh* mesh = easy3d::SurfaceMeshIO::load(file_path);
    auto face_normal = mesh->face_property<easy3d::vec3>("f:normal");
    auto face_label = mesh->face_property<int>("f:label", -1);
    for(auto f : mesh->faces()) {
        face_normal[f] = mesh->compute_face_normal(f);
    }
    int ct = 0;
    for(auto f : mesh->faces()) {
        if(face_label[f] == -1) {
            std::queue<easy3d::SurfaceMesh::Face> que;
            que.push(f);
            while(!que.empty()) {
                auto cur = que.front(); que.pop();
                if(face_label[cur] != -1)
                    continue;
                face_label[cur] = ct;
                for (auto h: mesh->halfedges(cur)) {
                    auto opp = mesh->face(mesh->opposite(h));
                    if (opp.is_valid() && face_label[opp] == -1) {
                        if (easy3d::dot(face_normal[cur], face_normal[opp]) > 0.86) {
                            que.push(opp);
                        }
                    }
                }
            }
            ct++;
        }
    }


    for(auto f : mesh->faces()) {
        if(face_label[f] != 0) {
            mesh->delete_face(f);
        }
    }
    for(auto v : mesh->vertices()) {
        int num = 0;
        for(auto f : mesh->faces(v)) {
            num++;
        }
        if(num == 0) {
            mesh->delete_vertex(v);
        }
    }

    std::vector<std::vector<size_t>> sites;
    std::vector<std::vector<size_t>> faces;
    faces.resize(ct);
    sites.resize(ct);
    for(auto f : mesh->faces()) {
        int id = face_label[f];
        faces[id].emplace_back(f.idx());
    }
    for(int i = 0; i < ct; i++) {
        int num = faces[i].size();
        std::set<size_t> s;
        for(int j = 0; j < num; j++) {
            easy3d::SurfaceMesh::Face f(faces[i][j]);
            for(auto v : mesh->vertices(f)) {
                s.insert(v.idx());
            }
        }
        sites[i].insert(sites[i].end(), s.begin(), s.end());
    }

    easy3d::SurfaceMesh* mesh0 = new easy3d::SurfaceMesh;
    std::map<int, int> mp;
    for(int i = 0; i < sites[0].size(); i++) {
        easy3d::SurfaceMesh::Vertex v(sites[0][i]);
        easy3d::vec3 pos = mesh->position(v);
        mesh0->add_vertex(pos);
        mp[sites[0][i]] = i;
    }
    for(int i = 0; i < faces[0].size(); i++) {
        easy3d::SurfaceMesh::Face f(faces[0][i]);
        std::vector<int>tmp;
        for(auto v : mesh->vertices(f)) {
            tmp.emplace_back(mp[v.idx()]);
        }
        mesh0->add_triangle(easy3d::SurfaceMesh::Vertex(tmp[0]),
                            easy3d::SurfaceMesh::Vertex(tmp[1]),
                            easy3d::SurfaceMesh::Vertex(tmp[2]));
    }
    easy3d::Viewer viewer("binary_segmentation");
    viewer.add_model(mesh0);
//    auto drawable = mesh->renderer()->get_triangles_drawable("faces");
//    std::string color_name = "f:color";
//    auto coloring = mesh->face_property<easy3d::vec3>(color_name, easy3d::vec3(0, 0, 0));
//    std::vector<easy3d::vec3> color_table(ct);
//    for(int i = 0; i < ct; i++) {
//        color_table[i] = easy3d::random_color();
//    }
//    for(auto f : mesh->faces()) {
//        coloring[f] = color_table[face_label[f]];
//    }
//    drawable->set_property_coloring(easy3d::State::FACE, color_name);
//    drawable->update();
//    viewer.update();
    viewer.run();
//

    return 0;
}