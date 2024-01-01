//
// Created by xiaowuga on 12/29/2023.
//
#include "common.h"
#include "io.h"

int main() {
    std::string file_path = "D:\\code\\defillet\\data\\bottle_fillet_remeshed.ply";
    std::vector<std::array<double, 3>> v_pos;
    std::vector<std::vector<size_t>> f_ind;
    read_ply_mesh(file_path,v_pos, f_ind);

    Surface_mesh cgal_mesh;
    int nb_v = v_pos.size();
    for(int i = 0; i < nb_v; i++) {
        cgal_mesh.add_vertex(Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]));
    }

    int nb_f = f_ind.size();
    for(int i = 0; i < nb_f; i++) {
        if(i < 10) {
            std::cout << f_ind[i][0] << ' ' << f_ind[i][1] << ' ' << f_ind[i][2] << std::endl;
        }
        cgal_mesh.add_face(Vertex_index(f_ind[i][0])
                           ,Vertex_index(f_ind[i][1]), Vertex_index(f_ind[i][2]));
    }
    std::cout << std::endl;
    int ct = 0;
    for (auto f : cgal_mesh.faces()) {
        if(ct ++ < 10) {
            for (auto fv : cgal_mesh.vertices_around_face(cgal_mesh.halfedge(f))) {
                std::cout << fv.idx() << " ";
            }
            std::cout << std::endl;
        }
    }
    return 0;
}