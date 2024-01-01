
#include<iostream>

#include "common.h"

using namespace std;

int main() {
    Triangulation T;
    std::vector<std::array<double, 3>> v_pos;
    std::vector<std::vector<size_t>> f_ind;
    std::string file_path = "D:\\code\\defillet\\data\\bottle_fillet_remeshed.ply";
    read_ply_mesh(file_path, v_pos, f_ind);
    std::vector<Point> cgal_points;
    my_points_convert_to_cgal_points(v_pos, cgal_points);
    auto start_time = std::chrono::high_resolution_clock::now();
    int nb_sites = cgal_points.size();

    std::map<Vertex_handle, int> mp;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T.insert(cgal_points[i]);
        mp[v] = i;
    }

    std::cout << T.number_of_vertices() <<std::endl;
    std::vector<Point> cgal_voronoi_vertices;
    std::vector<std::vector<size_t>> sites_regions(nb_sites);
    for(auto cell = T.finite_cells_begin(); cell != T.finite_cells_end(); cell++) {
        auto v = T.dual(cell);
        int vid = cgal_voronoi_vertices.size();
        cgal_voronoi_vertices.emplace_back(v);
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp[vh];
            sites_regions[id].emplace_back(vid);
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
    std::vector<Point> insides_voronoi_vertices;
    select_points_inside_bounding_box(cgal_points, cgal_voronoi_vertices, insides_voronoi_vertices);

    std::vector<std::array<double, 3>> my_voronoi_vertices;
    cgal_points_convert_to_my_points(insides_voronoi_vertices, my_voronoi_vertices);
    std::string out_path = "../data/voronoi11.ply";
    write_ply_points(out_path, my_voronoi_vertices);


    return 0;


}
