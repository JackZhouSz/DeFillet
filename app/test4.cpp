//
// Created by xiaowuga on 2024/10/31.
//

#include<easy3d/core/point_cloud.h>
#include<easy3d/fileio/point_cloud_io.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include <easy3d/util/file_system.h>
#include <easy3d/viewer/viewer.h>


void voronoi3d(const std::vector<easy3d::vec3>& sites
               , std::vector<easy3d::vec3>& vor_vertices
               , std::vector<std::vector<int>>& sites_nearby_vertices
               , std::vector<std::vector<int>>& vertices_nearby_sites) {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
    typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
    typedef Triangulation::Vertex_handle  Vertex_handle;
    typedef Triangulation::Cell_handle  Cell_handle;
    typedef K::Point_3          CGAL_Point;
    typedef K::Vector_3         CGAL_Vector_3;
    typedef K::Aff_transformation_3 Transformation;
    typedef CGAL::Surface_mesh<CGAL_Point> Surface_mesh;
    typedef Surface_mesh::Vertex_index Vertex_index;
    typedef Surface_mesh::Face_index Face_index;
    typedef Surface_mesh::Edge_index Edge_index;
    typedef Surface_mesh::Halfedge_index Halfedge_index;

    int nb_sites = sites.size();
    easy3d::Box3 box;

    for(auto p : sites) {box.grow(p);}
    Triangulation T;

    std::map<Vertex_handle, int> mp_sites;
    for(int i = 0; i < nb_sites; i++) {
        Vertex_handle v = T.insert(CGAL_Point(sites[i].x, sites[i].y, sites[i].z ));
        mp_sites[v] = i;
    }

    double diag = box.diagonal_length() * 0.1;
    double xmin = box.min_coord(0) - diag;
    double xmax = box.max_coord(0) + diag;
    double ymin = box.min_coord(1) - diag;
    double ymax = box.max_coord(1) + diag;
    double zmin = box.min_coord(2) - diag;
    double zmax = box.max_coord(2) + diag;
    T.insert(CGAL_Point(xmin, ymin, zmax));
    T.insert(CGAL_Point(xmax, ymin, zmax));
    T.insert(CGAL_Point(xmin, ymax, zmax));
    T.insert(CGAL_Point(xmax, ymax, zmax));
    T.insert(CGAL_Point(xmin, ymin, zmin));
    T.insert(CGAL_Point(xmax, ymin, zmin));
    T.insert(CGAL_Point(xmin, ymax, zmin));
    T.insert(CGAL_Point(xmax, ymax, zmin));

    vor_vertices.clear(); sites_nearby_vertices.clear(); vertices_nearby_sites.clear();
    sites_nearby_vertices.resize(nb_sites);
    for(auto cell = T.finite_cells_begin(); cell != T.finite_cells_end(); cell++) {
        auto v = T.dual(cell);
        int vid = vor_vertices.size();
        vor_vertices.emplace_back(easy3d::vec3(v.x(), v.y(), v.z()));
        std::vector<int> tmp;
        for(int i = 0; i < 4; i++) {
            Vertex_handle vh = cell->vertex(i);
            int id = mp_sites[vh];
            sites_nearby_vertices[id].emplace_back(vid);
            tmp.emplace_back(id);
        }
        vertices_nearby_sites.emplace_back(tmp);
    }

}
int main() {
    std::ofstream outfile("../res2.txt");
    std::string dir = "D:\\2_1\\";
    std::vector<std::string> files;
    easy3d::file_system::get_files(dir, files, false);
    std::vector<float> candidate_r = {7.5, 10, 12.5, 15};
    int kk = 10;
    double tot_err = 0;
    int dsa = 0;
    for(auto file : files) {
        dsa++;
        easy3d::PointCloud* cloud = easy3d::PointCloudIO::load(dir  + file);
        std::vector<easy3d::vec3> sites = cloud->points();
        std::vector<easy3d::vec3> vor_vertices;
        std::vector<std::vector<int>> sites_nearby_vertices;
        std::vector<std::vector<int>> vertices_nearby_sites;
        voronoi3d(sites, vor_vertices, sites_nearby_vertices, vertices_nearby_sites);
        double thr = cloud->bounding_box().diagonal_length() * 0.01;
        double minn = 1e9;
        int idx = 0;
        for(size_t i = 0; i < vor_vertices.size(); i++) {
            std::vector<double> tmp;
            bool flag = false;
            for(size_t j = 0; j < sites.size(); j++) {
                tmp.emplace_back((sites[j] - vor_vertices[i]).norm());
                if(j > 0 && tmp.back() - tmp[0] > thr) {
                    flag = true; break;
                }
            }
            if(flag) {
                continue;
            }
            double avg = 0;
            double var = 0;
            for(size_t j = 0; j < tmp.size(); j++) {
                avg += tmp[j];
            }
            avg /= tmp.size();
            for(size_t j = 0; j < tmp.size(); j++) {
                double a = tmp[j] - avg;
                var += a * a;
            }
            var /= tmp.size();
            if(var < minn) {
                minn = var;
                idx = i;
            }
        }

        float est_r = 0; minn = 1e9;
        for(int i = 0; i < sites.size(); i++) {
            est_r += (sites[i] - vor_vertices[idx]).norm();
        }
        est_r /= sites.size();
        float med = est_r;
        for(int i = 0; i < candidate_r.size(); i++) {
            double err = fabs(est_r - candidate_r[i]);
            if(err < minn) {
                minn = err;
                med = candidate_r[i];
            }
        }
        std::vector<std::pair<double,int>> ss;
        for(size_t i = 0; i < vor_vertices.size(); i++) {
            std::vector<double> tmp;
            bool flag = false;
            float var = 0.0;
            for(size_t j = 0; j < sites.size(); j++) {
                float val = fabs((sites[j] - vor_vertices[i]).norm() - med);
                // if( val > 1) {
                //     flag = true; break;
                // }
                var += val * val;
            }
            if(flag) {
                continue;
            }
            var /= sites.size();
            ss.emplace_back(std::make_pair(var, i));
        }
        std::cout << ss.size() << ' ' << kk << std::endl;
        std::sort(ss.begin(), ss.end());
        easy3d::vec3 center(0,0,0);
        for(int i = 0 ; i < kk; i++) {
            center += vor_vertices[ss[i].second];
        }
        center /= kk;
        est_r = 0;
        for(int i = 0; i < sites.size(); i++) {
            est_r += (sites[i] - center).norm();
        }
        est_r /= sites.size();
        tot_err += fabs(est_r - med);
        // std::cout << file << "   " << center.x <<  "   " << center.y << "   " << center.z << "   " << est_r << std::endl;
        // break;
    }
    tot_err /= dsa;
    std::cout << tot_err << std::endl;
    // outfile.close();
    // double avg_r = 0;
    // for(int i = 0; i < sites.size(); i++) {
    //     avg_r += (sites[i] - vor_vertices[idx]).norm();
    // }
    // avg_r /= sites.size();
    // std::cout << avg_r <<std::endl;

    return 0;
}

