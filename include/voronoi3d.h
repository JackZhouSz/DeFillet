//
// Created by xiaowuga on 2025/8/15.
//

#ifndef VORONOI3D_H
#define VORONOI3D_H

#include <easy3d/core/types.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>

namespace DeFillet {

    /**
     * @brief Compute the 3D Voronoi diagram from a set of input points.
     *
     * This function constructs a 3D Delaunay triangulation from the given set of sites,
     * derives the dual Voronoi diagram, and outputs the Voronoi vertices, radii, and connectivity.
     * It also ensures the Voronoi cells are bounded by inserting "infinite" points at a large distance.
     *
     * @param s    Input site coordinates (list of 3D points).
     * @param box  Input bounding box of all site coordinates
     * @param vv   Output Voronoi vertices (coordinates of Voronoi cell vertices).
     * @param vvr  Output radii of each Voronoi vertex (average distance to its generating sites).
     * @param snv  Output mapping: site index → list of Voronoi vertex indices associated with it.
     * @param vns  Output mapping: Voronoi vertex index → list of site indices that generated it.
     *
     */
    void voronoi3d(const std::vector<easy3d::vec3>& s
                   , const easy3d::Box3& box
                   , std::vector<easy3d::vec3>& vv
                   , std::vector<float>& vvr
                   , std::vector<std::vector<int>>& snv
                   , std::vector<std::vector<int>>& vns) {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Delaunay_triangulation_3<K> Triangulation;
        typedef Triangulation::Vertex_handle Vertex_handle;
        typedef K::Point_3 CGAL_Point;

        int nb_sites = s.size();

        Triangulation T;

        // Map each CGAL vertex handle to its original site index
        std::map<Vertex_handle, int> mp_sites;

        // Insert all input sites into the triangulation
        for (int i = 0; i < nb_sites; i++) {
            Vertex_handle v = T.insert(CGAL_Point(s[i].x, s[i].y, s[i].z));
            mp_sites[v] = i;
        }


        // Expand the bounding box by 10× its diagonal length to avoid infinite Voronoi cells
        double diag = box.diagonal_length() * 10;
        double xmin = box.min_coord(0) - diag;
        double xmax = box.max_coord(0) + diag;
        double ymin = box.min_coord(1) - diag;
        double ymax = box.max_coord(1) + diag;
        double zmin = box.min_coord(2) - diag;
        double zmax = box.max_coord(2) + diag;

        // Insert 8 corner points of the enlarged bounding box
        T.insert(CGAL_Point(xmin, ymin, zmax));
        T.insert(CGAL_Point(xmax, ymin, zmax));
        T.insert(CGAL_Point(xmin, ymax, zmax));
        T.insert(CGAL_Point(xmax, ymax, zmax));
        T.insert(CGAL_Point(xmin, ymin, zmin));
        T.insert(CGAL_Point(xmax, ymin, zmin));
        T.insert(CGAL_Point(xmin, ymax, zmin));
        T.insert(CGAL_Point(xmax, ymax, zmin));

        // Prepare output containers
        vv.clear();
        snv.clear();
        vns.clear();
        snv.resize(nb_sites);

        // Iterate over all finite cells in the triangulation
        for (auto cell = T.finite_cells_begin(); cell != T.finite_cells_end(); cell++) {

            // Dual of a Delaunay cell is a Voronoi vertex
            auto v = T.dual(cell);
            int vid = vv.size();

            vec3 pos(v.x(), v.y(), v.z());

            if(!box.contains(pos))
                continue;

            // Store Voronoi vertex position
            vv.emplace_back(pos);

            std::vector<int> tmp; // Sites that generate this Voronoi vertex
            float radius = 0;     // radius

            // For each vertex of the cell, link it to the Voronoi vertex
            for (int i = 0; i < 4; i++) {
                Vertex_handle vh = cell->vertex(i);
                int id = mp_sites[vh];
                snv[id].emplace_back(vid); // Site → Voronoi vertex mapping
                radius += (s[id] - vv.back()).norm(); // Distance for radius calc
                tmp.emplace_back(id); // Voronoi vertex → sites
            }

            // Store average distance as radius
            radius /= 4;
            vvr.emplace_back(radius);

            // Store list of site indices for this Voronoi vertex
            vns.emplace_back(tmp);
        }
    }
}

#endif //VORONOI3D_H
