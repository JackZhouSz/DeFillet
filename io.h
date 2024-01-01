//
// Created by xiaowuga on 12/17/2023.
//

#ifndef DEFILLET_IO_H
#define DEFILLET_IO_H
#include "common.h"

inline void read_ply_mesh(const std::string& file_path, std::vector<std::array<double, 3>>& v_pos
        , std::vector<std::vector<size_t>>& f_ind) {
    happly::PLYData ply_in(file_path);
    v_pos = ply_in.getVertexPositions();
    f_ind = ply_in.getFaceIndices<size_t>();
}

inline void write_ply_points(const std::string& file_path, std::vector<std::array<double, 3>>& points) {
    happly::PLYData plyOut;

    plyOut.addVertexPositions(points);

    plyOut.write(file_path, happly::DataFormat::ASCII);
}

#endif //DEFILLET_IO_H
