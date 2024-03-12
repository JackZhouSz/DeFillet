//
// Created by 13900K on 2024/3/12.
//



#include <easy3d/util/logging.h>
#include "viewer.h"

int main() {
    easy3d::logging::initialize();
//    const std::string file_name = resource::directory() + "/data/easy3d.ply";
    easy3d::ViewerImGui viewer("DeFillet");


    viewer.resize(800, 600);
    return viewer.run();
}