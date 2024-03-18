cmake_minimum_required(VERSION 3.20)
project(defillet)

set(CMAKE_CXX_STANDARD 17)

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message("asd")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Eigen3
find_package(Eigen3 CONFIG REQUIRED)

# CGAL
find_package(CGAL CONFIG REQUIRED)

#Libigl
find_package(libigl CONFIG REQUIRED)

#CLI11
find_package(CLI11 CONFIG REQUIRED)

# nlohmann_json
find_package(nlohmann_json CONFIG REQUIRED)

add_subdirectory(gcp)
add_subdirectory(geodesic)
add_subdirectory(lib-easy3d)

#add_subdirectory(lib-easy3d)
#add_executable(example1 example1.cpp defillet.cpp knn.cpp voronoi.cpp src/gcp.cpp optimize.cpp)
#target_link_libraries(example1 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
#        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#
#add_executable(example2 example2.cpp defillet.cpp knn.cpp voronoi.cpp src/gcp.cpp optimize.cpp)
#target_link_libraries(example2 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
#        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#
#
add_executable(example3 example3.cpp)
target_link_libraries(example3 PUBLIC easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#
#
add_executable(example4 example4.cpp src/MeshVoronoi3d.cpp)
target_link_libraries(example4 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)


set(DEFILLET_SOURCES src/MeshVoronoi3d.cpp src/fillet_seg.cpp
        src/viewer.cpp src/knn.cpp src/defillet.cpp src/surafce_mesh_segmenter.cpp)


add_executable(example5 example5.cpp ${DEFILLET_SOURCES})
target_include_directories(example5 PRIVATE src)
target_link_libraries(example5 PRIVATE gcp geodesic)
target_link_libraries(example5 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(example5 PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer 3rd_imgui)

add_executable(scoring src/scoring.cpp ${DEFILLET_SOURCES})
target_include_directories(scoring PRIVATE src)
target_link_libraries(scoring PRIVATE gcp geodesic)
target_link_libraries(scoring PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(scoring PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer 3rd_imgui)
target_link_libraries(scoring PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(gc src/gc.cpp ${DEFILLET_SOURCES})
target_include_directories(gc PRIVATE src)
target_link_libraries(gc PRIVATE gcp geodesic)
target_link_libraries(gc PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(gc PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer 3rd_imgui)
target_link_libraries(gc PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)