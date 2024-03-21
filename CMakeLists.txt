cmake_minimum_required(VERSION 3.20)
project(defillet)

set(CMAKE_CXX_STANDARD 17)

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
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

set(DEFILLET_SOURCES src/MeshVoronoi3d.cpp src/fillet_seg.cpp
        src/viewer.cpp src/knn.cpp src/defillet.cpp src/surafce_mesh_segmenter.cpp)


add_executable(gui app/gui.cpp ${DEFILLET_SOURCES})
target_include_directories(gui PRIVATE src)
target_link_libraries(gui PRIVATE gcp geodesic)
target_link_libraries(gui PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(gui PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui 3rd_imgui)

add_executable(scoring app/scoring.cpp ${DEFILLET_SOURCES})
target_include_directories(scoring PRIVATE src)
target_link_libraries(scoring PRIVATE gcp geodesic)
target_link_libraries(scoring PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(scoring PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui 3rd_imgui)
target_link_libraries(scoring PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(gc app/gc.cpp ${DEFILLET_SOURCES})
target_include_directories(gc PRIVATE src)
target_link_libraries(gc PRIVATE gcp geodesic)
target_link_libraries(gc PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(gc PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui 3rd_imgui)
target_link_libraries(gc PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(geo app/geo.cpp ${DEFILLET_SOURCES})
target_include_directories(geo PRIVATE src)
target_link_libraries(geo PRIVATE gcp geodesic)
target_link_libraries(geo PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(geo PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui 3rd_imgui)
target_link_libraries(geo PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(defillet app/defillet.cpp ${DEFILLET_SOURCES})
target_include_directories(defillet PRIVATE src)
target_link_libraries(defillet PRIVATE gcp geodesic)
target_link_libraries(defillet PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(defillet PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui 3rd_imgui)
target_link_libraries(defillet PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(decompose app/decompose.cpp src/surafce_mesh_segmenter.cpp)
target_include_directories(decompose PRIVATE src)
target_link_libraries(decompose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui)
target_link_libraries(decompose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(compose app/compose.cpp src/surafce_mesh_segmenter.cpp)
target_include_directories(compose PRIVATE src)
target_link_libraries(compose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui)
target_link_libraries(compose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)