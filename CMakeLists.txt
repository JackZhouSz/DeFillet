cmake_minimum_required(VERSION 3.20)
project(defillet)

set(CMAKE_CXX_STANDARD 17)

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_package(Easy3d CONFIG REQUIRED)

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

#imgui
find_package(imgui CONFIG REQUIRED)

#glfw
find_package(glfw3 CONFIG REQUIRED)

add_subdirectory(gcp)
add_subdirectory(geodesic)

#add_compile_definitions(SCALAR_COLOR_TEXTURE="${CMAKE_CURRENT_LIST_DIR}/lib-easy3d/resources/colormaps/rainbow.png")

set(source src/MeshVoronoi3d.cpp src/fillet_seg.cpp
        src/knn.cpp src/defillet.cpp
)


add_executable(DeFillet_GUI app/gui.cpp app/my_viewer.cpp app/my_picker.cpp)
target_include_directories(DeFillet_GUI PRIVATE app)
target_link_libraries(DeFillet_GUI PRIVATE gcp geodesic)
target_link_libraries(DeFillet_GUI PRIVATE CGAL::CGAL Eigen3::Eigen)
target_link_libraries(DeFillet_GUI PRIVATE glfw imgui::imgui easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui easy3d::fileio)
target_link_libraries(DeFillet_GUI PRIVATE OpenMP::OpenMP_CXX)


add_executable(scoring app/scoring.cpp ${source} src/surafce_mesh_segmenter.cpp)
target_include_directories(scoring PRIVATE src)
target_link_libraries(scoring PRIVATE gcp geodesic)
target_link_libraries(scoring PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(scoring PRIVATE easy3d::algo easy3d::util easy3d::fileio)
target_link_libraries(scoring PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(gc app/gc.cpp ${source})
target_include_directories(gc PRIVATE src)
target_link_libraries(gc PRIVATE gcp geodesic)
target_link_libraries(gc PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(gc PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
target_link_libraries(gc PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(geo app/geo.cpp ${source})
target_include_directories(geo PRIVATE src)
target_link_libraries(geo PRIVATE gcp geodesic)
target_link_libraries(geo PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(geo PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
target_link_libraries(geo PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(defillet app/defillet.cpp ${source})
target_include_directories(defillet PRIVATE src)
target_link_libraries(defillet PRIVATE gcp geodesic)
target_link_libraries(defillet PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(defillet PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
target_link_libraries(defillet PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(decompose app/decompose.cpp src/surafce_mesh_segmenter.cpp)
target_include_directories(decompose PRIVATE src)
target_link_libraries(decompose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
target_link_libraries(decompose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(compose app/compose.cpp src/surafce_mesh_segmenter.cpp)
target_include_directories(compose PRIVATE src)
target_link_libraries(compose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui)
target_link_libraries(compose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(scoring2 app/test.cpp ${source} src/surafce_mesh_segmenter.cpp src/fillet_seg_v2.cpp)
target_include_directories(scoring2 PRIVATE src)
target_link_libraries(scoring2 PRIVATE gcp geodesic)
target_link_libraries(scoring2 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(scoring2 PRIVATE easy3d::algo easy3d::util easy3d::fileio)
target_link_libraries(scoring2 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(eigen_value app/eigen_value.cpp ${source} src/surafce_mesh_segmenter.cpp src/fillet_seg_v2.cpp)
target_include_directories(scoring2 PRIVATE src)
target_link_libraries(eigen_value PRIVATE gcp geodesic)
target_link_libraries(eigen_value PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(eigen_value PRIVATE easy3d::algo easy3d::util easy3d::fileio)
target_link_libraries(eigen_value PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(test2 app/test2.cpp)
target_include_directories(test2 PRIVATE src)
#target_link_libraries(test2 PRIVATE gcp geodesic)
target_link_libraries(test2 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test2 PRIVATE easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer)
target_link_libraries(test2 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(test3 app/test3.cpp src/fillet_seg_v3.cpp src/knn.cpp)
target_include_directories(test3 PRIVATE src)
target_link_libraries(test2 PRIVATE gcp geodesic)
target_link_libraries(test3 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test3 PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer)
target_link_libraries(test3 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


find_package(assimp CONFIG REQUIRED)
add_executable(test_assimp app/test_assimp.cpp)
target_link_libraries(test_assimp PRIVATE easy3d::algo easy3d::util easy3d::fileio easy3d::viewer)
target_link_libraries(test_assimp PRIVATE assimp::assimp)


add_executable(vor_cor_sites app/vor_cor_sites.cpp ${source} src/surafce_mesh_segmenter.cpp src/fillet_seg_v3.cpp)
target_include_directories(vor_cor_sites PRIVATE src)
target_link_libraries(vor_cor_sites PRIVATE gcp geodesic)
target_link_libraries(vor_cor_sites PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(vor_cor_sites PRIVATE glfw easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer)
target_link_libraries(vor_cor_sites PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(test4 app/test4.cpp)
target_include_directories(test4 PRIVATE src)
target_link_libraries(test4 PRIVATE gcp geodesic)
target_link_libraries(test4 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test4 PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer)
target_link_libraries(test4 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(test5 app/test5.cpp)
target_include_directories(test5 PRIVATE src)
target_link_libraries(test5 PRIVATE gcp geodesic)
target_link_libraries(test5 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test5 PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
target_link_libraries(test5 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(test6 app/test6.cpp src/knn.cpp src/fillet_seg_v4.cpp src/knn4d.cpp)
target_include_directories(test6 PRIVATE src)
target_include_directories(test6 PRIVATE dkm)
#target_link_libraries(test5 PRIVATE gcp geodesic)
target_link_libraries(test6 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test6 PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
target_link_libraries(test6 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


add_executable(test7 app/test7.cpp src/knn.cpp src/fillet_seg_v5.cpp src/mean_shift.cpp)
target_include_directories(test7 PRIVATE src)
#target_include_directories(test7 PRIVATE dkm)
target_link_libraries(test7 PRIVATE gcp geodesic)
target_link_libraries(test7 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test7 PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
target_link_libraries(test7 PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)