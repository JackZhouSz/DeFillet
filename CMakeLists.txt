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

add_executable(fillet_seg app/fillet_seg.cpp src/knn4d.cpp src/knn.cpp src/fillet_seg_v8.cpp src/defillet_v2.cpp  src/voro_viewer.cpp src/surface_mesh_segmenter.cpp)
target_include_directories(fillet_seg PRIVATE src)
target_link_libraries(fillet_seg PRIVATE gcp geodesic)
target_link_libraries(fillet_seg PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(fillet_seg PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
target_link_libraries(fillet_seg PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)


#add_executable(read_clr app/read_clr.cpp)
#target_link_libraries(read_clr PUBLIC glfw imgui::imgui easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
#
#
#add_executable(eval app/eval.cpp src/knn.cpp)
#target_include_directories(eval PRIVATE src)
#target_link_libraries(eval PUBLIC easy3d::algo easy3d::util easy3d::fileio  OpenMP::OpenMP_CXX)
#
#
#add_executable(cur app/cur.cpp)
#target_include_directories(cur PRIVATE src)
#target_link_libraries(cur PUBLIC easy3d::algo easy3d::util easy3d::fileio  igl::igl_core OpenMP::OpenMP_CXX)
#
#
#add_executable(graph_cut_for_mesh graph_cut_for_mesh.cpp)
#target_include_directories(graph_cut_for_mesh PRIVATE src)
#target_link_libraries(graph_cut_for_mesh PUBLIC igl::igl_core gcp easy3d::fileio easy3d::core)


set(SEG_SRCS src/decompose.cpp)


add_executable(main_seg app/main_seg.cpp ${SEG_SRCS})
target_include_directories(main_seg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(main_seg PUBLIC Eigen3::Eigen igl::igl_core easy3d::core easy3d::fileio)



add_executable(test_decompose test/test_decompose.cpp ${SEG_SRCS})
target_include_directories(test_decompose PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_decompose PUBLIC Eigen3::Eigen igl::igl_core easy3d::core easy3d::fileio)


#find_package(coin CONFIG REQUIRED)
#
#
#add_executable(readiv app/readIV.cpp)
#target_link_libraries(readiv PRIVATE Coin::Coin)
#add_executable(defillet app/defillet.cpp src/defillet_v2.cpp)
#target_include_directories(defillet PRIVATE src)
#target_link_libraries(defillet PRIVATE  geodesic)
#target_link_libraries(defillet PRIVATE  Eigen3::Eigen igl::igl_core)
#target_link_libraries(defillet PRIVATE glfw imgui::imgui easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui easy3d::fileio)
#add_compile_definitions(SCALAR_COLOR_TEXTURE="${CMAKE_CURRENT_LIST_DIR}/lib-easy3d/resources/colormaps/rainbow.png")

#set(source src/MeshVoronoi3d.cpp src/fillet_seg.cpp
#        src/knn.cpp src/defillet1.cpp
#)
#
#
#add_executable(DeFillet_GUI app/gui.cpp app/my_viewer.cpp app/my_picker.cpp)
#target_include_directories(DeFillet_GUI PRIVATE app)
#target_link_libraries(DeFillet_GUI PRIVATE gcp geodesic)
#target_link_libraries(DeFillet_GUI PRIVATE CGAL::CGAL Eigen3::Eigen)
#target_link_libraries(DeFillet_GUI PRIVATE glfw imgui::imgui easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui easy3d::fileio)
#target_link_libraries(DeFillet_GUI PRIVATE OpenMP::OpenMP_CXX)
#
#
#add_executable(scoring app/scoring.cpp ${source} src/surafce_mesh_segmenter.cpp)
#target_include_directories(scoring PRIVATE src)
#target_link_libraries(scoring PRIVATE gcp geodesic)
#target_link_libraries(scoring PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
#target_link_libraries(scoring PRIVATE easy3d::algo easy3d::util easy3d::fileio)
#target_link_libraries(scoring PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#
#
#add_executable(gc app/gc.cpp ${source})
#target_include_directories(gc PRIVATE src)
#target_link_libraries(gc PRIVATE gcp geodesic)
#target_link_libraries(gc PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
#target_link_libraries(gc PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#target_link_libraries(gc PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#
#
#add_executable(geo app/geo.cpp ${source})
#target_include_directories(geo PRIVATE src)
#target_link_libraries(geo PRIVATE gcp geodesic)
#target_link_libraries(geo PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
#target_link_libraries(geo PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#target_link_libraries(geo PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#
#add_executable(defillet app/defillet1.cpp ${source})
#target_include_directories(defillet PRIVATE src)
#target_link_libraries(defillet PRIVATE gcp geodesic)
#target_link_libraries(defillet PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
#target_link_libraries(defillet PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#target_link_libraries(defillet PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#
#add_executable(decompose app/decompose.cpp src/surafce_mesh_segmenter.cpp)
#target_include_directories(decompose PRIVATE src)
#target_link_libraries(decompose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#target_link_libraries(decompose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
#
#
#add_executable(compose app/compose.cpp src/surafce_mesh_segmenter.cpp)
#target_include_directories(compose PRIVATE src)
#target_link_libraries(compose PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer easy3d::gui)
#target_link_libraries(compose PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)