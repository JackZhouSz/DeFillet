cmake_minimum_required(VERSION 3.20)
project(DeFillet)

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

#add_executable(fillet_seg app/fillet_seg.cpp src/fillet_seg_v8.cpp src/defillet_v2.cpp  src/voro_viewer.cpp src/surface_mesh_segmenter.cpp)
#target_include_directories(fillet_seg PRIVATE src)
#target_link_libraries(fillet_seg PRIVATE gcp geodesic)
#target_link_libraries(fillet_seg PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
#target_link_libraries(fillet_seg PRIVATE easy3d::gui easy3d::algo easy3d::util easy3d::fileio easy3d::viewer easy3d::renderer easy3d::viewer)
#target_link_libraries(fillet_seg PRIVATE CLI11::CLI11 nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)




set(SEG_SRCS src/decompose.cpp src/preliminary_fillet_selector.cpp)


#add_executable(main_seg app/main_seg.cpp ${SEG_SRCS})
#target_include_directories(main_seg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(main_seg PUBLIC Eigen3::Eigen igl::igl_core easy3d::core easy3d::fileio)




add_executable(test_seg test/test_seg.cpp src/fillet_detector.cpp src/utils.cpp src/knn4d.cpp src/axis_optimizer.cpp)
target_include_directories(test_seg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_seg PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test_seg PUBLIC easy3d::algo easy3d::core easy3d::fileio)
target_link_libraries(test_seg PUBLIC OpenMP::OpenMP_CXX)
