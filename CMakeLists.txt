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

add_subdirectory(gcp)
add_subdirectory(geodesic)





set(SEG_SRCS src/fillet_detector.cpp src/utils.cpp src/knn4d.cpp)


add_executable(test_seg test/test_seg.cpp ${SEG_SRCS})
target_include_directories(test_seg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_seg PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(test_seg PUBLIC easy3d::algo easy3d::core easy3d::fileio)
target_link_libraries(test_seg PUBLIC OpenMP::OpenMP_CXX CLI11::CLI11 nlohmann_json::nlohmann_json)
target_link_libraries(test_seg PUBLIC gcp)
