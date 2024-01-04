cmake_minimum_required(VERSION 3.26)
project(defillet)

set(CMAKE_CXX_STANDARD 14)


# Eigen3
find_package(Eigen3 CONFIG REQUIRED)
#target_link_libraries(main PRIVATE Eigen3::Eigen)

# CGAL
find_package(CGAL CONFIG REQUIRED)
#target_link_libraries(main PRIVATE CGAL::CGAL)


find_package(nlohmann_json CONFIG REQUIRED)

find_package(OpenCV CONFIG REQUIRED)

#set(Easy3D_DIR "D:\\software\\Easy3D-2.5.2\\build")
find_package(Easy3D CONFIG REQUIRED)

find_package(libigl CONFIG REQUIRED)

add_subdirectory(gcp)
add_subdirectory(geodesic)

add_executable(defillet main.cpp knn.cpp voronoi.cpp defillet.cpp)
target_link_libraries(defillet PUBLIC CGAL::CGAL
        nlohmann_json::nlohmann_json easy3d::renderer
        easy3d::util easy3d::viewer ${OpenCV_LIBS} gcp geodesic Eigen3::Eigen igl::igl_core)

add_executable(test test.cpp)
target_link_libraries(test PUBLIC Eigen3::Eigen CGAL::CGAL)

add_executable(test_knn test_knn.cpp knn.cpp)
target_link_libraries(test_knn PUBLIC Eigen3::Eigen CGAL::CGAL)

add_executable(test_geodesic test_geodesic.cpp)
target_link_libraries(test_geodesic PUBLIC geodesic)

add_executable(test_cgal_surface_mesh test_cgal_surface_mesh.cpp)
target_link_libraries(test_cgal_surface_mesh PUBLIC CGAL::CGAL easy3d::renderer
        easy3d::util easy3d::viewer)

add_executable(test_cg_solver test_cg_solver.cpp)
target_link_libraries(test_cg_solver PUBLIC Eigen3::Eigen)


