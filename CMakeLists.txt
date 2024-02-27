cmake_minimum_required(VERSION 3.20)
project(defillet)

set(CMAKE_CXX_STANDARD 14)


# Eigen3
find_package(Eigen3 CONFIG REQUIRED)

# CGAL
find_package(CGAL CONFIG REQUIRED)

#Libigl
find_package(libigl CONFIG REQUIRED)

add_subdirectory(gcp)
add_subdirectory(lib-easy3d)
add_subdirectory(geodesic)

add_executable(example1 example1.cpp defillet.cpp knn.cpp voronoi.cpp gcp.cpp optimize.cpp)
target_link_libraries(example1 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#add_executable(defillet main.cpp knn.cpp voronoi.cpp defillet1.cpp)
#target_link_libraries(defillet PUBLIC CGAL::CGAL easy3d::renderer
#        easy3d::util easy3d::viewer ${OpenCV_LIBS} gcp geodesic Eigen3::Eigen igl::igl_core)
#
#add_executable(test test.cpp knn.cpp voronoi.cpp defillet1.cpp)
#target_link_libraries(test PUBLIC CGAL::CGAL easy3d::renderer
#        easy3d::util easy3d::viewer ${OpenCV_LIBS} gcp geodesic Eigen3::Eigen igl::igl_core)
#
#add_executable(test_knn test_knn.cpp knn.cpp)
#target_link_libraries(test_knn PUBLIC Eigen3::Eigen CGAL::CGAL)
#
#add_executable(test_geodesic test_geodesic.cpp)
#target_link_libraries(test_geodesic PUBLIC geodesic)
#
#add_executable(test_cgal_surface_mesh test_cgal_surface_mesh.cpp)
#target_link_libraries(test_cgal_surface_mesh PUBLIC CGAL::CGAL easy3d::renderer
#        easy3d::util easy3d::viewer)
#
#add_executable(test_cg_solver test_cg_solver.cpp)
#target_link_libraries(test_cg_solver PUBLIC Eigen3::Eigen)
#
#add_executable(test_optimize test_optimize.cpp knn.cpp voronoi.cpp defillet1.cpp)
#target_link_libraries(test_optimize PUBLIC CGAL::CGAL Eigen3::Eigen igl::igl_core easy3d::renderer
#        easy3d::util easy3d::viewer ${OpenCV_LIBS} gcp geodesic)
#
#add_executable(test_KroneckerProduct test_KroneckerProduct.cpp)
#target_link_libraries(test_KroneckerProduct PUBLIC Eigen3::Eigen)
#
#
#add_executable(test_igl_cot_matrix test_igl_cot_matrix.cpp)
#target_link_libraries(test_igl_cot_matrix PUBLIC Eigen3::Eigen igl::igl_core)


