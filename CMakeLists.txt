cmake_minimum_required(VERSION 3.20)
project(defillet)

set(CMAKE_CXX_STANDARD 17)



# Eigen3
find_package(Eigen3 CONFIG REQUIRED)

# CGAL
find_package(CGAL CONFIG REQUIRED)

#Libigl
find_package(libigl CONFIG REQUIRED)


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
#add_executable(example3 example3.cpp defillet.cpp knn.cpp voronoi.cpp src/gcp.cpp optimize.cpp)
#target_link_libraries(example3 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
#        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
#
#
#add_executable(example4 example4.cpp src/MeshVoronoi3d.cpp src/fillet_seg.cpp)
#target_link_libraries(example4 PUBLIC gcp geodesic CGAL::CGAL Eigen3::Eigen
#        igl::igl_core easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)


set(DEFILLET_SOURCES src/MeshVoronoi3d.cpp src/fillet_seg.cpp src/gcp.cpp src/viewer.cpp)
add_executable(example5 example5.cpp ${DEFILLET_SOURCES})
target_include_directories(example5 PRIVATE src)
target_link_libraries(example5 PRIVATE gcp geodesic)
target_link_libraries(example5 PRIVATE CGAL::CGAL Eigen3::Eigen igl::igl_core)
target_link_libraries(example5 PRIVATE easy3d::algo easy3d::util easy3d::viewer easy3d::renderer 3rd_imgui)


#add_executable(segment segment.cpp)
#target_link_libraries(segment PUBLIC gcp easy3d::algo easy3d::util easy3d::viewer easy3d::renderer)
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


