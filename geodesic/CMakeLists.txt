cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME geodesic)
project(${PROJECT_NAME})

set(HEADER_FILES
        BaseModel.h
        Chen_Han.h
        DistanceApproach.h
        EdgePoint.h
        ExactDGPMethod.h
        FacePoint.h
        ICH_WindowFiltering.h
        Parameters.h
        Point3D.h
        RichModel.h
        VectorOperations.h
        Xin_Wang.h)

set(SRC_FILES
        BaseModel.cpp
        Chen_Han.cpp
        DistanceApproach.cpp
        ExactDGPMethod.cpp
        ICH_WindowFiltering.cpp
        Point3D.cpp
        RichModel.cpp
        VectorOperations.cpp
        Xin_Wang.cpp)

add_library(${PROJECT_NAME}
          ${HEAD_FILES} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
