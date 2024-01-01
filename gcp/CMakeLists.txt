cmake_minimum_required(VERSION 3.12)

if(DEFINED ENV{VCPKG_DIR})
    set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_DIR}/scripts/buildsystems/vcpkg.cmake
  CACHE STRING "Vcpkg toolchain file")
endif()

set(PROJECT_NAME gcp)
project(${PROJECT_NAME})

set(HEADER_FILES
        block.h
        energy.h
        GCoptimization.h
        graph.h
        LinkedBlockList.h
        )

set(SRC_FILES        
        GCoptimization.cpp    
        graph.cpp             
        LinkedBlockList.cpp   
        maxflow.cpp
        )

add_library(${PROJECT_NAME}
        ${HEADER_FILES}
        ${SRC_FILES}
        )

target_include_directories(${PROJECT_NAME} PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)