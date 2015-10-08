cmake_minimum_required(VERSION 2.8)

project(elwave-2d-cg CXX)

set(MFEM_DIR "/u/artemyev/projects/mfem/")
include_directories("${MFEM_DIR}")
link_directories("${MFEM_DIR}/build")

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/*.hpp")

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} mfem rt)

