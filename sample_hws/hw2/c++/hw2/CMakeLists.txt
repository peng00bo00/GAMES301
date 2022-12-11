cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerOptions.cmake")
project(hw2 LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Eigen 3.4.0
include(eigen)
# Libigl
include(libigl)
# sciplot
include(matplot)

# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)

find_package ( OpenMP )

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${sources})

set_target_properties(${PROJECT_NAME} PROPERTIES 
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    OUTPUT_NAME "hw2")

target_link_libraries(${PROJECT_NAME} PRIVATE igl::glfw igl::imgui matplot OpenMP::OpenMP_CXX)
