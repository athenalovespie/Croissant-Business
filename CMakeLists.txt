cmake_minimum_required(VERSION 3.20)

project(EigenHelloWorld VERSION 1.0 LANGUAGES CXX)

# Use a modern, beginner-friendly version of C++.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Download Eigen automatically the first time the project is configured.
include(FetchContent)
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Eigen)

add_executable(hello_eigen src/main.cpp)
target_link_libraries(hello_eigen PRIVATE Eigen3::Eigen)

