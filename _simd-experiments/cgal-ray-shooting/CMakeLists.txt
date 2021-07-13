cmake_minimum_required(VERSION 3.13)
project(ray-bbox-intersections-fanout4)
set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED COMPONENTS Core)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native")

add_executable(ray_shooting
        ray_shooting.cpp)

target_link_libraries(ray_shooting
        CGAL::CGAL CGAL::CGAL_Core)