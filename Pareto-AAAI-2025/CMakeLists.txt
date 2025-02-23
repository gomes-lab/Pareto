cmake_minimum_required(VERSION 3.14)

project(HydroDams)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost REQUIRED)

find_package(OpenMP REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(hclust-cpp-master)

set(GCC_OPT_FLAGS "-O3 -ftree-vectorize -ffast-math -pthread")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_OPT_FLAGS}")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
endif()

add_executable(HydroDams HydroDams.cpp ParetoDP.cpp Network.cpp Solution.cpp Node.h Consts.h Config.cpp IOHandler.cpp DivideAndConquer.cpp DPOriginal.cpp DPBase.cpp DPTransforms.cpp Transform.cpp DPBruteForce.cpp SolutionCompressor.cpp hclust-cpp-master/fastcluster.cpp)
add_executable(HydroDamsUnion HydroDamsUnion.cpp ParetoUnion.cpp Solution.cpp Node.h Consts.h Config.cpp DivideAndConquer.cpp hclust-cpp-master/fastcluster.cpp)
add_executable(HydroDamsBruteForce HydroDamsBruteForce.cpp ParetoDP.cpp Network.cpp Solution.cpp Node.h Consts.h Config.cpp IOHandler.cpp DivideAndConquer.cpp DPOriginal.cpp DPBase.cpp DPTransforms.cpp Transform.cpp DPBruteForce.cpp SolutionCompressor.cpp hclust-cpp-master/fastcluster.cpp)

target_link_libraries(HydroDams PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(HydroDamsUnion PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(HydroDamsBruteForce PRIVATE OpenMP::OpenMP_CXX)
