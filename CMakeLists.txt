cmake_minimum_required(VERSION 3.0.0)
# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif () 

project(loopy VERSION 0.1.0)

include(CTest)
enable_testing()
set (CMAKE_CXX_STANDARD 11)
find_package(Boost COMPONENTS program_options thread REQUIRED)

include_directories(include)

add_executable(loopy 
    src/main.cpp
    src/bn2cn.cpp
    src/directed_graph.cpp
    src/factor.cpp
    src/graph.cpp
    src/generator.cpp
    src/ipe2u.cpp
    src/loopy2u.cpp
    src/merlin.cpp
    src/program_options.cpp
    src/utils.cpp
    src/mmap2u.cpp
    include/algorithm.h
    include/base.h
    include/bayes_net.h
    include/bn2cn.h
    include/credal_net.h
    include/directed_graph.h
    include/enum.h
    include/factor.h
    include/generator.h
    include/graph.h
    include/index.h
    include/indexed_heap.h
    include/ipe2u.h
    include/loopy2u.h
    include/merlin.h
    include/mmap2u.h
    include/program_options.h
    include/set.h
    include/utils.h
    include/variable_set.h
    include/variable.h
    include/vector.h    
)

set(LINK_LIBS dl ${Boost_LIBRARIES})

target_link_libraries(${PROJECT_NAME} ${LINK_LIBS} ${LINK_LIBS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
