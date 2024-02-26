cmake_minimum_required(VERSION 3.27)
project(heat2d)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Add VTK build directory
list(APPEND CMAKE_PREFIX_PATH "~/Develop/vtk-build")

find_package(
    VTK REQUIRED COMPONENTS
    RenderingCore
    IOExodus
    IOParallel
    IOXML
    CommonColor
    CommonCore
    CommonDataModel
    FiltersSources
    FiltersGeneral
    InteractionStyle
    RenderingFreeType
    RenderingOpenGL2
    RenderingGL2PSOpenGL2
    RenderingAnnotation
)

message(STATUS "VTK libraries: ${VTK_LIBRARIES}")

add_executable(heat2d heat2d.cpp)

target_link_libraries(heat2d ${MPI_LIBRARIES})
target_link_libraries(heat2d ${VTK_LIBRARIES})
