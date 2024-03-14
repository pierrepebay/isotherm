cmake_minimum_required(VERSION 3.27)
project(heat2d)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Add VTK build directory
list(APPEND CMAKE_PREFIX_PATH "~/Develop/vtk-build")

# Make VTK optional
find_package(
    VTK QUIET COMPONENTS
    RenderingCore
    IOExodus
    IOParallel
    IOXML
    CommonColor
    CommonCore
    CommonDataModel
    FiltersCore
    FiltersSources
    FiltersGeneral
    InteractionStyle
    RenderingFreeType
    RenderingOpenGL2
    RenderingGL2PSOpenGL2
    RenderingAnnotation
)

if(VTK_FOUND)
    message(STATUS "VTK libraries: ${VTK_LIBRARIES}")
    add_definitions(-DVTK_AVAILABLE)
endif()

add_executable(heat2d heat2d.cpp)
add_executable(heat2d_parallel heat2d_parallel.cpp)

target_link_libraries(heat2d ${MPI_LIBRARIES})
target_link_libraries(heat2d_parallel ${MPI_LIBRARIES})

if(VTK_FOUND)
    target_link_libraries(heat2d ${VTK_LIBRARIES})
    target_link_libraries(heat2d_parallel ${VTK_LIBRARIES})

    vtk_module_autoinit(
        TARGETS heat2d
        MODULES ${VTK_LIBRARIES}
    )

    vtk_module_autoinit(
        TARGETS heat2d_parallel
        MODULES ${VTK_LIBRARIES}
    )
endif()
