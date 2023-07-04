cmake_minimum_required(VERSION 3.16)
project(gravomg)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# OpenMP
find_package(OpenMP REQUIRED)

# Add your project files
SET(SRC_FILES
    src/multigrid_solver.cpp
    src/utility.cpp
    src/sampling.cpp
    src/sig21/SSP_qslim_optimal_collapse_edge_callbacks.cpp
    src/sig21/cotmatrix_dense.cpp
    src/sig21/mqwf_dense.cpp
    src/sig21/SSP_decimate.cpp
    src/sig21/get_post_faces.cpp
    src/sig21/quasi_conformal_error.cpp
    src/sig21/sort_vec.cpp
    src/sig21/SSP_collapse_edge.cpp
    src/sig21/SSP_vertexRemoval.cpp
    src/sig21/SSP_qslim.cpp
    src/sig21/query_fine_to_coarse.cpp
    src/sig21/remove_vector_element.cpp
    src/sig21/SSP_vertexRemoval_optimal_collapse_edge_callbacks.cpp
    src/sig21/get_prolong.cpp
    src/sig21/mg_precompute.cpp
    src/sig21/mg_precompute_block.cpp
    src/sig21/always_try_never_care.cpp
    src/sig21/joint_lscm.cpp
    src/sig21/intersect_ordered.cpp
    src/sig21/get_collapse_onering_faces.cpp
    src/sig21/remove_unreferenced_lessF.cpp
    src/sig21/vector_mod.cpp
    src/sig21/compute_vertex_quadrics.cpp
    src/sig21/normalize_unit_area.cpp
    src/sig21/SSP_midpoint.cpp
    src/sig21/vector_area_matrix_size.cpp
    src/sig21/remove_row.cpp
    src/sig21/min_quad_with_fixed_mg.cpp
    src/sig21/compute_barycentric.cpp
    src/sig21/mg_VCycle.cpp
)

SET(HEADER_FILES
    include/gravomg/multigrid_solver.h
    include/gravomg/utility.h
    include/gravomg/sampling.h
)

add_library(gravomg ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(gravomg PUBLIC igl::core OpenMP::OpenMP_CXX)
target_include_directories(gravomg PUBLIC include src/sig21)

# Pardiso support
# Find the path to the MKL include directory
# If you have installed MKL in some other location, please add this to the list of paths to search
find_path(MKL_INCLUDE_DIR mkl_pardiso.h PATHS "/usr/include/mkl")

# Check if the MKL include directory was found
if(MKL_INCLUDE_DIR)
  message("Found MKL include directory: ${MKL_INCLUDE_DIR}")
  target_include_directories(gravomg PUBLIC ${MKL_INCLUDE_DIR})
  add_definitions(-DPARDISO_ENABLED=1)
else()
  message("MKL include directory not found, disabling Pardiso.")
endif()
