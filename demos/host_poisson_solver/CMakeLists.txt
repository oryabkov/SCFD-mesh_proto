
file(COPY test.msh DESTINATION .)

#find_library(GMSH_LIB_DIR Gmsh PATH_SUFFIXES lib DOC "directory where gmsh lib files are located")
#find_path(GMSH_INCLUDE_DIR "Gmsh.h" PATH_SUFFIXES include gmsh include/gmsh DOC "directory where gmsh header files are located")
#if(NOT GMSH_LIB_DIR)
#  message(FATAL_ERROR "gmsh lib not found; stop")
#endif(NOT GMSH_LIB_DIR)
#if(NOT GMSH_INCLUDE_DIR)
#  message(FATAL_ERROR "gmsh headers not found; stop")
#endif(NOT GMSH_INCLUDE_DIR)

set(GMSH_LIB_DIR "${GMSH_INSTALL_ROOT_DIR}/lib")
set(GMSH_INCLUDE_DIR "${GMSH_INSTALL_ROOT_DIR}/include")

include_directories(${GMSH_INCLUDE_DIR})
#message(${GMSH_LIB_DIR})
link_directories(${GMSH_LIB_DIR})

add_executable(host_poisson_solver.bin host_poisson_solver.cpp)
target_link_libraries(host_poisson_solver.bin gmsh)

