
file(COPY test_box3d.msh DESTINATION .)
file(COPY test_box3d.geo DESTINATION .)
file(COPY test_box3d_small_mesh.msh DESTINATION .)
file(COPY test_box3d_small_mesh.geo DESTINATION .)
file(COPY test_box3d_period_small_mesh.msh DESTINATION .)
file(COPY test_box3d_period_small_mesh.geo DESTINATION .)

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

add_executable(test_gmsh_mesh_wrap.bin test_gmsh_mesh_wrap.cpp)
target_link_libraries(test_gmsh_mesh_wrap.bin gmsh)
target_link_libraries(test_gmsh_mesh_wrap.bin gtest_main)

add_executable(test_host_mesh_gmsh_wrap.bin test_host_mesh_gmsh_wrap.cpp)
target_link_libraries(test_host_mesh_gmsh_wrap.bin gmsh)
target_link_libraries(test_host_mesh_gmsh_wrap.bin gtest_main)

add_executable(test_gmsh_cpp_api_geom.bin test_gmsh_cpp_api_geom.cpp)
target_link_libraries(test_gmsh_cpp_api_geom.bin gmsh)
target_link_libraries(test_gmsh_cpp_api_geom.bin gtest_main)

add_test(NAME test_gmsh_mesh_wrap COMMAND test_gmsh_mesh_wrap.bin)
add_test(NAME test_host_mesh_gmsh_wrap COMMAND test_host_mesh_gmsh_wrap.bin)
add_test(NAME test_gmsh_cpp_api_geom COMMAND test_gmsh_cpp_api_geom.bin)
