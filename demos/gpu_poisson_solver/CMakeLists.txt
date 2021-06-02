# Copyright © 2016-2021 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

# This file is part of SimpleCFD.

# SimpleCFD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 only of the License.

# SimpleCFD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SimpleCFD.  If not, see <http://www.gnu.org/licenses/>.

set(GMSH_LIB_DIR "${GMSH_INSTALL_ROOT_DIR}/lib")
set(GMSH_INCLUDE_DIR "${GMSH_INSTALL_ROOT_DIR}/include")

include_directories(${GMSH_INCLUDE_DIR})
#message(${GMSH_LIB_DIR})
link_directories(${GMSH_LIB_DIR})

function(generate_target add_double)
  set(CURRENT_TARGET_NAME "gpu_poisson_solver")
  set(CURRENT_TARGET_FLAGS "")
  if(add_double)
    set(CURRENT_TARGET_NAME "${CURRENT_TARGET_NAME}_double")
    set(CURRENT_TARGET_FLAGS "${CURRENT_TARGET_FLAGS} -DGPU_POISSON_SOLVER_SCALAR_TYPE=double")
  endif(add_double)

  set(CURRENT_TARGET_NAME "${CURRENT_TARGET_NAME}.bin")

  add_executable(${CURRENT_TARGET_NAME} gpu_poisson_solver.cu)
  if(NOT (CURRENT_TARGET_FLAGS STREQUAL ""))
    set_target_properties(${CURRENT_TARGET_NAME}
                          PROPERTIES COMPILE_FLAGS ${CURRENT_TARGET_FLAGS}
  			    )
  endif()

  target_link_libraries(${CURRENT_TARGET_NAME} gmsh)

endfunction(generate_target)

generate_target(0)
generate_target(1)
