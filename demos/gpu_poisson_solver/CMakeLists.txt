# Copyright © 2016,2017 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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

cmake_minimum_required (VERSION 2.6)
project (unstructured_mesh_framework/gpu_poisson_solver_sample)

set (GPU_POISSON_SOLVER_VERSION_MAJOR 0)
set (GPU_POISSON_SOLVER_VERSION_MINOR 1)

#set (SOURCEFILES gpu_poisson_solver.cu  )
set (HEADERFILES gpu_poisson_solver_config.h)

include(FindCUDA)

include_directories(../common)
include_directories(common-local)

function(generate_target add_tet add_double)
  set(CURRENT_TARGET_NAME "gpu_poisson_solver")
  set(CURRENT_TARGET_FLAGS "")
  if(add_tet)
    set(CURRENT_TARGET_NAME "${CURRENT_TARGET_NAME}_tet")
    set(CURRENT_TARGET_FLAGS "${CURRENT_TARGET_FLAGS} -DGPU_POISSON_SOLVER_MAX_FACES_N=4 -DGPU_POISSON_SOLVER_MAX_VERT_N=4")
  else(add_tet)
    set(CURRENT_TARGET_FLAGS "${CURRENT_TARGET_FLAGS} -DGPU_POISSON_SOLVER_MAX_FACES_N=6 -DGPU_POISSON_SOLVER_MAX_VERT_N=8")
  endif(add_tet)
  if(add_double)
    set(CURRENT_TARGET_NAME "${CURRENT_TARGET_NAME}_double")
    set(CURRENT_TARGET_FLAGS "${CURRENT_TARGET_FLAGS} -DGPU_POISSON_SOLVER_SCALAR_TYPE=double")
  endif(add_double)

  cuda_add_executable(${CURRENT_TARGET_NAME} gpu_poisson_solver.cu
                      OPTIONS ${CURRENT_TARGET_FLAGS})
  if(NOT (CURRENT_TARGET_FLAGS STREQUAL ""))
    set_target_properties(${CURRENT_TARGET_NAME}
                          PROPERTIES COMPILE_FLAGS ${CURRENT_TARGET_FLAGS}
  			    )
  endif()

endfunction(generate_target)

generate_target(0 0)
generate_target(1 0)
generate_target(0 1)
generate_target(1 1)
