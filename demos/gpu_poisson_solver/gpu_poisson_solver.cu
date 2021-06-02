// Copyright © 2016,2017 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

// This file is part of SimpleCFD.

// SimpleCFD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 only of the License.

// SimpleCFD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with SimpleCFD.  If not, see <http://www.gnu.org/licenses/>.

//TODO include cuda explicitly

#include "gpu_poisson_solver_config.h"
#include <string>
#include <cmath>
#include <scfd/utils/log_std.h>
#include <scfd/utils/main_try_catch_macro.h>
#include <scfd/utils/init_cuda.h>
#include <scfd/utils/constant_data.h>
#include <scfd/static_vec/vec.h>
#include <scfd/memory/cuda.h>
#include <scfd/communication/linear_partitioner.h>
//TODO add serial_map
#include <scfd/communication/serial_map.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>
#include <scfd/mesh/host_mesh.h>
#include <scfd/mesh/device_mesh.h>
#include <scfd/mesh/device_mesh_impl.h>
#include "gmsh_pos_output.h"

using real = GPU_POISSON_SOLVER_SCALAR_TYPE;
using ordinal = int;
static const ordinal dim = 3;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,dim,ordinal>;
using host_mesh_t = scfd::mesh::host_mesh<gmsh_wrap_t>;
using vec_t = scfd::static_vec::vec<real,dim>;
using log_t = scfd::utils::log_std;
using mem_t = scfd::memory::cuda_device;
using device_mesh = scfd::mesh::device_mesh<real,mem_t,dim,ordinal>;

typedef t_serial_map            t_map;
typedef t_tensor0_field_tml<real>   t_vars;

DEFINE_CONSTANT_BUFFER(t_gpu_mesh, mesh)

__device__ int  get_elem_type(int i)    //gmsh element type
{
    if (mesh().is_homogeneous) return mesh().homogeneous_elem_type; else return mesh().elem_type(i,0);
}

//these service function related to mesh reference information could be placed in gpu_mesh
__device__ int  get_elem_faces_n(int elem_type)
{
    //ISSUE to make faster access (perhaps, not tested yet) we can make constant array
    if (elem_type == 4) return 4;
    if (elem_type == 5) return 6;
    if (elem_type == 6) return 5;
    if (elem_type == 7) return 5;
    //TODO others, error
    return 0;   //just to remove warning
}

__device__ int  get_elem_vert_n(int elem_type)
{
    //ISSUE to make faster access (perhaps, not tested yet) we can make constant array
        if (elem_type == 4) return 4;
        if (elem_type == 5) return 8;
        if (elem_type == 6) return 6;
        if (elem_type == 7) return 5;
        //TODO others, error
        return 0;   //just to remove warning
}

//returns p0 reflected with respect to plane with normal norm and point p1 on it
__device__ t_vec reflect_point(const t_vec &norm, const t_vec &p1, const t_vec &p0)
{
    real d = -scalar_prod(norm,p1);
    //could be of any sign
    real    dest = std::abs(scalar_prod(norm,p0) + d);
    return p0 + norm*(real(2.f)*dest);
}

__global__ void ker_poisson_iteration(t_vars vars_old, t_vars vars_new, int bnd1, int bnd2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= 0)&&(i < mesh().n_cv))) return;

        int elem_type = get_elem_type(i);
    real    numerator(0.f), denominator(0.f);
    for (int j = 0;j < get_elem_faces_n(elem_type);++j) {
        int nb = mesh().Neighbour(i,j);
        t_vec   nb_center;
        real    dist, var_nb;
        if (nb != CUDA_EMPTY_IDX) {
            mesh().center_neighbour.getv(i,j,nb_center);
            var_nb = vars_old(nb);
        } else {
            nb_center = reflect_point(mesh().Norm.getv(i,j), mesh().center_faces.getv(i,j), mesh().center.getv(i));
            if (mesh().Boundary(i,j) == bnd1) {
                //dirichle 0. value
                var_nb = -vars_old(i);
            } else if (mesh().Boundary(i,j) == bnd2) {
                //dirichle 1. value
                var_nb = real(2.f)*real(1.f)-vars_old(i);
            } else {
                //neumann
                var_nb = vars_old(i);
            }
        }
        dist = scalar_prod(mesh().Norm.getv(i,j), nb_center - mesh().center.getv(i));

        numerator += mesh().faces_S(i,j)*var_nb/dist;
        denominator += mesh().faces_S(i,j)/dist;
    }
    vars_new(i) = numerator/denominator;
}

//B := A
__global__ void ker_fill_zero(t_vars A)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= 0)&&(i < mesh().n_cv))) return;
        A(i) = real(0.f);
}

//B := A
__global__ void ker_assign(t_vars B, t_vars A)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= 0)&&(i < mesh().n_cv))) return;
        B(i) = A(i);
}

int main(int argc, char **args)
{
    LogStd          log;
    t_cpu_mesh      cpu_mesh;
        t_map           map;
    t_gpu_mesh      gpu_mesh;
    t_vars          vars0, vars1;
    int         device_number;
    int         bnd1, bnd2, iters_num;
    dim3            dimBlock, dimGrid;

    USE_MAIN_TRY_CATCH(log)

    log.set_verbosity(1);

    //process args
    if (argc < 5) {
        printf("usage: gpu_poisson_solver device_number bnd1_id bnd2_id iters_num\n");
        printf("    mesh is read from 'mesh.dat' file\n");
        printf("example: ./gpu_poisson_solver 0 6 166 1000\n");
        return 1;
    }
    device_number = atoi(args[1]);
    bnd1 = atoi(args[2]);
    bnd2 = atoi(args[3]);
    iters_num = atoi(args[4]);

    MAIN_TRY("reading mesh")
        if (!cpu_mesh.read("mesh.dat")) throw std::runtime_error("failed to read mesh from mesh.dat");
        //init map object
        map = t_serial_map(cpu_mesh.cv.size());
        //TODO add 0th order stencil through add_stencil_element()
        map.complete();
    MAIN_CATCH(2)

    MAIN_TRY("init CUDA")
        if (!InitCUDA(device_number)) throw std::runtime_error("InitCUDA failed");
    MAIN_CATCH(3)

    MAIN_TRY("allocate memory for mesh in device")
        gpu_mesh.init(map);
    MAIN_CATCH(4)

    MAIN_TRY("copy mesh data to device")
        init_gpu_mesh(map, gpu_mesh, cpu_mesh);
        //copy info about gpu mesh to gpu constant buffer
                COPY_TO_CONSTANT_BUFFER(mesh, gpu_mesh);
    MAIN_CATCH(5)
    
    dimBlock = dim3(128);
    dimGrid = dim3((gpu_mesh.n_cv+128)/dimBlock.x);

    MAIN_TRY("allocating variables array")
        //TODO we could make cool init using MAP concept, like init init_local methods
        vars0.init(map.max_loc_ind() - map.min_loc_ind() + 1, map.min_loc_ind());
        vars1.init(map.max_loc_ind() - map.min_loc_ind() + 1, map.min_loc_ind());
        MAIN_CATCH(6)

        MAIN_TRY("iterate poisson equation")
        ker_fill_zero<<<dimGrid, dimBlock>>>(vars0);
        ker_fill_zero<<<dimGrid, dimBlock>>>(vars1);
            for (int i = 0;i < iters_num;++i) {
            log.info_f("iteration %d", i);
            //put result in vars1
            ker_poisson_iteration<<<dimGrid, dimBlock>>>(vars0, vars1, bnd1, bnd2);
            //vars0 := vars1
                        ker_assign<<<dimGrid, dimBlock>>>(vars0, vars1);
        }
        MAIN_CATCH(7)
        
        MAIN_TRY("writing pos output to result.pos")
            write_out_pos_scalar_file("result.pos", "poisson_phi", cpu_mesh, map, vars0);
        MAIN_CATCH(8)

        //NOTE memory deallocation done in destrcutor automatically

    return 0;
}
