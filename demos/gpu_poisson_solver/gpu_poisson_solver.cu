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

//#define SCFD_MAIN_TRY_CATCH_DISABLE_CATCH 1

#include "gpu_poisson_solver_config.h"
#include <string>
#include <cmath>
#include <scfd/utils/log_std.h>
#include <scfd/utils/main_try_catch_macro.h>
#include <scfd/utils/init_cuda.h>
#include <scfd/utils/constant_data.h>
#include <scfd/static_vec/vec.h>
#include <scfd/memory/cuda.h>
#include <scfd/for_each/cuda.h>
#include <scfd/for_each/cuda_impl.cuh>
#include <scfd/communication/linear_partitioner.h>
//TODO add serial_map
#include <scfd/communication/serial_map.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>
#include <scfd/mesh/host_mesh.h>
#include <scfd/mesh/device_mesh.h>
#include <scfd/mesh/device_mesh_impl.h>
#include "gmsh_pos_output.h"
#include "map_mock.h"

using real = GPU_POISSON_SOLVER_SCALAR_TYPE;
using ordinal = int;
static const ordinal dim = 3;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,dim,ordinal>;
using host_mesh_t = scfd::mesh::host_mesh<gmsh_wrap_t>;
using vec_t = scfd::static_vec::vec<real,dim>;
using log_t = scfd::utils::log_std;
using mem_t = scfd::memory::cuda_device;
using device_mesh_t = scfd::mesh::device_mesh<real,mem_t,dim,ordinal>;
using host_real_vector_t = std::vector<real>;
using for_each_t = scfd::for_each::cuda<>;

SCFD_DEVICE_MESH_INSTANTIATE(real,mem_t,dim,ordinal)

typedef scfd::communication::serial_map  map_t;
typedef scfd::arrays::array<real,mem_t>  vars_t;

DEFINE_CONSTANT_BUFFER(device_mesh_t, mesh)

__device__ int  get_elem_type(int i)    //gmsh element type
{
    if (mesh().is_homogeneous) return mesh().homogeneous_elem_type; else return mesh().elems_types(i);
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
__device__ vec_t reflect_point(const vec_t &norm, const vec_t &p1, const vec_t &p0)
{
    real d = -scalar_prod(norm,p1);
    //could be of any sign
    real    dest = std::abs(scalar_prod(norm,p0) + d);
    return p0 + norm*(real(2.f)*dest);
}

struct bnd_cond_data_t
{
    ordinal     dirichle_bnds_num;
    ordinal     dirichle_bnds[10];
    real        dirichle_vals[10];
};

DEFINE_CONSTANT_BUFFER(bnd_cond_data_t, bnd_cond_data)

struct force_data_t
{
    /// amplitude
    real    a;
    /// space multipliers
    real    omega[dim];
    /// space phase shifts
    real    phi[dim];
};

DEFINE_CONSTANT_BUFFER(force_data_t, force_data)

__global__ void ker_poisson_iteration(vars_t vars_old, vars_t vars_new, int bnd1, int bnd2)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= mesh().own_elems_range.i0)&&(i < mesh().own_elems_range.i0 + mesh().own_elems_range.n))) return;

    int     elem_type = get_elem_type(i);
    real    numerator(0.f), denominator(0.f);
    for (int j = 0;j < get_elem_faces_n(elem_type);++j) {
        int nb = mesh().elems_neighbours0(i,j);
        vec_t   nb_center;
        real    dist, var_nb;
        if (nb != device_mesh_t::special_id) {
            mesh().elems_neighbours0_centers.get_vec(nb_center,i,j);
            var_nb = vars_old(nb);
        } else {
            nb_center = reflect_point(mesh().elems_faces_norms.get_vec(i,j), mesh().elems_faces_centers.get_vec(i,j), mesh().elems_centers.get_vec(i));
            //real x = mesh().elems_faces_centers(i,j,0);
            ordinal  bnd_group_id = mesh().elems_faces_group_ids(i,j);
            bool     dirichle_found = false;
            real     dirichle_val;
            for (ordinal bnd_i = 0;bnd_i < bnd_cond_data().dirichle_bnds_num;++bnd_i)
            {
                if (bnd_group_id == bnd_cond_data().dirichle_bnds[bnd_i])
                {
                    dirichle_found = true;
                    dirichle_val = bnd_cond_data().dirichle_vals[bnd_i];
                    break;
                }
            }
            if (dirichle_found) {
                //dirichle value
                var_nb = real(2.f)*dirichle_val-vars_old(i);
            } else {
                //neumann
                var_nb = vars_old(i);
            }
            //var_nb = real(2.f)*x-vars_old(i);
        }
        dist = scalar_prod(mesh().elems_faces_norms.get_vec(i,j), nb_center - mesh().elems_centers.get_vec(i));

        numerator += mesh().elems_faces_areas(i,j)*var_nb/dist;
        denominator += mesh().elems_faces_areas(i,j)/dist;
    }
    vars_new(i) = numerator/denominator;
    //real x = mesh().elems_centers(i,0);
    //vars_new(i) = x;
}

//B := A
__global__ void ker_fill_zero(vars_t A)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= mesh().own_elems_range.i0)&&(i < mesh().own_elems_range.i0 + mesh().own_elems_range.n))) return;
    A(i) = real(0.f);
}

//B := A
__global__ void ker_assign(vars_t B, vars_t A)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (!((i >= mesh().own_elems_range.i0)&&(i < mesh().own_elems_range.i0 + mesh().own_elems_range.n))) return;
    B(i) = A(i);
}

void read_bnd_data(const std::string &fn, bnd_cond_data_t &res, std::set<ordinal> &periodic_bnds)
{
    std::ifstream   f(fn.c_str());
    ordinal         bnds_n;
    if (!(f >> bnds_n))
        throw std::runtime_error("read_bnd_data: failed to read from file " + fn);
    res.dirichle_bnds_num = 0;
    for (ordinal bnd_i = 0;bnd_i < bnds_n;++bnd_i)
    {
        ordinal         bnd_id;
        std::string     bnd_type;
        real            bnd_val;
        if (!(f >> bnd_id >> bnd_type >> bnd_val))
            throw std::runtime_error("read_bnd_data: failed to read from file " + fn);

        if (bnd_type == "D")
        {
            res.dirichle_bnds[res.dirichle_bnds_num] = bnd_id;
            res.dirichle_vals[res.dirichle_bnds_num] = bnd_val;
            ++res.dirichle_bnds_num;
            if (res.dirichle_bnds_num > 10)
                throw std::runtime_error("read_bnd_data: res.dirichle_bnds_num > 10");
        } 
        else if (bnd_type == "P")
        {
            periodic_bnds.insert(bnd_id);
        }
        else
        {
            throw std::runtime_error("read_bnd_data: unknown boundary type " + bnd_type);
        }
    }
}

force_data_t read_force_data(const std::string &fn)
{
    force_data_t res;

    std::ifstream   f(fn.c_str());
    if (!(f >> res.a))
        throw std::runtime_error("read_force_data: failed to read from file " + fn);
    for (int j = 0;j < dim;++j)
    {
        if (!(f >> res.omega[j] >> res.phi[j]))
            throw std::runtime_error("read_force_data: failed to read from file " + fn);
    }

    return res;
}

int main(int argc, char **args)
{
    std::string         mesh_fn, bnd_fn, force_fn;
    int                 device_number;
    ordinal             iters_num;
    std::set<ordinal>   periodic_bnds;

    log_t               log;
    auto                part = std::make_shared<partitioner_t>();
    auto                host_mesh = std::make_shared<host_mesh_t>();

    auto                map = std::make_shared<map_t>();
    device_mesh_t       gpu_mesh;
    host_real_vector_t  vars_host;
    vars_t              vars0, vars1;
    dim3                dimBlock, dimGrid;

    USE_MAIN_TRY_CATCH(log)

    log.set_verbosity(1);

    //process args
    if (argc < 6)
    {
        printf("Usage: ./gpu_poisson_solver.bin DEVICE_NUMBER MESH_FN BNDS_FN FORCE_FN ITERS_NUM\n");
        printf("Example: ./gpu_poisson_solver.bin 0 test.msh 5 27 1000\n");
        return 1;
    }
    device_number = atoi(args[1]);
    mesh_fn = args[2];
    bnd_fn = args[3];
    force_fn = args[4];
    iters_num = atoi(args[5]);

    MAIN_TRY("reading boundary data from from " + bnd_fn)
    bnd_cond_data_t bnd_cond_data_host;
    read_bnd_data(bnd_fn, bnd_cond_data_host, periodic_bnds);
    COPY_TO_CONSTANT_BUFFER(bnd_cond_data, bnd_cond_data_host);
    MAIN_CATCH(2)

    MAIN_TRY("reading force data from from " + force_fn)
    force_data_t force_data_host = read_force_data(force_fn);
    COPY_TO_CONSTANT_BUFFER(force_data, force_data_host);
    MAIN_CATCH(2)

    MAIN_TRY("reading mesh from " + mesh_fn)
    host_mesh->set_mesh_filename(mesh_fn);
    host_mesh->read(periodic_bnds);
    *part = partitioner_t(host_mesh->get_total_elems_num(), 1, 0);
    host_mesh->set_partitioner(part);
    host_mesh->enlarge_stencil(1);
    //init map object
    *map = map_t(host_mesh->get_total_elems_num());
    //TODO add 0th order stencil through add_stencil_element()
    map->complete();
    MAIN_CATCH(2)

    MAIN_TRY("init CUDA")
    //TODO
    //if (!InitCUDA(device_number)) throw std::runtime_error("InitCUDA failed");
    scfd::utils::init_cuda(-2, 0);
    MAIN_CATCH(3)

    MAIN_TRY("allocate memory for mesh in device and copy mesh data to device")
    for_each_t  for_each;
    gpu_mesh.params.has_elems_nodes_data = false;
    gpu_mesh.init_elems_data
    (
        *host_mesh, *map, map_mock(), map_mock(), for_each
    );
    COPY_TO_CONSTANT_BUFFER(mesh, gpu_mesh);
    MAIN_CATCH(4)
    
    //TODO
    dimBlock = dim3(128);
    dimGrid = dim3((gpu_mesh.own_elems_range.i1() + dimBlock.x)/dimBlock.x);

    MAIN_TRY("allocating variables array")
    //TODO we could make cool init using MAP concept, like init init_local methods
    //TODO
    vars0.init(map->max_loc_ind() - map->min_loc_ind() + 1, map->min_loc_ind());
    vars1.init(map->max_loc_ind() - map->min_loc_ind() + 1, map->min_loc_ind());
    vars_host = host_real_vector_t(map->get_total_size());
    MAIN_CATCH(5)

    MAIN_TRY("iterate poisson equation")
    ker_fill_zero<<<dimGrid, dimBlock>>>(vars0);
    ker_fill_zero<<<dimGrid, dimBlock>>>(vars1);
    for (int i = 0;i < iters_num;++i) 
    {
        log.info_f("iteration %d", i);
        //put result in vars1
        ker_poisson_iteration<<<dimGrid, dimBlock>>>(vars0, vars1, bnd1, bnd2);
        //vars0 := vars1
        ker_assign<<<dimGrid, dimBlock>>>(vars0, vars1);
    }
    MAIN_CATCH(6)

    MAIN_TRY("copy results to host")
    auto                    vars0_view = vars0.create_view(true);
    //auto                    vars0_view = gpu_mesh.elems_vols.create_view(true);
    for (ordinal i_ = 0;i_ < map->get_size();++i_) 
    {
        ordinal         i_glob = map->own_glob_ind(i_);
        ordinal         i = map->own_loc_ind(i_);
        vars_host[i_glob] = vars0_view(i);
        //vars_host[i_glob] = vars0_view(i,0);
    }
    vars0_view.release(false);
    MAIN_CATCH(7)
    
    MAIN_TRY("writing pos output to result.pos")
    write_out_pos_scalar_file("result.pos", "poisson_phi", *host_mesh, vars_host);
    MAIN_CATCH(8)

    //NOTE memory deallocation done in destrcutor automatically

    return 0;
}
