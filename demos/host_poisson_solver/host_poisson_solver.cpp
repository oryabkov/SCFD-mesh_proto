// Copyright © 2016-2021 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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

#include "host_poisson_solver_config.h"
#include <string>
#include <cmath>
#include <scfd/utils/log_std.h>
#include <scfd/utils/main_try_catch_macro.h>
#include <scfd/communication/linear_partitioner.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>
#include <scfd/mesh/host_mesh.h>
#include "gmsh_pos_output.h"

using real = HOST_POISSON_SOLVER_SCALAR_TYPE;
using ordinal = int;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,3,ordinal>;
using host_mesh_t = scfd::mesh::host_mesh<gmsh_wrap_t>;
using vec_t = scfd::static_vec::vec<real,3>;
using log_t = scfd::utils::log_std;

int bnd1, bnd2, iters_num;

//returns p0 reflected with respect to plane with normal norm and point p1 on it
vec_t reflect_point(const vec_t &norm, const vec_t &p1, const vec_t &p0)
{
    real d = -scalar_prod(norm,p1);
    //could be of any sign
    real    dest = std::abs(scalar_prod(norm,p0) + d);
    return p0 + norm*(real(2.f)*dest);
}


void    poisson_iteration(const host_mesh_t &host_mesh, const real *vars_old, real *vars_new)
{
    for (int i = 0;i < host_mesh.cv.size();++i) 
    {
        real    numerator(0.f), denominator(0.f);
        for (int j = 0;j < host_mesh.cv[i].faces_n;++j) 
        {
            int nb = host_mesh.cv[i].neighbours[j];
            vec_t   nb_center;
            real    dist, var_nb;
            if (nb != -1) 
            {
                nb_center = host_mesh.cv[nb].center;
                var_nb = vars_old[nb];
            } 
            else 
            {
                nb_center = reflect_point(host_mesh.cv[i].norms[j], host_mesh.cv[i].face_centers[j], host_mesh.cv[i].center);
                if (host_mesh.cv[i].boundaries[j] == bnd1) 
                {
                    //dirichle 0. value
                    var_nb = -vars_old[i];
                } 
                else if (host_mesh.cv[i].boundaries[j] == bnd2) 
                {
                    //dirichle 1. value
                    var_nb = real(2.f)*real(1.f)-vars_old[i];
                } 
                else 
                {
                    //neumann
                    var_nb = vars_old[i];
                }
            }
            dist = scalar_prod(host_mesh.cv[i].norms[j], nb_center - host_mesh.cv[i].center);
            numerator += host_mesh.cv[i].S[j]*var_nb/dist;
            denominator += host_mesh.cv[i].S[j]/dist;
        }
                vars_new[i] = numerator/denominator;
    }
}

void    fill_zero(const host_mesh_t &host_mesh, real *A)
{
    for (int i = 0;i < host_mesh.cv.size();++i) 
    {
        A[i] = 0.f;
    }
}

//B := A
void    assign(const host_mesh_t &host_mesh, real *B, const real *A)
{
    for (int i = 0;i < host_mesh.cv.size();++i) 
    {
        B[i] = A[i];
    }
}

int main(int argc, char **args)
{
    log_t           log;
    host_mesh_t     host_mesh;
    real            *vars0, *vars1;

    USE_MAIN_TRY_CATCH(log)
    
    log.set_verbosity(1);
    
    //process args
    if (argc < 4) 
    {
        printf("usage: host_poisson_solver bnd1_id bnd2_id iters_num\n");
        printf("example: ./cpu_poisson_solver 6 166 1000\n");
        return 1;
    }
    bnd1 = atoi(args[1]);
    bnd2 = atoi(args[2]);
    iters_num = atoi(args[3]);

    MAIN_TRY("reading mesh from mesh.dat")
    if (!host_mesh.read("mesh.dat")) throw std::runtime_error("failed to read mesh from mesh.dat");
    MAIN_CATCH(2)

    MAIN_TRY("allocating variables array")
    vars0 = new real[host_mesh.cv.size()];
    vars1 = new real[host_mesh.cv.size()];
    MAIN_CATCH(3)

    MAIN_TRY("iterate poisson equation")
    fill_zero(host_mesh, vars0);
    fill_zero(host_mesh, vars1);
    for (int i = 0;i < iters_num;++i) 
    {
        log.info_f("iteration %d", i);
        //put result in vars1
        poisson_iteration(host_mesh, vars0, vars1);
        //vars0 := vars1
        assign(host_mesh, vars0, vars1);
    }
    MAIN_CATCH(4)

    MAIN_TRY("writing pos output to result.pos")
    write_out_pos_scalar_file("result.pos", "poisson_phi", host_mesh, vars0);
    MAIN_CATCH(5)

    MAIN_TRY("deallocating variables array")
    delete []vars0;
    delete []vars1;
    MAIN_CATCH(6)

    return 0;
}
