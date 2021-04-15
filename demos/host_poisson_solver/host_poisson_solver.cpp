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

#include "cpu_poisson_solver_config.h"
#include <string>
#include <cmath>
#include <utils/Log.h>
#include <utils/main_try_catch_macro.h>
#include <mesh/t_cpu_mesh_tml.h>
#include "gmsh_pos_output.h"

typedef CPU_POISSON_SOLVER_SCALAR_TYPE  real;
typedef t_cpu_mesh_tml<real>            t_cpu_mesh;
typedef t_vec_tml<real,3>               t_vec;

int bnd1, bnd2, iters_num;

//returns p0 reflected with respect to plane with normal norm and point p1 on it
t_vec reflect_point(const t_vec &norm, const t_vec &p1, const t_vec &p0)
{
    real d = -scalar_prod(norm,p1);
    //could be of any sign
    real    dest = std::abs(scalar_prod(norm,p0) + d);
    return p0 + norm*(real(2.f)*dest);
}


void    poisson_iteration(const t_cpu_mesh &cpu_mesh, const real *vars_old, real *vars_new)
{
    for (int i = 0;i < cpu_mesh.cv.size();++i) {
        real    numerator(0.f), denominator(0.f);
        for (int j = 0;j < cpu_mesh.cv[i].faces_n;++j) {
            int nb = cpu_mesh.cv[i].neighbours[j];
            t_vec   nb_center;
            real    dist, var_nb;
            if (nb != -1) {
                nb_center = cpu_mesh.cv[nb].center;
                var_nb = vars_old[nb];
            } else {
                nb_center = reflect_point(cpu_mesh.cv[i].norms[j], cpu_mesh.cv[i].face_centers[j], cpu_mesh.cv[i].center);
                if (cpu_mesh.cv[i].boundaries[j] == bnd1) {
                    //dirichle 0. value
                    var_nb = -vars_old[i];
                } else if (cpu_mesh.cv[i].boundaries[j] == bnd2) {
                    //dirichle 1. value
                    var_nb = real(2.f)*real(1.f)-vars_old[i];
                } else {
                    //neumann
                    var_nb = vars_old[i];
                }
            }
            dist = scalar_prod(cpu_mesh.cv[i].norms[j], nb_center - cpu_mesh.cv[i].center);
            numerator += cpu_mesh.cv[i].S[j]*var_nb/dist;
            denominator += cpu_mesh.cv[i].S[j]/dist;
        }
                vars_new[i] = numerator/denominator;
    }
}

void    fill_zero(const t_cpu_mesh &cpu_mesh, real *A)
{
    for (int i = 0;i < cpu_mesh.cv.size();++i) {
        A[i] = 0.f;
    }
}

//B := A
void    assign(const t_cpu_mesh &cpu_mesh, real *B, const real *A)
{
    for (int i = 0;i < cpu_mesh.cv.size();++i) {
        B[i] = A[i];
    }
}

int main(int argc, char **args)
{
    LogStd          log;
    t_cpu_mesh      cpu_mesh;
    real            *vars0, *vars1;

    USE_MAIN_TRY_CATCH(log)
    
    log.set_verbosity(1);
    
    //process args
    if (argc < 4) {
        printf("usage: cpu_poisson_solver bnd1_id bnd2_id iters_num\n");
        printf("example: ./cpu_poisson_solver 6 166 1000\n");
        return 1;
    }
    bnd1 = atoi(args[1]);
    bnd2 = atoi(args[2]);
    iters_num = atoi(args[3]);

    MAIN_TRY("reading mesh from mesh.dat")
        if (!cpu_mesh.read("mesh.dat")) throw std::runtime_error("failed to read mesh from mesh.dat");
    MAIN_CATCH(2)

    MAIN_TRY("allocating variables array")
        vars0 = new real[cpu_mesh.cv.size()];
        vars1 = new real[cpu_mesh.cv.size()];
        MAIN_CATCH(3)

        MAIN_TRY("iterate poisson equation")
            fill_zero(cpu_mesh, vars0);
            fill_zero(cpu_mesh, vars1);
            for (int i = 0;i < iters_num;++i) {
            log.info_f("iteration %d", i);
            //put result in vars1
            poisson_iteration(cpu_mesh, vars0, vars1);
            //vars0 := vars1
                        assign(cpu_mesh, vars0, vars1);
        }
        MAIN_CATCH(4)

        MAIN_TRY("writing pos output to result.pos")
            write_out_pos_scalar_file("result.pos", "poisson_phi", cpu_mesh, vars0);
        MAIN_CATCH(5)

        MAIN_TRY("deallocating variables array")
        delete []vars0;
        delete []vars1;
        MAIN_CATCH(6)

    return 0;
}
