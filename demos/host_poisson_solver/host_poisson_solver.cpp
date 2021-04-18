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
#include <scfd/static_vec/vec.h>
#include <scfd/communication/linear_partitioner.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>
#include <scfd/mesh/host_mesh.h>
#include "gmsh_pos_output.h"

using real = HOST_POISSON_SOLVER_SCALAR_TYPE;
using ordinal = int;
static const ordinal dim = 3;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,dim,ordinal>;
using host_mesh_t = scfd::mesh::host_mesh<gmsh_wrap_t>;
using vec_t = scfd::static_vec::vec<real,dim>;
using log_t = scfd::utils::log_std;
using real_vector_t = std::vector<real>;
using vec_vector_t = std::vector<vec_t>;
template<class T>
class elems_faces_vector_t 
{
public:
    elems_faces_vector_t()
    {
        /// elems_num() must work for empty array
        max_faces_per_elem_ = 1;
    }
    elems_faces_vector_t(ordinal elems_n,ordinal max_faces_per_elem)
    {
        max_faces_per_elem_ = max_faces_per_elem;
        vals_.resize(elems_n*max_faces_per_elem);
    }

    ordinal     elems_num()const
    {
        return vals_.size()/max_faces_per_elem;
    }
    ordinal     max_faces_per_elem()const
    {
        return max_faces_per_elem_;
    }

    T          &operator()(ordinal elem_id, ordinal face_i)
    {
        return vals_[max_faces_per_elem_*elem_id + face_i];
    }
    const T    &operator()(ordinal elem_id, ordinal face_i)const
    {
        return vals_[max_faces_per_elem_*elem_id + face_i];
    }

private:
    ordinal             max_faces_per_elem_;
    std::vector<T>      vals_;
};
using real_elems_faces_vector_t = elems_faces_vector_t<real>;
using vec_elems_faces_vector_t = elems_faces_vector_t<vec_t>;

std::string mesh_fn;
ordinal     bnd1, bnd2, iters_num;

//returns p0 reflected with respect to plane with normal norm and point p1 on it
vec_t reflect_point(const vec_t &norm, const vec_t &p1, const vec_t &p0)
{
    real d = -scalar_prod(norm,p1);
    //could be of any sign
    real    dest = std::abs(scalar_prod(norm,p0) + d);
    return p0 + norm*(real(2.f)*dest);
}


void    poisson_iteration
(
    const host_mesh_t &host_mesh, const vec_vector_t &centers,
    const real_vector_t &vols, const real_elems_faces_vector_t &face_areas,
    const vec_elems_faces_vector_t &norms, const vec_elems_faces_vector_t &face_centers,
    const real_vector_t &vars_old, real_vector_t &vars_new
)
{
    ordinal neibs[host_mesh.get_elems_max_faces_num()];
    for (int i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        real    numerator(0.f), denominator(0.f);
        host_mesh.get_elem_neighbours0(i, neibs);
        for (int j = 0;j < host_mesh.get_elem_faces_num(i);++j) 
        {
            int nb = neibs[j];
            vec_t   nb_center;
            real    dist, var_nb;
            if (nb != host_mesh_t::special_id) 
            {
                nb_center = centers[nb];
                var_nb = vars_old[nb];
            } 
            else 
            {
                nb_center = reflect_point(norms(i,j), face_centers(i,j), centers[i]);
                if (host_mesh.get_elem_face_group_id(i,j) == bnd1) 
                {
                    //dirichle 0. value
                    var_nb = -vars_old[i];
                } 
                else if (host_mesh.get_elem_face_group_id(i,j) == bnd2) 
                {
                    //dirichle 1. value
                    var_nb = real(2)*real(1)-vars_old[i];
                } 
                else 
                {
                    //neumann
                    var_nb = vars_old[i];
                }
            }
            dist = scalar_prod(norms(i,j), nb_center - centers[i]);
            numerator += face_areas(i,j)*var_nb/dist;
            denominator += face_areas(i,j)/dist;
        }
        vars_new[i] = numerator/denominator;
    }
}

real_vector_t   allocate_real_vector(const host_mesh_t &host_mesh)
{
    return real_vector_t(host_mesh.get_total_elems_num());
}

vec_vector_t    allocate_vec_vector(const host_mesh_t &host_mesh)
{
    return vec_vector_t(host_mesh.get_total_elems_num());
}

real_elems_faces_vector_t   allocate_real_elems_faces_vector(const host_mesh_t &host_mesh)
{
    return real_elems_faces_vector_t(host_mesh.get_total_elems_num(),host_mesh.get_elems_max_faces_num());
}

vec_elems_faces_vector_t   allocate_vec_elems_faces_vector(const host_mesh_t &host_mesh)
{
    return vec_elems_faces_vector_t(host_mesh.get_total_elems_num(),host_mesh.get_elems_max_faces_num());
}

void    fill_zero(const host_mesh_t &host_mesh, real_vector_t &A)
{
    for (int i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        A[i] = 0.f;
    }
}

//B := A
void    assign(const host_mesh_t &host_mesh, real_vector_t &B, const real_vector_t &A)
{
    for (ordinal i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        B[i] = A[i];
    }
}

void    calc_centers(const host_mesh_t &host_mesh, vec_vector_t &centers)
{
    for (ordinal i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        vec_t       center = vec_t::make_zero();

        ordinal     elem_nodes[host_mesh.get_elems_max_nodes_num()];
        ordinal     nodes_n;
        host_mesh.get_elem_nodes(i, elem_nodes, &nodes_n);
        for (ordinal j = 0;j < nodes_n;++j) 
        {
            vec_t   vertex = host_mesh.get_node_coords(elem_nodes[j]);
            center += vertex;
        }
        center /= nodes_n;

        centers[i] = center;
    }
}

real get_tetrahedron_volume(real tet_verteces[4][dim])
{
    for (int vertex_j = 0;vertex_j < 3;vertex_j++)
    {
        for (int k = 0;k < dim;k++)
        {
              tet_verteces[vertex_j][k] -= tet_verteces[3][k];  
        }
    }
    //we now take tet_verteces[3][*] equals (0,0,0)
    real tet_volume;
    tet_volume = std::abs(+tet_verteces[0][0]*(tet_verteces[1][1]*tet_verteces[2][2]-tet_verteces[2][1]*tet_verteces[1][2])
                          -tet_verteces[0][1]*(tet_verteces[1][0]*tet_verteces[2][2]-tet_verteces[2][0]*tet_verteces[1][2])
                          +tet_verteces[0][2]*(tet_verteces[1][0]*tet_verteces[2][1]-tet_verteces[2][0]*tet_verteces[1][1]))/real(6.f);

    return tet_volume;
}

void    calc_vols(const host_mesh_t &host_mesh, real_vector_t &vols)
{
    for (ordinal i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        ordinal     elem_nodes[host_mesh.get_elems_max_nodes_num()];
        ordinal     nodes_n;
        host_mesh.get_elem_nodes(i, elem_nodes, &nodes_n);
        
        auto elem_type = host_mesh.get_elem_type(i);
        if (elem_type == 4) 
        {
            real m[4][dim];
            for (ordinal i1 = 0;i1 < 4;i1++)
            {
                vec_t   vertex = host_mesh.get_node_coords(elem_nodes[i1]);
                for (ordinal i2 = 0;i2 < dim;i2++)
                {
                    m[i1][i2] =  vertex[i2];
                }
            }              
            vols[i] = get_tetrahedron_volume(m);
        }
        else 
            throw std::logic_error("calc_vols::not supported element type");
    }
}

//check orientation of normals
vec_t check_orientation_of_a_normal_vector(const vec_t &normal_vector, const vec_t &vertex_center_vector)
{
    if (scalar_prod(normal_vector, vertex_center_vector) > real(0))
    {
        return normal_vector.inverted();
    }      
    else
    {
        return normal_vector;
    }
}

//gets normal vector to the face defined by three verteces and uses center point to find outer orientation.
vec_t construct_normal_vector_for_triangle_face(const vec_t &vertex_1, const vec_t &vertex_2, const vec_t &vertex_3, const vec_t &center)
{
    vec_t vector_2_1 = vertex_2-vertex_1;
    vec_t vector_3_1 = vertex_3-vertex_1;
    vec_t normal_vector = vector_prod(vector_2_1, vector_3_1);
    vec_t vector_1_center = center-vertex_1;
    normal_vector = check_orientation_of_a_normal_vector(normal_vector, vector_1_center);
    return normal_vector;
}     

void    calc_face_areas_and_norms
(
    const host_mesh_t &host_mesh, const vec_vector_t &centers,
    real_elems_faces_vector_t &face_areas, vec_elems_faces_vector_t &norms
)
{
    const auto &ref = host_mesh.mesh_elem_reference();

    for (ordinal i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        vec_t  center = centers[i];

        ordinal     elem_nodes[host_mesh.get_elems_max_nodes_num()];
        host_mesh.get_elem_nodes(i, elem_nodes);

        auto elem_type = host_mesh.get_elem_type(i);

        for (ordinal face_i = 0;face_i < ref.get_faces_n(elem_type);++face_i)
        {
            if (ref.get_face_verts_n(elem_type,face_i)==3)
            {
                // This is shit, i must find a way to return mesh().vertexes(i,vert_j) as vec_t!
                vec_t  verteces[ref.get_face_verts_n(elem_type,face_i)];
                
                for (ordinal j = 0;j < ref.get_face_verts_n(elem_type,face_i);++j) 
                {
                    verteces[j] = host_mesh.get_node_coords(elem_nodes[ref.get_face_vert_i(elem_type,face_i,j)]);
                }

                vec_t normal_vector = construct_normal_vector_for_triangle_face(verteces[0], verteces[1], verteces[2], center);
                real face_area = normal_vector.norm2();
                
                face_areas(i,face_i) = face_area*real(0.5);
                norms(i,face_i) = normal_vector/face_area;
            }
            else
                throw std::runtime_error("calc_face_areas_and_norms::not supported face type");
        }
    }
}

void    calc_face_centers(const host_mesh_t &host_mesh, vec_elems_faces_vector_t &face_centers)
{
    const auto &ref = host_mesh.mesh_elem_reference();

    for (ordinal i = 0;i < host_mesh.get_total_elems_num();++i) 
    {
        ordinal     elem_nodes[host_mesh.get_elems_max_nodes_num()];
        host_mesh.get_elem_nodes(i, elem_nodes);

        auto elem_type = host_mesh.get_elem_type(i);

        for (ordinal face_i = 0;face_i < ref.get_faces_n(elem_type);++face_i)
        {
            vec_t       center = vec_t::make_zero();

            for (ordinal j = 0;j < ref.get_face_verts_n(elem_type,face_i);++j) 
            {
                vec_t   vertex = host_mesh.get_node_coords(elem_nodes[ref.get_face_vert_i(elem_type,face_i,j)]);
                center += vertex;
            }
            center /= ref.get_face_verts_n(elem_type,face_i);

            face_centers(i,face_i) = center;
        }
    }
}

int main(int argc, char **args)
{
    log_t           log;
    auto            part = std::make_shared<partitioner_t>();
    auto            host_mesh = std::make_shared<host_mesh_t>();

    vec_vector_t                centers;
    real_vector_t               vols; 
    real_elems_faces_vector_t   face_areas;
    vec_elems_faces_vector_t    norms; 
    vec_elems_faces_vector_t    face_centers;

    real_vector_t   vars0, vars1;

    USE_MAIN_TRY_CATCH(log)
    
    log.set_verbosity(1);
    
    //process args
    if (argc < 5) 
    {
        printf("Usage: host_poisson_solver MESH_FN BND1_ID BND2_ID ITERS_NUM\n");
        printf("Example: ./host_poisson_solver mesh.msh 6 166 1000\n");
        return 1;
    }

    mesh_fn = args[1];
    bnd1 = atoi(args[2]);
    bnd2 = atoi(args[3]);
    iters_num = atoi(args[4]);

    MAIN_TRY("reading mesh from " + mesh_fn)
    host_mesh->set_mesh_filename(mesh_fn);
    host_mesh->read();
    *part = partitioner_t(host_mesh->get_total_elems_num(), 1, 0);
    host_mesh->set_partitioner(part);
    host_mesh->enlarge_stencil(1);
    MAIN_CATCH(2)

    MAIN_TRY("allocating supplementary arrays")
    centers = allocate_vec_vector(*host_mesh);
    vols = allocate_real_vector(*host_mesh);
    face_areas = allocate_real_elems_faces_vector(*host_mesh);
    norms = allocate_vec_elems_faces_vector(*host_mesh);
    face_centers = allocate_vec_elems_faces_vector(*host_mesh);
    MAIN_CATCH(3)

    MAIN_TRY("calculating supplementary value")
    calc_centers(*host_mesh, centers);
    calc_vols(*host_mesh, vols);
    calc_face_areas_and_norms(*host_mesh, centers, face_areas, norms);
    calc_face_centers(*host_mesh, face_centers);
    MAIN_CATCH(4)

    MAIN_TRY("allocating variables array")
    vars0 = allocate_real_vector(*host_mesh);
    vars1 = allocate_real_vector(*host_mesh);
    MAIN_CATCH(5)

    MAIN_TRY("iterate poisson equation")
    fill_zero(*host_mesh, vars0);
    fill_zero(*host_mesh, vars1);
    for (int i = 0;i < iters_num;++i) 
    {
        log.info_f("iteration %d", i);
        //put result in vars1
        poisson_iteration
        (
            *host_mesh, centers, vols, face_areas, norms, face_centers, vars0, vars1
        );
        //vars0 := vars1
        assign(*host_mesh, vars0, vars1);
    }
    MAIN_CATCH(6)

    MAIN_TRY("writing pos output to result.pos")
    write_out_pos_scalar_file("result.pos", "poisson_phi", *host_mesh, vars0);
    MAIN_CATCH(7)

    MAIN_TRY("deallocating variables array")
    MAIN_CATCH(8)

    return 0;
}
