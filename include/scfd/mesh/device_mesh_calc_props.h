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

#ifndef __SCFD_MESH_DEVICE_MESH_CALC_PROPS_H__
#define __SCFD_MESH_DEVICE_MESH_CALC_PROPS_H__

#include "mesh_deform_shepard_config.h"
#include <set>
#include <utils/device_tag.h>
#include <utils/scalar_traits.h>
#include <utils/constant_data.h>
#include <tensor_field/t_tensor_field_tml.h>
#include <for_each/for_each_1d.h>
#ifdef MESH_DEFORM_SHEPARD_CUDA_NODE
#include <for_each/for_each_1d_cuda_impl.cuh>
#endif
#ifdef MESH_DEFORM_SHEPARD_OPENMP_NODE
#include <for_each/for_each_1d_openmp_impl.h>
#endif
#include <for_each/for_each_storage_types.h>
#include <for_each/for_each_func_macro.h>
#include <mesh/t_mesh_elem_reference_tml.h>
#include <mesh/t_gpu_mesh_tml.h>
#include "mesh_deform_shepard.h"

namespace scfd
{
namespace mesh
{
namespace detail
{

typedef scalar_traits<real>                                             sc_tr;
typedef t_vec_tml<real,dim>                                             t_vec;
typedef t_mesh_elem_reference_tml<real>                                 t_elem_reference;
typedef for_each_1d<for_each_type>                                      t_for_each_1d;

DEFINE_CONSTANT_BUFFER(t_gpu_mesh, mesh)
DEFINE_CONSTANT_BUFFER(t_elem_reference, elem_ref)

//TODO just temporal -> use mesh_reference instead
__DEVICE_TAG__ int      get_elem_faces_n(int elem_type)
{
    //ISSUE to make faster access (perhaps, not tested yet) we can make constant array
    if (elem_type == 4) return 4;
    if (elem_type == 5) return 6;
    if (elem_type == 6) return 5;
    if (elem_type == 7) return 5;
    //TODO others, error
    return 0;       //just to remove warning
}

__DEVICE_TAG__ int      get_elem_vert_n(int elem_type)
{
    //ISSUE to make faster access (perhaps, not tested yet) we can make constant array
    if (elem_type == 4) return 4;
    if (elem_type == 5) return 8;
    if (elem_type == 6) return 6;
    if (elem_type == 7) return 5;
    //TODO others, error
    return 0;       //just to remove warning
}

template<class T,class Memory,int Dim,class Ord>
struct device_mesh_funcs
{

    struct calc_center
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            int elem_type = mesh().get_elem_type(i);
            int max_vertex_in_the_element=elem_ref().get_verts_n(elem_type);
            t_vec  sum_of_verteces_coords(real(0.f), real(0.f), real(0.f));
            for(int vert_j=0;vert_j<max_vertex_in_the_element;vert_j++)
            {
                for (int k_spacial_dimension = 0;k_spacial_dimension < dim;k_spacial_dimension++)
                {

                    //t_tensor2_field_tml<T,max_vert_n,dim,strg>      vertexes;
                    //where is "operator()(int i, int j, int k)" ???
                    sum_of_verteces_coords[k_spacial_dimension]+=mesh().vertexes(i,vert_j,k_spacial_dimension);
                    //DEBUG:printf("%e ", sum_of_verteces_coords[k_spacial_dimension]);
                }
            }
            for (int k_spacial_dimension = 0;k_spacial_dimension < dim;k_spacial_dimension++)
            {
                   mesh().center(i,k_spacial_dimension)=sum_of_verteces_coords[k_spacial_dimension]/((real) max_vertex_in_the_element);
            }


        }
    };

    struct calc_center_faces
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            int elem_type = mesh().get_elem_type(i);
            for (int face_i = 0;face_i < elem_ref().get_faces_n(elem_type);++face_i)
            {
                int     face_verts_n = elem_ref().get_face_verts_n(elem_type,face_i);
                t_vec  c(real(0.f), real(0.f), real(0.f));
                for (int face_vert_i = 0;face_vert_i < face_verts_n;++face_vert_i) 
                {
                    int vert_i = elem_ref().get_face_vert_i(elem_type,face_i,face_vert_i);
                    for (int k = 0;k < dim;++k) 
                        c[k] += mesh().vertexes(i, vert_i, k);
                }
                for (int k = 0;k < dim;++k) 
                    c[k] /= real(face_verts_n);
                for (int k = 0;k < dim;++k) 
                    mesh().center_faces(i,face_i,k) = c[k];
            }
        }
    };


    //check orientation of normals
    __DEVICE_TAG__ t_vec check_orientation_of_a_normal_vector(const t_vec &normal_vector, const t_vec &vertex_center_vector)
    {
        if(scalar_prod(normal_vector, vertex_center_vector)>0.0)
        {
            return(vector_invert(normal_vector));
        }      
        else
        {
            return normal_vector;
        }

    }

    //return vector as vertex_1-vertex_2
    __DEVICE_TAG__ t_vec bulid_vector_from_two_verteces(const t_vec &vertex_1, const t_vec &vertex_2)
    {
        return(vertex_1-vertex_2);
    }

    //gets normal vector to the face defined by three verteces and uses center point to find outer orientation.
    __DEVICE_TAG__ t_vec construct_normal_vector_for_triangle_face(const t_vec &vertex_1, const t_vec &vertex_2, const t_vec &vertex_3, const t_vec &center)
    {
        t_vec vector_2_1 = bulid_vector_from_two_verteces(vertex_2,vertex_1);
        t_vec vector_3_1 = bulid_vector_from_two_verteces(vertex_3,vertex_1);
        t_vec normal_vector = vector_prod(vector_2_1, vector_3_1);
        t_vec vector_1_center = bulid_vector_from_two_verteces(center, vertex_1);
        normal_vector = check_orientation_of_a_normal_vector(normal_vector, vector_1_center);
        return normal_vector;
    }       

    //calculates normals to faces and araes of faces
    struct calc_norm
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            int elem_type=mesh().get_elem_type(i);
            for (int face_j=0;face_j<elem_ref().get_faces_n(elem_type);face_j++)
            {
                int face_verteces_number=elem_ref().get_face_verts_n(elem_type,face_j);
                if(face_verteces_number==3)
                {
                    // This is shit, i must find a way to return mesh().vertexes(i,vert_j) as t_vec!
                    t_vec  verteces[3]={(real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f))};
                    
                    t_vec  center = (real(0.f), real(0.f), real(0.f));
                    
                    for (int face_vert_j = 0;face_vert_j < face_verteces_number;face_vert_j++) { 
                        int vert_j = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                        for (int k = 0;k < dim;++k){ 
                            verteces[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                        }

                    }
                    for (int k = 0;k < dim;++k){ 
                        center[k] = mesh().center(i,k);
                    }

                    t_vec normal_vector = construct_normal_vector_for_triangle_face(verteces[0], verteces[1], verteces[2], center);
                    real face_area = normal_vector.norm2();
                    

                    for (int k = 0;k < dim;++k){ 
                        mesh().Norm(i,face_j,k)=normal_vector[k]/face_area;
                    }
                    mesh().faces_S(i,face_j)=face_area*real(0.5);

                }
                else if(face_verteces_number==4)
                {
                    t_vec face_center(real(0.f), real(0.f), real(0.f));

                    for (int k = 0;k < dim;++k)
                    { 
                        face_center[k]=mesh().center_faces(i,face_j,k);
                    }

                    t_vec  verteces[4] =
                        {
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f))
                        };
                    
                    t_vec  center = (real(0.f), real(0.f), real(0.f)); 

                    for (int face_vert_j = 0;face_vert_j < face_verteces_number;face_vert_j++) 
                    { 
                        int vert_j = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                        for (int k = 0;k < dim;++k)
                        { 
                            verteces[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                        }
                    }
                    for (int k = 0;k < dim;++k)
                    { 
                        center[k] = mesh().center(i,k);
                    }

                    t_vec normal_vector(real(0.f), real(0.f), real(0.f));
                    real face_area = real(0.0);
                    //loop over triangles
                    for (int k = 0;k < dim;++k)
                    { 
                        mesh().Norm(i,face_j,k)=real(0.0);
                    }

                    for(int triangle_number=0;triangle_number<face_verteces_number;triangle_number++)
                    {
                        int vertex_index_1=triangle_number;
                        int vertex_index_2=vertex_index_1+1;
                        vertex_index_2=(vertex_index_2)%(face_verteces_number);

                        normal_vector = construct_normal_vector_for_triangle_face(verteces[vertex_index_1], verteces[vertex_index_2], face_center, center);
                        real normal_vector_norm=normal_vector.norm2();
                        //for (int k = 0;k < dim;++k){ 
                        //        normal_vector[k]=normal_vector[k]/normal_vector_norm; //using equal wight for normals
                        //}
                        face_area+=normal_vector_norm;
                        for (int k = 0;k < dim;++k)
                        { 
                               // mesh().Norm(i,face_j,k)+=real(0.25)*normal_vector[k]; // using equal wight for normals. 1/face_verteces_number wight
                            mesh().Norm(i,face_j,k)+=normal_vector[k]; // wight 0.5 is relative to triangle area
                        }
                    
                    }
                    for (int k = 0;k < dim;++k)
                    { 
                        mesh().Norm(i,face_j,k)/=face_area; // used for triangle area wighted
                    } 
                    mesh().faces_S(i,face_j)=face_area*real(0.5);        
                }
                else
                {
                    //are there any elements with 5 verteces per face?
                }

            }
        }
    };

    struct calc_faces_S
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            //TODO      
            //Do we need it?!?          
        }
    };


    __DEVICE_TAG__ real get_tetrahedron_volume(real tet_verteces[4][dim])
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
        tet_volume = 
            sc_tr::abs
            (
                +tet_verteces[0][0]*(tet_verteces[1][1]*tet_verteces[2][2]-tet_verteces[2][1]*tet_verteces[1][2])
                -tet_verteces[0][1]*(tet_verteces[1][0]*tet_verteces[2][2]-tet_verteces[2][0]*tet_verteces[1][2])
                +tet_verteces[0][2]*(tet_verteces[1][0]*tet_verteces[2][1]-tet_verteces[2][0]*tet_verteces[1][1])
            )/real(6.f);

        return tet_volume;
    }


    __DEVICE_TAG__ real volume_of_4_tetraheda_from_4_verteces_and_center
    (
        t_vec face_verteces[4], const t_vec &face_center, const t_vec &element_center
    )
    {
        real tet_volumes=0.0;
        real m[4][dim];

        for(int tet_number=0;tet_number<4;tet_number++)
        {
            int vertex_index_1=tet_number;
            int vertex_index_2=vertex_index_1+1;
            vertex_index_2=(vertex_index_2)%(4);                

            for(int k=0;k<dim;k++)
            {
                m[0][k] = face_verteces[vertex_index_1][k];
                m[1][k] = face_verteces[vertex_index_2][k];
                m[2][k]=face_center[k];
                m[3][k]=element_center[k];
            }
            tet_volumes+=get_tetrahedron_volume(m);
        }

        return tet_volumes;
    }


    struct calc_vol
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            int elem_type = mesh().get_elem_type(i);
            if (elem_type == 4) 
            {
                /*
                real m[dim][dim]; //[number of verteces-1][dimension]
                for (int i1 = 0;i1 < dim;i1++)
                for (int i2 = 0;i2 < dim;i2++)
                    m[i1][i2] =  mesh().vertexes(i,i1+1,i2);
                
                t_vec   v0;
                mesh().vertexes.getv(i,0,v0);
                for (int i1 = 0;i1 < dim;i1++)
                for (int i2 = 0;i2 < dim;i2++)
                    m[i1][i2] -= v0[i2];
                mesh().Vol(i,0) = sc_tr::abs(  +m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
                                -m[0][1]*(m[1][0]*m[2][2]-m[2][0]*m[1][2])
                                +m[0][2]*(m[1][0]*m[2][1]-m[2][0]*m[1][1]))/real(6.f);
                */
                real m[4][dim];
                for (int i1 = 0;i1 < 4;i1++)
                {
                    for (int i2 = 0;i2 < dim;i2++)
                    {
                        m[i1][i2] =  mesh().vertexes(i,i1,i2);
                    }
                }              
                mesh().Vol(i,0) = get_tetrahedron_volume(m);

            }
            if (elem_type == 5) 
            {
                //hexahedron
                real hexahedron_volume = real(0.f);
                t_vec element_center(real(0.f), real(0.f), real(0.f));
                for (int k = 0;k < dim;++k)
                { 
                    element_center[k]=mesh().center(i,k);
                }
                int number_of_faces = elem_ref().get_faces_n(elem_type);
                for (int face_j=0;face_j<number_of_faces;face_j++)
                {
                    int face_verteces_number=elem_ref().get_face_verts_n(elem_type,face_j);
                    t_vec face_center(real(0.f), real(0.f), real(0.f));
                    t_vec  verteces[4]=
                        {
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f))
                        };

                    for (int k = 0;k < dim;++k)
                    { 
                        face_center[k]=mesh().center_faces(i,face_j,k);
                    }                                
                    for (int face_vert_j = 0;face_vert_j < face_verteces_number;face_vert_j++) 
                    { 
                        int vert_j = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                        for (int k = 0;k < dim;++k)
                        { 
                            verteces[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                        }
                    }                                
                    hexahedron_volume += volume_of_4_tetraheda_from_4_verteces_and_center(verteces, face_center, element_center);

                }
                mesh().Vol(i,0) = hexahedron_volume;

            }
            if (elem_type == 6) 
            {
                //prism
                real prism_volume = real(0.f);
                t_vec element_center(real(0.f), real(0.f), real(0.f));
                for (int k = 0;k < dim;++k){ 
                    element_center[k]=mesh().center(i,k);
                }
                int number_of_faces = elem_ref().get_faces_n(elem_type);
                for (int face_j=0;face_j<number_of_faces;face_j++)
                {
                    int face_verteces_number=elem_ref().get_face_verts_n(elem_type,face_j);
                    t_vec face_center(real(0.f), real(0.f), real(0.f));
                    t_vec  verteces[4]={(real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f)),
                            (real(0.f), real(0.f), real(0.f))};

                    for (int k = 0;k < dim;++k)
                    { 
                        face_center[k]=mesh().center_faces(i,face_j,k);
                    }                                
                    for (int face_vert_j = 0;face_vert_j < face_verteces_number;face_vert_j++) 
                    { 
                        int vert_j = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                        for (int k = 0;k < dim;++k)
                        { 
                            verteces[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                        }
                    }
                    if(face_verteces_number==4)
                    {
                        prism_volume += volume_of_4_tetraheda_from_4_verteces_and_center(verteces, face_center, element_center);
                    }
                    else
                    { //face_verteces_number==3
                        real m[4][dim];
                        for (int face_vert_j = 0;face_vert_j < 3;face_vert_j++) 
                        { 
                            int vert_j = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                            for (int k = 0;k < dim;++k)
                            { 
                                m[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                            }
                        }
                        for (int k = 0;k < dim;++k)
                        {
                            m[3][k] = element_center[k];
                        }

                        prism_volume += get_tetrahedron_volume(m);
                    }

                }
                mesh().Vol(i,0) = prism_volume;
                

            }
            if (elem_type == 7) 
            {
                //pyramid
                int number_of_faces = elem_ref().get_faces_n(elem_type);
                int verteces_indexes[4];
                int quad_face_number=-1;
                t_vec face_center(real(0.f), real(0.f), real(0.f));

                for (int face_j=0;face_j<number_of_faces;face_j++)
                {
                    int face_verteces_number=elem_ref().get_face_verts_n(elem_type,face_j);
                    if(face_verteces_number==4)
                    {
                        quad_face_number=face_j;
                        for (int k = 0;k < dim;++k)
                        { 
                            face_center[k]=mesh().center_faces(i,face_j,k);
                        } 
                        for (int face_vert_j = 0;face_vert_j < face_verteces_number;face_vert_j++) 
                        { 
                            verteces_indexes[face_vert_j] = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j);
                        }
                    }
                }
                int vertex_pinacle=-1;
                for (int face_j=0;face_j<number_of_faces;face_j++)
                {
                    int face_verteces_number=elem_ref().get_face_verts_n(elem_type,face_j); 
                    if(face_verteces_number==3)
                    {
                        for (int face_vert_j = 0;face_vert_j < 3;face_vert_j++) 
                        { 
                            int vertex_index = elem_ref().get_face_vert_i(elem_type,face_j,face_vert_j); 
                            for(int k=0;k<4;k++)
                            {
                                if(verteces_indexes[k]!=vertex_index)
                                {
                                    vertex_pinacle=vertex_index;
                                    break;
                                }
                            }
                        }
                    }
                    if(vertex_pinacle!=-1)
                        break;

                }
                t_vec  vertex_pinacle_coords(real(0.f), real(0.f), real(0.f));
                t_vec  verteces[4]=
                    {
                        (real(0.f), real(0.f), real(0.f)),
                        (real(0.f), real(0.f), real(0.f)),
                        (real(0.f), real(0.f), real(0.f)),
                        (real(0.f), real(0.f), real(0.f))
                    };
                
                for (int k = 0;k < dim;++k)
                { 
                    vertex_pinacle_coords[k] = mesh().vertexes(i, vertex_pinacle, k);
                }

                for (int face_vert_j = 0;face_vert_j < 4;face_vert_j++) 
                { 
                    int vert_j = elem_ref().get_face_vert_i(elem_type,quad_face_number,face_vert_j);
                    for (int k = 0;k < dim;++k)
                    { 
                        verteces[face_vert_j][k] = mesh().vertexes(i, vert_j, k);
                    }
                }


                mesh().Vol(i,0) = volume_of_4_tetraheda_from_4_verteces_and_center(verteces, face_center, vertex_pinacle_coords);



            }
        }
    };

    struct update_center_neighbour
    {
        __DEVICE_TAG__ void operator()(const int &i)const
        {
            int elem_type = mesh().get_elem_type(i);
            for (int face_i = 0;face_i < elem_ref().get_faces_n(elem_type);++face_i) 
            {
                int     nb = mesh().Neighbour(i, face_i);
                if (nb != CUDA_EMPTY_IDX) 
                {
                    t_vec   nb_center;
                    mesh().center.getv(nb, nb_center);
                    mesh().center_neighbour.setv(i, face_i, nb_center);
                }
            }
        }
    };

};

void init_mesh_deform_shepard
(
    const t_map_nodes &map_nodes, const t_cpu_mesh &cpu_mesh, 
    const std::set<int> &bnd_ids,
    t_vec_field &deform_vec,
    t_idx_field &buffer_idx,
    t_vec_field &deform_vec_buffer,
    t_idx_field &deform_idx_buffer
)
{
    int bnd_nodes_n = 0;
    for (int i_loc = map_nodes.min_loc_ind();i_loc <= map_nodes.max_loc_ind();++i_loc) 
    {
        if (!map_nodes.check_loc_has_loc_ind(i_loc)) continue;
        int i_glob = map_nodes.loc2glob(i_loc);
        if (bnd_ids.find(cpu_mesh.nodes[i_glob].bnd_id) == bnd_ids.end()) continue;
        bnd_nodes_n++;
    }
    deform_vec_buffer.init(bnd_nodes_n);
    deform_idx_buffer.init(bnd_nodes_n);
    deform_vec.init(map_nodes.max_loc_ind() - map_nodes.min_loc_ind() + 1, map_nodes.min_loc_ind());
    buffer_idx.init(map_nodes.max_loc_ind() - map_nodes.min_loc_ind() + 1, map_nodes.min_loc_ind());

    int bnd_nodes_i = 0;
    typename t_idx_field::view_type buffer_idx_view(buffer_idx, false),
                                    deform_idx_buffer_view(deform_idx_buffer, false);
    for (int i_loc = map_nodes.min_loc_ind();i_loc <= map_nodes.max_loc_ind();++i_loc) 
    {
        if (!map_nodes.check_loc_has_loc_ind(i_loc)) continue;
        int i_glob = map_nodes.loc2glob(i_loc);
        
        if (bnd_ids.find(cpu_mesh.nodes[i_glob].bnd_id) == bnd_ids.end()) 
        {
            buffer_idx_view(i_loc) = -1;
            continue;
        }
        buffer_idx_view(i_loc) = bnd_nodes_i;
        deform_idx_buffer_view(bnd_nodes_i) = i_loc;
        bnd_nodes_i++;
    }
    buffer_idx_view.release(true);
    deform_idx_buffer_view.release(true);

    t_elem_reference        elem_ref_;
    COPY_TO_CONSTANT_BUFFER(elem_ref, elem_ref_);
}

//deform_vec is input field (contains deforamtion in boundary nodes) but also is used as buffer
void mesh_deform_shepard
(
    const t_vec_field &deform_vec,
    const t_gpu_mesh &gpu_mesh, 
    const t_idx_field &buffer_idx,
    const t_vec_field &deform_vec_buffer,
    const t_idx_field &deform_idx_buffer
)
{
    COPY_TO_CONSTANT_BUFFER(mesh, gpu_mesh);

    t_for_each_1d           for_each_1d;
#ifdef MESH_DEFORM_SHEPARD_CUDA_NODE
    for_each_1d.block_size = 256;
#endif
#ifdef MESH_DEFORM_SHEPARD_OPENMP_NODE
    for_each_1d.threads_num = -1;
#endif

    //copy boundary deformations to separate buffer
    //put new coords to vertex array
    //update geometry features
    for_each_1d( calc_center(), 0, gpu_mesh.n_cv );
    for_each_1d( calc_center_faces(), 0, gpu_mesh.n_cv );
    for_each_1d( calc_norm(), 0, gpu_mesh.n_cv );
    //for_each_1d( calc_faces_S(), 0, gpu_mesh.n_cv );
    for_each_1d( calc_vol(), 0, gpu_mesh.n_cv );
    //TODO we need to sync all updated geometry features between processors 
    //(for those one which are stored not only for own elements, like element centers)
    for_each_1d( update_center_neighbour(), 0, gpu_mesh.n_cv );
}

void copy_deform_2_cpu_mesh
(
    const t_gpu_mesh    &gpu_mesh,
    const t_map_elems   &map_elems,
    const t_map_nodes   &map_nodes,
    t_cpu_mesh          &cpu_mesh
)
{
    t_tensor2_field_view_tml<real,vert_max_n,dim,tf_storage>         vertexes_view(gpu_mesh.vertexes, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        for (int j = 0;j < cpu_mesh.cv[i_glob].vert_n;++j)
            vertexes_view.getv(i_loc, j, cpu_mesh.cv[i_glob].vertexes[j]);
    }
    vertexes_view.release(false);

    //ISSUE do we need to copy stencil element's centers as well??
    t_tensor1_field_view_tml<real,dim,tf_storage>                    center_view(gpu_mesh.center, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        center_view.getv(i_loc, cpu_mesh.cv[i_glob].center);
    }
    center_view.release(false);

    t_tensor2_field_view_tml<real,faces_max_n,dim,tf_storage>        center_faces_view(gpu_mesh.center_faces, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            center_faces_view.getv(i_loc, j, cpu_mesh.cv[i_glob].face_centers[j]);
    }
    center_faces_view.release(false);

    t_tensor2_field_view_tml<real,faces_max_n,dim,tf_storage>        norm_view(gpu_mesh.Norm, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            norm_view.getv(i_loc, j, cpu_mesh.cv[i_glob].norms[j]);
    }
    norm_view.release(false);

    t_tensor1_field_view_tml<real,faces_max_n,tf_storage>            faces_S_view(gpu_mesh.faces_S, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            cpu_mesh.cv[i_glob].S[j] = faces_S_view(i_loc, j);
    }
    faces_S_view.release(false);

    t_tensor1_field_view_tml<real,1,tf_storage>                      vol_view(gpu_mesh.Vol, true);
    for (int i_ = 0;i_ < map_elems.get_size();++i_) 
    {
        int i_glob = map_elems.own_glob_ind(i_),
            i_loc = map_elems.own_loc_ind(i_);
        cpu_mesh.cv[i_glob].vol = vol_view(i_loc, 0);
    }
    vol_view.release(false);

    //TODO nodes data
}

}  /// namespace detail
}  /// namespace mesh
}  /// namespace scfd

#endif
