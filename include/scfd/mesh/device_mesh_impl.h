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

#ifndef __SCFD_MESH_DEVICE_MESH_IMPL_H__
#define __SCFD_MESH_DEVICE_MESH_IMPL_H__

#include "device_mesh.h"
#include "device_mesh_funcs.h"
#include "gmsh_mesh_elem_reference.h"

//TODO do something with possible simultanious several types INSTANTIATE (mesh and elem_ref names conflict)
//maybe create some kind of DEFINE_TEMPLATE_CONSTANT_BUFFER with explicit template params argumnet
#define SCFD_DEVICE_MESH_INSTANTIATE(T,MEMORY,DIM,ORD)                            \
  namespace scfd { namespace mesh { namespace detail {                            \
    DEFINE_CONSTANT_BUFFER((device_mesh<T,MEMORY,DIM,ORD>), mesh)                 \
    DEFINE_CONSTANT_BUFFER((gmsh_mesh_elem_reference<T>), elem_ref)               \
  } } }                                                                           \
  template class scfd::mesh::device_mesh<T,MEMORY,DIM,ORD>;

namespace scfd
{
namespace mesh
{

template<class T,class Memory,int Dim,class Ord>
template<class BasicMesh,class MapElems,class MapFaces,class MapNodes,class ForEach>
void    device_mesh<T,Memory,Dim,Ord>::init_elems_data
(
    const host_mesh<BasicMesh> &cpu_mesh,
    const MapElems &map_e, const MapFaces &map_f, const MapNodes &map_n,
    const ForEach &for_each
)
{
    using vec_t = static_vec::vec<T,dim>;
    using host_mesh_t = host_mesh<BasicMesh>;
    using host_ordinal = host_mesh_t::ordinal_type;
    using elem_reference_t = gmsh_mesh_elem_reference<T>;

    //n_cv = _n_cv; n_cv_all = _n_cv_all; i0 = _i0;
    own_elems_range.n = map_e.get_size();
    own_elems_range.i0 = 0;
    elems_range.n = map_e.max_loc_ind() - map_e.min_loc_ind() + 1;
    elems_range.i0 = map_e.min_loc_ind();

    ordinal_type    max_faces_n = cpu_mesh.get_elems_max_faces_num(),
                    max_nodes_n = cpu_mesh.get_elems_max_nodes_num(),
                    max_prim_nodes_n = cpu_mesh.get_elems_max_prim_nodes_num();

    elems_types.init(own_elems_range.n, own_elems_range.i0);
    elems_centers.init(elems_range.n, elems_range.i0);  //ISSUE
    elems_neighbours0_centers.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_faces_centers.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_vertexes.init(elems_range.n,max_nodes_n,elems_range.i0,0);
    elems_neighbours0.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_neighbours0_loc_face_i.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_faces_group_ids.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_group_ids.init(own_elems_range.n,own_elems_range.i0);
    elems_faces_norms.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_faces_areas.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
    elems_vols.init(elems_range.n, elems_range.i0);
    elems_prim_nodes_ids.init(own_elems_range.n,max_prim_nodes_n);
    elems_nodes_ids.init(own_elems_range.n,max_nodes_n);

    //TODO why only local?
    auto                    elem_type_view = elems_types.create_view(false);
    for(Ord i_ = 0;i_ < map_e.get_size();++i_) 
    {
        int     i_glob = map_e.own_glob_ind(i_),
                i = map_e.own_loc_ind(i_);
        elem_type_view(i) = cpu_mesh.get_elem_type(i_glob);
    }
    elem_type_view.release();

    is_homogeneous = cpu_mesh.is_homogeneous();
    homogeneous_elem_type = cpu_mesh.homogeneous_elem_type();

    /*auto                    center_view = elems_centers.create_view(false);
    for (Ord i = center_view.begin();i < center_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        center_view.setv(i, cpu_mesh.cv[i_glob].elems_centers);
    }
    center_view.release();*/

    /*auto        center_neighbour_view = elems_neighbours0_centers.create_view(false);
    for (Ord i = center_neighbour_view.begin();i < center_neighbour_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            if (cpu_mesh.cv[i_glob].neighbours[j] != -1) 
                center_neighbour_view.setv(i,j,cpu_mesh.cv[cpu_mesh.cv[i_glob].neighbours[j]].elems_centers);
    }
    center_neighbour_view.release();*/

    /*auto        center_faces_view = elems_faces_centers.create_view(false);
    for (Ord i = center_faces_view.begin();i < center_faces_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            center_faces_view.setv(i,j,cpu_mesh.cv[i_glob].face_centers[j]);
    }
    center_faces_view.release();*/

    auto         vertexes_view = elems_vertexes.create_view(false);
    for (Ord i = map_e.min_loc_ind();i <= map_e.max_loc_ind();++i) 
    {
        if (!map_e.check_loc_has_loc_ind(i)) continue;
        int i_glob = map_e.loc2glob(i);

        host_ordinal     elem_nodes[cpu_mesh.get_elems_max_nodes_num()];
        host_ordinal     nodes_n;
        cpu_mesh.get_elem_nodes(i, elem_nodes, &nodes_n);

        for (Ord j = 0;j < nodes_n;++j)
        {
            vec_t   vertex = cpu_mesh.get_node_coords(elem_nodes[j]);
            vertexes_view.setv(i,j,vertex);
        }
    }
    vertexes_view.release();

    /*auto                      vol_view = elems_vols.create_view(false);
    for (Ord i = vol_view.begin();i < vol_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        vol_view(i,0) = cpu_mesh.cv[i_glob].vol;
    }
    vol_view.release();*/

    /*auto            faces_S_view = elems_faces_areas.create_view(false);
    for (Ord i = faces_S_view.begin();i < faces_S_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            faces_S_view(i,j) = cpu_mesh.cv[i_glob].S[j];
    }
    faces_S_view.release();*/

    //TODO add cv index-accessors

    /*auto        norm_view = elems_faces_norms.create_view(false);
    for(Ord i = norm_view.begin();i < norm_view.end();i++) {
        int i_glob = map_e.loc2glob(i);
        for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
            norm_view.setv(i,j,cpu_mesh.cv[i_glob].norms[j]);
    }
    norm_view.release();*/

    auto          neighbours_view = elems_neighbours0.create_view(false);
    for (Ord i_ = 0;i_ < map_e.get_size();++i_) 
    {
        host_ordinal    i_glob = map_e.own_glob_ind(i_);
        Ord             i = map_e.own_loc_ind(i_);
        host_ordinal    neibs[cpu_mesh.get_elems_max_faces_num()];
        cpu_mesh.get_elem_neighbours0(i_glob, neibs);
        for (Ord j = 0;j < cpu_mesh.get_elem_faces_num(i_glob);++j) 
        {
            if (neibs[j] != host_mesh_t::special_id) 
                neighbours_view(i,j) = map_e.glob2loc(neibs[j]); 
            else 
                neighbours_view(i,j) = special_id;
        }
    }
    neighbours_view.release();

    auto          neighbours_loc_iface_view = elems_neighbours0_loc_face_i.create_view(false);
    for (Ord i_ = 0;i_ < map_e.get_size();++i_) 
    {
        host_ordinal    i_glob = map_e.own_glob_ind(i_);
        Ord             i = map_e.own_loc_ind(i_);
        Ord             loc_face_i[cpu_mesh.get_elems_max_faces_num()];
        cpu_mesh.get_elem_neighbours0_loc_face_i(i_glob, loc_face_i);
        for (Ord j = 0;j < cpu_mesh.get_elem_faces_num(i_glob);++j) 
        {
            neighbours_loc_iface_view(i,j) = loc_face_i[j];
        }
    }
    neighbours_loc_iface_view.release();

    auto          boundaries_view = elems_faces_group_ids.create_view(false);
    for (Ord i_ = 0;i_ < map_e.get_size();++i_) 
    {
        host_ordinal    i_glob = map_e.own_glob_ind(i_);
        Ord             i = map_e.own_loc_ind(i_);
        for (Ord j = 0;j < cpu_mesh.get_elem_faces_num(i_glob);++j)
        { 
            boundaries_view(i,j) = cpu_mesh.get_elem_face_group_id(i_glob,j);
        }
    }
    boundaries_view.release();

    auto                    vol_id_view = elems_group_ids.create_view(false);
    for (Ord i_ = 0;i_ < map_e.get_size();++i_) 
    {
        host_ordinal    i_glob = map_e.own_glob_ind(i_);
        Ord             i = map_e.own_loc_ind(i_);
        vol_id_view(i) = cpu_mesh.get_elem_group_id(i_glob);
    }
    vol_id_view.release();

    auto   elem_node_ids_view = elems_prim_nodes_ids.create_view(false);
    for(Ord i_ = 0;i_ < map_e.get_size();i_++) 
    {
        int     i_glob = map_e.own_glob_ind(i_),
            i_loc = map_e.own_loc_ind(i_);
        for (Ord vert_i = 0;vert_i < cpu_mesh.cv[i_glob].vert_n;++vert_i) 
        {
            elem_node_ids_view(i_loc,vert_i) = map_n.glob2loc( cpu_mesh.cv_2_node_ids[i_glob].ids[vert_i] );
        }
    }
    elem_node_ids_view.release();

    elem_reference_t        elem_ref_;
    COPY_TO_CONSTANT_BUFFER(elem_ref, elem_ref_);
    COPY_TO_CONSTANT_BUFFER(mesh, *this);

    //copy boundary deformations to separate buffer
    //put new coords to vertex array
    //update geometry features
    for_each( device_mesh_funcs_t::calc_center(), 0, gpu_mesh.n_cv );
    for_each( device_mesh_funcs_t::calc_center_faces(), 0, gpu_mesh.n_cv );
    for_each( device_mesh_funcs_t::calc_norm(), 0, gpu_mesh.n_cv );
    //for_each_1d( calc_faces_S(), 0, gpu_mesh.n_cv );
    for_each( device_mesh_funcs_t::calc_vol(), 0, gpu_mesh.n_cv );
    //TODO we need to sync all updated geometry features between processors 
    //(for those one which are stored not only for own elements, like element centers)
    for_each( device_mesh_funcs_t::update_center_neighbour(), 0, gpu_mesh.n_cv );
}

template<class T,class Memory,int Dim,class Ord>
template<class BasicMesh,class MapElems,class MapFaces,class MapNodes,class ForEach>
void    device_mesh<T,Memory,Dim,Ord>::init_nodes_data
(
    const host_mesh<BasicMesh> &cpu_mesh,
    const MapElems &map_e, const MapFaces &map_f, const MapNodes &map_n,
    const ForEach &for_each
)
{
    own_nodes_range.n = map_n.get_size();
    own_nodes_range.i0 = 0;
    nodes_range.n = map_n.max_loc_ind() - map_n.min_loc_ind() + 1;
    nodes_range.i0 = map_n.min_loc_ind();

    nodes_coords.init(own_nodes_range.n,own_nodes_range.i0);
    node_vol_id.init(own_nodes_range.n,own_nodes_range.i0);
    //node_bnd_id.init(n_nodes);

    auto                 node_coords_view = nodes_coords.create_view(false);
    for(Ord i_ = 0;i_ < map_n.get_size();++i_) 
    {
        int     i_glob = map_n.own_glob_ind(i_),
            i_loc = map_n.own_loc_ind(i_);
        node_coords_view.setv(i_loc, cpu_mesh.nodes[i_glob].c);
    }
    node_coords_view.release();

    auto                   node_vol_id_view = node_vol_id.create_view(false);
    for(Ord i_ = 0;i_ < map_n.get_size();++i_) 
    {
        int     i_glob = map_n.own_glob_ind(i_),
            i_loc = map_n.own_loc_ind(i_);
        node_vol_id_view(i_loc) = cpu_mesh.nodes[i_glob].vol_id;
    }
    node_vol_id_view.release();

    /*auto                   node_bnd_id_view = node_bnd_id.create_view(false);
    for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
        int     i_glob = map_n.own_glob_ind(i_),
            i_loc = map_n.own_loc_ind(i_);
        node_bnd_id_view(i_loc) = cpu_mesh.nodes[i_glob].bnd_id;
    }
    node_bnd_id_view.release();*/

    node_2_elem_graph_refs.init(own_nodes_range.n);
    node_2_elem_graph_sz = 0;
    for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
        int     i_node_glob = map_n.own_glob_ind(i_);
        node_2_elem_graph_sz += cpu_mesh.node_2_cv_ids_ref[i_node_glob].second - cpu_mesh.node_2_cv_ids_ref[i_node_glob].first;
    }
    node_2_elem_graph_elem_ids.init(node_2_elem_graph_sz);
    node_2_elem_graph_node_ids.init(node_2_elem_graph_sz);

    auto            node_2_elem_graph_refs_view = node_2_elem_graph_refs.create_view(false);
    auto            node_2_elem_graph_elem_ids_view = node_2_elem_graph_elem_ids.create_view(false);
    auto            node_2_elem_graph_node_ids_view = node_2_elem_graph_node_ids.create_view(false);
    Ord             curr_graph_loc_idx = 0;
    for(Ord i_ = 0;i_ < map_n.get_size();++i_) 
    {
        int     i_node_glob = map_n.own_glob_ind(i_),
            i_node_loc = map_n.own_loc_ind(i_);
        node_2_elem_graph_refs_view(i_node_loc, 0) = curr_graph_loc_idx;
        for (int ii = cpu_mesh.node_2_cv_ids_ref[i_node_glob].first;ii < cpu_mesh.node_2_cv_ids_ref[i_node_glob].second;++ii) 
        {
            node_2_elem_graph_elem_ids_view(curr_graph_loc_idx) = map_e.glob2loc( cpu_mesh.node_2_cv_ids_data[ii] );
            node_2_elem_graph_node_ids_view(curr_graph_loc_idx) = i_node_loc;
            curr_graph_loc_idx++;
        }
        node_2_elem_graph_refs_view(i_node_loc, 1) = curr_graph_loc_idx;
    }
    node_2_elem_graph_refs_view.release();
    node_2_elem_graph_elem_ids_view.release();
    node_2_elem_graph_node_ids_view.release();
}

}  /// namespace mesh
}  /// namespace scfd

#endif
