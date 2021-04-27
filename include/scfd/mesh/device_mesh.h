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

#ifndef __SCFD_MESH_DEVICE_MESH_H__
#define __SCFD_MESH_DEVICE_MESH_H__

#include <scfd/arrays/tensorN_array.h>

//TODO rename these macroses somehow
#ifndef  GPU_MESH_DEFAULT_MAX_FACES_N
#define  GPU_MESH_DEFAULT_MAX_FACES_N 4
#endif
#ifndef  GPU_MESH_DEFAULT_MAX_VERT_N
#define  GPU_MESH_DEFAULT_MAX_VERT_N  4
#endif

//TODO
#define CUDA_EMPTY_IDX  -1000000000

namespace scfd
{
namespace mesh
{

template<class Ord>
struct index_range_descr
{
    Ord i0, n;
};

//map_n (nodes map) and map_e (element map) during all calls must be the same 
//ISSUE perhaps, better to make copy of maps inside (but what to do with gpu in this case?)

//ISSUE should not we make own flag?? (because this structure is copied to const memory)

//TODO for now Ord = SCFD_ARRAYS_ORDINAL_TYPE is only supported (arrays don't have Ord template parameter for now)

template<class T,class Memory,int Dim = 3,class Ord = SCFD_ARRAYS_ORDINAL_TYPE>
struct device_mesh
{
    using scalar_type = T;
    using memory_type = Memory;
    static const int dim = Dim;
    using ordinal_type = Ord;
    //using glob_ordinal_type = int;
    using elem_type_ordinal_type = int;
    using index_range_descr_type = index_range_descr<Ord>;

    using namespace arrays;

    //ISSUE neither n_cv
    //but i0, n_cv_all somehow does

    ordinal_type        max_faces_n, max_prim_nodes_n;

    /// Elements data part

    index_range_descr_type                          elems_range, 
                                                    own_elems_range;

    //vars buffers CONTAINS space for elements in [i0, i0+n_cv_all)
    //Ord                                             i0, n_cv_all;
    //gpu OWN elements in [0,n_cv)
    //Ord                                             n_cv;
    //if is_homogeneous == true then all elements in mesh have the same elem_type
    bool                                            is_homogeneous;
    elem_type_ordinal_type                          homogeneous_elem_type;  //valid only if is_homogeneous == true
    tensor1_array<elem_type_ordinal_type,Memory,1>  elems_types;              //valid only if is_homogeneous == false
    tensor1_array<T,Memory,Dim>                     elems_centers;
    tensor2_array<T,Memory,dyn_dim,Dim>             elems_neighbours0_centers;
    tensor2_array<T,Memory,dyn_dim,Dim>             elems_faces_centers;
    //TODO replace with var array
    tensor2_array<T,Memory,dyn_dim,Dim>             elems_vertexes;
    tensor1_array<Ord,Memory,dyn_dim>               elems_neighbours0;
    tensor1_array<Ord,Memory,dyn_dim>               elems_neighbours0_loc_face_i;
    tensor1_array<Ord,Memory,dyn_dim>               elems_faces_group_ids;
    tensor1_array<Ord,Memory,1>                     elems_group_ids;
    tensor2_array<T,Memory,dyn_dim,Dim>             elems_faces_norms;
    tensor1_array<T,Memory,dyn_dim>                 elems_faces_areas;
    tensor1_array<T,Memory,1>                       elems_vols;
    //elements to nodes graphs
    tensor1_array<Ord,Memory,dyn_dim>               elems_prim_nodes_ids;
    //TODO replace with var array
    tensor1_array<Ord,Memory,dyn_dim>               elems_nodes_ids;

    /// Nodes data part

    index_range_descr_type                          nodes_range, 
                                                    own_nodes_range;
    //Ord                                             i0_nodes, n_nodes_all;
    //Ord                                             n_nodes;
    tensor1_array<T,Memory,Dim>                     nodes_coords;
    tensor0_array<Ord,Memory>                       nodes_group_ids;
    //tensor0_array<Ord,Memory>                       node_bnd_id;

    //nodes to elements graph part
    Ord                                             node_2_elem_graph_sz;
    tensor1_array<Ord,Memory,2>                     node_2_elem_graph_refs;
    tensor0_array<Ord,Memory>                       node_2_elem_graph_elem_ids;
    tensor0_array<Ord,Memory>                       node_2_elem_graph_node_ids;

    /// Faces data part

    index_range_descr_type                          faces_range, 
                                                    own_faces_range;

    __DEVICE_TAG__ elem_type_ordinal_type  get_elem_type(Ord i)const
    {
        if (is_homogeneous) return homogeneous_elem_type; else return elems_types(i,0);
    }

    //void  init(int _n_cv, int _n_cv_all, int _i0, int _maxQ_real)
    template<class MapElems>
    void    init_elems(const MapElems &map_e)
    {
        
    }
    template<class MapNodes>
    void    init_nodes(const MapNodes &map_n)
    {
        own_nodes_range.n = map_n.get_size();
        own_nodes_range.i0 = 0;
        nodes_range.n = map_n.max_loc_ind() - map_n.min_loc_ind() + 1;
        nodes_range.i0 = map_n.min_loc_ind();

        nodes_coords.init(own_nodes_range.n,own_nodes_range.i0);
        node_vol_id.init(own_nodes_range.n,own_nodes_range.i0);
        //node_bnd_id.init(n_nodes);
    }
    template<class MapElems,class MapNodes,class BasicMesh>
    void    init_node_2_elem_graph(const MapElems &map_e, const MapNodes &map_n, host_mesh<BasicMesh> &cpu_mesh)
    {
        node_2_elem_graph_refs.init(own_nodes_range.n);
        node_2_elem_graph_sz = 0;
        for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_node_glob = map_n.own_glob_ind(i_);
            node_2_elem_graph_sz += cpu_mesh.node_2_cv_ids_ref[i_node_glob].second - cpu_mesh.node_2_cv_ids_ref[i_node_glob].first;
        }
        node_2_elem_graph_elem_ids.init(node_2_elem_graph_sz);
        node_2_elem_graph_node_ids.init(node_2_elem_graph_sz);
    }

    //TODO fix index calculation type (enumerate through map, not through view range)
    template<class MapElems,class BasicMesh>
    void    init_elems_data(const MapElems &map_e, host_mesh<BasicMesh> &cpu_mesh)
    {
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
        elems_vertexes.init(own_elems_range.n,max_nodes_n,own_elems_range.i0,0);
        elems_neighbours0.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
        elems_neighbours0_loc_face_i.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
        elems_faces_group_ids.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
        elems_group_ids.init(own_elems_range.n,own_elems_range.i0);
        elems_faces_norms.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
        elems_faces_areas.init(own_elems_range.n,max_faces_n,own_elems_range.i0,0);
        elems_vols.init(elems_range.n, elems_range.i0);
        elems_prim_nodes_ids.init(own_elems_range.n,max_prim_nodes_n);
        elems_nodes_ids.init(own_elems_range.n,max_nodes_n);

        auto                    center_view = elems_centers.create_view(false);
        for (Ord i = center_view.begin();i < center_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            center_view.setv(i, cpu_mesh.cv[i_glob].elems_centers);
        }
        center_view.release();

        auto        center_neighbour_view = elems_neighbours0_centers.create_view(false);
        for (Ord i = center_neighbour_view.begin();i < center_neighbour_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) 
                    center_neighbour_view.setv(i,j,cpu_mesh.cv[cpu_mesh.cv[i_glob].neighbours[j]].elems_centers);
        }
        center_neighbour_view.release();

        auto        center_faces_view = elems_faces_centers.create_view(false);
        for (Ord i = center_faces_view.begin();i < center_faces_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                center_faces_view.setv(i,j,cpu_mesh.cv[i_glob].face_centers[j]);
        }
        center_faces_view.release();

        auto         vertexes_view = elems_vertexes.create_view(false);
        for (Ord i = vertexes_view.begin();i < vertexes_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].vert_n;++j)
                vertexes_view.setv(i,j,cpu_mesh.cv[i_glob].elems_vertexes[j]);
        }
        vertexes_view.release();

        auto                      vol_view = elems_vols.create_view(false);
        for (Ord i = vol_view.begin();i < vol_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            vol_view(i,0) = cpu_mesh.cv[i_glob].vol;
        }
        vol_view.release();

        auto            faces_S_view = elems_faces_areas.create_view(false);
        for (Ord i = faces_S_view.begin();i < faces_S_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                faces_S_view(i,j) = cpu_mesh.cv[i_glob].S[j];
        }
        faces_S_view.release();

        //TODO add cv index-accessors

        auto        norm_view = elems_faces_norms.create_view(false);
        for(Ord i = norm_view.begin();i < norm_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                norm_view.setv(i,j,cpu_mesh.cv[i_glob].norms[j]);
        }
        norm_view.release();

        auto          neighbours_view = elems_neighbours0.create_view(false);
        for(Ord i = neighbours_view.begin();i < neighbours_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) neighbours_view(i,j) = map_e.glob2loc(cpu_mesh.cv[i_glob].neighbours[j]); else neighbours_view(i,j) = CUDA_EMPTY_IDX;
        }
        neighbours_view.release();

        auto          neighbours_loc_iface_view = elems_neighbours0_loc_face_i.create_view(false);
        for(Ord i = neighbours_loc_iface_view.begin();i < neighbours_loc_iface_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                neighbours_loc_iface_view(i,j) = cpu_mesh.cv[i_glob].neighbours_loc_iface[j];
        }
        neighbours_loc_iface_view.release();

        auto          boundaries_view = elems_faces_group_ids.create_view(false);
        for(Ord i = boundaries_view.begin();i < boundaries_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                boundaries_view(i,j) = cpu_mesh.cv[i_glob].boundaries[j];
        }
        boundaries_view.release();

        auto                    vol_id_view = elems_group_ids.create_view(false);
        for(Ord i = vol_id_view.begin();i < vol_id_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            vol_id_view(i,0) = cpu_mesh.cv[i_glob].vol_id;
        }
        vol_id_view.release();

        auto                    elem_type_view = elems_types.create_view(false);
        for(Ord i = elem_type_view.begin();i < elem_type_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            elem_type_view(i,0) = cpu_mesh.cv[i_glob].elems_types;
        }
        elem_type_view.release();

        is_homogeneous = cpu_mesh.is_homogeneous;
        homogeneous_elem_type = cpu_mesh.homogeneous_elem_type;
    }
    template<class MapNodes,class BasicMesh>
    void    init_nodes_data(const MapNodes &map_n, host_mesh<BasicMesh> &cpu_mesh)
    {
        auto                 node_coords_view = nodes_coords.create_view(false);
        for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_glob = map_n.own_glob_ind(i_),
                i_loc = map_n.own_loc_ind(i_);
            node_coords_view.setv(i_loc, cpu_mesh.nodes[i_glob].c);
        }
        node_coords_view.release();

        auto                   node_vol_id_view = node_vol_id.create_view(false);
        for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
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
    }
    template<class MapElems,class MapNodes,class BasicMesh>
    void    init_elem_node_ids_data(const MapElems &map_e, const MapNodes &map_n, host_mesh<BasicMesh> &cpu_mesh)
    {
        auto   elem_node_ids_view = elems_prim_nodes_ids.create_view(false);
        for(Ord i_ = 0;i_ < map_e.get_size();i_++) {
            int     i_glob = map_e.own_glob_ind(i_),
                i_loc = map_e.own_loc_ind(i_);
            for (Ord vert_i = 0;vert_i < cpu_mesh.cv[i_glob].vert_n;++vert_i) {
                elem_node_ids_view(i_loc,vert_i) = map_n.glob2loc( cpu_mesh.cv_2_node_ids[i_glob].ids[vert_i] );
            }
        }
        elem_node_ids_view.release();
    }
    template<class MapElems,class MapNodes,class BasicMesh>
    void    init_node_2_elem_graph_data(const MapElems &map_e, const MapNodes &map_n, host_mesh<BasicMesh> &cpu_mesh)
    {
        auto            node_2_elem_graph_refs_view = node_2_elem_graph_refs.create_view(false);
        auto              node_2_elem_graph_elem_ids_view = node_2_elem_graph_elem_ids.create_view(false);
        auto              node_2_elem_graph_node_ids_view = node_2_elem_graph_node_ids.create_view(false);
        Ord                                             curr_graph_loc_idx = 0;
        for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_node_glob = map_n.own_glob_ind(i_),
                i_node_loc = map_n.own_loc_ind(i_);
            node_2_elem_graph_refs_view(i_node_loc, 0) = curr_graph_loc_idx;
            for (int ii = cpu_mesh.node_2_cv_ids_ref[i_node_glob].first;ii < cpu_mesh.node_2_cv_ids_ref[i_node_glob].second;++ii) {
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

    //TODO change bane
    template<class MapElems,class BasicMesh>
    void    dump_elems_geom_data(const MapElems &map_e, host_mesh<BasicMesh> &cpu_mesh)const
    {
        auto                    center_view = elems_centers.create_view(true);
        for (Ord i = center_view.begin();i < center_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            center_view.getv(i, cpu_mesh.cv[i_glob].elems_centers);
        }
        center_view.release(false);

        auto        center_neighbour_view = elems_neighbours0_centers.create_view(true);
        for (Ord i = center_neighbour_view.begin();i < center_neighbour_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) 
                    center_neighbour_view.getv(i,j,cpu_mesh.cv[cpu_mesh.cv[i_glob].neighbours[j]].elems_centers);
        }
        center_neighbour_view.release(false);

        auto        center_faces_view = elems_faces_centers.create_view(true);
        for (Ord i = center_faces_view.begin();i < center_faces_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                center_faces_view.getv(i,j,cpu_mesh.cv[i_glob].face_centers[j]);
        }
        center_faces_view.release(false);

        auto         vertexes_view = elems_vertexes.create_view(true);
        for (Ord i = vertexes_view.begin();i < vertexes_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].vert_n;++j)
                vertexes_view.getv(i,j,cpu_mesh.cv[i_glob].elems_vertexes[j]);
        }
        vertexes_view.release(false);

        auto                      vol_view = elems_vols.create_view(true);
        for (Ord i = vol_view.begin();i < vol_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            cpu_mesh.cv[i_glob].vol = vol_view(i,0);
        }
        vol_view.release(false);

        auto            faces_S_view = elems_faces_areas.create_view(true);
        for (Ord i = faces_S_view.begin();i < faces_S_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                cpu_mesh.cv[i_glob].S[j] = faces_S_view(i,j);
        }
        faces_S_view.release(false);

        //TODO add cv index-accessors

        auto        norm_view = elems_faces_norms.create_view(true);
        for(Ord i = norm_view.begin();i < norm_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (Ord j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                norm_view.getv(i,j,cpu_mesh.cv[i_glob].norms[j]);
        }
        norm_view.release(false);
    }
    template<class MapNodes,class BasicMesh>
    void    dump_nodes_geom_data(const MapNodes &map_n, host_mesh<BasicMesh> &cpu_mesh)const
    {
        auto                    node_coords_view = nodes_coords.create_view(true);
        for(Ord i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_glob = map_n.own_glob_ind(i_),
                i_loc = map_n.own_loc_ind(i_);
            node_coords_view.getv(i_loc, cpu_mesh.nodes[i_glob].c);
        }
        node_coords_view.release(false);
    }

    ~device_mesh()
    {
    }
};

}  /// namespace mesh
}  /// namespace scfd

#endif
