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

#include <tensor_field/t_tensor_field_tml.h>

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

//GPU-oriented classes

//man_n (nodes map) and map_e (element map) during all calls must be the same 
//ISSUE perhaps, better to make cope of maps inside (but what to do with gpu in this case?)

//TODO should be renamed to device_mesh

//ISSUE should not we make own flag?? (because this structure is copied to const memory)
template<class T,int dim = 3, int max_faces_n = GPU_MESH_DEFAULT_MAX_FACES_N, int max_vert_n = GPU_MESH_DEFAULT_MAX_VERT_N, t_tensor_field_storage strg = TFS_DEVICE>
struct device_mesh
{
    //ISSUE neither n_cv
    //but i0, n_cv_all somehow does

    //elements data part

    //vars buffers CONTAINS space for elements in [i0, i0+n_cv_all)
    int                                             i0, n_cv_all;
    //gpu OWN elements in [0,n_cv)
    int                                             n_cv;
    //if is_homogeneous == true then all elements in mesh have the same elem_type
    bool                                            is_homogeneous;
    int                                             homogeneous_elem_type;  //valid only if is_homogeneous == true
    t_tensor1_field_tml<int,1,strg>                 elem_type;              //valid only if is_homogeneous == false
    t_tensor1_field_tml<T,dim,strg>                 center;
    t_tensor2_field_tml<T,max_faces_n,dim,strg>     center_neighbour;
    t_tensor2_field_tml<T,max_faces_n,dim,strg>     center_faces;
    t_tensor2_field_tml<T,max_vert_n,dim,strg>      vertexes;
    t_tensor1_field_tml<int,max_faces_n,strg>       Neighbour;
    t_tensor1_field_tml<int,max_faces_n,strg>       Neighbour_loc_iface;
    t_tensor1_field_tml<int,max_faces_n,strg>       Boundary;
    t_tensor1_field_tml<int,1,strg>                 Volume_id;
    t_tensor2_field_tml<T,max_faces_n,dim,strg>     Norm;
    t_tensor1_field_tml<T,max_faces_n,strg>         faces_S;
    t_tensor1_field_tml<T,1,strg>                   Vol;

    //nodes data part
    int                                             i0_nodes, n_nodes_all;
    int                                             n_nodes;
    t_tensor1_field_tml<T,dim,strg>                 node_coords;
    t_tensor0_field_tml<int,strg>                   node_vol_id;
    t_tensor0_field_tml<int,strg>                   node_bnd_id;

    //elements to nodes graph part
    t_tensor1_field_tml<int,max_vert_n,strg>        elem_node_ids;

    //nodes to elements graph part
    int                                             node_2_elem_graph_sz;
    t_tensor1_field_tml<int,2,strg>                 node_2_elem_graph_refs;
    t_tensor0_field_tml<int,strg>                   node_2_elem_graph_elem_ids;
    t_tensor0_field_tml<int,strg>                   node_2_elem_graph_node_ids;

    __DEVICE_TAG__ int      get_elem_type(int i)const
    {
        if (is_homogeneous) return homogeneous_elem_type; else return elem_type(i,0);
    }

    //void  init(int _n_cv, int _n_cv_all, int _i0, int _maxQ_real)
    template<class MAP_ELEMS>
    void    init_elems(const MAP_ELEMS &map_e)
    {
        //n_cv = _n_cv; n_cv_all = _n_cv_all; i0 = _i0;
        n_cv = map_e.get_size();
        n_cv_all = map_e.max_loc_ind() - map_e.min_loc_ind() + 1;
        i0 = map_e.min_loc_ind();

        elem_type.init(n_cv);
        center.init(n_cv_all, i0);  //ISSUE
        center_neighbour.init(n_cv);
        center_faces.init(n_cv);
        vertexes.init(n_cv);
        Neighbour.init(n_cv);
        Neighbour_loc_iface.init(n_cv);
        Boundary.init(n_cv);
        Volume_id.init(n_cv);
        Norm.init(n_cv);
        faces_S.init(n_cv);
        Vol.init(n_cv_all, i0);
    }
    template<class MAP_NODES>
    void    init_nodes(const MAP_NODES &map_n)
    {
        n_nodes = map_n.get_size();
        n_nodes_all = map_n.max_loc_ind() - map_n.min_loc_ind() + 1;
        i0_nodes = map_n.min_loc_ind();

        node_coords.init(n_nodes);
        node_vol_id.init(n_nodes);
        node_bnd_id.init(n_nodes);
    }
    //this called strictly after init_elems()
    void    init_elem_node_ids()
    {
        elem_node_ids.init(n_cv);
    }
    template<class MAP_ELEMS,class MAP_NODES,class CPU_MESH>
    void    init_node_2_elem_graph(const MAP_ELEMS &map_e, const MAP_NODES &map_n, CPU_MESH &cpu_mesh)
    {
        node_2_elem_graph_refs.init(n_nodes);
        node_2_elem_graph_sz = 0;
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_node_glob = map_n.own_glob_ind(i_);
            node_2_elem_graph_sz += cpu_mesh.node_2_cv_ids_ref[i_node_glob].second - cpu_mesh.node_2_cv_ids_ref[i_node_glob].first;
        }
        node_2_elem_graph_elem_ids.init(node_2_elem_graph_sz);
        node_2_elem_graph_node_ids.init(node_2_elem_graph_sz);
    }

    //TODO fix index calculation type (enumerate through map, not through view range)
    template<class MAP_ELEMS,class CPU_MESH>
    void    init_elems_data(const MAP_ELEMS &map_e, CPU_MESH &cpu_mesh)
    {
        t_tensor1_field_view_tml<T,dim,strg>                    center_view(center, false);
        for (int i = center_view.begin();i < center_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            center_view.setv(i, cpu_mesh.cv[i_glob].center);
        }
        center_view.release();

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        center_neighbour_view(center_neighbour, false);
        for (int i = center_neighbour_view.begin();i < center_neighbour_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) 
                    center_neighbour_view.setv(i,j,cpu_mesh.cv[cpu_mesh.cv[i_glob].neighbours[j]].center);
        }
        center_neighbour_view.release();

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        center_faces_view(center_faces, false);
        for (int i = center_faces_view.begin();i < center_faces_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                center_faces_view.setv(i,j,cpu_mesh.cv[i_glob].face_centers[j]);
        }
        center_faces_view.release();

        t_tensor2_field_view_tml<T,max_vert_n,dim,strg>         vertexes_view(vertexes, false);
        for (int i = vertexes_view.begin();i < vertexes_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].vert_n;++j)
                vertexes_view.setv(i,j,cpu_mesh.cv[i_glob].vertexes[j]);
        }
        vertexes_view.release();

        t_tensor1_field_view_tml<T,1,strg>                      vol_view(Vol, false);
        for (int i = vol_view.begin();i < vol_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            vol_view(i,0) = cpu_mesh.cv[i_glob].vol;
        }
        vol_view.release();

        t_tensor1_field_view_tml<T,max_faces_n,strg>            faces_S_view(faces_S, false);
        for (int i = faces_S_view.begin();i < faces_S_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                faces_S_view(i,j) = cpu_mesh.cv[i_glob].S[j];
        }
        faces_S_view.release();

        //TODO add cv index-accessors

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        norm_view(Norm, false);
        for(int i = norm_view.begin();i < norm_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                norm_view.setv(i,j,cpu_mesh.cv[i_glob].norms[j]);
        }
        norm_view.release();

        t_tensor1_field_view_tml<int,max_faces_n,strg>          neighbours_view(Neighbour, false);
        for(int i = neighbours_view.begin();i < neighbours_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) neighbours_view(i,j) = map_e.glob2loc(cpu_mesh.cv[i_glob].neighbours[j]); else neighbours_view(i,j) = CUDA_EMPTY_IDX;
        }
        neighbours_view.release();

        t_tensor1_field_view_tml<int,max_faces_n,strg>          neighbours_loc_iface_view(Neighbour_loc_iface, false);
        for(int i = neighbours_loc_iface_view.begin();i < neighbours_loc_iface_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                neighbours_loc_iface_view(i,j) = cpu_mesh.cv[i_glob].neighbours_loc_iface[j];
        }
        neighbours_loc_iface_view.release();

        t_tensor1_field_view_tml<int,max_faces_n,strg>          boundaries_view(Boundary, false);
        for(int i = boundaries_view.begin();i < boundaries_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                boundaries_view(i,j) = cpu_mesh.cv[i_glob].boundaries[j];
        }
        boundaries_view.release();

        t_tensor1_field_view_tml<int,1,strg>                    vol_id_view(Volume_id, false);
        for(int i = vol_id_view.begin();i < vol_id_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            vol_id_view(i,0) = cpu_mesh.cv[i_glob].vol_id;
        }
        vol_id_view.release();

        t_tensor1_field_view_tml<int,1,strg>                    elem_type_view(elem_type, false);
        for(int i = elem_type_view.begin();i < elem_type_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            elem_type_view(i,0) = cpu_mesh.cv[i_glob].elem_type;
        }
        elem_type_view.release();

        is_homogeneous = cpu_mesh.is_homogeneous;
        homogeneous_elem_type = cpu_mesh.homogeneous_elem_type;
    }
    template<class MAP_NODES,class CPU_MESH>
    void    init_nodes_data(const MAP_NODES &map_n, CPU_MESH &cpu_mesh)
    {
        t_tensor1_field_view_tml<T,dim,strg>                 node_coords_view(node_coords, false);
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_glob = map_n.own_glob_ind(i_),
                i_loc = map_n.own_loc_ind(i_);
            node_coords_view.setv(i_loc, cpu_mesh.nodes[i_glob].c);
        }
        node_coords_view.release();

        t_tensor0_field_view_tml<int,strg>                   node_vol_id_view(node_vol_id, false);
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_glob = map_n.own_glob_ind(i_),
                i_loc = map_n.own_loc_ind(i_);
            node_vol_id_view(i_loc) = cpu_mesh.nodes[i_glob].vol_id;
        }
        node_vol_id_view.release();

        t_tensor0_field_view_tml<int,strg>                   node_bnd_id_view(node_bnd_id, false);
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
            int     i_glob = map_n.own_glob_ind(i_),
                i_loc = map_n.own_loc_ind(i_);
            node_bnd_id_view(i_loc) = cpu_mesh.nodes[i_glob].bnd_id;
        }
        node_bnd_id_view.release();
    }
    template<class MAP_ELEMS,class MAP_NODES,class CPU_MESH>
    void    init_elem_node_ids_data(const MAP_ELEMS &map_e, const MAP_NODES &map_n, CPU_MESH &cpu_mesh)
    {
        t_tensor1_field_view_tml<int,max_vert_n,strg>   elem_node_ids_view(elem_node_ids, false);
        for(int i_ = 0;i_ < map_e.get_size();i_++) {
            int     i_glob = map_e.own_glob_ind(i_),
                i_loc = map_e.own_loc_ind(i_);
            for (int vert_i = 0;vert_i < cpu_mesh.cv[i_glob].vert_n;++vert_i) {
                elem_node_ids_view(i_loc,vert_i) = map_n.glob2loc( cpu_mesh.cv_2_node_ids[i_glob].ids[vert_i] );
            }
        }
        elem_node_ids_view.release();
    }
    template<class MAP_ELEMS,class MAP_NODES,class CPU_MESH>
    void    init_node_2_elem_graph_data(const MAP_ELEMS &map_e, const MAP_NODES &map_n, CPU_MESH &cpu_mesh)
    {
        t_tensor1_field_view_tml<int,2,strg>            node_2_elem_graph_refs_view(node_2_elem_graph_refs, false);
        t_tensor0_field_view_tml<int,strg>              node_2_elem_graph_elem_ids_view(node_2_elem_graph_elem_ids, false);
        t_tensor0_field_view_tml<int,strg>              node_2_elem_graph_node_ids_view(node_2_elem_graph_node_ids, false);
        int                                             curr_graph_loc_idx = 0;
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
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
    template<class MAP_ELEMS,class CPU_MESH>
    void    dump_elems_geom_data(const MAP_ELEMS &map_e, CPU_MESH &cpu_mesh)const
    {
        t_tensor1_field_view_tml<T,dim,strg>                    center_view(center, true);
        for (int i = center_view.begin();i < center_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            center_view.getv(i, cpu_mesh.cv[i_glob].center);
        }
        center_view.release(false);

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        center_neighbour_view(center_neighbour, true);
        for (int i = center_neighbour_view.begin();i < center_neighbour_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                if (cpu_mesh.cv[i_glob].neighbours[j] != -1) 
                    center_neighbour_view.getv(i,j,cpu_mesh.cv[cpu_mesh.cv[i_glob].neighbours[j]].center);
        }
        center_neighbour_view.release(false);

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        center_faces_view(center_faces, true);
        for (int i = center_faces_view.begin();i < center_faces_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                center_faces_view.getv(i,j,cpu_mesh.cv[i_glob].face_centers[j]);
        }
        center_faces_view.release(false);

        t_tensor2_field_view_tml<T,max_vert_n,dim,strg>         vertexes_view(vertexes, true);
        for (int i = vertexes_view.begin();i < vertexes_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].vert_n;++j)
                vertexes_view.getv(i,j,cpu_mesh.cv[i_glob].vertexes[j]);
        }
        vertexes_view.release(false);

        t_tensor1_field_view_tml<T,1,strg>                      vol_view(Vol, true);
        for (int i = vol_view.begin();i < vol_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            cpu_mesh.cv[i_glob].vol = vol_view(i,0);
        }
        vol_view.release(false);

        t_tensor1_field_view_tml<T,max_faces_n,strg>            faces_S_view(faces_S, true);
        for (int i = faces_S_view.begin();i < faces_S_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                cpu_mesh.cv[i_glob].S[j] = faces_S_view(i,j);
        }
        faces_S_view.release(false);

        //TODO add cv index-accessors

        t_tensor2_field_view_tml<T,max_faces_n,dim,strg>        norm_view(Norm, true);
        for(int i = norm_view.begin();i < norm_view.end();i++) {
            int i_glob = map_e.loc2glob(i);
            for (int j = 0;j < cpu_mesh.cv[i_glob].faces_n;++j)
                norm_view.getv(i,j,cpu_mesh.cv[i_glob].norms[j]);
        }
        norm_view.release(false);
    }
    template<class MAP_NODES,class CPU_MESH>
    void    dump_nodes_geom_data(const MAP_NODES &map_n, CPU_MESH &cpu_mesh)const
    {
        t_tensor1_field_view_tml<T,dim,strg>                    node_coords_view(node_coords, true);
        for(int i_ = 0;i_ < map_n.get_size();++i_) {
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
