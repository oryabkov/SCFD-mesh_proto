// Copyright © 2016 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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

#ifndef __SCFD_NODES_PARTITIONER_H__
#define __SCFD_NODES_PARTITIONER_H__

#include <vector>
#include <map>
#include <algorithm>

namespace scfd
{
namespace communication
{

//WARNING i think we will reqiure MPI here!

//supposed to satisfy PARTITIONER concept

struct nodes_partitioner
{
    int                     total_size;
    int                     my_rank;
    bool                    is_complete;
    //first int is global index of node, second int is rank
    std::map<int,int>       ranks;
    //supposed to be sorted
    std::vector<int>        own_glob_indices;
    std::map<int,int>       own_glob_indices_2_ind;

    nodes_partitioner() {}
    //map_elem could be at 'before complete' stage
    //mesh must have all nodes incident(to owned by map_elen) and 2nd order incident (to owned by map_elem) elements nodes-elem graph
    template<class CPU_MESH, class MAP>
    nodes_partitioner(int comm_size, int _my_rank, const CPU_MESH &mesh, const MAP &map_elem) : my_rank(_my_rank), is_complete(false)
    {
        total_size = mesh.nodes.size();
        for (int i = 0;i < map_elem.get_size();++i) {
            int elem_glob_i = map_elem.own_glob_ind(i);
            for (int vert_i = 0;vert_i < mesh.cv_2_node_ids[elem_glob_i].vert_n;++vert_i) {
                int node_i = mesh.cv_2_node_ids[elem_glob_i].ids[vert_i];
                if (ranks.find(node_i) != ranks.end()) continue;
                std::pair<int,int> ref_range = mesh.node_2_cv_ids_ref[node_i];
                int min_elem_id = -1;
                for (int j = ref_range.first;j < ref_range.second;++j) {
                    int nb_elem_glob_i = mesh.node_2_cv_ids_data[j];
                    if ((min_elem_id == -1)||(nb_elem_glob_i < min_elem_id)) min_elem_id = nb_elem_glob_i;
                }
                if (map_elem.check_glob_owned(min_elem_id)) {
                    ranks[node_i] = my_rank;
                    own_glob_indices.push_back(node_i);
                }
            }
        }
        std::sort(own_glob_indices.begin(), own_glob_indices.end());
        for (int i = 0;i < own_glob_indices.size();++i) {
            own_glob_indices_2_ind[ own_glob_indices[i] ] = i;
        }
        //my_rank = part.my_rank;
        //is_complete = part.is_complete;
    }

    //'read' before construction part:
    int     get_own_rank()const { return my_rank; }
    //i_glob here could be any global index both before and after construction part
    bool    check_glob_owned(int i_glob)const 
    {
        std::map<int,int>::const_iterator       it = ranks.find(i_glob);
        if (it == ranks.end()) return false;
        return it->second == get_own_rank();
    }
    int     get_total_size()const { return total_size; }
    //get_size is not part of PARTITIONER concept
    int     get_size()const { return own_glob_indices.size(); }
    int     own_glob_ind(int i)const
    {
        assert((i >= 0)&&(i < get_size()));
        return own_glob_indices[i];
    }
    int     own_glob_ind_2_ind(int i_glob)const
    {
        std::map<int,int>::const_iterator       it = own_glob_indices_2_ind.find(i_glob);
        assert(it != own_glob_indices_2_ind.end());
        return it->second;
    }

    //'construction' part:
    void    add_stencil_element(int i_glob)
    {
        if (check_glob_owned(i_glob)) return;
        ranks[i_glob] = -1;
    }
    void    complete()
    {
        is_complete = true;
        //TODO i think the only way is to use mpi communications here
    }

    //'read' after construction part:
    //i_glob is ether index owner by calling process, ether index from stencil, otherwise behavoiur is  undefined
    //returns rank of process that owns i_glob index
    int     get_rank(int i_glob)const
    {
        assert(is_complete);
        std::map<int,int>::const_iterator it = ranks.find(i_glob);
        assert(it != ranks.end());
        if (it->second == -1) throw std::logic_error("nodes_partitioner:: not realized yet!!");
        return it->second;
    }
    //for (rank == get_own_rank()) result and behavoiur coincides with check_glob_owned(i_glob)
    //for rank != get_own_rank():
    //i_glob is ether index owner by calling process, ether index from stencil, otherwise behavoiur is  undefined
    //returns, if index i_glob is owned by rank processor
    bool    check_glob_owned(int i_glob, int rank)const 
    { 
        assert(is_complete);
        if (rank == get_own_rank()) {
            return check_glob_owned(i_glob);
        } else {
            return get_rank(i_glob) == rank; 
        }
    }
};

}  /// namespace communication
}  /// namespace scfd

#endif
