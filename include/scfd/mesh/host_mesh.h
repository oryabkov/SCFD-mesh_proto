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

#ifndef __SCFD_MESH_HOST_MESH_H__
#define __SCFD_MESH_HOST_MESH_H__

#include <memory>
#include <algorithm>
#include "detail/face_key.h"

//TODO building face_key for element face fully repeated 3 times - move into separate method
//TODO walkthrough part + stencil repeats several times - mb, create some method with functor parameter?
//TODO graphs (elems_to_faces_graph_ and faces_to_elems_graph_) are build for more elements then intended 
//(last ghost level) including incorrect ones (outer faces). Not a problem if no one will use them but some 
//inconsistency presents.

namespace scfd
{
namespace mesh
{

//ISSUE what to use: inheritance or aglommeration?
template<class BasicMesh>
class host_mesh : public BasicMesh
{
public:
    using parent_type = BasicMesh;
    using basic_mesh_type = BasicMesh;
    using scalar_type = typename basic_mesh_type::scalar_type;
    using ordinal_type = typename basic_mesh_type::ordinal_type;
    using elem_type_ordinal_type = typename basic_mesh_type::elem_type_ordinal_type;
    using partitioner_type = typename basic_mesh_type::partitioner_type;

public:
    /// Need default contructable BasicMesh
    host_mesh()
    {
        /// Initialization
    }
    /// No 'empty' state
    /*host_mesh(std::shared_ptr<const BasicMesh>  basic_mesh) : 
        basic_mesh_(basic_mesh)
    {
        /// Initialization
    }*/

    /// See gmsh_mesh_wrap.h for PartElems description
    void enlarge_stencil(ordinal_type ghost_level)
    {
        parent_type::enlarge_stencil(ghost_level);
        build_faces(ghost_level);
    }

    /// Suppose we need to duplicate BasicMesh interface here?
    /// In case of inheritance we get it at once

    /// Faces interface
    ordinal_type get_elem_faces_num(ordinal_type i)const
    {
        return elems_to_faces_graph_.get_range_size(i);
    }
    void         get_elem_faces(ordinal_type i, ordinal_type *faces)const
    {
        auto it_range = elems_to_faces_graph_.get_range(i);
        ordinal_type j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            faces[j] = *it;
        }
    }
    /// Actually it can return either 1 or 2
    ordinal_type get_face_elems_num(ordinal_type i)const
    {
        //TODO some check in debug mode for result (1 or 2)
        return faces_to_elems_graph_.get_range_size(i);
    }
    void         get_face_elems(ordinal_type i,ordinal_type elems[2])const
    {
        auto it_range = faces_to_elems_graph_.get_range(i);
        ordinal_type j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            elems[j] = it->first;
        }
    }

    /// Neighbours0 interface
    /// No need to return number of neighbours0 - use get_elem_faces_num
    void get_elem_neighbours0(ordinal_type i, ordinal_type *elems)const
    {
        auto it_range = elems_to_neighbours0_graph_.get_range(i);
        ordinal_type j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            elems[j] = it->first;
        }
    }

private:
    using face_key_t = detail::face_key<ordinal_type>;
    using face_key_equal_func = detail::face_key_equal_func<ordinal_type>;
    using face_key_less_func = detail::face_key_less_func<ordinal_type>;

    using elems_to_faces_graph_t = detail::ranges_sparse_arr<ordinal_type,ordinal_type>;
    /// Here pair's first is id of incident element, second - local face index inside this element
    using faces_to_elems_graph_t = detail::ranges_sparse_arr<std::pair<ordinal_type,ordinal_type>,ordinal_type>;
    /// Here pair's first is id of incident element, second - local face index inside this element
    using elems_to_neighbours0_graph_t = detail::ranges_sparse_arr<std::pair<ordinal_type,ordinal_type>,ordinal_type>;

    //std::shared_ptr<const BasicMesh>  basic_mesh_;

private:
    elems_to_faces_graph_t          elems_to_faces_graph_;
    faces_to_elems_graph_t          faces_to_elems_graph_;
    elems_to_neighbours0_graph_t    elems_to_neighbours0_graph_;

protected:

    //TODO virtual?
    void build_faces(ordinal_type ghost_level)
    {
        std::set<ordinal_type> stencil_ids;
        calc_elems_stencil(ghost_level,stencil_ids);
        const auto &part = *parent_type::get_partitioner();

        /// Create faces dict

        std::map<face_key_t,ordinal_type,face_key_less_func>    faces;
        /// Process own elements
        for (ordinal_type i = 0;i < part.get_size();++i)
        {
            ordinal_type  elem_id = part.own_glob_ind(i);
            build_faces_for_elem(elem_id,faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            build_faces_for_elem(elem_id,faces);
        }

        /// Create graphs elems_to_faces_graph_ and faces_to_elems_graph_ based on dict
    
        elems_to_faces_graph_.reserve(part.get_size()+stencil_ids.size());
        elems_to_neighbours0_graph_.reserve(part.get_size()+stencil_ids.size());
        faces_to_elems_graph_.reserve(faces.size());

        /// Estmate sizes of graphs elems_to_faces_graph_ and faces_to_elems_graph_
        /// Process own elements
        for (ordinal_type i = 0;i < part.get_size();++i)
        {
            ordinal_type  elem_id = part.own_glob_ind(i);
            reserve_graphs_for_elem(elem_id, faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            reserve_graphs_for_elem(elem_id, faces);
        }

        /// Complete graphs structures
        elems_to_faces_graph_.complete_structure();
        elems_to_neighbours0_graph_.complete_structure();
        faces_to_elems_graph_.complete_structure();

        /// Fill actual graphs elems_to_faces_graph_ and faces_to_elems_graph_
        /// Process own elements
        for (ordinal_type i = 0;i < part.get_size();++i)
        {
            ordinal_type  elem_id = part.own_glob_ind(i);
            fill_graphs_for_elem(elem_id, faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            fill_graphs_for_elem(elem_id, faces);
        }

        /// Fill actual graph elems_to_neighbours0_graph_
        /// Process own elements
        for (ordinal_type i = 0;i < part.get_size();++i)
        {
            ordinal_type  elem_id = part.own_glob_ind(i);
            fill_neib_graph_for_elem(elem_id);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            fill_neib_graph_for_elem(elem_id);
        }
    }
    void build_faces_for_elem(ordinal_type elem_id, std::map<face_key_t,ordinal_type,face_key_less_func> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        //TODO in general it is not good idea to take it each time (some cached value?)
        ordinal_type glob_max_faces_num = parent_type::get_elems_glob_max_faces_num();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        //ISSUE are there any performance problmes with local allocation here?
        for (ordinal_type j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            ordinal_type face_id = elem_id*glob_max_faces_num + j;
            face_key_t face_key(*this,elem_id,j);
            auto face_it = faces.find(face_key);
            if (face_it == faces.end())
            {
                faces[face_key] = face_id;
            }
            else
            {
                /// Between two poosible ids for face we choose the minimal one
                /// (i.e. one that was taken from element with minimal id)
                face_it->second = std::min(face_it->second,face_id);
            }
        }
    }
    void reserve_graphs_for_elem(ordinal_type elem_id, const std::map<face_key_t,ordinal_type,face_key_less_func> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        for (ordinal_type j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            face_key_t face_key(*this,elem_id,j);
            auto face_it = faces.find(face_key);
            if (face_it == faces.end())
                throw std::logic_error("host_mesh::reserve_graphs_for_elem: no face found!");
            ordinal_type face_id = face_it->second;
            elems_to_faces_graph_.inc_max_range_size(elem_id,1);
            elems_to_neighbours0_graph_.inc_max_range_size(elem_id,1);
            faces_to_elems_graph_.inc_max_range_size(face_id,1);
        }
    }
    void fill_graphs_for_elem(ordinal_type elem_id, const std::map<face_key_t,ordinal_type,face_key_less_func> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        for (ordinal_type j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            face_key_t face_key(*this,elem_id,j);
            auto face_it = faces.find(face_key);
            if (face_it == faces.end())
                throw std::logic_error("host_mesh::reserve_graphs_for_elem: no face found!");
            ordinal_type face_id = face_it->second;
            elems_to_faces_graph_.add_to_range(elem_id, face_id);            
            faces_to_elems_graph_.add_to_range(face_id,std::pair<ordinal_type,ordinal_type>(elem_id,j));
        }
    }
    void fill_neib_graph_for_elem(ordinal_type elem_id)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        auto it_range = elems_to_faces_graph_.get_range(elem_id);
        ordinal_type loc_face_i = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++loc_face_i)
        {
            ordinal_type face_id = *it;
            if (get_face_elems_num(face_id) != 2) 
            {
                //TODO -1 -> to some special value
                elems_to_neighbours0_graph_.add_to_range
                (
                    elem_id, std::pair<ordinal_type,ordinal_type>(-1,-1)
                );
                continue;
            }
            ordinal_type elems[2];
            ordinal_type neib0_id;
            get_face_elems(face_id,elems);
            /// We now know there are two elements for this face
            for (ordinal_type j = 0;j < 2;++j)
            {
                if (elems[j] == elem_id) continue;
                neib0_id = elems[j];
            }
            elems_to_neighbours0_graph_.add_to_range
            (
                elem_id, std::pair<ordinal_type,ordinal_type>(neib0_id,loc_face_i)
            );
        }
    }
    void calc_elems_stencil(ordinal_type ghost_level,std::set<ordinal_type> &stencil_ids)
    {
        const auto &part = *parent_type::get_partitioner();
        stencil_ids.clear();
        std::set<ordinal_type> curr_ids,next_ids;
        for (ordinal_type i = 0;i < part.get_size();++i)
        {
            ordinal_type  elem_id = part.own_glob_ind(i);
            calc_elems_stencil_for_elem(elem_id,stencil_ids,next_ids);
        }
        for (ordinal_type level = 1;level < ghost_level;++level)
        {
            std::swap(curr_ids,next_ids);
            next_ids.clear();
            for (auto elem_id : curr_ids)
            {
                calc_elems_stencil_for_elem(elem_id,stencil_ids,next_ids);
            }
        }
    }
    void calc_elems_stencil_for_elem
    (
        ordinal_type elem_id, std::set<ordinal_type> &stencil_ids, std::set<ordinal_type> &next_ids
    )
    {
        const auto &part = *parent_type::get_partitioner();
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        ordinal_type nodes_n;
        ordinal_type nodes[parent_type::get_elems_max_prim_nodes_num()];
        parent_type::get_elem_prim_nodes(elem_id, &nodes_n, nodes);
        for (ordinal_type node_i = 0;node_i < nodes_n;++node_i)
        {
            ordinal_type incid_elems_n;
            ordinal_type incid_elems[parent_type::get_nodes_max_incident_elems_num()];            
            parent_type::get_node_incident_elems(nodes[node_i],incid_elems,&incid_elems_n);
            for (ordinal_type incid_elems_i = 0;incid_elems_i < incid_elems_n;++incid_elems_i)
            {
                ordinal_type incid_elems_id = incid_elems[incid_elems_i];
                /// Simple cutoff (also not very effective)
                if (incid_elems_id == elem_id) continue;
                if (part.check_glob_owned(incid_elems_id)) continue;
                if (stencil_ids.find(incid_elems_id) != stencil_ids.end()) continue;
                stencil_ids.insert(incid_elems_id);
                next_ids.insert(incid_elems_id);
            }
        }
    }

};

}  /// namespace mesh
}  /// namespace scfd

#endif