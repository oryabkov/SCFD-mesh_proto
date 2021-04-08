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

#include <memory>
#include <algorithm>

#ifndef __SCFD_MESH_HOST_MESH_H__
#define __SCFD_MESH_HOST_MESH_H__

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
    void enlarge_stencil(Ord ghost_level)
    {
        parent_type::enlarge_stencil(part, ghost_level);

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
        Ord j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            faces[j] = *it;
        }
    }
    /// Actually it can return either 1 or 2
    void         get_face_elems_num(ordinal_type i)const
    {
        //TODO some check in debug mode for result (1 or 2)
        return faces_to_elems_graph_.get_range_size(i);
    }
    void         get_face_elems(ordinal_type i,ordinal_type elems[2])const
    {
        auto it_range = faces_to_elems_graph_.get_range(i);
        Ord j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            faces[j] = *it;
        }
    }

    /// Neighbours0 interface
    /// No need to return number of neighbours0 - use get_elem_faces_num
    ordinal_type get_elem_neighbours0(ordinal_type i, ordinal_type *faces)const
    {

    }

private:
    struct face_key_t
    {
        //TODO temporal solution (max 4 nodes) but will be enough for most cases
        Ord     nodes_n_;
        Ord     sorted_prim_nodes_[4];

        face_key_t() = default;
        face_key_t(Ord nodes_n, Ord prim_nodes[4])
        {
            nodes_n_ = nodes_n;
            for (Ord j = 0;j < nodes_n_;++j)
                sorted_prim_nodes_[j] = prim_nodes[j];
            std::sort(sorted_prim_nodes_,sorted_prim_nodes_+nodes_n_);
        }

        Ord     nodes_n()const
        {
            return nodes_n_;
        }
        Ord     sorted_prim_node(Ord j)
        {
            return sorted_prim_nodes_[j];
        }
    };
    struct face_key_equal_func
    {
        bool operator()(const face_key_t &f1, const face_key_t &f2)const
        {
            if (f1.nodes_n() != f2.nodes_n()) 
                return false;
            for (int i = 0; i < f1.nodes_n(); i++) 
            {
                if (f1.sorted_prim_node(i) != f2.sorted_prim_node(i)) 
                    return false;
            }
            return true;
        }
    };
    struct face_key_less_func
    {
        bool operator()(const face_key_t &f1, const face_key_t &f2) const
        {
            if (f1.nodes_n() < f2.nodes_n()) return true;
            if (f1.nodes_n() > f2.nodes_n()) return false;
            for (int i = 0; i < f1.nodes_n(); i++) 
            {
                if (f1.sorted_prim_node(i) < f2.sorted_prim_node(i)) return true;
                if (f1.sorted_prim_node(i) > f2.sorted_prim_node(i)) return false;
            }
            return false;
        }
    }; 

    using elems_to_faces_graph_t = detail::ranges_sparse_arr<Ord,Ord>;
    /// Here pair's first is id of incident element, second - local face index inside this element
    using faces_to_elems_graph_t = detail::ranges_sparse_arr<std::pair<Ord,Ord>,Ord>;
    /// Here pair's first is id of incident element, second - local face index inside this element
    using elems_to_neighbours0_graph_t = detail::ranges_sparse_arr<std::pair<Ord,Ord>,Ord>;

    //std::shared_ptr<const BasicMesh>  basic_mesh_;

private:
    elems_to_faces_graph_t          elems_to_faces_graph_;
    faces_to_elems_graph_t          faces_to_elems_graph_;
    elems_to_neighbours0_graph_t    elems_to_neighbours0_graph_;

protected:

    //TODO virtual?
    void build_faces(Ord ghost_level)
    {
        std::set<Ord> stencil_ids;
        calc_elems_stencil(ghost_level,stencil_ids);
        const auto &part = *parent_type::get_partitioner();

        /// Create faces dict

        std::map<face_key_t,Ord>    faces;
        /// Process own elements
        for (Ord i = 0;i < part.get_size();++i)
        {
            Ord  elem_id = part.own_glob_ind(i);
            build_faces_for_elem(elem_id,faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            build_faces_for_elem(elem_id,faces);
        }

        /// Create graphs elems_to_faces_graph_ and faces_to_elems_graph_ based on dict
    
        elems_to_faces_graph_.reserve(faces.size());
        faces_to_elems_graph_.reserve(part.get_size()+stencil_ids.size());

        /// Estmate sizes of graphs elems_to_faces_graph_ and faces_to_elems_graph_
        /// Process own elements
        for (Ord i = 0;i < part.get_size();++i)
        {
            Ord  elem_id = part.own_glob_ind(i);
            reserve_graphs_for_elem(elem_id, faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            reserve_graphs_for_elem(elem_id, faces);
        }

        /// Complete graphs structures
        elems_to_faces_graph_.complete_structure();
        faces_to_elems_graph_.complete_structure();

        /// Fill actual graphs elems_to_faces_graph_ and faces_to_elems_graph_
        /// Process own elements
        for (Ord i = 0;i < part.get_size();++i)
        {
            Ord  elem_id = part.own_glob_ind(i);
            fill_graphs_for_elem(elem_id, faces);
        }
        /// Process stencil elements
        for (auto elem_id : stencil_ids)
        {
            fill_graphs_for_elem(elem_id, faces);
        }

    }
    void build_faces_for_elem(Ord elem_id, std::map<face_key_t,Ord> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        //TODO in general it is not good idea to take it each time (some cached value?)
        Ord glob_max_faces_num = parent_type::get_elems_glob_max_faces_num();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        //ISSUE are there any performance problmes with local allocation here?
        Ord nodes[parent_type::get_elems_max_prim_nodes_num()];
        parent_type::get_elem_prim_nodes(elem_id, nullptr, nodes);
        for (Ord j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            Ord face_id = elem_id*glob_max_faces_num + j;
            //TODO temporal solution (max 4 nodes) but will be enough for most cases
            Ord face_nodes[4];
            for (Ord face_vert_i = 0;face_vert_i < ref.get_face_verts_n(elem_type,j);++face_vert_i)
            {
                face_nodes[face_vert_i] = nodes[ref.get_face_vert_i(elem_type,j,face_vert_i)];
            }
            face_key_t face_key(ref.get_face_verts_n(elem_type,j), face_nodes);
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
    void reserve_graphs_for_elem(Ord elem_id, const std::map<face_key_t,Ord> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        Ord nodes[parent_type::get_elems_max_prim_nodes_num()];
        parent_type::get_elem_prim_nodes(elem_id, nullptr, nodes);
        for (Ord j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            Ord face_nodes[4];
            for (Ord face_vert_i = 0;face_vert_i < ref.get_face_verts_n(elem_type,j);++face_vert_i)
            {
                face_nodes[face_vert_i] = nodes[ref.get_face_vert_i(elem_type,j,face_vert_i)];
            }
            face_key_t face_key(ref.get_face_verts_n(elem_type,j), face_nodes);
            auto face_it = faces.find(face_key);
            if (face_it == faces.end())
                throw std::logic_error("host_mesh::reserve_graphs_for_elem: no face found!");
            Ord face_id = face_it->second;
            elems_to_faces_graph_.inc_max_range_size(elem_id,1);
            faces_to_elems_graph_.inc_max_range_size(face_id,1);
        }
    }
    void fill_graphs_for_elem(Ord elem_id, const std::map<face_key_t,Ord> &faces)
    {
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        Ord nodes[parent_type::get_elems_max_prim_nodes_num()];
        parent_type::get_elem_prim_nodes(elem_id, nullptr, nodes);
        for (Ord j = 0;j < ref.get_faces_n(elem_type);++j)
        {
            Ord face_nodes[4];
            for (Ord face_vert_i = 0;face_vert_i < ref.get_face_verts_n(elem_type,j);++face_vert_i)
            {
                face_nodes[face_vert_i] = nodes[ref.get_face_vert_i(elem_type,j,face_vert_i)];
            }
            face_key_t face_key(ref.get_face_verts_n(elem_type,j), face_nodes);
            auto face_it = faces.find(face_key);
            if (face_it == faces.end())
                throw std::logic_error("host_mesh::reserve_graphs_for_elem: no face found!");
            Ord face_id = face_it->second;
            elems_to_faces_graph_.add_to_range(elem_id, face_id);            
            faces_to_elems_graph_.add_to_range(face_id,std::pair<Ord,Ord>(elem_id,j));
        }
    }
    void calc_elems_stencil(Ord ghost_level,std::set<Ord> &stencil_ids)
    {
        const auto &part = *parent_type::get_partitioner();
        stencil_ids.clear();
        std::set<Ord> curr_ids,next_ids;
        for (Ord i = 0;i < part.get_size();++i)
        {
            Ord  elem_id = part.own_glob_ind(i);
            calc_elems_stencil_for_elem(elem_id,stencil_ids,next_ids);
        }
        for (Ord level = 1;level < ghost_level;++level)
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
        Ord elem_id, std::set<Ord> &stencil_ids, std::set<Ord> &next_ids
    )
    {
        const auto &part = *parent_type::get_partitioner();
        const auto &ref = parent_type::mesh_elem_reference();
        elem_type_ordinal_type  elem_type = parent_type::get_elem_type(elem_id);
        Ord nodes_n;
        Ord nodes[parent_type::get_elems_max_prim_nodes_num()];
        parent_type::get_elem_prim_nodes(elem_id, &nodes_n, nodes);
        for (Ord node_i = 0;node_i < nodes_n;++node_i)
        {
            Ord incid_elems_n;
            Ord incid_elems[parent_type::get_nodes_max_incident_elems_num()];            
            parent_type::get_node_incident_elems(nodes[node_i],incid_elems,&incid_elems_n);
            for (Ord incid_elems_i = 0;incid_elems_i < incid_elems_n;++incid_elems_i)
            {
                Ord incid_elems_id = incid_elems[incid_elems_i];
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