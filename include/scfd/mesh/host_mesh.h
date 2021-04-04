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
    template<class PartElems>
    void read(const PartElems &part, Ord ghost_level = 1)
    {
        parent_type::read(part, ghost_level);

        /// Build faces
        std::map<face_key_t,Ord>    faces;
    }

    /// Suppose we need to duplicate BasicMesh interface here?
    /// In case of inheritance we get it at once

    /// Faces interface
    ordinal_type get_elem_faces_num(ordinal_type i)const;
    void         get_elem_faces(ordinal_type i, ordinal_type *faces)const;
    /// Actually it can return either 1 or 2
    void         get_face_elems_num(ordinal_type i)const;
    void         get_face_elems(ordinal_type i,ordinal_type elems[2])const;

    /// Neighbours0 interface
    /// No need to return number of neighbours0 - use get_elem_faces_num
    ordinal_type get_elem_neighbours0(ordinal_type i, ordinal_type *faces)const;

private:
    struct face_key_t
    {
        //TODO temporal solution (max 4 nodes) but will be enough for most cases
        Ord     nodes_n_;
        Ord     sorted_prim_nodes_[4];

        face_key() = default;
        face_key(Ord nodes_n, Ord prim_nodes[4])
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


    //std::shared_ptr<const BasicMesh>  basic_mesh_;


private:


};

}  /// namespace mesh
}  /// namespace scfd

#endif