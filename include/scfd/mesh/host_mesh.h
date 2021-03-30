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

    /// Suppose we need to duplicate BasicMesh interface here?
    /// In case of inheritance we get it at once

    /// Nodes to elements graph access interface
    ordinal_type get_node_incident_elems_num(ordinal_type i)const;
    ordinal_type get_nodes_max_incident_elems_num()const;
    //TODO here i'm not sure about external storage for result; mb, return internal array point
    void         get_node_incident_elems(ordinal_type i,ordinal_type *elems)const;

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
    //std::shared_ptr<const BasicMesh>  basic_mesh_;

private:


};

}  /// namespace mesh
}  /// namespace scfd

#endif