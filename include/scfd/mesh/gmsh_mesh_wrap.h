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
#include <stdexcept>
#include <map>
#include "Gmsh.h"
#include "GModel.h"
#include "MElement.h"
#include "MQuadrangle.h"
#include "MTriangle.h"
#include "MPrism.h"
#include "MTetrahedron.h"

#ifndef __SCFD_MESH_GMSH_MESH_WRAP_H__
#define __SCFD_MESH_GMSH_MESH_WRAP_H__

namespace scfd
{
namespace mesh
{

class file_read_error : public std::runtime_error
{
public:
    file_read_error(const std::string &pref, const std::string &fn) : 
        std::runtime_error(pref + "Error while reading file " + fn)
    {
    }
};

/// Note that gmsh internally supports only doubles, so this interface also performs types conversions.
/// Besides that elements shift is performed (1st 3d element in gmsh is not likely to have index '0')
template<class T,int dim = 3,class Ord = int>
class gmsh_mesh_wrap
{
public:
    using scalar_type = T;
    using ordinal_type = Ord;
    using elem_type_ordinal_type = int;

public:
    /// No 'empty' state
    gmsh_mesh_wrap(const std::string &fn)
    {
        /// Initialization

        /// Read mesh
        g_model_ = std::make_shared<GModel>();
        if (!g_model_->readMSH(fn)) 
        {
            throw file_read_error("gmsh_mesh_wrap::gmsh_mesh_wrap()",fn);
        }

        /// Calculate elements_index_shift

        /// Init elements_tags

    }

    elem_type_ordinal_type get_elem_type(Ord i)const
    {

    }
    Ord get_elem_tag(Ord i)const
    {

    }
    Ord get_elems_max_prim_nodes_num()const
    {
                
    }
    Ord get_elems_max_nodes_num()const
    {

    }
    Ord get_elem_prim_nodes_num(Ord i)const
    {

    }
    void get_elem_prim_nodes(Ord i, Ord *nodes)const
    {

    }
    Ord get_elem_nodes_num(Ord i)const
    {

    }
    void get_elem_nodes(Ord i, Ord *nodes)const
    {

    }
    void get_node_coords(Ord i,T *coords)const
    {
        
    }

private:
    using elem_type_ord_t = elem_type_ordinal_type;

private:
    std::shared_ptr<GModel>         g_model_;
    Ord                             elements_index_shift;
    std::map<Ord,Ord>               elements_tags;
};

}  /// namespace mesh
}  /// namespace scfd

#endif