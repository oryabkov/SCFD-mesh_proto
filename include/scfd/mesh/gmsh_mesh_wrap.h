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
#include <limits>
#include <gmsh/GmshGlobal.h>
#include <gmsh/GModel.h>
#include <gmsh/MElement.h>
#include <gmsh/MQuadrangle.h>
#include <gmsh/MTriangle.h>
#include <gmsh/MPrism.h>
#include <gmsh/MTetrahedron.h>

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
    gmsh_mesh_wrap()
    {
        /// Initialization
        g_model_ = std::make_shared<GModel>();
    }

    void set_mesh_filename(const std::string &fn)
    {
        fn_ = fn;
    }
    void read()
    {
        /// Read mesh
        if (!g_model_->readMSH(fn_)) 
        {
            throw file_read_error("gmsh_mesh_wrap::read(): ",fn);
        }

        std::vector<GEntity*> entities;
        g_model_->getEntities(entities, dim);

        /// Build elements_index_shift_ and check tags contigiousness
        Ord elems_n = g_model_->getNumMeshElements(dim),
            min_elem_tag = std::numeric_limits<int>::max(),
            max_elem_tag = std::numeric_limits<int>::min();
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                min_elem_tag = std::min(min_elem_tag,s->getNum());
                max_elem_tag = std::max(max_elem_tag,s->getNum());
            }
        }
        if (max_elem_tag+1-min_elem_tag != elems_n)
            throw 
                std::logic_error
                (
                    "gmsh_mesh_wrap::read(): case of non-contigious tags is not supported yet"
                );
        elements_index_shift_ = min_elem_tag;

        /// Build entities tags map and nodes to elements graph
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                Ord elem_i = elem_tag_to_elem_i(s->getNum());
                elements_group_ids_[elem_i] = e->tag();
            }
        }
    }
    template<class MapElems>
    void read(const MapElems &map, Ord ghost_level = 1)
    {
        
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
    /// Here theoretically types convesion could be done, so extrnal space is used
    
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
    std::string                     fn_;
    /// Stick to old private C++ API
    std::shared_ptr<GModel>         g_model_;
    Ord                             elements_index_shift_;
    std::map<Ord,Ord>               elements_group_ids_;

    /// Converts internal gmsh tag into 'visible' element index
    Ord elem_tag_to_elem_i(Ord elem_tag)const
    {
        return elem_tag - elements_index_shift_;
    }
    Ord elem_i_to_elem_tag(Ord elem_i)const
    {
        return elem_i + elements_index_shift_;
    }
};

}  /// namespace mesh
}  /// namespace scfd

#endif