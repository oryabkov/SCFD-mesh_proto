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

#ifndef __SCFD_MESH_GMSH_MESH_WRAP_H__
#define __SCFD_MESH_GMSH_MESH_WRAP_H__

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
#include "detail/ranges_sparse_arr.h"
#include "gmsh_mesh_elem_reference.h"

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
    using mesh_elem_reference_type = gmsh_mesh_elem_reference<T>;

public:
    /// No 'empty' state
    gmsh_mesh_wrap()
    {
        /// Initialization
        g_model_ = std::make_shared<GModel>();
    }

    const mesh_elem_reference_type &mesh_elem_reference()const
    {
        return mesh_elem_reference_;
    }

    /// Following reads are gmsh specific (not part of BasicMesh concept)
    void set_mesh_filename(const std::string &fn)
    {
        fn_ = fn;
    }
    /// only one call to either read() or read_parted() method is allowed during lifetime
    void read()
    {
        /// Read mesh
        if (!g_model_->readMSH(fn_)) 
        {
            throw file_read_error("gmsh_mesh_wrap::read(): ",fn);
        }

        std::vector<GEntity*> entities;
        g_model_->getEntities(entities, dim);

        /// Calc elements_index_shift_, elems_max_prim_nodes_num_, elems_max_nodes_num_
        /// and check for elements tags contigiousness
        Ord elems_n = g_model_->getNumMeshElements(dim),
            min_elem_tag = std::numeric_limits<int>::max(),
            max_elem_tag = std::numeric_limits<int>::min();
        elems_max_prim_nodes_num_ = 0;
        elems_max_nodes_num_ = 0;
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                min_elem_tag = std::min(min_elem_tag,s->getNum());
                max_elem_tag = std::max(max_elem_tag,s->getNum());
                elems_max_prim_nodes_num_ = std::max(elems_max_prim_nodes_num_,s->getNumPrimaryVertices());
                elems_max_nodes_num_ = std::max(elems_max_nodes_num_,s->getNumVertices());
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
        //TODO would be usefull to call it but don't know how to efficently take nodes number from gmsh
        //nodes_to_elems_graph_.reserve(nodes_num?)

        /// First build group ids map and estimate nodes to elements graph sizes
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                Ord elem_id = elem_tag_to_elem_id(s->getNum());
                elements_group_ids_[elem_id] = e->tag();
                for (Ord elem_vert_i = 0; elem_vert_i < s->getNumVertices(); ++elem_vert_i) 
                {
                    nodes_to_elems_graph_.inc_max_range_size(s->getVertex(elem_vert_i)->getNum(),1);
                }
            }
        }
        /// Complete graph structure
        nodes_to_elems_graph_.complete_structure();
        /// Fill actual incidence data
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                Ord elem_id = elem_tag_to_elem_id(s->getNum());
                for (Ord elem_vert_i = 0; elem_vert_i < s->getNumVertices(); ++elem_vert_i)
                {
                    //nodes_to_elems_graph_.inc_max_range_size(s->getVertex(elem_vert_i)->getNum(),1);
                    nodes_to_elems_graph_.add_to_range
                    (
                        s->getVertex(elem_vert_i)->getNum(),
                        std::pair<Ord,Ord>(elem_id,elem_vert_i)
                    );
                }
            }
        }
        
    }
    void read_parted()
    {
        //TODO
    }

    /// PartElems satisfies Partitioner concept without own_glob_ind_2_ind, so Map can be used here
    /// This read method is part of the BasicMesh concept
    /// This read method can only be used after read() or read_parted() calls, its purpose is
    /// to ensure that all elements defined by part together with ghost_level number of 2nd (nodal) 
    /// neighbours are read and accesible through interface calls.
    template<class PartElems>
    void read(const PartElems &part, Ord ghost_level = 1)
    {
        //For now do nothing, because we read whole mesh data in read() method
    }

    /// Elements interface

    elem_type_ordinal_type get_elem_type(Ord i)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        return s->getType();
    }
    Ord get_elem_group_id(Ord i)const
    {
        return elements_group_ids_[i];
    }
    Ord get_elems_max_prim_nodes_num()const
    {
        return elems_max_prim_nodes_num_;
    }
    Ord get_elem_prim_nodes_num(Ord i)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        return s->getNumPrimaryVertices();
    }
    /// Here theoretically types convesion could be done, so extrnal space is used
    void get_elem_prim_nodes(Ord i, Ord &prim_nodes_num, Ord *nodes)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        prim_nodes_num = s->getNumPrimaryVertices();
        /// TODO we explicitly use here that primary vertices goes 1st (used in mesh_prepare)
        /// Chrch is this always true.
        for (Ord j = 0;j < s->getNumPrimaryVertices();++j)
        {   
            nodes[j] = s->getVertex(j)->getNum();
        }
    }
    Ord get_elems_max_nodes_num()const
    {
        return elems_max_nodes_num_;
    }
    Ord get_elem_nodes_num(Ord i)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        return s->getNumVertices();
    }
    void get_elem_nodes(Ord i, Ord &nodes_num, Ord *nodes)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        nodes_num = s->getNumVertices();
        /// TODO we explicitly use here that primary vertices goes 1st (used in mesh_prepare)
        /// Chrch is this always true.
        for (Ord j = 0;j < s->getNumVertices();++j)
        {   
            nodes[j] = s->getVertex(j)->getNum();
        }
    }

    /// Nodes interface

    void get_node_coords(Ord i,T *coords)const
    {
        GVertex *v = g_model_->getVertexByTag(node_id_to_node_tag(i));
        if (dim >= 1) coords[0] = static_cast<T>(v->x());
        if (dim >= 2) coords[1] = static_cast<T>(v->y());
        if (dim >= 3) coords[2] = static_cast<T>(v->z());
    }
    /// Nodes to elements graph access interface
    Ord get_node_incident_elems_num(Ord i)const
    {
        return nodes_to_elems_graph_.get_range_size(i);
    }
    Ord get_nodes_max_incident_elems_num()const
    {
        return nodes_to_elems_graph_.get_max_ranges_size();
    }
    //TODO here i'm not sure about external storage for result; mb, return internal array point
    void get_node_incident_elems(Ord i,Ord *elems)const
    {
        auto it_range = get_range(i);
        Ord j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            elems[j] = it->first;
        }
    }

private:
    using elem_type_ord_t = elem_type_ordinal_type;
    /// Here pair's first is id of incident element, second - local node index inside this element
    using nodes_to_elems_graph_t = detail::ranges_sparse_arr<std::pair<Ord,Ord>,Ord>;

private:
    mesh_elem_reference_type        mesh_elem_reference_;

    std::string                     fn_;
    /// Stick to old private C++ API
    std::shared_ptr<GModel>         g_model_;
    Ord                             elements_index_shift_;
    Ord                             elems_max_prim_nodes_num_,
                                    elems_max_nodes_num_;
    nodes_to_elems_graph_t          nodes_to_elems_graph_;
    std::map<Ord,Ord>               elements_group_ids_;

    /// Converts internal gmsh tag into 'visible' element index
    Ord elem_tag_to_elem_id(Ord elem_tag)const
    {
        return elem_tag - elements_index_shift_;
    }
    Ord elem_id_to_elem_tag(Ord elem_id)const
    {
        return elem_id + elements_index_shift_;
    }
    //For now seems we can use non-contigious node ids
    Ord node_tag_to_node_id(Ord node_tag)const
    {
        return node_tag;
    }
    Ord node_id_to_node_tag(Ord node_id)const
    {
        return node_id;
    }
};

}  /// namespace mesh
}  /// namespace scfd

#endif