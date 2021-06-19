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
#include <scfd/static_vec/vec.h>
#include <scfd/static_mat/mat.h>
#include "detail/ranges_sparse_arr.h"
#include "detail/face_key.h"
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

/// PartElems satisfies Partitioner concept without own_glob_ind_2_ind, so Map can be used here
/// Note that gmsh internally supports only doubles, so this interface also performs types conversions.
/// Besides that elements shift is performed (1st 3d element in gmsh is not likely to have index '0')
template<class T,class PartElems,int Dim = 3,class Ord = int>
class gmsh_mesh_wrap
{
public:
    using scalar_type = T;
    using ordinal_type = Ord;
    using elem_type_ordinal_type = int;
    using partitioner_type = PartElems;
    using mesh_elem_reference_type = gmsh_mesh_elem_reference<T>;
    
    static const int          dim = Dim;
    static const ordinal_type special_id = std::numeric_limits<ordinal_type>::max();

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
            throw file_read_error("gmsh_mesh_wrap::read(): ",fn_);
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
        bool is_first = true;
        is_homogeneous_ = true;
        for (auto e : entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                if (is_first)
                {
                    homogeneous_elem_type_ = s->getTypeForMSH();
                }
                else
                {
                    if (homogeneous_elem_type_ != s->getTypeForMSH()) is_homogeneous_ = false;
                } 

                min_elem_tag = std::min(min_elem_tag,Ord(s->getNum()));
                max_elem_tag = std::max(max_elem_tag,Ord(s->getNum()));
                elems_max_prim_nodes_num_ = std::max(elems_max_prim_nodes_num_,Ord(s->getNumPrimaryVertices()));
                elems_max_nodes_num_ = std::max(elems_max_nodes_num_,Ord(s->getNumVertices()));
                is_first = false;
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
        /// Also calcs elems_max_faces_num_
        elems_max_faces_num_ = 0;
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
                elems_max_faces_num_ = std::max(elems_max_faces_num_,s->getNumFaces());
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

        std::vector<GEntity*> bnd_entities;
        g_model_->getEntities(bnd_entities, dim-1);
        for (auto e : bnd_entities)
        {
            for (Ord j = 0; j < e->getNumMeshElements(); ++j)
            {
                MElement *s = e->getMeshElement(j);
                //TODO temporal solution (max 4 nodes) but will be enough for most cases
                Ord nodes[4];
                if (s->getNumPrimaryVertices() > 4)
                    throw 
                        std::logic_error
                        (
                            "gmsh_mesh_wrap::read(): s->getNumPrimaryVertices() > 4 for bnd elements is not supported"
                        );
                for (Ord vert_i = 0;vert_i < s->getNumPrimaryVertices();++vert_i)
                {   
                    nodes[vert_i] = node_tag_to_node_id(s->getVertex(vert_i)->getNum());
                }
                face_key_t face_key(s->getNumPrimaryVertices(), nodes);
                bnd_faces_group_ids_[face_key] = e->tag();
            }
        }

        
    }
    //ISSUE is it part of BasicMesh concept?
    Ord get_total_elems_num()const
    {
        return g_model_->getNumMeshElements(dim);
    }
    //TODO think read_parted is not suited here; different wrap class is requered for parted gmsh read
    /*void read_parted()
    {
        //TODO
    }*/
    void set_partitioner(const std::shared_ptr<PartElems> &elems_partitioner)
    {
        elems_partitioner_ = elems_partitioner;
    }

    std::shared_ptr<PartElems> get_partitioner()const
    {
        return elems_partitioner_;
    }

    /// This read method is part of the BasicMesh concept
    /// This read method can only be used after read() or read_parted() calls, its purpose is
    /// to ensure that all elements defined by part together with ghost_level number of 2nd (nodal) 
    /// neighbours are read and accesible through interface calls.
    void enlarge_stencil(Ord ghost_level)
    {
        //For this wrap type do nothing, because we read whole mesh data in read() method
    }

    /// Elements interface

    elem_type_ordinal_type get_elem_type(Ord i)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        return s->getTypeForMSH();
    }
    bool                   is_homogeneous()const
    {
        return is_homogeneous_;
    }
    elem_type_ordinal_type homogeneous_elem_type()const
    {
        return homogeneous_elem_type_;
    }
    Ord get_elem_group_id(Ord i)const
    {
        return elements_group_ids_.at(i);
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
    void get_elem_prim_nodes(Ord i, Ord *nodes, Ord *prim_nodes_num = nullptr)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        if (prim_nodes_num)
            *prim_nodes_num = s->getNumPrimaryVertices();
        /// TODO we explicitly use here that primary vertices goes 1st (used in mesh_prepare)
        /// Chrch is this always true.
        for (Ord j = 0;j < s->getNumPrimaryVertices();++j)
        {   
            nodes[j] = node_tag_to_node_id(s->getVertex(j)->getNum());
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
    void get_elem_nodes(Ord i, Ord *nodes, Ord *nodes_num = nullptr)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        if (nodes_num) *nodes_num = s->getNumVertices();
        /// TODO we explicitly use here that primary vertices goes 1st (used in mesh_prepare)
        /// Chrch is this always true.
        for (Ord j = 0;j < s->getNumVertices();++j)
        {   
            nodes[j] = node_tag_to_node_id(s->getVertex(j)->getNum());
        }
    }
    Ord get_elem_node(Ord i, Ord j)const
    {
        MElement *s = g_model_->getMeshElementByTag(elem_id_to_elem_tag(i));
        return node_tag_to_node_id(s->getVertex(j)->getNum());
    }

    /// Nodes interface

    void get_node_coords(Ord i,T *coords)const
    {
        MVertex *v = g_model_->getMeshVertexByTag(node_id_to_node_tag(i));
        if (dim >= 1) coords[0] = static_cast<T>(v->x());
        if (dim >= 2) coords[1] = static_cast<T>(v->y());
        if (dim >= 3) coords[2] = static_cast<T>(v->z());
    }
    Ord get_node_group_id(Ord i)const
    {
        //TODO
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
    void get_node_incident_elems(Ord i,Ord *elems, Ord *elems_num = nullptr)const
    {
        if (elems_num)
            *elems_num = get_node_incident_elems_num(i);
        auto it_range = nodes_to_elems_graph_.get_range(i);
        Ord j = 0;
        for (auto it = it_range.first;it != it_range.second;++it,++j)
        {
            elems[j] = it->first;
        }
    }

    /// Parts of face interface on this level

    /// Maximum faces per element across all mesh (perhaps, only for local part)
    Ord get_elems_max_faces_num()const
    {
        return elems_max_faces_num_;
    }
    /// Maximum faces per element across all mesh (not only local part)
    Ord get_elems_glob_max_faces_num()const
    {
        return elems_max_faces_num_;
    }
    //ISSUE not sure about this type of interface for boundary groups; mb move it on lever of host_mesh?
    //TODO temporal solution (max 4 nodes) but will be enough for most cases
    bool check_face_has_group_id(Ord nodes_n, Ord prim_nodes[4])const
    {
        face_key_t  face_key(nodes_n, prim_nodes);
        return bnd_faces_group_ids_.find(face_key) != bnd_faces_group_ids_.end();
    }
    //TODO temporal solution (max 4 nodes) but will be enough for most cases
    Ord get_face_group_id(Ord nodes_n, Ord prim_nodes[4])const
    {
        face_key_t  face_key(nodes_n, prim_nodes);
        auto it = bnd_faces_group_ids_.find(face_key);
        if (it == bnd_faces_group_ids_.end())
            throw std::logic_error("gmsh_mesh_wrap::get_face_group_id: no face exists");
        return it->second;
    }
    bool check_elem_face_has_group_id(Ord elem_id, Ord face_i)const
    {
        face_key_t face_key(*this,elem_id,face_i);
        return bnd_faces_group_ids_.find(face_key) != bnd_faces_group_ids_.end();
    }
    Ord get_elem_face_group_id(Ord elem_id, Ord face_i)const
    {
        face_key_t face_key(*this,elem_id,face_i);
        auto it = bnd_faces_group_ids_.find(face_key);
        if (it == bnd_faces_group_ids_.end()) return special_id;
        return it->second;
    }

private:
    using elem_type_ord_t = elem_type_ordinal_type;
    using face_key_t = detail::face_key<Ord>;
    using face_key_less_func = detail::face_key_less_func<Ord>;
    /// Here pair's first is id of incident element, second - local node index inside this element
    using nodes_to_elems_graph_t = detail::ranges_sparse_arr<std::pair<Ord,Ord>,Ord>;

    using nodes_virt_master_ids_arr_t =  detail::sparse_arr<Ord>;
    //using faces_virt_master_ids_arr_t =  detail::sparse_arr<Ord>;

    using vec_t = static_vec::vec<scalar_type,dim>;
    using mat_t = static_mat::mat<scalar_type,dim,dim>;

private:
    mesh_elem_reference_type        mesh_elem_reference_;

    std::string                     fn_;
    /// Stick to old private C++ API
    std::shared_ptr<GModel>         g_model_;
    std::shared_ptr<PartElems>      elems_partitioner_;
    bool                            is_homogeneous_;
    elem_type_ordinal_type          homogeneous_elem_type_;
    Ord                             elements_index_shift_;
    Ord                             elems_max_prim_nodes_num_,
                                    elems_max_nodes_num_;
    nodes_to_elems_graph_t          nodes_to_elems_graph_;
    //TODO turn into sparse_arr?
    std::map<Ord,Ord>               elements_group_ids_;
    Ord                             elems_max_faces_num_;

    std::map<face_key_t,Ord,face_key_less_func>        bnd_faces_group_ids_;

    nodes_virt_master_ids_arr_t     nodes_virt_master_ids_arr_;

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

    //TODO O(n^2) algo is basically used here - use kd-tree instead of ref_all_nodes
    /// a and b pair gives affine transfrom
    /// graph is symmetric by construction
    void add_virt_nodes_graph_connections
    (
        const std::set<Ord> &nodes,
        const mat_t &a, const vec_t &b,
        const std::set<Ord> &ref_all_nodes,
        std::map<Ord,std::set<Ord>> &graph
    )
    {
        for (auto node_id : nodes)
        {
            vec_t  c, c_pair;
            get_node_coords(node_id,c.d);
            c_pair = a*c + b;
            bool    node_pair_id_found = false;
            Ord     node_pair_id;
            for (auto node1_id : ref_all_nodes)
            {
                vec_t  c1, c_diff;
                get_node_coords(node1_id,c1.d);
                c_diff = c - c1;
                bool is_found = true;
                for (int j = 0;j < dim;++j)
                    if (!(std::abs(c_diff[j]) <= std::numeric_limits<T>::epsilon()))
                        is_found = false;
                if (is_found)
                {
                    node_pair_id_found = true;
                    node_pair_id = node1_id;
                    break;
                }
            }
            if (!node_pair_id_found)
                throw 
                    std::runtime_error
                    (
                        "gmsh_mesh_wrap::add_virt_nodes_graph_connections: node pair for node " + 
                        std::to_string(node_id) + " was not found; perhaps periodicity in mesh is broken"
                    );
            graph[node_id].insert(node_pair_id);
            graph[node_pair_id].insert(node_id);
        }
    }
    /// Here we assert that f dimension is dim-1 (i.e. it is considered to be 'face')
    std::vector<GEntity*> get_face_bound_g_entities(int curr_dim, GEntity* f)
    {
        if (curr_dim == dim-1)
            return std::vector<GEntity*>({f});
        else
        {
            switch (curr_dim)
            {
                case 0: 
                    return f->vertices();
                case 1:
                    return f->edges();
                default:
                    throw std::logic_error("gmsh_mesh_wrap::get_face_bound_g_entities: dim > 2");
            };
        }
    }
    void build_virt_nodes(const std::set<Ord> &periodic_g_faces_tags)
    {
        std::set<Ord> master_periodic_g_faces_tags,
                      subordinates_periodic_g_faces_tags;

        /// Check for periodic_g_faces_tags logic integrity

        for (auto g_face_tag : periodic_g_faces_tags)
        {
            GEntity *f = g_model_->getEntityByTag(dim-1, g_face_tag);
            if (f->getMeshMaster()->tag() == f->tag())
            {
                master_periodic_g_faces_tags.insert(f->tag());
            }
            else 
            {
                subordinates_periodic_g_faces_tags.insert(f->tag());
            }
        }
        for (auto g_face_tag : periodic_g_faces_tags)
        {
            GEntity *f = g_model_->getEntityByTag(dim-1, g_face_tag);
            if (f->getMeshMaster()->tag() == f->tag()) continue;
            if (master_periodic_g_faces_tags.find(f->getMeshMaster()->tag()) == master_periodic_g_faces_tags.end())
                throw 
                    std::runtime_error
                    (
                        "gmsh_mesh_wrap::build_virt_nodes: face " + std::to_string(f->tag()) + " is present in"
                        "periodic_g_faces_tags while its master " + std::to_string(f->getMeshMaster()->tag()) + 
                        " is not"
                    );
            /// ISSUE is it possible that on master surface is referenced by two other surfaces??
            master_periodic_g_faces_tags.erase(f->getMeshMaster()->tag());
        }
        if (!master_periodic_g_faces_tags.empty())
            throw 
                std::runtime_error
                (
                    "gmsh_mesh_wrap::build_virt_nodes: master face " + 
                    std::to_string(*master_periodic_g_faces_tags.begin()) + 
                    " is present in periodic_g_faces_tags while there is not pair for it"
                );

        /// Build nodes virtual connectivity graph
        /// First, build reference nodes set (must be replaced with kd-tree in future)
        std::set<Ord> ref_all_nodes;
        for (int curr_dim = dim-1;dim >= 0;--dim)
        {
            std::vector<GEntity *> entities;
            g_model_->getEntities(entities, curr_dim);
            for (auto entity : entities)
            {
                for (Ord vertex_i = 0;vertex_i < entity->getNumMeshVertices();++vertex_i)
                {
                    ref_all_nodes.insert(node_tag_to_node_id(entity->getMeshVertex(vertex_i)->getNum()));
                }
            }
        }
        /// Second, build graph itself
        std::map<Ord,std::set<Ord>> graph;
        for (auto g_face_tag : subordinates_periodic_g_faces_tags)
        {
            GEntity *f = g_model_->getEntityByTag(dim-1, g_face_tag);

            auto at = f->affineTransform;

            mat_t a(at[0].at[1],at[2],
                    at[4].at[5],at[6],
                    at[8].at[9],at[10]);
            vec_t b(at[3],at[7],at[11]);

            for (int curr_dim = dim-1;dim >= 0;--dim)
            {
                std::vector<GEntity*> entities = get_face_bound_g_entities(curr_dim, f);
                for (auto entity : entities)
                {
                    std::set<Ord> nodes;
                    for (Ord vertex_i = 0;vertex_i < entity->getNumMeshVertices();++vertex_i)
                    {
                        nodes.insert(node_tag_to_node_id(entity->getMeshVertex(vertex_i)->getNum()));
                    }
                    add_virt_nodes_graph_connections
                    (
                        nodes, a, b, ref_all_nodes, graph
                    );
                }
            }
        }

        //nodes_virt_master_ids_arr_
    }
};

}  /// namespace mesh
}  /// namespace scfd

#endif