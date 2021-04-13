// Copyright © 2016-2021 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

// This file is part of SCFD.

// SCFD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 only of the License.

// SCFD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with SCFD.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include "gtest/gtest.h"
#include <scfd/communication/linear_partitioner.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>

using real = double;
using ordinal = int;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,3,ordinal>;

TEST(GMSHMeshWrapTest, BasicRead) 
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        gmsh_wrap = std::make_shared<gmsh_wrap_t>();
        gmsh_wrap->set_mesh_filename("test.msh");
        gmsh_wrap->read();
        *part = partitioner_t(gmsh_wrap->get_total_elems_num(), 1, 0);
        gmsh_wrap->set_partitioner(part);

        ASSERT_EQ(gmsh_wrap->get_total_elems_num(), 1138);
        ASSERT_EQ(gmsh_wrap->get_elems_max_faces_num(), 4);
        ASSERT_EQ(gmsh_wrap->get_elems_glob_max_faces_num(), 4);
        ASSERT_EQ(gmsh_wrap->get_elems_max_nodes_num(), 4);
        ASSERT_EQ(gmsh_wrap->get_elem_prim_nodes_num(0), 4);
        ordinal prim_nodes_num; 
        ordinal nodes[4];
        gmsh_wrap->get_elem_prim_nodes(0, &prim_nodes_num, nodes);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 160);
        ASSERT_EQ(nodes[1], 293);
        ASSERT_EQ(nodes[2], 189);
        ASSERT_EQ(nodes[3], 298);
        //ASSERT_EQ(nodes[3], 299);
        gmsh_wrap->get_elem_prim_nodes(1137, &prim_nodes_num, nodes);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 165);
        ASSERT_EQ(nodes[1], 307);
        ASSERT_EQ(nodes[2], 72);
        ASSERT_EQ(nodes[3], 92);
        gmsh_wrap->get_elem_prim_nodes(1704-621, &prim_nodes_num, nodes);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 62);
        ASSERT_EQ(nodes[1], 140);
        ASSERT_EQ(nodes[2], 172);
        ASSERT_EQ(nodes[3], 61);
        
        ASSERT_EQ(gmsh_wrap->get_elem_type(0),MSH_TET_4);
        ASSERT_EQ(gmsh_wrap->get_elem_type(500),MSH_TET_4);
        ASSERT_EQ(gmsh_wrap->get_elem_type(1137),MSH_TET_4);
        ASSERT_EQ(gmsh_wrap->get_elem_group_id(0),1);
        ASSERT_EQ(gmsh_wrap->get_elem_group_id(500),1);
        ASSERT_EQ(gmsh_wrap->get_elem_group_id(1137),1);
        //test neigbours0 1704-621 1340-621
        //boundary triangles of element 1704-621
        //triangle 220 -> surface 18
        //triangle 301 -> surface 22
        //We know that this is tetrahedron
        std::set<ordinal>   elem_1704_face_groups,
                            elem_1704_face_groups_ref = {18,22};
        for (ordinal j = 0;j < 4;++j)
        {
            if (gmsh_wrap->check_elem_face_has_group_id(1704-621,j))
                elem_1704_face_groups.insert(gmsh_wrap->get_elem_face_group_id(1704-621,j));
        }
        ASSERT_EQ(elem_1704_face_groups.size(),2);
        ASSERT_EQ(elem_1704_face_groups,elem_1704_face_groups_ref);

    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}
