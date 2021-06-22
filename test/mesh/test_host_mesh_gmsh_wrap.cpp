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
#include <scfd/mesh/host_mesh.h>

using real = double;
using ordinal = int;
using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<real,partitioner_t,3,ordinal>;
using host_mesh_t = scfd::mesh::host_mesh<gmsh_wrap_t>;

TEST(TestHostMeshGMSHWrap, BasicRead) 
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        host_mesh = std::make_shared<host_mesh_t>();
        host_mesh->set_mesh_filename("test.msh");
        host_mesh->read();
        *part = partitioner_t(host_mesh->get_total_elems_num(), 1, 0);
        host_mesh->set_partitioner(part);
        host_mesh->enlarge_stencil(1);

        ASSERT_EQ(host_mesh->get_total_elems_num(), 1138);
        ASSERT_EQ(host_mesh->get_elems_max_faces_num(), 4);
        ASSERT_EQ(host_mesh->get_elems_glob_max_faces_num(), 4);
        ASSERT_EQ(host_mesh->get_elems_max_nodes_num(), 4);
        ASSERT_EQ(host_mesh->get_elem_prim_nodes_num(0), 4);
        ordinal prim_nodes_num; 
        ordinal nodes[4];
        host_mesh->get_elem_prim_nodes(0, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 160);
        ASSERT_EQ(nodes[1], 293);
        ASSERT_EQ(nodes[2], 189);
        ASSERT_EQ(nodes[3], 298);
        //ASSERT_EQ(nodes[3], 299);
        host_mesh->get_elem_prim_nodes(1137, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 165);
        ASSERT_EQ(nodes[1], 307);
        ASSERT_EQ(nodes[2], 72);
        ASSERT_EQ(nodes[3], 92);
        host_mesh->get_elem_prim_nodes(1704-621, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 62);
        ASSERT_EQ(nodes[1], 140);
        ASSERT_EQ(nodes[2], 172);
        ASSERT_EQ(nodes[3], 61);
        
        ASSERT_EQ(host_mesh->get_elem_type(0),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(500),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(1137),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_group_id(0),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(500),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(1137),1);
        //test neigbours0 1704-621 1340-621
        std::set<ordinal>   elem_1704_neibs0_set;
        ordinal             elem_1704_neibs0[host_mesh->get_elem_faces_num(1704-621)];
        host_mesh->get_elem_neighbours0(1704-621, elem_1704_neibs0);
        for (ordinal j = 0;j < host_mesh->get_elem_faces_num(1704-621);++j)
        {
            elem_1704_neibs0_set.insert(elem_1704_neibs0[j]);
        }
        ASSERT_TRUE(elem_1704_neibs0_set.find(1340-621) != elem_1704_neibs0_set.end());
        //boundary triangles of element 1704-621
        //triangle 220 -> surface 18
        //triangle 301 -> surface 22
        //We know that this is tetrahedron
        std::set<ordinal>   elem_1704_face_groups;
        for (ordinal j = 0;j < 4;++j)
        {
            if (host_mesh->check_elem_face_has_group_id(1704-621,j))
                elem_1704_face_groups.insert(host_mesh->get_elem_face_group_id(1704-621,j));
        }
        ASSERT_EQ(elem_1704_face_groups.size(),2);
        ASSERT_EQ(elem_1704_face_groups,std::set<ordinal>({18,22}));

    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

TEST(TestHostMeshGMSHWrap, BasicReadPeriodic1)
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        host_mesh = std::make_shared<host_mesh_t>();
        host_mesh->set_mesh_filename("test_box3d_period_small_mesh.msh");
        /// No periodic surfaces is set (so, mesh is not periodic)
        host_mesh->read();
        *part = partitioner_t(host_mesh->get_total_elems_num(), 1, 0);
        host_mesh->set_partitioner(part);
        host_mesh->enlarge_stencil(1);

        ASSERT_EQ(host_mesh->get_total_nodes_num(), 14);
        for (ordinal i = 1;i <= 14;++i)
        {
            ASSERT_EQ(host_mesh->get_node_virt_master_id(i), i);
        }

        ordinal prim_nodes_num; 
        ordinal nodes[4];
        host_mesh->get_elem_prim_nodes(65-45, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 13);
        ASSERT_EQ(nodes[1], 8);
        ASSERT_EQ(nodes[2], 14);
        ASSERT_EQ(nodes[3], 12);

        ordinal prim_virt_nodes_num; 
        ordinal virt_nodes[4];
        host_mesh->get_elem_prim_virt_nodes(65-45, virt_nodes, &prim_virt_nodes_num);
        ASSERT_EQ(prim_virt_nodes_num, 4);
        ASSERT_EQ(virt_nodes[0], 13);
        ASSERT_EQ(virt_nodes[1], 8);
        ASSERT_EQ(virt_nodes[2], 14);
        ASSERT_EQ(virt_nodes[3], 12);

        ordinal elems_num0; 
        ordinal elems0[4];
        ASSERT_EQ(host_mesh->get_node_incident_elems_num(8),4);
        host_mesh->get_node_incident_elems(8,elems0,&elems_num0);
        ASSERT_EQ(elems_num0, 4);
        ///TODO order is in fact not garanteed here
        ASSERT_EQ(elems0[0], 50-45);
        ASSERT_EQ(elems0[1], 51-45);
        ASSERT_EQ(elems0[2], 56-45);
        ASSERT_EQ(elems0[3], 65-45);

        ordinal elems_num; 
        ordinal elems[4];
        ASSERT_EQ(host_mesh->get_virt_node_incident_elems_num(8),4);
        host_mesh->get_virt_node_incident_elems(8,elems,&elems_num);
        ASSERT_EQ(elems_num, 4);
        ///TODO order is in fact not garanteed here
        ASSERT_EQ(elems[0], 50-45);
        ASSERT_EQ(elems[1], 51-45);
        ASSERT_EQ(elems[2], 56-45);
        ASSERT_EQ(elems[3], 65-45);

        ASSERT_EQ(host_mesh->get_elem_type(54-45),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(64-45),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(67-45),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(50-45),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_type(60-45),MSH_TET_4);
        ASSERT_EQ(host_mesh->get_elem_group_id(54-45),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(64-45),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(67-45),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(50-45),1);
        ASSERT_EQ(host_mesh->get_elem_group_id(60-45),1);
        //test neigbours0 54-45 
        std::set<ordinal>   elem_54_neibs0_set;
        ordinal             elem_54_neibs0[host_mesh->get_elem_faces_num(54-45)];
        host_mesh->get_elem_neighbours0(54-45, elem_54_neibs0);
        for (ordinal j = 0;j < host_mesh->get_elem_faces_num(54-45);++j)
        {
            elem_54_neibs0_set.insert(elem_54_neibs0[j]);
        }
        std::set<ordinal>   elem_54_neibs0_ref_set({64-45,67-45});
        ASSERT_TRUE(elem_54_neibs0_set == elem_54_neibs0_ref_set);
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

