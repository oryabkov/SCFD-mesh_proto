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

TEST(TestGMSHMeshWrap, BasicRead) 
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
        gmsh_wrap->get_elem_prim_nodes(0, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 160);
        ASSERT_EQ(nodes[1], 293);
        ASSERT_EQ(nodes[2], 189);
        ASSERT_EQ(nodes[3], 298);
        //ASSERT_EQ(nodes[3], 299);
        gmsh_wrap->get_elem_prim_nodes(1137, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 165);
        ASSERT_EQ(nodes[1], 307);
        ASSERT_EQ(nodes[2], 72);
        ASSERT_EQ(nodes[3], 92);
        gmsh_wrap->get_elem_prim_nodes(1704-621, nodes, &prim_nodes_num);
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
        std::set<ordinal>   elem_1704_face_groups;
        for (ordinal j = 0;j < 4;++j)
        {
            if (gmsh_wrap->check_elem_face_has_group_id(1704-621,j))
                elem_1704_face_groups.insert(gmsh_wrap->get_elem_face_group_id(1704-621,j));
        }
        ASSERT_EQ(elem_1704_face_groups.size(),2);
        ASSERT_EQ(elem_1704_face_groups,std::set<ordinal>({18,22}));

        real coords[3];
        gmsh_wrap->get_node_coords(160,coords);
        std::cout << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
        //1 0.8852084354774035 0.4191415113762025

    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

TEST(TestGMSHMeshWrap, BasicReadPeriodic1)
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        gmsh_wrap = std::make_shared<gmsh_wrap_t>();
        gmsh_wrap->set_mesh_filename("test_box3d_period_small_mesh.msh");
        /// No periodic surfaces is set (so, mesh is not periodic)
        gmsh_wrap->read();
        *part = partitioner_t(gmsh_wrap->get_total_elems_num(), 1, 0);
        gmsh_wrap->set_partitioner(part);

        ASSERT_EQ(gmsh_wrap->get_total_nodes_num(), 14);
        for (ordinal i = 1;i <= 14;++i)
        {
            ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(i), i);
        }

        ordinal prim_nodes_num; 
        ordinal nodes[4];
        gmsh_wrap->get_elem_prim_nodes(65-45, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 13);
        ASSERT_EQ(nodes[1], 8);
        ASSERT_EQ(nodes[2], 14);
        ASSERT_EQ(nodes[3], 12);

        ordinal prim_virt_nodes_num; 
        ordinal virt_nodes[4];
        gmsh_wrap->get_elem_prim_virt_nodes(65-45, virt_nodes, &prim_virt_nodes_num);
        ASSERT_EQ(prim_virt_nodes_num, 4);
        ASSERT_EQ(virt_nodes[0], 13);
        ASSERT_EQ(virt_nodes[1], 8);
        ASSERT_EQ(virt_nodes[2], 14);
        ASSERT_EQ(virt_nodes[3], 12);

        ordinal elems_num; 
        ordinal elems[4];
        ASSERT_EQ(gmsh_wrap->get_virt_node_incident_elems_num(8),4);
        gmsh_wrap->get_virt_node_incident_elems(8,elems,&elems_num);
        ASSERT_EQ(elems_num, 4);
        ///TODO order is in fact not garanteed here
        ASSERT_EQ(elems[0], 50-45);
        ASSERT_EQ(elems[1], 51-45);
        ASSERT_EQ(elems[2], 56-45);
        ASSERT_EQ(elems[3], 65-45);
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

TEST(TestGMSHMeshWrap, BasicReadPeriodic2)
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        gmsh_wrap = std::make_shared<gmsh_wrap_t>();
        gmsh_wrap->set_mesh_filename("test_box3d_period_small_mesh.msh");
        /// No periodic surfaces is set (so, mesh is not periodic)
        gmsh_wrap->read(std::set<ordinal>({5,27}));
        *part = partitioner_t(gmsh_wrap->get_total_elems_num(), 1, 0);
        gmsh_wrap->set_partitioner(part);

        ASSERT_EQ(gmsh_wrap->get_total_nodes_num(), 14);
        /*for (ordinal i = 1;i <= 14;++i)
        {
            ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(i), i);
        }*/
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(1 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(2 ), 2 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(3 ), 3 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(4 ), 4 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(5 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(6 ), 2 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(7 ), 4 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(8 ), 3 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(9 ), 9 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(10), 10);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(11), 11);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(12), 12);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(13), 13);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(14), 9 );

        ordinal prim_nodes_num; 
        ordinal nodes[4];
        gmsh_wrap->get_elem_prim_nodes(65-45, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 13);
        ASSERT_EQ(nodes[1], 8);
        ASSERT_EQ(nodes[2], 14);
        ASSERT_EQ(nodes[3], 12);

        ordinal prim_virt_nodes_num; 
        ordinal virt_nodes[4];
        gmsh_wrap->get_elem_prim_virt_nodes(65-45, virt_nodes, &prim_virt_nodes_num);
        ASSERT_EQ(prim_virt_nodes_num, 4);
        ASSERT_EQ(virt_nodes[0], 13);
        ASSERT_EQ(virt_nodes[1], 3);
        ASSERT_EQ(virt_nodes[2], 9);
        ASSERT_EQ(virt_nodes[3], 12);
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

/// Fully periodic case 
TEST(TestGMSHMeshWrap, BasicReadPeriodic3)
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        gmsh_wrap = std::make_shared<gmsh_wrap_t>();
        gmsh_wrap->set_mesh_filename("test_box3d_period_small_mesh.msh");
        /// No periodic surfaces is set (so, mesh is not periodic)
        gmsh_wrap->read(std::set<ordinal>({5,14,18,22,26,27}));
        *part = partitioner_t(gmsh_wrap->get_total_elems_num(), 1, 0);
        gmsh_wrap->set_partitioner(part);

        ASSERT_EQ(gmsh_wrap->get_total_nodes_num(), 14);
        /*for (ordinal i = 1;i <= 14;++i)
        {
            ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(i), i);
        }*/
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(1 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(2 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(3 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(4 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(5 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(6 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(7 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(8 ), 1 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(9 ), 9 );
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(10), 10);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(11), 11);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(12), 10);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(13), 11);
        ASSERT_EQ(gmsh_wrap->get_node_virt_master_id(14), 9 );

        ordinal prim_nodes_num; 
        ordinal nodes[4];
        gmsh_wrap->get_elem_prim_nodes(65-45, nodes, &prim_nodes_num);
        ASSERT_EQ(prim_nodes_num, 4);
        ASSERT_EQ(nodes[0], 13);
        ASSERT_EQ(nodes[1], 8);
        ASSERT_EQ(nodes[2], 14);
        ASSERT_EQ(nodes[3], 12);

        ordinal prim_virt_nodes_num; 
        ordinal virt_nodes[4];
        gmsh_wrap->get_elem_prim_virt_nodes(65-45, virt_nodes, &prim_virt_nodes_num);
        ASSERT_EQ(prim_virt_nodes_num, 4);
        ASSERT_EQ(virt_nodes[0], 11);
        ASSERT_EQ(virt_nodes[1], 1);
        ASSERT_EQ(virt_nodes[2], 9);
        ASSERT_EQ(virt_nodes[3], 10);
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}
