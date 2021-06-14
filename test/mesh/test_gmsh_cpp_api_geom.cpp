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
#include <stdexcept>
#include <set>
#include "gtest/gtest.h"
#include <gmsh/GmshGlobal.h>
#include <gmsh/GModel.h>

/*TEST(TestGMSHCPPAPIGeom, BasicReadGEO) 
{
    try 
    {
        //auto  g_model_ = std::make_shared<GModel>();;

        std::string fn = "test.geo";
        /// Read mesh
        if (!GModel::readGEO(fn)) 
        {
            throw std::runtime_error("gmsh_mesh_wrap::read(): " + fn);
        }

        auto g_model_ = GModel::current();

        //ASSERT_EQ(g_model_->getNumRegions(), 1);
        ASSERT_EQ(g_model_->getNumFaces(), 6);
        ASSERT_EQ(g_model_->getNumEdges(), 12);
        ASSERT_EQ(g_model_->getNumVertices(), 8);
         

        
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}*/

TEST(TestGMSHCPPAPIGeom, BasicReadMSH) 
{
    try 
    {
        auto  g_model_ = std::make_shared<GModel>();;

        std::string fn = "test.msh";
        /// Read mesh
        if (!g_model_->readMSH(fn)) 
        {
            throw std::runtime_error("gmsh_mesh_wrap::read(): " + fn);
        }

        ASSERT_EQ(g_model_->getNumRegions(), 1);
        ASSERT_EQ(g_model_->getNumFaces(), 6);
        ASSERT_EQ(g_model_->getNumEdges(), 12);
        ASSERT_EQ(g_model_->getNumVertices(), 8);
         
        const GRegion*  r = *(g_model_->firstRegion());
        ASSERT_EQ(r->faces().size(),6);
        std::set<int> faces_tags;
        for (auto face : r->faces())
        {
            faces_tags.insert(face->tag());
        }
        ASSERT_EQ(faces_tags, std::set<int>({5,14,18,22,26,27}));

        ASSERT_EQ(r->edges().size(),12);
        ASSERT_EQ(r->vertices().size(),8);
        //TODO add tags checks

        const GFace*  f = g_model_->getFaceByTag(5);

        ASSERT_EQ(f->edges().size(),4);
        std::set<int> edges_tags;
        for (auto edge : f->edges())
        {
            edges_tags.insert(edge->tag());
        }
        ASSERT_EQ(edges_tags, std::set<int>({1,2,3,4}));

        ASSERT_EQ(f->vertices().size(),4);
        std::set<int> vert_tags;
        for (auto vert : f->vertices())
        {
            vert_tags.insert(vert->tag());
        }
        ASSERT_EQ(vert_tags, std::set<int>({1,2,3,4}));

        const GEdge*  e = g_model_->getEdgeByTag(7);
        ASSERT_EQ(e->vertices().size(),2);
        std::set<int> edge_vert_tags;
        for (auto vert : e->vertices())
        {
            edge_vert_tags.insert(vert->tag());
        }
        ASSERT_EQ(edge_vert_tags, std::set<int>({5,6}));        

    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}
