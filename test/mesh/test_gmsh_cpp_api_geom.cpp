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

        //std::string fn = "test.msh";
        std::string fn = "test_box3d_small_mesh.msh";
        /// Read mesh
        if (!g_model_->readMSH(fn)) 
        {
            throw std::runtime_error("gmsh_mesh_wrap::read(): " + fn);
        }

        ASSERT_EQ(g_model_->getNumRegions(), 1);
        ASSERT_EQ(g_model_->getNumFaces(), 6);
        ASSERT_EQ(g_model_->getNumEdges(), 12);
        ASSERT_EQ(g_model_->getNumVertices(), 8);
         
        GRegion*  r = *(g_model_->firstRegion());
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

        //const GFace*  f = g_model_->getFaceByTag(5);
        //NOTE const contradicts with getMeshVertex later here
        GFace*  f = g_model_->getFaceByTag(5);

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

        ASSERT_EQ(f->getNumMeshVertices(),1);
        ASSERT_EQ(f->getMeshVertex(0)->getNum(),9);

        ASSERT_EQ(r->getNumMeshVertices(),0);
    
        GVertex*  v = g_model_->getVertexByTag(14);
        ASSERT_EQ(v->getNumMeshVertices(),1);
        ASSERT_EQ(v->getMeshVertex(0)->getNum(),8);
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}

TEST(TestGMSHCPPAPIGeom, BasicReadMSHPeriod) 
{
    try 
    {
        auto  g_model_ = std::make_shared<GModel>();;

        //std::string fn = "test.msh";
        std::string fn = "test_box3d_period_small_mesh.msh";
        /// Read mesh
        if (!g_model_->readMSH(fn)) 
        {
            throw std::runtime_error("gmsh_mesh_wrap::read(): " + fn);
        }

        ASSERT_EQ(g_model_->getNumRegions(), 1);
        ASSERT_EQ(g_model_->getNumFaces(), 6);
        ASSERT_EQ(g_model_->getNumEdges(), 12);
        ASSERT_EQ(g_model_->getNumVertices(), 8);
         
        GRegion*  r = *(g_model_->firstRegion());
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

        //const GFace*  f = g_model_->getFaceByTag(5);
        //NOTE const contradicts with getMeshVertex later here
        GFace*  f = g_model_->getFaceByTag(5);

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

        ASSERT_EQ(f->getNumMeshVertices(),1);
        ASSERT_EQ(f->getMeshVertex(0)->getNum(),9);

        ASSERT_EQ(r->getNumMeshVertices(),0);
    
        GVertex*  v = g_model_->getVertexByTag(14);
        ASSERT_EQ(v->getNumMeshVertices(),1);
        ASSERT_EQ(v->getMeshVertex(0)->getNum(),8);

        GFace   *f5 = g_model_->getFaceByTag(5),
                *f14 = g_model_->getFaceByTag(14),
                *f18 = g_model_->getFaceByTag(18),
                *f22 = g_model_->getFaceByTag(22),
                *f26 = g_model_->getFaceByTag(26),
                *f27 = g_model_->getFaceByTag(27);

        ASSERT_EQ(f18->getMeshMaster()->tag(),26);
        ASSERT_EQ(f22->getMeshMaster()->tag(),14);
        ASSERT_EQ(f27->getMeshMaster()->tag(),5);

        ASSERT_EQ(f26->getMeshMaster()->tag(),26);
        ASSERT_EQ(f14->getMeshMaster()->tag(),14);
        ASSERT_EQ(f5->getMeshMaster()->tag(),5);

        auto at = f18->affineTransform;

        /*std::cout << at[0] << " " << at[1] << " " << at[2] << " " << at[3] << std::endl
                  << at[4] << " " << at[5] << " " << at[6] << " " << at[7] << std::endl
                  << at[8] << " " << at[9] << " " << at[10] << " " << at[11] << std::endl
                  << at[12] << " " << at[13] << " " << at[14] << " " << at[15] << std::endl;

        std::cout << std::endl;*/

        EXPECT_DOUBLE_EQ(at[0], 1); EXPECT_DOUBLE_EQ(at[1], 0); EXPECT_DOUBLE_EQ(at[2], 0); EXPECT_DOUBLE_EQ(at[3], 1);
        EXPECT_DOUBLE_EQ(at[4], 0); EXPECT_DOUBLE_EQ(at[5], 1); EXPECT_DOUBLE_EQ(at[6], 0); EXPECT_DOUBLE_EQ(at[7], 0);
        EXPECT_DOUBLE_EQ(at[8], 0); EXPECT_DOUBLE_EQ(at[9], 0); EXPECT_DOUBLE_EQ(at[10],1); EXPECT_DOUBLE_EQ(at[11],0);
        EXPECT_DOUBLE_EQ(at[12],0); EXPECT_DOUBLE_EQ(at[13],0); EXPECT_DOUBLE_EQ(at[14],0); EXPECT_DOUBLE_EQ(at[15],1);

        at = f22->affineTransform;

        /*std::cout << at[0] << " " << at[1] << " " << at[2] << " " << at[3] << std::endl
                  << at[4] << " " << at[5] << " " << at[6] << " " << at[7] << std::endl
                  << at[8] << " " << at[9] << " " << at[10] << " " << at[11] << std::endl
                  << at[12] << " " << at[13] << " " << at[14] << " " << at[15] << std::endl;

        std::cout << std::endl;*/

        EXPECT_DOUBLE_EQ(at[0], 1); EXPECT_DOUBLE_EQ(at[1], 0); EXPECT_DOUBLE_EQ(at[2], 0); EXPECT_DOUBLE_EQ(at[3], 0);
        EXPECT_DOUBLE_EQ(at[4], 0); EXPECT_DOUBLE_EQ(at[5], 1); EXPECT_DOUBLE_EQ(at[6], 0); EXPECT_DOUBLE_EQ(at[7], 1);
        EXPECT_DOUBLE_EQ(at[8], 0); EXPECT_DOUBLE_EQ(at[9], 0); EXPECT_DOUBLE_EQ(at[10],1); EXPECT_DOUBLE_EQ(at[11],0);
        EXPECT_DOUBLE_EQ(at[12],0); EXPECT_DOUBLE_EQ(at[13],0); EXPECT_DOUBLE_EQ(at[14],0); EXPECT_DOUBLE_EQ(at[15],1);

        at = f27->affineTransform;

        /*std::cout << at[0] << " " << at[1] << " " << at[2] << " " << at[3] << std::endl
                  << at[4] << " " << at[5] << " " << at[6] << " " << at[7] << std::endl
                  << at[8] << " " << at[9] << " " << at[10] << " " << at[11] << std::endl
                  << at[12] << " " << at[13] << " " << at[14] << " " << at[15] << std::endl;

        std::cout << std::endl;*/

        EXPECT_DOUBLE_EQ(at[0], 1); EXPECT_DOUBLE_EQ(at[1], 0); EXPECT_DOUBLE_EQ(at[2], 0); EXPECT_DOUBLE_EQ(at[3], 0);
        EXPECT_DOUBLE_EQ(at[4], 0); EXPECT_DOUBLE_EQ(at[5], 1); EXPECT_DOUBLE_EQ(at[6], 0); EXPECT_DOUBLE_EQ(at[7], 0);
        EXPECT_DOUBLE_EQ(at[8], 0); EXPECT_DOUBLE_EQ(at[9], 0); EXPECT_DOUBLE_EQ(at[10],1); EXPECT_DOUBLE_EQ(at[11],1);
        EXPECT_DOUBLE_EQ(at[12],0); EXPECT_DOUBLE_EQ(at[13],0); EXPECT_DOUBLE_EQ(at[14],0); EXPECT_DOUBLE_EQ(at[15],1);
        
        //{18} = {26} Translate {1, 0, 0};
        //{22} = {14} Translate {0, 1, 0};
        //{27} = {5} Translate {0, 0, 1};
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}
