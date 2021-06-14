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

#include <gmsh/GmshGlobal.h>
#include <gmsh/GModel.h>

TEST(TestGMSHCPPAPIGeom, BasicReadGEO) 
{
    try 
    {
        auto  g_model_ = std::make_shared<GModel>();;

        std::string fn = "test.geo";
        /// Read mesh
        if (!g_model_->readMSH(fn)) 
        {
            throw file_read_error("gmsh_mesh_wrap::read(): ",fn);
        }

        ASSERT_EQ(g_model_->getNumRegions(), 1138);

        
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
        FAIL();
    }    
}
