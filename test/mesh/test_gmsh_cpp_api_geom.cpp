
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
