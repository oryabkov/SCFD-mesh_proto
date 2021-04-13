
#include <iostream>
#include <scfd/communication/linear_partitioner.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>

using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<double,partitioner_t,3>;

int main(int argc, char const *args[])
{
    try 
    {
        auto        part = std::make_shared<partitioner_t>();
        auto        gmsh_wrap = std::make_shared<gmsh_wrap_t>();
        gmsh_wrap->set_mesh_filename("test.msh");
        gmsh_wrap->read();
        *part = partitioner_t(gmsh_wrap->get_total_elems_num(), 1, 0);
        gmsh_wrap->set_partitioner(part);

        std::cout << "gmsh_wrap->get_total_elems_num() = " << gmsh_wrap->get_total_elems_num() << std::endl;
        std::cout << "expected number " << 1138 << std::endl;
    } 
    catch(const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "exit" << std::endl;
    }

    return 0;
}