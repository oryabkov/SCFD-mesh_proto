
#include <scfd/communication/linear_partitioner.h>
#include <scfd/mesh/gmsh_mesh_wrap.h>

using partitioner_t = scfd::communication::linear_partitioner;
using gmsh_wrap_t = scfd::mesh::gmsh_mesh_wrap<double,partitioner_t,3>;

int main(int argc, char const *args[])
{
    
    return 0;
}