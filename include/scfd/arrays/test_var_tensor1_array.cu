
#include <scfd/utils/init_cuda.h>
#include <scfd/memory/host.h>
#include <scfd/memory/cuda.h>
#include <scfd/arrays/var_tensor1_array.h>

using namespace scfd;

int main(int argc, char const *args[])
{
    arrays::var_tensor1_array<double,memory::cuda_host>     host_arr;
    arrays::var_tensor1_array<double,memory::cuda_device>   cuda_arr;


    
    return 0;
}