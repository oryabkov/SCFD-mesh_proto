// Copyright © 2016 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

// This file is part of SimpleCFD.

// SimpleCFD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 only of the License.

// SimpleCFD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with SimpleCFD.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __SCFD_MPI_COMMUNICATOR_H__
#define __SCFD_MPI_COMMUNICATOR_H__

#include <stdexcept>
#include <mpi.h>

namespace scfd
{
namespace communication
{

class mpi_communicator
{
public:
    int size()const
    {
        int res;
        if (MPI_Comm_size(MPI_COMM_WORLD, &res) != MPI_SUCCESS) throw std::runtime_error("mpi_communicator::size: MPI_Comm_size failed");
        return res;
    }
    int my_rank()const
    {
        int res;
        if (MPI_Comm_rank(MPI_COMM_WORLD, &res) != MPI_SUCCESS) throw std::runtime_error("mpi_communicator::my_rank: MPI_Comm_rank failed");
        return res;
    }

    template<class T>
    T   reduce_max(const T &local_val)const
    {
        throw std::logic_error("mpi_communicator::reduce_max: not supported for given type");
    }
    template<class T>
    T   reduce_sum(const T &local_val)const
    {
        throw std::logic_error("mpi_communicator::reduce_sum: not supported for given type");
    }
    //TODO MPI handle
};

template<>
float   mpi_communicator::reduce_max<float>(const float &local_val)const
{
    float       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
    return res_val;
}

template<>
double   mpi_communicator::reduce_max<double>(const double &local_val)const
{
    double       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    return res_val;
}

template<>
int      mpi_communicator::reduce_max<int>(const int &local_val)const
{
    int       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    return res_val;
}

template<>
float   mpi_communicator::reduce_sum<float>(const float &local_val)const
{
    float       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
    return res_val;
}

template<>
double   mpi_communicator::reduce_sum<double>(const double &local_val)const
{
    double       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return res_val;
}

template<>
int      mpi_communicator::reduce_sum<int>(const int &local_val)const
{
    int       res_val;
    MPI_Allreduce( &local_val, &res_val, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    return res_val;
}

}  /// namespace communication
}  /// namespace scfd

#endif
