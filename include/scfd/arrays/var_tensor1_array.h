// Copyright © 2016-2020 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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

#ifndef __SCFD_ARRAYS_VAR_TENSOR1_ARRAY_H__
#define __SCFD_ARRAYS_VAR_TENSOR1_ARRAY_H__

#include <scfd/utils/device_tag.h>
#include "array.h"

namespace scfd
{
namespace arrays
{

template<class T, class Memory>
class var_tensor1_array
{
public:
    using value_type = T;
    using pointer_type = T*;
    using ordinal_type = arrays::ordinal_type;
    using memory_type = Memory;

public:
    void pre_init(ordinal_type size0)
    {
        try 
        {
            offsets_.init(size0);
            lens_.init(size0);
        }
        catch (...)
        {
            offsets_.free();
            lens_.free();
            std::throw_with_nested( std::runtime_error("var_tensor1_array::pre_init: failed memory allocation") );
        }
    }
    __DEVICE_TAG__ void set_var_tensor_dim(ordinal_type i0,ordinal_type tensor_dim)
    {
        lens_(i0) = tensor_dim;
    }
    void init()
    {
        if (lens_.size() == 0) return;
        /// Calc offsets_ using prefix sum
        /// TODO general (GPU/CPU) implementation
        ordinal_type curr_offset = 0;
        for (ordinal_type i = 0;i < lens_.size();curr_offset += lens_(i++))
        {
            offsets_(i) = curr_offset;
        }
        /// allocate array for values
        values_.init(curr_offset);
    }
    void free()
    {
        offsets_.free();
        lens_.free();
        values_.free();
    }
    __DEVICE_TAG__ T &operator()(ordinal_type i0,ordinal_type i1)const
    {
        ordinal_type    offset = offsets_(i0);
        return values_(offset);
    }  
    
private:
    array<ordinal_type,Memory>  offsets_, lens_;
    array<T,Memory>             values_;

};

} /// namespace arrays
} /// namespace scfd

#endif
