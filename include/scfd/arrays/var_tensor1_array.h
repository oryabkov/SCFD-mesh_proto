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
    template<class T1, class MemoryIn, class MemoryOut>
    friend void copy
    (
        const var_tensor1_array<T1,MemoryIn> &arr_in,
        var_tensor1_array<T1,MemoryOut> &arr_out
    );
public:
    using value_type = T;
    using pointer_type = T*;
    using ordinal_type = arrays::ordinal_type;
    using memory_type = Memory;

public:
    void pre_init(ordinal_type size0,ordinal_type i0_0)
    {
        try 
        {
            offsets_.init(size0,i0_0);
            lens_.init(size0,i0_0);
        }
        catch (...)
        {
            offsets_.free();
            lens_.free();
            std::throw_with_nested( std::runtime_error("arrays::var_tensor1_array::pre_init: failed memory allocation") );
        }
        /// TODO general (GPU/CPU) implementation
        for (ordinal_type i = 0;i < lens_.size();++i)
        {
            lens_(i) = 0;
        }
    }
    __DEVICE_TAG__ void set_var_tensor_dim(ordinal_type i0,ordinal_type tensor_dim)
    {
        lens_(i0) = tensor_dim;
    }
    void init(bool create_keys_array)
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

        if (create_keys_array)
        {
            //TODO create create_keys_array
        }
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

    /// Direct access to internal arrays; TODO some kind of const_arrays would be usefull here

    __DEVICE_TAG__ array<ordinal_type,Memory> offsets()const
    {
        return offsets_;
    }
    __DEVICE_TAG__ array<ordinal_type,Memory> lens()const
    {
        return lens_;
    }
    __DEVICE_TAG__ array<T,Memory>            values()const
    {
        return values_;
    }
    __DEVICE_TAG__ array<ordinal_type,Memory> keys()const
    {
        return keys_;
    }
    
private:
    array<ordinal_type,Memory>  offsets_, lens_;
    array<T,Memory>             values_;
    /// This one is optional and has the same dimension as values_
    array<ordinal_type,Memory>  keys_;

};

template<class T, class MemoryIn, class MemoryOut>
void copy
(
    const var_tensor1_array<T,MemoryIn> &arr_in,
    var_tensor1_array<T,MemoryOut> &arr_out
)
{
    //TODO check and allow case for equal sizes
    if (!arr_out.is_free() && !arr_out.is_own())
    {
        throw 
            std::runtime_error
            (
                "arrays::copy: try to copy to non-empty non-owning array"
            );
    }

    arr_out.free();

}

} /// namespace arrays
} /// namespace scfd

#endif
