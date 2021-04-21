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
    
private:
    array<ordinal_type,Memory>  offsets_, lens_;
    array<T,Memory>             values_;

};

} /// namespace arrays
} /// namespace scfd

#endif
