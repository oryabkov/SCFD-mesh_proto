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

#ifndef __SCFD_STATIC_MAT_BOOL_ARRAY_H__
#define __SCFD_STATIC_MAT_BOOL_ARRAY_H__

namespace scfd
{
namespace static_mat
{
namespace detail
{

template< bool ... b> struct bool_array{};
template< bool ... b> struct check_all_are_true: std::is_same< bool_array<b...>, bool_array<(b||true)...> >{};

}
}
}

#endif
