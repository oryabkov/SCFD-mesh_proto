// Copyright © 2016-2021 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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

#ifndef __SCFD_MESH_SPARSE_ARR_H__
#define __SCFD_MESH_SPARSE_ARR_H__

#include <map>

namespace scfd
{
namespace mesh
{
namespace detail
{

template<class T,class Ord = int>
struct sparse_arr
{
    Ord                  n;
    std::vector<T>       objs_;
    std::map<Ord,Ord>    glob_ind_2_vec_ind;

    //it's real total mesh size, not number of actually read elements (which could be zero, for example)
    Ord             size()const 
    { 
        return n; 
    }
    Ord             find_elem__(Ord glob_i)const
    {
        std::map<Ord,Ord>::const_iterator       it = glob_ind_2_vec_ind.find(glob_i);
        assert(it != glob_ind_2_vec_ind.end());
        return it->second;
    }
    T            &operator[](Ord glob_i)
    {
        //ISSUE check if we have glob_i?? i think do it in DEBUG mode
        return objs_[ find_elem__(glob_i) ];
    }
    const T      &operator[](Ord glob_i)const
    {
        //ISSUE check if we have glob_i?? i think do it in DEBUG mode
        return objs_[ find_elem__(glob_i) ];
    }
    //TODO rename has_elem and add_elem into has and and
    bool            has_elem(Ord glob_i)const
    {
        return (glob_ind_2_vec_ind.find(glob_i) != glob_ind_2_vec_ind.end());
    }
    void            add_elem(Ord glob_i, const T &e)
    {
        //ISSUE check if we already have glob_i?? i think do it in DEBUG mode
        glob_ind_2_vec_ind[glob_i] = objs_.size();
        objs_.push_back(e);
    }
};

}  /// namespace detail
}  /// namespace mesh
}  /// namespace scfd

#endif
