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

#ifndef __SCFD_PARTED_MAP_H__
#define __SCFD_PARTED_MAP_H__

#include <stdio.h>
#include <cassert>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <boost/foreach.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/find.hpp>

namespace communication
{

//supposed to satisfy MAP concept
//PARTITIONER must satisfy PARTITIONER concept
template<class PARTITIONER>
struct parted_map
{
    PARTITIONER             part;
    bool                    r_stencil_only;
    std::map<int,int>       stencil_glob_2_loc;
    std::vector<int>        l_stencil_loc_2_glob;
    std::vector<int>        r_stencil_loc_2_glob;

    //creates uninitialized map
    parted_map() {}
    //dummy constructor
    parted_map(const PARTITIONER &_part,bool _r_stencil_only = false) : part(_part), r_stencil_only(_r_stencil_only)
    {                
    }

    //'read' before construction part:
    int     get_own_rank()const { return part.get_own_rank(); }
    bool    check_glob_owned(int i)const { return part.check_glob_owned(i); }
    int     get_total_size()const { return part.get_total_size(); }
    int     get_size()const { return part.get_size(); }
    int     own_glob_ind(int i)const { return part.own_glob_ind(i); }

    //'construction' part:
    //t_simple_map just ignores this information
    //i is global index here
    void    add_stencil_element(int i)
    {
        if (check_glob_owned(i)) return;
        part.add_stencil_element(i);
        //before complete() stencil_glob_2_loc simply store all stencil elements global indices
        stencil_glob_2_loc[i] = 0;
    }
    void    complete()
    {
        part.complete();
        int median_elem = 0;
        if (get_size() > 0) {
            if (!r_stencil_only) median_elem = own_glob_ind(0); else median_elem = 0;
        }
        typedef std::pair<const int,int> pair_int_int;
        BOOST_FOREACH(pair_int_int &e, stencil_glob_2_loc) {
            if (e.first < median_elem) {
                //to the left
                l_stencil_loc_2_glob.push_back(e.first);
            } else {
                //to the right
                r_stencil_loc_2_glob.push_back(e.first);
            } 
        }
        boost::sort(l_stencil_loc_2_glob);
        boost::sort(r_stencil_loc_2_glob);
        for (int i = 0;i < l_stencil_loc_2_glob.size();++i) {
            int     glob_ind = l_stencil_loc_2_glob[i],
                loc_ind = -l_stencil_loc_2_glob.size() + i;
            stencil_glob_2_loc[glob_ind] = loc_ind;
        }
        for (int i = 0;i < r_stencil_loc_2_glob.size();++i) {
            int     glob_ind = r_stencil_loc_2_glob[i],
                loc_ind = get_size() + i;
            stencil_glob_2_loc[glob_ind] = loc_ind;
        }
    }

    //'read' after construction part:
    int     get_rank(int i)const { return part.get_rank(i); }
    bool    check_glob_owned(int i, int rank)const { return part.check_glob_owned(i, rank); }
    int     loc2glob(int i_loc)const
    {
        if (i_loc < 0) {
            return l_stencil_loc_2_glob[i_loc + l_stencil_loc_2_glob.size()];
        } else if (i_loc >= get_size()) {
            return r_stencil_loc_2_glob[i_loc - get_size()];
        } else {
            return part.own_glob_ind(i_loc);
        }
    }
    int     glob2loc(int i_glob)const
    {
        if (check_glob_owned(i_glob))
            return part.own_glob_ind_2_ind(i_glob);
        else
            return stencil_glob_2_loc.at(i_glob);
            //return stencil_glob_2_loc[i_glob];
    }
    int     own_loc_ind(int i)const
    {
        return i;
    }
    int     min_loc_ind()const
    {
        return -l_stencil_loc_2_glob.size();
    }
    int     max_loc_ind()const
    {
        return r_stencil_loc_2_glob.size()-1 +  get_size();
    }
    int     min_own_loc_ind()const
    {
        return 0;
    }
    int     max_own_loc_ind()const
    {
        return get_size()-1;
    }
    bool    check_glob_has_loc_ind(int i_glob)const
    {
        if (check_glob_owned(i_glob)) return true;
        return stencil_glob_2_loc.find(i_glob) != stencil_glob_2_loc.end();
    }
    bool    check_loc_has_loc_ind(int i_loc)const
    {
        return (i_loc >= min_loc_ind())&&(i_loc <= max_loc_ind());
    }
    bool  is_loc_glob_ind_order_preserv()const
    {
        return true;
    }
};

}

#endif
