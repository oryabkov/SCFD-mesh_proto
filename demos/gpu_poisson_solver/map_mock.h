#ifndef __MAP_MOCK_H__
#define __MAP_MOCK_H__

struct map_mock
{
    //dummy constructor
    map_mock()
    {
    }

    //'read' before construction part:
    int     get_own_rank()const { throw std::runtime_error("map_mock call!"); }
    bool    check_glob_owned(int i)const { throw std::runtime_error("map_mock call!"); }
    int     get_total_size()const { throw std::runtime_error("map_mock call!"); }
    int     get_size()const { throw std::runtime_error("map_mock call!"); }
    //in this map 'i' in own_glob_ind and own_loc_ind simply coincides with local index (not true for general MAP)
    int     own_glob_ind(int i)const { throw std::runtime_error("map_mock call!"); }

    //'construction' part:
    //t_simple_map just ignores this information
    void    add_stencil_element(int i) { throw std::runtime_error("map_mock call!"); }
    void    complete() { throw std::runtime_error("map_mock call!"); }

    //'read' after construction part:
    int     get_rank(int i)const { throw std::runtime_error("map_mock call!"); }
    bool    check_glob_owned(int i, int rank)const { throw std::runtime_error("map_mock call!"); }
    int     loc2glob(int i_loc)const { throw std::runtime_error("map_mock call!"); }
    int     glob2loc(int i_glob)const { throw std::runtime_error("map_mock call!"); }
    int     own_loc_ind(int i)const { throw std::runtime_error("map_mock call!"); }
    int     min_loc_ind()const { throw std::runtime_error("map_mock call!"); }
    int     max_loc_ind()const { throw std::runtime_error("map_mock call!"); }
    int     min_own_loc_ind()const { throw std::runtime_error("map_mock call!"); }
    int     max_own_loc_ind()const { throw std::runtime_error("map_mock call!"); }
    bool    check_glob_has_loc_ind(int i_glob)const { throw std::runtime_error("map_mock call!"); }
    bool    check_loc_has_loc_ind(int i_loc)const { throw std::runtime_error("map_mock call!"); }
    bool    is_loc_glob_ind_order_preserv()const { throw std::runtime_error("map_mock call!"); }
};

#endif
