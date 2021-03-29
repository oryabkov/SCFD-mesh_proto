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

#ifndef __SCFD_MESH_DEVICE_MESH_PARAMS_H__
#define __SCFD_MESH_DEVICE_MESH_PARAMS_H__

namespace scfd
{
namespace mesh
{

struct device_mesh_params
{
    /// Elements data

    bool has_elems_centers_data = true;
    bool has_elems_centers_neighbour_data = true;
    bool has_elems_faces_centers_data = true;
    bool has_elems_vertexes_data = false;
    bool has_elems_neighbours0_data = true;
    bool has_elems_neighbours0_loc_face_i_data = true;
    bool has_elems_faces_boundaries_tags_data = true;
    bool has_elems_volumes_tags_data = true;
    bool has_elems_faces_norms_data = true;
    bool has_elems_faces_areas_data = true;
    bool has_elems_vols_data = true;

    /// Nodes data 

    //TODO

    /// Faces data

    //TODO
};

}  /// namespace mesh
}  /// namespace scfd

#endif
