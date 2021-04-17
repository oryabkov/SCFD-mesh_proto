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

#ifndef __GMSH_POS_OUTPUT__
#define __GMSH_POS_OUTPUT__

#include <cstdio>
#include <scfd/mesh/host_mesh.h>

template<class T,class Part,class Ord>
void write_out_pos_scalar_file
( 
    const char f_name[], const char v_name[], 
    const scfd::mesh::host_mesh<scfd::mesh::gmsh_mesh_wrap<T,Part,3,Ord>> &mesh, 
    const std::vector<T> &data
)
{
    using vec_t = scfd::static_vec::vec<T,3>;

    FILE *stream;
    //stream = fopen( f_name, "a" );
    stream = fopen( f_name, "w" );

    fprintf( stream, "View");
    fprintf( stream, " '");
    fprintf( stream, v_name);
    fprintf( stream, "' {\n");
    fprintf( stream, "TIME{0};\n");
    fflush(stream);

    for(int i = 0;i < mesh.get_total_elems_num();i++) 
    {
        T   par = data[i];

        if (mesh.get_elem_type(i) == 4) fprintf( stream, "SS(");
        if (mesh.get_elem_type(i) == 5) fprintf( stream, "SH(");
        if (mesh.get_elem_type(i) == 6) fprintf( stream, "SI(");
        if (mesh.get_elem_type(i) == 7) fprintf( stream, "SY(");

        Ord     elem_nodes[mesh.get_elems_max_nodes_num()];
        Ord     nodes_n;
        mesh.get_elem_nodes(i, elem_nodes, &nodes_n);

        for (int j = 0;j < nodes_n;++j) 
        {
            vec_t   vertex = mesh.get_node_coords(elem_nodes[j]);

            fprintf( stream, "%f,%f,%f", vertex[0],vertex[1],vertex[2] );
            if (j != nodes_n-1) fprintf( stream, "," );
        }

        fprintf( stream, "){" );

        for (int j = 0;j < nodes_n;++j) 
        {
            fprintf( stream,"%e", par );
            if (j != nodes_n-1) fprintf( stream, "," );
        }
        fprintf( stream, "};\n");
    }
    fflush(stream);

    fprintf( stream, "};\n");
    fflush(stream);

    fclose(stream);
}

#endif