// Copyright © 2016,2017 Ryabkov Oleg Igorevich, Evstigneev Nikolay Mikhaylovitch

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
#include <mesh/t_cpu_mesh_tml.h>

template<class T>
void write_out_pos_scalar_file( const char f_name[], const char v_name[], const t_cpu_mesh_tml<T> &mesh, T *data)
{
    FILE *stream;
    //stream = fopen( f_name, "a" );
    stream = fopen( f_name, "w" );

    fprintf( stream, "View");
    fprintf( stream, " '");
    fprintf( stream, v_name);
    fprintf( stream, "' {\n");
    fprintf( stream, "TIME{0};\n");
    fflush(stream);

    for(int i = 0;i < mesh.cv.size();i++) {
        T   par = data[i];

        if (mesh.cv[i].elem_type == 4) fprintf( stream, "SS(");
        if (mesh.cv[i].elem_type == 5) fprintf( stream, "SH(");
        if (mesh.cv[i].elem_type == 6) fprintf( stream, "SI(");
        if (mesh.cv[i].elem_type == 7) fprintf( stream, "SY(");

        for (int j = 0;j < mesh.cv[i].vert_n;++j) {
            fprintf( stream, "%f,   %f, %f",
                mesh.cv[i].vertexes[j][0],mesh.cv[i].vertexes[j][1],mesh.cv[i].vertexes[j][2] );
                    if (j != mesh.cv[i].vert_n-1) fprintf( stream, ",   " );
        }

        fprintf( stream, "){" );

                for (int j = 0;j < mesh.cv[i].vert_n;++j) {
            fprintf( stream,"%e", par );
                        if (j != mesh.cv[i].vert_n-1) fprintf( stream, ",   " );
        }
        fprintf( stream, "};\n");
    }
    fflush(stream);

    fprintf( stream, "};\n");
    fflush(stream);

    fclose(stream);
}

#endif