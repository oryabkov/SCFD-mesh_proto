
Two levels of mesh abstaction:

BasicMesh
Basic level of abstarction (nodes, element to nodes, nodes to element), elem types, group ids,
some misc staff.

BasicMesh has following issues:
What about elements and nodes enumeration?
Must they be contigious? Must they start from zero?
For now, elements are shrinked to zero, non-contigious case is logic_error.
At the same time nodes id are taken as is from gmsh file.
Not sure about it.

host_mesh<cBasicMesh>
Wrap of BasicMesh abstraction
Most of faces interface (faces to elements, elements to faces, 0 order neighbours).

/// Assertion: ghost_level >= 1
enlarge_stencil(ordinal_type ghost_level)
Verifies that 2nd order (through nodes) neighbours of level ghost_level are read.
It does not work 'backwards' (i.e. call of enlarge_stencil(2) after enlarge_stencil(3)
does nothing).

There is following convention about data and connectivity graphs build.

Define own elements + 2nd order neigbours of level n as ELEMENTS_LEVEL(n).

Say, we have performed enlarge_stencil(n). 
After that we have:
 *For all elements ELEMENTS_LEVEL(n) we have elements to nodes connectivity data. 
 *For all nodes incident to ELEMENTS_LEVEL(n) we have their coords data.
 *For all nodes incident to ELEMENTS_LEVEL(n-1) we have their nodes to elements connectivity data. 
 *For all faces incident to ELEMENTS_LEVEL(n-1) we face their faces to elements connectivity data. 
 *For all elements ELEMENTS_LEVEL(n-1) we have elements to faces connectivity data. 

 Any other request will have undefined behaviour. They can cause error, return correct result, 
 return incorrect result, etc.

 About virtual nodes and faces
Virtual nodes and faces are special mechanism to support periodic meshes (and perhaps, more complicated 
topological structures, like 'torus quaters').
Only nodes and faces can be considered 'virtual' for now (because in there is no point to create
virtual elements, at least for now).
Virtual node or face is a group of physical nodes or faces 'merged' together. During calculations these nodes 
or faces can be dealt as one node or face (for example, used to gather one matrix element for implicit scheme).
Physical nodes or faces are basically geometrical entities. Each physical node has its own distinct coordinates
and each physical face consists of distinct physical nodes. 
Virtual faces are induced from virtual nodes.
Virtual nodes and faces are stored through 'master' node or 'master' face id. Basically all nodes and faces that 
belong to the group of specific virtual node or face are assigned a special id which is virtual node or face id.
There is strict agreement that 'master' node or 'master' face id is id of one of its physical nodes or faces.
This agreement makes some operations easier. For example, when we need to store some graph (for example, graph
of incidence from nodes to elements), we need to create only one special 'virtual' case of this graph - from
phisical entities to virtual entities (from phisical nodes to virtual elements, in example). There is no need
to create separate virtual entity to phisical or virtual entity graph because virtual entities ids are subset 
of physical entities ids.
