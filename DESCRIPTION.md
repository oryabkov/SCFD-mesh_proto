
Two levels of mesh abstaction:

BasicMesh
Basic level of abstarction (nodes, element to nodes, nodes to element), elem types, group ids,
some misc staff.

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

