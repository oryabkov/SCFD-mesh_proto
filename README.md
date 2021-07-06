
  #ABOUT 
This is prototype version for (presumably unstructured) mesh framework in SCFD.
Previuos version (SimpleCFD/samples_n_protos/unstructured_mesh_framework) is 
taken as basis for this new one. Also, most recent versions of 'communication' 
(SimpleCFD/samples_n_protos/communication) and 'for_each' 
(SimpleCFD/samples_n_protos/for_each) are supposed to be used.
Main goals are the following:
* Getting rid of internal mesh format (both file fromat and CPU maintain structures) 
as something excessive. There will be no converters. It is supposed that solvers will
read externals mesh files at once.
* Separation between device mesh and host mesh should be increased. Just some basic 
features should be requiered from mesh class (basically nodes + elements, perhaps tags 
also will be needed). Plan is to use some sort of mesh traits class.
* Periodic connections (boundaries) must be supported.
* High-order elements must be supported.
* Faces fileds astraction must be supported.
* Parallel reading for MPI version is preferrable but limited by current mesh frameworks
* Calculation of all geometric properties (normals, surfaces, volumes, etc) is now will 
be perfomed on device side (not on some sort of conversion preface and not on pure cpu 
preface).
* Ideally pure 2d case must be supported (not sure if it is possible for now).
* Most properties must be done in run-time level when it is not degrading performance
(previously most of the properties like max_faces were done compile-time because of old arrays).
* Maybe, structured meshes?

Initially we are oriented on gmhs usage (more familiar and as we know capabilities are best).
Only substantial drawback is parallel reading which is not directly supported. Next closest 
candidate is VTK, but for now many drawbacks are faced (for example, absence of tags spit during 
gmsh save, Salome does not support this format).

For now, this is only prototype and no direct relation with SCFD is supported (for example, not copied in 'all' repo of SCFD). When main features are outlined new module (perhaps, called 'mesh')
will be created and this repo will be depricated.

  #Abstractions
We need some low-level host mesh wrap (elements,nodes,tags). Because of tags organization in gmsh
we cannot directly use GModel + some traits, distinct class is needed. Also types conversion will be performed there.

We also need some higher level host mesh abstraction that will create more internal information
(faces, neigbours, nodes to elements graph, etc...).