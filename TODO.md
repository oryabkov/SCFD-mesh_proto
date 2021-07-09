
1. Add test for face_2_elem loc_face_index
2. Some interface functions are missing (like node 2 elem inelem_node_ind) + tests.
3. BasicMesh concept
4. host_mesh interface description (is it some concept?)
5. loc_face_i rename to inelem_face_ind (in neib. interfaces)
6. move functors structures into detail folders (device_mesh)
7. in sparse_arr make preallocation for min/max indexes - use it for graphs in mesh classes
8. what about more intellectual hash in sparse_arr? (switch between vector/map or some more intellectual?)
9. is it possible to optimize faces build? mb, preface with writing into msh file?
   parallel version?