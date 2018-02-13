Two steps:
- Generate the L09 level mesh by opening the file MS4_RANS_wall.geo in Gmsh [1], clicking on Mesh/Define/2D and saving the .msh file.
- Generate L08 -> L00 level meshes by removing every other line from the L09 level mesh (post-processing on .msh) and writing to a new mesh file.

The resulting values of first element height off the wall should be:

L00 2.584926326078569e-04
L01 1.016293036155189e-04
L02 4.532225455902131e-05
L03 2.143284369136409e-05
L04 1.042578002302168e-05
L05 5.142187573629273e-06
L06 2.553657618145825e-06
L07 1.272499050089030e-06
L08 6.351709303078160e-07
L09 3.173163011862635e-07

References:
[1] http://gmsh.info/
