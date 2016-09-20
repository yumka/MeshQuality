__author__ = 'ymelniko'

from mesh_quality_routines import get_mesh_quality

# Quality assessment of unstructured tetrahedral mesh
# filename: ,vtu filename
# quality_metric: metric name, should match name in Paraview/vtk
# qmin: lower quality bound
# qmax: upper quality bound
# n_int: number of subranges for quality output
# Examples of metrics names:
#'Condition' = condition number of the matrix used for transforming element into ideal tetrahedra
#'MinAngle' = minumum dihedral angle of element

#filename = 'box_0.vtu'

filename='mesh1_3d.e'

quality_metric = 'Condition'
get_mesh_quality(filename, quality_metric,qmin=1,qmax=5, n_int=5)

