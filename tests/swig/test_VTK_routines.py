import os
import numpy as np
import vtk
from vtk.util import numpy_support as npvtk 
import cdmath

M1 = Mesh(0.0, 1.0, 10, 0., 1., 5)

field1 = Field("test field 1", CELLS, M1, 1)
for j in range(field1.getNumberOfComponents()):
    for i in range(field1.getNumberOfElements()):
        field1[i, j] = i + j

fileNameVTK = "2D_structured_cell_field"
field1.writeVTK(fileNameVTK)

M2 = Mesh("meshSquare.med")
field2 = Field("test field 2", NODES, M2, 1)
for j in range(field2.getNumberOfComponents()):
    for i in range(field2.getNumberOfElements()):
        field2[i, j] = i + j

fileNameVTK = "2D_unstructured_node_field"
field2.writeVTK(fileNameVTK)

M3 = Mesh("meshCube.med")
field3 = Field("test field 3", NODES, M3, 1)
for j in range(field3.getNumberOfComponents()):
    for i in range(field3.getNumberOfElements()):
        field3[i, j] = i + j

fileNameVTK = "3D_unstructured_node_field"
field3.writeVTK(fileNameVTK)

M4 = Mesh("meshSphere.med")
field4 = Field("test field 4", NODES, M4, 1)
for j in range(field4.getNumberOfComponents()):
    for i in range(field4.getNumberOfElements()):
        field4[i, j] = i + j

fileNameVTK = "Sphere_unstructured_node_field"
field4.writeVTK(fileNameVTK)

M5 = Mesh(0.0, 1.0, 4, 0.0, 1.0, 4, 0.0, 1.0, 4)
field5 = Field("test field 5", CELLS, M5, 1)
for j in range(field5.getNumberOfComponents()):
    for i in range(field5.getNumberOfElements()):
        field5[i, j] = i + j

fileNameVTK = "Cube_structured_cell_field"
field5.writeVTK(fileNameVTK)


