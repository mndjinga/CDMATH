import cdmath
import FiniteElements2DWithCDMATH
import FiniteElements3DWithCDMATH
import matplotlib.pyplot as plt
from math import log

### 2D triangle mesh
nbMeshes=5
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/2DTriangles/'
mesh_name='triangleMeshSquare'
i=0
for filename in ['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5']:
    error_tab[i], mesh_size_tab[i] =FiniteElements2DWithCDMATH.solve(mesh_path+filename)
    error_tab[i]=log(error_tab[i])
    mesh_size_tab[i] = log(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.savefig(mesh_name+".png")

### 3D tetrahedra mesh
nbMeshes=6
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/3DTetrahedra/'
mesh_name='meshCubeTetrahedra'
i=0
for filename in ['meshCubeTetrahedra_1','meshCubeTetrahedra_2','meshCubeTetrahedra_3','meshCubeTetrahedra_4','meshCubeTetrahedra_5','meshCubeTetrahedra6']:
    error_tab[i], mesh_size_tab[i] =FiniteElements3DWithCDMATH.solve(mesh_path+filename)
    error_tab[i]=log(error_tab[i])
    mesh_size_tab[i] = log(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.savefig(mesh_name+".png")

