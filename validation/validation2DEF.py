import cdmath
import FiniteElements2DWithCDMATH
import matplotlib.pyplot as plt
from math import log10

#### 2D FE triangle mesh
nbMeshes=5
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/2DTriangles/'
mesh_name='meshSquareTrianglesFE'
i=0
for filename in ['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5']:
    error_tab[i], mesh_size_tab[i] =FiniteElements2DWithCDMATH.solve(mesh_path+filename)
    error_tab[i]=log10(error_tab[i])
    mesh_size_tab[i] = log10(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.title('Convergence of finite elements for Laplace operator on a 2D triangular mesh')
plt.savefig(mesh_name+".png")

