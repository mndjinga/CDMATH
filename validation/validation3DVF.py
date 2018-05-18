import cdmath
import FiniteVolumes3DWithCDMATH
import matplotlib.pyplot as plt
from math import log10

#### 3D FV tetrahedra mesh
nbMeshes=6
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/3DTetrahedra/'
mesh_name='meshCubeTetrahedraFV'
i=0
for filename in ['meshCubeTetrahedra_1','meshCubeTetrahedra_2','meshCubeTetrahedra_3','meshCubeTetrahedra_4','meshCubeTetrahedra_5','meshCubeTetrahedra_6']:
    error_tab[i], mesh_size_tab[i] =FiniteVolumes3DWithCDMATH.solve_file(mesh_path+filename)
    error_tab[i]=log10(error_tab[i])
    mesh_size_tab[i] = log10(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.title('Convergence of finite volumes for Laplace operator on a 3D tetrahedral mesh')
plt.savefig(mesh_name+".png")

### 3D FV rectangular mesh
nbMeshes=4
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_name='meshCubeRectangles'
i=0
for nx in [11,51,101]:
    my_mesh=cdmath.Mesh(0,1,nx,0,1,nx,0,1,nx)
    error_tab[i], mesh_size_tab[i] =FiniteVolumes3DWithCDMATH.solve(my_mesh,str(nx)+'x'+str(nx)+'x'+str(nx))
    error_tab[i]=log10(error_tab[i])
    mesh_size_tab[i] = log10(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.title('Convergence of finite volumes for Laplace operator on a 3D rectangular mesh')
plt.savefig(mesh_name+".png")

