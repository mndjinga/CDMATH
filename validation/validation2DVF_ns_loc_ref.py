import cdmath
import FiniteVolumes2DWithCDMATH
import matplotlib.pyplot as plt
from math import log10

##### 2D FV refined squares mesh
nbMeshes=7
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/2DLocRefinedSquares/'
mesh_name='meshSquareRefSquaresFV'
i=0
for filename in ['meshLocRefSquares_1','meshLocRefSquares_2','meshLocRefSquares_3','meshLocRefSquares_4','meshLocRefSquares_5','meshLocRefSquares_6','meshLocRefSquares_7']:
    error_tab[i], mesh_size_tab[i] =FiniteVolumes2DWithCDMATH.solve_file(mesh_path+filename)
    error_tab[i]=log10(error_tab[i])
    mesh_size_tab[i] = log10(mesh_size_tab[i])
    i=i+1
    
plt.plot(mesh_size_tab, error_tab)
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.title('Convergence of finite volumes for Laplace operator on a 2D refined squares mesh')
plt.savefig(mesh_name+".png")
