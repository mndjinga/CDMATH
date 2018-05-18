import cdmath
import FiniteVolumes2DWithCDMATH
import matplotlib.pyplot as plt
import numpy as np
from math import log10

##### 2D FV triangle mesh
nbMeshes=5
error_tab=[0]*nbMeshes
mesh_size_tab=[0]*nbMeshes
mesh_path='../validation/2DTriangles/'
mesh_name='meshSquareTrianglesFV'
i=0
# Computation of the numerical error
for filename in ['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5']:
    error_tab[i], mesh_size_tab[i] =FiniteVolumes2DWithCDMATH.solve_file(mesh_path+filename)
    error_tab[i]=log10(error_tab[i])
    mesh_size_tab[i] = log10(mesh_size_tab[i])
    i=i+1
    
# Least square linear regression
# Find the best a,b such that f(x)=ax+b best approximates the convergence curve
# The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
a1=np.dot(mesh_size_tab,mesh_size_tab)
a2=np.sum(mesh_size_tab)
a3=nbMeshes
b1=np.dot(error_tab,mesh_size_tab)   
b2=np.sum(error_tab)

det=a1*a3-a2*a2
a=( a3*b1-a2*b2)/det
b=(-a2*b1+a1*b2)/det

print "FV on 2D triangle mesh :scheme order is ", -a

# Plot of figures
plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+str(a))
plt.legend()
plt.xlabel('log(number of nodes)')
plt.ylabel('log(error)')
plt.title('Convergence of finite volumes for Laplace operator on a 2D triangular mesh')
plt.savefig(mesh_name+".png")
