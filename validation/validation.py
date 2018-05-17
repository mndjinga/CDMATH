import cdmath
import FiniteElements2DWithCDMATH

error_tab=[0]*5
mesh_size_tab=[0]*5
i=0
mesh_path='../validation/2DTriangles/'

for filename in ['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5']:
    error_tab[i] =FiniteElements2DWithCDMATH.solve(mesh_path+filename)
    mesh_size_tab[i]=cdmath.Mesh(mesh_path+filename+".med").getNumberOfNodes()
