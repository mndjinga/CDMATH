import FiniteElements2DWithCDMATH

error_tab=[0]*5
mesh_size_tab=[0]*5
i=0
for filename in [triangleMeshSquare_1,triangleMeshSquare_2,triangleMeshSquare_3,triangleMeshSquare_4,triangleMeshSquare_5]:
    error_tab[i] =FiniteElements2DWithCDMATH(filename)
    mesh_size_tab[i]=cdmath.Mesh(filename+"med").getNumberOfNodes()
