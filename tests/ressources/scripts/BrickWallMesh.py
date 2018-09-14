# -*-coding:utf-8 -*

#### Maillage d'un rectangle par une grille décalée (mur de briques)
### input : xmin, xmax, nx, ymin, ymax, ny
### output : squareWithBrickWall.vtu squareWithBrickWall.med

import MEDCoupling as mc
import math
import MEDLoader as ML

xmin=0.
xmax=1.
ymin=0.
ymax=1.
nx=15
ny=15

dx = (xmax-xmin)/nx
dy=(ymax-ymin)/ny

print "nx=",nx,"ny=",ny, "ncells=",nx*ny

# Building the initial rectangular cell, centered at 0,0
d = mc.DataArrayDouble(4,2)
d[0,0] = -dx/2
d[0,1] =  dy/2
d[1,0] =  dx/2
d[1,1] =  dy/2
d[2,0] =  dx/2
d[2,1] = -dy/2
d[3,0] = -dx/2
d[3,1] = -dy/2

d.setInfoOnComponents(["X [m]","Y [m]"])

print "Uniform array ?", d.magnitude().isUniform(0.5*math.sqrt(dx*dx+dy*dy),1e-10)

# translation of the first cell
translationToPerform = [[(0.5*(1+j%2)+i)*dx,(0.5+j)*dy] for i in range(nx) for j in range(ny)]
    
ds = len(translationToPerform)*[None]
for pos,t in enumerate(translationToPerform):
                 ds[pos] = d[:]         # Perform a deep copy of d and place it at position 'pos' in ds
                 ds[pos] += t             # Adding a vector to a set of coordinates does a translation
                 pass
# Identifying duplicate tuples
d2 = mc.DataArrayDouble.Aggregate(ds)
oldNbOfTuples = d2.getNumberOfTuples()
c,cI = d2.findCommonTuples(1e-10)
tmp = c[cI[0]:cI[0+1]]
print tmp
a = cI.deltaShiftIndex()
b = a - 1
myNewNbOfTuples = oldNbOfTuples - sum(b.getValues())
o2n, newNbOfTuples = mc.DataArrayInt.ConvertIndexArrayToO2N(oldNbOfTuples,c,cI)
print "Have I got the right number of tuples?"
print "myNewNbOfTuples = %d, newNbOfTuples = %d" % (myNewNbOfTuples, newNbOfTuples)
assert(myNewNbOfTuples == newNbOfTuples)
print "Old number of tuple was ", oldNbOfTuples

# Extracting the unique set of tuples
d3 = d2.renumberAndReduce(o2n, newNbOfTuples)
n2o = o2n.invertArrayO2N2N2O(newNbOfTuples)
d3_bis = d2[n2o]
print "Are d3 and d3_bis equal ? %s" % (str(d3.isEqual(d3_bis, 1e-12)))
# Now build an unstructured mesh representing the final pattern
mesh = mc.MEDCouplingUMesh("squareWithBrickWall",2)
mesh.setCoords(d3)
print "Mesh dimension is", mesh.getMeshDimension()
print "Spatial dimension is", mesh.getCoords().getNumberOfComponents()
mesh.allocateCells(nx*ny)
for i in xrange(nx*ny):
        cell_connec = o2n[4*i:4*(i+1)]
        mesh.insertNextCell(mc.NORM_POLYGON, cell_connec.getValues())
        pass

# Crée des polygones pour rendre conforme les mailles
mesh.conformize2D(1e-10)
# Check that everything is coherent (will throw if not)
mesh.checkConsistencyLight()

# Crée les éléments 1D pour pouvoir imposer les conditions aux limites
mesh_1d = mesh.computeSkin()

# Identifie les segments de chaque côté pour créer les groupes
tol = 1e-10

# PB: getCellsInBoundingBox renvoie aussi les segments qui touchent la bounding box
# => On boucle sur les coordonnées des barycentres

barycenters = mesh_1d.computeIsoBarycenterOfNodesPerCell()
ids_left = []
ids_right = []
ids_bottom = []
ids_top = []
for i, coord in enumerate(barycenters):
    x, y = coord
    if abs(y-ymin) < tol :
      ids_bottom.append(i)
    elif abs(y-ymax) < tol :
      ids_top.append(i)
    elif abs(x-xmax) < tol or abs(x-xmax-dx/4) < tol or abs(x-xmax-dx/2) < tol:
      ids_right.append(i)
    elif abs(x-xmin) < tol or abs(x-xmin-dx/4) < tol or abs(x-xmin-dx/2) < tol:
      ids_left.append(i)
    else:
        raise ValueError("Pb with boundary construction : barycenter does not belong to any border group")
	
arr_left = mc.DataArrayInt(ids_left)
arr_right = mc.DataArrayInt(ids_right)
arr_bottom = mc.DataArrayInt(ids_bottom)
arr_top = mc.DataArrayInt(ids_top)

arr_left.setName("Left")
arr_right.setName("Right")
arr_bottom.setName("Bottom")
arr_top.setName("Top")

# Trie les cellules par type conformément à la convention MED fichier
o2n = mesh.sortCellsInMEDFileFrmt()
meshMEDFile=ML.MEDFileUMesh.New()
# Ecrit le maillage 2D
meshMEDFile.setMeshAtLevel(0,mesh)
# Ecrit le maillage 1D
meshMEDFile.setMeshAtLevel(-1,mesh_1d)
# Ecrit les groupes
meshMEDFile.addGroup(-1, arr_left)
meshMEDFile.addGroup(-1, arr_right)
meshMEDFile.addGroup(-1, arr_bottom)
meshMEDFile.addGroup(-1, arr_top)

filename = "squareWithBrickWall"+".med"
# Write the result into a VTU file that can be read with ParaView
mesh.writeVTK("squareWithBrickWall.vtu")
# Write the result into a MED file that can be read with Salomé
meshMEDFile.write(filename,2) # 2 stands for write from scratch
