# -*-coding:utf-8 -*

#### Maillage d'un cercle par une grille de carrés
### input : xcenter, ycenter, Radius, n
### output : CircleWithSquares.vtu

import MEDCoupling as mc
import math
import MEDLoader as ML

xcenter=0.
ycenter=0.
Radius=1.
n=17

xmin=-Radius
xmax=Radius
ymin=-Radius
ymax=Radius

dx = (xmax-xmin)/n
dy=(ymax-ymin)/n

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
translationToPerform = []
for i in range(n) :
    for j in range(n):
        if (xcenter-xmin-(0.5+i)*dx)**2+(ycenter-ymin-(0.5+j)*dy)**2<Radius*Radius :
            translationToPerform.append([xmin+(0.5+i)*dx,ymin+(0.5+j)*dy] )

ncells= len(translationToPerform) 
print "n=",n,"ncells=",ncells
  
ds = ncells*[None]
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
mesh = mc.MEDCouplingUMesh("diskWithSquares",2)
mesh.setCoords(d3)
print "Mesh dimension is", mesh.getMeshDimension()
print "Spatial dimension is", mesh.getCoords().getNumberOfComponents()
mesh.allocateCells(ncells)
for i in xrange(ncells):
        cell_connec = o2n[4*i:4*(i+1)]
        mesh.insertNextCell(mc.NORM_POLYGON, cell_connec.getValues())
        pass

# Check that everything is coherent (will throw if not)
mesh.checkConsistencyLight()

# Crée les éléments 1D pour pouvoir imposer les conditions aux limites
mesh_1d = mesh.computeSkin()

# Trie les cellules par type conformément à la convention MED fichier
o2n = mesh.sortCellsInMEDFileFrmt()
meshMEDFile=ML.MEDFileUMesh.New()
# Ecrit le maillage 2D
meshMEDFile.setMeshAtLevel(0,mesh)
# Ecrit le maillage 1D
meshMEDFile.setMeshAtLevel(-1,mesh_1d)
# Ecrit les groupes
arr_circle = mc.DataArrayInt(range(mesh_1d.getNumberOfCells()))
arr_circle.setName("Circle")
meshMEDFile.addGroup(-1, arr_circle)

filename = "diskWithSquares"+".med"
# Write the result into a VTU file that can be read with ParaView
mesh.writeVTK("diskWithSquares.vtu")
# Write the result into a MED file that can be read with Salomé
meshMEDFile.write(filename,2) # 2 stands for write from scratch
