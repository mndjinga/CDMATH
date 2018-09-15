# -*- coding: UTF8 -*-

import MEDCoupling as MC
import MEDLoader as ML

## Create a 2D grid mesh given
# @param nb_segs_x the number of segments in the side of the grid
# @param x_start the x coordinate of the first point
# @param mesh_name the name of the mesh
def createMesh(nb_segs_x, x_start=0, mesh_name="Mesh"):
  mesh_dim = 2
  nb_nodes = nb_segs_x+1
  dx = (1.-x_start)/nb_segs_x
  mesh = MC.MEDCouplingIMesh(mesh_name, 2, [nb_nodes, nb_nodes], [x_start, 0.], [dx, dx])
  return mesh

def createCheckerboardMesh(nb_segs_x, mesh_name=""):
  
  if not mesh_name:
    mesh_name = "checkerboard_%ix%i"%(nb_segs_x, nb_segs_x)
  
  # First mesh
  mesh_1 = createMesh(nb_segs_x, 0, mesh_name)
  
  amr = MC.MEDCouplingCartesianAMRMesh(mesh_1)
  
  # Exemple sur une grille 4x4
  ## 1ere diagonale
  #amr.addPatch([(0, 1), (0, 1)], [2, 2])
  #amr.addPatch([(1, 2), (1, 2)], [2, 2])
  #amr.addPatch([(2, 3), (2, 3)], [2, 2])
  #amr.addPatch([(3, 4), (3, 4)], [2, 2])
  
  ## en-dessous de la diagonale, en bas à droite
  #amr.addPatch([(0, 1), (2, 3)], [2, 2])
  #amr.addPatch([(1, 2), (3, 4)], [2, 2])
  
  ## au-dessus de la diagonale, en haut à gauche
  #amr.addPatch([(2, 3), (0, 1)], [2, 2])
  #amr.addPatch([(3, 4), (1, 2)], [2, 2])
  
  # généralise avec une double boucle
  for i in range(0, nb_segs_x, 1):
    for j in range(0, nb_segs_x/2):
      if i%2 == 0:
        amr.addPatch([(i, i+1), (2*j, 2*j+1)], [2, 2])
      else:
        amr.addPatch([(i, i+1), (1+2*j, 2*j+2)], [2, 2])
  
  # Crée un seul maillage avec tous les rafinements
  mesh = amr.buildUnstructured()
  mesh.setName(mesh_name)
  # Crée des polygones pour rendre conforme les mailles
  mesh.conformize2D(1e-10)
  # Merge les noeuds confondus
  arr, areNodesMerged, newNbOfNodes = mesh.mergeNodes(1e-10)

  # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
  mesh_1d = mesh.computeSkin()

  # Identifie les segments de chaque côté pour créer les groupes
  tol = 1e-10
  #tol2 = 0
  #arr_left = mesh_1d.getCellsInBoundingBox([0-tol, tol, -tol, 1+tol], tol2)
  #arr_right = mesh_1d.getCellsInBoundingBox([1-tol, 1+tol, -tol, 1+tol], tol2)
  #arr_bottom = mesh_1d.getCellsInBoundingBox([0-tol, 1+tol, -tol, tol], tol2)
  #arr_top = mesh_1d.getCellsInBoundingBox([0-tol, 1+tol, 1-tol, 1+tol], tol2)
  
  # PB: getCellsInBoundingBox renvoie aussi les segments qui touchent la bounding box
  # => On boucle sur les coordonnées des barycentres

  barycenters = mesh_1d.computeIsoBarycenterOfNodesPerCell()
  ids_left = []
  ids_right = []
  ids_bottom = []
  ids_top = []
  for i, coord in enumerate(barycenters):
    x, y = coord
    if abs(x) < tol:
      ids_left.append(i)
    elif abs(x-1) < tol:
      ids_right.append(i)
    elif abs(y) < tol:
      ids_bottom.append(i)
    elif abs(y-1) < tol:
      ids_top.append(i)

  arr_left = MC.DataArrayInt(ids_left)
  arr_right = MC.DataArrayInt(ids_right)
  arr_bottom = MC.DataArrayInt(ids_bottom)
  arr_top = MC.DataArrayInt(ids_top)
  
  arr_left.setName("left")
  arr_right.setName("right")
  arr_bottom.setName("bottom")
  arr_top.setName("top")

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
  filename = mesh_name+".med"
  meshMEDFile.write(filename,2) # 2 stands for write from scratch
  
  return meshMEDFile

if __name__ == '__main__':
  createCheckerboardMesh(4)
  createCheckerboardMesh(8)
  createCheckerboardMesh(16)
  createCheckerboardMesh(32)
  createCheckerboardMesh(64)
  createCheckerboardMesh(128)
  createCheckerboardMesh(256)