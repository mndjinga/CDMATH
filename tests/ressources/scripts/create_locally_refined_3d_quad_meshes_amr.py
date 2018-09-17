# -*- coding: UTF8 -*-

import MEDCoupling as MC
import MEDLoader as ML

## Create a 2D grid mesh given
# @param nb_segs_x the number of segments in the side of the grid
# @param x_start the x coordinate of the first point
# @param mesh_name the name of the mesh
def createMesh(nb_segs_x, x_start=0, mesh_name="Mesh"):
  mesh_dim = 3
  nb_nodes = nb_segs_x+1
  dx = (1.-x_start)/nb_segs_x
  mesh = MC.MEDCouplingIMesh(mesh_name, mesh_dim, [nb_nodes]*mesh_dim, [x_start, 0., 0.], [dx]*mesh_dim)
  return mesh

def createLocallyRefinedMesh(nb_segs_x, mesh_name):
  
  # First mesh
  mesh_1 = createMesh(nb_segs_x, 0, mesh_name)
  
  amr = MC.MEDCouplingCartesianAMRMesh(mesh_1)
  
  # 1er raffinement
  amr.addPatch([(nb_segs_x/2,nb_segs_x),(0,nb_segs_x/2),(0,nb_segs_x/2)],[2,2,2])
  # 2eme raffinement
  amr[0].addPatch([(nb_segs_x/2,nb_segs_x),(0,nb_segs_x/2),(0,nb_segs_x/2)],[2,2,2])
  
  # Crée un seul maillage avec tous les rafinements
  mesh = amr.buildUnstructured()
  mesh.setName(mesh_name)
  # Merge les noeuds confondus (à faire avant le conformize2D)
  arr, areNodesMerged, newNbOfNodes = mesh.mergeNodes(1e-10)
  # Crée des polyèdres pour rendre conforme les mailles
  mesh.convertAllToPoly()
  mesh.conformize3D(1e-10)
  mesh.unPolyze()

  # Crée les éléments 1D pour pouvoir imposer les conditions aux limites
  mesh_2d = mesh.computeSkin()

  # Trie les cellules par type conformément à la convention MED fichier
  o2n = mesh.sortCellsInMEDFileFrmt()
  o2n = mesh_2d.sortCellsInMEDFileFrmt()

  # Identifie les faces de chaque côté pour créer les groupes
  tol = 1e-10

  barycenters = mesh_2d.computeIsoBarycenterOfNodesPerCell()
  ids_left = []
  ids_right = []
  ids_bottom = []
  ids_top = []
  ids_front = []
  ids_back = []
  for i, coord in enumerate(barycenters):
    x, y, z = coord
    if abs(x) < tol:
      ids_left.append(i)
    elif abs(x-1) < tol:
      ids_right.append(i)
    elif abs(y) < tol:
      ids_bottom.append(i)
    elif abs(y-1) < tol:
      ids_top.append(i)
    elif abs(z) < tol:
      ids_back.append(i)
    elif abs(z-1) < tol:
      ids_front.append(i)

  arr_left = MC.DataArrayInt(ids_left)
  arr_right = MC.DataArrayInt(ids_right)
  arr_bottom = MC.DataArrayInt(ids_bottom)
  arr_top = MC.DataArrayInt(ids_top)
  arr_back = MC.DataArrayInt(ids_back)
  arr_front = MC.DataArrayInt(ids_front)

  arr_left.setName("Left")
  arr_right.setName("Right")
  arr_bottom.setName("Bottom")
  arr_top.setName("Top")
  arr_back.setName("Back")
  arr_front.setName("Front")

  meshMEDFile=ML.MEDFileUMesh.New()
  # Ecrit le maillage 3D
  meshMEDFile.setMeshAtLevel(0,mesh)
  # Ecrit le maillage 2D
  meshMEDFile.setMeshAtLevel(-1,mesh_2d)
  # Ecrit les groupes
  meshMEDFile.addGroup(-1, arr_left)
  meshMEDFile.addGroup(-1, arr_right)
  meshMEDFile.addGroup(-1, arr_bottom)
  meshMEDFile.addGroup(-1, arr_top)
  meshMEDFile.addGroup(-1, arr_back)
  meshMEDFile.addGroup(-1, arr_front)
  filename = mesh_name+".med"
  meshMEDFile.write(filename,2) # 2 stands for write from scratch
  
  return meshMEDFile

if __name__ == '__main__':
  createLocallyRefinedMesh(4, "cubeWithLocRefCubes_1")
  createLocallyRefinedMesh(8, "cubeWithLocRefCubes_2")
  createLocallyRefinedMesh(16, "cubeWithLocRefCubes_3")
  createLocallyRefinedMesh(32, "cubeWithLocRefCubes_4")
