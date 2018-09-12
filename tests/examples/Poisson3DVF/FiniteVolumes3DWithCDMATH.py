# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Poisson -\triangle u = f sur unu cube avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2016
# Description : Utilisation de la méthode des volumes finis avec champs u et f discrétisés aux cellules d'un maillage quelconque
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant CDMATH
#================================================================================================================================

import cdmath
from math import sin, pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
import PV_routines
import VTK_routines

import sys

if len(sys.argv) >1 :#non rectangular mesh
    my_mesh = cdmath.Mesh(sys.argv[1])
else :   #rectangular mesh
# Maillage du domaine cubique [0,1]x[0,1]x[0,1], définition des bords
#====================================================================================
    xmin=0
    xmax=1
    ymin=0
    ymax=1
    zmin=0 
    zmax=1
    
    nx=21
    ny=21
    nz=21
    
    my_mesh = cdmath.Mesh(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz)#Ajouter ,zmin,zmax,nz si calcul 3D
eps=1e-6
my_mesh.setGroupAtPlan(0,0,eps,"DirichletBorder")#Bord GAUCHE
my_mesh.setGroupAtPlan(1,0,eps,"DirichletBorder")#Bord DROIT
my_mesh.setGroupAtPlan(0,1,eps,"DirichletBorder")#Bord BAS
my_mesh.setGroupAtPlan(1,1,eps,"DirichletBorder")#Bord HAUT
my_mesh.setGroupAtPlan(0,2,eps,"DirichletBorder")#Bord AVANT 
my_mesh.setGroupAtPlan(1,2,eps,"DirichletBorder")#Bord ARRIERE 

nbCells = my_mesh.getNumberOfCells()

print("Mesh building done")
print("nb of cells ", nbCells)

#Discrétisation du second membre et extraction du nb max de voisins d'une cellule
#================================================================================
my_RHSfield = cdmath.Field("RHS field", cdmath.CELLS, my_mesh, 1)
maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
#parcours des cellules pour discrétisation du second membre et extraction du nb max de voisins d'une cellule
for i in range(nbCells): 
	Ci = my_mesh.getCell(i)
	x = Ci.x()
	y = Ci.y()
	z = Ci.z() 
	my_RHSfield[i]=3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)#mettre la fonction definie au second membre de l edp
	# compute maximum number of neighbours
	maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)

# sauvegarde sur le disque dur du second membre discrétisé dans un fichier paraview
my_RHSfield.writeVTK("FiniteVolumes3DRHSField")

print("Right hand side discretisation done")
print("Max nb of neighbours=", maxNbNeighbours)

# Construction de la matrice et du vecteur second membre du système linéaire
#===========================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbCells)
#Parcours des cellules du domaine
for i in range(nbCells):
	RHS[i]=my_RHSfield[i] #la valeur moyenne du second membre f dans la cellule i
	Ci=my_mesh.getCell(i)
	for j in range(Ci.getNumberOfFaces()):# parcours des faces voisinnes
		Fj=my_mesh.getFace(Ci.getFaceId(j))
		if not Fj.isBorder():
			k=Fj.getCellId(0)
			if k==i :
				k=Fj.getCellId(1)
			Ck=my_mesh.getCell(k)
			distance=Ci.getBarryCenter().distance(Ck.getBarryCenter())
			coeff=Fj.getMeasure()/Ci.getMeasure()/distance
			Rigidite.setValue(i,k,-coeff) # terme extradiagonal
		else:
			coeff=Fj.getMeasure()/Ci.getMeasure()/Ci.getBarryCenter().distance(Fj.getBarryCenter())
		Rigidite.addValue(i,i,coeff) # terme diagonal

print("Linear system matrix building done")

# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"GMRES","ILU")
SolSyst=LS.solve()

print "Preconditioner used : ", LS.getNameOfPc()
print "Number of iterations used : ", LS.getNumberOfIter()
print "Final residual : ", LS.getResidu()
print("Linear system solved")

# Création du champ résultat
#===========================
my_ResultField = cdmath.Field("ResultField", cdmath.CELLS, my_mesh, 1)
for i in range(nbCells):
    my_ResultField[i]=SolSyst[i];
#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_ResultField.writeVTK("FiniteVolumes3DResultField")

#Postprocessing 
#==============
# save 2D picture
PV_routines.Save_PV_data_to_picture_file("FiniteVolumes3DResultField"+'_0.vtu',"ResultField",'CELLS',"FiniteVolumes3DResultField")

# extract diagonal values
resolution=100
curv_abs=np.linspace(0,sqrt(2),resolution+1)
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,1,0],[1,0,0], resolution)
plt.legend()
plt.xlabel('Position on diagonal line')
plt.ylabel('Value on diagonal line')
if len(sys.argv) >1 :
    plt.title('Plot over diagonal line for finite Volumes \n for Laplace operator on a 3D mesh '+my_mesh.getName())
    plt.plot(curv_abs, diag_data, label= str(nbCells)+ ' cells mesh')
    plt.savefig("FiniteVolumes3DResultField_"+str(nbCells)+ '_cells'+"_PlotOverDiagonalLine.png")
else :   
    plt.title('Plot over diagonal line for finite Volumes \n for Laplace operator on a 3D rectangular grid')
    plt.plot(curv_abs, diag_data, label= str(nx) +'x'+str(ny)+ ' cells mesh')
    plt.savefig("FiniteVolumes3DResultField_"+str(nx) +'x'+str(ny)+ '_cells'+"_PlotOverDiagonalLine.png")

print("Numerical solution of 3D poisson equation using finite volumes done")

#Calcul de l'erreur commise par rapport à la solution exacte
#===========================================================
#The following formulas use the fact that the exact solution is equal to the right hand side divided by 3*pi*pi
max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(3*pi*pi)
max_sol_num=my_ResultField.max()
min_sol_num=my_ResultField.min()
erreur_abs=0
for i in range(nbCells) :
	if erreur_abs < abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i]) :
		erreur_abs = abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i])

print("Absolute error = max(| exact solution - numerical solution |) = ",erreur_abs )
print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
print ("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)

assert erreur_abs/max_abs_sol_exacte <1.