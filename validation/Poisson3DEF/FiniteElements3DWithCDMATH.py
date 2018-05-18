# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Poisson 3D -\triangle u = f avec conditions aux limites de Dirichlet u=0
# Author      : Michaël Ndjinga, Sédrick Kameni
# Copyright   : CEA Saclay 2017
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u et f discrétisés aux noeuds d'un maillage tétraédrique
#		Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#================================================================================================================================

import cdmath
from math import sin, pi

def solve(filename):
    #Préprocessing optionnel: création du fichier my_mesh.med contenant la géométrie et le maillage du domaine de calcul à partir de commandes python (import salome)
    
    #Chargement du maillage tétraédrique du domaine cubique [0,1]x[0,1]x[0,1], définition des bords
    #==============================================================================================
    my_mesh = cdmath.Mesh(filename+".med")
    if(not my_mesh.isTetrahedral()) :
        raise ValueError("Wrong cell types : mesh is not made of tetrahedra")
    eps=1e-6
    my_mesh.setGroupAtPlan(0.,0,eps,"DirichletBorder")#Bord GAUCHE
    my_mesh.setGroupAtPlan(1.,0,eps,"DirichletBorder")#Bord DROIT
    my_mesh.setGroupAtPlan(0.,1,eps,"DirichletBorder")#Bord BAS
    my_mesh.setGroupAtPlan(1.,1,eps,"DirichletBorder")#Bord HAUT
    my_mesh.setGroupAtPlan(0.,2,eps,"DirichletBorder")#Bord AVANT
    my_mesh.setGroupAtPlan(1.,2,eps,"DirichletBorder")#Bord ARRIERE
    
    nbNodes = my_mesh.getNumberOfNodes()
    nbCells = my_mesh.getNumberOfCells()
    
    print("Mesh building done")
    print("nb of nodes=", nbNodes)
    print("nb of cells=", nbCells)
    
    #Discrétisation du second membre et détermination des noeuds intérieurs
    #======================================================================
    my_RHSfield = cdmath.Field("RHS field", cdmath.NODES, my_mesh, 1)
    B = cdmath.Field("EXA_SOL field", cdmath.NODES, my_mesh, 1)
    nbInteriorNodes = 0
    nbBoundaryNodes = 0
    maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
    interiorNodes=[]
    boundaryNodes=[]
    
    #parcours des noeuds pour discrétisation du second membre et extraction 1) des noeuds intérieur 2) des noeuds frontière 3) du nb max voisins d'un noeud
    for i in range(nbNodes):
        Ni=my_mesh.getNode(i)
        x = Ni.x()
        y = Ni.y()
        z = Ni.z()
    
        my_RHSfield[i]=3*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)#mettre la fonction definie au second membre de l'edp
        B[i]=sin(pi*x)*sin(pi*y)*sin(pi*z)
        if my_mesh.isBorderNode(i): # Détection des noeuds frontière
            boundaryNodes.append(i)
            nbBoundaryNodes=nbBoundaryNodes+1
        else: # Détection des noeuds intérieurs
            interiorNodes.append(i)
            nbInteriorNodes=nbInteriorNodes+1
            maxNbNeighbours= max(1+1*Ni.getNumberOfCells(),maxNbNeighbours) # need a function Ni.getNumberOfNeighbourNodes();
    
    # sauvegarde sur le disque dur du second membre discrétisé dans un fichier paraview
    my_RHSfield.writeVTK("FiniteElements3DRHSField"+filename) 
    B.writeVTK("FiniteElements3DEXSOLField"+filename) 
    
    print("Right hand side discretisation done")
    print("nb of interior nodes=", nbInteriorNodes)
    print("nb of boundary nodes=", nbBoundaryNodes)
    print("Max nb of neighbours=", maxNbNeighbours)
    
    # Construction de la matrice de rigidité et du vecteur second membre du système linéaire
    #=======================================================================================
    Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours) # warning : third argument is maximum number of non zero coefficients per line of the matrix
    RHS=cdmath.Vector(nbInteriorNodes)
    
    # Vecteurs gradient de la fonction de forme associée à chaque noeud d'un hexaèdre
    GradShapeFunc0=cdmath.Vector(3)
    GradShapeFunc1=cdmath.Vector(3)
    GradShapeFunc2=cdmath.Vector(3)
    GradShapeFunc3=cdmath.Vector(3)
    
    #On parcourt les triangles du domaine
    for i in range(nbCells):
    
        Ci=my_mesh.getCell(i)
    
        #Contribution à la matrice de rigidité
        nodeId0=Ci.getNodeId(0)
        nodeId1=Ci.getNodeId(1)
        nodeId2=Ci.getNodeId(2)
        nodeId3=Ci.getNodeId(3)
        N0=my_mesh.getNode(nodeId0)
        N1=my_mesh.getNode(nodeId1)
        N2=my_mesh.getNode(nodeId2)
        N3=my_mesh.getNode(nodeId3)
    
        #Formule des gradients voir EF P1 -> calcul déterminants
        GradShapeFunc0[0]= (N2.y()*N3.z()-N2.z()*N3.y()-N1.y()*N3.z()+N3.y()*N1.z()+N1.y()*N2.z()-N2.y()*N1.z())/6
        GradShapeFunc0[1]=-(N2.x()*N3.z()-N2.z()*N3.x()-N1.x()*N3.z()+N3.x()*N1.z()+N1.x()*N2.z()-N2.x()*N1.z())/6
        GradShapeFunc0[2]=(N2.x()*N3.y()-N2.y()*N3.x()-N1.x()*N3.y()+N3.x()*N1.y()+N1.x()*N2.y()-N2.x()*N1.y())/6
        GradShapeFunc1[0]=- (N2.y()*N3.z()-N2.z()*N3.y()-N0.y()*N3.z()+N3.y()*N0.z()+N0.y()*N2.z()-N2.y()*N0.z())/6
        GradShapeFunc1[1]=(N2.x()*N3.z()-N2.z()*N3.x()-N0.x()*N3.z()+N3.x()*N0.z()+N0.x()*N2.z()-N2.x()*N0.z())/6
        GradShapeFunc1[2]=-(N2.x()*N3.y()-N2.y()*N3.x()-N0.x()*N3.y()+N3.x()*N0.y()+N0.x()*N2.y()-N2.x()*N0.y())/6
        GradShapeFunc2[0]= -(N0.y()*N3.z()-N0.z()*N3.y()-N1.y()*N3.z()+N3.y()*N1.z()+N1.y()*N0.z()-N0.y()*N1.z())/6
        GradShapeFunc2[1]=(N0.x()*N3.z()-N0.z()*N3.x()-N1.x()*N3.z()+N3.x()*N1.z()+N1.x()*N0.z()-N0.x()*N1.z())/6
        GradShapeFunc2[2]= -(N0.x()*N3.y()-N0.y()*N3.x()-N1.x()*N3.y()+N3.x()*N1.y()+N1.x()*N0.y()-N0.x()*N1.y())/6
        GradShapeFunc3[0]=-(N2.y()*N0.z()-N2.z()*N0.y()-N1.y()*N0.z()+N0.y()*N1.z()+N1.y()*N2.z()-N2.y()*N1.z())/6
        GradShapeFunc3[1]=(N2.x()*N0.z()-N2.z()*N0.x()-N1.x()*N0.z()+N0.x()*N1.z()+N1.x()*N2.z()-N2.x()*N1.z())/6
        GradShapeFunc3[2]=-(N2.x()*N0.y()-N2.y()*N0.x()-N1.x()*N0.y()+N0.x()*N1.y()+N1.x()*N2.y()-N2.x()*N1.y())/6
        
        #Création d'un tableau (numéro du noeud, gradient de la fonction de forme)
        GradShapeFuncs={nodeId0 : GradShapeFunc0}
        GradShapeFuncs[nodeId1]=GradShapeFunc1
        GradShapeFuncs[nodeId2]=GradShapeFunc2
        GradShapeFuncs[nodeId3]=GradShapeFunc3
    
        # Remplissage de  la matrice de rigidité et du second membre
        for j in [nodeId0,nodeId1,nodeId2,nodeId3] :
            if boundaryNodes.count(j)==0 : #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
                j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
                #Ajout de la contribution de la cellule triangulaire i au second membre du noeud j 
                RHS[j_int]=Ci.getMeasure()/4*my_RHSfield[j]+RHS[j_int] # intégrale dans le triangle du produit f x fonction de base
                #Contribution de la cellule triangulaire i à la ligne j_int du système linéaire
                for k in [nodeId0,nodeId1,nodeId2,nodeId3] :
                    if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
                        k_int=interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
                        Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
                    #else: si condition limite non nulle au bord, ajouter la contribution du bord au second membre de la cellule j
    
    print("Linear system matrix building done")
    
    # Résolution du système linéaire
    #=================================
    LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"CG","ILU")#,"ILU" Remplacer CG par CHOLESKY pour solveur direct
    SolSyst=LS.solve()
    print "Preconditioner used : ", LS.getNameOfPc()
    print "Number of iterations used : ", LS.getNumberOfIter()
    print("Linear system solved")
    
    # Création du champ résultat
    #===========================
    my_ResultField = cdmath.Field("Result field", cdmath.NODES, my_mesh, 1)
    for j in range(nbInteriorNodes):
        my_ResultField[interiorNodes[j]]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs
    
    for j in range(nbBoundaryNodes):
        my_ResultField[boundaryNodes[j]]=0;#remplissage des valeurs pour les noeuds frontière (condition limite)
    #sauvegarde sur le disque dur du résultat dans un fichier paraview
    my_ResultField.writeVTK("FiniteElements3DResultField"+filename)
    
    print("Numerical solution of 3D poisson equation using finite elements done")
    
    
    #Calcul de l'erreur commise par rapport à la solution exacte
    #===========================================================
    #The following formulas use the fact that the exact solution is equal the right hand side divided by 3*pi*pi
    max_abs_sol_exacte=max(my_RHSfield.max(),-my_RHSfield.min())/(3*pi*pi)
    max_sol_num=my_ResultField.max()
    min_sol_num=my_ResultField.min()
    erreur_abs=0
    for i in range(nbNodes) :
        if erreur_abs < abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i]) :
            erreur_abs = abs(my_RHSfield[i]/(3*pi*pi) - my_ResultField[i])
    
    print("Relative error = max(| exact solution - numerical solution |)/max(| exact solution |) = ",erreur_abs/max_abs_sol_exacte)
    print ("Maximum numerical solution = ", max_sol_num, " Minimum numerical solution = ", min_sol_num)
    
    return erreur_abs/max_abs_sol_exacte, my_mesh.getNumberOfNodes()

if __name__ == """__main__""":
    solve("meshCube")

