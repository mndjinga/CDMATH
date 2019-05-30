#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF du système des ondes 2D sans terme source
#                \partial_t p + c^2 \div q = 0
#                \partial_t q +    \grad p = 0
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Propagation d'une onde de choc
#               Utilisation du schéma upwind explicite ou implicite sur un maillage général
#               Initialisation par une surpression sphérique
#               Conditions aux limites périodiques
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sin, cos, pi, sqrt
import cdmath
import PV_routines
import VTK_routines
import sys

p0=155e5#reference pressure in a pressurised nuclear vessel
c0=700.#reference sound speed for water at 155 bars, 600K
rho0=p0/c0*c0#reference density
precision=1e-5

def initial_conditions_shock(my_mesh):
    print "Initial data : Spherical wave"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    rayon = 0.15
    xcentre = 0.5
    ycentre = 0.5
    zcentre = 0.5


    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()

        velocity_field[i,0] = 0
        velocity_field[i,1] = 0
        velocity_field[i,2] = 0

        valX = (x - xcentre) * (x - xcentre)
        valY = (y - ycentre) * (y - ycentre)

        if(dim==2):
            val =  sqrt(valX + valY)
        if(dim==3):
            z = my_mesh.getCell(i).z()
            valZ = (z - zcentre) * (z - zcentre)
            val =  sqrt(valX + valY + valZ)

        if val < rayon:
            pressure_field[i] = p0
            pass
        else:
            pressure_field[i] = p0/2
            pass
        pass

    return pressure_field, velocity_field

def jacobianMatrices(normal):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)
    absA=cdmath.Matrix(dim+1,dim+1)

    absA[0,0]=c0
    for i in range(dim):
        A[i+1,0]=normal[i]
        A[0,i+1]=c0*c0*normal[i]
        for j in range(dim):
            absA[i+1,j+1]=c0*normal[i]*normal[j]
    
    return A, absA
    
def Flux(U, normal):
    dim=normal.size()

    result=cdmath.Vector(dim+1)
    for i in range(dim):
        result[0]  +=normal[i]*U[i+1]
        result[i+1] =normal[i]*U[0]
        
    result[0]=c0*c0*result[0]
    
    return result
    
def computeFluxes(U, SumFluxes):
    my_mesh =U.getMesh();
    nbCells = my_mesh.getNumberOfCells();
    dim=my_mesh.getMeshDimension();
    nbComp=U.getNumberOfComponents();
    Fcourant=cdmath.Vector(nbComp)
    Fautre=cdmath.Vector(nbComp)
    Ucourant=cdmath.Vector(nbComp)
    Uautre=cdmath.Vector(nbComp)
    normal=cdmath.Vector(dim)
    sumFluxCourant=cdmath.Vector(nbComp)

    for j in range(nbCells):#On parcourt les cellules         
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();
        for i in range(nbComp) :
            Ucourant[i]=U[j,i];
            sumFluxCourant[i]=0;

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

            cellAutre = -1;
            if ( not Fk.isBorder()) :
                # hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if (Fk.getCellsId()[0] == j) :
                    # hypothese verifiée 
                    cellAutre = Fk.getCellsId()[1];
                elif(Fk.getCellsId()[1] == j) :
                    # hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                else :
                    raise ValueError("computeFluxes: problem with mesh, unknown cel number");
                
                for i in range(nbComp):
                    Uautre[i]=U[cellAutre,i]
            else :
                if(Fk.getGroupName() != "Neumann"):#Wall boundary condition unless Wall/Neumann specified explicitly
                    for i in range(nbComp):
                        Uautre[i]=Ucourant[i]
                    qn=0# normal momentum
                    for i in range(dim):
                        qn+=Ucourant[i+1]*normal[i]
                    for i in range(dim):
                        Uautre[i+1]-=2*qn*normal[i]
                elif(Fk.getGroupName() == "Neumann"):
                    for i in range(nbComp):
                        Uautre[i]=Ucourant[i]
                else:
                    print Fk.getGroupName()
                    raise ValueError("computeFluxes: Unknown boundary condition name");
            
            Fcourant=Flux(Ucourant,normal);
            Fautre  =Flux(Uautre,  normal);

            A, absA=jacobianMatrices( normal);
            sumFluxCourant = sumFluxCourant + (Fcourant+Fautre +absA*(Ucourant-Uautre))*Fk.getMeasure()*0.5
                
        #On divise par le volume de la cellule la contribution des flux au snd membre
        for i in range(nbComp):
            SumFluxes[j,i]=sumFluxCourant[i]/Cj.getMeasure()


def WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    meshName=my_mesh.getName()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    SumFluxes = cdmath.Field("Fluxes", cdmath.CELLS, my_mesh, dim+1)

    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_shock(my_mesh)

    U              = cdmath.Field("Conservative vector", cdmath.CELLS, my_mesh, dim+1)
    for i in range(nbCells):
        U[i,0] =       pressure_field[i]
        U[i,1] =  rho0*velocity_field[i,0]
        U[i,2] =  rho0*velocity_field[i,1]

    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity");

    dx_min=my_mesh.minRatioVolSurf()

    dt = cfl * dx_min / c0
    
    print("Starting computation of the linear wave system with an explicit UPWIND scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        computeFluxes(U,SumFluxes);

        SumFluxes*=dt;
        maxVector=SumFluxes.normMax()
        isStationary= maxVector[0]/p0<precision and maxVector[1]/rho0<precision and maxVector[2]/rho0<precision;
        U-=SumFluxes;
    
        time=time+dt;
        it=it+1;
 
         #Sauvegardes
        if(it%output_freq==0 or it>=ntmax or isStationary or time >=tmax):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
            print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0

            for k in range(nbCells):
                pressure_field[k]=U[k,0]
                velocity_field[k,0]=U[k,1]/rho0
                if(dim>1):
                    velocity_field[k,1]=U[k,2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=U[k,3]/rho0

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity",False);

    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
        diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem"+str(dim)+"DUpwind"+meshName+"_velocity_Stat")
        
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"


def solve(my_mesh,filename,resolution):
    print("Resolution of the Wave system with Wall boundary conditions:")

    # Problem data
    tmax = 1.
    ntmax = 100
    cfl = 1./my_mesh.getSpaceDimension()
    output_freq = 100

    WaveSystemVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution)

def solve_file( filename,resolution):
    my_mesh = cdmath.Mesh(filename+".med")
    solve(my_mesh, filename,resolution)
    
if __name__ == """__main__""":
    if len(sys.argv) >1 :
        filename=sys.argv[1]
        my_mesh = cdmath.Mesh(filename)
        solve(my_mesh,filename,100)
    else :
        raise ValueError("WaveSystemUpwind.py expects a mesh file name")
