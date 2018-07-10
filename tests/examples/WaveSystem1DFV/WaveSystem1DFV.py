#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import cdmath
import PV_routines
import VTK_routines

rho0=1000#reference density
c0=1500#reference sound speed
p0=rho0*c0*c0#reference pressure
precision=1e-5

def initial_conditions_wave_system(my_mesh):
    nbCells = my_mesh.getNumberOfCells()
    dim=1
    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, dim)
    U              = cdmath.Field("Conservative vector", cdmath.CELLS, my_mesh, dim+1)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()

        
        pressure_field[i] = p0
        if(x>0.5):
            velocity_field[i,0] =   1
        else:    
            velocity_field[i,0] =  -1
            
        U[i,0] =   p0
        U[i,1] =  rho0*velocity_field[i,0]
        
    return U, pressure_field, velocity_field

def jacobianMatrices():
    A=cdmath.Matrix(2,2)
    absA=cdmath.Matrix(2,2)

    absA[0,0]=c0
    absA[1,1]=c0
    A[1,0]=1
    A[0,1]=c0*c0
    
    return A, absA
    
def Flux(U):

    result=cdmath.Vector(2)
    result[0] = U[1]
    result[1] = U[0]
        
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
                else :
                    # hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                
                for i in range(nbComp):
                    Uautre[i]=U[cellAutre,i]
            else :
                if(Fk.getGroupName() == "Wall" or Fk.getGroupName() == "Paroi" or Fk.getGroupName() == "Haut" or Fk.getGroupName() == "Bas" or Fk.getGroupName() == "Gauche" or Fk.getGroupName() == "Droite"):#Wall boundary condition unless Neumannspecified explicitly
                    Uautre=Ucourant;
                    qn=0# normal momentum
                    for i in range(dim):
                        qn+=Ucourant[i+1]*normal[i]
                    for i in range(dim):
                        Uautre[i+1]-=2*qn*normal[i]
                elif(Fk.getGroupName() == "Neumann"):
                    Uautre=Ucourant;
                else:
                    print Fk.getGroupName()
                    raise ValueError("computeFluxes: Unknown boundary condition name");
            
            Fcourant=Flux(Ucourant);
            Fautre  =Flux(Uautre);

            A, absA=jacobianMatrices();
            
            sumFluxCourant = sumFluxCourant + (Fcourant+Fautre +absA*(Ucourant-Uautre))*Fk.getMeasure()*0.5
 
        #On divise par le volume de la cellule la contribution des flux au snd membre
        for i in range(nbComp):
            SumFluxes[j,i]=sumFluxCourant[i]/Cj.getMeasure();


def WaveSystem1DVF(ntmax, tmax, cfl, my_mesh, output_freq, resolution):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    SumFluxes = cdmath.Field("Fluxes", cdmath.CELLS, my_mesh, dim+1)

    # Initial conditions #
    print("Construction of the initial condition …")
    U, pressure_field, velocity_field = initial_conditions_wave_system(my_mesh)

    #keep the initial data in memory to measure the error later
    Uinitial=U
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("velocity");

    dx_min=my_mesh.minRatioSurfVol()

    dt = cfl * dx_min / c0
    
    print("Starting computation of the linear wave system with an UPWIND scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        computeFluxes(U,SumFluxes);
        
        SumFluxes*=dt;
        maxVector=SumFluxes.normMax()
        isStationary= maxVector[0]/p0<precision and maxVector[1]/rho0<precision
        U-=SumFluxes;
    
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
            print "|| Un+1 - Un || : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 

            for k in range(nbCells):
                pressure_field[k]=U[k,0]/p0
                velocity_field[k,0]=U[k,1]/rho0

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("velocity",False);
    
    print "|| Un+1 - Un || : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        for k in range(nbCells):
            pressure_field[k]=U[k,0]/p0
            velocity_field[k,0]=U[k,1]/rho0

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("velocity_Stat");

        maxVector=(Uinitial-U).normMax()
        error_p=maxVector[0]/p0
        error_u=maxVector[1]/rho0
        
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh,resolution):
    print("RESOLUTION OF THE 1D Wave system:")

    # Problem data
    tmax = 1.
    ntmax = 100
    cfl = 0.95
    output_freq = 10

    WaveSystem1DVF(ntmax, tmax, cfl, my_mesh, output_freq,resolution)

if __name__ == """__main__""":

    M=cdmath.Mesh(0,1,50)
    
    xinf=0
    xsup=1
 
    M=cdmath.Mesh(xinf,xsup,50)
    M.setGroupAtPlan(xsup,0,precision,"Neumann");
    M.setGroupAtPlan(xinf,0,precision,"Neumann");

    solve(M,100)
