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
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()

    if(dim!=2):
        raise ValueError("initial_conditions_wave_system: Mesh dimension should be 2")

    pressure_field = cdmath.Field("Pressure",            cdmath.CELLS, my_mesh, 1)
    velocity_field = cdmath.Field("Velocity",            cdmath.CELLS, my_mesh, 3)
    U              = cdmath.Field("Conservative vector", cdmath.CELLS, my_mesh, dim+1)

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()

        pressure_field[i] = p0
        velocity_field[i,0] =  sin(pi*x)*cos(pi*y)
        velocity_field[i,1] = -sin(pi*y)*cos(pi*x)
        velocity_field[i,2] = 0

        U[i,0] =   p0
        U[i,1] =  rho0*sin(pi*x)*cos(pi*y)
        U[i,2] = -rho0*sin(pi*y)*cos(pi*x)
        
    return U, pressure_field, velocity_field

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
                if(Fk.getGroupName() == "" or Fk.getGroupName() == "Wall" or Fk.getGroupName() == "Paroi" or Fk.getGroupName() == "Haut" or Fk.getGroupName() == "Bas" or Fk.getGroupName() == "Gauche" or Fk.getGroupName() == "Droite"):
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


def WaveSystem2DVF(ntmax, tmax, cfl, my_mesh, output_freq, outputFileName,resolution):
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
    Uinitial = cdmath.Field("Initial vector", cdmath.CELLS, my_mesh, dim+1)
    for k in range(nbCells):
        for i in range(dim+1):
            Uinitial[k,i]=U[k,i]
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem2DFV"+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem2DFV"+"_velocity");

    dx_min=my_mesh.minRatioSurfVol()

    dt = cfl * dx_min / c0
    
    print("Starting computation of the linear wave system with an UPWIND scheme …")
    
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
        if(it%output_freq==0):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
            print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0

            for k in range(nbCells):
                pressure_field[k]=U[k,0]
                velocity_field[k,0]=U[k,1]/rho0
                if(dim>1):
                    velocity_field[k,1]=U[k,2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=U[k,3]/rho0
            print "total pressure ",pressure_field.integral()[0]
            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem2DFV"+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem2DFV"+"_velocity",False);

            maxVector=(Uinitial-U).normMax()
            error_p=maxVector[0]/p0
            error_u=sqrt(maxVector[1]*maxVector[1]+maxVector[2]*maxVector[2])/rho0

            print "max(|Pnum-Pexact|/p0)= ", error_p, "max(||Qnum-Qexact|/rho0)= ", error_u
            print
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        for k in range(nbCells):
            pressure_field[k]=U[k,0]
            velocity_field[k,0]=U[k,1]/rho0
            if(dim>1):
                velocity_field[k,1]=U[k,2]/rho0
                if(dim>2):
                    velocity_field[k,2]=U[k,3]/rho0

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem2DFV"+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem2DFV"+"_velocity_Stat");

        maxVector=(Uinitial-U).normMax()
        error_p=maxVector[0]/p0
        error_u=sqrt(maxVector[1]*maxVector[1]+maxVector[2]*maxVector[2])/rho0

        print "max(|Pnum-Pexact|/p0)= ", error_p, "max(||Qnum-Qexact|/rho0)= ", error_u
        
        #Postprocessing : Extraction of the diagonal data
        diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
        diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        #Postprocessing : save 2D picture
        PV_routines.Save_PV_data_to_picture_file("WaveSystem2DFV"+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem2DFV"+"_pressure_Stat")
        PV_routines.Save_PV_data_to_picture_file("WaveSystem2DFV"+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem2DFV"+"_velocity_Stat")
        
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh,filename,resolution):
    print("Resolution of the 2D Wave system:")

    # Problem data
    tmax = 1.
    ntmax = 10000
    cfl = 0.45
    output_freq = 100

    WaveSystem2DVF(ntmax, tmax, cfl, my_mesh, output_freq, filename,resolution)

def solve_file( filename,resolution):
    my_mesh = cdmath.Mesh(filename+".med")
    solve(my_mesh, filename,resolution)
    
if __name__ == """__main__""":
    M=cdmath.Mesh("meshSquare.med")
    solve(M,'SquaresWithTrianglesCells',100)

    xinf=0
    xsup=1
    yinf=0
    ysup=1
    M=cdmath.Mesh(xinf,xsup,15,yinf,ysup,15)
    
    M.setGroupAtPlan(xsup,0,precision,"Wall");
    M.setGroupAtPlan(xinf,0,precision,"Wall");
    M.setGroupAtPlan(ysup,1,precision,"Wall");
    M.setGroupAtPlan(yinf,1,precision,"Wall");

    solve(M,'SquaresWithSquareCells',100)
