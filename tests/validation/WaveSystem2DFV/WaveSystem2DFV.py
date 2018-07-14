#!/usr/bin/env python
# -*-coding:utf-8 -*

from math import sin, cos, pi, sqrt
import time
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

    for i in range(nbCells):
        x = my_mesh.getCell(i).x()
        y = my_mesh.getCell(i).y()

        pressure_field[i] = p0
        velocity_field[i,0] =  sin(pi*x)*cos(pi*y)
        velocity_field[i,1] = -sin(pi*y)*cos(pi*x)
        velocity_field[i,2] = 0
        
    return pressure_field, velocity_field

def jacobianMatrices(normal, coeff):
    dim=normal.size()
    A=cdmath.Matrix(dim+1,dim+1)
    absA=cdmath.Matrix(dim+1,dim+1)

    absA[0,0]=c0*coeff
    for i in range(dim):
        A[i+1,0]=normal[i]*coeff
        A[0,i+1]=c0*c0*normal[i]*coeff
        for j in range(dim):
            absA[i+1,j+1]=c0*normal[i]*normal[j]*coeff
    
    return (A-absA)/2
    
    
def computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt):
    nbCells = my_mesh.getNumberOfCells()
    dim=my_mesh.getMeshDimension()
    nbComp=dim+1
    normal=cdmath.Vector(dim)

    implMat=cdmath.SparseMatrix(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp*nbCells)
    #implMat=cdmath.SparseMatrixPetsc(nbCells*nbComp,nbCells*nbComp,(nbVoisinsMax+1)*nbComp)

    idMoinsJacCL=cdmath.Matrix(nbComp)
    
    for j in range(nbCells):#On parcourt les cellules
        Cj = my_mesh.getCell(j)
        nbFaces = Cj.getNumberOfFaces();

        for k in range(nbFaces) :
            indexFace = Cj.getFacesId()[k];
            Fk = my_mesh.getFace(indexFace);
            for i in range(dim) :
                normal[i] = Cj.getNormalVector(k, i);#normale sortante

            Am=jacobianMatrices( normal,dt*Fk.getMeasure()/Cj.getMeasure());

            cellAutre =-1
            if ( not Fk.isBorder()) :
                # hypothese: La cellule d'index indexC1 est la cellule courante index j */
                if (Fk.getCellsId()[0] == j) :
                    # hypothese verifiée 
                    cellAutre = Fk.getCellsId()[1];
                elif(Fk.getCellsId()[1] == j) :
                    # hypothese non verifiée 
                    cellAutre = Fk.getCellsId()[0];
                else :
                    raise ValueError("computeFluxes: problem with mesh, unknown cel number")
                    
                implMat.addValue(j*nbComp,cellAutre*nbComp,Am)
                implMat.addValue(j*nbComp,        j*nbComp,Am*(-1.))
            else  :
                if(Fk.getGroupName() == "" or Fk.getGroupName() == "Wall" or Fk.getGroupName() == "Paroi" or Fk.getGroupName() != "Neumann"):#Wall boundary condition unless Neumannspecified explicitly
                    v=cdmath.Vector(dim+1)
                    for i in range(dim) :
                        v[i+1]=normal[i]
                    idMoinsJacCL=v.tensProduct(v)*2
                    
                    implMat.addValue(j*nbComp,j*nbComp,Am*(-1.)*idMoinsJacCL)
                    
                elif(Fk.getGroupName() != "Neumann"):#Nothing to do for Neumann boundary condition
                    print Fk.getGroupName()
                    raise ValueError("computeFluxes: Unknown boundary condition name");
                
    return implMat

def WaveSystem2DVF(ntmax, tmax, cfl, my_mesh, output_freq,resolution):
    dim=my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    
    dt = 0.
    time = 0.
    it=0;
    isStationary=False;
    
    nbVoisinsMax=10;
    iterGMRESMax=50
    isImplicit=False
    
    #iteration vectors
    Un=cdmath.Vector(nbCells*(dim+1))
    dUn=cdmath.Vector(nbCells*(dim+1))
    
    # Initial conditions #
    print("Construction of the initial condition …")
    pressure_field, velocity_field = initial_conditions_wave_system(my_mesh)
    initial_pressure, initial_velocity = initial_conditions_wave_system(my_mesh)

    for k in range(nbCells):
        Un[k*(dim+1)+0] =     initial_pressure[k]
        Un[k*(dim+1)+1] =rho0*initial_velocity[k,0]
        Un[k*(dim+1)+2] =rho0*initial_velocity[k,1]
            
    #sauvegarde de la donnée initiale
    pressure_field.setTime(time,it);
    pressure_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_pressure");
    velocity_field.setTime(time,it);
    velocity_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_velocity");

    total_pressure_initial=pressure_field.integral()[0]#For conservation test later
    
    dx_min=my_mesh.minRatioSurfVol()

    dt = cfl * dx_min / c0

    divMat=computeDivergenceMatrix(my_mesh,nbVoisinsMax,dt)
    if(isImplicit):
        #Adding the identity matrix on the diagonal
        for j in range(nbCells*nbComp):
            divMat.addValue(j,j,1)
        LS=cdmath.LinearSolver(divMat,Un,iterGMRESMax, precision, "GMRES","ILU")
    
    print("Starting computation of the linear wave system with an UPWIND scheme …")
    
    # Starting time loop
    while (it<ntmax and time <= tmax and not isStationary):
        if(isImplicit):
            dUn=Un.deepCopy()
            LS.setSndMember(Un)
            Un=LS.solve();
            cvgceLS=LS.getStatus();
            iterGMRES=LS.getNumberOfIter();
            if(not cvgceLS):
                print "Linear system did not converge ", iterGMRES, " GMRES iterations"
                raise ValueError("Pas de convergence du système linéaire");
            dUn-=Un
        else:
            dUn=divMat*Un
            Un-=dUn
        
        maxVector=dUn.maxVector(3)
        isStationary= maxVector[0]/p0<precision and maxVector[1]/rho0<precision and maxVector[2]/rho0<precision;
    
        time=time+dt;
        it=it+1;
    
        #Sauvegardes
        if(it%output_freq==0):
            print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
            print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
            if(isImplicit):
                print "Linear system converged in ", iterGMRES, " GMRES iterations"

            delta_press=0
            delta_velx=0
            delta_vely=0
            for k in range(nbCells):
                pressure_field[k]=Un[k*(dim+1)+0]
                velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
                if(dim>1):
                    velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                    if(dim>2):
                        velocity_field[k,2]=Un[k*(dim+1)+3]/rho0
                if (abs(initial_pressure[k]-pressure_field[k])>delta_press):
                    delta_press=abs(initial_pressure[k]-pressure_field[k])
                if (abs(initial_velocity[k,0]-velocity_field[k,0])>delta_velx):
                    delta_velx=abs(initial_velocity[k,0]-velocity_field[k,0])
                if (abs(initial_velocity[k,1]-velocity_field[k,1])>delta_vely):
                    delta_vely=abs(initial_velocity[k,1]-velocity_field[k,1])
            delta_vel=sqrt(delta_velx*delta_velx+delta_vely*delta_vely)

            pressure_field.setTime(time,it);
            pressure_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_pressure",False);
            velocity_field.setTime(time,it);
            velocity_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_velocity",False);

            print "Ecart au stationnaire exact : error_p= ",delta_press/p0," error_ux= ",delta_velx/rho0," error_uy= ",delta_vely/rho0
            print
    print"-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt)
    print "Variation temporelle relative : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 ,", velocity y", maxVector[2]/rho0
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        assert abs(total_pressure_initial-pressure_field.integral()[0])/p0<precision
        
        delta_press=0
        delta_velx=0
        delta_vely=0
        for k in range(nbCells):
            pressure_field[k]=Un[k*(dim+1)+0]
            velocity_field[k,0]=Un[k*(dim+1)+1]/rho0
            if(dim>1):
                velocity_field[k,1]=Un[k*(dim+1)+2]/rho0
                if(dim>2):
                    velocity_field[k,2]=Un[k*(dim+1)+3]/rho0
            if (abs(initial_pressure[k]-pressure_field[k])>delta_press):
                delta_press=abs(initial_pressure[k]-pressure_field[k])
            if (abs(initial_velocity[k,0]-velocity_field[k,0])>delta_velx):
                delta_velx=abs(initial_velocity[k,0]-velocity_field[k,0])
            if (abs(initial_velocity[k,1]-velocity_field[k,1])>delta_vely):
                delta_vely=abs(initial_velocity[k,1]-velocity_field[k,1])
        delta_vel=sqrt(delta_velx*delta_velx+delta_vely*delta_vely)

        pressure_field.setTime(time,0);
        pressure_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeVTK("WaveSystem2DFV"+str(nbCells)+"_velocity_Stat");

        #Postprocessing : Extraction of the diagonal data
        #diag_data_press=VTK_routines.Extract_field_data_over_line_to_numpyArray(pressure_field,[0,1,0],[1,0,0], resolution)    
        #diag_data_vel  =VTK_routines.Extract_field_data_over_line_to_numpyArray(velocity_field,[0,1,0],[1,0,0], resolution)    
        #Postprocessing : save 2D picture
        #PV_routines.Save_PV_data_to_picture_file("WaveSystem2DFV"+str(nbCells)+"_pressure_Stat"+'_0.vtu',"Pressure",'CELLS',"WaveSystem2DFV"+str(nbCells)+"_pressure_Stat")
        #PV_routines.Save_PV_data_to_picture_file("WaveSystem2DFV"+str(nbCells)+"_velocity_Stat"+'_0.vtu',"Velocity",'CELLS',"WaveSystem2DFV"+str(nbCells)+"_velocity_Stat")
        
        return delta_press/p0, delta_velx/rho0, nbCells, #diag_data_press, diag_data_vel
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")


def solve(my_mesh,filename,resolution):
    start = time.time()
    print("Resolution of the 2D Wave system:")

    # Problem data
    tmax = 1000.
    ntmax = 10000
    cfl = 0.45
    output_freq = 100

    error_p, error_u, nbCells = WaveSystem2DVF(ntmax, tmax, cfl, my_mesh, output_freq, resolution)
    end = time.time()

    return error_p, error_u, nbCells, end - start

def solve_file( filename,resolution):
    my_mesh = cdmath.Mesh(filename+".med")

    return solve(my_mesh, filename,resolution)
    

if __name__ == """__main__""":
    M=cdmath.Mesh(0,1,20,0,1,20)
    
    M.setGroupAtPlan(xsup,0,precision,"Wall");
    M.setGroupAtPlan(xinf,0,precision,"Wall");
    M.setGroupAtPlan(ysup,1,precision,"Wall");
    M.setGroupAtPlan(yinf,1,precision,"Wall");

    solve(M,"SquareRegularSquares",100)
