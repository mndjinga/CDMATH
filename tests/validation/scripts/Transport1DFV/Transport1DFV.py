#!/usr/bin/env python
# -*-coding:utf-8 -*

import cdmath
import time, json

test_desc={}

rho0=1000#reference density
c0=1500#reference sound speed
p0=rho0*c0*c0#reference pressure
precision=1e-5

def initial_conditions_wave_system(my_mesh):
    test_desc["Initial_data"]="Schock"
    dim     = my_mesh.getMeshDimension()
    nbCells = my_mesh.getNumberOfCells()
    if(dim!=1):
        raise ValueError("initial_conditions_wave_system: Mesh dimension should be 1")
        
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
    test_desc["Numerical_method_name"]="Upwind"
    A=cdmath.Matrix(2,2)
    absA=cdmath.Matrix(2,2)

    absA[0,0]=c0
    absA[1,1]=c0
    A[1,0]=1
    A[0,1]=c0*c0
    
    return A, absA
    
def Flux(U):

    result=cdmath.Vector(2)
    result[0] = c0*c0*U[1]
    result[1] = U[0]
            
    return result
    
def numericalFlux(Uj,Ujp1,absA):

    Fj   = Flux(Uj)
    Fjp1 = Flux(Ujp1)
            
    return Fj+Fjp1 +absA*(Uj-Ujp1)
        
def computeFluxes(U, SumFluxes):
    my_mesh =U.getMesh();
    nbCells = my_mesh.getNumberOfCells();
    dim=my_mesh.getMeshDimension();
    nbComp=U.getNumberOfComponents();
    Fj=cdmath.Vector(nbComp)
    Fjp1=cdmath.Vector(nbComp)
    Fjm1=cdmath.Vector(nbComp)
    Uj=cdmath.Vector(nbComp)
    Ujp1=cdmath.Vector(nbComp)
    Ujm1=cdmath.Vector(nbComp)
    normal=cdmath.Vector(dim)
    sumFluxCourant=cdmath.Vector(nbComp)

    A, absA=jacobianMatrices();

    for j in range(nbCells):#On parcourt les cellules
        for i in range(nbComp) :
            Uj[i]=U[j,i];
            sumFluxCourant[i]=0;

        test_desc["Boundary_conditions"]="Neumann"
        if ( j==0) :
            for i in range(nbComp) :
                Ujp1[i]=U[j+1,i];
                Ujm1[i]=U[j  ,i];
        elif ( j==nbCells-1) :
            for i in range(nbComp) :
                Ujp1[i]=U[j  ,i];
                Ujm1[i]=U[j-1,i];
        else :
            for i in range(nbComp) :
                Ujp1[i]=U[j+1,i];
                Ujm1[i]=U[j-1,i];
            
        Fr=numericalFlux(Uj,Ujp1,absA)
        Fl=numericalFlux(Ujm1,Uj,absA)

        sumFluxCourant = (Fr - Fl)*0.5/Cj.getMeasure()
            
            Fcourant=Flux(Ucourant,normal);
            Fautre  =Flux(Uautre,  normal);

            A, absA=jacobianMatrices( normal);
            
            sumFluxCourant = sumFluxCourant + (Fcourant+Fautre +absA*(Ucourant-Uautre))*Fk.getMeasure()*0.5

        #On divise par le volume de la cellule la contribution des flux au snd membre
        for i in range(nbComp):
            SumFluxes[j,i]=sumFluxCourant[i];


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
            print
    
            for k in range(nbCells):
                pressure_field[k]=U[k,0]/p0
                velocity_field[k,0]=U[k,1]/rho0

            pressure_field.setTime(time,it);
            pressure_field.writeCSV("pressure");
            velocity_field.setTime(time,it);
            velocity_field.writeCSV("velocity");
    
    print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
    print "|| Un+1 - Un || : pressure ", maxVector[0]/p0 ,", velocity x", maxVector[1]/rho0 
    print

    if(it>=ntmax):
        print "Nombre de pas de temps maximum ntmax= ", ntmax, " atteint"
        raise ValueError("Maximum number of time steps reached : Stationary state not found !!!!!!!")
    elif(isStationary):
        print "Régime stationnaire atteint au pas de temps ", it, ", t= ", time
        for k in range(nbCells):
            pressure_field[k]=U[k,0]/p0
            velocity_field[k,0]=U[k,1]/rho0

        pressure_field.setTime(time,0);
        pressure_field.writeCSV("pressure_Stat");
        velocity_field.setTime(time,0);
        velocity_field.writeCSV("velocity_Stat");
        
    else:
        print "Temps maximum Tmax= ", tmax, " atteint"
        raise ValueError("Maximum time reached : Stationary state not found !!!!!!!")

    return time, it

def solve(my_mesh,resolution):
    start = time.time()
    test_name="Resolution of the Wave system in dimension " +str( my_mesh.getSpaceDimension())+" on "+str(my_mesh.getNumberOfCells())+ " cells"
    test_name_comment="Classical characteristic based scheme"
    test_model="Wave system"
    test_method="Upwind"
    test_initial_data="Schock"
    test_bc="Neumann"

    print test_name
    print "Numerical method : ", test_method
    print "Initial data : ", test_initial_data
    print "Boundary conditions : ",test_bc
    print "Mesh name : ",meshName , ", ", my_mesh.getNumberOfCells(), " cells"
 
    # Problem data
    tmax = 1.
    ntmax = 100
    cfl = 0.95
    output_freq = 10

    t_final, ndt_final = WaveSystem1DVF(ntmax, tmax, cfl, my_mesh, output_freq,resolution)
    end = time.time()

    test_desc["Global_name"]=test_name
    test_desc["Global_comment"]=test_name_comment
    test_desc["PDE_model"]=test_model
    test_desc["PDE_is_stationary"]=False
    test_desc["PDE_search_for_stationary_solution"]=False
    test_desc["Numerical_method_name"]=test_method
    test_desc["Numerical_method_space_discretization"]="Finite volumes"
    test_desc["Numerical_method_time_discretization"]="Implicit"
    test_desc["Space_dimension"]=my_mesh.getSpaceDimension()
    test_desc["Mesh_dimension"]=my_mesh.getMeshDimension()
    test_desc["Mesh_is_unstructured"]=False
    test_desc["Mesh_cell_type"]="Quadrangles"
    test_desc["Mesh_number_of_elements"]=my_mesh.getNumberOfCells()
    test_desc["Mesh_max_number_of_neighbours"]=2
    test_desc["Geometry"]="Interval"
    test_desc["Boundary_conditions"]=test_bc
    test_desc["Initial_data"]=test_initial_data
    test_desc["Part_of_mesh_convergence_analysis"]=True
    test_desc["Numerical_parameter_cfl"]=cfl
    test_desc["Simulation_parameter_maximum_time_step"]=ntmax
    test_desc["Simulation_parameter_maximum_time"]=tmax
    test_desc["Simulation_output_frequency"]=output_freq
    test_desc["Simulation_final_time_after_run"]=t_final
    test_desc["Simulation_final_number_of_time_steps_after_run"]=ndt_final
    test_desc["Computational_time_taken_by_run"]=end-start
    test_desc["Part_of_mesh_convergence_analysis"]=True

    with open('Transport'+str(my_mesh.getMeshDimension())+'DFV_'+meshName+ "Cells.json", 'w') as outfile:  
        json.dump(test_desc, outfile)

if __name__ == """__main__""":

    xinf=0
    xsup=1
 
    M=cdmath.Mesh(xinf,xsup,50)

    solve(M,100)
