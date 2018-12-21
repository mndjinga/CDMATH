import cdmath
import 1DTransportEquationUpwindExplicit
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

    
def test_validation1DTransportEquationUpwindExplicit(cfl,isSmooth):
    start = time.time()
    #### 1D regular grid
    meshList=[10,100,500,1000]
    meshType="1D regular grid"
    testColor="Green"
    nbMeshes=len(meshList)
    mesh_size_tab=[0]*nbMeshes
    mesh_name='RegularGrid'

    a=0  ;  b=1
    error_u_tab=[0]*nbMeshes
    sol_u=[0]*nbMeshes
    total_var_u=[0]*nbMeshes
    min_u=[0]*nbMeshes
    max_u=[0]*nbMeshes
    time_tab=[0]*nbMeshes

    plt.close('all')
    i=0

    # Storing of numerical errors, mesh sizes and solution
    for nx in meshList:
        min_u[i], max_u[i], sol_u[i], total_var_u[i], error_u_tab[i], time_tab[i] = 1DTransportEquationUpwindExplicit.solve_file(nx,cfl,a,b,isSmooth)
        assert max_u[i]>0. and max_u[i]<1.
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    # Plot of solution
    for i in range(nbMeshes):
        plt.plot(np.linspace(a,b,meshList[i]), sol_u[i],   label= str(mesh_size_tab[i]) + ' cells')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Plot of the numerical solution of the transport equation \n with explicit upwind scheme on a 1D regular grid')
    plt.savefig(mesh_name+"_1DTransportEquationUpwind"+"_PlotOfSolution.png")    
    plt.close()

    # Plot of maximal value
    plt.close()
    plt.plot(mesh_size_tab, max_u, label='Maximum value')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max |u|')
    plt.title('Maximum velocity norm of the transport equation \n with explicit upwind scheme on a 1D regular grid')
    plt.savefig(mesh_name+"_1DTransportEquationUpwind"+"_MaxSolution.png")
    
    # Plot of total variation
    plt.close()
    plt.plot(mesh_size_tab, total_var_u, label='Total variation')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Var(u)')
    plt.title('Total variation for the transport equation \n with explicit upwind scheme on a 1D regular grid')
    plt.savefig(mesh_name+"_1DTransportEquationUpwind"+"_TotalVariation.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=log10(mesh_size_tab[i])
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation1DTransportEquationUpwind() : Make sure you use distinct meshes and at least two meshes'

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "Explicit Upwind scheme for Transport Equation on 1D regular grid : scheme order is ", -au
    
    # Plot of convergence curve
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='log(|error|)')
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(|error u|)')
    plt.title('Convergence of finite volumes for the transport equation \n with explicit upwind scheme on a 1D regular grid')
    plt.savefig(mesh_name+"1DTransportEquationUpwind_"+"_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(Number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the transport equation \n with explicit upwind scheme on a 1D regular grid')
    plt.savefig(mesh_name+"1DTransportEquationUpwind_"+"_ComputationalTime.png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Transport_Equation"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    convergence_synthesis["Numerical_method_name"]="Upwind scheme"
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]="Explicit"
    convergence_synthesis["Initial_data"]="sine"
    convergence_synthesis["Boundary_conditions"]="Periodic"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="1D regular grid"
    convergence_synthesis["Numerical_ersolution"]=max_u
    convergence_synthesis["Final_time"]=t_final  
    convergence_synthesis["Final_time_step"]=ndt_final  
    convergence_synthesis["Scheme_order"]=-au
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_1DTransportEquationUpwindExplicit_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >2 :
        cfl = float(sys.argv[1])
        isSmooth = float(sys.argv[2])
        test_validation1DTransportEquationUpwindExplicit(cfl,isSmooth)
    else :
        test_validation1DTransportEquationUpwindExplicit(0.99,True)

