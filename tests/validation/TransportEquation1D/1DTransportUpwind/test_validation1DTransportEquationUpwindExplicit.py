import cdmath
import 1DTransportEquationUpwindExplicit
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt
import sys
import time, json

    
def test_validation1DTransportEquationUpwindExplicit(cfl):
    start = time.time()
    #### 2D triangular meshes
    meshList=[10,100,500,1000]
    meshType="1D regular grid"
    testColor="Green"
    nbMeshes=len(meshList)
    mesh_size_tab=[0]*nbMeshes
    mesh_name='RegularGrid'

    error_u_tab=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    max_u=[0]*nbMeshes

    plt.close('all')
    i=0

    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], t_final[i], ndt_final[i], max_vel[i], diag_data_press[i], diag_data_vel[i], time_tab[i], cond_number[i] = 1DTransportEquationUpwindExplicit.solve_file(mesh_path+filename, mesh_name, resolution,scaling,meshType,testColor,cfl,"Periodic")
        assert max_u[i]>0.94 and max_u[i]<1.1
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    end = time.time()

    # Plot over diagonal line
    for i in range(nbMeshes):
        if(scaling==0):
            plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells - no scaling')
        else:
            plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells - with scaling')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for the stationary wave system \n with centered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_PlotOverDiagonalLine.png")    
    plt.close()

    # Plot of number of maximal unknown
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm - no scaling')
    else:
        plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm - with scaling')
    plt.legend()
    plt.xlabel('Number of cells')
    plt.ylabel('Max velocity norm')
    plt.title('Maximum velocity norm for the stationary Wave System \n with centered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_MaxVelNorm.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=0.5*log10(mesh_size_tab[i])
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemBrickWallFVCentered() : Make sure you use distinct meshes and at least two meshes'

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "Explicit Upwind scheme for Transport Equation on 1D regular grid : scheme order is ", -au
    
    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='|error on stationary pressure| ')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(|error p|)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with centered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Pressure_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ConvergenceCurve.png")
    
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|) - no scaling')
    else:
        plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|) - with scaling')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(|error u|)')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with centered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    if(scaling==0):
        plt.plot(mesh_size_tab, time_tab, label='log(cpu time) - no scaling')
    else:
        plt.plot(mesh_size_tab, time_tab, label='log(cpu time) - with scaling')
    plt.legend()
    plt.xlabel('1/2 log(Number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the stationary Wave System \n with centered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"2DWaveSystemCentered_"+"scaling"+str(scaling)+"_ComputationalTime.png")

    plt.close('all')

    convergence_synthesis={}

    convergence_synthesis["PDE_model"]="Wave system"
    convergence_synthesis["PDE_is_stationary"]=False
    convergence_synthesis["PDE_search_for_stationary_solution"]=True
    convergence_synthesis["Numerical_method_name"]="Upwind scheme"
    convergence_synthesis["Numerical_method_space_discretization"]="Finite volumes"
    convergence_synthesis["Numerical_method_time_discretization"]="Implicit"
    convergence_synthesis["Initial_data"]="Constant pressure, divergence free velocity"
    convergence_synthesis["Boundary_conditions"]="Periodic"
    convergence_synthesis["Numerical_parameter_cfl"]=cfl
    convergence_synthesis["Space_dimension"]=2
    convergence_synthesis["Mesh_dimension"]=2
    convergence_synthesis["Mesh_names"]=meshList
    convergence_synthesis["Mesh_type"]=meshType
    convergence_synthesis["Mesh_path"]=mesh_path
    convergence_synthesis["Mesh_description"]=mesh_name
    convergence_synthesis["Mesh_sizes"]=mesh_size_tab
    convergence_synthesis["Mesh_cell_type"]="BrickWall"
    convergence_synthesis["Numerical_error_velocity"]=error_u_tab
    convergence_synthesis["Numerical_error_pressure"]=error_p_tab
    convergence_synthesis["Max_vel_norm"]=max_vel
    convergence_synthesis["Final_time"]=t_final  
    convergence_synthesis["Final_time_step"]=ndt_final  
    convergence_synthesis["Scheme_order"]=-au
    convergence_synthesis["Scheme_order_vel"]=-au
    convergence_synthesis["Scaling_preconditioner"]=scaling
    convergence_synthesis["Condition_numbers"]=cond_number
    convergence_synthesis["Test_color"]=testColor
    convergence_synthesis["Computational_time"]=end-start

    with open('Convergence_1DTransportEquationUpwindExplicit_'+mesh_name+'.json', 'w') as outfile:  
        json.dump(convergence_synthesis, outfile)

if __name__ == """__main__""":
    if len(sys.argv) >1 :
        cfl = float(sys.argv[1])
        test_validation1DTransportEquationUpwindExplicit(cfl)
    else :
        test_validation1DTransportEquationUpwindExplicit(50)

