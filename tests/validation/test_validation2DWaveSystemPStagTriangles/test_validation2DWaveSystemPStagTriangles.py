import cdmath
import WaveSystemPStag
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt

    
def test_validation2DWaveSystemPStagTriangles():
    #### 2D triangular meshes
    meshList=['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5','triangleMeshSquare_6']
    nbMeshes=len(meshList)
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../../ressources/2DTriangles/'
    mesh_name='meshSquareWithTriangles'
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)

    error_p_tab_noscaling=[0]*nbMeshes
    error_u_tab_noscaling=[0]*nbMeshes
    diag_data_press_noscaling=[0]*nbMeshes
    diag_data_vel_noscaling=[0]*nbMeshes
    time_tab_noscaling=[0]*nbMeshes
    t_final_noscaling=[0]*nbMeshes
    ndt_final_noscaling=[0]*nbMeshes
    max_vel_noscaling=[0]*nbMeshes
    cond_number_noscaling=[0]*nbMeshes

    error_p_tab_scaling=[0]*nbMeshes
    error_u_tab_scaling=[0]*nbMeshes
    diag_data_press_scaling=[0]*nbMeshes
    diag_data_vel_scaling=[0]*nbMeshes
    time_tab_scaling=[0]*nbMeshes
    t_final_scaling=[0]*nbMeshes
    ndt_final_scaling=[0]*nbMeshes
    max_vel_scaling=[0]*nbMeshes
    cond_number_scaling=[0]*nbMeshes

    plt.close('all')
    i=0

    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_p_tab_noscaling[i], error_u_tab_noscaling[i], mesh_size_tab[i], t_final_noscaling[i], ndt_final_noscaling[i], max_vel_noscaling[i], diag_data_press_noscaling[i], diag_data_vel_noscaling[i], time_tab_noscaling[i], cond_number_noscaling[i] = WaveSystemPStag.solve_file(mesh_path+filename, mesh_name, resolution,0)
        print max_vel_noscaling[i]
        assert max_vel_noscaling[i]>0.94 and max_vel_noscaling[i]<1.1
        error_p_tab_noscaling[i]=log10(error_p_tab_noscaling[i])
        error_u_tab_noscaling[i]=log10(error_u_tab_noscaling[i])
        time_tab_noscaling[i]=log10(time_tab_noscaling[i])

        error_p_tab_scaling[i],   error_u_tab_scaling[i],   mesh_size_tab[i],   t_final_scaling[i], ndt_final_scaling[i],   max_vel_scaling[i],   diag_data_press_scaling[i],   diag_data_vel_scaling[i],   time_tab_scaling[i],   cond_number_scaling[i] = WaveSystemPStag.solve_file(mesh_path+filename, mesh_name, resolution,2)
        print max_vel_scaling[i]
        assert max_vel_scaling[i]>0.94 and max_vel_scaling[i]<1.1
        error_p_tab_scaling[i]=log10(error_p_tab_scaling[i])
        error_u_tab_scaling[i]=log10(error_u_tab_scaling[i])
        time_tab_scaling[i]=log10(time_tab_scaling[i])

        i=i+1
    
    # Plot over diagonal line
    for i in range(nbMeshes):
        plt.plot(curv_abs, diag_data_press_noscaling[i], label= str(mesh_size_tab[i]) + ' cells - no scaling')
        plt.plot(curv_abs, diag_data_press_scaling[i], label= str(mesh_size_tab[i]) + ' cells - with scaling')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Pressure on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+'_Pressure_2DWaveSystemSquaresPStag_'+"PlotOverDiagonalLine.png")
    plt.close()

    plt.clf()
    for i in range(nbMeshes):
        plt.plot(curv_abs, diag_data_vel_noscaling[i],   label= str(mesh_size_tab[i]) + ' cells - no scaling')
        plt.plot(curv_abs, diag_data_vel_scaling[i],   label= str(mesh_size_tab[i]) + ' cells - with scaling')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for the stationary wave system \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemSquaresPStag_"+"PlotOverDiagonalLine.png")    
    plt.close()

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemSquaresFVPStag() : Make sure you use distinct meshes and at least two meshes'

    b1u_noscaling=np.dot(error_u_tab_noscaling,mesh_size_tab)   
    b2u_noscaling=np.sum(error_u_tab_noscaling)
    au_noscaling=( a3*b1u_noscaling-a2*b2u_noscaling)/det
    bu_noscaling=(-a2*b1u_noscaling+a1*b2u_noscaling)/det
    
    print "FVPStag on 2D triangular meshes : scheme order for velocity without scaling is ", -au_noscaling
    
    b1u_scaling=np.dot(error_u_tab_scaling,mesh_size_tab)   
    b2u_scaling=np.sum(error_u_tab_scaling)
    au_scaling=( a3*b1u_scaling-a2*b2u_scaling)/det
    bu_scaling=(-a2*b1u_scaling+a1*b2u_scaling)/det
    
    print "FVPStag on 2D triangular meshes : scheme order for velocity with scaling is ", -au_scaling
    
    # Plot of number of time steps
    plt.close()
    plt.plot(mesh_size_tab, ndt_final_noscaling, label='Number of time step to reach stationary regime - no scaling')
    plt.plot(mesh_size_tab, ndt_final_scaling, label='Number of time step to reach stationary regime - with scaling')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max time steps for stationary regime')
    plt.title('Number of times steps required for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquarePStags_"+"TimeSteps.png")
    
    # Plot of number of stationary time
    plt.close()
    plt.plot(mesh_size_tab, t_final_noscaling, label='Time where stationary regime is reached - no scaling')
    plt.plot(mesh_size_tab, t_final_scaling, label='Time where stationary regime is reached - with scaling')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max time for stationary regime')
    plt.title('Simulated time for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquaresPStag_"+"TimeFinal.png")
    
    # Plot of number of maximal velocity norm
    plt.close()
    plt.plot(mesh_size_tab, max_vel_noscaling, label='Maximum velocity norm - no scaling')
    plt.plot(mesh_size_tab, max_vel_scaling, label='Maximum velocity norm - with scaling')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max velocity norm')
    plt.title('Maximum velocity norm for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquaresPStag_"+"MaxVelNorm.png")
    
    # Plot of condition number 
    plt.close()
    plt.plot(mesh_size_tab, cond_number_noscaling, label='Maximum velocity norm - no scaling')
    plt.plot(mesh_size_tab, cond_number_scaling, label='Maximum velocity norm - with scaling')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Condition number')
    plt.title('Condition number for the stationary Wave System \n with PStagggered scheme on 2D square meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquaresPStag_"+"condition_number.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=log10(mesh_size_tab[i])
        
    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab_noscaling, label='|error on stationary pressure| - no scaling')
    plt.plot(mesh_size_tab, error_p_tab_scaling, label='|error on stationary pressure| - with scaling')
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('|error p|')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Pressure_2DWaveSystemSquaresPStag_"+"ConvergenceCurve.png")
    
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab_noscaling, label='log(|error on stationary velocity|) - no scaling')
    plt.plot(mesh_size_tab, error_u_tab_scaling, label='log(|error on stationary velocity|) - with scaling')
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('|error u|')
    plt.title('Convergence of finite volumes for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemSquaresPStag_"+"ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab_noscaling, label='log(cpu time) - no scaling')
    plt.plot(mesh_size_tab, time_tab_scaling, label='log(cpu time) - with scaling')
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes for the stationary Wave System \n with PStagggered scheme on 2D triangular meshes')
    plt.savefig(mesh_name+"2DWaveSystemSquaresPStag_ComputationalTimeSquares.png")

    plt.close('all')

if __name__ == """__main__""":
    test_validation2DWaveSystemPStagTriangles()
