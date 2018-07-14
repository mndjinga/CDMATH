import cdmath
import WaveSystem2DFV
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt

    
def test_validation2DWaveSystemSquaresFV():
    #### 2D square mesh
    meshList=[15,31,51]#,151]
    nbMeshes=len(meshList)
    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_name='meshSquareWithSquaresFV'
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    time_tab=[0]*nbMeshes
    t_final=[0]*nbMeshes
    ndt_final=[0]*nbMeshes
    max_vel=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0

    #plt.figure('pressure')
    #plt.figure('velocity')
    # Storing of numerical errors, mesh sizes and diagonal values
    for nx in meshList:
        my_mesh=cdmath.Mesh(0,1,nx,0,1,nx)
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], t_final[i], ndt_final[i], max_vel[i], time_tab[i] =WaveSystem2DFV.solve(my_mesh,str(nx)+'x'+str(nx), resolution)
        #plt.show('pressure')
        #plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells')
        #plt.close('pressure')
        #plt.show('velocity')
        #plt.plot(curv_abs, diag_data_vel[i],   label= str(mesh_size_tab[i]) + ' cells')
        #plt.close('velocity')
        error_u_tab[i]=log10(error_u_tab[i])
        time_tab[i]=log10(time_tab[i])
        i=i+1
    
    ## Plot over diagonal line
    #plt.show('pressure')
    #plt.legend()
    #plt.xlabel('Position on diagonal line')
    #plt.ylabel('Pressure on diagonal line')
    #plt.title('Plot over diagonal line for stationary wave system \n on 2D square meshes')
    #plt.savefig(mesh_name+'_Pressure_2DWaveSystemSquares_'+"PlotOverDiagonalLine.png")

    #plt.close('pressure')

    #plt.show('velocity')
    #plt.legend()
    #plt.xlabel('Position on diagonal line')
    #plt.ylabel('Velocity on diagonal line')
    #plt.title('Plot over diagonal line for the stationary wave system \n on 2D square meshes')
    #plt.savefig(mesh_name+"_Velocity_2DWaveSystemSquares_"+"PlotOverDiagonalLine.png")
    
    #plt.close('velocity')

    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemSquaresFV() : Make sure you use distinct meshes and at least two meshes'

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "FV on 2D square meshes : scheme order for velocity is ", -au
    
    # Plot of number of time steps
    plt.close()
    plt.plot(mesh_size_tab, ndt_final, label='Number of time step to reach stationary regime')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max time steps for stationary regime')
    plt.title('Number of times steps required \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquares_"+"TimeSteps.png")
    
    # Plot of number of stationary time
    plt.close()
    plt.plot(mesh_size_tab, t_final, label='Time where stationary regime is reached')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max time for stationary regime')
    plt.title('Simulated time  \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquares_"+"TimeFinal.png")
    
    # Plot of number of maximal velocity norm
    plt.close()
    plt.plot(mesh_size_tab, max_vel, label='Maximum velocity norm')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('Max velocity norm')
    plt.title('Maximum velocity norm  \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_2DWaveSystemSquares_"+"TimeFinal.png")
    
    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='|error on stationary pressure|')
    plt.legend()
    plt.xlabel('number of cells')
    plt.ylabel('|error p|')
    plt.title('Convergence of finite volumes \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_Pressure_2DWaveSystemSquares_"+"ConvergenceCurve.png")
    
    for i in range(nbMeshes):
        mesh_size_tab[i]=log10(mesh_size_tab[i])
        
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|)')
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('|log(error u)|')
    plt.title('Convergence of finite volumes \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_Velocity_2DWaveSystemSquares_"+"ConvergenceCurve.png")
    
    # Plot of computational time
    plt.close()
    plt.plot(mesh_size_tab, time_tab, label='log(cpu time)')
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(cpu time)')
    plt.title('Computational time of finite volumes \n for the stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"2DWaveSystemSquares_ComputationalTimeSquares.png")

    plt.close('all')

if __name__ == """__main__""":
    test_validation2DWaveSystemSquaresFV()
