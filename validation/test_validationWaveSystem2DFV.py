import cdmath
import WaveSystem2DFV
import matplotlib.pyplot as plt
import numpy as np
from math import log10, sqrt

def test_validation2DWaveSystemTrianglesFV():
    #### 2D triangle mesh
    meshList=['triangleMeshSquare_1','triangleMeshSquare_2','triangleMeshSquare_3','triangleMeshSquare_4','triangleMeshSquare_5']
    nbMeshes=len(meshList)
    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../validation/2DTriangles/'
    mesh_name='meshSquareWithTrianglesFV'
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    # Storing of numerical errors, mesh sizes and diagonal values
    for filename in meshList:
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], diag_data_press[i], diag_data_vel[i] =WaveSystem2DFV.solve_file(mesh_path+filename, resolution)
        plt.figure('pressure')
        plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells')
        plt.figure('velocity')
        plt.plot(curv_abs, diag_data_vel[i], label= str(mesh_size_tab[i]) + ' cells')
        error_p_tab[i]=log10(error_p_tab[i])
        error_u_tab[i]=log10(error_u_tab[i])
        mesh_size_tab[i] = log10(mesh_size_tab[i])
        i=i+1
    
    # Plot over diagonal line
    plt.figure('pressure')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Pressure on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system on 2D triangular meshes')
    plt.savefig(mesh_name+'_Pressure_'+"PlotOverDiagonalLine.png")

    plt.figure('velocity')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_"+"PlotOverDiagonalLine.png")
    
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemTrianglesFV() : Make sure you use distinct meshes and at least two meshes'

    b1p=np.dot(error_p_tab,mesh_size_tab)   
    b2p=np.sum(error_p_tab)
    
    ap=( a3*b1p-a2*b2p)/det
    bp=(-a2*b1p+a1*b2p)/det
    
    print "Stationary Wave System on 2D triangle meshes : pressure scheme order is ", -ap
    assert abs(ap+1)<0.1

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "Stationary Wave System on 2D triangle meshes : scheme order is ", -au
    assert abs(au+1)<0.1

    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='log(|error on stationary pressure|)')
    plt.plot(mesh_size_tab, ap*np.array(mesh_size_tab)+bp,label='least square slope : '+'%.3f' % ap)
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(|error|)')
    plt.title('Convergence of finite volumes for stationary Wave System on 2D triangular meshes')
    plt.savefig(mesh_name+"_Pressure_"+"ConvergenceCurve.png")
    
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|)')
    plt.plot(mesh_size_tab, au*np.array(mesh_size_tab)+bu,label='least square slope : '+'%.3f' % au)
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(|error|)')
    plt.title('Convergence of finite volumes for stationary Wave System on 2D triangular meshes')
    plt.savefig(mesh_name+"_Velocity_"+"ConvergenceCurve.png")
    
def test_validation2DWaveSystemSquaresFV():
    #### 2D square mesh
    meshList=[11,51,151,201]
    nbMeshes=len(meshList)
    error_p_tab=[0]*nbMeshes
    error_u_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_name='meshSquareWithSquaresFV'
    diag_data_press=[0]*nbMeshes
    diag_data_vel=[0]*nbMeshes
    resolution=100
    curv_abs=np.linspace(0,sqrt(2),resolution+1)
    plt.close('all')
    i=0
    # Storing of numerical errors, mesh sizes and diagonal values
    for nx in meshList:
        my_mesh=cdmath.Mesh(0,1,nx,0,1,nx)
        error_p_tab[i], error_u_tab[i], mesh_size_tab[i], diag_data_press[i], diag_data_vel[i] =WaveSystem2DFV.solve(my_mesh,str(nx)+'x'+str(nx), resolution)
        plt.figure('pressure')
        plt.plot(curv_abs, diag_data_press[i], label= str(mesh_size_tab[i]) + ' cells')
        plt.figure('velocity')
        plt.plot(curv_abs, diag_data_vel[i], label= str(mesh_size_tab[i]) + ' cells')
        error_p_tab[i]=log10(error_p_tab[i])
        error_u_tab[i]=log10(error_u_tab[i])
        mesh_size_tab[i] = log10(mesh_size_tab[i])
        i=i+1
    
    # Plot over diagonal line
    plt.figure('pressure')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Pressure on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system on a 2D square mesh')
    plt.savefig(mesh_name+'_Pressure_'+"PlotOverDiagonalLine.png")

    plt.figure('velocity')
    plt.legend()
    plt.xlabel('Position on diagonal line')
    plt.ylabel('Velocity on diagonal line')
    plt.title('Plot over diagonal line for stationary wave system on a 2D square mesh')
    plt.savefig(mesh_name+"_Velocity_"+"PlotOverDiagonalLine.png")
    
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DWaveSystemSquaresFV() : Make sure you use distinct meshes and at least two meshes'

    b1p=np.dot(error_p_tab,mesh_size_tab)   
    b2p=np.sum(error_p_tab)
    
    ap=( a3*b1p-a2*b2p)/det
    bp=(-a2*b1p+a1*b2p)/det
    
    print "Stationary Wave System on 2D square mesh : pressure scheme order is ", -ap
    assert abs(ap+1)<0.1

    b1u=np.dot(error_u_tab,mesh_size_tab)   
    b2u=np.sum(error_u_tab)
    
    au=( a3*b1u-a2*b2u)/det
    bu=(-a2*b1u+a1*b2u)/det
    
    print "Stationary Wave System on 2D square mesh : scheme order is ", -au
    assert abs(au+1)<0.1

    # Plot of convergence curves
    plt.close()
    plt.plot(mesh_size_tab, error_p_tab, label='log(|error on stationary pressure|)')
    plt.plot(mesh_size_tab, ap*np.array(mesh_size_tab)+bp,label='least square slope : '+'%.3f' % ap)
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(|error|)')
    plt.title('Convergence of finite volumes for stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_Pressure_"+"ConvergenceCurve.png")
    
    plt.close()
    plt.plot(mesh_size_tab, error_u_tab, label='log(|error on stationary velocity|)')
    plt.plot(mesh_size_tab, au*np.array(mesh_size_tab)+bu,label='least square slope : '+'%.3f' % au)
    plt.legend()
    plt.xlabel('log(number of cells)')
    plt.ylabel('log(|error|)')
    plt.title('Convergence of finite volumes for stationary Wave System on 2D square meshes')
    plt.savefig(mesh_name+"_Velocity_"+"ConvergenceCurve.png")
    
if __name__ == """__main__""":
    test_validation2DWaveSystemTrianglesFV()
    test_validation2DWaveSystemSquaresFV()
