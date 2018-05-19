import cdmath
import FiniteVolumes2DWithCDMATH
import matplotlib.pyplot as plt
import numpy as np
from math import log10

def test_validation2DVF_ns_loc_ref():
    ##### 2D FV refined squares mesh
    nbMeshes=7
    error_tab=[0]*nbMeshes
    mesh_size_tab=[0]*nbMeshes
    mesh_path='../validation/2DLocRefinedSquares/'
    mesh_name='meshSquareWithLocRefSquaresFV'
    i=0
    # Computation of the numerical error
    for filename in ['meshLocRefSquares_1','meshLocRefSquares_2','meshLocRefSquares_3','meshLocRefSquares_4','meshLocRefSquares_5','meshLocRefSquares_6','meshLocRefSquares_7']:
        error_tab[i], mesh_size_tab[i] =FiniteVolumes2DWithCDMATH.solve_file(mesh_path+filename)
        error_tab[i]=log10(error_tab[i])
        mesh_size_tab[i] = log10(mesh_size_tab[i])
        i=i+1
        
    # Least square linear regression
    # Find the best a,b such that f(x)=ax+b best approximates the convergence curve
    # The vector X=(a,b) solves a symmetric linear system AX=B with A=(a1,a2\\a2,a3), B=(b1,b2)
    a1=np.dot(mesh_size_tab,mesh_size_tab)
    a2=np.sum(mesh_size_tab)
    a3=nbMeshes
    b1=np.dot(error_tab,mesh_size_tab)   
    b2=np.sum(error_tab)
    
    det=a1*a3-a2*a2
    assert det!=0, 'test_validation2DVF_ns_loc_ref() : Make sure you use distinct meshes'
    a=( a3*b1-a2*b2)/det
    b=(-a2*b1+a1*b2)/det
    
    print "FV on 2D refined squares mesh : scheme order is ", -a
    #assert abs(a-0.23)<0.1 #The scheme is not converging
    
    # Plot of figures
    plt.plot(mesh_size_tab, error_tab, label='log(|numerical-exact|)')
    plt.plot(mesh_size_tab, a*np.array(mesh_size_tab)+b,label='least square slope : '+'%.3f' % a)
    plt.legend()
    plt.xlabel('log(number of nodes)')
    plt.ylabel('log(error)')
    plt.title('Convergence of finite volumes for Laplace operator on a 2D refined squares mesh')
    plt.savefig(mesh_name+".png")

if __name__ == """__main__""":
    test_validation2DVF_ns_loc_ref()
