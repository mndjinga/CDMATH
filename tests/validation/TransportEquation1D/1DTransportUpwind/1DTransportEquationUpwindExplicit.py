#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du transport 1D \partial_t u + c \partial_x u = 0 avec conditions aux limites périodiques
# Author      : Michaël Ndjinga, Katia Ait Ameur
# Copyright   : CEA Saclay 2018
# Description : Utilisation de la méthode des volumes P0 avec champs u discrétisé aux cellules d'un maillage 1D régulier
#		        Création et sauvegarde du champ résultant et des figures
#================================================================================================================================


from math import sin, pi, ceil
import numpy as np
import matplotlib.pyplot as plt
import sys

def solve(nx,cfl,a,b):
    start = time.time()
    ##################### Simulation parameters
    dx = (b - a) / nx #space step

    c = 0.25 # advection velocity
    tmax = (b-a)/c # runs the simulation for 0 <= t <= tMax
    dt = cfl * dx / c
    ntmax = ceil(tmax/dt)

    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    ########################## Initial data
    
    u_initial = [ 0.5*(1+sin(4*pi*xi-pi*.5))*int(xi<0.5)*int(0<xi) + int(0.6<xi)*int(xi<0.85)  for xi in x];# to be used with a=0, b=1

    max_initial=max(u_initial)
    min_initial=min(u_initial)
    total_var_initial = np.sum([abs(u_initial[i] - u_initial[(i-1)%nx]) for i in range(nx)])

    time = 0.
    it = 0
    output_freq = 10

    u=u_initial

    ########################### Postprocessing initialisation
    # Picture frame
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('u')
    plt.xlim(a,b)
    plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.1*(max_initial-min_initial) )
    plt.title('Upwind scheme for transport equation')
    line1, = plt.plot(x, u, label='u') #new picture for video # Returns a tuple of line objects, thus the comma

    writer.grab_frame()
    plt.savefig("TransportEquation_UpwindScheme_"+str(nx)+"Cells_ResultField_"+str(it)+".png")

    print("Saving initial data at T=0")
    np.savetxt("TransportEquation_UpwindScheme_"+str(nx)+"Cells_ResultField_0.txt", u, delimiter="\n")

    ############################# Time loop
    while (it < ntmax and time <= tmax):
        for i in reversed(range(nx)):
            u[i] = u[i] - c * dt / dx * (u[i] - u[(i-1)%nx])

        time += dt
        it += 1

        if cfl<1 :
        # Postprocessing
        line1.set_ydata(u)
        writer.grab_frame()
        if (it % output_freq == 0):
            print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
            np.savetxt( "TransportEquation_UpwindScheme_"+str(nx)+"Cells_ResultField_"+str(it)+".txt", u, delimiter="\n")
            plt.savefig("TransportEquation_UpwindScheme_"+str(nx)+"Cells_ResultField_"+str(it)+".png")
            #plt.show()
            pass
        pass

    assert max(u) <= max_initial
    assert min(u) >= min_initial
    assert np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]) <= total_var_initial

    print("Simulation of transport equation with upwind scheme done.")
    
    end = time.time()

    #return min, max, sol, total variation, l1 error and elapsed time
    return min(u), max(u), u, np.sum([abs(u[i] - u[(i-1)%nx]) for i in range(nx)]), dx*np.sum([abs(u[i] - u_initial[i]) for i in range(nx)]), end-start


if __name__ == """__main__""":
    if len(sys.argv) >2 :
        nx = int(sys.argv[1])
        cfl = float(sys.argv[2])
        solve(nx,cfl,0,1)
    else :
        nx = 50 # number of cells
        cfl = 0.99 # c*dt/dx <= CFL
        solve(nx,cfl,0,1)
    
