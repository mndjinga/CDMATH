#!/usr/bin/env python
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Résolution VF de l'équation du Burgers 1D \partial_t u + 1/2 \partial_x u^2 = 0 
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2019
# Description : Utilisation du schéma de Lax-Friedrichs explicite sur un maillage 1D régulier
#               Condition limite Dirichlet à gauche, Neumann à droite
#		        Création et sauvegarde du champ résultant et des figures
#               Génération d'une video sauvegardée dans un fichier .mp4
#================================================================================================================================


from math import sin, cos, pi
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import sys

def Burgers1D():
    ##################### Simulation parameters
    a = 0.0 # space domain :  a <= x <= b
    b = 1.3
    nx=1300
    dx = (b - a) / nx #space step

    tmax = 0.9 # runs the simulation for 0 <= t <= tMax
    cfl=0.4

    x=[a+0.5*dx + i*dx for i in range(nx)]   # array of cell center (1D mesh)
    
    u_bc_L=1.0
    u_bc_R=0.0

    ########################## Initial data
    xa=0.1
    xb=0.5
    u_initial = [ (xi<xa)+(xa<=xi)*(xi<=xb)*(np.cos(np.pi*(xa-xi)/(xa-xb))+1.0)*0.5  for xi in x];# to be used with a=0, b=1
    u         = [ (xi<xa)+(xa<=xi)*(xi<=xb)*(np.cos(np.pi*(xa-xi)/(xa-xb))+1.0)*0.5  for xi in x];# to be used with a=0, b=1

    max_initial=max(u_initial)
    min_initial=min(u_initial)

    time = 0.
    it = 0
    output_freq = 10

    # Video settings
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title="Lax-Friedrichs explicit scheme for Burgers equation", artist = "CEA Saclay", comment="Non linear wave propagation")
    writer=FFMpegWriter(fps=output_freq, metadata=metadata, codec='h264')
    with writer.saving(plt.figure(), "1DBurgersEquation_Lax-FriedrichsExplicit"+".mp4", ntmax):
        ########################### Postprocessing initialisation
        # Picture frame
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('u')
        plt.xlim(a,b)
        plt.ylim( min_initial - 0.1*(max_initial-min_initial), max_initial +  0.1*(max_initial-min_initial) )
        plt.title('Lax-Friedrichs explicit scheme for Burgers equation')
        line1, = plt.plot(x, u, label='u') #new picture for video # Returns a tuple of line objects, thus the comma
    
        print("Starting time loop")
        print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
        np.savetxt("BurgersEquation_Lax-FriedrichsExplicit_ResultField_0"+".txt", u, delimiter="\n")
        writer.grab_frame()
        plt.savefig("BurgersEquation_Lax-FriedrichsExplicit_ResultField_0"+".png")

        ############################# Time loop
        while (it < ntmax and time <= tmax):
            upw = max(abs(u))
            dt = cfl * dx / upw
            for i in reversed(range(nx)):
                u[i] = u[i] - c * dt / dx * (u[i] - u[(i-1)%nx])
    
            time += dt
            it += 1

            # Postprocessing
            line1.set_ydata(u)
            writer.grab_frame()
            if (it % output_freq == 0):
                print("-- Iter: " + str(it) + ", Time: " + str(time) + ", dt: " + str(dt))
                np.savetxt( "BurgersEquation_Lax-FriedrichsExplicit_ResultField_"+str(it)+".txt", u, delimiter="\n")
                plt.savefig("BurgersEquation_Lax-FriedrichsExplicit_ResultField_"+str(it)+".png")
                #plt.show()
    
    print("Simulation of Burgers equation with explicit Lax-Friedrichs scheme done.")
    
    return

if __name__ == """__main__""":
    Burgers1D()
