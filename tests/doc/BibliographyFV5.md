## Some bibliographical remarks

- Order 1 convergence on triangular meshes   
  *R. Herbin, An error estimate for a four point finite volume scheme for the convection-diffusion equation on a triangular mesh, Num. Meth. P.D.E., 165-173, 1995.*


- On triangular meshes, the FV5 scheme order is 2 provided 
    - the center of the circumscribed circle is used instead of the center of mass in each cell
    - the Delaunay conditions are satisfied (no neighboring cell is included in the circumscribed circle of an arbitrary cell)


- Non convergence on highly deformed meshes  
  *K. Domelevo, P. Omnes, A finite volume method for the Laplace equation on almost arbitrary 2D grids, Mathematical Modelling and Numerical Analysis, 2005*


- The scheme is order 1 if the mesh is conforming except on a line  
  *J. Droniou, C. Le Potier, Construction and Convergence Study of Schemes Preserving the Elliptic Local Maximum Principle, SIAM Journal on Numerical Analysis, 2011*


- *J. Droniou 2017* (ordre 2 si condition de type Delaunay)  
  A mixed finite element method for a sixth order elliptic problem  
  The gradient discretisation method for optimal control problems, with super-convergence for non-conforming finite elements and mixed-hybrid mimetic finite differences  
  Wsp-approximation properties of elliptic projectors on polynomial spaces, with application to the error analysis of a Hybrid High-Order discretisation of Leray-Lions problems  
  A Hybrid High-Order method for Leray-Lions elliptic equations on general meshes


- It is possible to converge with order 1 on the gradient, but only order 1 on the function ie there is no equivalent of the Aubin-Nitsche lemma in the finite volume context  
  *P. Omnes, Error estimates for a finite volume method for the Laplace equation in dimension one through discrete Green functions. International Journal on Finite Volumes 6(1), 18p., electronic only, 2009*
