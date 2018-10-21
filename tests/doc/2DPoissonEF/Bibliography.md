- The minimum angle condition (1968)  
  $\alpha_K\geq\alpha_0>0,\,\forall K\in \mathcal{T}$ where $\alpha_K$ is the minimal angle of $K$  
  is sufficient for the convergence of the quadratic finite element method
  $$
  \exists C,\quad || u - u_h||_1 \leq C \frac{h^2}{\sin(\alpha_0)}
  $$ 
  as $h\to 0$.

  It is equivalent to
  $$
  \exists C,\forall h\forall K\in\mathcal{T}_h,\quad |K|\geq C h^2
  $$
  *J. Brandts, S. Korotov, M. Křı́žek, On the equivalence of ball conditions for simplicial finite elements in R d , Appl. Math. Lett. 22 (2009), 1210–1212*

- The maximum angle condition (1976)  
  $\gamma_K\leq\gamma_0<\pi,\,\forall K\in \mathcal{T}$ where $\gamma_K$ is the maximal angle of $K$  
  is sufficient for the convergence of the quadratic finite element method
  $$
  \exists C,\quad || u - u_h||_1 \leq C h ||u||_2}
  $$ 
  as $h\to 0$.

- The maximum angle condition is not necessary for convergence  
  *A. Hannukainen, S. Korotov, M. Křı́žek, Maximum angle condition is not necessary for convergence of the finite element method,  M. Numer. Math. (2012) 120: 79*

- Nonobtuse simplicial partitions yield maximum principle thanks to a diagonally dominant stiffness matrices  
  *J. Karátson, S. Korotov, M. Křı́žek, On discrete maximum principles for nonlinear elliptic problems, Math. Comput. Simulation 76 (2007), 99–108*  
  *J. Brandts, S. Korotov, M. Křı́žek, The discrete maximum principle for linear simplicial finite element approximations of a reaction-diffusion problem, Linear Algebra Appl. 429 (2008), 2344–2357*
