/**
  @page step_72 The step-72 tutorial program  
* This tutorial depends on   step-71  ,   step-15  .
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Motivation">Motivation</a><a href="#Motivation">Motivation</a>
        <li><a href="#ComputingtheJacobianfromtheresidual"> Computing the Jacobian from the residual </a><a href="#ComputingtheJacobianfromtheresidual"> Computing the Jacobian from the residual </a>
        <li><a href="#ComputingtheJacobianandtheresidualfromtheenergyfunctional"> Computing the Jacobian and the residual from the energy functional </a><a href="#ComputingtheJacobianandtheresidualfromtheenergyfunctional"> Computing the Jacobian and the residual from the energy functional </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMinimalSurfaceProblemParameterscodeclass">The <code>MinimalSurfaceProblemParameters</code> class</a><a href="#ThecodeMinimalSurfaceProblemParameterscodeclass">The <code>MinimalSurfaceProblemParameters</code> class</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a>
        <li><a href="#Boundarycondition">Boundary condition</a><a href="#Boundarycondition">Boundary condition</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#MinimalSurfaceProblemMinimalSurfaceProblem">MinimalSurfaceProblem::MinimalSurfaceProblem</a><a href="#MinimalSurfaceProblemMinimalSurfaceProblem">MinimalSurfaceProblem::MinimalSurfaceProblem</a>
        <li><a href="#MinimalSurfaceProblemsetup_system">MinimalSurfaceProblem::setup_system</a><a href="#MinimalSurfaceProblemsetup_system">MinimalSurfaceProblem::setup_system</a>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
      <ul>
        <li><a href="#Manualassembly">Manual assembly</a><a href="#Manualassembly">Manual assembly</a>
        <li><a href="#Assemblyviadifferentiationoftheresidualvector">Assembly via differentiation of the residual vector</a><a href="#Assemblyviadifferentiationoftheresidualvector">Assembly via differentiation of the residual vector</a>
        <li><a href="#Assemblyviadifferentiationoftheenergyfunctional">Assembly via differentiation of the energy functional</a><a href="#Assemblyviadifferentiationoftheenergyfunctional">Assembly via differentiation of the energy functional</a>
      </ul>
        <li><a href="#MinimalSurfaceProblemsolve">MinimalSurfaceProblem::solve</a><a href="#MinimalSurfaceProblemsolve">MinimalSurfaceProblem::solve</a>
        <li><a href="#MinimalSurfaceProblemrefine_mesh">MinimalSurfaceProblem::refine_mesh</a><a href="#MinimalSurfaceProblemrefine_mesh">MinimalSurfaceProblem::refine_mesh</a>
        <li><a href="#MinimalSurfaceProblemset_boundary_values">MinimalSurfaceProblem::set_boundary_values</a><a href="#MinimalSurfaceProblemset_boundary_values">MinimalSurfaceProblem::set_boundary_values</a>
        <li><a href="#MinimalSurfaceProblemcompute_residual">MinimalSurfaceProblem::compute_residual</a><a href="#MinimalSurfaceProblemcompute_residual">MinimalSurfaceProblem::compute_residual</a>
        <li><a href="#MinimalSurfaceProblemdetermine_step_length">MinimalSurfaceProblem::determine_step_length</a><a href="#MinimalSurfaceProblemdetermine_step_length">MinimalSurfaceProblem::determine_step_length</a>
        <li><a href="#MinimalSurfaceProblemoutput_results">MinimalSurfaceProblem::output_results</a><a href="#MinimalSurfaceProblemoutput_results">MinimalSurfaceProblem::output_results</a>
        <li><a href="#MinimalSurfaceProblemrun">MinimalSurfaceProblem::run</a><a href="#MinimalSurfaceProblemrun">MinimalSurfaceProblem::run</a>
        <li><a href="#Themainfunction">The main function</a><a href="#Themainfunction">The main function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
*   <br>  
* <i>This program was contributed by Jean-Paul Pelteret and Wolfgang Bangerth.
* 
Wolfgang Bangerth's work is partially supported by National Science
Foundation grants OCI-1148116, OAC-1835673, DMS-1821210, and EAR-1925595;
and by the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-1550901 and The University of California-Davis.
</i>
* 

* <a name="Introduction"></a><h1>Introduction</h1>
* 

* <a name="Motivation"></a><h3>Motivation</h3>
* 

* This program solves the same problem as   step-15  , that is, it solvesfor the[minimal surface equation](https://en.wikipedia.org/wiki/Minimal_surface) 
* @f{align*}
    F(u) \dealcoloneq
* 
-\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.
  @f}
* 
* Among the issues we had identified there (see the<a href="step_15#extensions">Possibilities for extensions</a> section)was that when wanting to usea Newton iteration, we needed to compute the derivative of theresidual of the equation with regard to the solution   $u$   (here,because the right hand side is zero, the residual is simply the lefthand side). For the equation we have here, this is cumbersome but notimpossible
* 
*  -  but one can easily imagine much more complicatedequations where just implementing the residual itself correctly is achallenge, let alone doing so for the derivative necessary to computethe Jacobian matrix. We will address this issue in this program: Usingthe automatic differentiation techniques discussed in great detail in  step-71  , we will come up with a way how we only have to implement theresidual and get the Jacobian for free.
* In fact, we can even go one step further. While in   step-15   we havejust taken the equation as a given, the minimal surface equation isactually the product of minimizing an energy. Specifically,the minimal surface equations are the Euler-Lagrange equations thatcorrespond to minimizing the energy  @f[
    E(u) = \int_\Omega \Psi \left( u \right)
  @f]where theenergy density* is given by  @f[
    \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}}.
  @f]This is the same as saying that we seek to find the stationary point ofthe variation of the energy functional  @f[
    \min\limits_{u} E \left( u \right)
      \quad \rightarrow \quad
      \delta E \left( u, \varphi \right) \dealcoloneq
      \left(\varphi, F(u)\right) = 0
      \qquad
      \forall \varphi,
  @f]as this is where the equilibrium solution to the boundary value problem lies.
* The key point then is that, maybe, we don't even need to implement theresidual, but that implementing the simpler energy density   $\Psi(u)$  might actually be enough.
* Our goal then is this: Whenusing a Newton iteration, we need to repeatedly solve thelinear partial differential equation 
* @f{align*}
    F'(u^{n},\delta u^{n}) &=- F(u^{n})
  @f}
* so that we can compute the update 
* @f{align*}
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}
* with the solution   $\delta u^{n}$   of the Newton step. As discussed in   step-15  ,we can compute the derivative   $F'(u,\delta u)$   by hand andobtain  @f[
  F'(u,\delta u)
  =
* 


* 
* 


* 
* 
- \nabla \cdot \left( \frac{1}{\left(1+|\nabla u|^{2}\right)^{\frac{1}{2}}}\nabla
  \delta u \right) +
  \nabla \cdot \left( \frac{\nabla u \cdot
  \nabla \delta u}{\left(1+|\nabla u|^{2}\right)^{\frac{3}{2}}} \nabla u
  \right).
  @f]
* So here then is what this program is about: It is about techniquesthat can help us with computing   $F'(u,\delta u)$   without having toimplement it explicitly, either by providing an implementation of  $F(u)$   or an implementation of   $E(u)$  . More precisely, we willimplement three different approaches and compare them in terms ofrun-time but also
* 
*  -  maybe more importantly
* 
*  -  how much human effortit takes to implement them:
* 
*  - The method used in   step-15   to form the Jacobian matrix.
* 
*  - Computing the Jacobian matrix from an implementation of the  residual   $F(u)$  , using automatic differentiation.
* 
*  - Computing both the residual and Jacobian matrix from an  implementation of the energy functional   $E(u)$  , also using automatic  differentiation.
* For the first of these methods, there are no conceptual changescompared to   step-15  .
* 

* <a name="ComputingtheJacobianfromtheresidual"></a><h3> Computing the Jacobian from the residual </h3>
* 

* For the second method, let us outline how we will approach the issueusing automatic differentiationto compute the linearization of the residual vector. To this end, letus change notation for a moment and denote by   $F(U)$   not the residualof the differential equation, but in fact theresidual vector*
* 
*  - i.e., thediscrete residual*. We do so because that is what we
*actually* do when we discretize the problem on a given mesh: We solvethe problem   $F(U)=0$   where   $U$   is the vector of unknowns.
* More precisely, the   $i$  th component of the residual is given by@f[
  F(U)_i \dealcoloneq
  \int\limits_{\Omega}\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,
@f]where   $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$  . Given this, thecontribution for cell   $K$   is@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,
@f]Its first order Taylor expansion is given as@f[
  F(U + \delta U)_i^K
  \approx F(U)_i^K
  + \sum_{j}^{n_{\textrm{dofs}}} \left[ \frac{\partial F(U)_i^K}{\partial
  U_j} \delta U_j \right],
@f]and consequently we can compute the contribution of cell   $K$   to theJacobian matrix   $J$   as   $J(U)_{ij}^K = \frac{\partial F(U)_i^K}{\partial U_j}$  . Theimportant point here is that on cell   $K$  , we can express@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+\left|
  \sum_{j'}^{n_\textrm{dofs}} U_{j'} \nabla \varphi_{j'}\right|^{2}}}
  \left(\sum_{j''}^{n_\textrm{dofs}} U_{j''} \nabla \varphi_{j''}\right)\right] \, dV.
@f]For clarity, we have used   $j'$   and   $j''$   as counting indices to makeclear that they are distinct from each other and from   $j$   above.Because in this formula,   $F(U)$   only depends on the coefficients  $U_j$  , we can compute the derivative   $J(U)_{ij}^K$   as a matrix viaautomatic differentiation of   $F(U)_i^K$  . By the same argument as wealways use, it is clear that   $F(U)^K$   does not actually depend on
*all* unknowns   $U_j$  , but only on those unknowns for which   $j$   is ashape function that lives on cell   $K$  , and so in practice, we restrict  $F(U)^K$   and   $J(U)^K$   to that part of the vector and matrix thatcorresponds to thelocal* DoF indices, and then distribute from thelocal cell   $K$   to the global objects.
* Using all of these realizations, the approach will then be toimplement   $F(U)^K$   in the program and let the automatic differentiationmachinery compute the derivatives   $J(U)^K$   from that.
* 

* <a name="ComputingtheJacobianandtheresidualfromtheenergyfunctional"></a><h3> Computing the Jacobian and the residual from the energy functional </h3>
* 

* For the final implementation of the assembly process, we will move a levelhigher than the residual: our entire linear system will be determineddirectly from the energy functional that governs the physics of thisboundary value problem. We can take advantage of the fact that we cancalculate the total energy in the domain directly from the localcontributions, i.e.,@f[
  E \left( U \right) \dealcoloneq \int\limits_{\Omega} \Psi \left( u
  \right) \, dV .
@f]In the discrete setting, this means that on each finite element we have@f[
   E \left( U \right)^K
    \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
    \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
    \mathbf{x}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times W_{q}}_{\text{JxW(q)}} .
@f]If we implement the cell energy, which depends on the field solution,we can compute its first (discrete) variation@f[
  F(U)^K_i
    = \frac{\partial E(U)^K}{\partial U_i}
@f]and, thereafter, its second (discrete) variation@f[
  J(U)^K_{ij}
    = \frac{\partial^{2}  E(U)^K}{\partial U_i \partial U_j}.
@f]So, from the cell contribution to the total energy function, we may expectto have the approximate residual and tangent contributions generatedfor us as long as we can provide an implementation of the local energy  $E(U)^K$  . Again, due to the design of theautomatic differentiation variables used in this tutorial, in practicethese approximations for the contributions to the residual vector andtangent matrix are actually accurate to machine precision.
* 

*  <a name="CommProg"></a> <h1> The commented program</h1>
*  The majority of this tutorial is an exact replica of   step-15  . So, in the interest of brevity and maintaining a focus on the changes implemented here, we will only document what's new and simply indicate which sections of code are a repetition of what has come before.
* 

* 
*  
*  
*  <a name="Includefiles"></a>  <h3>Include files</h3>
* 

* 
*  There are a few new header files that have been included in this tutorial. The first is the one that provides the declaration of the ParameterAcceptor class.
* 

* 
* @code
 #include <deal.II/base/quadrature_lib.h>
 #include <deal.II/base/function.h>
 #include <deal.II/base/parameter_acceptor.h>
 #include <deal.II/base/timer.h>
 #include <deal.II/base/utilities.h>
* 
 @endcode
* 
*  This is the second, which is an all-inclusive header that will allow us to incorporate the automatic differentiation (AD) functionality within this code.
* 

* 
* @code
 #include <deal.II/differentiation/ad.h>
* 
 #include <deal.II/lac/vector.h>
 #include <deal.II/lac/full_matrix.h>
 #include <deal.II/lac/sparse_matrix.h>
 #include <deal.II/lac/dynamic_sparsity_pattern.h>
 #include <deal.II/lac/solver_cg.h>
 #include <deal.II/lac/precondition.h>
 #include <deal.II/lac/affine_constraints.h>
* 
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/grid_refinement.h>
* 
 #include <deal.II/dofs/dof_handler.h>
 #include <deal.II/dofs/dof_tools.h>
* 
 #include <deal.II/fe/fe_values.h>
 #include <deal.II/fe/fe_values_extractors.h>
 #include <deal.II/fe/fe_q.h>
* 
 @endcode
* 
*  And the next three provide some multi-threading capability using the generic   MeshWorker::mesh_loop()   framework.
* 

* 
* @code
 #include <deal.II/meshworker/copy_data.h>
 #include <deal.II/meshworker/mesh_loop.h>
 #include <deal.II/meshworker/scratch_data.h>
* 
 #include <deal.II/numerics/vector_tools.h>
 #include <deal.II/numerics/matrix_tools.h>
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/numerics/error_estimator.h>
* 
 
 #include <fstream>
 #include <iostream>
* 
 #include <deal.II/numerics/solution_transfer.h>
* 
 @endcode
* 
*  We then open a namespace for this program and import everything from the dealii namespace into it, as in previous programs:
* 

* 
* @code
 namespace Step72
 {
   using namespace dealii;
* 
 @endcode
* 
*   <a name="ThecodeMinimalSurfaceProblemParameterscodeclass"></a>  <h3>The <code>MinimalSurfaceProblemParameters</code> class</h3>
* 

* 
*  In this tutorial we will implement three different approaches for assembling the linear system. One mirrors the hand implementation originally provided in   step-15  , while the other two use the Sacado automatic differentiation library that is provided as a part of the Trilinos framework.   
*   To facilitate switching between the three implementations, we have this really basic parameters class that has only two options that are configurable.
* 

* 
* @code
   class MinimalSurfaceProblemParameters : public ParameterAcceptor
   {
   public:
     MinimalSurfaceProblemParameters();
* 
 @endcode
* 
*  Selection for the formulation and corresponding AD framework to be used:
* 

* 
* 
*  -  formulation = 0 : Unassisted implementation (full hand linearization).
* 

* 
* 
*  -  formulation = 1 : Automated linearization of the finite element residual.
* 

* 
* 
*  -  formulation = 2 : Automated computation of finite element residual and linearization using a variational formulation.
* 

* 
* @code
     unsigned int formulation = 0;
* 
 @endcode
* 
*  The maximum acceptable tolerance for the linear system residual. We will see that the assembly time becomes appreciable once we use the AD framework, so we have increased the tolerance selected in   step-15   by one order of magnitude. This way, the computations do not take too long to complete.
* 

* 
* @code
     double tolerance = 1e-2;
   };
* 
 
   MinimalSurfaceProblemParameters::MinimalSurfaceProblemParameters()
     : ParameterAcceptor("Minimal Surface Problem/")
   {
     add_parameter(
       "Formulation", formulation, "", this->prm, Patterns::Integer(0, 2));
     add_parameter("Tolerance", tolerance, "", this->prm, Patterns::Double(0.0));
   }
* 
 
* 
 @endcode
* 
*   <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a>  <h3>The <code>MinimalSurfaceProblem</code> class template</h3>
* 

* 
*  The class template is essentially the same as in   step-15  . The only functional changes to the class are that:
* 

* 
* 
*  - the run() function now takes in two arguments: one to choose which assembly approach is to be adopted, and one for the tolerance for the permissible final residual is, and
* 

* 
* 
*  - there are now three different assembly functions that implement the three methods of assembling the linear system. We'll provide details on these later on.
* 

* 
*  

* 
* @code
   template <int dim>
   class MinimalSurfaceProblem
   {
   public:
     MinimalSurfaceProblem();
* 
     void run(const int formulation, const double tolerance);
* 
   private:
     void   setup_system(const bool initial_step);
     void   assemble_system_unassisted();
     void   assemble_system_with_residual_linearization();
     void   assemble_system_using_energy_functional();
     void   solve();
     void   refine_mesh();
     void   set_boundary_values();
     double compute_residual(const double alpha) const;
     double determine_step_length() const;
     void   output_results(const unsigned int refinement_cycle) const;
* 
     Triangulation<dim> triangulation;
* 
     DoFHandler<dim> dof_handler;
     FE_Q<dim>       fe;
     QGauss<dim>     quadrature_formula;
* 
     AffineConstraints<double> hanging_node_constraints;
* 
     SparsityPattern      sparsity_pattern;
     SparseMatrix<double> system_matrix;
* 
     Vector<double> current_solution;
     Vector<double> newton_update;
     Vector<double> system_rhs;
   };
* 
 @endcode
* 
*   <a name="Boundarycondition"></a>  <h3>Boundary condition</h3>
* 

* 
*  There are no changes to the boundary conditions applied to the problem.
* 

* 
* @code
   template <int dim>
   class BoundaryValues : public Function<dim>
   {
   public:
     virtual double value(const Point<dim> & p,
                          const unsigned int component = 0) const override;
   };
* 
 
   template <int dim>
   double BoundaryValues<dim>::value(const Point<dim> &p,
                                     const unsigned int  /*component*/ ) const
   {
     return std::sin(2 numbers::PI (p[0] + p[1]));
   }
* 
 
 @endcode
* 
*   <a name="ThecodeMinimalSurfaceProblemcodeclassimplementation"></a>  <h3>The <code>MinimalSurfaceProblem</code> class implementation</h3>
* 

* 
*   <a name="MinimalSurfaceProblemMinimalSurfaceProblem"></a>  <h4>MinimalSurfaceProblem::MinimalSurfaceProblem</h4>
* 

* 
*  There have been no changes made to the class constructor.
* 

* 
* @code
   template <int dim>
   MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
     : dof_handler(triangulation)
     , fe(2)
     , quadrature_formula(fe.degree + 1)
   {}
* 
 
 @endcode
* 
*   <a name="MinimalSurfaceProblemsetup_system"></a>  <h4>MinimalSurfaceProblem::setup_system</h4>
* 

* 
*  There have been no changes made to the function that sets up the class data structures, namely the DoFHandler, the hanging node constraints applied to the problem, and the linear system.
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
   {
     if (initial_step)
       {
         dof_handler.distribute_dofs(fe);
         current_solution.reinit(dof_handler.n_dofs());
* 
         hanging_node_constraints.clear();
         DoFTools::make_hanging_node_constraints(dof_handler,
                                                 hanging_node_constraints);
         hanging_node_constraints.close();
       }
* 
     newton_update.reinit(dof_handler.n_dofs());
     system_rhs.reinit(dof_handler.n_dofs());
* 
     DynamicSparsityPattern dsp(dof_handler.n_dofs());
     DoFTools::make_sparsity_pattern(dof_handler, dsp);
* 
     hanging_node_constraints.condense(dsp);
* 
     sparsity_pattern.copy_from(dsp);
     system_matrix.reinit(sparsity_pattern);
   }
* 
 @endcode
* 
*   <a name="Assemblingthelinearsystem"></a>  <h4>Assembling the linear system</h4>
* 

* 
*   <a name="Manualassembly"></a>  <h5>Manual assembly</h5>
* 

* 
*  The assembly functions are the interesting contributions to this tutorial. The assemble_system_unassisted() method implements exactly the same assembly function as is detailed in   step-15  , but in this instance we use the   MeshWorker::mesh_loop()   function to multithread the assembly process. The reason for doing this is quite simple: When using automatic differentiation, we know that there is to be some additional computational overhead incurred. In order to mitigate this performance loss, we'd like to take advantage of as many (easily available) computational resources as possible. The   MeshWorker::mesh_loop()   concept makes this a relatively straightforward task. At the same time, for the purposes of fair comparison, we need to do the same to the implementation that uses no assistance when computing the residual or its linearization. (The   MeshWorker::mesh_loop()   function is first discussed in   step-12   and   step-16  , if you'd like to read up on it.)   
*   The steps required to implement the multithreading are the same between the three functions, so we'll use the assemble_system_unassisted() function as an opportunity to focus on the multithreading itself.
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::assemble_system_unassisted()
   {
     system_matrix = 0;
     system_rhs    = 0;
* 
     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
* 
 @endcode
* 
*  The   MeshWorker::mesh_loop()   expects that we provide two exemplar data structures. The first, `ScratchData`, is to store all large data that is to be reused between threads. The `CopyData` will hold the contributions to the linear system that come from each cell. These independent matrix-vector pairs must be accumulated into the global linear system sequentially. Since we don't need anything on top of what the   MeshWorker::ScratchData   and   MeshWorker::CopyData   classes already provide, we use these exact class definitions for our problem. Note that we only require a single instance of a local matrix, local right-hand side vector, and cell degree of freedom index vector
* 
*  -  the   MeshWorker::CopyData   therefore has `1` for all three of its template arguments.
* 

* 
* @code
     using ScratchData = MeshWorker::ScratchData<dim>;
     using CopyData    = MeshWorker::CopyData<1, 1, 1>;
* 
 @endcode
* 
*  We also need to know what type of iterator we'll be working with during assembly. For simplicity, we just ask the compiler to work this out for us using the decltype() specifier, knowing that we'll be iterating over active cells owned by the   @p dof_handler.  
* 

* 
* @code
     using CellIteratorType = decltype(dof_handler.begin_active());
* 
 @endcode
* 
*  Here we initialize the exemplar data structures. Since we know that we need to compute the shape function gradients, weighted Jacobian, and the position of the quadrate points in real space, we pass these flags into the class constructor.
* 

* 
* @code
     const ScratchData sample_scratch_data(fe,
                                           quadrature_formula,
                                           update_gradients |
                                             update_quadrature_points |
                                             update_JxW_values);
     const CopyData    sample_copy_data(dofs_per_cell);
* 
 @endcode
* 
*  Now we define a lambda function that will perform the assembly on a single cell. The three arguments are those that will be expected by   MeshWorker::mesh_loop(),   due to the arguments that we'll pass to that final call. We also capture the   @p this   pointer, which means that we'll have access to "this" (i.e., the current `MinimalSurfaceProblem<dim>`) class instance, and its private member data (since the lambda function is defined within a MinimalSurfaceProblem<dim> method).     
*   At the top of the function, we initialize the data structures that are dependent on the cell for which the work is being performed. Observe that the reinitialization call actually returns an instance to an FEValues object that is initialized and stored within (and, therefore, reused by) the `scratch_data` object.     
*   Similarly, we get aliases to the local matrix, local RHS vector, and local cell DoF indices from the `copy_data` instance that   MeshWorker::mesh_loop()   provides. We then initialize the cell DoF indices, knowing that the local matrix and vector are already correctly sized.
* 

* 
* @code
     const auto cell_worker = [this](const CellIteratorType &cell,
                                     ScratchData &           scratch_data,
                                     CopyData &              copy_data) {
       const auto &fe_values = scratch_data.reinit(cell);
* 
       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
       std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
       cell->get_dof_indices(local_dof_indices);
* 
 @endcode
* 
*  For Newton's method, we require the gradient of the solution at the point about which the problem is being linearized.       
*   Once we have that, we can perform assembly for this cell in the usual way.  One minor difference to   step-15   is that we've used the (rather convenient) range-based loops to iterate over all quadrature points and degrees-of-freedom.
* 

* 
* @code
       std::vector<Tensor<1, dim>> old_solution_gradients(
         fe_values.n_quadrature_points);
       fe_values.get_function_gradients(current_solution,
                                        old_solution_gradients);
* 
       for (const unsigned int q : fe_values.quadrature_point_indices())
         {
           const double coeff =
             1.0 / std::sqrt(1.0 + old_solution_gradients[q]
                                     old_solution_gradients[q]);
* 
           for (const unsigned int i : fe_values.dof_indices())
             {
               for (const unsigned int j : fe_values.dof_indices())
                 cell_matrix(i, j) +=
                   (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
                      coeff                         //   a_n
                      fe_values.shape_grad(j, q))   //   \nabla \phi_j)
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 
-                                //
* 


* 
* 
-
                     (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
                      coeff coeff coeff         //   a_n^3
                      (fe_values.shape_grad(j, q)   //   (\nabla \phi_j
                         old_solution_gradients[q]) //      \nabla u_n)
                      old_solution_gradients[q]))   //   \nabla u_n)))
                    fe_values.JxW(q));              // dx
* 
               cell_rhs(i)
* 
-= (fe_values.shape_grad(i, q)  // \nabla \phi_i
                               coeff                     // a_n
                               old_solution_gradients[q] // u_n
                               fe_values.JxW(q));        // dx
             }
         }
     };
* 
 @endcode
* 
*  The second lambda function that   MeshWorker::mesh_loop()   requires is one that performs the task of accumulating the local contributions in the global linear system. That is precisely what this one does, and the details of the implementation have been seen before. The primary point to recognize is that the local contributions are stored in the `copy_data` instance that is passed into this function. This `copy_data` has been filled with data during   @a   some call to the `cell_worker`.
* 

* 
* @code
     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
       const std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
* 
       for (unsigned int i = 0; i < dofs_per_cell; ++i)
         {
           for (unsigned int j = 0; j < dofs_per_cell; ++j)
             system_matrix.add(local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i, j));
* 
           system_rhs(local_dof_indices[i]) += cell_rhs(i);
         }
     };
* 
 @endcode
* 
*  We have all of the required functions definitions in place, so now we call the   MeshWorker::mesh_loop()   to perform the actual assembly.  We pass a flag as the last parameter which states that we only want to perform the assembly on the cells. Internally,   MeshWorker::mesh_loop()   then distributes the available work to different threads, making efficient use of the multiple cores almost all of today's processors have to offer.
* 

* 
* @code
     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                           cell_worker,
                           copier,
                           sample_scratch_data,
                           sample_copy_data,
                           MeshWorker::assemble_own_cells);
* 
 @endcode
* 
*  And finally, as is done in   step-15  , we remove hanging nodes from the system and apply zero boundary values to the linear system that defines the Newton updates   $\delta u^n$  .
* 

* 
* @code
     hanging_node_constraints.condense(system_matrix);
     hanging_node_constraints.condense(system_rhs);
* 
     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              0,
                                              Functions::ZeroFunction<dim>(),
                                              boundary_values);
     MatrixTools::apply_boundary_values(boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
   }
* 
 @endcode
* 
*   <a name="Assemblyviadifferentiationoftheresidualvector"></a>  <h5>Assembly via differentiation of the residual vector</h5>
* 

* 
*  As outlined in the introduction, what we need to do for this second approach is implement the local contributions   $F(U)^K$   from cell   $K$   to the residual vector, and then let the AD machinery deal with how to compute the derivatives   $J(U)_{ij}^K=\frac{\partial F(U)^K_i}{\partial U_j}$   from it.   
*   For the following, recall that @f[
 F(U)_i^K \dealcoloneq
 \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
 u|^{2}}} \nabla u \right] \, dV ,
 @f] where   $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$  .   
*   Let us see how this is implemented in practice:
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::assemble_system_with_residual_linearization()
   {
     system_matrix = 0;
     system_rhs    = 0;
* 
     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
* 
     using ScratchData      = MeshWorker::ScratchData<dim>;
     using CopyData         = MeshWorker::CopyData<1, 1, 1>;
     using CellIteratorType = decltype(dof_handler.begin_active());
* 
     const ScratchData sample_scratch_data(fe,
                                           quadrature_formula,
                                           update_gradients |
                                             update_quadrature_points |
                                             update_JxW_values);
     const CopyData    sample_copy_data(dofs_per_cell);
* 
 @endcode
* 
*  We'll define up front the AD data structures that we'll be using, utilizing the techniques shown in   step-71  . In this case, we choose the helper class that will automatically compute the linearization of the finite element residual using Sacado forward automatic differentiation types. These number types can be used to compute first derivatives only. This is exactly what we want, because we know that we'll only be linearizing the residual, which means that we only need to compute first-order derivatives. The return values from the calculations are to be of type `double`.     
*   We also need an extractor to retrieve some data related to the field solution to the problem.
* 

* 
* @code
     using ADHelper = Differentiation::AD::ResidualLinearization<
       Differentiation::AD::NumberTypes::sacado_dfad,
       double>;
     using ADNumberType = typename ADHelper::ad_type;
* 
     const FEValuesExtractors::Scalar u_fe(0);
* 
 @endcode
* 
*  With this, let us define the lambda function that will be used to compute the cell contributions to the Jacobian matrix and the right hand side:
* 

* 
* @code
     const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                            ScratchData &           scratch_data,
                                            CopyData &              copy_data) {
       const auto &       fe_values     = scratch_data.reinit(cell);
       const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();
* 
       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
       std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
       cell->get_dof_indices(local_dof_indices);
* 
 @endcode
* 
*  We'll now create and initialize an instance of the AD helper class. To do this, we need to specify how many independent variables and dependent variables there are. The independent variables will be the number of local degrees of freedom that our solution vector has, i.e., the number   $j$   in the per-element representation of the discretized solution vector   $u (\mathbf{x})|_K = \sum\limits_{j} U^K_i \varphi_j(\mathbf{x})$   that indicates how many solution coefficients are associated with each finite element. In deal.II, this equals   FiniteElement::dofs_per_cell.   The number of dependent variables will be the number of entries in the local residual vector that we will be forming. In this particular problem (like many others that employ the [standard Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method)) the number of local solution coefficients matches the number of local residual equations.
* 

* 
* @code
       const unsigned int n_independent_variables = local_dof_indices.size();
       const unsigned int n_dependent_variables   = dofs_per_cell;
       ADHelper ad_helper(n_independent_variables, n_dependent_variables);
* 
 @endcode
* 
*  Next we inform the helper of the values of the solution, i.e., the actual values for   $U_j$   about which we wish to linearize. As this is done on each element individually, we have to extract the solution coefficients from the global solution vector. In other words, we define all of those coefficients   $U_j$   where   $j$   is a local degree of freedom as the independent variables that enter the computation of the vector   $F(U)^{K}$   (the dependent function).       
*   Then we get the complete set of degree of freedom values as represented by auto-differentiable numbers. The operations performed with these variables are tracked by the AD library from this point until the object goes out of scope. So it is   <em>  precisely these variables  </em>   with respect to which we will compute derivatives of the residual entries.
* 

* 
* @code
       ad_helper.register_dof_values(current_solution, local_dof_indices);
* 
       const std::vector<ADNumberType> &dof_values_ad =
         ad_helper.get_sensitive_dof_values();
* 
 @endcode
* 
*  Then we do some problem specific tasks, the first being to compute all values, (spatial) gradients, and the like based on "sensitive" AD degree of freedom values. In this instance we are retrieving the solution gradients at each quadrature point. Observe that the solution gradients are now sensitive to the values of the degrees of freedom as they use the   @p ADNumberType   as the scalar type and the   @p dof_values_ad   vector provides the local DoF values.
* 

* 
* @code
       std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
         fe_values.n_quadrature_points);
       fe_values[u_fe].get_function_gradients_from_local_dof_values(
         dof_values_ad, old_solution_gradients);
* 
 @endcode
* 
*  The next variable that we declare will store the cell residual vector contributions. This is rather self-explanatory, save for one <b>very important</b> detail: Note that each entry in the vector is hand-initialized with a value of zero. This is a   <em>  highly recommended  </em>   practice, as some AD libraries appear not to safely initialize the internal data structures of these number types. Not doing so could lead to some very hard to understand or detect bugs (appreciate that the author of this program mentions this out of, generally bad, experience). So out of an abundance of caution it's worthwhile zeroing the initial value explicitly. After that, apart from a sign change the residual assembly looks much the same as we saw for the cell RHS vector before: We loop over all quadrature points, ensure that the coefficient now encodes its dependence on the (sensitive) finite element DoF values by using the correct `ADNumberType`, and finally we assemble the components of the residual vector. For complete clarity, the finite element shape functions (and their gradients, etc.) as well as the "JxW" values remain scalar valued, but the   @p coeff   and the    @p old_solution_gradients   at each quadrature point are computed in terms of the independent variables.
* 

* 
* @code
       std::vector<ADNumberType> residual_ad(n_dependent_variables,
                                             ADNumberType(0.0));
       for (const unsigned int q : fe_values.quadrature_point_indices())
         {
           const ADNumberType coeff =
             1.0 / std::sqrt(1.0 + old_solution_gradients[q]
                                     old_solution_gradients[q]);
* 
           for (const unsigned int i : fe_values.dof_indices())
             {
               residual_ad[i] += (fe_values.shape_grad(i, q)   // \nabla \phi_i
                                  coeff                      // a_n
                                  old_solution_gradients[q]) // u_n
                                 fe_values.JxW(q);           // dx
             }
         }
* 
 @endcode
* 
*  Once we have the full cell residual vector computed, we can register it with the helper class.       
*   Thereafter, we compute the residual values (basically, extracting the real values from what we already computed) and their Jacobian (the linearization of each residual component with respect to all cell DoFs) at the evaluation point. For the purposes of assembly into the global linear system, we have to respect the sign difference between the residual and the RHS contribution: For Newton's method, the right hand side vector needs to be equal to thenegative* residual vector.
* 

* 
* @code
       ad_helper.register_residual_vector(residual_ad);
* 
       ad_helper.compute_residual(cell_rhs);
       cell_rhs=
* 
-1.0;
* 
       ad_helper.compute_linearization(cell_matrix);
     };
* 
 @endcode
* 
*  The remainder of the function equals what we had previously:
* 

* 
* @code
     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
       const std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
* 
       for (unsigned int i = 0; i < dofs_per_cell; ++i)
         {
           for (unsigned int j = 0; j < dofs_per_cell; ++j)
             system_matrix.add(local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i, j));
* 
           system_rhs(local_dof_indices[i]) += cell_rhs(i);
         }
     };
* 
     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                           cell_worker,
                           copier,
                           sample_scratch_data,
                           sample_copy_data,
                           MeshWorker::assemble_own_cells);
* 
     hanging_node_constraints.condense(system_matrix);
     hanging_node_constraints.condense(system_rhs);
* 
     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              0,
                                              Functions::ZeroFunction<dim>(),
                                              boundary_values);
     MatrixTools::apply_boundary_values(boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
   }
* 
 @endcode
* 
*   <a name="Assemblyviadifferentiationoftheenergyfunctional"></a>  <h5>Assembly via differentiation of the energy functional</h5>
* 

* 
*  In this third approach, we compute residual and Jacobian as first and second derivatives of the local energy functional @f[
 E\left( U \right)^K
 \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
 \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
 \mathbf{X}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times
 W_{q}}_{\text{JxW(q)}}
 @f] with the energy density given by @f[
 \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}} .
 @f]   
*   Let us again see how this is done:
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::assemble_system_using_energy_functional()
   {
     system_matrix = 0;
     system_rhs    = 0;
* 
     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
* 
     using ScratchData      = MeshWorker::ScratchData<dim>;
     using CopyData         = MeshWorker::CopyData<1, 1, 1>;
     using CellIteratorType = decltype(dof_handler.begin_active());
* 
     const ScratchData sample_scratch_data(fe,
                                           quadrature_formula,
                                           update_gradients |
                                             update_quadrature_points |
                                             update_JxW_values);
     const CopyData    sample_copy_data(dofs_per_cell);
* 
 @endcode
* 
*  In this implementation of the assembly process, we choose the helper class that will automatically compute both the residual and its linearization from the cell contribution to an energy functional using nested Sacado forward automatic differentiation types. The selected number types can be used to compute both first and second derivatives. We require this, as the residual defined as the sensitivity of the potential energy with respect to the DoF values (i.e. its gradient). We'll then need to linearize the residual, implying that second derivatives of the potential energy must be computed. You might want to compare this with the definition of `ADHelper` used int previous function, where we used   `Differentiation::AD::ResidualLinearization<Differentiation::AD::NumberTypes::sacado_dfad,double>`.  
* 

* 
* @code
     using ADHelper = Differentiation::AD::EnergyFunctional<
       Differentiation::AD::NumberTypes::sacado_dfad_dfad,
       double>;
     using ADNumberType = typename ADHelper::ad_type;
* 
     const FEValuesExtractors::Scalar u_fe(0);
* 
 @endcode
* 
*  Let us then again define the lambda function that does the integration on a cell.     
*   To initialize an instance of the helper class, we now only require that the number of independent variables (that is, the number of degrees of freedom associated with the element solution vector) are known up front. This is because the second-derivative matrix that results from an energy functional is necessarily square (and also, incidentally, symmetric).
* 

* 
* @code
     const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
                                            ScratchData &           scratch_data,
                                            CopyData &              copy_data) {
       const auto &fe_values = scratch_data.reinit(cell);
* 
       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
       std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
       cell->get_dof_indices(local_dof_indices);
* 
       const unsigned int n_independent_variables = local_dof_indices.size();
       ADHelper           ad_helper(n_independent_variables);
* 
 @endcode
* 
*  Once more, we register all cell DoFs values with the helper, followed by extracting the "sensitive" variant of these values that are to be used in subsequent operations that must be differentiated
* 
*  -  one of those being the calculation of the solution gradients.
* 

* 
* @code
       ad_helper.register_dof_values(current_solution, local_dof_indices);
* 
       const std::vector<ADNumberType> &dof_values_ad =
         ad_helper.get_sensitive_dof_values();
* 
       std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
         fe_values.n_quadrature_points);
       fe_values[u_fe].get_function_gradients_from_local_dof_values(
         dof_values_ad, old_solution_gradients);
* 
 @endcode
* 
*  We next create a variable that stores the cell total energy. Once more we emphasize that we explicitly zero-initialize this value, thereby ensuring the integrity of the data for this starting value.       
*   The aim for our approach is then to compute the cell total energy, which is the sum of the internal (due to right hand side functions, typically linear in   $U$  ) and external energies. In this particular case, we have no external energies (e.g., from source terms or Neumann boundary conditions), so we'll focus on the internal energy part.       
*   In fact, computing   $E(U)^K$   is almost trivial, requiring only the following lines:
* 

* 
* @code
       ADNumberType energy_ad = ADNumberType(0.0);
       for (const unsigned int q : fe_values.quadrature_point_indices())
         {
           const ADNumberType psi = std::sqrt(1.0 + old_solution_gradients[q]
                                                      old_solution_gradients[q]);
* 
           energy_ad += psi fe_values.JxW(q);
         }
* 
 @endcode
* 
*  After we've computed the total energy on this cell, we'll register it with the helper.  Based on that, we may now compute the desired quantities, namely the residual values and their Jacobian at the evaluation point. As before, the Newton right hand side needs to be the negative of the residual:
* 

* 
* @code
       ad_helper.register_energy_functional(energy_ad);
* 
       ad_helper.compute_residual(cell_rhs);
       cell_rhs=
* 
-1.0;
* 
       ad_helper.compute_linearization(cell_matrix);
     };
* 
 @endcode
* 
*  As in the previous two functions, the remainder of the function is as before:
* 

* 
* @code
     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
       const std::vector<types::global_dof_index> &local_dof_indices =
         copy_data.local_dof_indices[0];
* 
       for (unsigned int i = 0; i < dofs_per_cell; ++i)
         {
           for (unsigned int j = 0; j < dofs_per_cell; ++j)
             system_matrix.add(local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i, j));
* 
           system_rhs(local_dof_indices[i]) += cell_rhs(i);
         }
     };
* 
     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
                           cell_worker,
                           copier,
                           sample_scratch_data,
                           sample_copy_data,
                           MeshWorker::assemble_own_cells);
* 
     hanging_node_constraints.condense(system_matrix);
     hanging_node_constraints.condense(system_rhs);
* 
     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              0,
                                              Functions::ZeroFunction<dim>(),
                                              boundary_values);
     MatrixTools::apply_boundary_values(boundary_values,
                                        system_matrix,
                                        newton_update,
                                        system_rhs);
   }
* 
 
 @endcode
* 
*   <a name="MinimalSurfaceProblemsolve"></a>  <h4>MinimalSurfaceProblem::solve</h4>
* 

* 
*  The solve function is the same as is used in   step-15  .
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::solve()
   {
     SolverControl            solver_control(system_rhs.size(),
                                  system_rhs.l2_norm() 1e-6);
     SolverCG<Vector<double>> solver(solver_control);
* 
     PreconditionSSOR<SparseMatrix<double>> preconditioner;
     preconditioner.initialize(system_matrix, 1.2);
* 
     solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
* 
     hanging_node_constraints.distribute(newton_update);
* 
     const double alpha = determine_step_length();
     current_solution.add(alpha, newton_update);
   }
* 
 
 @endcode
* 
*   <a name="MinimalSurfaceProblemrefine_mesh"></a>  <h4>MinimalSurfaceProblem::refine_mesh</h4>
* 

* 
*  Nothing has changed since   step-15   with respect to the mesh refinement procedure and transfer of the solution between adapted meshes.
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::refine_mesh()
   {
     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
* 
     KellyErrorEstimator<dim>::estimate(
       dof_handler,
       QGauss<dim
* 
- 1>(fe.degree + 1),
       std::map<types::boundary_id, const Function<dim>>(),
       current_solution,
       estimated_error_per_cell);
* 
     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                     estimated_error_per_cell,
                                                     0.3,
                                                     0.03);
* 
     triangulation.prepare_coarsening_and_refinement();
     SolutionTransfer<dim> solution_transfer(dof_handler);
     solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
     triangulation.execute_coarsening_and_refinement();
* 
     dof_handler.distribute_dofs(fe);
* 
     Vector<double> tmp(dof_handler.n_dofs());
     solution_transfer.interpolate(current_solution, tmp);
     current_solution = tmp;
* 
     hanging_node_constraints.clear();
     DoFTools::make_hanging_node_constraints(dof_handler,
                                             hanging_node_constraints);
     hanging_node_constraints.close();
* 
     set_boundary_values();
* 
     setup_system(false);
   }
* 
 
* 
 @endcode
* 
*   <a name="MinimalSurfaceProblemset_boundary_values"></a>  <h4>MinimalSurfaceProblem::set_boundary_values</h4>
* 

* 
*  The choice of boundary conditions remains identical to   step-15  ...
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::set_boundary_values()
   {
     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              0,
                                              BoundaryValues<dim>(),
                                              boundary_values);
     for (auto &boundary_value : boundary_values)
       current_solution(boundary_value.first) = boundary_value.second;
* 
     hanging_node_constraints.distribute(current_solution);
   }
* 
 
 @endcode
* 
*   <a name="MinimalSurfaceProblemcompute_residual"></a>  <h4>MinimalSurfaceProblem::compute_residual</h4>
* 

* 
*  ... as does the function used to compute the residual during the solution iteration procedure. One could replace this by differentiation of the energy functional if one really wanted, but for simplicity we here simply copy what we already had in   step-15  .
* 

* 
* @code
   template <int dim>
   double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
   {
     Vector<double> residual(dof_handler.n_dofs());
* 
     Vector<double> evaluation_point(dof_handler.n_dofs());
     evaluation_point = current_solution;
     evaluation_point.add(alpha, newton_update);
* 
     const QGauss<dim> quadrature_formula(fe.degree + 1);
     FEValues<dim>     fe_values(fe,
                             quadrature_formula,
                             update_gradients | update_quadrature_points |
                               update_JxW_values);
* 
     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
     const unsigned int n_q_points    = quadrature_formula.size();
* 
     Vector<double>              cell_residual(dofs_per_cell);
     std::vector<Tensor<1, dim>> gradients(n_q_points);
* 
     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
* 
     for (const auto &cell : dof_handler.active_cell_iterators())
       {
         cell_residual = 0;
         fe_values.reinit(cell);
* 
         fe_values.get_function_gradients(evaluation_point, gradients);
* 
         for (unsigned int q = 0; q < n_q_points; ++q)
           {
             const double coeff =
               1.0 / std::sqrt(1.0 + gradients[q] gradients[q]);
* 
             for (unsigned int i = 0; i < dofs_per_cell; ++i)
               cell_residual(i)
* 
-= (fe_values.shape_grad(i, q) // \nabla \phi_i
                                    coeff                    // a_n
                                    gradients[q]             // u_n
                                    fe_values.JxW(q));       // dx
           }
* 
         cell->get_dof_indices(local_dof_indices);
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
           residual(local_dof_indices[i]) += cell_residual(i);
       }
* 
     hanging_node_constraints.condense(residual);
* 
     for (types::global_dof_index i :
          DoFTools::extract_boundary_dofs(dof_handler))
       residual(i) = 0;
* 
     return residual.l2_norm();
   }
* 
 
* 
 @endcode
* 
*   <a name="MinimalSurfaceProblemdetermine_step_length"></a>  <h4>MinimalSurfaceProblem::determine_step_length</h4>
* 

* 
*  The choice of step length (or, under-relaxation factor) for the nonlinear iterations procedure remains fixed at the value chosen and discussed in   step-15  .
* 

* 
* @code
   template <int dim>
   double MinimalSurfaceProblem<dim>::determine_step_length() const
   {
     return 0.1;
   }
* 
 
* 
 @endcode
* 
*   <a name="MinimalSurfaceProblemoutput_results"></a>  <h4>MinimalSurfaceProblem::output_results</h4>
* 

* 
*  This last function to be called from `run()` outputs the current solution (and the Newton update) in graphical form as a VTU file. It is entirely the same as what has been used in previous tutorials.
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::output_results(
     const unsigned int refinement_cycle) const
   {
     DataOut<dim> data_out;
* 
     data_out.attach_dof_handler(dof_handler);
     data_out.add_data_vector(current_solution, "solution");
     data_out.add_data_vector(newton_update, "update");
     data_out.build_patches();
* 
     const std::string filename =
       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
     std::ofstream output(filename);
     data_out.write_vtu(output);
   }
* 
 
 @endcode
* 
*   <a name="MinimalSurfaceProblemrun"></a>  <h4>MinimalSurfaceProblem::run</h4>
* 

* 
*  In the run function, most remains the same as was first implemented in   step-15  . The only observable changes are that we can now choose (via the parameter file) what the final acceptable tolerance for the system residual is, and that we can choose which method of assembly we wish to utilize. To make the second choice clear, we output to the console some message which indicates the selection. Since we're interested in comparing the time taken to assemble for each of the three methods, we've also added a timer that keeps a track of how much time is spent during assembly. We also track the time taken to solve the linear system, so that we can contrast those numbers to the part of the code which would normally take the longest time to execute.
* 

* 
* @code
   template <int dim>
   void MinimalSurfaceProblem<dim>::run(const int    formulation,
                                        const double tolerance)
   {
     std::cout << "******** Assembly approach*******" << std::endl;
     const std::array<std::string, 3> method_descriptions = {
       {"Unassisted implementation (full hand linearization).",
        "Automated linearization of the finite element residual.",
        "Automated computation of finite element residual and linearization using a variational formulation."}};
     AssertIndexRange(formulation, method_descriptions.size());
     std::cout << method_descriptions[formulation] << std::endl << std::endl;
* 
 
     TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);
* 
     GridGenerator::hyper_ball(triangulation);
     triangulation.refine_global(2);
* 
     setup_system( /*first time=*/ true);
     set_boundary_values();
* 
     double       last_residual_norm = std::numeric_limits<double>::max();
     unsigned int refinement_cycle   = 0;
     do
       {
         std::cout << "Mesh refinement step " << refinement_cycle << std::endl;
* 
         if (refinement_cycle != 0)
           refine_mesh();
* 
         std::cout << "  Initial residual: " << compute_residual(0) << std::endl;
* 
         for (unsigned int inner_iteration = 0; inner_iteration < 5;
              ++inner_iteration)
           {
             {
               TimerOutput::Scope t(timer, "Assemble");
* 
               if (formulation == 0)
                 assemble_system_unassisted();
               else if (formulation == 1)
                 assemble_system_with_residual_linearization();
               else if (formulation == 2)
                 assemble_system_using_energy_functional();
               else
                 AssertThrow(false, ExcNotImplemented());
             }
* 
             last_residual_norm = system_rhs.l2_norm();
* 
             {
               TimerOutput::Scope t(timer, "Solve");
               solve();
             }
* 
 
             std::cout << "  Residual: " << compute_residual(0) << std::endl;
           }
* 
         output_results(refinement_cycle);
* 
         ++refinement_cycle;
         std::cout << std::endl;
       }
     while (last_residual_norm > tolerance);
   }
 } // namespace Step72
* 
 @endcode
* 
*   <a name="Themainfunction"></a>  <h4>The main function</h4>
* 

* 
*  Finally the main function. This follows the scheme of most other main functions, with two obvious exceptions:
* 

* 
* 
*  - We call   Utilities::MPI::MPI_InitFinalize   in order to set up (via a hidden default parameter) the number of threads using the execution of multithreaded tasks.
* 

* 
* 
*  - We also have a few lines dedicates to reading in or initializing the user-defined parameters that will be considered during the execution of the program.
* 

* 
* @code
 int main(int argc, charargv[])
 {
   try
     {
       using namespace Step72;
* 
       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
* 
       std::string prm_file;
       if (argc > 1)
         prm_file = argv[1];
       else
         prm_file = "parameters.prm";
* 
       const MinimalSurfaceProblemParameters parameters;
       ParameterAcceptor::initialize(prm_file);
* 
       MinimalSurfaceProblem<2> minimal_surface_problem_2d;
       minimal_surface_problem_2d.run(parameters.formulation,
                                      parameters.tolerance);
     }
   catch (std::exception &exc)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Exception on processing: " << std::endl
                 << exc.what() << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
* 
       return 1;
     }
   catch (...)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Unknown exception!" << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       return 1;
     }
   return 0;
 }
 @endcode
* <a name="Results"></a><h1>Results</h1>
* 

* Since there was no change to the physics of the problem that has first been analyzedin   step-15  , there is nothing to report about that. The only outwardly noticeabledifference between them is that, by default, this program will only run 9 meshrefinement steps (as opposed to   step-15  , which executes 11 refinements).This will be observable in the simulation status that appears between theheader text that prints which assembly method is being used, and the finaltimings. (All timings reported below were obtained in release mode.)
* @code
Mesh refinement step 0
  Initial residual: 1.53143
  Residual: 1.08746
  Residual: 0.966748
  Residual: 0.859602
  Residual: 0.766462
  Residual: 0.685475
* 
...
* 
Mesh refinement step 9
  Initial residual: 0.00924594
  Residual: 0.00831928
  Residual: 0.0074859
  Residual: 0.0067363
  Residual: 0.00606197
  Residual: 0.00545529
@endcode
* 
* So what is interesting for us to compare is how long the assembly process takesfor the three different implementations, and to put that into some greater context.Below is the output for the hand linearization (as computed on a circa 2012four core, eight thread laptop
* 
*  -  but we're only really interested in therelative time between the different implementations):
* @code
******** Assembly approach*******
Unassisted implementation (full hand linearization).
* 
...
* 
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      35.1s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |      1.56s |       4.5% |
| Solve                           |        50 |      30.8s |        88% |
+---------------------------------+-----------+------------+------------+
@endcode
* And for the implementation that linearizes the residual in an automatedmanner using the Sacado dynamic forward AD number type:
* @code
******** Assembly approach*******
Automated linearization of the finite element residual.
* 
...
* 
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      40.1s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |       8.8s |        22% |
| Solve                           |        50 |      28.6s |        71% |
+---------------------------------+-----------+------------+------------+
@endcode
* And, lastly, for the implementation that computes both the residual andits linearization directly from an energy functional (using nested Sacadodynamic forward AD numbers):
* @code
******** Assembly approach*******
Automated computation of finite element residual and linearization using a variational formulation.
* 
...
* 
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      48.8s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |      16.7s |        34% |
| Solve                           |        50 |      29.3s |        60% |
+---------------------------------+-----------+------------+------------+
@endcode
* 
* It's evident that the more work that is passed off to the automatic differentiationframework to perform, the more time is spent during the assembly process. Accumulatedover all refinement steps, using one level of automatic differentiation resultedin   $5.65 \times$   more computational time being spent in the assembly stage whencompared to unassisted assembly, while assembling the discrete linear system took  $10.7 \times$   longer when deriving directly from the energy functional.Unsurprisingly, the overall time spent solving the linear system remained unchanged.This means that the proportion of time spent in the solve phase to the assembly phaseshifted significantly as the number of times automated differentiation was performedat the finite element level. For many, this might mean that leveraging higher-orderdifferentiation (at the finite element level) in production code leads to anunacceptable overhead, but it may still be useful during the prototyping phase.A good compromise between the two may, therefore, be the automated linearizationof the finite element residual, which offers a lot of convenience at a measurable,but perhaps not unacceptable, cost. Alternatively, one could considernot re-building the Newton matrix in every step
* 
*  -  a topic that isexplored in substantial depth in   step-77  .
* Of course, in practice the actual overhead is very much dependent on the problem being evaluated(e.g., how many components there are in the solution, what the nature of the functionbeing differentiated is, etc.). So the exact results presented here should beinterpreted within the context of this scalar problem alone, and when it comes toother problems, some preliminary investigation by the user is certainly warranted.
* 

* <a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>
* 

* Like   step-71  , there are a few items related to automatic differentiation that couldbe evaluated further:
* 
*  - The use of other AD frameworks should be investigated, with the outlook that  alternative implementations may provide performance benefits.
* 
*  - It is also worth evaluating AD number types other than those that have been  hard-coded into this tutorial. With regard to twice differentiable types  employed at the finite-element level, mixed differentiation modes ("RAD-FAD")  should in principle be more computationally efficient than the single  mode ("FAD-FAD") types employed here. The reason that the RAD-FAD type was not  selected by default is that, at the time of writing, there remain some  bugs in its implementation within the Sacado library that lead to memory leaks.  This is documented in the   @ref auto_symb_diff   module.
* 
*  - It might be the case that using reduced precision types (i.e., `float`) as the  scalar types for the AD numbers could render a reduction in computational  expense during assembly. Using `float` as the data type for the  matrix and the residual is not unreasonable, given that the Newton  update is only meant to get us closer to the solution, but not  actuallyto* the solution; as a consequence, it makes sense to  consider using reduced-precision data types for computing these  updates, and then accumulating these updates in a solution vector  that uses the full `double` precision accuracy.
* 
*  - One further method of possibly reducing resources during assembly is to frame  the AD implementations as a constitutive model. This would be similar to the  approach adopted in   step-71  , and pushes the starting point for the automatic  differentiation one level higher up the chain of computations. This, in turn,  means that less operations are tracked by the AD library, thereby reducing the  cost of differentiating (even though one would perform the differentiation at  each cell quadrature point).
* 
*  -   step-77   is yet another variation of   step-15   that addresses a very  different part of the problem: Line search and whether it is  necessary to re-build the Newton matrix in every nonlinear  iteration. Given that the results above show that using automatic  differentiation comes at a cost, the techniques in   step-77   have the  potential to offset these costs to some degree. It would therefore  be quite interesting to combine the current program with the  techniques in   step-77  . For production codes, this would certainly be  the way to go.
* 

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-72.cc"  
*
*/
