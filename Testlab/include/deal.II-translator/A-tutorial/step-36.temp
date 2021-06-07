/**
  @page step_36 The step-36 tutorial program  
* This tutorial depends on   step-4  .
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#EigenvaluesandDirichletboundaryconditions">Eigenvalues and Dirichlet boundary conditions</a><a href="#EigenvaluesandDirichletboundaryconditions">Eigenvalues and Dirichlet boundary conditions</a>
        <li><a href="#Implementationdetails">Implementation details</a><a href="#Implementationdetails">Implementation details</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeEigenvalueProblemcodeclasstemplate">The <code>EigenvalueProblem</code> class template</a><a href="#ThecodeEigenvalueProblemcodeclasstemplate">The <code>EigenvalueProblem</code> class template</a>
        <li><a href="#ImplementationofthecodeEigenvalueProblemcodeclass">Implementation of the <code>EigenvalueProblem</code> class</a><a href="#ImplementationofthecodeEigenvalueProblemcodeclass">Implementation of the <code>EigenvalueProblem</code> class</a>
      <ul>
        <li><a href="#EigenvalueProblemEigenvalueProblem">EigenvalueProblem::EigenvalueProblem</a><a href="#EigenvalueProblemEigenvalueProblem">EigenvalueProblem::EigenvalueProblem</a>
        <li><a href="#EigenvalueProblemmake_grid_and_dofs">EigenvalueProblem::make_grid_and_dofs</a><a href="#EigenvalueProblemmake_grid_and_dofs">EigenvalueProblem::make_grid_and_dofs</a>
        <li><a href="#EigenvalueProblemassemble_system">EigenvalueProblem::assemble_system</a><a href="#EigenvalueProblemassemble_system">EigenvalueProblem::assemble_system</a>
        <li><a href="#EigenvalueProblemsolve">EigenvalueProblem::solve</a><a href="#EigenvalueProblemsolve">EigenvalueProblem::solve</a>
        <li><a href="#EigenvalueProblemoutput_results">EigenvalueProblem::output_results</a><a href="#EigenvalueProblemoutput_results">EigenvalueProblem::output_results</a>
        <li><a href="#EigenvalueProblemrun">EigenvalueProblem::run</a><a href="#EigenvalueProblemrun">EigenvalueProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Runningtheproblem">Running the problem</a><a href="#Runningtheproblem">Running the problem</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
*   <br>  
* <i>This program was contributed by Toby D. Young and Wolfgang
Bangerth.  </i>
* <a name="Preamble"></a><a name="Preamble"></a><h1>Preamble</h1>
* 

* The problem we want to solve in this example is an eigenspectrumproblem. Eigenvalue problems appear in a wide context of problems, forexample in the computation of electromagnetic standing waves incavities, vibration modes of drum membranes, or oscillations of lakesand estuaries. One of the most enigmatic applications is probably thecomputation of stationary or quasi-static wave functions in quantummechanics. The latter application is what we would like to investigatehere, though the general techniques outlined in this program are ofcourse equally applicable to the other applications above.
* Eigenspectrum problems have the general form
* @f{align*}
	L \Psi &= \varepsilon \Psi \qquad &&\text{in}\ \Omega\quad,
	\\
	\Psi &= 0 &&\text{on}\ \partial\Omega\quad,
@f}
* where the Dirichlet boundary condition on   $\Psi=\Psi(\mathbf x)$   could also bereplaced by Neumann or Robin conditions;   $L$   is an operator that generallyalso contains differential operators.
* Under suitable conditions, the above equations have a set of solutions  $\Psi_\ell,\varepsilon_\ell$  ,   $\ell\in {\cal I}$  , where   $\cal I$   canbe a finite or infinite set (and in the latter case it may be a discrete orsometimes at least in part a continuous set). In either case, let us note thatthere isno longer just a single solution, but a set of solutions (the variouseigenfunctions and corresponding eigenvalues) that we want tocompute. The problem of numerically finding all eigenvalues(eigenfunctions) of such eigenvalue problems is a formidablechallenge. In fact, if the set   $\cal I$   is infinite, the challenge isof course intractable.  Most of the time however we are really onlyinterested in a small subset of these values (functions); andfortunately, the interface to the SLEPc library that we will use forthis tutorial program allows us to select which portion of theeigenspectrum and how many solutions we want to solve for.
* In this program, the eigenspectrum solvers we use are classes providedby deal.II that wrap around the linear algebra implementation of the<a href="http://www.grycap.upv.es/slepc/" target="_top">SLEPc</a>library; SLEPc itself builds on the <a
href="http://www.mcs.anl.gov/petsc/" target="_top">PETSc</a> libraryfor linear algebra contents.
* <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1>
* 

* The basic equation of stationary quantum mechanics is theSchrödinger equation which models the motion of particles in anexternal potential   $V(\mathbf x)$  . The particle is described by a wavefunction   $\Psi(\mathbf x)$   that satisfies a relation of the(nondimensionalized) form
* @f{align*} [-\Delta + V(\mathbf x)]
\Psi(\mathbf x) &= \varepsilon \Psi(\mathbf x) \qquad &&\text{in}\
\Omega\quad, \\ \Psi &= 0 &&\text{on}\ \partial\Omega\quad.
@f}
* As a consequence, this particle can only exist in a certain number ofeigenstates that correspond to the energy eigenvalues  $\varepsilon_\ell$   admitted as solutions of this equation. Theorthodox (Copenhagen) interpretation of quantum mechanics posits that, if aparticle has energy   $\varepsilon_\ell$   then the probability of findingit at location   $\mathbf x$   is proportional to   $|\Psi_\ell(\mathbf
x)|^2$   where   $\Psi_\ell$   is the eigenfunction that corresponds to thiseigenvalue.
* In order to numerically find solutions to this equation, i.e. a set ofpairs of eigenvalues/eigenfunctions, we use the usual finite elementapproach of multiplying the equation from the left with test functions,integrating by parts, and searching for solutions in finitedimensional spaces by approximating   $\Psi(\mathbf
x)\approx\Psi_h(\mathbf x)=\sum_{j}\phi_j(\mathbf x)\tilde\psi_j$  ,where   $\tilde\psi$   is a vector of expansion coefficients. We thenimmediately arrive at the following equation that discretizes thecontinuous eigenvalue problem: @f[ \sum_j [(\nabla\phi_i,
\nabla\phi_j)+(V(\mathbf x)\phi_i,\phi_j)] \tilde{\psi}_j =
\varepsilon_h \sum_j (\phi_i, \phi_j) \tilde{\psi}_j\quad.  @f] Inmatrix and vector notation, this equation then reads: @f[ A
\tilde{\Psi} = \varepsilon_h M \tilde{\Psi} \quad, @f] where   $A$   isthe stiffness matrix arising from the differential operator   $L$  , and  $M$   is the mass matrix. The solution to the eigenvalue problem is aneigenspectrum   $\varepsilon_{h,\ell}$  , with associated eigenfunctions  $\Psi_\ell=\sum_j \phi_j\tilde{\psi}_j$  .
* 

* <a name="EigenvaluesandDirichletboundaryconditions"></a><h3>Eigenvalues and Dirichlet boundary conditions</h3>
* 

* In this program, we use Dirichlet boundary conditions for the wavefunction   $\Psi$  . What this means, from the perspective of a finiteelement code, is that only the interior degrees of freedom are realdegrees of <i>freedom</i>: the ones on the boundary are not free butare forced to have a zero value, after all. On the other hand, thefinite element method gains much of its power and simplicity fromthe fact that we just do the same thing on every cell, withouthaving to think too much about where a cell is, whether it boundson a less refined cell and consequently has a hanging node, or isadjacent to the boundary. All such checks would make the assemblyof finite element linear systems unbearably difficult to write andeven more so to read.
* Consequently, of course, when you distribute degrees of freedom withyour DoFHandler object, you don't care whether some of the degreesof freedom you enumerate are at a Dirichlet boundary. They all getnumbers. We just have to take care of these degrees of freedom at alater time when we apply boundary values. There are two basic waysof doing this (either using   MatrixTools::apply_boundary_values()  <i>after</i> assembling the linear system, or using  AffineConstraints::distribute_local_to_global()   <i>during</i> assembly;see the   @ref constraints   "constraints module" for more information),but both result in the same: a linear system that has a totalnumber of rows equal to the number of <i>all</i> degrees of freedom,including those that lie on the boundary. However, degrees offreedom that are constrained by Dirichlet conditions are separatedfrom the rest of the linear system by zeroing out the correspondingrow and column, putting a single positive entry on the diagonal,and the corresponding Dirichlet value on the right hand side.
* If you assume for a moment that we had renumbered degrees of freedomin such a way that all of those on the Dirichlet boundary come last,then the linear system we would get when solving a regular PDE witha right hand side would look like this:
* @f{align*}
  \begin{pmatrix}
    A_i & 0 \\ 0 & D_b
  \end{pmatrix}
  \begin{pmatrix}
    U_i \\ U_b
  \end{pmatrix}
  =
  \begin{pmatrix}
    F_i \\ F_b
  \end{pmatrix}.
@f}
* Here, subscripts   $i$   and   $b$   correspond to interior and boundarydegrees of freedom, respectively. The interior degrees of freedomsatisfy the linear system   $A_i U_i=F_i$   which yields the correctsolution in the interior, and boundary values are determined by  $U_b = D_b^{-1} F_b$   where   $D_b$   is a diagonal matrix that resultsfrom the process of eliminating boundary degrees of freedom, and  $F_b$   is chosen in such a way that   $U_{b,j}=D_{b,jj}^{-1} F_{b,j}$  has the correct boundary values for every boundary degree of freedom  $j$  . (For the curious, the entries of thematrix   $D_b$   result from adding modified local contributions to theglobal matrix where for the local matrices the diagonal elements, if non-zero,are set to their absolute value; otherwise, they are set to the average ofabsolute values of the diagonal. This process guarantees that the entriesof   $D_b$   are positive and of a size comparable to the rest of the diagonalentries, ensuring that the resulting matrix does not incur unreasonablelosses of accuracy due to roundoff involving matrix entries of drasticallydifferent size. The actual values that end up on the diagonal are difficultto predict and you should treat them as arbitrary and unpredictable, butpositive.)
* For "regular" linear systems, this all leads to the correct solution.On the other hand, for eigenvalue problems, this is not so trivial.There, eliminating boundary values affects both matrices  $A$   and   $M$   that we will solve with in the current tutorial program.After elimination of boundary values, we then receive an eigenvalueproblem that can be partitioned like this:
* @f{align*}
  \begin{pmatrix}
    A_i & 0 \\ 0 & D_A
  \end{pmatrix}
  \begin{pmatrix}
    \tilde\Psi_i \\ \tilde\Psi_b
  \end{pmatrix}
  =
  \epsilon_h
  \begin{pmatrix}
    M_i & 0 \\ 0 & D_M
  \end{pmatrix}
  \begin{pmatrix}
    \tilde\Psi_i \\ \tilde\Psi_b
  \end{pmatrix}.
@f}
* This form makes it clear that there are two sets of eigenvalues:the ones we care about, and spurious eigenvalues from theseparated problem@f[
  D_A \tilde \Psi_b = \epsilon_h D_M \Psi_b.
@f]These eigenvalues are spurious since they result from an eigenvaluesystem that operates only on boundary nodes
* 
*  -  nodes that are notreal degrees of <i>freedom</i>.Of course, since the two matrices   $D_A,D_M$   are diagonal, we canexactly quantify these spurious eigenvalues: they are  $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$   (where the indices  $j$   corresponds exactly to the degrees of freedom that are constrainedby Dirichlet boundary values).
* So how does one deal with them? The fist part is to recognize when oureigenvalue solver finds one of them. To this end, the program computesand prints an interval within which these eigenvalues lie, by computingthe minimum and maximum of the expression   $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$  over all constrained degrees of freedom. In the program below, thisalready suffices: we find that this interval lies outside the set ofsmallest eigenvalues and corresponding eigenfunctions we are interestedin and compute, so there is nothing we need to do here.
* On the other hand, it may happen that we find that one of the eigenvalueswe compute in this program happens to be in this interval, and in thatcase we would not know immediately whether it is a spurious or a trueeigenvalue. In that case, one could simply scale the diagonal elements ofeither matrix after computing the two matrices,thus shifting them away from the frequency of interest in the eigen-spectrum.This can be done by using the following code, making sure that all spuriouseigenvalues are exactly equal to   $1.234\cdot 10^5$  :
* @code
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          stiffness_matrix.set(i, i, 1.234e5);
          mass_matrix.set(i, i, 1);
        }
@endcode
* However, this strategy is not pursued here as the spurious eigenvalueswe get from our program as-is happen to be greater than the lowestfive that we will calculate and are interested in.
* 

* <a name="Implementationdetails"></a><h3>Implementation details</h3>
* 

* The program below is essentially just a slightly modified version of  step-4  . The things that are different are the following:
*   <ul>  
*   <li>  The main class (named   <code>EigenvalueProblem</code>  ) now nolonger has a single solution vector, but a whole set of vectors forthe various eigenfunctions we want to compute. Moreover, the  <code>main</code>   function, which has the top-level control overeverything here, initializes and finalizes the interface to SLEPc andPETSc simultaneously via   <code>SlepcInitialize</code>   and  <code>SlepFinalize</code>  .  </li>  
*   <li>  We use PETSc matrices and vectors as in   step-17   and  step-18   since that is what the SLEPc eigenvalue solversrequire.  </li>  
*   <li>  The function   <code>EigenvalueProblem::solve</code>   is entirelydifferent from anything seen so far in the tutorial, as it does notjust solve a linear system but actually solves the eigenvalue problem.It is built on the SLEPc library, and more immediately on the deal.IISLEPc wrappers in the class   SLEPcWrappers::SolverKrylovSchur.</li>  
*   <li>  We use the ParameterHandler class to describe a few inputparameters, such as the exact form of the potential   $V({\mathbf
x})$  , the number of global refinement steps of the mesh,or the number of eigenvalues we want to solve for. We could go muchfurther with this but stop at making only a few of the things that onecould select at run time actual input file parameters. In order to seewhat could be done in this regard, take a look at   @ref step_29  "  step-29  " and   step-33  .  </li>  
*   <li>  We use the FunctionParser class to make the potential   $V(\mathbf
x)$   a run-time parameter that can be specified in the input file as aformula.  </li>  
*   </ul>  
* The rest of the program follows in a pretty straightforward way from  step-4  .
* 

*  <a name="CommProg"></a> <h1> The commented program</h1>
*   <a name="Includefiles"></a>  <h3>Include files</h3>
* 

* 
*  As mentioned in the introduction, this program is essentially only a slightly revised version of   step-4  . As a consequence, most of the following include files are as used there, or at least as used already in previous tutorial programs:
* 

* 
* @code
 #include <deal.II/base/logstream.h>
 #include <deal.II/base/quadrature_lib.h>
 #include <deal.II/base/function.h>
 #include <deal.II/base/function_parser.h>
 #include <deal.II/base/parameter_handler.h>
 #include <deal.II/base/utilities.h>
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/dofs/dof_handler.h>
 #include <deal.II/dofs/dof_tools.h>
 #include <deal.II/fe/fe_q.h>
 #include <deal.II/fe/fe_values.h>
 #include <deal.II/numerics/vector_tools.h>
 #include <deal.II/numerics/matrix_tools.h>
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/lac/affine_constraints.h>
 #include <deal.II/lac/full_matrix.h>
* 
 @endcode
* 
*  IndexSet is used to set the size of each   PETScWrappers::MPI::Vector:  
* 

* 
* @code
 #include <deal.II/base/index_set.h>
* 
 @endcode
* 
*  PETSc appears here because SLEPc depends on this library:
* 

* 
* @code
 #include <deal.II/lac/petsc_sparse_matrix.h>
 #include <deal.II/lac/petsc_vector.h>
* 
 @endcode
* 
*  And then we need to actually import the interfaces for solvers that SLEPc provides:
* 

* 
* @code
 #include <deal.II/lac/slepc_solver.h>
* 
 @endcode
* 
*  We also need some standard C++:
* 

* 
* @code
 #include <fstream>
 #include <iostream>
* 
 @endcode
* 
*  Finally, as in previous programs, we import all the deal.II class and function names into the namespace into which everything in this program will go:
* 

* 
* @code
 namespace Step36
 {
   using namespace dealii;
* 
 @endcode
* 
*   <a name="ThecodeEigenvalueProblemcodeclasstemplate"></a>  <h3>The <code>EigenvalueProblem</code> class template</h3>
* 

* 
*  Following is the class declaration for the main class template. It looks pretty much exactly like what has already been shown in   step-4  :
* 

* 
* @code
   template <int dim>
   class EigenvalueProblem
   {
   public:
     EigenvalueProblem(const std::string &prm_file);
     void run();
* 
   private:
     void         make_grid_and_dofs();
     void         assemble_system();
     unsigned int solve();
     void         output_results() const;
* 
     Triangulation<dim> triangulation;
     FE_Q<dim>          fe;
     DoFHandler<dim>    dof_handler;
* 
 @endcode
* 
*  With these exceptions: For our eigenvalue problem, we need both a stiffness matrix for the left hand side as well as a mass matrix for the right hand side. We also need not just one solution function, but a whole set of these for the eigenfunctions we want to compute, along with the corresponding eigenvalues:
* 

* 
* @code
     PETScWrappers::SparseMatrix             stiffness_matrix, mass_matrix;
     std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
     std::vector<double>                     eigenvalues;
* 
 @endcode
* 
*  And then we need an object that will store several run-time parameters that we will specify in an input file:
* 

* 
* @code
     ParameterHandler parameters;
* 
 @endcode
* 
*  Finally, we will have an object that contains "constraints" on our degrees of freedom. This could include hanging node constraints if we had adaptively refined meshes (which we don't have in the current program). Here, we will store the constraints for boundary nodes   $U_i=0$  .
* 

* 
* @code
     AffineConstraints<double> constraints;
   };
* 
 @endcode
* 
*   <a name="ImplementationofthecodeEigenvalueProblemcodeclass"></a>  <h3>Implementation of the <code>EigenvalueProblem</code> class</h3>
* 

* 
*   <a name="EigenvalueProblemEigenvalueProblem"></a>  <h4>EigenvalueProblem::EigenvalueProblem</h4>
* 

* 
*  First up, the constructor. The main new part is handling the run-time input parameters. We need to declare their existence first, and then read their values from the input file whose name is specified as an argument to this function:
* 

* 
* @code
   template <int dim>
   EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file)
     : fe(1)
     , dof_handler(triangulation)
   {
 @endcode
* 
*  TODO investigate why the minimum number of refinement steps required to obtain the correct eigenvalue degeneracies is 6
* 

* 
* @code
     parameters.declare_entry(
       "Global mesh refinement steps",
       "5",
       Patterns::Integer(0, 20),
       "The number of times the 1-cell coarse mesh should "
       "be refined globally for our computations.");
     parameters.declare_entry("Number of eigenvalues/eigenfunctions",
                              "5",
                              Patterns::Integer(0, 100),
                              "The number of eigenvalues/eigenfunctions "
                              "to be computed.");
     parameters.declare_entry("Potential",
                              "0",
                              Patterns::Anything(),
                              "A functional description of the potential.");
* 
     parameters.parse_input(prm_file);
   }
* 
 
 @endcode
* 
*   <a name="EigenvalueProblemmake_grid_and_dofs"></a>  <h4>EigenvalueProblem::make_grid_and_dofs</h4>
* 

* 
*  The next function creates a mesh on the domain   $[-1,1]^d$  , refines it as many times as the input file calls for, and then attaches a DoFHandler to it and initializes the matrices and vectors to their correct sizes. We also build the constraints that correspond to the boundary values   $u|_{\partial\Omega}=0$  .   
*   For the matrices, we use the PETSc wrappers. These have the ability to allocate memory as necessary as non-zero entries are added. This seems inefficient: we could as well first compute the sparsity pattern, initialize the matrices with it, and as we then insert entries we can be sure that we do not need to re-allocate memory and free the one used previously. One way to do that would be to use code like this:   <div class=CodeFragmentInTutorialComment>


* 
* 
@code
   DynamicSparsityPattern
      dsp (dof_handler.n_dofs(),
           dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, dsp);
   dsp.compress ();
   stiffness_matrix.reinit (dsp);
   mass_matrix.reinit (dsp);
 @endcode
* 
 </div>   instead of the two   <code>reinit()</code>   calls for the stiffness and mass matrices below.   
*   This doesn't quite work, unfortunately. The code above may lead to a few entries in the non-zero pattern to which we only ever write zero entries; most notably, this holds true for off-diagonal entries for those rows and columns that belong to boundary nodes. This shouldn't be a problem, but for whatever reason, PETSc's ILU preconditioner, which we use to solve linear systems in the eigenvalue solver, doesn't like these extra entries and aborts with an error message.   
*   In the absence of any obvious way to avoid this, we simply settle for the second best option, which is have PETSc allocate memory as necessary. That said, since this is not a time critical part, this whole affair is of no further importance.
* 

* 
* @code
   template <int dim>
   void EigenvalueProblem<dim>::make_grid_and_dofs()
   {
     GridGenerator::hyper_cube(triangulation,
* 
-1, 1);
     triangulation.refine_global(
       parameters.get_integer("Global mesh refinement steps"));
     dof_handler.distribute_dofs(fe);
* 
     DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
     constraints.close();
* 
     stiffness_matrix.reinit(dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
     mass_matrix.reinit(dof_handler.n_dofs(),
                        dof_handler.n_dofs(),
                        dof_handler.max_couplings_between_dofs());
* 
 @endcode
* 
*  The next step is to take care of the eigenspectrum. In this case, the outputs are eigenvalues and eigenfunctions, so we set the size of the list of eigenfunctions and eigenvalues to be as large as we asked for in the input file. When using a   PETScWrappers::MPI::Vector,   the Vector is initialized using an IndexSet. IndexSet is used not only to resize the   PETScWrappers::MPI::Vector   but it also associates an index in the   PETScWrappers::MPI::Vector   with a degree of freedom (see   step-40   for a more detailed explanation). The function complete_index_set() creates an IndexSet where every valid index is part of the set. Note that this program can only be run sequentially and will throw an exception if used in parallel.
* 

* 
* @code
     IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
     eigenfunctions.resize(
       parameters.get_integer("Number of eigenvalues/eigenfunctions"));
     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
       eigenfunctions[i].reinit(eigenfunction_index_set, MPI_COMM_WORLD);
* 
     eigenvalues.resize(eigenfunctions.size());
   }
* 
 
 @endcode
* 
*   <a name="EigenvalueProblemassemble_system"></a>  <h4>EigenvalueProblem::assemble_system</h4>
* 

* 
*  Here, we assemble the global stiffness and mass matrices from local contributions   $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf x) \cdot
 \nabla\varphi_j(\mathbf x) + V(\mathbf x)\varphi_i(\mathbf
 x)\varphi_j(\mathbf x)$   and   $M^K_{ij} = \int_K \varphi_i(\mathbf
 x)\varphi_j(\mathbf x)$   respectively. This function should be immediately familiar if you've seen previous tutorial programs. The only thing new would be setting up an object that described the potential   $V(\mathbf x)$   using the expression that we got from the input file. We then need to evaluate this object at the quadrature points on each cell. If you've seen how to evaluate function objects (see, for example the coefficient in   step-5  ), the code here will also look rather familiar.
* 

* 
* @code
   template <int dim>
   void EigenvalueProblem<dim>::assemble_system()
   {
     QGauss<dim> quadrature_formula(fe.degree + 1);
* 
     FEValues<dim> fe_values(fe,
                             quadrature_formula,
                             update_values | update_gradients |
                               update_quadrature_points | update_JxW_values);
* 
     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
     const unsigned int n_q_points    = quadrature_formula.size();
* 
     FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
* 
     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
* 
     FunctionParser<dim> potential;
     potential.initialize(FunctionParser<dim>::default_variable_names(),
                          parameters.get("Potential"),
                          typename FunctionParser<dim>::ConstMap());
* 
     std::vector<double> potential_values(n_q_points);
     for (const auto &cell : dof_handler.active_cell_iterators())
       {
         fe_values.reinit(cell);
         cell_stiffness_matrix = 0;
         cell_mass_matrix      = 0;
* 
         potential.value_list(fe_values.get_quadrature_points(),
                              potential_values);
* 
         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
           for (unsigned int i = 0; i < dofs_per_cell; ++i)
             for (unsigned int j = 0; j < dofs_per_cell; ++j)
               {
                 cell_stiffness_matrix(i, j) +=           
                   (fe_values.shape_grad(i, q_point)    
                      fe_values.shape_grad(j, q_point)    
                    +                                     
                    potential_values[q_point]           
                      fe_values.shape_value(i, q_point) 
                      fe_values.shape_value(j, q_point)   
                    )                                   
                   fe_values.JxW(q_point);                
* 
                 cell_mass_matrix(i, j) +=              
                   (fe_values.shape_value(i, q_point) 
                    fe_values.shape_value(j, q_point)   
                    )                                 
                   fe_values.JxW(q_point);              
               }
* 
 @endcode
* 
*  Now that we have the local matrix contributions, we transfer them into the global objects and take care of zero boundary constraints:
* 

* 
* @code
         cell->get_dof_indices(local_dof_indices);
* 
         constraints.distribute_local_to_global(cell_stiffness_matrix,
                                                local_dof_indices,
                                                stiffness_matrix);
         constraints.distribute_local_to_global(cell_mass_matrix,
                                                local_dof_indices,
                                                mass_matrix);
       }
* 
 @endcode
* 
*  At the end of the function, we tell PETSc that the matrices have now been fully assembled and that the sparse matrix representation can now be compressed as no more entries will be added:
* 

* 
* @code
     stiffness_matrix.compress(VectorOperation::add);
     mass_matrix.compress(VectorOperation::add);
* 
 
 @endcode
* 
*  Before leaving the function, we calculate spurious eigenvalues, introduced to the system by zero Dirichlet constraints. As discussed in the introduction, the use of Dirichlet boundary conditions coupled with the fact that the degrees of freedom located at the boundary of the domain remain part of the linear system we solve, introduces a number of spurious eigenvalues. Below, we output the interval within which they all lie to ensure that we can ignore them should they show up in our computations.
* 

* 
* @code
     double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
            max_spurious_eigenvalue =
* 
-std::numeric_limits<double>::max();
* 
     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
       if (constraints.is_constrained(i))
         {
           const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i);
           min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
           max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
         }
* 
     std::cout << "   Spurious eigenvalues are all in the interval "
               << "[" << min_spurious_eigenvalue << ","
               << max_spurious_eigenvalue << "]" << std::endl;
   }
* 
 
 @endcode
* 
*   <a name="EigenvalueProblemsolve"></a>  <h4>EigenvalueProblem::solve</h4>
* 

* 
*  This is the key new functionality of the program. Now that the system is set up, here is a good time to actually solve the problem: As with other examples this is done using a "solve" routine. Essentially, it works as in other programs: you set up a SolverControl object that describes the accuracy to which we want to solve the linear systems, and then we select the kind of solver we want. Here we choose the Krylov-Schur solver of SLEPc, a pretty fast and robust choice for this kind of problem:
* 

* 
* @code
   template <int dim>
   unsigned int EigenvalueProblem<dim>::solve()
   {
 @endcode
* 
*  We start here, as we normally do, by assigning convergence control we want:
* 

* 
* @code
     SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9);
     SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control);
* 
 @endcode
* 
*  Before we actually solve for the eigenfunctions and
* 
*  - alues, we have to also select which set of eigenvalues to solve for. Lets select those eigenvalues and corresponding eigenfunctions with the smallest real part (in fact, the problem we solve here is symmetric and so the eigenvalues are purely real). After that, we can actually let SLEPc do its work:
* 

* 
* @code
     eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL);
* 
     eigensolver.set_problem_type(EPS_GHEP);
* 
     eigensolver.solve(stiffness_matrix,
                       mass_matrix,
                       eigenvalues,
                       eigenfunctions,
                       eigenfunctions.size());
* 
 @endcode
* 
*  The output of the call above is a set of vectors and values. In eigenvalue problems, the eigenfunctions are only determined up to a constant that can be fixed pretty arbitrarily. Knowing nothing about the origin of the eigenvalue problem, SLEPc has no other choice than to normalize the eigenvectors to one in the   $l_2$   (vector) norm. Unfortunately this norm has little to do with any norm we may be interested from a eigenfunction perspective: the   $L_2(\Omega)$   norm, or maybe the   $L_\infty(\Omega)$   norm.     
*   Let us choose the latter and rescale eigenfunctions so that they have   $\|\phi_i(\mathbf x)\|_{L^\infty(\Omega)}=1$   instead of   $\|\Phi\|_{l_2}=1$   (where   $\phi_i$   is the   $i$  th eigen<i>function</i> and   $\Phi_i$   the corresponding vector of nodal values). For the   $Q_1$   elements chosen here, we know that the maximum of the function   $\phi_i(\mathbf x)$   is attained at one of the nodes, so   $\max_{\mathbf
 x}\phi_i(\mathbf x)=\max_j (\Phi_i)_j$  , making the normalization in the   $L_\infty$   norm trivial. Note that this doesn't work as easily if we had chosen   $Q_k$   elements with   $k>1$  : there, the maximum of a function does not necessarily have to be attained at a node, and so   $\max_{\mathbf x}\phi_i(\mathbf x)\ge\max_j (\Phi_i)_j$   (although the equality is usually nearly true).
* 

* 
* @code
     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
       eigenfunctions[i] /= eigenfunctions[i].linfty_norm();
* 
 @endcode
* 
*  Finally return the number of iterations it took to converge:
* 

* 
* @code
     return solver_control.last_step();
   }
* 
 
 @endcode
* 
*   <a name="EigenvalueProblemoutput_results"></a>  <h4>EigenvalueProblem::output_results</h4>
* 

* 
*  This is the last significant function of this program. It uses the DataOut class to generate graphical output from the eigenfunctions for later visualization. It works as in many of the other tutorial programs.   
*   The whole collection of functions is then output as a single VTK file.
* 

* 
* @code
   template <int dim>
   void EigenvalueProblem<dim>::output_results() const
   {
     DataOut<dim> data_out;
* 
     data_out.attach_dof_handler(dof_handler);
* 
     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
       data_out.add_data_vector(eigenfunctions[i],
                                std::string("eigenfunction_") +
                                  Utilities::int_to_string(i));
* 
 @endcode
* 
*  The only thing worth discussing may be that because the potential is specified as a function expression in the input file, it would be nice to also have it as a graphical representation along with the eigenfunctions. The process to achieve this is relatively straightforward: we build an object that represents   $V(\mathbf x)$   and then we interpolate this continuous function onto the finite element space. The result we also attach to the DataOut object for visualization.
* 

* 
* @code
     Vector<double> projected_potential(dof_handler.n_dofs());
     {
       FunctionParser<dim> potential;
       potential.initialize(FunctionParser<dim>::default_variable_names(),
                            parameters.get("Potential"),
                            typename FunctionParser<dim>::ConstMap());
       VectorTools::interpolate(dof_handler, potential, projected_potential);
     }
     data_out.add_data_vector(projected_potential, "interpolated_potential");
* 
     data_out.build_patches();
* 
     std::ofstream output("eigenvectors.vtk");
     data_out.write_vtk(output);
   }
* 
 
 @endcode
* 
*   <a name="EigenvalueProblemrun"></a>  <h4>EigenvalueProblem::run</h4>
* 

* 
*  This is the function which has the top-level control over everything. It is almost exactly the same as in   step-4  :
* 

* 
* @code
   template <int dim>
   void EigenvalueProblem<dim>::run()
   {
     make_grid_and_dofs();
* 
     std::cout << "   Number of active cells:       "
               << triangulation.n_active_cells() << std::endl
               << "   Number of degrees of freedom: " << dof_handler.n_dofs()
               << std::endl;
* 
     assemble_system();
* 
     const unsigned int n_iterations = solve();
     std::cout << "   Solver converged in " << n_iterations << " iterations."
               << std::endl;
* 
     output_results();
* 
     std::cout << std::endl;
     for (unsigned int i = 0; i < eigenvalues.size(); ++i)
       std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i]
                 << std::endl;
   }
 } // namespace Step36
* 
 @endcode
* 
*   <a name="Thecodemaincodefunction"></a>  <h3>The <code>main</code> function</h3>
* 

* 
* @code
 int main(int argc, char*argv)
 {
   try
     {
       using namespace dealii;
       using namespace Step36;
* 
       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
* 
 
 @endcode
* 
*  This program can only be run in serial. Otherwise, throw an exception.
* 

* 
* @code
       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                   ExcMessage(
                     "This program can only be run in serial, use ./step-36"));
* 
       EigenvalueProblem<2> problem("step-36.prm");
       problem.run();
     }
* 
 @endcode
* 
*  All the while, we are watching out if any exceptions should have been generated. If that is so, we panic...
* 

* 
* @code
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
* 
 @endcode
* 
*  If no exceptions are thrown, then we tell the program to stop monkeying around and exit nicely:
* 

* 
* @code
   std::cout << std::endl << "   Job done." << std::endl;
* 
   return 0;
 }
 @endcode
* <a name="Results"></a><h1>Results</h1>
* 

* <a name="Runningtheproblem"></a><h3>Running the problem</h3>
* 

* The problem's input is parameterized by an input file   <code>\step-36.prm</code>  which could, for example, contain the following text:
* @code
set Global mesh refinement steps         = 5
set Number of eigenvalues/eigenfunctions = 5
set Potential                            = 0
@endcode
* 
* Here, the potential is zero inside the domain, and we know that theeigenvalues are given by   $\lambda_{(mn)}=\frac{\pi^2}{4}(m^2+n^2)$   where  $m,n\in{\mathbb N^+}$  . Eigenfunctions are sines and cosines with   $m$   and   $n$  periods in   $x$   and   $y$   directions. This matches the output our programgenerates:
* @code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       1024
   Number of degrees of freedom: 1089
   Solver converged in 67 iterations.
* 
      Eigenvalue 0 : 4.93877
      Eigenvalue 1 : 12.3707
      Eigenvalue 2 : 12.3707
      Eigenvalue 3 : 19.8027
      Eigenvalue 4 : 24.837
* 
   Job done.  @endcode These eigenvalues are exactly the ones thatcorrespond to pairs   $(m,n)=(1,1)$  ,   $(1,2)$   and   $(2,1)$  ,   $(2,2)$  , and  $(3,1)$  . A visualization of the corresponding eigenfunctions wouldlook like this:
*   <table width="80%">
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.0.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.1.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.2.png" alt=""></td>
</tr>
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.3.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.4.png" alt=""></td>
<td></td>
</tr>
</table>  
* <a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>
* 

* It is always worth playing a few games in the playground! So here goeswith a few suggestions:
*   <ul>  
*   <li>   The potential used above (called the <i>infinite well</i> becauseit is a flat potential surrounded by infinitely high walls) isinteresting because it allows for analytically known solutions. Apartfrom that, it is rather boring, however. That said, it is trivial toplay around with the potential by just setting it to somethingdifferent in the input file. For example, let us assume that we wantedto work with the following potential in2d:@f[
  V(x,y) = \left\{
       \begin{array}{ll}
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
-100 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy>0
         \\
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
-5 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy\le 0
         \\
         0 & \text{otherwise}
      \end{array} \right.\quad.
@f]In other words, the potential is
* 
*  - 00 in two sectors of a circle of radius0.75,
* 
*  -  in the other two sectors, and zero outside the circle. We can achievethis by using the following in the input file:
* @code
set Potential = if (x^2 + y^2 < 0.75^2, if (x*y > 0,
* 
-100,
* 
-5), 0)
@endcode
* If in addition we also increase the mesh refinement by one level, we get thefollowing results:
* @code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       4096
   Number of degrees of freedom: 4225
* 
   Eigenvalue 0 :
* 
-74.2562
   Eigenvalue 1 :
* 
-72.7322
   Eigenvalue 2 :
* 
-42.7406
   Eigenvalue 3 :
* 
-42.2232
   Eigenvalue 4 :
* 
-37.0744
@endcode
* 
* The output file also contains an interpolated version of the potential, whichlooks like this (note that as expected the lowest few eigenmodes haveprobability densities   $|\Psi(\mathbf x)|^2$   that are significant only where thepotential is the lowest, i.e. in the top right and bottom left sector of innercircle of the potential):
*   <img src="https://www.dealii.org/images/steps/developer/step-36.mod.potential.png" alt="">  
* The first five eigenfunctions are now like this:
*   <table width="80%">
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.0.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.1.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.2.png" alt=""></td>
</tr>
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.3.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.4.png" alt=""></td>
<td></td>
</tr>
</table>  
*   <li>   In our derivation of the problem we have assumed that theparticle is confined to a domain   $\Omega$   and that at the boundary ofthis domain its probability   $|\Psi|^2$   of being is zero. This isequivalent to solving the eigenvalue problem on all of   ${\mathbb R}^d$  and assuming that the energy potential is finite only inside a region  $\Omega$   and infinite outside. It is relatively easy to show that  $|\Psi(\mathbf x)|^2$   at all locations   $\mathbf x$   where   $V(\mathbf
x)=\infty$  . So the question is what happens if our potential is not ofthis form, i.e. there is no bounded domain outside of which thepotential is infinite? In that case, it may be worth to just considera very large domain at the boundary of which   $V(\mathbf x)$   is atleast very large, if not infinite. Play around with a few cases likethis and explore how the spectrum and eigenfunctions change as we makethe computational region larger and larger.
*   <li>   What happens if we investigate the simple harmonic oscillatorproblem   $V(\mathbf x)=c|\mathbf x|^2$  ? This potential is exactly ofthe form discussed in the previous paragraph and has hyper sphericalsymmetry. One may want to use a large spherical domain with a largeouter radius, to approximate the whole-space problem (say, by invoking  GridGenerator::hyper_ball).  
*   <li>   The plots above show the wave function   $\Psi(\mathbf x)$  , but thephysical quantity of interest is actually the probability density  $|\Psi(\mathbf x)|^2$   for the particle to be at location   $\mathbf x$  .Some visualization programs can compute derived quantities from the data inan input file, but we can also do so right away when creating the outputfile. The facility to do that is the DataPostprocessor class that canbe used in conjunction with the DataOut class. Examples of how thiscan be done can be found in   step-29   and  step-33  .
*   <li>   What happens if the particle in the box has %internal degrees offreedom? For example, if the particle were a spin-  $1/2$   particle? Inthat case, we may want to start solving a vector-valued probleminstead.
*   <li>   Our implementation of the deal.II library here uses thePETScWrappers and SLEPcWrappers and is suitable for running on serialmachine architecture. However, for larger grids and with a largernumber of degrees-of-freedom, we may want to run our application onparallel architectures. A parallel implementation of the above codecan be particularly useful here since the generalized eigenspectrumproblem is somewhat more expensive to solve than the standard problemsconsidered in most of the earlier tutorials. Fortunately, modifying the aboveprogram to be MPI compliant is a relatively straightforwardprocedure. A sketch of how this can be done can be found in   @ref
step_17   "  step-17  ".
*   <li>   Finally, there are alternatives to using the SLEPc eigenvaluesolvers. deal.II has interfaces to one of them, ARPACK (see <a
href="../../external-libs/arpack.html">the ARPACK configuration page</a> forsetup instructions), implemented in the ArpackSolver class. Here is a short andquick overview of what one would need to change to use it, provided you have aworking installation of ARPACK and deal.II has been configured properly for it(see the deal.II <a href="../../readme.html" target="body">README</a> file):
* First, in order to use the ARPACK interfaces, we can go back to using standarddeal.II matrices and vectors, so we start by replacing the PETSc and SLEPcheaders
* @code
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
@endcode
* with these:
* @code
#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
@endcode
* ARPACK allows complex eigenvalues, so we will also need
* @code
#include <complex>
@endcode
* 
* Secondly, we switch back to the deal.II matrix and vector definitions in themain class:
* @code
    SparsityPattern                     sparsity_pattern;
    SparseMatrix<double>                stiffness_matrix, mass_matrix;
    std::vector<Vector<double> >        eigenfunctions;
    std::vector<std::complex<double>>   eigenvalues;
@endcode
* and initialize them as usual in   <code>make_grid_and_dofs()</code>  :
* @code
    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
* 
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    constraints.condense (sparsity_pattern);
    sparsity_pattern.compress();
* 
    stiffness_matrix.reinit (sparsity_pattern);
    mass_matrix.reinit (sparsity_pattern);
@endcode
* 
* For solving the eigenvalue problem with ARPACK, we finally need to modify  <code>solve()</code>  :
* @code
  template <int dim>
  unsigned int EigenvalueProblem<dim>::solve ()
  {
    SolverControl solver_control (dof_handler.n_dofs(), 1e-9);
* 
    SparseDirectUMFPACK inverse;
    inverse.initialize (stiffness_matrix);
* 
    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors);
* 
    ArpackSolver eigensolver (solver_control, additional_data);
    eigensolver.solve (stiffness_matrix,
                       mass_matrix,
                       inverse,
                       eigenvalues,
                       eigenfunctions,
                       eigenvalues.size());
* 
    for (unsigned int i=0; i<eigenfunctions.size(); ++i)
      eigenfunctions[i] /= eigenfunctions[i].linfty_norm ();
* 
    return solver_control.last_step ();
  }
@endcode
* Note how we have used an exact decomposition (using SparseDirectUMFPACK) as apreconditioner to ARPACK.  </ul>  
* 

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-36.cc"  
* */
