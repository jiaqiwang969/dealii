/**
 * @page step_40 The step-40 tutorial program
 * This tutorial depends on   step-6  . @htmlonly <table class="tutorial"
 * width="50%"> <tr><th colspan="2"><b><small>Table of
 * contents</small></b><b><small>Table of contents</small></b></th></tr>
 * <tr><td width="50%" valign="top">
 * <ol>
 * <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro"
 * class=bold>Introduction</a>
 * <ul>
 * <li><a href="#Thetestcase">The testcase</a><a href="#Thetestcase">The
 * testcase</a>
 * </ul>
 * <li> <a href="#CommProg" class=bold>The commented program</a><a
 * href="#CommProg" class=bold>The commented program</a>
 * <ul>
 * <li><a href="#Includefiles">Include files</a><a
 * href="#Includefiles">Include files</a>
 * <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The
 * <code>LaplaceProblem</code> class template</a><a
 * href="#ThecodeLaplaceProblemcodeclasstemplate">The
 * <code>LaplaceProblem</code> class template</a>
 * <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The
 * <code>LaplaceProblem</code> class implementation</a><a
 * href="#ThecodeLaplaceProblemcodeclassimplementation">The
 * <code>LaplaceProblem</code> class implementation</a>
 * <ul>
 * <li><a href="#Constructor">Constructor</a><a
 * href="#Constructor">Constructor</a>
 * <li><a
 * href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a><a
 * href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
 * <li><a
 * href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a><a
 * href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
 * <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a><a
 * href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
 * <li><a href="#LaplaceProblemrefine_grid">LaplaceProblem::refine_grid</a><a
 * href="#LaplaceProblemrefine_grid">LaplaceProblem::refine_grid</a>
 * <li><a
 * href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a><a
 * href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
 * <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a><a
 * href="#LaplaceProblemrun">LaplaceProblem::run</a>
 * <li><a href="#main">main()</a><a href="#main">main()</a>
 * </ul>
 * </ul>
 * </ol></td><td width="50%" valign="top"><ol>
 * <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results"
 * class=bold>Results</a>
 * <ul>
 * <li><a href="#Possibilitiesforextensions">Possibilities for
 * extensions</a><a href="#Possibilitiesforextensions">Possibilities for
 * extensions</a>
 * </ul>
 * <li> <a href="#PlainProg" class=bold>The plain program</a><a
 * href="#PlainProg" class=bold>The plain program</a>
 * </ol> </td> </tr> </table>
 * @endhtmlonly
 * <br>
 * <i>This program was contributed by Timo Heister, Martin Kronbichler and
 * Wolfgang Bangerth.
 * <br>
 * This material is based upon work partly supported by the National Science
 * Foundation under Award No. EAR-0426271 and The California Institute of
 * Technology. Any opinions, findings, and conclusions or recommendations
 * expressed in this publication are those of the author and do not
 * necessarily reflect the views of the National Science Foundation or of The
 * California Institute of Technology. </i>
 *
 *
 * @note   As a prerequisite of this program, you need to have both PETSc and
 * thep4est library installed. The installation of deal.IItogether with these
 * two additional libraries is described in the <a href="../../readme.html"
 * target="body">README</a> file. Note also thatto work properly, this program
 * needs access to the Hyprepreconditioner package implementing algebraic
 * multigrid; it can beinstalled as part of PETSc but has to be explicitly
 * enabled duringPETSc configuration; see the page linked to from the
 * installationinstructions for PETSc.
 *
 *  <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1>
 *
 *  @dealiiVideoLecture{41.5,41.75} Given today's computers, most finite
 * element computations can be done ona single machine. The majority of
 * previous tutorial programs thereforeshows only this, possibly splitting up
 * work among a number ofprocessors that, however, can all access the same,
 * shared memoryspace. That said, there are problems that are simply too big
 * for asingle machine and in that case the problem has to be split up in
 * asuitable way among multiple machines each of which contributes itspart to
 * the whole. A simple way to do that was shown in   step-17   and  step-18  ,
 * where we show how a program can use <a
 * Consequently, we need to approach the problem differently: to scale tovery
 * large problems each processor can only store its own little pieceof the
 * Triangulation and DoFHandler objects. deal.II implements such ascheme in
 * the   parallel::distributed   namespace and the classestherein. It builds
 * on an external library, <a href="http://www.p4est.org/">p4est</a> (a play
 * on the expression<i>parallel forest</i> that describes the parallel storage
 * of ahierarchically constructed mesh as a forest of quad- oroct-trees). You
 * need to <a href="../../external-libs/p4est.html">install and configure
 * p4est</a>but apart from that all of its workings are hidden under the
 * surfaceof deal.II. In essence, what the
 * parallel::distributed::Triangulation   class andcode inside the DoFHandler
 * class do is to splitthe global mesh so that every processor only stores a
 * small bit it"owns" along with one layer of "ghost" cells that surround the
 * ones itowns. What happens in the rest of the domain on which we want to
 * solvethe partial differential equation is unknown to each processor and
 * canonly be inferred through communication with other machines if
 * suchinformation is needed. This implies that we also have to think
 * aboutproblems in a different way than we did in, for example,   step-17
 * and  step-18  : no processor can have the entire solution vector
 * forpostprocessing, for example, and every part of a program has to
 * beparallelized because no processor has all the information necessaryfor
 * sequential operations.
 * A general overview of how this parallelization happens is described inthe   @ref distributed   documentation module. You should read it for atop-level overview before reading through the source code of thisprogram. A concise discussion of many terms we will use in the programis also provided in the   @ref distributed_paper   "Distributed Computing paper".
 * It is probably worthwhile reading it for background information on
 * howthings work internally in this program.
 *
 *  <a name="Thetestcase"></a><h3>The testcase</h3>
 *
 *  This program essentially re-solves what we already do in  step-6  , i.e.
 * it solves the Laplace equation
 * @f{align*}
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -\Delta u &= f \qquad &&\text{in}\ \Omega=[0,1]^2, \\
 * u &= 0 \qquad &&\text{on}\ \partial\Omega.
 * @f}
 * The difference of course is now that we want to do so on a mesh thatmay
 * have a billion cells, with a billion or so degrees offreedom. There is no
 * doubt that doing so is completely silly for sucha simple problem, but the
 * point of a tutorial program is, after all,not to do something useful but to
 * show how useful programs can beimplemented using deal.II. Be that as it
 * may, to make things at leasta tiny bit interesting, we choose the right
 * hand side as adiscontinuous function,
 * @f{align*}
 * f(x,y)
 * =
 * \left\{
 * \begin{array}{ll}
 *  1 & \text{if}\ y > \frac 12 + \frac 14 \sin(4\pi x), \\
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
 * -1 & \text{otherwise},
 * \end{array}
 * \right.
 * @f}
 * so that the solution has a singularity along the sinusoidal linesnaking its
 * way through the domain. As a consequence, mesh refinementwill be
 * concentrated along this line. You can see this in the meshpicture shown
 * below in the results section.
 * Rather than
 * continuing here and giving a long introduction, let us gostraight to the
 * program code. If you have read through   step-6   and the  @ref distributed
 * documentation module, most of things that are goingto happen should be
 * familiar to you already. In fact, comparing the twoprograms you will notice
 * that the additional effort necessary to make thingswork in %parallel is
 * almost insignificant: the two programs have about thesame number of lines
 * of code (though   step-6   spends more space on dealing withcoefficients
 * and output). In either case, the comments below will only be onthe things
 * that set   step-40   apart from   step-6   and that aren't already
 * coveredin the   @ref distributed   documentation module.
 *
 *
 * @note   This program will be able to compute on as many processors as you
 * wantto throw at it, and for as large a problem as you have the memory and
 * patienceto solve. However, there <i>is</i> a limit: the number of unknowns
 * can notexceed the largest number that can be stored with an object of type
 * types::global_dof_index.   By default, this is an alias for
 * <code>unsignedint</code>, which on most machines today is a 32-bit integer,
 * limiting you tosome 4 billion (in reality, since this program uses PETSc,
 * you will be limitedto half that as PETSc uses signed integers). However,
 * this can be changedduring configuration to use 64-bit integers, see the
 * ReadMe file. This willgive problem sizes you are unlikely to exceed anytime
 * soon.
 *
 *  <a name="CommProg"></a> <h1> The commented program</h1> <a
 * name="Includefiles"></a>  <h3>Include files</h3>
 *
 *
 * Most of the include files we need for this program have already been
 * discussed in previous programs. In particular, all of the following should
 * already be familiar friends:
 *
 *
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 *
 * #include <deal.II/lac/generic_linear_algebra.h>
 *
 * @endcode
 *
 * This program can use either PETSc or Trilinos for its parallel algebra
 * needs. By default, if deal.II has been configured with PETSc, it will use
 * PETSc. Otherwise, the following few lines will check that deal.II has been
 * configured with Trilinos and take that.
 *
 *
 * But there may be cases where you want to use Trilinos, even though deal.II
 * hasalso* been configured with PETSc, for example to compare the performance
 * of these two libraries. To do this, add the following \#define to the
 * source code:   <div class=CodeFragmentInTutorialComment>
 *
 *
 * @code
 * #define FORCE_USE_OF_TRILINOS
 * @endcode
 *
 * </div>
 *
 *
 * Using this logic, the following lines will then import either the PETSc or
 * Trilinos wrappers into the namespace `LA` (for "linear algebra). In the
 * former case, we are also defining the macro `USE_PETSC_LA` so that we can
 * detect if we are using PETSc (see solve() for an example where this is
 * necessary).
 *
 *
 * @code
 * namespace LA
 * {
 * #if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
 * !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
 * using namespace dealii::LinearAlgebraPETSc;
 * #  define USE_PETSC_LA
 * #elif defined(DEAL_II_WITH_TRILINOS)
 * using namespace dealii::LinearAlgebraTrilinos;
 * #else
 * #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
 * #endif
 * } // namespace LA
 *
 *
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 *
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 *
 * @endcode
 *
 * The following, however, will be new or be used in new roles. Let's walk
 * through them. The first of these will provide the tools of the
 * Utilities::System   namespace that we will use to query things like the
 * number of processors associated with the current MPI universe, or the
 * number within this universe the processor this job runs on is:
 *
 *
 * @code
 * #include <deal.II/base/utilities.h>
 * @endcode
 *
 * The next one provides a class, ConditionOStream that allows us to write
 * code that would output things to a stream (such as   <code>std::cout</code>
 * on every processor but throws the text away on all but one of them. We
 * could achieve the same by simply putting an   <code>if</code>   statement
 * in front of each place where we may generate output, but this doesn't make
 * the code any prettier. In addition, the condition whether this processor
 * should or should not produce output to the screen is the same every time
 *
 *  -  and consequently it should be simple enough to put it into the statements that generate output itself.
 *
 *
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * @endcode
 *
 * After these
 * preliminaries, here is where it becomes more interesting. As mentioned in
 * the   @ref distributed   module, one of the fundamental truths of solving
 * problems on large numbers of processors is that there is no way for any
 * processor to store everything (e.g. information about all cells in the
 * mesh, all degrees of freedom, or the values of all elements of the solution
 * vector). Rather, every processor will <i>own</i> a few of each of these
 * and, if necessary, may <i>know</i> about a few more, for example the ones
 * that are located on cells adjacent to the ones this processor owns itself.
 * We typically call the latter <i>ghost cells</i>, <i>ghost nodes</i> or
 * <i>ghost elements of a vector</i>. The point of this discussion here is
 * that we need to have a way to indicate which elements a particular
 * processor owns or need to know of. This is the realm of the IndexSet class:
 * if there are a total of   $N$   cells, degrees of freedom, or vector
 * elements, associated with (non-negative) integral indices   $[0,N)$  , then
 * both the set of elements the current processor owns as well as the
 * (possibly larger) set of indices it needs to know about are subsets of the
 * set   $[0,N)$  . IndexSet is a class that stores subsets of this set in an
 * efficient format:
 *
 *
 * @code
 * #include <deal.II/base/index_set.h>
 * @endcode
 *
 * The next header file is necessary for a single function,
 * SparsityTools::distribute_sparsity_pattern.   The role of this function
 * will be explained below.
 *
 *
 * @code
 * #include <deal.II/lac/sparsity_tools.h>
 * @endcode
 *
 * The final two, new header files provide the class
 * parallel::distributed::Triangulation   that provides meshes distributed
 * across a potentially very large number of processors, while the second
 * provides the namespace   parallel::distributed::GridRefinement   that
 * offers functions that can adaptively refine such distributed meshes:
 *
 *
 * @code
 * #include <deal.II/distributed/tria.h>
 * #include <deal.II/distributed/grid_refinement.h>
 *
 * #include <fstream>
 * #include <iostream>
 *
 * namespace Step40
 * {
 * using namespace dealii;
 *
 * @endcode
 *
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a>  <h3>The
 * <code>LaplaceProblem</code> class template</h3>
 *
 *
 * Next let's declare the main class of this program. Its structure is almost
 * exactly that of the   step-6   tutorial program. The only significant
 * differences are:
 *
 *
 *
 *  - The   <code>mpi_communicator</code>   variable that describes the set of processors we want this code to run on. In practice, this will be MPI_COMM_WORLD, i.e. all processors the batch scheduling system has assigned to this particular job.
 *
 *
 *
 *  - The presence of the   <code>pcout</code>   variable of type ConditionOStream.
 *
 *
 *
 *  - The obvious use of   parallel::distributed::Triangulation   instead of Triangulation.
 *
 *
 *
 *  - The presence of two IndexSet objects that denote which sets of degrees of freedom (and associated elements of solution and right hand side vectors) we own on the current processor and which we need (as ghost elements) for the algorithms in this program to work.
 *
 *
 *
 *  - The fact that all matrices and vectors are now distributed. We use either the PETSc or Trilinos wrapper classes so that we can use one of the sophisticated preconditioners offered by Hypre (with PETSc) or ML (with Trilinos). Note that as part of this class, we store a solution vector that does not only contain the degrees of freedom the current processor owns, but also (as ghost elements) all those vector elements that correspond to "locally relevant" degrees of freedom (i.e. all those that live on locally owned cells or the layer of ghost cells that surround it).
 *
 *
 * @code
 * template <int dim>
 * class LaplaceProblem
 * {
 * public:
 *   LaplaceProblem();
 *
 *   void run();
 *
 * private:
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void refine_grid();
 *   void output_results(const unsigned int cycle) const;
 *
 *   MPI_Comm mpi_communicator;
 *
 *   parallel::distributed::Triangulation<dim> triangulation;
 *
 *   FE_Q<dim>       fe;
 *   DoFHandler<dim> dof_handler;
 *
 *   IndexSet locally_owned_dofs;
 *   IndexSet locally_relevant_dofs;
 *
 *   AffineConstraints<double> constraints;
 *
 *   LA::MPI::SparseMatrix system_matrix;
 *   LA::MPI::Vector       locally_relevant_solution;
 *   LA::MPI::Vector       system_rhs;
 *
 *   ConditionalOStream pcout;
 *   TimerOutput        computing_timer;
 * };
 *
 *
 * @endcode
 *
 * <a name="ThecodeLaplaceProblemcodeclassimplementation"></a>  <h3>The
 * <code>LaplaceProblem</code> class implementation</h3>
 *
 *
 * <a name="Constructor"></a>  <h4>Constructor</h4>
 *
 *
 * Constructors and destructors are rather trivial. In addition to what we do
 * in   step-6  , we set the set of processors we want to work on to all
 * machines available (MPI_COMM_WORLD); ask the triangulation to ensure that
 * the mesh remains smooth and free to refined islands, for example; and
 * initialize the   <code>pcout</code>   variable to only allow processor zero
 * to output anything. The final piece is to initialize a timer that we use to
 * determine how much compute time the different parts of the program take:
 *
 *
 * @code
 * template <int dim>
 * LaplaceProblem<dim>::LaplaceProblem()
 *   : mpi_communicator(MPI_COMM_WORLD)
 *   , triangulation(mpi_communicator,
 *                   typename Triangulation<dim>::MeshSmoothing(
 *                     Triangulation<dim>::smoothing_on_refinement |
 *                     Triangulation<dim>::smoothing_on_coarsening))
 *   , fe(2)
 *   , dof_handler(triangulation)
 *   , pcout(std::cout,
 *           (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *   , computing_timer(mpi_communicator,
 *                     pcout,
 *                     TimerOutput::summary,
 *                     TimerOutput::wall_times)
 * {}
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemsetup_system"></a>
 * <h4>LaplaceProblem::setup_system</h4>
 *
 *
 * The following function is, arguably, the most interesting one in the entire
 * program since it goes to the heart of what distinguishes %parallel
 * step-40   from sequential   step-6  . At the top we do what we always do:
 * tell the DoFHandler object to distribute degrees of freedom. Since the
 * triangulation we use here is distributed, the DoFHandler object is smart
 * enough to recognize that on each processor it can only distribute degrees
 * of freedom on cells it owns; this is followed by an exchange step in which
 * processors tell each other about degrees of freedom on ghost cell. The
 * result is a DoFHandler that knows about the degrees of freedom on locally
 * owned cells and ghost cells (i.e. cells adjacent to locally owned cells)
 * but nothing about cells that are further away, consistent with the basic
 * philosophy of distributed computing that no processor can know everything.
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::setup_system()
 * {
 *   TimerOutput::Scope t(computing_timer, "setup");
 *
 *   dof_handler.distribute_dofs(fe);
 *
 * @endcode
 *
 * The next two lines extract some information we will need later on, namely
 * two index sets that provide information about which degrees of freedom are
 * owned by the current processor (this information will be used to initialize
 * solution and right hand side vectors, and the system matrix, indicating
 * which elements to store on the current processor and which to expect to be
 * stored somewhere else); and an index set that indicates which degrees of
 * freedom are locally relevant (i.e. live on cells that the current processor
 * owns or on the layer of ghost cells around the locally owned cells; we need
 * all of these degrees of freedom, for example, to estimate the error on the
 * local cells).
 *
 *
 * @code
 *   locally_owned_dofs = dof_handler.locally_owned_dofs();
 *   DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *
 * @endcode
 *
 * Next, let us initialize the solution and right hand side vectors. As
 * mentioned above, the solution vector we seek does not only store elements
 * we own, but also ghost entries; on the other hand, the right hand side
 * vector only needs to have the entries the current processor owns since all
 * we will ever do is write into it, never read from it on locally owned cells
 * (of course the linear solvers will read from it, but they do not care about
 * the geometric location of degrees of freedom).
 *
 *
 * @code
 *   locally_relevant_solution.reinit(locally_owned_dofs,
 *                                    locally_relevant_dofs,
 *                                    mpi_communicator);
 *   system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *
 * @endcode
 *
 * The next step is to compute hanging node and boundary value constraints,
 * which we combine into a single object storing all constraints.
 * As with
 * all other things in %parallel, the mantra must be that no processor can
 * store all information about the entire universe. As a consequence, we need
 * to tell the AffineConstraints object for which degrees of freedom it can
 * store constraints and for which it may not expect any information to store.
 * In our case, as explained in the   @ref distributed   module, the degrees
 * of freedom we need to care about on each processor are the locally relevant
 * ones, so we pass this to the   AffineConstraints::reinit   function. As a
 * side note, if you forget to pass this argument, the AffineConstraints class
 * will allocate an array with length equal to the largest DoF index it has
 * seen so far. For processors with high MPI process number, this may be very
 * large
 *
 *  -  maybe on the order of billions. The program would then allocate more memory than for likely all other operations combined for this single array.
 *
 *
 * @code
 *   constraints.clear();
 *   constraints.reinit(locally_relevant_dofs);
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<dim>(),
 *                                            constraints);
 *   constraints.close();
 *
 * @endcode
 *
 * The last part of this function deals with initializing the matrix with accompanying sparsity pattern. As in previous tutorial programs, we use the DynamicSparsityPattern as an intermediate with which we then initialize the system matrix. To do so we have to tell the sparsity pattern its size but as above there is no way the resulting object will be able to store even a single pointer for each global degree of freedom; the best we can hope for is that it stores information about each locally relevant degree of freedom, i.e. all those that we may ever touch in the process of assembling the matrix (the   @ref distributed_paper   "distributed computing paper"
 * has a long discussion why one really needs the locally relevant, and not
 * the small set of locally active degrees of freedom in this context). So we
 * tell the sparsity pattern its size and what DoFs to store anything for and
 * then ask   DoFTools::make_sparsity_pattern   to fill it (this function
 * ignores all cells that are not locally owned, mimicking what we will do
 * below in the assembly process). After this, we call a function that
 * exchanges entries in these sparsity pattern between processors so that in
 * the end each processor really knows about all the entries that will exist
 * in that part of the finite element matrix that it will own. The final step
 * is to initialize the matrix with the sparsity pattern.
 *
 *
 * @code
 *   DynamicSparsityPattern dsp(locally_relevant_dofs);
 *
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 *   SparsityTools::distribute_sparsity_pattern(dsp,
 *                                              dof_handler.locally_owned_dofs(),
 *                                              mpi_communicator,
 *                                              locally_relevant_dofs);
 *
 *   system_matrix.reinit(locally_owned_dofs,
 *                        locally_owned_dofs,
 *                        dsp,
 *                        mpi_communicator);
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemassemble_system"></a>
 * <h4>LaplaceProblem::assemble_system</h4>
 *
 *
 * The function that then assembles the linear system is comparatively boring,
 * being almost exactly what we've seen before. The points to watch out for
 * are:
 *
 *
 *
 *  - Assembly must only loop over locally owned cells. There are multiple ways to test that; for example, we could compare a cell's subdomain_id against information from the triangulation as in <code>cell->subdomain_id() == triangulation.locally_owned_subdomain()</code>, or skip all cells for which the condition <code>cell->is_ghost() || cell->is_artificial()</code> is true. The simplest way, however, is to simply ask the cell whether it is owned by the local processor.
 *
 *
 *
 *  - Copying local contributions into the global matrix must include distributing constraints and boundary values. In other words, we cannot (as we did in   step-6  ) first copy every local contribution into the global matrix and only in a later step take care of hanging node constraints and boundary values. The reason is, as discussed in   step-17  , that the parallel vector classes do not provide access to arbitrary elements of the matrix once they have been assembled into it
 *
 *  -  in parts because they may simply no longer reside on the current processor but have instead been shipped to a different machine.
 *
 *
 *
 *  - The way we compute the right hand side (given the formula stated in the introduction) may not be the most elegant but will do for a program whose focus lies somewhere entirely different.
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::assemble_system()
 * {
 *   TimerOutput::Scope t(computing_timer, "assembly");
 *
 *   const QGauss<dim> quadrature_formula(fe.degree + 1);
 *
 *   FEValues<dim> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients |
 *                             update_quadrature_points | update_JxW_values);
 *
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *   const unsigned int n_q_points    = quadrature_formula.size();
 *
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 *
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 *
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     if (cell->is_locally_owned())
 *       {
 *         cell_matrix = 0.;
 *         cell_rhs    = 0.;
 *
 *         fe_values.reinit(cell);
 *
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           {
 *             const double rhs_value =
 *               (fe_values.quadrature_point(q_point)[1] >
 *                    0.5 +
 *                      0.25 std::sin(4.0 numbers::PI
 *                                      fe_values.quadrature_point(q_point)[0]) ?
 *                  1. :
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
 * -1.);
 *
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                   cell_matrix(i, j) += fe_values.shape_grad(i, q_point)
 *                                        fe_values.shape_grad(j, q_point)
 *                                        fe_values.JxW(q_point);
 *
 *                 cell_rhs(i) += rhs_value
 *                                fe_values.shape_value(i, q_point)
 *                                fe_values.JxW(q_point);
 *               }
 *           }
 *
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(cell_matrix,
 *                                                cell_rhs,
 *                                                local_dof_indices,
 *                                                system_matrix,
 *                                                system_rhs);
 *       }
 *
 * @endcode
 *
 *
 *
 * @code
 *   system_matrix.compress(VectorOperation::add);
 *   system_rhs.compress(VectorOperation::add);
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemsolve"></a>  <h4>LaplaceProblem::solve</h4>
 *
 *
 * Even though solving linear systems on potentially tens of thousands of
 * processors is by far not a trivial job, the function that does this is
 *
 *  -  at least at the outside
 *
 *  -  relatively simple. Most of the parts you've seen before. There are really only two things worth mentioning:
 *
 *
 *
 *  - Solvers and preconditioners are built on the deal.II wrappers of PETSc and Trilinos functionality. It is relatively well known that the primary bottleneck of massively %parallel linear solvers is not actually the communication between processors, but the fact that it is difficult to produce preconditioners that scale well to large numbers of processors. Over the second half of the first decade of the 21st century, it has become clear that algebraic multigrid (AMG) methods turn out to be extremely efficient in this context, and we will use one of them
 *
 *  -  either the BoomerAMG implementation of the Hypre package that can be interfaced to through PETSc, or a preconditioner provided by ML, which is part of Trilinos
 *
 *  -  for the current program. The rest of the solver itself is boilerplate and has been shown before. Since the linear system is symmetric and positive definite, we can use the CG method as the outer solver.
 *
 *
 *
 *  - Ultimately, we want a vector that stores not only the elements of the solution for degrees of freedom the current processor owns, but also all other locally relevant degrees of freedom. On the other hand, the solver itself needs a vector that is uniquely split between processors, without any overlap. We therefore create a vector at the beginning of this function that has these properties, use it to solve the linear system, and only assign it to the vector we want at the very end. This last step ensures that all ghost elements are also copied as necessary.
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::solve()
 * {
 *   TimerOutput::Scope t(computing_timer, "solve");
 *   LA::MPI::Vector    completely_distributed_solution(locally_owned_dofs,
 *                                                   mpi_communicator);
 *
 *   SolverControl solver_control(dof_handler.n_dofs(), 1e-12);
 *
 * #ifdef USE_PETSC_LA
 *   LA::SolverCG solver(solver_control, mpi_communicator);
 * #else
 *   LA::SolverCG solver(solver_control);
 * #endif
 *
 *   LA::MPI::PreconditionAMG preconditioner;
 *
 *   LA::MPI::PreconditionAMG::AdditionalData data;
 *
 * #ifdef USE_PETSC_LA
 *   data.symmetric_operator = true;
 * #else
 *   // Trilinos defaults are good
 * #endif
 *   preconditioner.initialize(system_matrix, data);
 *
 *   solver.solve(system_matrix,
 *                completely_distributed_solution,
 *                system_rhs,
 *                preconditioner);
 *
 *   pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *         << std::endl;
 *
 *   constraints.distribute(completely_distributed_solution);
 *
 *   locally_relevant_solution = completely_distributed_solution;
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemrefine_grid"></a>
 * <h4>LaplaceProblem::refine_grid</h4>
 *
 *
 * The function that estimates the error and refines the grid is again almost
 * exactly like the one in   step-6  . The only difference is that the
 * function that flags cells to be refined is now in namespace
 * parallel::distributed::GridRefinement
 *
 *  -  a namespace that has functions that can communicate between all involved processors and determine global thresholds to use in deciding which cells to refine and which to coarsen.
 * Note that we didn't have to do anything special about the
 * KellyErrorEstimator class: we just give it a vector with as many elements
 * as the local triangulation has cells (locally owned cells, ghost cells, and
 * artificial ones), but it only fills those entries that correspond to cells
 * that are locally owned.
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::refine_grid()
 * {
 *   TimerOutput::Scope t(computing_timer, "refine");
 *
 *   Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 *   KellyErrorEstimator<dim>::estimate(
 *     dof_handler,
 *     QGauss<dim
 *
 * - 1>(fe.degree + 1),
 *     std::map<types::boundary_id, const Function<dim>>(),
 *     locally_relevant_solution,
 *     estimated_error_per_cell);
 *   parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
 *     triangulation, estimated_error_per_cell, 0.3, 0.03);
 *   triangulation.execute_coarsening_and_refinement();
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemoutput_results"></a>
 * <h4>LaplaceProblem::output_results</h4>
 *
 *
 * Compared to the corresponding function in   step-6  , the one here is a tad
 * more complicated. There are two reasons: the first one is that we do not
 * just want to output the solution but also for each cell which processor
 * owns it (i.e. which "subdomain" it is in). Secondly, as discussed at length
 * in   step-17   and   step-18  , generating graphical data can be a
 * bottleneck in parallelizing. In   step-18  , we have moved this step out of
 * the actual computation but shifted it into a separate program that later
 * combined the output from various processors into a single file. But this
 * doesn't scale: if the number of processors is large, this may mean that the
 * step of combining data on a single processor later becomes the longest
 * running part of the program, or it may produce a file that's so large that
 * it can't be visualized any more. We here follow a more sensible approach,
 * namely creating individual files for each MPI process and leaving it to the
 * visualization program to make sense of that. To start, the top of the
 * function looks like it usually does. In addition to attaching the solution
 * vector (the one that has entries for all locally relevant, not only the
 * locally owned, elements), we attach a data vector that stores, for each
 * cell, the subdomain the cell belongs to. This is slightly tricky, because
 * of course not every processor knows about every cell. The vector we attach
 * therefore has an entry for every cell that the current processor has in its
 * mesh (locally owned ones, ghost cells, and artificial cells), but the
 * DataOut class will ignore all entries that correspond to cells that are not
 * owned by the current processor. As a consequence, it doesn't actually
 * matter what values we write into these vector entries: we simply fill the
 * entire vector with the number of the current MPI process (i.e. the
 * subdomain_id of the current process); this correctly sets the values we
 * care for, i.e. the entries that correspond to locally owned cells, while
 * providing the wrong value for all other elements
 *
 *  -  but these are then ignored anyway.
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
 * {
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(locally_relevant_solution, "u");
 *
 *   Vector<float> subdomain(triangulation.n_active_cells());
 *   for (unsigned int i = 0; i < subdomain.size(); ++i)
 *     subdomain(i) = triangulation.locally_owned_subdomain();
 *   data_out.add_data_vector(subdomain, "subdomain");
 *
 *   data_out.build_patches();
 *
 * @endcode
 *
 * The next step is to write this data to disk. We write up to 8 VTU files in
 * parallel with the help of MPI-IO. Additionally a PVTU record is generated,
 * which groups the written VTU files.
 *
 *
 * @code
 *   data_out.write_vtu_with_pvtu_record(
 *     "./", "solution", cycle, mpi_communicator, 2, 8);
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="LaplaceProblemrun"></a>  <h4>LaplaceProblem::run</h4>
 *
 *
 * The function that controls the overall behavior of the program is again
 * like the one in   step-6  . The minor difference are the use of
 * <code>pcout</code> instead of <code>std::cout</code>   for output to the
 * console (see also   step-17  ) and that we only generate graphical output
 * if at most 32 processors are involved. Without this limit, it would be just
 * too easy for people carelessly running this program without reading it
 * first to bring down the cluster interconnect and fill any file system
 * available :-) A functional difference to   step-6   is the use of a square
 * domain and that we start with a slightly finer mesh (5 global refinement
 * cycles)
 *
 *  -  there just isn't much of a point showing a massively %parallel program starting on 4 cells (although admittedly the point is only slightly stronger starting on 1024).
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::run()
 * {
 *   pcout << "Running with "
 * #ifdef USE_PETSC_LA
 *         << "PETSc"
 * #else
 *         << "Trilinos"
 * #endif
 *         << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
 *         << " MPI rank(s)..." << std::endl;
 *
 *   const unsigned int n_cycles = 8;
 *   for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *     {
 *       pcout << "Cycle " << cycle << ':' << std::endl;
 *
 *       if (cycle == 0)
 *         {
 *           GridGenerator::hyper_cube(triangulation);
 *           triangulation.refine_global(5);
 *         }
 *       else
 *         refine_grid();
 *
 *       setup_system();
 *
 *       pcout << "   Number of active cells:       "
 *             << triangulation.n_global_active_cells() << std::endl
 *             << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 *
 *       assemble_system();
 *       solve();
 *
 *       if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
 *         {
 *           TimerOutput::Scope t(computing_timer, "output");
 *           output_results(cycle);
 *         }
 *
 *       computing_timer.print_summary();
 *       computing_timer.reset();
 *
 *       pcout << std::endl;
 *     }
 * }
 * } // namespace Step40
 *
 *
 *
 * @endcode
 *
 * <a name="main"></a>  <h4>main()</h4>
 *
 *
 * The final function,   <code>main()</code>  , again has the same structure
 * as in all other programs, in particular   step-6  . Like the other programs
 * that use MPI, we have to initialize and finalize MPI, which is done using
 * the helper object   Utilities::MPI::MPI_InitFinalize.   The constructor of
 * that class also initializes libraries that depend on MPI, such as p4est,
 * PETSc, SLEPc, and Zoltan (though the last two are not used in this
 * tutorial). The order here is important: we cannot use any of these
 * libraries until they are initialized, so it does not make sense to do
 * anything before creating an instance of   Utilities::MPI::MPI_InitFinalize.
 *
 *
 * After the solver finishes, the LaplaceProblem destructor will run followed
 * by   Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize().   This order is
 * also important:   Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize()
 * calls   <code>PetscFinalize</code>   (and finalization functions for other
 * libraries), which will delete any in-use PETSc objects. This must be done
 * after we destruct the Laplace solver to avoid double deletion errors.
 * Fortunately, due to the order of destructor call rules of C++, we do not
 * need to worry about any of this: everything happens in the correct order
 * (i.e., the reverse of the order of construction). The last function called
 * by   Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize()   is
 * <code>MPI_Finalize</code>  : i.e., once this object is destructed the
 * program should exit since MPI will no longer be available.
 *
 *
 * @code
 * int main(int argc, charargv[])
 * {
 * try
 *   {
 *     using namespace dealii;
 *     using namespace Step40;
 *
 *     Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 *
 *     LaplaceProblem<2> laplace_problem_2d;
 *     laplace_problem_2d.run();
 *   }
 * catch (std::exception &exc)
 *   {
 *     std::cerr << std::endl
 *               << std::endl
 *               << "----------------------------------------------------"
 *               << std::endl;
 *     std::cerr << "Exception on processing: " << std::endl
 *               << exc.what() << std::endl
 *               << "Aborting!" << std::endl
 *               << "----------------------------------------------------"
 *               << std::endl;
 *
 *     return 1;
 *   }
 * catch (...)
 *   {
 *     std::cerr << std::endl
 *               << std::endl
 *               << "----------------------------------------------------"
 *               << std::endl;
 *     std::cerr << "Unknown exception!" << std::endl
 *               << "Aborting!" << std::endl
 *               << "----------------------------------------------------"
 *               << std::endl;
 *     return 1;
 *   }
 *
 * return 0;
 * }
 * @endcode
 * <a name="Results"></a><h1>Results</h1>
 *
 *  When you run the program, on a single processor or with your local
 * MPIinstallation on a few, you should get output like this:
 * @code
 * Cycle 0:
 * Number of active cells:       1024
 * Number of degrees of freedom: 4225
 * Solved in 10 iterations.
 *
 *
 *
 *
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |     0.176s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | assembly                        |         1 |    0.0209s |        12% |
 * | output                          |         1 |    0.0189s |        11% |
 * | setup                           |         1 |    0.0299s |        17% |
 * | solve                           |         1 |    0.0419s |        24% |
 * +---------------------------------+-----------+------------+------------+
 *
 *
 *
 *
 * Cycle 1:
 * Number of active cells:       1954
 * Number of degrees of freedom: 8399
 * Solved in 10 iterations.
 *
 *
 *
 *
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |     0.327s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | assembly                        |         1 |    0.0368s |        11% |
 * | output                          |         1 |    0.0208s |       6.4% |
 * | refine                          |         1 |     0.157s |        48% |
 * | setup                           |         1 |    0.0452s |        14% |
 * | solve                           |         1 |    0.0668s |        20% |
 * +---------------------------------+-----------+------------+------------+
 *
 *
 *
 *
 * Cycle 2:
 * Number of active cells:       3664
 * Number of degrees of freedom: 16183
 * Solved in 11 iterations.
 *
 * ...
 * @endcode
 *
 * The exact numbers differ, depending on how many processors we use;this is
 * due to the fact that the preconditioner depends on thepartitioning of the
 * problem, the solution then differs in the last fewdigits, and consequently
 * the mesh refinement differs slightly.The primary thing to notice here,
 * though, is that the number ofiterations does not increase with the size of
 * the problem. Thisguarantees that we can efficiently solve even the largest
 * problems. When run on a sufficiently large number of machines (say a
 * fewthousand), this program can relatively easily solve problems with
 * wellover one billion unknowns in less than a minute. On the other hand,such
 * big problems can no longer be visualized, so we also ran theprogram on only
 * 16 processors. Here are a mesh, along with itspartitioning onto the 16
 * processors, and the corresponding solution: <table width="100%"> <tr> <td>
 * <img src="https://www.dealii.org/images/steps/developer/step-40.mesh.png"
 * alt=""> </td> <td> <img
 * src="https://www.dealii.org/images/steps/developer/step-40.solution.png"
 * alt=""> </td> </tr> </table>
 * The mesh on the left has a mere 7,069 cells. This is of course aproblem we would easily have been able to solve already on a singleprocessor using   step-6  , but the point of the program was to show howto write a program that scales to many more machines. For example,here are two graphs that show how the run time of a large number of partsof the program scales on problems with around 52 and 375 million degrees offreedom if we take more and more processors (these and the next couple ofgraphs are taken from an earlier version of the  @ref distributed_paper   "Distributed Computing paper"
 * ; updated graphs showingdata of runs on even larger numbers of processors,
 * and a lotmore interpretation can be found in the final version of the
 * paper): <table width="100%"> <tr> <td> <img
 * src="https://www.dealii.org/images/steps/developer/step-40.strong2.png"
 * alt=""> </td> <td> <img
 * src="https://www.dealii.org/images/steps/developer/step-40.strong.png"
 * alt=""> </td> </tr> </table>
 * While the strong scaling graphs above show that we can solve a problem
 * offixed size faster and faster if we take more and more processors, the
 * moreinteresting question may be how big problems can become so that they
 * can stillbe solved within a reasonable time on a machine of a particular
 * size. We showthis in the following two graphs for 256 and 4096 processors:
 * <table width="100%"> <tr> <td> <img
 * src="https://www.dealii.org/images/steps/developer/step-40.256.png" alt="">
 * </td> <td> <img
 * src="https://www.dealii.org/images/steps/developer/step-40.4096.png"
 * alt=""> </td> </tr> </table>
 * What these graphs show is that all parts of the program scale linearly withthe number of degrees of freedom. This time, lines are wobbly at the left asthe size of local problems is too small. For more discussions of these resultswe refer to the   @ref distributed_paper   "Distributed Computing paper".
 * So how large are the largest problems one can solve? At the time of
 * writingthis problem, thelimiting factor is that the program uses the
 * BoomerAMG algebraicmultigrid method from the <a
 * href="http://acts.nersc.gov/hypre/" target="_top">Hypre package</a> asa
 * preconditioner, which unfortunately uses signed 32-bit integers toindex the
 * elements of a %distributed matrix. This limits the size ofproblems to
 * $2^{31}-1=2,147,483,647$   degrees of freedom. From the graphsabove it is
 * obvious that the scalability would extend beyond thisnumber, and one could
 * expect that given more than the 4,096 machinesshown above would also
 * further reduce the compute time. That said, onecan certainly expect that
 * this limit will eventually be lifted by thehypre developers. On the other
 * hand, this does not mean that deal.II cannot solve biggerproblems. Indeed,
 * step-37   shows how one can solve problems that are notjust a little, but
 * very substantially larger than anything we have shownhere.
 *
 *
 * <a name="extensions"></a><a
 * name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>
 *
 *  In a sense, this program is the ultimate solver for the Laplaceequation:
 * it can essentially solve the equation to whatever accuracyyou want, if only
 * you have enough processors available. Since theLaplace equation by itself
 * is not terribly interesting at this levelof accuracy, the more interesting
 * possibilities for extensiontherefore concern not so much this program but
 * what comes beyondit. For example, several of the other programs in this
 * tutorial havesignificant run times, especially in 3d. It would therefore
 * beinteresting to use the techniques explained here to extend otherprograms
 * to support parallel distributed computations. We have donethis for
 * step-31   in the   step-32   tutorial program, but the same wouldapply to,
 * for example,   step-23   and   step-25   for hyperbolic timedependent
 * problems,   step-33   for gas dynamics, or   step-35   for theNavier-Stokes
 * equations. Maybe equally interesting is the problem of postprocessing.
 * Asmentioned above, we only show pictures of the solution and the meshfor 16
 * processors because 4,096 processors solving 1 billion unknownswould produce
 * graphical output on the order of several 10gigabyte. Currently, no program
 * is able to visualize this amount ofdata in any reasonable way unless it
 * also runs on at least severalhundred processors. There are, however,
 * approaches where visualizationprograms directly communicate with solvers on
 * each processor with eachvisualization process rendering the part of the
 * scene computed by thesolver on this processor. Implementing such an
 * interface would allowto quickly visualize things that are otherwise not
 * amenable tographical display.
 *
* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-40.cc"
 *
 */
