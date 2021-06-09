/**
@page step_40 The step-40 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The <code>LaplaceProblem</code> class template</a>
        <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The <code>LaplaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblemrefine_grid">LaplaceProblem::refine_grid</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
        <li><a href="#main">main()</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-40/doc/intro.dox

 <br> 

<i>This program was contributed by Timo Heister, Martin Kronbichler and Wolfgang
Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation or of The
California Institute of Technology.
</i>




 @note  作为这个程序的前提条件，你需要同时安装PETSc和p4est库。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II与这两个附加库的安装。还要注意的是，为了正常工作，本程序需要访问实现代数多网格的Hypre预处理程序包；它可以作为PETSc的一部分安装，但必须在配置PETSc时明确启用；参见PETSc安装说明中的链接页面。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{41.5,41.75} 

鉴于今天的计算机，大多数有限元计算可以在一台机器上完成。因此，以前的大多数教程程序只显示了这一点，可能是在一些处理器之间进行分工，但这些处理器都可以访问相同的共享内存空间。也就是说，有些问题对于单台机器来说实在是太大了，在这种情况下，必须以适当的方式将问题分割给多台机器，每台机器都为整体贡献自己的一部分。在第17步和第18步中展示了一个简单的方法，我们展示了一个程序如何使用<a
href="http://www.mpi-forum.org/" target="_top">MPI</a>来并行组装线性系统，存储它，解决它，并计算误差估计。所有这些操作的扩展都是相对微不足道的（关于操作 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"），但是有一个明显的缺点：为了使这个实现适度简单，每个MPI处理器都必须保留自己的整个Triangulation和DoFHandler对象的副本。因此，虽然我们可以怀疑（有充分的理由）上面列出的操作可以扩展到成千上万的计算机和数十亿个单元和数十亿个自由度的问题规模，但在每一个最后的处理器上为这成千上万的计算机所解决的整个问题建立一个大的网格显然是不能扩展的：这将需要永远，也许更重要的是没有一台机器会有足够的内存来存储一个有十亿个单元的网格（至少在写这篇文章时没有）。在现实中，像第17步和第18步这样的程序不可能在超过100或200个处理器上运行，即使在那里，存储Triangulation和DoFHandler对象也会消耗每台机器上的绝大部分内存。

因此，我们需要以不同的方式来处理这个问题：为了扩展到非常大的问题，每个处理器只能存储自己的一小块三角形和DoFHandler对象。deal.II在 parallel::distributed 命名空间和其中的类中实现了这样一个方案。它建立在一个外部库上，<a
href="http://www.p4est.org/">p4est</a>（对表达式<i>parallel forest</i>的发挥，描述了将分层构造的网格作为四叉树或八叉树的森林进行并行存储）。你需要<a
href="../../external-libs/p4est.html">install and configure p4est</a>，但除此之外，它的所有工作原理都隐藏在deal.II的表面之下。

本质上， parallel::distributed::Triangulation 类和DoFHandler类中的代码所做的是分割全局网格，使每个处理器只存储其 "拥有 "的一小部分，以及围绕其拥有的单元的一层 "幽灵 "单元。在我们想要解决偏微分方程的领域的其余部分发生了什么，对每个处理器来说都是未知的，如果需要这些信息，只能通过与其他机器的交流来推断。这意味着我们还必须以不同于例如第17步和第18步的方式来思考问题：例如，没有一个处理器可以拥有用于后处理的整个解矢量，程序的每一部分都必须被并行化，因为没有一个处理器拥有顺序操作所需的所有信息。

在 @ref distributed 文档模块中描述了这种并行化如何发生的一般概述。在阅读本程序的源代码之前，你应该先阅读它，以获得一个顶层的概述。在 @ref distributed_paper "分布式计算论文 "中也提供了关于我们将在程序中使用的许多术语的简明讨论。也许值得一读，以了解本程序内部如何工作的背景信息。




<a name="Thetestcase"></a><h3>The testcase</h3>


这个程序基本上重新解决了我们在步骤6中已经做的事情，即它解决了拉普拉斯方程

@f{align*}


  -\Delta u &= f \qquad &&\text{in}\ \Omega=[0,1]^2, \\
  u &= 0 \qquad &&\text{on}\ \partial\Omega.


@f}

当然不同的是，现在我们要在一个可能有十亿个单元，有十亿个左右自由度的网格上这样做。毫无疑问，对于这样一个简单的问题，这样做是完全愚蠢的，但毕竟一个教程程序的重点不是做一些有用的东西，而是展示如何使用deal.II来实现有用的程序。尽管如此，为了使事情至少有一点点有趣，我们选择右侧为一个不连续的函数。

@f{align*}
  f(x,y)
  =
  \left\{
  \begin{array}{ll}
    1 & \text{if}\ y > \frac 12 + \frac 14 \sin(4\pi x), \\


    -1 & \text{otherwise},
  \end{array}
  \right.


@f}

使得解沿着蜿蜒穿过域的正弦线有一个奇点。因此，网格的细化将集中在这条线上。你可以在下面结果部分的网格图中看到这一点。

与其在这里继续做冗长的介绍，不如让我们直接进入程序代码。如果你已经读完了步骤6和 @ref distributed 文档模块，大部分将要发生的事情你应该已经熟悉了。事实上，比较这两个程序，你会发现在%parallel中工作所需的额外努力几乎是微不足道的：这两个程序的代码行数差不多（尽管步骤6在处理系数和输出方面花费了更多的空间）。在任何情况下，下面的评论将只针对使step-40与step-6不同的事情，而且在 @ref distributed 文档模块中还没有涵盖。




 @note  这个程序将能够在你想扔给它的多少个处理器上进行计算，以及你有多少内存和耐心来解决多大的问题。然而，<i>is</i>有一个限制：未知数的数量不能超过可以用类型 types::global_dof_index. 的对象存储的最大数量。默认情况下，这是<code>unsigned int</code>的别名，在今天大多数机器上是一个32位的整数，限制了你大约40亿（实际上，由于这个程序使用PETSc，你将被限制在一半，因为PETSc使用有符号整数）。然而，这可以在配置过程中改变为使用64位整数，见ReadMe文件。这将使问题的大小在短期内不太可能超过。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * Most of the include files we need for this program have already been
 * discussed in previous programs. In particular, all of the following should
 * already be familiar friends:
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
 * This program can use either PETSc or Trilinos for its parallel
 * algebra needs. By default, if deal.II has been configured with
 * PETSc, it will use PETSc. Otherwise, the following few lines will
 * check that deal.II has been configured with Trilinos and take that.
 * 

 * 
 * But there may be cases where you want to use Trilinos, even though
 * deal.II has *also* been configured with PETSc, for example to
 * compare the performance of these two libraries. To do this,
 * add the following \#define to the source code:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * #define FORCE_USE_OF_TRILINOS
 * @endcode
 * </div>
 * 

 * 
 * Using this logic, the following lines will then import either the
 * PETSc or Trilinos wrappers into the namespace `LA` (for "linear
 * algebra). In the former case, we are also defining the macro
 * `USE_PETSC_LA` so that we can detect if we are using PETSc (see
 * solve() for an example where this is necessary).
 * 
 * @code
 * namespace LA
 * {
 * #if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
 *   !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
 *   using namespace dealii::LinearAlgebraPETSc;
 * #  define USE_PETSC_LA
 * #elif defined(DEAL_II_WITH_TRILINOS)
 *   using namespace dealii::LinearAlgebraTrilinos;
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
 * Utilities::System namespace that we will use to query things like the
 * number of processors associated with the current MPI universe, or the
 * number within this universe the processor this job runs on is:
 * 
 * @code
 * #include <deal.II/base/utilities.h>
 * @endcode
 * 
 * The next one provides a class, ConditionOStream that allows us to write
 * code that would output things to a stream (such as <code>std::cout</code>
 * on every processor but throws the text away on all but one of them. We
 * could achieve the same by simply putting an <code>if</code> statement in
 * front of each place where we may generate output, but this doesn't make the
 * code any prettier. In addition, the condition whether this processor should
 * or should not produce output to the screen is the same every time -- and
 * consequently it should be simple enough to put it into the statements that
 * generate output itself.
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * @endcode
 * 
 * After these preliminaries, here is where it becomes more interesting. As
 * mentioned in the @ref distributed module, one of the fundamental truths of
 * solving problems on large numbers of processors is that there is no way for
 * any processor to store everything (e.g. information about all cells in the
 * mesh, all degrees of freedom, or the values of all elements of the solution
 * vector). Rather, every processor will <i>own</i> a few of each of these
 * and, if necessary, may <i>know</i> about a few more, for example the ones
 * that are located on cells adjacent to the ones this processor owns
 * itself. We typically call the latter <i>ghost cells</i>, <i>ghost nodes</i>
 * or <i>ghost elements of a vector</i>. The point of this discussion here is
 * that we need to have a way to indicate which elements a particular
 * processor owns or need to know of. This is the realm of the IndexSet class:
 * if there are a total of $N$ cells, degrees of freedom, or vector elements,
 * associated with (non-negative) integral indices $[0,N)$, then both the set
 * of elements the current processor owns as well as the (possibly larger) set
 * of indices it needs to know about are subsets of the set $[0,N)$. IndexSet
 * is a class that stores subsets of this set in an efficient format:
 * 
 * @code
 * #include <deal.II/base/index_set.h>
 * @endcode
 * 
 * The next header file is necessary for a single function,
 * SparsityTools::distribute_sparsity_pattern. The role of this function will
 * be explained below.
 * 
 * @code
 * #include <deal.II/lac/sparsity_tools.h>
 * @endcode
 * 
 * The final two, new header files provide the class
 * parallel::distributed::Triangulation that provides meshes distributed
 * across a potentially very large number of processors, while the second
 * provides the namespace parallel::distributed::GridRefinement that offers
 * functions that can adaptively refine such distributed meshes:
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
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceProblem</code> class template</h3>
 * 

 * 
 * Next let's declare the main class of this program. Its structure is
 * almost exactly that of the step-6 tutorial program. The only significant
 * differences are:
 * - The <code>mpi_communicator</code> variable that
 * describes the set of processors we want this code to run on. In practice,
 * this will be MPI_COMM_WORLD, i.e. all processors the batch scheduling
 * system has assigned to this particular job.
 * - The presence of the <code>pcout</code> variable of type ConditionOStream.
 * - The obvious use of parallel::distributed::Triangulation instead of
 * Triangulation.
 * - The presence of two IndexSet objects that denote which sets of degrees of
 * freedom (and associated elements of solution and right hand side vectors)
 * we own on the current processor and which we need (as ghost elements) for
 * the algorithms in this program to work.
 * - The fact that all matrices and vectors are now distributed. We use
 * either the PETSc or Trilinos wrapper classes so that we can use one of
 * the sophisticated preconditioners offered by Hypre (with PETSc) or ML
 * (with Trilinos). Note that as part of this class, we store a solution
 * vector that does not only contain the degrees of freedom the current
 * processor owns, but also (as ghost elements) all those vector elements
 * that correspond to "locally relevant" degrees of freedom (i.e. all
 * those that live on locally owned cells or the layer of ghost cells that
 * surround it).
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     MPI_Comm mpi_communicator;
 * 
 *     parallel::distributed::Triangulation<dim> triangulation;
 * 
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     LA::MPI::SparseMatrix system_matrix;
 *     LA::MPI::Vector       locally_relevant_solution;
 *     LA::MPI::Vector       system_rhs;
 * 
 *     ConditionalOStream pcout;
 *     TimerOutput        computing_timer;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>LaplaceProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * Constructors and destructors are rather trivial. In addition to what we
 * do in step-6, we set the set of processors we want to work on to all
 * machines available (MPI_COMM_WORLD); ask the triangulation to ensure that
 * the mesh remains smooth and free to refined islands, for example; and
 * initialize the <code>pcout</code> variable to only allow processor zero
 * to output anything. The final piece is to initialize a timer that we
 * use to determine how much compute time the different parts of the program
 * take:
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem()
 *     : mpi_communicator(MPI_COMM_WORLD)
 *     , triangulation(mpi_communicator,
 *                     typename Triangulation<dim>::MeshSmoothing(
 *                       Triangulation<dim>::smoothing_on_refinement |
 *                       Triangulation<dim>::smoothing_on_coarsening))
 *     , fe(2)
 *     , dof_handler(triangulation)
 *     , pcout(std::cout,
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *     , computing_timer(mpi_communicator,
 *                       pcout,
 *                       TimerOutput::summary,
 *                       TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * The following function is, arguably, the most interesting one in the
 * entire program since it goes to the heart of what distinguishes %parallel
 * step-40 from sequential step-6.
 *   

 * 
 * At the top we do what we always do: tell the DoFHandler object to
 * distribute degrees of freedom. Since the triangulation we use here is
 * distributed, the DoFHandler object is smart enough to recognize that on
 * each processor it can only distribute degrees of freedom on cells it
 * owns; this is followed by an exchange step in which processors tell each
 * other about degrees of freedom on ghost cell. The result is a DoFHandler
 * that knows about the degrees of freedom on locally owned cells and ghost
 * cells (i.e. cells adjacent to locally owned cells) but nothing about
 * cells that are further away, consistent with the basic philosophy of
 * distributed computing that no processor can know everything.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "setup");
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * The next two lines extract some information we will need later on,
 * namely two index sets that provide information about which degrees of
 * freedom are owned by the current processor (this information will be
 * used to initialize solution and right hand side vectors, and the system
 * matrix, indicating which elements to store on the current processor and
 * which to expect to be stored somewhere else); and an index set that
 * indicates which degrees of freedom are locally relevant (i.e. live on
 * cells that the current processor owns or on the layer of ghost cells
 * around the locally owned cells; we need all of these degrees of
 * freedom, for example, to estimate the error on the local cells).
 * 
 * @code
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 * @endcode
 * 
 * Next, let us initialize the solution and right hand side vectors. As
 * mentioned above, the solution vector we seek does not only store
 * elements we own, but also ghost entries; on the other hand, the right
 * hand side vector only needs to have the entries the current processor
 * owns since all we will ever do is write into it, never read from it on
 * locally owned cells (of course the linear solvers will read from it,
 * but they do not care about the geometric location of degrees of
 * freedom).
 * 
 * @code
 *     locally_relevant_solution.reinit(locally_owned_dofs,
 *                                      locally_relevant_dofs,
 *                                      mpi_communicator);
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 * 
 * @endcode
 * 
 * The next step is to compute hanging node and boundary value
 * constraints, which we combine into a single object storing all
 * constraints.
 *     

 * 
 * As with all other things in %parallel, the mantra must be that no
 * processor can store all information about the entire universe. As a
 * consequence, we need to tell the AffineConstraints object for which
 * degrees of freedom it can store constraints and for which it may not
 * expect any information to store. In our case, as explained in the
 * @ref distributed module, the degrees of freedom we need to care about on
 * each processor are the locally relevant ones, so we pass this to the
 * AffineConstraints::reinit function. As a side note, if you forget to
 * pass this argument, the AffineConstraints class will allocate an array
 * with length equal to the largest DoF index it has seen so far. For
 * processors with high MPI process number, this may be very large --
 * maybe on the order of billions. The program would then allocate more
 * memory than for likely all other operations combined for this single
 * array.
 * 
 * @code
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(),
 *                                              constraints);
 *     constraints.close();
 * 
 * @endcode
 * 
 * The last part of this function deals with initializing the matrix with
 * accompanying sparsity pattern. As in previous tutorial programs, we use
 * the DynamicSparsityPattern as an intermediate with which we
 * then initialize the system matrix. To do so we have to tell the sparsity
 * pattern its size but as above there is no way the resulting object will
 * be able to store even a single pointer for each global degree of
 * freedom; the best we can hope for is that it stores information about
 * each locally relevant degree of freedom, i.e. all those that we may
 * ever touch in the process of assembling the matrix (the
 * @ref distributed_paper "distributed computing paper" has a long
 * discussion why one really needs the locally relevant, and not the small
 * set of locally active degrees of freedom in this context).
 *     

 * 
 * So we tell the sparsity pattern its size and what DoFs to store
 * anything for and then ask DoFTools::make_sparsity_pattern to fill it
 * (this function ignores all cells that are not locally owned, mimicking
 * what we will do below in the assembly process). After this, we call a
 * function that exchanges entries in these sparsity pattern between
 * processors so that in the end each processor really knows about all the
 * entries that will exist in that part of the finite element matrix that
 * it will own. The final step is to initialize the matrix with the
 * sparsity pattern.
 * 
 * @code
 *     DynamicSparsityPattern dsp(locally_relevant_dofs);
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 *     SparsityTools::distribute_sparsity_pattern(dsp,
 *                                                dof_handler.locally_owned_dofs(),
 *                                                mpi_communicator,
 *                                                locally_relevant_dofs);
 * 
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          dsp,
 *                          mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system</h4>
 * 

 * 
 * The function that then assembles the linear system is comparatively
 * boring, being almost exactly what we've seen before. The points to watch
 * out for are:
 * - Assembly must only loop over locally owned cells. There
 * are multiple ways to test that; for example, we could compare a cell's
 * subdomain_id against information from the triangulation as in
 * <code>cell->subdomain_id() ==
 * triangulation.locally_owned_subdomain()</code>, or skip all cells for
 * which the condition <code>cell->is_ghost() ||
 * cell->is_artificial()</code> is true. The simplest way, however, is to
 * simply ask the cell whether it is owned by the local processor.
 * - Copying local contributions into the global matrix must include
 * distributing constraints and boundary values. In other words, we cannot
 * (as we did in step-6) first copy every local contribution into the global
 * matrix and only in a later step take care of hanging node constraints and
 * boundary values. The reason is, as discussed in step-17, that the
 * parallel vector classes do not provide access to arbitrary elements of
 * the matrix once they have been assembled into it -- in parts because they
 * may simply no longer reside on the current processor but have instead
 * been shipped to a different machine.
 * - The way we compute the right hand side (given the
 * formula stated in the introduction) may not be the most elegant but will
 * do for a program whose focus lies somewhere entirely different.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "assembly");
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix = 0.;
 *           cell_rhs    = 0.;
 * 
 *           fe_values.reinit(cell);
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             {
 *               const double rhs_value =
 *                 (fe_values.quadrature_point(q_point)[1] >
 *                      0.5 +
 *                        0.25 * std::sin(4.0 * numbers::PI *
 *                                        fe_values.quadrature_point(q_point)[0]) ?
 *                    1. :
 *                    -1.);
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
 *                                          fe_values.shape_grad(j, q_point) *
 *                                          fe_values.JxW(q_point);
 * 
 *                   cell_rhs(i) += rhs_value *                         
 *                                  fe_values.shape_value(i, q_point) * 
 *                                  fe_values.JxW(q_point);
 *                 }
 *             }
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints.distribute_local_to_global(cell_matrix,
 *                                                  cell_rhs,
 *                                                  local_dof_indices,
 *                                                  system_matrix,
 *                                                  system_rhs);
 *         }
 * 
 * @endcode
 * 
 * Notice that the assembling above is just a local operation. So, to
 * form the "global" linear system, a synchronization between all
 * processors is needed. This could be done by invoking the function
 * compress(). See @ref GlossCompress "Compressing distributed objects"
 * for more information on what is compress() designed to do.
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * Even though solving linear systems on potentially tens of thousands of
 * processors is by far not a trivial job, the function that does this is --
 * at least at the outside -- relatively simple. Most of the parts you've
 * seen before. There are really only two things worth mentioning:
 * - Solvers and preconditioners are built on the deal.II wrappers of PETSc
 * and Trilinos functionality. It is relatively well known that the
 * primary bottleneck of massively %parallel linear solvers is not
 * actually the communication between processors, but the fact that it is
 * difficult to produce preconditioners that scale well to large numbers
 * of processors. Over the second half of the first decade of the 21st
 * century, it has become clear that algebraic multigrid (AMG) methods
 * turn out to be extremely efficient in this context, and we will use one
 * of them -- either the BoomerAMG implementation of the Hypre package
 * that can be interfaced to through PETSc, or a preconditioner provided
 * by ML, which is part of Trilinos -- for the current program. The rest
 * of the solver itself is boilerplate and has been shown before. Since
 * the linear system is symmetric and positive definite, we can use the CG
 * method as the outer solver.
 * - Ultimately, we want a vector that stores not only the elements
 * of the solution for degrees of freedom the current processor owns, but
 * also all other locally relevant degrees of freedom. On the other hand,
 * the solver itself needs a vector that is uniquely split between
 * processors, without any overlap. We therefore create a vector at the
 * beginning of this function that has these properties, use it to solve the
 * linear system, and only assign it to the vector we want at the very
 * end. This last step ensures that all ghost elements are also copied as
 * necessary.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     TimerOutput::Scope t(computing_timer, "solve");
 *     LA::MPI::Vector    completely_distributed_solution(locally_owned_dofs,
 *                                                     mpi_communicator);
 * 
 *     SolverControl solver_control(dof_handler.n_dofs(), 1e-12);
 * 
 * #ifdef USE_PETSC_LA
 *     LA::SolverCG solver(solver_control, mpi_communicator);
 * #else
 *     LA::SolverCG solver(solver_control);
 * #endif
 * 
 *     LA::MPI::PreconditionAMG preconditioner;
 * 
 *     LA::MPI::PreconditionAMG::AdditionalData data;
 * 
 * #ifdef USE_PETSC_LA
 *     data.symmetric_operator = true;
 * #else
 *     /* Trilinos defaults are good */
 * #endif
 *     preconditioner.initialize(system_matrix, data);
 * 
 *     solver.solve(system_matrix,
 *                  completely_distributed_solution,
 *                  system_rhs,
 *                  preconditioner);
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *           << std::endl;
 * 
 *     constraints.distribute(completely_distributed_solution);
 * 
 *     locally_relevant_solution = completely_distributed_solution;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrefine_grid"></a> 
 * <h4>LaplaceProblem::refine_grid</h4>
 * 

 * 
 * The function that estimates the error and refines the grid is again
 * almost exactly like the one in step-6. The only difference is that the
 * function that flags cells to be refined is now in namespace
 * parallel::distributed::GridRefinement -- a namespace that has functions
 * that can communicate between all involved processors and determine global
 * thresholds to use in deciding which cells to refine and which to coarsen.
 *   

 * 
 * Note that we didn't have to do anything special about the
 * KellyErrorEstimator class: we just give it a vector with as many elements
 * as the local triangulation has cells (locally owned cells, ghost cells,
 * and artificial ones), but it only fills those entries that correspond to
 * cells that are locally owned.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::refine_grid()
 *   {
 *     TimerOutput::Scope t(computing_timer, "refine");
 * 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       locally_relevant_solution,
 *       estimated_error_per_cell);
 *     parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
 *       triangulation, estimated_error_per_cell, 0.3, 0.03);
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemoutput_results"></a> 
 * <h4>LaplaceProblem::output_results</h4>
 * 

 * 
 * Compared to the corresponding function in step-6, the one here is a tad
 * more complicated. There are two reasons: the first one is that we do not
 * just want to output the solution but also for each cell which processor
 * owns it (i.e. which "subdomain" it is in). Secondly, as discussed at
 * length in step-17 and step-18, generating graphical data can be a
 * bottleneck in parallelizing. In step-18, we have moved this step out of
 * the actual computation but shifted it into a separate program that later
 * combined the output from various processors into a single file. But this
 * doesn't scale: if the number of processors is large, this may mean that
 * the step of combining data on a single processor later becomes the
 * longest running part of the program, or it may produce a file that's so
 * large that it can't be visualized any more. We here follow a more
 * sensible approach, namely creating individual files for each MPI process
 * and leaving it to the visualization program to make sense of that.
 *   

 * 
 * To start, the top of the function looks like it usually does. In addition
 * to attaching the solution vector (the one that has entries for all locally
 * relevant, not only the locally owned, elements), we attach a data vector
 * that stores, for each cell, the subdomain the cell belongs to. This is
 * slightly tricky, because of course not every processor knows about every
 * cell. The vector we attach therefore has an entry for every cell that the
 * current processor has in its mesh (locally owned ones, ghost cells, and
 * artificial cells), but the DataOut class will ignore all entries that
 * correspond to cells that are not owned by the current processor. As a
 * consequence, it doesn't actually matter what values we write into these
 * vector entries: we simply fill the entire vector with the number of the
 * current MPI process (i.e. the subdomain_id of the current process); this
 * correctly sets the values we care for, i.e. the entries that correspond
 * to locally owned cells, while providing the wrong value for all other
 * elements -- but these are then ignored anyway.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(locally_relevant_solution, "u");
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells());
 *     for (unsigned int i = 0; i < subdomain.size(); ++i)
 *       subdomain(i) = triangulation.locally_owned_subdomain();
 *     data_out.add_data_vector(subdomain, "subdomain");
 * 
 *     data_out.build_patches();
 * 
 * @endcode
 * 
 * The next step is to write this data to disk. We write up to 8 VTU files
 * in parallel with the help of MPI-IO. Additionally a PVTU record is
 * generated, which groups the written VTU files.
 * 
 * @code
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", cycle, mpi_communicator, 2, 8);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * The function that controls the overall behavior of the program is again
 * like the one in step-6. The minor difference are the use of
 * <code>pcout</code> instead of <code>std::cout</code> for output to the
 * console (see also step-17) and that we only generate graphical output if
 * at most 32 processors are involved. Without this limit, it would be just
 * too easy for people carelessly running this program without reading it
 * first to bring down the cluster interconnect and fill any file system
 * available :-)
 *   

 * 
 * A functional difference to step-6 is the use of a square domain and that
 * we start with a slightly finer mesh (5 global refinement cycles) -- there
 * just isn't much of a point showing a massively %parallel program starting
 * on 4 cells (although admittedly the point is only slightly stronger
 * starting on 1024).
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     pcout << "Running with "
 * #ifdef USE_PETSC_LA
 *           << "PETSc"
 * #else
 *           << "Trilinos"
 * #endif
 *           << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator)
 *           << " MPI rank(s)..." << std::endl;
 * 
 *     const unsigned int n_cycles = 8;
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation);
 *             triangulation.refine_global(5);
 *           }
 *         else
 *           refine_grid();
 * 
 *         setup_system();
 * 
 *         pcout << "   Number of active cells:       "
 *               << triangulation.n_global_active_cells() << std::endl
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 *         assemble_system();
 *         solve();
 * 
 *         if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
 *           {
 *             TimerOutput::Scope t(computing_timer, "output");
 *             output_results(cycle);
 *           }
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 * 
 *         pcout << std::endl;
 *       }
 *   }
 * } // namespace Step40
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="main"></a> 
 * <h4>main()</h4>
 * 

 * 
 * The final function, <code>main()</code>, again has the same structure as in
 * all other programs, in particular step-6. Like the other programs that use
 * MPI, we have to initialize and finalize MPI, which is done using the helper
 * object Utilities::MPI::MPI_InitFinalize. The constructor of that class also
 * initializes libraries that depend on MPI, such as p4est, PETSc, SLEPc, and
 * Zoltan (though the last two are not used in this tutorial). The order here
 * is important: we cannot use any of these libraries until they are
 * initialized, so it does not make sense to do anything before creating an
 * instance of Utilities::MPI::MPI_InitFinalize.
 * 

 * 
 * After the solver finishes, the LaplaceProblem destructor will run followed
 * by Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize(). This order is
 * also important: Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize() calls
 * <code>PetscFinalize</code> (and finalization functions for other
 * libraries), which will delete any in-use PETSc objects. This must be done
 * after we destruct the Laplace solver to avoid double deletion
 * errors. Fortunately, due to the order of destructor call rules of C++, we
 * do not need to worry about any of this: everything happens in the correct
 * order (i.e., the reverse of the order of construction). The last function
 * called by Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize() is
 * <code>MPI_Finalize</code>: i.e., once this object is destructed the program
 * should exit since MPI will no longer be available.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step40;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       LaplaceProblem<2> laplace_problem_2d;
 *       laplace_problem_2d.run();
 *     }
 *   catch (std::exception &exc)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Exception on processing: " << std::endl
 *                 << exc.what() << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 * 
 *       return 1;
 *     }
 *   catch (...)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Unknown exception!" << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       return 1;
 *     }
 * 
 *   return 0;
 * }
 * @endcode
examples/step-40/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你在单个处理器上或在几个本地MPI安装上运行该程序时，你应该得到这样的输出。

@code
Cycle 0:
   Number of active cells:       1024
   Number of degrees of freedom: 4225
   Solved in 10 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.176s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| assembly                        |         1 |    0.0209s |        12% |
| output                          |         1 |    0.0189s |        11% |
| setup                           |         1 |    0.0299s |        17% |
| solve                           |         1 |    0.0419s |        24% |
+---------------------------------+-----------+------------+------------+



Cycle 1:
   Number of active cells:       1954
   Number of degrees of freedom: 8399
   Solved in 10 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.327s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| assembly                        |         1 |    0.0368s |        11% |
| output                          |         1 |    0.0208s |       6.4% |
| refine                          |         1 |     0.157s |        48% |
| setup                           |         1 |    0.0452s |        14% |
| solve                           |         1 |    0.0668s |        20% |
+---------------------------------+-----------+------------+------------+



Cycle 2:
   Number of active cells:       3664
   Number of degrees of freedom: 16183
   Solved in 11 iterations.


...
@endcode



确切的数字是不同的，这取决于我们使用多少个处理器；这是由于预处理程序取决于问题的分区，然后解决方案在最后几位上有所不同，因此，网格细化也略有不同。不过，这里最值得注意的是，迭代次数并不随问题的大小而增加。这保证了我们甚至可以有效地解决最大的问题。

当在足够多的机器上运行时（比如说几千台），这个程序可以相对容易地在不到一分钟的时间内解决有远超过10亿个未知数的问题。另一方面，这样的大问题已经不能被视觉化，所以我们也只在16个处理器上运行该程序。下面是一个网格，以及它在16个处理器上的划分，还有相应的解决方案。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.mesh.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.solution.png" alt="">
</td>
</tr>
</table> 

左边的网格仅有7,069个单元。当然，这个问题我们在单台处理器上使用step-6就已经很容易解决了，但是这个程序的重点是展示如何编写一个可以扩展到更多机器的程序。例如，这里有两张图，显示了如果我们采取越来越多的处理器，程序的大量部分的运行时间是如何在大约5200万和37500万自由度的问题上扩展的（这些和接下来的几张图取自 @ref distributed_paper "分布式计算论文 "的早期版本；显示在更大数量的处理器上运行数据的更新图，以及更多的解释可以在该论文的最终版本中找到）。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.strong2.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.strong.png" alt="">
</td>
</tr>
</table> 

可以清楚地看到，这个程序可以很好地扩展到非常多的处理器。关于我们认为的 "可扩展 "程序的讨论，见 @ref GlossParallelScaling "本词汇表条目"）。曲线，特别是线性求解器，在图形的右端变得有点摇摆不定，因为每个处理器要做的事情太少，无法抵消通信成本（在上面两个例子中，每个处理器要解决的整个问题的部分，在使用4,096个处理器时，只有13,000和90,000个自由度；一个好的经验法则是，如果每个处理器至少有100,000个未知数，并行程序就会运行良好）。

虽然上面的强扩展图显示，如果我们采取越来越多的处理器，我们可以越来越快地解决一个固定大小的问题，但更有趣的问题可能是，问题可以变得多大，以便在一个特定大小的机器上仍然可以在合理的时间内解决它们。我们在下面两张256和4096处理器的图中展示了这一点。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.256.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.4096.png" alt="">
</td>
</tr>
</table> 

这些图显示的是，程序的所有部分都随着自由度数的增加而线性扩展。这一次，由于局部问题的规模太小，线条在左边摇摆不定。关于这些结果的更多讨论，我们参考了 @ref distributed_paper "分布式计算论文"。

那么，一个人能够解决的最大问题是多大？在写这个问题的时候，限制因素是程序使用<a
href="http://acts.nersc.gov/hypre/" target="_top">Hypre package</a>中的BoomerAMG代数多网格方法作为预处理程序，不幸的是，它使用有符号的32位整数来索引%分布式矩阵的元素。这将问题的大小限制在 $2^{31}-1=2,147,483,647$ 个自由度。从上面的图中可以看出，可扩展性会超过这个数字，而且可以预期，给定超过上面显示的4096台机器也会进一步减少计算时间。也就是说，人们当然可以期待，这个限制最终会被hybre的开发者解除。

另一方面，这并不意味着deal.II不能解决更大的问题。事实上，step-37展示了如何解决不仅仅是一点点，而是大大超过我们在这里所展示的任何问题的问题。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


从某种意义上说，这个程序是拉普拉斯方程的终极解算器：只要你有足够的处理器，它基本上可以把方程解到你想要的精度。由于拉普拉斯方程本身在这种精度水平上并不十分有趣，因此，更有趣的扩展可能性不在于这个程序，而在于它之后的内容。例如，本教程中的其他几个程序都有相当长的运行时间，特别是在3D中。因此，使用这里解释的技术来扩展其他程序以支持并行的分布式计算将是有趣的。我们在step-32教程程序中对step-31做了这样的处理，但同样的做法也适用于，例如，用于双曲时间相关问题的step-23和step-25，用于气体动力学的step-33，或用于纳维-斯托克斯方程的step-35。

也许同样有趣的是后处理的问题。如上所述，我们只展示了16个处理器的解决方案和网格的图片，因为4,096个处理器解决10亿个未知数会产生几10G的图形输出。目前，除非在至少几百个处理器上运行，否则没有任何程序能够以任何合理的方式将如此大量的数据可视化。然而，有一些方法，可视化程序直接与每个处理器上的求解器进行通信，每个可视化进程渲染这个处理器上的求解器所计算的场景部分。实现这样的接口将允许快速可视化那些在其他方面不适合用图形显示的东西。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-40.cc"
*/
