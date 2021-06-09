/**
@page step_17 The step-17 tutorial program
This tutorial depends on step-8.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Overview">Overview</a>
        <li><a href="#ParallelizingsoftwarewithMPI">Parallelizing software with MPI</a>
        <li><a href="#Whatthisprogramdoes">What this program does</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeElasticProblemcodeclasstemplate">The <code>ElasticProblem</code> class template</a>
        <li><a href="#Righthandsidevalues">Right hand side values</a>
        <li><a href="#ThecodeElasticProblemcodeclassimplementation">The <code>ElasticProblem</code> class implementation</a>
      <ul>
        <li><a href="#ElasticProblemElasticProblem">ElasticProblem::ElasticProblem</a>
        <li><a href="#ElasticProblemsetup_system">ElasticProblem::setup_system</a>
        <li><a href="#ElasticProblemassemble_system">ElasticProblem::assemble_system</a>
        <li><a href="#ElasticProblemsolve">ElasticProblem::solve</a>
        <li><a href="#ElasticProblemrefine_grid">ElasticProblem::refine_grid</a>
        <li><a href="#ElasticProblemoutput_results">ElasticProblem::output_results</a>
        <li><a href="#ElasticProblemrun">ElasticProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-17/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Overview"></a><h3>Overview</h3>


这个程序没有引入任何新的数学思想；事实上，它所做的只是做与step-8已经做的完全相同的计算，但它以一种不同的方式来做：我们没有使用deal.II自己的线性代数类，而是在deal.II提供的类之上建立一切，这些类包裹着<a
href="http://www.mcs.anl.gov/petsc/" target="_top">PETSc</a>库的线性代数实现。由于PETSc允许将矩阵和向量分布在MPI网络中的几台计算机上，因此产生的代码甚至能够以%并行方式解决问题。如果你不知道PETSc是什么，那么这将是一个快速浏览其主页的好时机。

作为这个程序的先决条件，你需要安装PETSc，如果你想在一个集群上以%并行方式运行，你还需要<a
href="http://www-users.cs.umn.edu/~karypis/metis/index.html"
target="_top">METIS</a>来划分网格。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II和这两个附加库的安装。

现在，关于细节：如前所述，该程序不计算任何新的东西，所以对有限元类等的使用与以前完全相同。与以前的程序不同的是，我们用几乎所有的类 <code>Vector</code> and <code>SparseMatrix</code> 代替了它们的近似值 <code>PETScWrappers::MPI::Vector</code> 和 <code>PETScWrappers::MPI::SparseMatrix</code> ，它们存储数据的方式使MPI网络中的每个处理器只存储矩阵或矢量的一部分。更具体地说，每个处理器将只存储与它 "拥有 "的自由度相对应的矩阵的那些行。对于向量，它们要么只存储与处理器拥有的自由度相对应的元素（这是右手边所需要的），要么也存储一些额外的元素，以确保每个处理器都能访问处理器拥有的单元（所谓 @ref GlossLocallyActiveDof "本地活动的自由度"）或邻近单元（所谓 @ref GlossLocallyRelevantDof "本地相关自由度"）上的解组件。

来自PETScWrapper命名空间的类所提供的接口与deal.II线性代数类的接口非常相似，但它们不是自己实现这一功能，而是简单地传递给它们相应的PETSc函数。因此，包装器只是用来给PETSc一个更现代的、面向对象的接口，并使PETSc和deal.II对象的使用尽可能地互换。使用PETSc的主要意义在于它可以在%并行状态下运行。我们将利用这一点，将域划分为与MPI网络中的进程一样多的块（"子域"）。同时，PETSc还提供了假的MPI存根，所以如果PETSc的配置中没有MPI，你可以在一台机器上运行这个程序。




<a name="ParallelizingsoftwarewithMPI"></a><h3>Parallelizing software with MPI</h3>


开发软件以通过MPI在%parallel中运行，需要改变一下思维方式，因为我们通常必须分割所有的数据结构，使每个处理器只存储整个问题的一部分。因此，你通常不能在每个处理器上访问一个解决方案向量的所有组成部分 -- 每个处理器可能根本没有足够的内存来容纳整个解决方案向量。由于数据被分割或 "分布 "在各个处理器上，我们把MPI使用的编程模型称为 "分布式内存计算"（与 "共享内存计算 "相反，后者意味着多个处理器可以访问一个内存空间中的所有数据，例如，当一台机器的多个核心在一个共同任务上工作时）。分布式内存计算的一些基本原理在 @ref distributed "使用分布式内存的多处理器并行计算 "文档模块中讨论，该模块本身是 @ref Parallel "并行计算 "模块的一个子模块。

一般来说，为了真正能够扩展到大量的处理器，我们需要在可用的处理器之间分割出<i>every</i>数据结构，其大小随着整个问题的大小而扩展。关于程序 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"）。这包括，例如，三角形、矩阵和所有全局向量（解决方案，右手边）。如果不拆分所有这些对象，其中一个对象将被复制到所有的处理器上，如果问题大小（和可用的处理器数量）变得很大，最终会简单地变得太大。另一方面，在每个处理器上保留大小与整个问题大小无关的对象是完全可以的。例如，可执行文件的每个副本将创建自己的有限元对象，或者我们在汇编中使用的局部矩阵）。)

在当前的程序中（以及相关的第18步），我们不会走得这么远，而是对MPI的使用做一个比较温和的介绍。更具体地说，我们要并行化的数据结构只有矩阵和向量。然而，我们并没有拆分Triangulation和DoFHandler类：每个进程仍然拥有这些对象的完整副本，而且所有进程都拥有其他进程所拥有的确切副本。然后，我们只需在每个处理器上的三角形的每个副本中，标记哪个处理器拥有哪些单元。这个过程被称为将网格 "分割 "为 @ref GlossSubdomainId "子域"。

对于较大的问题，必须在每个处理器上存储<i>entire</i>网格，显然会产生一个瓶颈。分割网格是稍微的，虽然没有多复杂（从用户的角度来看，虽然它<i>much</i>下更复杂）来实现，我们将展示如何在step-40和其他一些程序中这样做。在讨论这个程序的某个功能如何工作的过程中，我们会多次评论它不会扩展到大型问题，以及为什么不会。所有这些问题都将在第18步，特别是第40步中得到解决，它可以扩展到非常多的进程。

从哲学上讲，MPI的运作方式如下。你通常通过以下方式运行一个程序

@code
  mpirun -np 32 ./step-17
@endcode

这意味着在（比如）32个处理器上运行它。如果你是在一个集群系统上，你通常需要<i>schedule</i>程序在32个处理器可用时运行；这将在你的集群的文档中描述。但是在系统内部，每当这些处理器可用时，通常会执行上述相同的调用）。)这样做的目的是，MPI系统将启动32个<i>copies</i>的 <code>step-17</code> 的可执行文件。(这些正在运行的可执行文件中的每一个的MPI术语是，你有32个 @ref GlossMPIProcess "MPI进程"。)这可能发生在不同的机器上，甚至不能从对方的内存空间中读取，也可能发生在同一台机器上，但最终的结果是一样的：这32个副本中的每一个都将以操作系统分配给它的一些内存运行，而且它不能直接读取其他31个副本的内存。为了在一个共同的任务中进行协作，这32个副本就必须<i>communicate</i>相互协作。MPI是<i>Message Passing Interface</i>的缩写，通过允许程序<i>send messages</i>来实现这一点。你可以把它看作是邮件服务：你可以把一封写给特定地址的信放入邮件，它将被送达。但这是你能控制事物的程度。如果你想让收信人对信的内容做些什么，例如把你想要的数据从那边返回给你，那么需要发生两件事。(i)接收方需要实际去检查他们的邮箱里是否有东西，(ii)如果有的话，做出适当的反应，比如说发送数据回来。如果你等待这个返回信息，但原来的接收者却心不在焉，没有注意，那么你就不走运了：你只需要等待，直到你在那边的请求将被解决。在某些情况下，错误会导致原始接收者永远不检查你的邮件，在这种情况下，你将永远等待--这被称为<i>deadlock</i>。(  @dealiiVideoLectureSeeAlso{39,41,41.25,41.5}) 

在实践中，人们通常不在发送和接收单个消息的层面上编程，而是使用更高层次的操作。例如，在程序中，我们将使用函数调用，从每个处理器获取一个数字，将它们全部相加，然后将总和返回给所有处理器。在内部，这是用单个消息实现的，但对用户来说，这是透明的。我们称这种操作为<i>collectives</i>，因为<i>all</i>处理器参与其中。集合体允许我们编写程序，其中不是每个可执行文件的副本都在做完全不同的事情（这将是难以置信的编程难度），但实质上所有副本都在为自己做同样的事情（尽管是在不同的数据上），通过相同的代码块运行；然后他们通过集合体进行数据通信；然后再回到为自己做事情，通过相同的数据块运行。这是能够编写程序的关键部分，也是确保程序能够在任何数量的处理器上运行的关键部分，因为我们不需要为每个参与的处理器编写不同的代码。

这并不是说程序从来都是以不同的处理器在其可执行文件的副本中运行不同的代码块的方式来编写的。程序内部也经常以其他方式而不是通过集合体进行通信。但是在实践中，%并行有限元代码几乎总是遵循这样的方案：程序的每个副本在同一时间运行相同的代码块，中间穿插着所有处理器相互交流的阶段）。)

在现实中，即使是调用MPI集体函数的水平也太低了。相反，下面的程序根本不会包含对MPI的任何直接调用，而只包含对deal.II的用户隐藏这种通信的函数。这样做的好处是，你不需要学习MPI的细节和相当复杂的函数调用。也就是说，你确实必须理解上文所述的MPI背后的一般哲学。




<a name="Whatthisprogramdoes"></a><h3>What this program does</h3>


然后，这个程序演示的技术是。

- 如何使用PETSc封装类；这在本程序的主类的声明中已经可以看到，  <code>ElasticProblem</code>  。

- 如何将网格划分为子域；这发生在 <code>ElasticProblem::setup_system()</code> 函数。

- 如何对运行在MPI网络上的作业进行并行化操作；在这里，这是在很多地方都要注意的，最明显的是在 <code>ElasticProblem::assemble_system()</code> 函数中。

- 如何处理只存储向量项子集的向量，对于这些向量，我们必须确保它们在当前处理器上存储我们需要的东西。例如见 <code>ElasticProblem::solve()</code> and <code>ElasticProblem::refine_grid()</code> 函数。

- 如何处理同时在多个处理器上运行的程序的状态输出。这是通过程序中的 <code>pcout</code> 变量完成的，在构造函数中初始化。

由于这一切只能用实际的代码来证明，让我们直接进入代码，不再多说。


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
 * First the usual assortment of header files we have already used in previous
 * example programs:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/multithread_info.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/sparsity_tools.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * And here come the things that we need particularly for this example program
 * and that weren't in step-8. First, we replace the standard output
 * <code>std::cout</code> by a new stream <code>pcout</code> which is used in
 * parallel computations for generating output only on one of the MPI
 * processes.
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * @endcode
 * 
 * We are going to query the number of processes and the number of the present
 * process by calling the respective functions in the Utilities::MPI
 * namespace.
 * 
 * @code
 * #include <deal.II/base/mpi.h>
 * @endcode
 * 
 * Then, we are going to replace all linear algebra components that involve
 * the (global) linear system by classes that wrap interfaces similar to our
 * own linear algebra classes around what PETSc offers (PETSc is a library
 * written in C, and deal.II comes with wrapper classes that provide the PETSc
 * functionality with an interface that is similar to the interface we already
 * had for our own linear algebra classes). In particular, we need vectors and
 * matrices that are distributed across several
 * @ref GlossMPIProcess "processes" in MPI programs (and
 * simply map to sequential, local vectors and matrices if there is only a
 * single process, i.e., if you are running on only one machine, and without
 * MPI support):
 * 
 * @code
 * #include <deal.II/lac/petsc_vector.h>
 * #include <deal.II/lac/petsc_sparse_matrix.h>
 * @endcode
 * 
 * Then we also need interfaces for solvers and preconditioners that PETSc
 * provides:
 * 
 * @code
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/petsc_precondition.h>
 * @endcode
 * 
 * And in addition, we need some algorithms for partitioning our meshes so
 * that they can be efficiently distributed across an MPI network. The
 * partitioning algorithm is implemented in the <code>GridTools</code>
 * namespace, and we need an additional include file for a function in
 * <code>DoFRenumbering</code> that allows to sort the indices associated with
 * degrees of freedom so that they are numbered according to the subdomain
 * they are associated with:
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * 
 * @endcode
 * 
 * And this is simply C++ again:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step17
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclasstemplate"></a> 
 * <h3>The <code>ElasticProblem</code> class template</h3>
 * 

 * 
 * The first real part of the program is the declaration of the main
 * class.  As mentioned in the introduction, almost all of this has
 * been copied verbatim from step-8, so we only comment on the few
 * differences between the two tutorials.  There is one (cosmetic)
 * change in that we let <code>solve</code> return a value, namely
 * the number of iterations it took to converge, so that we can
 * output this to the screen at the appropriate place.
 * 
 * @code
 *   template <int dim>
 *   class ElasticProblem
 *   {
 *   public:
 *     ElasticProblem();
 *     void run();
 * 
 *   private:
 *     void         setup_system();
 *     void         assemble_system();
 *     unsigned int solve();
 *     void         refine_grid();
 *     void         output_results(const unsigned int cycle) const;
 * 
 * @endcode
 * 
 * The first change is that we have to declare a variable that
 * indicates the @ref GlossMPICommunicator "MPI communicator" over
 * which we are supposed to distribute our computations.
 * 
 * @code
 *     MPI_Comm mpi_communicator;
 * 
 * @endcode
 * 
 * Then we have two variables that tell us where in the parallel
 * world we are. The first of the following variables,
 * <code>n_mpi_processes</code>, tells us how many MPI processes
 * there exist in total, while the second one,
 * <code>this_mpi_process</code>, indicates which is the number of
 * the present process within this space of processes (in MPI
 * language, this corresponds to the @ref GlossMPIRank "rank" of
 * the process). The latter will have a unique value for each
 * process between zero and (less than)
 * <code>n_mpi_processes</code>. If this program is run on a
 * single machine without MPI support, then their values are
 * <code>1</code> and <code>0</code>, respectively.
 * 
 * @code
 *     const unsigned int n_mpi_processes;
 *     const unsigned int this_mpi_process;
 * 
 * @endcode
 * 
 * Next up is a stream-like variable <code>pcout</code>. It is, in essence,
 * just something we use for convenience: in a parallel program,
 * if each process outputs status information, then there quickly
 * is a lot of clutter. Rather, we would want to only have one
 * @ref GlossMPIProcess "process" output everything once, for
 * example the one with @ref GlossMPIRank "rank" zero. At the same
 * time, it seems silly to prefix <i>every</i> place where we
 * create output with an <code>if (my_rank==0)</code> condition.
 *     

 * 
 * To make this simpler, the ConditionalOStream class does exactly
 * this under the hood: it acts as if it were a stream, but only
 * forwards to a real, underlying stream if a flag is set. By
 * setting this condition to <code>this_mpi_process==0</code>
 * (where <code>this_mpi_process</code> corresponds to the rank of
 * an MPI process), we make sure that output is only generated
 * from the first process and that we don't get the same lines of
 * output over and over again, once per process. Thus, we can use
 * <code>pcout</code> everywhere and in every process, but on all
 * but one process nothing will ever happen to the information
 * that is piped into the object via
 * <code>operator&lt;&lt;</code>.
 * 
 * @code
 *     ConditionalOStream pcout;
 * 
 * @endcode
 * 
 * The remainder of the list of member variables is fundamentally the
 * same as in step-8. However, we change the declarations of matrix
 * and vector types to use parallel PETSc objects instead. Note that
 * we do not use a separate sparsity pattern, since PETSc manages this
 * internally as part of its matrix data structures.
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 *     PETScWrappers::MPI::SparseMatrix system_matrix;
 * 
 *     PETScWrappers::MPI::Vector solution;
 *     PETScWrappers::MPI::Vector system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsidevalues"></a> 
 * <h3>Right hand side values</h3>
 * 

 * 
 * The following is taken from step-8 without change:
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override
 *     {
 *       Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
 *       Assert(dim >= 2, ExcInternalError());
 * 
 *       Point<dim> point_1, point_2;
 *       point_1(0) = 0.5;
 *       point_2(0) = -0.5;
 * 
 *       if (((p - point_1).norm_square() < 0.2 * 0.2) ||
 *           ((p - point_2).norm_square() < 0.2 * 0.2))
 *         values(0) = 1;
 *       else
 *         values(0) = 0;
 * 
 *       if (p.square() < 0.2 * 0.2)
 *         values(1) = 1;
 *       else
 *         values(1) = 0;
 *     }
 * 
 *     virtual void
 *     vector_value_list(const std::vector<Point<dim>> &points,
 *                       std::vector<Vector<double>> &  value_list) const override
 *     {
 *       const unsigned int n_points = points.size();
 * 
 *       Assert(value_list.size() == n_points,
 *              ExcDimensionMismatch(value_list.size(), n_points));
 * 
 *       for (unsigned int p = 0; p < n_points; ++p)
 *         RightHandSide<dim>::vector_value(points[p], value_list[p]);
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclassimplementation"></a> 
 * <h3>The <code>ElasticProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="ElasticProblemElasticProblem"></a> 
 * <h4>ElasticProblem::ElasticProblem</h4>
 * 

 * 
 * The first step in the actual implementation is the constructor of
 * the main class. Apart from initializing the same member variables
 * that we already had in step-8, we here initialize the MPI
 * communicator variable we shall use with the global MPI
 * communicator linking all processes together (in more complex
 * applications, one could here use a communicator object that only
 * links a subset of all processes), and call the Utilities::MPI
 * helper functions to determine the number of processes and where
 * the present one fits into this picture. In addition, we make sure
 * that output is only generated by the (globally) first process.
 * We do so by passing the stream we want to output to
 * (<code>std::cout</code>) and a true/false flag as arguments where
 * the latter is determined by testing whether the process currently
 * executing the constructor call is the first in the MPI universe.
 * 
 * @code
 *   template <int dim>
 *   ElasticProblem<dim>::ElasticProblem()
 *     : mpi_communicator(MPI_COMM_WORLD)
 *     , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
 *     , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
 *     , pcout(std::cout, (this_mpi_process == 0))
 *     , fe(FE_Q<dim>(1), dim)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsetup_system"></a> 
 * <h4>ElasticProblem::setup_system</h4>
 * 

 * 
 * Next, the function in which we set up the various variables
 * for the global linear system to be solved needs to be implemented.
 *   

 * 
 * However, before we proceed with this, there is one thing to do for a
 * parallel program: we need to determine which MPI process is
 * responsible for each of the cells. Splitting cells among
 * processes, commonly called "partitioning the mesh", is done by
 * assigning a @ref GlossSubdomainId "subdomain id" to each cell. We
 * do so by calling into the METIS library that does this in a very
 * efficient way, trying to minimize the number of nodes on the
 * interfaces between subdomains. Rather than trying to call METIS
 * directly, we do this by calling the
 * GridTools::partition_triangulation() function that does this at a
 * much higher level of programming.
 *   

 * 
 * @note As mentioned in the introduction, we could avoid this manual
 * partitioning step if we used the parallel::shared::Triangulation
 * class for the triangulation object instead (as we do in step-18).
 * That class does, in essence, everything a regular triangulation
 * does, but it then also automatically partitions the mesh after
 * every mesh creation or refinement operation.
 *   

 * 
 * Following partitioning, we need to enumerate all degrees of
 * freedom as usual.  However, we would like to enumerate the
 * degrees of freedom in a way so that all degrees of freedom
 * associated with cells in subdomain zero (which resides on process
 * zero) come before all DoFs associated with cells on subdomain
 * one, before those on cells on process two, and so on. We need
 * this since we have to split the global vectors for right hand
 * side and solution, as well as the matrix into contiguous chunks
 * of rows that live on each of the processors, and we will want to
 * do this in a way that requires minimal communication. This
 * particular enumeration can be obtained by re-ordering degrees of
 * freedom indices using DoFRenumbering::subdomain_wise().
 *   

 * 
 * The final step of this initial setup is that we get ourselves an
 * IndexSet that indicates the subset of the global number of unknowns
 * this
 * process is responsible for. (Note that a degree of freedom is not
 * necessarily owned by the process that owns a cell just because
 * the degree of freedom lives on this cell: some degrees of freedom
 * live on interfaces between subdomains, and are consequently only owned by
 * one of the processes adjacent to this interface.)
 *   

 * 
 * Before we move on, let us recall a fact already discussed in the
 * introduction: The triangulation we use here is replicated across
 * all processes, and each process has a complete copy of the entire
 * triangulation, with all cells. Partitioning only provides a way
 * to identify which cells out of all each process "owns", but it
 * knows everything about all of them. Likewise, the DoFHandler
 * object knows everything about every cell, in particular the
 * degrees of freedom that live on each cell, whether it is one that
 * the current process owns or not. This can not scale to large
 * problems because eventually just storing the entire mesh, and
 * everything that is associated with it, on every process will
 * become infeasible if the problem is large enough. On the other
 * hand, if we split the triangulation into parts so that every
 * process stores only those cells it "owns" but nothing else (or,
 * at least a sufficiently small fraction of everything else), then
 * we can solve large problems if only we throw a large enough
 * number of MPI processes at them. This is what we are going to in
 * step-40, for example, using the
 * parallel::distributed::Triangulation class.  On the other hand,
 * most of the rest of what we demonstrate in the current program
 * will actually continue to work whether we have the entire
 * triangulation available, or only a piece of it.
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::setup_system()
 *   {
 *     GridTools::partition_triangulation(n_mpi_processes, triangulation);
 * 
 *     dof_handler.distribute_dofs(fe);
 *     DoFRenumbering::subdomain_wise(dof_handler);
 * 
 * @endcode
 * 
 * We need to initialize the objects denoting hanging node constraints for
 * the present grid. As with the triangulation and DoFHandler objects, we
 * will simply store <i>all</i> constraints on each process; again, this
 * will not scale, but we show in step-40 how one can work around this by
 * only storing on each MPI process the constraints for degrees of freedom
 * that actually matter on this particular process.
 * 
 * @code
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 * @endcode
 * 
 * Now we create the sparsity pattern for the system matrix. Note that we
 * again compute and store all entries and not only the ones relevant
 * to this process (see step-18 or step-40 for a more efficient way to
 * handle this).
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     hanging_node_constraints,
 *                                     false);
 * 
 * @endcode
 * 
 * Now we determine the set of locally owned DoFs and use that to
 * initialize parallel vectors and matrix. Since the matrix and vectors
 * need to work in parallel, we have to pass them an MPI communication
 * object, as well as information about the partitioning contained in the
 * IndexSet @p locally_owned_dofs.  The IndexSet contains information
 * about the global size (the <i>total</i> number of degrees of freedom)
 * and also what subset of rows is to be stored locally.  Note that the
 * system matrix needs that partitioning information for the rows and
 * columns. For square matrices, as it is the case here, the columns
 * should be partitioned in the same way as the rows, but in the case of
 * rectangular matrices one has to partition the columns in the same way
 * as vectors are partitioned with which the matrix is multiplied, while
 * rows have to partitioned in the same way as destination vectors of
 * matrix-vector multiplications:
 * 
 * @code
 *     const std::vector<IndexSet> locally_owned_dofs_per_proc =
 *       DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
 *     const IndexSet locally_owned_dofs =
 *       locally_owned_dofs_per_proc[this_mpi_process];
 * 
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          dsp,
 *                          mpi_communicator);
 * 
 *     solution.reinit(locally_owned_dofs, mpi_communicator);
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemassemble_system"></a> 
 * <h4>ElasticProblem::assemble_system</h4>
 * 

 * 
 * We now assemble the matrix and right hand side of the
 * problem. There are some things worth mentioning before we go into
 * detail. First, we will be assembling the system in parallel,
 * i.e., each process will be responsible for assembling on cells
 * that belong to this particular process. Note that the degrees of
 * freedom are split in a way such that all DoFs in the interior of
 * cells and between cells belonging to the same subdomain belong to
 * the process that <code>owns</code> the cell. However, even then
 * we sometimes need to assemble on a cell with a neighbor that
 * belongs to a different process, and in these cases when we add up
 * the local contributions into the global matrix or right hand side
 * vector, we have to transfer these entries to the process that
 * owns these elements. Fortunately, we don't have to do this by
 * hand: PETSc does all this for us by caching these elements
 * locally, and sending them to the other processes as necessary
 * when we call the <code>compress()</code> functions on the matrix
 * and vector at the end of this function.
 *   

 * 
 * The second point is that once we have handed over matrix and
 * vector contributions to PETSc, it is a) hard, and b) very
 * inefficient to get them back for modifications. This is not only
 * the fault of PETSc, it is also a consequence of the distributed
 * nature of this program: if an entry resides on another processor,
 * then it is necessarily expensive to get it. The consequence of
 * this is that we should not try to first assemble the matrix and
 * right hand side as if there were no hanging node constraints and
 * boundary values, and then eliminate these in a second step
 * (using, for example, AffineConstraints::condense()). Rather, we
 * should try to eliminate hanging node constraints before handing
 * these entries over to PETSc. This is easy: instead of copying
 * elements by hand into the global matrix (as we do in step-4), we
 * use the AffineConstraints::distribute_local_to_global() functions
 * to take care of hanging nodes at the same time. We also already
 * did this in step-6. The second step, elimination of boundary
 * nodes, could also be done this way by putting the boundary values
 * into the same AffineConstraints object as hanging nodes (see the
 * way it is done in step-6, for example); however, it is not
 * strictly necessary to do this here because eliminating boundary
 * values can be done with only the data stored on each process
 * itself, and consequently we use the approach used before in
 * step-4, i.e., via MatrixTools::apply_boundary_values().
 *   

 * 
 * All of this said, here is the actual implementation starting with
 * the general setup of helper variables.  (Note that we still use
 * the deal.II full matrix and vector types for the local systems as
 * these are small and need not be shared across processes.)
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim>   quadrature_formula(fe.degree + 1);
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
 *     std::vector<double> lambda_values(n_q_points);
 *     std::vector<double> mu_values(n_q_points);
 * 
 *     Functions::ConstantFunction<dim> lambda(1.), mu(1.);
 * 
 *     RightHandSide<dim>          right_hand_side;
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));
 * 
 * 
 * @endcode
 * 
 * The next thing is the loop over all elements. Note that we do
 * not have to do <i>all</i> the work on every process: our job
 * here is only to assemble the system on cells that actually
 * belong to this MPI process, all other cells will be taken care
 * of by other processes. This is what the if-clause immediately
 * after the for-loop takes care of: it queries the subdomain
 * identifier of each cell, which is a number associated with each
 * cell that tells us about the owner process. In more generality,
 * the subdomain id is used to split a domain into several parts
 * (we do this above, at the beginning of
 * <code>setup_system()</code>), and which allows to identify
 * which subdomain a cell is living on. In this application, we
 * have each process handle exactly one subdomain, so we identify
 * the terms <code>subdomain</code> and <code>MPI process</code>.
 *     

 * 
 * Apart from this, assembling the local system is relatively uneventful
 * if you have understood how this is done in step-8. As mentioned above,
 * distributing local contributions into the global matrix
 * and right hand sides also takes care of hanging node constraints in the
 * same way as is done in step-6.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->subdomain_id() == this_mpi_process)
 *         {
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 *           fe_values.reinit(cell);
 * 
 *           lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
 *           mu.value_list(fe_values.get_quadrature_points(), mu_values);
 * 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const unsigned int component_i =
 *                 fe.system_to_component_index(i).first;
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const unsigned int component_j =
 *                     fe.system_to_component_index(j).first;
 * 
 *                   for (unsigned int q_point = 0; q_point < n_q_points;
 *                        ++q_point)
 *                     {
 *                       cell_matrix(i, j) +=
 *                         ((fe_values.shape_grad(i, q_point)[component_i] *
 *                           fe_values.shape_grad(j, q_point)[component_j] *
 *                           lambda_values[q_point]) +
 *                          (fe_values.shape_grad(i, q_point)[component_j] *
 *                           fe_values.shape_grad(j, q_point)[component_i] *
 *                           mu_values[q_point]) +
 *                          ((component_i == component_j) ?
 *                             (fe_values.shape_grad(i, q_point) *
 *                              fe_values.shape_grad(j, q_point) *
 *                              mu_values[q_point]) :
 *                             0)) *
 *                         fe_values.JxW(q_point);
 *                     }
 *                 }
 *             }
 * 
 *           right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
 *                                             rhs_values);
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const unsigned int component_i =
 *                 fe.system_to_component_index(i).first;
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *                 cell_rhs(i) += fe_values.shape_value(i, q_point) *
 *                                rhs_values[q_point](component_i) *
 *                                fe_values.JxW(q_point);
 *             }
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix,
 *                                                               cell_rhs,
 *                                                               local_dof_indices,
 *                                                               system_matrix,
 *                                                               system_rhs);
 *         }
 * 
 * @endcode
 * 
 * The next step is to "compress" the vector and the system matrix. This
 * means that each process sends the additions that were made to those
 * entries of the matrix and vector that the process did not own itself to
 * the process that owns them. After receiving these additions from other
 * processes, each process then adds them to the values it already
 * has. These additions are combining the integral contributions of shape
 * functions living on several cells just as in a serial computation, with
 * the difference that the cells are assigned to different processes.
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 * 
 * @endcode
 * 
 * The global matrix and right hand side vectors have now been
 * formed. We still have to apply boundary values, in the same way as we
 * did, for example, in step-3, step-4, and a number of other programs.
 *     

 * 
 * The last argument to the call to
 * MatrixTools::apply_boundary_values() below allows for some
 * optimizations. It controls whether we should also delete
 * entries (i.e., set them to zero) in the matrix columns
 * corresponding to boundary nodes, or to keep them (and passing
 * <code>true</code> means: yes, do eliminate the columns). If we
 * do eliminate columns, then the resulting matrix will be
 * symmetric again if it was before; if we don't, then it
 * won't. The solution of the resulting system should be the same,
 * though. The only reason why we may want to make the system
 * symmetric again is that we would like to use the CG method,
 * which only works with symmetric matrices. The reason why we may
 * <i>not</i> want to make the matrix symmetric is because this
 * would require us to write into column entries that actually
 * reside on other processes, i.e., it involves communicating
 * data. This is always expensive.
 *     

 * 
 * Experience tells us that CG also works (and works almost as
 * well) if we don't remove the columns associated with boundary
 * nodes, which can be explained by the special structure of this
 * particular non-symmetry. To avoid the expense of communication,
 * we therefore do not eliminate the entries in the affected
 * columns.
 * 
 * @code
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(dim),
 *                                              boundary_values);
 *     MatrixTools::apply_boundary_values(
 *       boundary_values, system_matrix, solution, system_rhs, false);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsolve"></a> 
 * <h4>ElasticProblem::solve</h4>
 * 

 * 
 * Having assembled the linear system, we next need to solve
 * it. PETSc offers a variety of sequential and parallel solvers,
 * for which we have written wrappers that have almost the same
 * interface as is used for the deal.II solvers used in all previous
 * example programs. The following code should therefore look rather
 * familiar.
 *   

 * 
 * At the top of the function, we set up a convergence monitor, and
 * assign it the accuracy to which we would like to solve the linear
 * system. Next, we create an actual solver object using PETSc's CG
 * solver which also works with parallel (distributed) vectors and
 * matrices. And finally a preconditioner; we choose to use a block
 * Jacobi preconditioner which works by computing an incomplete LU
 * decomposition on each diagonal block of the matrix.  (In other
 * words, each MPI process computes an ILU from the rows it stores
 * by throwing away columns that correspond to row indices not
 * stored locally; this yields a square matrix block from which we
 * can compute an ILU. That means that if you run the program with
 * only one process, then you will use an ILU(0) as a
 * preconditioner, while if it is run on many processes, then we
 * will have a number of blocks on the diagonal and the
 * preconditioner is the ILU(0) of each of these blocks. In the
 * extreme case of one degree of freedom per processor, this
 * preconditioner is then simply a Jacobi preconditioner since the
 * diagonal matrix blocks consist of only a single entry. Such a
 * preconditioner is relatively easy to compute because it does not
 * require any kind of communication between processors, but it is
 * in general not very efficient for large numbers of processors.)
 *   

 * 
 * Following this kind of setup, we then solve the linear system:
 * 
 * @code
 *   template <int dim>
 *   unsigned int ElasticProblem<dim>::solve()
 *   {
 *     SolverControl solver_control(solution.size(), 1e-8 * system_rhs.l2_norm());
 *     PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
 * 
 *     PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 * @endcode
 * 
 * The next step is to distribute hanging node constraints. This is a
 * little tricky, since to fill in the value of a constrained node you
 * need access to the values of the nodes to which it is constrained (for
 * example, for a Q1 element in 2d, we need access to the two nodes on the
 * big side of a hanging node face, to compute the value of the
 * constrained node in the middle).
 *     

 * 
 * The problem is that we have built our vectors (in
 * <code>setup_system()</code>) in such a way that every process
 * is responsible for storing only those elements of the solution
 * vector that correspond to the degrees of freedom this process
 * "owns". There are, however, cases where in order to compute the
 * value of the vector entry for a constrained degree of freedom
 * on one process, we need to access vector entries that are
 * stored on other processes.  PETSc (and, for that matter, the
 * MPI model on which it is built) does not allow to simply query
 * vector entries stored on other processes, so what we do here is
 * to get a copy of the "distributed" vector where we store all
 * elements locally. This is simple, since the deal.II wrappers
 * have a conversion constructor for the deal.II Vector
 * class. (This conversion of course requires communication, but
 * in essence every process only needs to send its data to every
 * other process once in bulk, rather than having to respond to
 * queries for individual elements):
 * 
 * @code
 *     Vector<double> localized_solution(solution);
 * 
 * @endcode
 * 
 * Of course, as in previous discussions, it is clear that such a
 * step cannot scale very far if you wanted to solve large
 * problems on large numbers of processes, because every process
 * now stores <i>all elements</i> of the solution vector. (We will
 * show how to do this better in step-40.)  On the other hand,
 * distributing hanging node constraints is simple on this local
 * copy, using the usual function
 * AffineConstraints::distributed(). In particular, we can compute
 * the values of <i>all</i> constrained degrees of freedom,
 * whether the current process owns them or not:
 * 
 * @code
 *     hanging_node_constraints.distribute(localized_solution);
 * 
 * @endcode
 * 
 * Then transfer everything back into the global vector. The
 * following operation copies those elements of the localized
 * solution that we store locally in the distributed solution, and
 * does not touch the other ones. Since we do the same operation
 * on all processors, we end up with a distributed vector (i.e., a
 * vector that on every process only stores the vector entries
 * corresponding to degrees of freedom that are owned by this
 * process) that has all the constrained nodes fixed.
 *     

 * 
 * We end the function by returning the number of iterations it
 * took to converge, to allow for some output.
 * 
 * @code
 *     solution = localized_solution;
 * 
 *     return solver_control.last_step();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrefine_grid"></a> 
 * <h4>ElasticProblem::refine_grid</h4>
 * 

 * 
 * Using some kind of refinement indicator, the mesh can be
 * refined. The problem is basically the same as with distributing
 * hanging node constraints: in order to compute the error indicator
 * (even if we were just interested in the indicator on the cells
 * the current process owns), we need access to more elements of the
 * solution vector than just those the current processor stores. To
 * make this happen, we do essentially what we did in
 * <code>solve()</code> already, namely get a <i>complete</i> copy
 * of the solution vector onto every process, and use that to
 * compute. This is in itself expensive as explained above and it
 * is particular unnecessary since we had just created and then
 * destroyed such a vector in <code>solve()</code>, but efficiency
 * is not the point of this program and so let us opt for a design
 * in which every function is as self-contained as possible.
 *   

 * 
 * Once we have such a "localized" vector that contains <i>all</i>
 * elements of the solution vector, we can compute the indicators
 * for the cells that belong to the present process. In fact, we
 * could of course compute <i>all</i> refinement indicators since
 * our Triangulation and DoFHandler objects store information about
 * all cells, and since we have a complete copy of the solution
 * vector. But in the interest in showing how to operate in
 * %parallel, let us demonstrate how one would operate if one were
 * to only compute <i>some</i> error indicators and then exchange
 * the remaining ones with the other processes. (Ultimately, each
 * process needs a complete set of refinement indicators because
 * every process needs to refine their mesh, and needs to refine it
 * in exactly the same way as all of the other processes.)
 *   

 * 
 * So, to do all of this, we need to:
 * - First, get a local copy of the distributed solution vector.
 * - Second, create a vector to store the refinement indicators.
 * - Third, let the KellyErrorEstimator compute refinement
 * indicators for all cells belonging to the present
 * subdomain/process. The last argument of the call indicates
 * which subdomain we are interested in. The three arguments
 * before it are various other default arguments that one usually
 * does not need (and does not state values for, but rather uses the
 * defaults), but which we have to state here explicitly since we
 * want to modify the value of a following argument (i.e., the one
 * indicating the subdomain).
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::refine_grid()
 *   {
 *     const Vector<double> localized_solution(solution);
 * 
 *     Vector<float> local_error_per_cell(triangulation.n_active_cells());
 *     KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                        QGauss<dim - 1>(fe.degree + 1),
 *                                        {},
 *                                        localized_solution,
 *                                        local_error_per_cell,
 *                                        ComponentMask(),
 *                                        nullptr,
 *                                        MultithreadInfo::n_threads(),
 *                                        this_mpi_process);
 * 
 * @endcode
 * 
 * Now all processes have computed error indicators for their own
 * cells and stored them in the respective elements of the
 * <code>local_error_per_cell</code> vector. The elements of this
 * vector for cells not owned by the present process are
 * zero. However, since all processes have a copy of the entire
 * triangulation and need to keep these copies in sync, they need
 * the values of refinement indicators for all cells of the
 * triangulation. Thus, we need to distribute our results. We do
 * this by creating a distributed vector where each process has
 * its share and sets the elements it has computed. Consequently,
 * when you view this vector as one that lives across all
 * processes, then every element of this vector has been set
 * once. We can then assign this parallel vector to a local,
 * non-parallel vector on each process, making <i>all</i> error
 * indicators available on every process.
 *     

 * 
 * So in the first step, we need to set up a parallel vector. For
 * simplicity, every process will own a chunk with as many
 * elements as this process owns cells, so that the first chunk of
 * elements is stored with process zero, the next chunk with
 * process one, and so on. It is important to remark, however,
 * that these elements are not necessarily the ones we will write
 * to. This is a consequence of the order in which cells are arranged,
 * i.e., the order in which the elements of the vector correspond
 * to cells is not ordered according to the subdomain these cells
 * belong to. In other words, if on this process we compute
 * indicators for cells of a certain subdomain, we may write the
 * results to more or less random elements of the distributed
 * vector; in particular, they may not necessarily lie within the
 * chunk of vector we own on the present process. They will
 * subsequently have to be copied into another process' memory
 * space, an operation that PETSc does for us when we call the
 * <code>compress()</code> function. This inefficiency could be
 * avoided with some more code, but we refrain from it since it is
 * not a major factor in the program's total runtime.
 *     

 * 
 * So here is how we do it: count how many cells belong to this
 * process, set up a distributed vector with that many elements to
 * be stored locally, copy over the elements we computed
 * locally, and finally compress the result. In fact, we really only copy
 * the elements that are nonzero, so we may miss a few that we
 * computed to zero, but this won't hurt since the original values
 * of the vector are zero anyway.
 * 
 * @code
 *     const unsigned int n_local_cells =
 *       GridTools::count_cells_with_subdomain_association(triangulation,
 *                                                         this_mpi_process);
 *     PETScWrappers::MPI::Vector distributed_all_errors(
 *       mpi_communicator, triangulation.n_active_cells(), n_local_cells);
 * 
 *     for (unsigned int i = 0; i < local_error_per_cell.size(); ++i)
 *       if (local_error_per_cell(i) != 0)
 *         distributed_all_errors(i) = local_error_per_cell(i);
 *     distributed_all_errors.compress(VectorOperation::insert);
 * 
 * 
 * @endcode
 * 
 * So now we have this distributed vector that contains the
 * refinement indicators for all cells. To use it, we need to
 * obtain a local copy and then use it to mark cells for
 * refinement or coarsening, and actually do the refinement and
 * coarsening. It is important to recognize that <i>every</i>
 * process does this to its own copy of the triangulation, and
 * does it in exactly the same way.
 * 
 * @code
 *     const Vector<float> localized_all_errors(distributed_all_errors);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     localized_all_errors,
 *                                                     0.3,
 *                                                     0.03);
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemoutput_results"></a> 
 * <h4>ElasticProblem::output_results</h4>
 * 

 * 
 * The final function of significant interest is the one that
 * creates graphical output. This works the same way as in step-8,
 * with two small differences. Before discussing these, let us state
 * the general philosophy this function will work: we intend for all
 * of the data to be generated on a single process, and subsequently
 * written to a file. This is, as many other parts of this program
 * already discussed, not something that will scale. Previously, we
 * had argued that we will get into trouble with triangulations,
 * DoFHandlers, and copies of the solution vector where every
 * process has to store all of the data, and that there will come to
 * be a point where each process simply doesn't have enough memory
 * to store that much data. Here, the situation is different: it's
 * not only the memory, but also the run time that's a problem. If
 * one process is responsible for processing <i>all</i> of the data
 * while all of the other processes do nothing, then this one
 * function will eventually come to dominate the overall run time of
 * the program.  In particular, the time this function takes is
 * going to be proportional to the overall size of the problem
 * (counted in the number of cells, or the number of degrees of
 * freedom), independent of the number of processes we throw at it.
 *   

 * 
 * Such situations need to be avoided, and we will show in step-18
 * and step-40 how to address this issue. For the current problem,
 * the solution is to have each process generate output data only
 * for its own local cells, and write them to separate files, one
 * file per process. This is how step-18 operates. Alternatively,
 * one could simply leave everything in a set of independent files
 * and let the visualization software read all of them (possibly
 * also using multiple processors) and create a single visualization
 * out of all of them; this is the path step-40, step-32, and all
 * other parallel programs developed later on take.
 *   

 * 
 * More specifically for the current function, all processes call
 * this function, but not all of them need to do the work associated
 * with generating output. In fact, they shouldn't, since we would
 * try to write to the same file multiple times at once. So we let
 * only the first process do this, and all the other ones idle
 * around during this time (or start their work for the next
 * iteration, or simply yield their CPUs to other jobs that happen
 * to run at the same time). The second thing is that we not only
 * output the solution vector, but also a vector that indicates
 * which subdomain each cell belongs to. This will make for some
 * nice pictures of partitioned domains.
 *   

 * 
 * To implement this, process zero needs a complete set of solution
 * components in a local vector. Just as with the previous function,
 * the efficient way to do this would be to re-use the vector
 * already created in the <code>solve()</code> function, but to keep
 * things more self-contained, we simply re-create one here from the
 * distributed solution vector.
 *   

 * 
 * An important thing to realize is that we do this localization operation
 * on all processes, not only the one that actually needs the data. This
 * can't be avoided, however, with the simplified communication model of MPI
 * we use for vectors in this tutorial program: MPI does not have a way to
 * query data on another process, both sides have to initiate a
 * communication at the same time. So even though most of the processes do
 * not need the localized solution, we have to place the statement
 * converting the distributed into a localized vector so that all processes
 * execute it.
 *   

 * 
 * (Part of this work could in fact be avoided. What we do is
 * send the local parts of all processes to all other processes. What we
 * would really need to do is to initiate an operation on all processes
 * where each process simply sends its local chunk of data to process
 * zero, since this is the only one that actually needs it, i.e., we need
 * something like a gather operation. PETSc can do this, but for
 * simplicity's sake we don't attempt to make use of this here. We don't,
 * since what we do is not very expensive in the grand scheme of things:
 * it is one vector communication among all processes, which has to be
 * compared to the number of communications we have to do when solving the
 * linear system, setting up the block-ILU for the preconditioner, and
 * other operations.)
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     const Vector<double> localized_solution(solution);
 * 
 * @endcode
 * 
 * This being done, process zero goes ahead with setting up the
 * output file as in step-8, and attaching the (localized)
 * solution vector to the output object.
 * 
 * @code
 *     if (this_mpi_process == 0)
 *       {
 *         std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
 * 
 *         DataOut<dim> data_out;
 *         data_out.attach_dof_handler(dof_handler);
 * 
 *         std::vector<std::string> solution_names;
 *         switch (dim)
 *           {
 *             case 1:
 *               solution_names.emplace_back("displacement");
 *               break;
 *             case 2:
 *               solution_names.emplace_back("x_displacement");
 *               solution_names.emplace_back("y_displacement");
 *               break;
 *             case 3:
 *               solution_names.emplace_back("x_displacement");
 *               solution_names.emplace_back("y_displacement");
 *               solution_names.emplace_back("z_displacement");
 *               break;
 *             default:
 *               Assert(false, ExcInternalError());
 *           }
 * 
 *         data_out.add_data_vector(localized_solution, solution_names);
 * 
 * @endcode
 * 
 * The only other thing we do here is that we also output one
 * value per cell indicating which subdomain (i.e., MPI
 * process) it belongs to. This requires some conversion work,
 * since the data the library provides us with is not the one
 * the output class expects, but this is not difficult. First,
 * set up a vector of integers, one per cell, that is then
 * filled by the subdomain id of each cell.
 *         

 * 
 * The elements of this vector are then converted to a
 * floating point vector in a second step, and this vector is
 * added to the DataOut object, which then goes off creating
 * output in VTK format:
 * 
 * @code
 *         std::vector<unsigned int> partition_int(triangulation.n_active_cells());
 *         GridTools::get_subdomain_association(triangulation, partition_int);
 * 
 *         const Vector<double> partitioning(partition_int.begin(),
 *                                           partition_int.end());
 * 
 *         data_out.add_data_vector(partitioning, "partitioning");
 * 
 *         data_out.build_patches();
 *         data_out.write_vtk(output);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrun"></a> 
 * <h4>ElasticProblem::run</h4>
 * 

 * 
 * Lastly, here is the driver function. It is almost completely
 * unchanged from step-8, with the exception that we replace
 * <code>std::cout</code> by the <code>pcout</code> stream. Apart
 * from this, the only other cosmetic change is that we output how
 * many degrees of freedom there are per process, and how many
 * iterations it took for the linear solver to converge:
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 10; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, -1, 1);
 *             triangulation.refine_global(3);
 *           }
 *         else
 *           refine_grid();
 * 
 *         pcout << "   Number of active cells:       "
 *               << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (by partition:";
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p)
 *           pcout << (p == 0 ? ' ' : '+')
 *                 << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
 *                                                                     p));
 *         pcout << ")" << std::endl;
 * 
 *         assemble_system();
 *         const unsigned int n_iterations = solve();
 * 
 *         pcout << "   Solver converged in " << n_iterations << " iterations."
 *               << std::endl;
 * 
 *         output_results(cycle);
 *       }
 *   }
 * } // namespace Step17
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The <code>main()</code> works the same way as most of the main
 * functions in the other example programs, i.e., it delegates work to
 * the <code>run</code> function of a managing object, and only wraps
 * everything into some code to catch exceptions:
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step17;
 * 
 * @endcode
 * 
 * Here is the only real difference: MPI and PETSc both require that we
 * initialize these libraries at the beginning of the program, and
 * un-initialize them at the end. The class MPI_InitFinalize takes care
 * of all of that. The trailing argument `1` means that we do want to
 * run each MPI process with a single thread, a prerequisite with the
 * PETSc parallel linear algebra.
 * 
 * @code
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       ElasticProblem<2> elastic_problem;
 *       elastic_problem.run();
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
examples/step-17/doc/results.dox



<a name="Results"></a><h1>Results</h1>



如果上述程序在单处理器机器上编译和运行，它产生的结果应该与我们已经通过步骤8得到的结果非常相似。然而，如果我们在多核机器或计算机集群上运行它，就会变得更加有趣。运行MPI程序的最基本方法是使用一个命令行，如

@code
  mpirun -np 32 ./step-17
@endcode

以32个处理器运行step-17可执行文件。

如果你在一个集群上工作，那么中间通常有一个步骤，你需要设置一个作业脚本，并将该脚本提交给调度器。只要调度器能够为你的工作分配32个未使用的处理器，它就会执行这个脚本。如何编写这样的作业脚本因集群而异，你应该找到你的集群的文档来看看如何做。在我的系统上，我必须使用 <code>qsub</code> 这个命令，加上一大堆的选项来并行运行一个作业）。)

无论是直接还是通过调度器，如果你在8个处理器上运行这个程序，你应该得到如下输出。

@code
Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 162 (by partition: 22+22+20+20+18+16+20+24)
   Solver converged in 23 iterations.
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 302 (by partition: 38+42+36+34+44+44+36+28)
   Solver converged in 35 iterations.
Cycle 2:
   Number of active cells:       238
   Number of degrees of freedom: 570 (by partition: 68+80+66+74+58+68+78+78)
   Solver converged in 46 iterations.
Cycle 3:
   Number of active cells:       454
   Number of degrees of freedom: 1046 (by partition: 120+134+124+130+154+138+122+124)
   Solver converged in 55 iterations.
Cycle 4:
   Number of active cells:       868
   Number of degrees of freedom: 1926 (by partition: 232+276+214+248+230+224+234+268)
   Solver converged in 77 iterations.
Cycle 5:
   Number of active cells:       1654
   Number of degrees of freedom: 3550 (by partition: 418+466+432+470+442+474+424+424)
   Solver converged in 93 iterations.
Cycle 6:
   Number of active cells:       3136
   Number of degrees of freedom: 6702 (by partition: 838+796+828+892+866+798+878+806)
   Solver converged in 127 iterations.
Cycle 7:
   Number of active cells:       5962
   Number of degrees of freedom: 12446 (by partition: 1586+1484+1652+1552+1556+1576+1560+1480)
   Solver converged in 158 iterations.
Cycle 8:
   Number of active cells:       11320
   Number of degrees of freedom: 23586 (by partition: 2988+2924+2890+2868+2864+3042+2932+3078)
   Solver converged in 225 iterations.
Cycle 9:
   Number of active cells:       21424
   Number of degrees of freedom: 43986 (by partition: 5470+5376+5642+5450+5630+5470+5416+5532)
   Solver converged in 282 iterations.
Cycle 10:
   Number of active cells:       40696
   Number of degrees of freedom: 83754 (by partition: 10660+10606+10364+10258+10354+10322+10586+10604)
   Solver converged in 392 iterations.
Cycle 11:
   Number of active cells:       76978
   Number of degrees of freedom: 156490 (by partition: 19516+20148+19390+19390+19336+19450+19730+19530)
   Solver converged in 509 iterations.
Cycle 12:
   Number of active cells:       146206
   Number of degrees of freedom: 297994 (by partition: 37462+37780+37000+37060+37232+37328+36860+37272)
   Solver converged in 705 iterations.
Cycle 13:
   Number of active cells:       276184
   Number of degrees of freedom: 558766 (by partition: 69206+69404+69882+71266+70348+69616+69796+69248)
   Solver converged in 945 iterations.
Cycle 14:
   Number of active cells:       523000
   Number of degrees of freedom: 1060258 (by partition: 132928+132296+131626+132172+132170+133588+132252+133226)
   Solver converged in 1282 iterations.
Cycle 15:
   Number of active cells:       987394
   Number of degrees of freedom: 1994226 (by partition: 253276+249068+247430+248402+248496+251380+248272+247902)
   Solver converged in 1760 iterations.
Cycle 16:
   Number of active cells:       1867477
   Number of degrees of freedom: 3771884 (by partition: 468452+474204+470818+470884+469960+
471186+470686+475694)
   Solver converged in 2251 iterations.
@endcode

(这次运行比examples/目录中的代码多用了几个细化周期。该运行还使用了2004年的METIS版本，产生了不同的分区；因此，你今天得到的数字略有不同）。)

可以看出，我们可以很容易地达到近400万个未知数。事实上，这段代码在8个进程中的运行时间不到7分钟，直到（包括）第14个周期，14分钟包括倒数第二步。(这些数字与该代码最初编写的时间有关，即2004年。)虽然我失去了最后一步的时间信息，但你会明白的。所有这些都是在通过运行 <code>make release</code> 启用发布模式之后，并且由于上述程序注释中所述的原因，关闭了图形输出的生成。(  @dealiiVideoLectureSeeAlso{18})  我做的最大的2D计算大约有710万个未知数，在32个进程上完成。花了大约40分钟。毫不奇怪，一个人能够走多远的限制因素是他有多少内存，因为每个进程都必须持有整个网格和DoFHandler对象，尽管矩阵和向量被分割开来。对于7.1M的计算，每个未知数的内存消耗约为600字节，这并不坏，但我们必须考虑到这是针对每个未知数的，无论我们是否在本地存储矩阵和向量条目。




下面是在程序的第12个周期中产生的一些输出，即大约有30万个未知数。

 <table align="center" style="width:80%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-ux.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-uy.png" alt="" width="100%"></td>
  </tr>
</table> 

正如人们所希望的那样，这里显示的X位移（左）和Y位移（右）与我们在第8步中已经看到的密切相关。正如第22步所示，我们也可以制作一个位移场的矢量图，而不是把它绘制成两个独立的标量场。不过，可能更有趣的是，在这一步看一下网格和分区。

 <table align="center" width="80%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-grid.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-partition.png" alt="" width="100%"></td>
  </tr>
</table> 

同样，网格（左边）显示了与之前看到的相同的细化模式。右图显示了8个过程中的域的划分，每个过程用不同的颜色表示。图片显示，在网格单元较小的地方，子域较小，考虑到分区算法试图平衡每个子域的单元数，这一事实是需要预期的；这种平衡在上图所示的输出中也很容易识别，每个子域的度数大致相同。




值得注意的是，如果我们用不同的进程数来运行同一个程序，我们可能会得到稍微不同的输出：不同的网格，不同的未知数和迭代收敛的次数。其原因是，虽然矩阵和右手边是相同的，与使用的进程数无关，但预处理程序不是：它对每个处理器的 <em>  矩阵块分别执行ILU(0)  </em>  。因此，随着进程数的增加，它作为预处理程序的有效性会降低，这使得迭代次数增加。由于不同的预处理程序会导致计算出的解有细微的变化，这将导致细化时标记的网格单元略有不同，在后续步骤中的差异也更大。不过，解决方案看起来总是非常相似的。




最后，这里是3D模拟的一些结果。你可以通过改变以下内容来重复这些结果

@code
        ElasticProblem<2> elastic_problem;
@endcode

至

@code
        ElasticProblem<3> elastic_problem;
@endcode

在主函数中。如果你再并行运行该程序，你会得到与此类似的东西（这是针对一个有16个进程的工作）。

@code
Cycle 0:
   Number of active cells:       512
   Number of degrees of freedom: 2187 (by partition: 114+156+150+114+114+210+105+102+120+120+96+123+141+183+156+183)
   Solver converged in 27 iterations.
Cycle 1:
   Number of active cells:       1604
   Number of degrees of freedom: 6549 (by partition: 393+291+342+354+414+417+570+366+444+288+543+525+345+387+489+381)
   Solver converged in 42 iterations.
Cycle 2:
   Number of active cells:       4992
   Number of degrees of freedom: 19167 (by partition: 1428+1266+1095+1005+1455+1257+1410+1041+1320+1380+1080+1050+963+1005+1188+1224)
   Solver converged in 65 iterations.
Cycle 3:
   Number of active cells:       15485
   Number of degrees of freedom: 56760 (by partition: 3099+3714+3384+3147+4332+3858+3615+3117+3027+3888+3942+3276+4149+3519+3030+3663)
   Solver converged in 96 iterations.
Cycle 4:
   Number of active cells:       48014
   Number of degrees of freedom: 168762 (by partition: 11043+10752+9846+10752+9918+10584+10545+11433+12393+11289+10488+9885+10056+9771+11031+8976)
   Solver converged in 132 iterations.
Cycle 5:
   Number of active cells:       148828
   Number of degrees of freedom: 492303 (by partition: 31359+30588+34638+32244+30984+28902+33297+31569+29778+29694+28482+28032+32283+30702+31491+28260)
   Solver converged in 179 iterations.
Cycle 6:
   Number of active cells:       461392
   Number of degrees of freedom: 1497951 (by partition: 103587+100827+97611+93726+93429+88074+95892+88296+96882+93000+87864+90915+92232+86931+98091+90594)
   Solver converged in 261 iterations.
@endcode






最后一步，达到150万个未知数，在8台双处理器机器（2003年可用的那种）上进行16个进程，需要大约55分钟。这个工作产生的图形输出相当大（第5周期已经打印了大约82MB的数据），所以我们要显示第4周期的输出。

 <table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.4-3d-partition.png" width="100%" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.4-3d-ux.png" alt="" width="100%"></td>
  </tr>
</table> 




左图显示的是将立方体划分为16个过程，而右图显示的是沿两个切面通过立方体的X位移。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


该程序在每个处理器上都保留一份三角形和DoFHandler对象的完整副本。它还创建了解决方案矢量的完整副本，并且只在一个处理器上创建输出。就并行化而言，所有这些显然是瓶颈。

在内部，在deal.II中，将分层和非结构化三角计算中使用的数据结构并行化是一个难点，我们又花了几年时间才实现了这一点。step-40教程程序和 @ref distributed 文档模块谈到了如何做这些步骤，以及从应用的角度来看需要什么。当前程序的一个明显的扩展是使用这个功能将计算完全分布到比这里使用的更多的处理器。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-17.cc"
*/
