/**
@page step_36 The step-36 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#EigenvaluesandDirichletboundaryconditions">Eigenvalues and Dirichlet boundary conditions</a>
        <li><a href="#Implementationdetails">Implementation details</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeEigenvalueProblemcodeclasstemplate">The <code>EigenvalueProblem</code> class template</a>
        <li><a href="#ImplementationofthecodeEigenvalueProblemcodeclass">Implementation of the <code>EigenvalueProblem</code> class</a>
      <ul>
        <li><a href="#EigenvalueProblemEigenvalueProblem">EigenvalueProblem::EigenvalueProblem</a>
        <li><a href="#EigenvalueProblemmake_grid_and_dofs">EigenvalueProblem::make_grid_and_dofs</a>
        <li><a href="#EigenvalueProblemassemble_system">EigenvalueProblem::assemble_system</a>
        <li><a href="#EigenvalueProblemsolve">EigenvalueProblem::solve</a>
        <li><a href="#EigenvalueProblemoutput_results">EigenvalueProblem::output_results</a>
        <li><a href="#EigenvalueProblemrun">EigenvalueProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Runningtheproblem">Running the problem</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-36/doc/intro.dox

 <br> 

<i>This program was contributed by Toby D. Young and Wolfgang
Bangerth.  </i>

<a name="Preamble"></a>

<a name="Preamble"></a><h1>Preamble</h1>


在这个例子中，我们要解决的问题是一个特征谱问题。特征值问题出现在广泛的问题背景中，例如在计算腔体中的电磁驻波、鼓膜的振动模式或湖泊和河口的振荡中。最神秘的应用之一可能是量子力学中静止或准静止波函数的计算。后一种应用是我们在此想研究的，尽管本方案中概述的一般技术当然同样适用于上述其他应用。

Eigenspectrum问题的一般形式是

@f{align*}
	L \Psi &= \varepsilon \Psi \qquad &&\text{in}\ \Omega\quad,
	\\
	\Psi &= 0 &&\text{on}\ \partial\Omega\quad,


@f}

其中 $\Psi=\Psi(\mathbf x)$ 上的迪里切特边界条件也可以用诺伊曼或罗宾条件代替； $L$ 是一个算子，一般也包含微分算子。

在适当的条件下，上述方程有一组解 $\Psi_\ell,\varepsilon_\ell$  ,  $\ell\in {\cal I}$  ，其中 $\cal I$ 可以是一个有限的或无限的集合（在后一种情况下，它可能是一个离散的或有时至少是部分连续的集合）。在这两种情况下，让我们注意到，我们要计算的不再只是一个单一的解决方案，而是一组解决方案（各种特征函数和相应的特征值）。从数字上寻找这类特征值问题的所有特征值（特征函数）的问题是一个艰巨的挑战。事实上，如果集合 $\cal I$ 是无限的，这个挑战当然是难以解决的。  然而，大多数时候，我们实际上只对这些值（函数）的一小部分感兴趣；幸运的是，我们将在这个教程程序中使用的SLEPc库的接口允许我们选择特征谱的哪一部分和我们想要解决的多少个解决方案。

在这个程序中，我们使用的eigenspectrum求解器是由deal.II提供的类，围绕<a href="http://www.grycap.upv.es/slepc/" target="_top">SLEPc</a>库的线性代数实现；SLEPc本身建立在<a
href="http://www.mcs.anl.gov/petsc/" target="_top">PETSc</a>库的线性代数内容之上。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


静止的量子力学的基本方程是薛定谔方程，它模拟粒子在外部势中的运动  $V(\mathbf x)$  。粒子由一个波函数 $\Psi(\mathbf x)$ 描述，它满足一个（非维度化）形式的关系

@f{align*} [-\Delta + V(\mathbf x)]
\Psi(\mathbf x) &= \varepsilon \Psi(\mathbf x) \qquad &&\text{in}\
\Omega\quad, \\ \Psi &= 0 &&\text{on}\ \partial\Omega\quad.


@f}

因此，这个粒子只能存在于一定数量的特征态中，这些特征态对应于作为该方程解的能量特征值 $\varepsilon_\ell$ 。量子力学的正统（哥本哈根）解释认为，如果一个粒子具有能量 $\varepsilon_\ell$ ，那么在位置 $\mathbf x$ 找到它的概率与 $|\Psi_\ell(\mathbf
x)|^2$ 成正比，其中 $\Psi_\ell$ 是对应于这个特征值的特征函数。

为了在数值上找到这个方程的解，即一组成对的特征值/特征函数，我们使用通常的有限元方法，将左边的方程与测试函数相乘，通过部分积分，在有限维空间中通过近似 $\Psi(\mathbf
x)\approx\Psi_h(\mathbf x)=\sum_{j}\phi_j(\mathbf x)\tilde\psi_j$ 寻找解，其中 $\tilde\psi$ 是一个扩展系数的矢量。然后，我们立即得出以下方程，将连续特征值问题离散化：@f[ \sum_j [(\nabla\phi_i,
\nabla\phi_j)+(V(\mathbf x)\phi_i,\phi_j)] \tilde{\psi}_j =
\varepsilon_h \sum_j (\phi_i, \phi_j) \tilde{\psi}_j\quad.  @f] 在矩阵和矢量符号中，这个方程然后读作。@f[ A
\tilde{\Psi} = \varepsilon_h M \tilde{\Psi} \quad, @f] 其中 $A$ 是由微分算子 $L$ 产生的刚度矩阵，而 $M$ 是质量矩阵。特征值问题的解决方案是一个特征谱  $\varepsilon_{h,\ell}$  ，以及相关的特征函数  $\Psi_\ell=\sum_j \phi_j\tilde{\psi}_j$  。




<a name="EigenvaluesandDirichletboundaryconditions"></a><h3>Eigenvalues and Dirichlet boundary conditions</h3>


在这个程序中，我们对波函数 $\Psi$ 使用了Dirichlet边界条件。从有限元代码的角度来看，这意味着只有内部的自由度是<i>freedom</i>的实数度：边界上的自由度不是自由的，但毕竟是被迫有一个零值。另一方面，有限元方法的力量和简单性主要来自于我们只是在每个单元上做同样的事情，而不需要过多考虑一个单元在哪里，它是否在一个不太精细的单元上，因此有一个悬挂的节点，或者与边界相邻。所有这些检查都会使有限元线性系统的组装变得难以忍受，甚至更难阅读。

因此，当然，当你用DoFHandler对象分配自由度时，你并不关心你列举的一些自由度是否处于迪里切特边界。他们都会得到数字。我们只需要在以后应用边界值的时候照顾到这些自由度。有两种基本的方法（要么使用 MatrixTools::apply_boundary_values() <i>after</i>装配线性系统，要么使用 AffineConstraints::distribute_local_to_global() <i>during</i>装配；更多信息见 @ref constraints "约束模块"），但两者的结果都一样：一个线性系统的总行数等于<i>all</i>自由度的数量，包括那些位于边界的自由度。然而，受迪里希特条件约束的自由度与线性系统的其他部分分开，方法是将相应的行和列清零，在对角线上放一个正条目，并在右手边放相应的迪里希特值。

如果你暂时假设我们将自由度重新编号，使迪里切特边界上的所有自由度排在最后，那么我们在求解一个有右手边的常规PDE时得到的线性系统将是这样的。

@f{align*}
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

这里，下标 $i$ 和 $b$ 分别对应于内部和边界自由度。内部自由度满足线性系统 $A_i U_i=F_i$ ，在内部产生正确的解，边界值由 $U_b = D_b^{-1} F_b$ 决定，其中 $D_b$ 是一个对角矩阵，由消除边界自由度的过程产生， $F_b$ 是以这样的方式选择的： $U_{b,j}=D_{b,jj}^{-1} F_{b,j}$ 对每个边界自由度 $j$ 都有正确的边界值。 对于好奇的人来说，矩阵 $D_b$ 的条目是将修改后的局部贡献加入全局矩阵的结果，对于局部矩阵，如果非零，对角线元素被设置为其绝对值；否则，它们被设置为对角线的绝对值的平均值。这个过程保证了 $D_b$ 的条目是正的，并且大小与对角线的其他条目相当，确保所产生的矩阵不会因为涉及到大小迥异的矩阵条目的舍入而产生不合理的精度损失。最终出现在对角线上的实际数值是很难预测的，你应该把它们当作任意的、不可预测的，但却是正的。)

对于 "常规 "的线性系统，这一切都导致了正确的解决方案。另一方面，对于特征值问题，这就不那么简单了。在那里，消除边界值会影响到我们在当前教程程序中要解决的矩阵 $A$ 和 $M$ 。消除边界值后，我们就会得到一个可以这样划分的特征值问题。

@f{align*}
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

这种形式清楚地表明，有两组特征值：我们关心的那些，以及来自分离问题的虚假特征值。

@f[
  D_A \tilde \Psi_b = \epsilon_h D_M \Psi_b.


@f]

这些特征值是虚假的，因为它们是由一个只在边界节点上操作的特征值系统产生的--这些节点不是<i>freedom</i>的实数度。当然，由于两个矩阵 $D_A,D_M$ 是对角线，我们可以准确地量化这些虚假的特征值：它们是 $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$ （其中指数 $j$ 正好对应于受迪里切特边界值约束的自由度）。

那么，如何处理这些问题呢？第一部分是识别我们的特征值求解器何时找到其中一个。为此，程序通过计算所有受限自由度上表达式 $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$ 的最小和最大，计算并打印出这些特征值所在的区间。在下面的程序中，这已经足够了：我们发现这个区间位于我们感兴趣并计算的最小特征值和相应的特征函数的集合之外，所以这里我们不需要做什么。

另一方面，我们可能会发现我们在这个程序中计算的一个特征值恰好在这个区间内，在这种情况下，我们不会立即知道它是一个虚假的还是一个真正的特征值。在这种情况下，我们可以在计算完两个矩阵后简单地缩放其中一个矩阵的对角线元素，从而将它们从特征谱的感兴趣的频率上移开。这可以通过使用以下代码来完成，确保所有假的特征值都正好等于 $1.234\cdot 10^5$  。

@code
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          stiffness_matrix.set(i, i, 1.234e5);
          mass_matrix.set(i, i, 1);
        }
@endcode

然而，这里没有采用这种策略，因为我们从程序中得到的虚假特征值恰好大于我们将计算的和感兴趣的最低的五个。




<a name="Implementationdetails"></a><h3>Implementation details</h3>


下面的程序实质上只是步骤4的一个稍加修改的版本。有所不同的是以下几点。

 <ul> 

 <li>  主类（名为  <code>EigenvalueProblem</code>  ）现在不再有单一的解向量，而是有一整套我们想要计算的各种特征函数的向量。此外， <code>main</code> 函数对这里的一切都有顶层控制，它通过 <code>SlepcInitialize</code> 和 <code>SlepFinalize</code> 同时初始化和最终确定SLEPc和PETSc的接口。 </li> 

 <li>  我们在步骤17和步骤18中使用PETSc矩阵和向量，因为这是SLEPc特征值求解器所要求的。 </li> 

 <li>  函数 <code>EigenvalueProblem::solve</code> 与教程中迄今为止所见的任何函数都完全不同，因为它不只是求解一个线性系统，而是实际求解特征值问题。它建立在SLEPc库上，更直接的是建立在类 SLEPcWrappers::SolverKrylovSchur.</li> 的deal.II SLEPc包装器上。

 <li>  我们使用ParameterHandler类来描述一些输入参数，如势的确切形式 $V({\mathbf
x})$  ，网格的全局细化步数，或我们要解决的特征值的数量。我们可以在这方面做得更多，但只限于在运行时选择一些实际的输入文件参数。为了看看在这方面可以做什么，看看 @ref step_29 "步骤-29 "和步骤-33。 </li> 

 <li>  我们使用FunctionParser类使潜在  $V(\mathbf
x)$  的运行时参数，可以在输入文件中指定为公式。 </li> 

 </ul> 

程序的其余部分以一种相当直接的方式从第4步开始。


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
 * As mentioned in the introduction, this program is essentially only a
 * slightly revised version of step-4. As a consequence, most of the following
 * include files are as used there, or at least as used already in previous
 * tutorial programs:
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/function_parser.h>
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/full_matrix.h>
 * 
 * @endcode
 * 
 * IndexSet is used to set the size of each PETScWrappers::MPI::Vector:
 * 
 * @code
 * #include <deal.II/base/index_set.h>
 * 
 * @endcode
 * 
 * PETSc appears here because SLEPc depends on this library:
 * 
 * @code
 * #include <deal.II/lac/petsc_sparse_matrix.h>
 * #include <deal.II/lac/petsc_vector.h>
 * 
 * @endcode
 * 
 * And then we need to actually import the interfaces for solvers that SLEPc
 * provides:
 * 
 * @code
 * #include <deal.II/lac/slepc_solver.h>
 * 
 * @endcode
 * 
 * We also need some standard C++:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * Finally, as in previous programs, we import all the deal.II class and
 * function names into the namespace into which everything in this program
 * will go:
 * 
 * @code
 * namespace Step36
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeEigenvalueProblemcodeclasstemplate"></a> 
 * <h3>The <code>EigenvalueProblem</code> class template</h3>
 * 

 * 
 * Following is the class declaration for the main class template. It looks
 * pretty much exactly like what has already been shown in step-4:
 * 
 * @code
 *   template <int dim>
 *   class EigenvalueProblem
 *   {
 *   public:
 *     EigenvalueProblem(const std::string &prm_file);
 *     void run();
 * 
 *   private:
 *     void         make_grid_and_dofs();
 *     void         assemble_system();
 *     unsigned int solve();
 *     void         output_results() const;
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 * @endcode
 * 
 * With these exceptions: For our eigenvalue problem, we need both a
 * stiffness matrix for the left hand side as well as a mass matrix for
 * the right hand side. We also need not just one solution function, but a
 * whole set of these for the eigenfunctions we want to compute, along
 * with the corresponding eigenvalues:
 * 
 * @code
 *     PETScWrappers::SparseMatrix             stiffness_matrix, mass_matrix;
 *     std::vector<PETScWrappers::MPI::Vector> eigenfunctions;
 *     std::vector<double>                     eigenvalues;
 * 
 * @endcode
 * 
 * And then we need an object that will store several run-time parameters
 * that we will specify in an input file:
 * 
 * @code
 *     ParameterHandler parameters;
 * 
 * @endcode
 * 
 * Finally, we will have an object that contains "constraints" on our
 * degrees of freedom. This could include hanging node constraints if we
 * had adaptively refined meshes (which we don't have in the current
 * program). Here, we will store the constraints for boundary nodes
 * $U_i=0$.
 * 
 * @code
 *     AffineConstraints<double> constraints;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeEigenvalueProblemcodeclass"></a> 
 * <h3>Implementation of the <code>EigenvalueProblem</code> class</h3>
 * 

 * 
 * 
 * <a name="EigenvalueProblemEigenvalueProblem"></a> 
 * <h4>EigenvalueProblem::EigenvalueProblem</h4>
 * 

 * 
 * First up, the constructor. The main new part is handling the run-time
 * input parameters. We need to declare their existence first, and then read
 * their values from the input file whose name is specified as an argument
 * to this function:
 * 
 * @code
 *   template <int dim>
 *   EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file)
 *     : fe(1)
 *     , dof_handler(triangulation)
 *   {
 * @endcode
 * 
 * TODO investigate why the minimum number of refinement steps required to
 * obtain the correct eigenvalue degeneracies is 6
 * 
 * @code
 *     parameters.declare_entry(
 *       "Global mesh refinement steps",
 *       "5",
 *       Patterns::Integer(0, 20),
 *       "The number of times the 1-cell coarse mesh should "
 *       "be refined globally for our computations.");
 *     parameters.declare_entry("Number of eigenvalues/eigenfunctions",
 *                              "5",
 *                              Patterns::Integer(0, 100),
 *                              "The number of eigenvalues/eigenfunctions "
 *                              "to be computed.");
 *     parameters.declare_entry("Potential",
 *                              "0",
 *                              Patterns::Anything(),
 *                              "A functional description of the potential.");
 * 
 *     parameters.parse_input(prm_file);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemmake_grid_and_dofs"></a> 
 * <h4>EigenvalueProblem::make_grid_and_dofs</h4>
 * 

 * 
 * The next function creates a mesh on the domain $[-1,1]^d$, refines it as
 * many times as the input file calls for, and then attaches a DoFHandler to
 * it and initializes the matrices and vectors to their correct sizes. We
 * also build the constraints that correspond to the boundary values
 * $u|_{\partial\Omega}=0$.
 *   

 * 
 * For the matrices, we use the PETSc wrappers. These have the ability to
 * allocate memory as necessary as non-zero entries are added. This seems
 * inefficient: we could as well first compute the sparsity pattern,
 * initialize the matrices with it, and as we then insert entries we can be
 * sure that we do not need to re-allocate memory and free the one used
 * previously. One way to do that would be to use code like this:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   DynamicSparsityPattern
 *      dsp (dof_handler.n_dofs(),
 *           dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern (dof_handler, dsp);
 *   dsp.compress ();
 *   stiffness_matrix.reinit (dsp);
 *   mass_matrix.reinit (dsp);
 * @endcode
 * </div>
 * instead of the two <code>reinit()</code> calls for the
 * stiffness and mass matrices below.
 *   

 * 
 * This doesn't quite work, unfortunately. The code above may lead to a few
 * entries in the non-zero pattern to which we only ever write zero entries;
 * most notably, this holds true for off-diagonal entries for those rows and
 * columns that belong to boundary nodes. This shouldn't be a problem, but
 * for whatever reason, PETSc's ILU preconditioner, which we use to solve
 * linear systems in the eigenvalue solver, doesn't like these extra entries
 * and aborts with an error message.
 *   

 * 
 * In the absence of any obvious way to avoid this, we simply settle for the
 * second best option, which is have PETSc allocate memory as
 * necessary. That said, since this is not a time critical part, this whole
 * affair is of no further importance.
 * 
 * @code
 *   template <int dim>
 *   void EigenvalueProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -1, 1);
 *     triangulation.refine_global(
 *       parameters.get_integer("Global mesh refinement steps"));
 *     dof_handler.distribute_dofs(fe);
 * 
 *     DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
 *     constraints.close();
 * 
 *     stiffness_matrix.reinit(dof_handler.n_dofs(),
 *                             dof_handler.n_dofs(),
 *                             dof_handler.max_couplings_between_dofs());
 *     mass_matrix.reinit(dof_handler.n_dofs(),
 *                        dof_handler.n_dofs(),
 *                        dof_handler.max_couplings_between_dofs());
 * 
 * @endcode
 * 
 * The next step is to take care of the eigenspectrum. In this case, the
 * outputs are eigenvalues and eigenfunctions, so we set the size of the
 * list of eigenfunctions and eigenvalues to be as large as we asked for
 * in the input file. When using a PETScWrappers::MPI::Vector, the Vector
 * is initialized using an IndexSet. IndexSet is used not only to resize the
 * PETScWrappers::MPI::Vector but it also associates an index in the
 * PETScWrappers::MPI::Vector with a degree of freedom (see step-40 for a
 * more detailed explanation). The function complete_index_set() creates
 * an IndexSet where every valid index is part of the set. Note that this
 * program can only be run sequentially and will throw an exception if used
 * in parallel.
 * 
 * @code
 *     IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs();
 *     eigenfunctions.resize(
 *       parameters.get_integer("Number of eigenvalues/eigenfunctions"));
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
 *       eigenfunctions[i].reinit(eigenfunction_index_set, MPI_COMM_WORLD);
 * 
 *     eigenvalues.resize(eigenfunctions.size());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemassemble_system"></a> 
 * <h4>EigenvalueProblem::assemble_system</h4>
 * 

 * 
 * Here, we assemble the global stiffness and mass matrices from local
 * contributions $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf x) \cdot
 * \nabla\varphi_j(\mathbf x) + V(\mathbf x)\varphi_i(\mathbf
 * x)\varphi_j(\mathbf x)$ and $M^K_{ij} = \int_K \varphi_i(\mathbf
 * x)\varphi_j(\mathbf x)$ respectively. This function should be immediately
 * familiar if you've seen previous tutorial programs. The only thing new
 * would be setting up an object that described the potential $V(\mathbf x)$
 * using the expression that we got from the input file. We then need to
 * evaluate this object at the quadrature points on each cell. If you've
 * seen how to evaluate function objects (see, for example the coefficient
 * in step-5), the code here will also look rather familiar.
 * 
 * @code
 *   template <int dim>
 *   void EigenvalueProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
 *     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     FunctionParser<dim> potential;
 *     potential.initialize(FunctionParser<dim>::default_variable_names(),
 *                          parameters.get("Potential"),
 *                          typename FunctionParser<dim>::ConstMap());
 * 
 *     std::vector<double> potential_values(n_q_points);
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         cell_stiffness_matrix = 0;
 *         cell_mass_matrix      = 0;
 * 
 *         potential.value_list(fe_values.get_quadrature_points(),
 *                              potential_values);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               {
 *                 cell_stiffness_matrix(i, j) +=           
 *                   (fe_values.shape_grad(i, q_point) *    
 *                      fe_values.shape_grad(j, q_point)    
 *                    +                                     
 *                    potential_values[q_point] *           
 *                      fe_values.shape_value(i, q_point) * 
 *                      fe_values.shape_value(j, q_point)   
 *                    ) *                                   
 *                   fe_values.JxW(q_point);                
 * 
 *                 cell_mass_matrix(i, j) +=              
 *                   (fe_values.shape_value(i, q_point) * 
 *                    fe_values.shape_value(j, q_point)   
 *                    ) *                                 
 *                   fe_values.JxW(q_point);              
 *               }
 * 
 * @endcode
 * 
 * Now that we have the local matrix contributions, we transfer them
 * into the global objects and take care of zero boundary constraints:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraints.distribute_local_to_global(cell_stiffness_matrix,
 *                                                local_dof_indices,
 *                                                stiffness_matrix);
 *         constraints.distribute_local_to_global(cell_mass_matrix,
 *                                                local_dof_indices,
 *                                                mass_matrix);
 *       }
 * 
 * @endcode
 * 
 * At the end of the function, we tell PETSc that the matrices have now
 * been fully assembled and that the sparse matrix representation can now
 * be compressed as no more entries will be added:
 * 
 * @code
 *     stiffness_matrix.compress(VectorOperation::add);
 *     mass_matrix.compress(VectorOperation::add);
 * 
 * 
 * @endcode
 * 
 * Before leaving the function, we calculate spurious eigenvalues,
 * introduced to the system by zero Dirichlet constraints. As
 * discussed in the introduction, the use of Dirichlet boundary
 * conditions coupled with the fact that the degrees of freedom
 * located at the boundary of the domain remain part of the linear
 * system we solve, introduces a number of spurious eigenvalues.
 * Below, we output the interval within which they all lie to
 * ensure that we can ignore them should they show up in our
 * computations.
 * 
 * @code
 *     double min_spurious_eigenvalue = std::numeric_limits<double>::max(),
 *            max_spurious_eigenvalue = -std::numeric_limits<double>::max();
 * 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *       if (constraints.is_constrained(i))
 *         {
 *           const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i);
 *           min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev);
 *           max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev);
 *         }
 * 
 *     std::cout << "   Spurious eigenvalues are all in the interval "
 *               << "[" << min_spurious_eigenvalue << ","
 *               << max_spurious_eigenvalue << "]" << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemsolve"></a> 
 * <h4>EigenvalueProblem::solve</h4>
 * 

 * 
 * This is the key new functionality of the program. Now that the system is
 * set up, here is a good time to actually solve the problem: As with other
 * examples this is done using a "solve" routine. Essentially, it works as
 * in other programs: you set up a SolverControl object that describes the
 * accuracy to which we want to solve the linear systems, and then we select
 * the kind of solver we want. Here we choose the Krylov-Schur solver of
 * SLEPc, a pretty fast and robust choice for this kind of problem:
 * 
 * @code
 *   template <int dim>
 *   unsigned int EigenvalueProblem<dim>::solve()
 *   {
 * @endcode
 * 
 * We start here, as we normally do, by assigning convergence control we
 * want:
 * 
 * @code
 *     SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9);
 *     SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control);
 * 
 * @endcode
 * 
 * Before we actually solve for the eigenfunctions and -values, we have to
 * also select which set of eigenvalues to solve for. Lets select those
 * eigenvalues and corresponding eigenfunctions with the smallest real
 * part (in fact, the problem we solve here is symmetric and so the
 * eigenvalues are purely real). After that, we can actually let SLEPc do
 * its work:
 * 
 * @code
 *     eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL);
 * 
 *     eigensolver.set_problem_type(EPS_GHEP);
 * 
 *     eigensolver.solve(stiffness_matrix,
 *                       mass_matrix,
 *                       eigenvalues,
 *                       eigenfunctions,
 *                       eigenfunctions.size());
 * 
 * @endcode
 * 
 * The output of the call above is a set of vectors and values. In
 * eigenvalue problems, the eigenfunctions are only determined up to a
 * constant that can be fixed pretty arbitrarily. Knowing nothing about
 * the origin of the eigenvalue problem, SLEPc has no other choice than to
 * normalize the eigenvectors to one in the $l_2$ (vector)
 * norm. Unfortunately this norm has little to do with any norm we may be
 * interested from a eigenfunction perspective: the $L_2(\Omega)$ norm, or
 * maybe the $L_\infty(\Omega)$ norm.
 *     

 * 
 * Let us choose the latter and rescale eigenfunctions so that they have
 * $\|\phi_i(\mathbf x)\|_{L^\infty(\Omega)}=1$ instead of
 * $\|\Phi\|_{l_2}=1$ (where $\phi_i$ is the $i$th eigen<i>function</i>
 * and $\Phi_i$ the corresponding vector of nodal values). For the $Q_1$
 * elements chosen here, we know that the maximum of the function
 * $\phi_i(\mathbf x)$ is attained at one of the nodes, so $\max_{\mathbf
 * x}\phi_i(\mathbf x)=\max_j (\Phi_i)_j$, making the normalization in the
 * $L_\infty$ norm trivial. Note that this doesn't work as easily if we
 * had chosen $Q_k$ elements with $k>1$: there, the maximum of a function
 * does not necessarily have to be attained at a node, and so
 * $\max_{\mathbf x}\phi_i(\mathbf x)\ge\max_j (\Phi_i)_j$ (although the
 * equality is usually nearly true).
 * 
 * @code
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
 *       eigenfunctions[i] /= eigenfunctions[i].linfty_norm();
 * 
 * @endcode
 * 
 * Finally return the number of iterations it took to converge:
 * 
 * @code
 *     return solver_control.last_step();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemoutput_results"></a> 
 * <h4>EigenvalueProblem::output_results</h4>
 * 

 * 
 * This is the last significant function of this program. It uses the
 * DataOut class to generate graphical output from the eigenfunctions for
 * later visualization. It works as in many of the other tutorial programs.
 *   

 * 
 * The whole collection of functions is then output as a single VTK file.
 * 
 * @code
 *   template <int dim>
 *   void EigenvalueProblem<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 * 
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i)
 *       data_out.add_data_vector(eigenfunctions[i],
 *                                std::string("eigenfunction_") +
 *                                  Utilities::int_to_string(i));
 * 
 * @endcode
 * 
 * The only thing worth discussing may be that because the potential is
 * specified as a function expression in the input file, it would be nice
 * to also have it as a graphical representation along with the
 * eigenfunctions. The process to achieve this is relatively
 * straightforward: we build an object that represents $V(\mathbf x)$ and
 * then we interpolate this continuous function onto the finite element
 * space. The result we also attach to the DataOut object for
 * visualization.
 * 
 * @code
 *     Vector<double> projected_potential(dof_handler.n_dofs());
 *     {
 *       FunctionParser<dim> potential;
 *       potential.initialize(FunctionParser<dim>::default_variable_names(),
 *                            parameters.get("Potential"),
 *                            typename FunctionParser<dim>::ConstMap());
 *       VectorTools::interpolate(dof_handler, potential, projected_potential);
 *     }
 *     data_out.add_data_vector(projected_potential, "interpolated_potential");
 * 
 *     data_out.build_patches();
 * 
 *     std::ofstream output("eigenvectors.vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemrun"></a> 
 * <h4>EigenvalueProblem::run</h4>
 * 

 * 
 * This is the function which has the top-level control over everything. It
 * is almost exactly the same as in step-4:
 * 
 * @code
 *   template <int dim>
 *   void EigenvalueProblem<dim>::run()
 *   {
 *     make_grid_and_dofs();
 * 
 *     std::cout << "   Number of active cells:       "
 *               << triangulation.n_active_cells() << std::endl
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 *     assemble_system();
 * 
 *     const unsigned int n_iterations = solve();
 *     std::cout << "   Solver converged in " << n_iterations << " iterations."
 *               << std::endl;
 * 
 *     output_results();
 * 
 *     std::cout << std::endl;
 *     for (unsigned int i = 0; i < eigenvalues.size(); ++i)
 *       std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i]
 *                 << std::endl;
 *   }
 * } // namespace Step36
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step36;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 * 
 * @endcode
 * 
 * This program can only be run in serial. Otherwise, throw an exception.
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
 *                   ExcMessage(
 *                     "This program can only be run in serial, use ./step-36"));
 * 
 *       EigenvalueProblem<2> problem("step-36.prm");
 *       problem.run();
 *     }
 * 
 * @endcode
 * 
 * All the while, we are watching out if any exceptions should have been
 * generated. If that is so, we panic...
 * 
 * @code
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
 * @endcode
 * 
 * If no exceptions are thrown, then we tell the program to stop monkeying
 * around and exit nicely:
 * 
 * @code
 *   std::cout << std::endl << "   Job done." << std::endl;
 * 
 *   return 0;
 * }
 * @endcode
examples/step-36/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Runningtheproblem"></a><h3>Running the problem</h3>


该问题的输入由一个输入文件 <code>\step-36.prm</code> 设定参数，例如，该文件可以包含以下文本。

@code
set Global mesh refinement steps         = 5
set Number of eigenvalues/eigenfunctions = 5
set Potential                            = 0
@endcode



这里，域内电势为零，我们知道特征值由 $\lambda_{(mn)}=\frac{\pi^2}{4}(m^2+n^2)$ 给出，其中 $m,n\in{\mathbb N^+}$  。特征函数是正弦和余弦，在 $m$ 和 $n$ 方向的周期为 $x$ 和 $y$ 。这与我们的程序产生的输出相匹配。

@code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       1024
   Number of degrees of freedom: 1089
   Solver converged in 67 iterations.


      Eigenvalue 0 : 4.93877
      Eigenvalue 1 : 12.3707
      Eigenvalue 2 : 12.3707
      Eigenvalue 3 : 19.8027
      Eigenvalue 4 : 24.837


   Job done.  @endcode 这些特征值正是对应于 $(m,n)=(1,1)$ 、 $(1,2)$ 和 $(2,1)$ 、 $(2,2)$ 和 $(3,1)$ 等对。相应的特征方程的可视化看起来是这样的。

 <table width="80%">
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

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


在操场上玩几个游戏总是值得的!所以这里有几个建议。

 <ul> 

 <li> 上面使用的势（称为<i>infinite well</i>，因为它是一个由无限高的墙包围的平坦势）很有趣，因为它允许有分析上已知的解决方案。然而，除此之外，它是相当无聊的。也就是说，通过在输入文件中设置不同的势来玩弄这个势是很容易的。例如，让我们假设我们想在2d中使用以下势。

@f[
  V(x,y) = \left\{
       \begin{array}{ll}


         -100 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy>0
         \\


         -5 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy\le 0
         \\
         0 & \text{otherwise}
      \end{array} \right.\quad.


@f]

换句话说，在半径为0.75的圆的两个扇面中，电位为-100，在另外两个扇面中为-5，而在圆外为零。我们可以通过在输入文件中使用以下内容来实现这一点。

@code
set Potential = if (x^2 + y^2 < 0.75^2, if (x*y > 0, -100, -5), 0)
@endcode

此外，如果我们还将网格细化程度提高一级，我们会得到以下结果。

@code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       4096
   Number of degrees of freedom: 4225


   Eigenvalue 0 : -74.2562
   Eigenvalue 1 : -72.7322
   Eigenvalue 2 : -42.7406
   Eigenvalue 3 : -42.2232
   Eigenvalue 4 : -37.0744
@endcode



输出文件还包含一个内插的势的版本，看起来像这样（注意，正如预期的那样，最低的几个特征模式的概率密度 $|\Psi(\mathbf x)|^2$ 只有在势最低的地方才是显著的，即在势的内圈的右上角和左下角部分）。

 <img src="https://www.dealii.org/images/steps/developer/step-36.mod.potential.png" alt=""> 

前五个特征函数现在是这样的。

 <table width="80%">
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

 <li> 在我们对问题的推导中，我们假设粒子被限制在一个域 $\Omega$ 中，并且在这个域的边界处，它的概率 $|\Psi|^2$ 为零。这相当于解决所有 ${\mathbb R}^d$ 上的特征值问题，并假设能量势只在 $\Omega$ 区域内是有限的，而在区域外是无限的。比较容易的是，在 $|\Psi(\mathbf x)|^2$ 的所有位置 $\mathbf x$ ， $V(\mathbf
x)=\infty$ 。那么问题来了，如果我们的势不是这种形式的，即没有一个势是无限的有界域，会发生什么？在这种情况下，可能值得只考虑一个非常大的边界域，其中 $V(\mathbf x)$ 即使不是无限的，也至少是非常大的。在这样的情况下玩一玩，探索一下当我们使计算区域越来越大时，频谱和特征函数如何变化。

 <li>  如果我们研究简单的谐波振荡器问题 $V(\mathbf x)=c|\mathbf x|^2$ 会怎样？这个势正是上一段所讨论的形式，具有超球面对称性。人们可能想用一个大的外半径的球面域，来近似于全空间问题（例如，通过引用 GridGenerator::hyper_ball). 

 <li>  上面的图显示了波函数  $\Psi(\mathbf x)$  ，但感兴趣的物理量实际上是粒子处于位置  $|\Psi(\mathbf x)|^2$  的概率密度。一些可视化程序可以从输入文件中的数据计算出衍生量，但我们也可以在创建输出文件时立即这样做。这样做的工具是DataPostprocessor类，可以和DataOut类一起使用。如何做到这一点的例子可以在步骤29和步骤33中找到。

 <li>  如果盒子里的粒子有%的内部自由度会怎样？例如，如果该粒子是一个自旋-  $1/2$  粒子？在这种情况下，我们可能要开始解决一个矢量值的问题，而不是。

 <li>  我们这里的deal.II库的实现使用PETScWrappers和SLEPcWrappers，适合在串行机器架构上运行。然而，对于更大的网格和更多的自由度，我们可能希望在并行架构上运行我们的应用程序。上述代码的并行实现在这里可能特别有用，因为广义的特征谱问题比前面大多数教程中考虑的标准问题的解决成本更高。幸运的是，修改上述程序使其符合MPI标准是一个相对简单的过程。关于如何做到这一点的简图可以在 @ref
step_17 "step-17 "中找到。

deal.II有与其中之一ARPACK（见<a
href="../../external-libs/arpack.html">the ARPACK configuration page</a>的设置说明）的接口，在ArpackSolver类中实现。下面是一个简短的、快速的概述，说明使用它需要改变什么，前提是你有一个工作的ARPACK安装，并且deal.II已经为它正确配置了（见deal.II <a href="../../readme.html" target="body">README</a>文件）。

首先，为了使用ARPACK接口，我们可以回到使用标准的deal.II矩阵和向量，所以我们首先替换PETSc和SLEPc头文件

@code
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
@endcode

与这些。

@code
#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
@endcode

ARPACK允许复杂的特征值，所以我们还需要

@code
#include <complex>
@endcode



其次，我们在主类中切换回deal.II矩阵和向量定义。

@code
    SparsityPattern                     sparsity_pattern;
    SparseMatrix<double>                stiffness_matrix, mass_matrix;
    std::vector<Vector<double> >        eigenfunctions;
    std::vector<std::complex<double>>   eigenvalues;
@endcode

并按照惯例对它们进行初始化  <code>make_grid_and_dofs()</code>  。

@code
    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());


    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    constraints.condense (sparsity_pattern);
    sparsity_pattern.compress();


    stiffness_matrix.reinit (sparsity_pattern);
    mass_matrix.reinit (sparsity_pattern);
@endcode



为了用ARPACK解决特征值问题，我们最后需要修改  <code>solve()</code>  。

@code
  template <int dim>
  unsigned int EigenvalueProblem<dim>::solve ()
  {
    SolverControl solver_control (dof_handler.n_dofs(), 1e-9);


    SparseDirectUMFPACK inverse;
    inverse.initialize (stiffness_matrix);


    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors);


    ArpackSolver eigensolver (solver_control, additional_data);
    eigensolver.solve (stiffness_matrix,
                       mass_matrix,
                       inverse,
                       eigenvalues,
                       eigenfunctions,
                       eigenvalues.size());


    for (unsigned int i=0; i<eigenfunctions.size(); ++i)
      eigenfunctions[i] /= eigenfunctions[i].linfty_norm ();


    return solver_control.last_step ();
  }
@endcode

请注意我们是如何使用精确分解（使用SparseDirectUMFPACK）作为ARPACK的预处理程序的。   </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-36.cc"
*/
