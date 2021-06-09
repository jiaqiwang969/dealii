/**
@page step_15 The step-15 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Foreword">Foreword</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Weakformulationoftheproblem">Weak formulation of the problem</a>
        <li><a href="#Questionsabouttheappropriatesolver"> Questions about the appropriate solver </a>
        <li><a href="#Choiceofsteplengthandglobalization"> Choice of step length and globalization </a>
        <li><a href="#Summaryofthealgorithmandtestcase"> Summary of the algorithm and testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a>
        <li><a href="#Boundarycondition">Boundary condition</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#MinimalSurfaceProblemMinimalSurfaceProblem">MinimalSurfaceProblem::MinimalSurfaceProblem</a>
        <li><a href="#MinimalSurfaceProblemsetup_system">MinimalSurfaceProblem::setup_system</a>
        <li><a href="#MinimalSurfaceProblemassemble_system">MinimalSurfaceProblem::assemble_system</a>
        <li><a href="#MinimalSurfaceProblemsolve">MinimalSurfaceProblem::solve</a>
        <li><a href="#MinimalSurfaceProblemrefine_mesh">MinimalSurfaceProblem::refine_mesh</a>
        <li><a href="#MinimalSurfaceProblemset_boundary_values">MinimalSurfaceProblem::set_boundary_values</a>
        <li><a href="#MinimalSurfaceProblemcompute_residual">MinimalSurfaceProblem::compute_residual</a>
        <li><a href="#MinimalSurfaceProblemdetermine_step_length">MinimalSurfaceProblem::determine_step_length</a>
        <li><a href="#MinimalSurfaceProblemoutput_results">MinimalSurfaceProblem::output_results</a>
        <li><a href="#MinimalSurfaceProblemrun">MinimalSurfaceProblem::run</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Steplengthcontrol"> Step length control </a>
        <li><a href="#Integratingmeshrefinementandnonlinearandlinearsolvers"> Integrating mesh refinement and nonlinear and linear solvers </a>
        <li><a href="#UsingautomaticdifferentiationtocomputetheJacobianmatrix"> Using automatic differentiation to compute the Jacobian matrix </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-15/doc/intro.dox

 <br> 

<i>
This program grew out of a student project by Sven Wetterauer at the
University of Heidelberg, Germany. Most of the work for this program
is by him.
</i> <br>  。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Foreword"></a><h3>Foreword</h3>


这个程序涉及到一个非线性椭圆偏微分方程的例子，即[最小表面方程](https://en.wikipedia.org/wiki/Minimal_surface)。你可以想象这个方程的解来描述一个肥皂膜所跨越的表面，这个肥皂膜被一个封闭的金属环所包围。我们想象金属丝不只是一个平面的环，实际上是弯曲的。肥皂膜的表面张力将使该表面变成最小表面。最小表面方程的解描述了这个形状，电线的垂直位移是一个边界条件。为了简单起见，我们在此假定表面可以写成图形 $u=u(x,y)$ ，尽管很明显，构建这样的情况并不难，即钢丝弯曲的方式使得表面只能局部构建成图形，而不能整体构建。

因为这个方程是非线性的，所以我们不能直接解决它。相反，我们必须使用牛顿的方法来迭代计算解决方案。

 @dealiiVideoLecture{31.5,31.55,31.6}  （  @dealiiVideoLectureSeeAlso{31.65,31.7})  ）。




<a name="Classicalformulation"></a><h3>Classical formulation</h3>


在经典的意义上，该问题以如下形式给出。


  @f{align*}


    -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.
  @f}



 $\Omega$ 是我们将导线的位置投影到 $x-y$ 空间得到的域。在这个例子中，我们选择 $\Omega$ 作为单位盘。

如上所述，我们用牛顿方法解决这个方程，在这个方法中，我们从 $(n-1)$ 个方程计算出 $n$ 个近似解，并使用阻尼参数 $\alpha^n$ 来获得更好的全局收敛行为。   @f{align*}
    F'(u^{n},\delta u^{n})&=- F(u^{n})
    \\
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}

与@f[
    F(u) \dealcoloneq -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)
  @f]

而 $F'(u,\delta u)$ 是F在 $\delta u$ 方向的导数。

@f[
  F'(u,\delta u)=\lim \limits_{\epsilon \rightarrow 0}{\frac{F(u+\epsilon \delta u)-
  F(u)}{\epsilon}}.


@f]



通过寻找 $F'(u,\delta u)$ 是什么，我们发现我们必须在每一个牛顿步骤中解决一个线性椭圆PDE， $\delta u^n$ 是的解。

  @f[


  - \nabla \cdot \left( \frac{1}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{1}{2}}}\nabla
  \delta u^{n} \right) +
  \nabla \cdot \left( \frac{\nabla u^{n} \cdot
  \nabla \delta u^{n}}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{3}{2}}} \nabla u^{n}
  \right)  =


  -\left( - \nabla \cdot \left( \frac{1}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{1}{2}}}
  \nabla u^{n} \right) \right)
  @f]



为了解决最小表面方程，我们必须重复解决这个方程，每一个牛顿步骤一次。为了解决这个问题，我们必须看一下这个问题的边界条件。假设 $u^{n}$ 已经有了正确的边界值，那么牛顿更新 $\delta u^{n}$ 的边界条件应该为零，这样才能在加入两者之后有正确的边界条件。  在第一个牛顿步骤中，我们从解 $u^{0}\equiv 0$ 开始，牛顿更新仍然要向解 $u^{1}$ 传递正确的边界条件。


综上所述，我们必须在第一步用边界条件 $\delta
u^{0}=g$ 解决上述PDE，并在接下来的所有步骤中用 $\delta u^{n}=0$ 解决。

 @note  在某种意义上，人们可能会说，如果程序已经实现了 $F(u)$ ，那么还必须实现 $F'(u,\delta)$ 就是重复的。像往常一样，重复会诱发错误，我们希望能避免它。虽然我们在这个程序中没有探讨这个问题，但我们将在下面的<a
  href="#extensions">Possibilities for extensions</a>部分的结尾处，特别是在步骤72中再来讨论这个问题。




<a name="Weakformulationoftheproblem"></a><h3>Weak formulation of the problem</h3>


从上面的强表述开始，我们通过将PDE的两边都乘以一个检验函数 $\varphi$ 并对两边进行部分积分，得到弱表述。   @f[
  \left( \nabla \varphi , \frac{1}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{1}{2}}}\nabla
  \delta u^{n} \right)-\left(\nabla \varphi ,\frac{\nabla u^{n} \cdot \nabla
  \delta u^{n}}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{3}{2}}}\nabla u^{n}  \right)
  = -\left(\nabla \varphi , \frac{1}{\left(1+|\nabla u^{n}|^{2}\right)^{\frac{1}{2}}} \nabla u^{n}
   \right).
  @f]

这里的解 $\delta u^{n}$ 是 $H^{1}(\Omega)$ 中的一个函数，受制于上面讨论的边界条件。将这一空间还原为具有基数 $\left\{
\varphi_{0},\dots , \varphi_{N-1}\right\}$ 的有限维空间，我们可以写出该解。

@f[
  \delta u^{n}=\sum_{j=0}^{N-1} \delta U_{j} \varphi_{j}.


@f]



使用基函数作为测试函数并定义 $a_{n} \dealcoloneq \frac{1}
{\sqrt{1+|\nabla u^{n}|^{2}}}$ ，我们可以重写弱的表述。

@f[
  \sum_{j=0}^{N-1}\left[ \left( \nabla \varphi_{i} , a_{n} \nabla \varphi_{j} \right) -
  \left(\nabla u^{n}\cdot \nabla \varphi_{i} , a_{n}^{3} \nabla u^{n} \cdot \nabla
  \varphi_{j} \right) \right] \cdot \delta U_{j}=-\left( \nabla \varphi_{i} , a_{n}
  \nabla u^{n}\right) \qquad \forall i=0,\dots ,N-1,


@f]



其中解  $\delta u^{n}$  是由系数  $\delta U^{n}_{j}$  给出的。这个线性方程组可以改写为。

@f[
  A^{n}\; \delta U^{n}=b^{n},


@f]



其中，矩阵 $A^{n}$ 的条目由以下方式给出。

@f[
  A^{n}_{ij} \dealcoloneq \left( \nabla \varphi_{i} , a_{n} \nabla \varphi_{j} \right) -
  \left(\nabla u^{n}\cdot \nabla \varphi_{i} , a_{n}^{3} \nabla u^{n} \cdot \nabla
  \varphi_{j} \right),


@f]



而右手边 $b^{n}$ 由以下公式给出。

@f[
  b^{n}_{i} \dealcoloneq -\left( \nabla \varphi_{i} , a_{n} \nabla u^{n}\right).


@f]






<a name="Questionsabouttheappropriatesolver"></a><h3> Questions about the appropriate solver </h3>


与上述牛顿步骤相对应的矩阵可以重新表述，以更好地显示其结构。稍微重写一下，我们可以得到它的形式是

@f[
  A_{ij}
  =
  \left(
    \nabla \varphi_i,
    B
    \nabla \varphi_j
  \right),


@f]

其中矩阵 $B$ （在 $d$ 空间维度中的大小为 $d \times d$ ）由以下表达式给出。

@f[
  B
  =
  a_n \left\{
   \mathbf I


   -
   a_n^2 [\nabla u_n] \otimes [\nabla u_n]
  \right\}
  =
  a_n \left\{
   \mathbf I


   -
  \frac{\nabla u_n}{\sqrt{1+|\nabla u^{n}|^{2}}} \otimes
  \frac{\nabla u_n}{\sqrt{1+|\nabla u^{n}|^{2}}}
  \right\}.


@f]

从这个表达式来看，显然 $B$ 是对称的，所以 $A$ 也是对称的。另一方面， $B$ 也是正定的，这使 $A$ 也具有同样的属性。这可以通过注意到向量 $v_1 =
\frac{\nabla u^n}{|\nabla u^n|}$ 是 $B$ 的特征向量，其特征值为 $\lambda_1=a_n \left(1-\frac{|\nabla u^n|^2}{1+|\nabla u^n|^2}\right) > 0$ ，而所有与 $v_1$ 相互垂直的向量 $v_2\ldots v_d$ 都是特征向量，其特征值为 $a_n$  。由于所有的特征值都是正的， $B$ 是正定的， $A$ 也是正定的。因此我们可以使用CG方法来解决牛顿步骤。矩阵 $A$ 是对称和正定的这一事实不应该令人惊讶。它是由取一个能量函数的导数而产生的算子的结果：最小表面方程只是最小化了一些非二次能量。因此，牛顿矩阵，作为标量能量的二阶导数矩阵，必须是对称的，因为与 $i$ th和 $j$ th自由度有关的导数显然应该相减。同样，如果能量函数是凸的，那么第二导数的矩阵必须是正定的，上面的直接计算只是重申了这一点。)

然而，值得注意的是，在 $\nabla u$ 变大的问题上，正定性会退化。换句话说，如果我们简单地将所有边界值乘以2，那么一阶 $u$ 和 $\nabla u$ 也将被乘以2，但结果是 $B$ 的最小特征值将变小，矩阵将变得更加涣散。更具体地说，对于 $|\nabla u^n|\rightarrow\infty$ ，我们有 $\lambda_1 \propto a_n \frac{1}{|\nabla u^n|^2}$ 而 $\lambda_2\ldots \lambda_d=a_n$ ；因此， $B$ 的条件数，也就是 $A$ 的条件数的一个乘法因子，会像 ${\cal O}(|\nabla u^n|^2)$ 那样增长）。用目前的程序很容易验证，确实将目前程序中使用的边界值乘以越来越大的数值，导致问题最终不再能用我们这里使用的简单预设条件的CG方法解决。




<a name="Choiceofsteplengthandglobalization"></a><h3> Choice of step length and globalization </h3>


如上所述，牛顿方法的工作原理是计算一个方向 $\delta u^n$ ，然后以一个步长 $0 < \alpha^n \le 1$ 执行更新 $u^{n+1} = u^{n}+\alpha^n
\delta u^{n}$  。一个常见的现象是，对于强非线性模型，如果我们总是选择 $\alpha^n=1$ ，牛顿方法就不会收敛，除非我们从一个足够接近非线性问题的解 $u$ 的初始猜测 $u^0$ 开始。在实践中，我们并不总是有这样的初始猜测，因此采取完整的牛顿步骤（即使用 $\alpha=1$ ）往往是行不通的。

因此，一个常见的策略是，当迭代 $u^n$ 离解 $u$ 还很远时，在前几步使用较小的步长，随着我们越来越接近，使用较大的 $\alpha^n$ 值，直到最后我们可以开始使用全步 $\alpha^n=1$ ，因为我们已经足够接近解。当然，问题是如何选择 $\alpha^n$ 。基本上有两种广泛使用的方法：直线搜索和信任区域方法。

在这个程序中，我们总是简单地选择步长等于0.1。这确保了在手头的测试案例中，我们确实得到了收敛，尽管很明显，由于没有最终恢复到全步长，我们放弃了使牛顿方法如此吸引人的快速、二次收敛性。很明显，如果这个程序是为了解决更多的现实问题，我们最终必须解决这个问题。我们将在<a href="#Results">results section</a>中对这个问题进行更多的评论，并在第77步中使用一个更好的方法。




<a name="Summaryofthealgorithmandtestcase"></a><h3> Summary of the algorithm and testcase </h3>


总的来说，我们这里的程序在许多方面与step-6并无不同。主类的布局基本上是相同的。另一方面， <code>run()</code> 函数中的驱动算法是不同的，工作原理如下。<ol>  <li>  从函数 $u^{0}\equiv 0$ 开始，以这样的方式修改它，使 $u^0$ 的值沿边界等于正确的边界值 $g$ （这发生在 <code>MinimalSurfaceProblem::set_boundary_values</code> ）。设 $n=0$  .   </li> 

 <li>  通过解决 $A^{n}\;\delta
  U^{n}=b^{n}$ 系统和 $\delta u^{n}=0$ 上的边界条件，计算牛顿更新。   </li> 

 <li>  计算一个步长  $\alpha^n$  。在这个程序中，我们总是设置  $\alpha^n=0.1$  。为了便于以后的扩展，这发生在一个自己的函数中，即 <code>MinimalSurfaceProblem::determine_step_length</code>  中。 (总是选择 $\alpha^n=0.1$ 的策略当然不是最优的--我们应该选择一个对给定搜索方向有效的步长--但这需要做一些工作来实现。最后，我们把这类事情留给外部软件包：step-77就是这样做的）。)   </li> 

 <li>  解的新近似值由  $u^{n+1}=u^{n}+\alpha^n \delta u^{n}$  给出。   </li> 

 <li>  如果 $n$ 是5的倍数，则细化网格，将解 $u^{n+1}$ 转移到新的网格中，并设置 $u^{n+1}$ 的值，以便沿边界有 $u^{n+1}|_{\partial\Gamma}=g$ （同样在 <code>MinimalSurfaceProblem::set_boundary_values</code> ）。请注意，这并不是自动保证的，即使在网格细化之前我们已经有了 $u^{n+1}|_{\partial\Gamma}=g$ ，因为网格细化会在网格中增加新的节点，当我们把旧的解决方案带到新的网格时，必须把旧的解决方案插值到新的节点。我们通过插值选择的数值可能接近于精确的边界条件，但一般来说，并不是正确的数值。   </li> 

 <li>  设置 $n\leftarrow n+1$ 并转到步骤2。   </li>   </ol> 。

我们解决的测试案例选择如下。我们寻求在单位盘 $\Omega=\{\mathbf x: \|\mathbf
x\|<1\}\subset {\mathbb R}^2$ 上找到最小曲面的解决方案，其中曲面沿边界达到 $u(x,y)|{\partial\Omega} = g(x,y) \dealcoloneq \sin(2 \pi (x+y))$ 的数值。


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
 * The first few files have already been covered in previous examples and will
 * thus not be further commented on.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * We will use adaptive mesh refinement between Newton iterations. To do so,
 * we need to be able to work with a solution on the new mesh, although it was
 * computed on the old one. The SolutionTransfer class transfers the solution
 * from the old to the new mesh:
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * @endcode
 * 
 * We then open a namespace for this program and import everything from the
 * dealii namespace into it, as in previous programs:
 * 
 * @code
 * namespace Step15
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class template</h3>
 * 

 * 
 * The class template is basically the same as in step-6.  Three additions
 * are made:
 * - There are two solution vectors, one for the Newton update
 * $\delta u^n$, and one for the current iterate $u^n$.
 * - The <code>setup_system</code> function takes an argument that denotes
 * whether this is the first time it is called or not. The difference is
 * that the first time around we need to distribute the degrees of freedom
 * and set the solution vector for $u^n$ to the correct size. The following
 * times, the function is called after we have already done these steps as
 * part of refining the mesh in <code>refine_mesh</code>.
 * - We then also need new functions: <code>set_boundary_values()</code>
 * takes care of setting the boundary values on the solution vector
 * correctly, as discussed at the end of the
 * introduction. <code>compute_residual()</code> is a function that computes
 * the norm of the nonlinear (discrete) residual. We use this function to
 * monitor convergence of the Newton iteration. The function takes a step
 * length $\alpha^n$ as argument to compute the residual of $u^n + \alpha^n
 * \; \delta u^n$. This is something one typically needs for step length
 * control, although we will not use this feature here. Finally,
 * <code>determine_step_length()</code> computes the step length $\alpha^n$
 * in each Newton iteration. As discussed in the introduction, we here use a
 * fixed step length and leave implementing a better strategy as an
 * exercise. (step-77 does this differently: It simply uses an
 * external package for the whole solution process, and a good
 * line search strategy is part of what that package provides.)
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class MinimalSurfaceProblem
 *   {
 *   public:
 *     MinimalSurfaceProblem();
 *     void run();
 * 
 *   private:
 *     void   setup_system(const bool initial_step);
 *     void   assemble_system();
 *     void   solve();
 *     void   refine_mesh();
 *     void   set_boundary_values();
 *     double compute_residual(const double alpha) const;
 *     double determine_step_length() const;
 *     void   output_results(const unsigned int refinement_cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     DoFHandler<dim> dof_handler;
 *     FE_Q<dim>       fe;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> current_solution;
 *     Vector<double> newton_update;
 *     Vector<double> system_rhs;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="Boundarycondition"></a> 
 * <h3>Boundary condition</h3>
 * 

 * 
 * The boundary condition is implemented just like in step-4.  It is chosen
 * as $g(x,y)=\sin(2 \pi (x+y))$:
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> &p,
 *                                     const unsigned int /*component*/) const
 *   {
 *     return std::sin(2 * numbers::PI * (p[0] + p[1]));
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="MinimalSurfaceProblemMinimalSurfaceProblem"></a> 
 * <h4>MinimalSurfaceProblem::MinimalSurfaceProblem</h4>
 * 

 * 
 * The constructor and destructor of the class are the same as in the first
 * few tutorials.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
 *     : dof_handler(triangulation)
 *     , fe(2)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemsetup_system"></a> 
 * <h4>MinimalSurfaceProblem::setup_system</h4>
 * 

 * 
 * As always in the setup-system function, we setup the variables of the
 * finite element method. There are same differences to step-6, because
 * there we start solving the PDE from scratch in every refinement cycle
 * whereas here we need to take the solution from the previous mesh onto the
 * current mesh. Consequently, we can't just reset solution vectors. The
 * argument passed to this function thus indicates whether we can
 * distributed degrees of freedom (plus compute constraints) and set the
 * solution vector to zero or whether this has happened elsewhere already
 * (specifically, in <code>refine_mesh()</code>).
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step)
 *   {
 *     if (initial_step)
 *       {
 *         dof_handler.distribute_dofs(fe);
 *         current_solution.reinit(dof_handler.n_dofs());
 * 
 *         hanging_node_constraints.clear();
 *         DoFTools::make_hanging_node_constraints(dof_handler,
 *                                                 hanging_node_constraints);
 *         hanging_node_constraints.close();
 *       }
 * 
 * 
 * @endcode
 * 
 * The remaining parts of the function are the same as in step-6.
 * 

 * 
 * 
 * @code
 *     newton_update.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 *     hanging_node_constraints.condense(dsp);
 * 
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemassemble_system"></a> 
 * <h4>MinimalSurfaceProblem::assemble_system</h4>
 * 

 * 
 * This function does the same as in the previous tutorials except that now,
 * of course, the matrix and right hand side functions depend on the
 * previous iteration's solution. As discussed in the introduction, we need
 * to use zero boundary values for the Newton updates; we compute them at
 * the end of this function.
 *   

 * 
 * The top of the function contains the usual boilerplate code, setting up
 * the objects that allow us to evaluate shape functions at quadrature
 * points and temporary storage locations for the local matrices and
 * vectors, as well as for the gradients of the previous solution at the
 * quadrature points. We then start the loop over all cells:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::assemble_system()
 *   {
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_gradients | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<Tensor<1, dim>> old_solution_gradients(n_q_points);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * For the assembly of the linear system, we have to obtain the values
 * of the previous solution's gradients at the quadrature
 * points. There is a standard way of doing this: the
 * FEValues::get_function_gradients function takes a vector that
 * represents a finite element field defined on a DoFHandler, and
 * evaluates the gradients of this field at the quadrature points of the
 * cell with which the FEValues object has last been reinitialized.
 * The values of the gradients at all quadrature points are then written
 * into the second argument:
 * 
 * @code
 *         fe_values.get_function_gradients(current_solution,
 *                                          old_solution_gradients);
 * 
 * @endcode
 * 
 * With this, we can then do the integration loop over all quadrature
 * points and shape functions.  Having just computed the gradients of
 * the old solution in the quadrature points, we are able to compute
 * the coefficients $a_{n}$ in these points.  The assembly of the
 * system itself then looks similar to what we always do with the
 * exception of the nonlinear terms, as does copying the results from
 * the local objects into the global ones:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double coeff =
 *               1.0 / std::sqrt(1 + old_solution_gradients[q] *
 *                                     old_solution_gradients[q]);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                   cell_matrix(i, j) +=
 *                     (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
 *                        * coeff                         //   * a_n
 *                        * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
 *                       -                                //  -
 *                       (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
 *                        * coeff * coeff * coeff         //   * a_n^3
 *                        * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
 *                           * old_solution_gradients[q]) //      * \nabla u_n)
 *                        * old_solution_gradients[q]))   //   * \nabla u_n)))
 *                      * fe_values.JxW(q));              // * dx
 * 
 *                 cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
 *                                 * coeff                     // * a_n
 *                                 * old_solution_gradients[q] // * u_n
 *                                 * fe_values.JxW(q));        // * dx
 *               }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               system_matrix.add(local_dof_indices[i],
 *                                 local_dof_indices[j],
 *                                 cell_matrix(i, j));
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *           }
 *       }
 * 
 * @endcode
 * 
 * Finally, we remove hanging nodes from the system and apply zero
 * boundary values to the linear system that defines the Newton updates
 * $\delta u^n$:
 * 
 * @code
 *     hanging_node_constraints.condense(system_matrix);
 *     hanging_node_constraints.condense(system_rhs);
 * 
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(),
 *                                              boundary_values);
 *     MatrixTools::apply_boundary_values(boundary_values,
 *                                        system_matrix,
 *                                        newton_update,
 *                                        system_rhs);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemsolve"></a> 
 * <h4>MinimalSurfaceProblem::solve</h4>
 * 

 * 
 * The solve function is the same as always. At the end of the solution
 * process we update the current solution by setting
 * $u^{n+1}=u^n+\alpha^n\;\delta u^n$.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(system_rhs.size(),
 *                                  system_rhs.l2_norm() * 1e-6);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     solver.solve(system_matrix, newton_update, system_rhs, preconditioner);
 * 
 *     hanging_node_constraints.distribute(newton_update);
 * 
 *     const double alpha = determine_step_length();
 *     current_solution.add(alpha, newton_update);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemrefine_mesh"></a> 
 * <h4>MinimalSurfaceProblem::refine_mesh</h4>
 * 

 * 
 * The first part of this function is the same as in step-6... However,
 * after refining the mesh we have to transfer the old solution to the new
 * one which we do with the help of the SolutionTransfer class. The process
 * is slightly convoluted, so let us describe it in detail:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::refine_mesh()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       current_solution,
 *       estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 * 
 * @endcode
 * 
 * Then we need an additional step: if, for example, you flag a cell that
 * is once more refined than its neighbor, and that neighbor is not
 * flagged for refinement, we would end up with a jump of two refinement
 * levels across a cell interface.  To avoid these situations, the library
 * will silently also have to refine the neighbor cell once. It does so by
 * calling the Triangulation::prepare_coarsening_and_refinement function
 * before actually doing the refinement and coarsening.  This function
 * flags a set of additional cells for refinement or coarsening, to
 * enforce rules like the one-hanging-node rule.  The cells that are
 * flagged for refinement and coarsening after calling this function are
 * exactly the ones that will actually be refined or coarsened. Usually,
 * you don't have to do this by hand
 * (Triangulation::execute_coarsening_and_refinement does this for
 * you). However, we need to initialize the SolutionTransfer class and it
 * needs to know the final set of cells that will be coarsened or refined
 * in order to store the data from the old mesh and transfer to the new
 * one. Thus, we call the function by hand:
 * 
 * @code
 *     triangulation.prepare_coarsening_and_refinement();
 * 
 * @endcode
 * 
 * With this out of the way, we initialize a SolutionTransfer object with
 * the present DoFHandler and attach the solution vector to it, followed
 * by doing the actual refinement and distribution of degrees of freedom
 * on the new mesh
 * 
 * @code
 *     SolutionTransfer<dim> solution_transfer(dof_handler);
 *     solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * Finally, we retrieve the old solution interpolated to the new
 * mesh. Since the SolutionTransfer function does not actually store the
 * values of the old solution, but rather indices, we need to preserve the
 * old solution vector until we have gotten the new interpolated
 * values. Thus, we have the new values written into a temporary vector,
 * and only afterwards write them into the solution vector object:
 * 
 * @code
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     solution_transfer.interpolate(current_solution, tmp);
 *     current_solution = tmp;
 * 
 * @endcode
 * 
 * On the new mesh, there are different hanging nodes, for which we have to
 * compute constraints again, after throwing away previous content of the
 * object. To be on the safe side, we should then also make sure that the
 * current solution's vector entries satisfy the hanging node constraints
 * (see the discussion in the documentation of the SolutionTransfer class
 * for why this is necessary). We could do this by calling
 * `hanging_node_constraints.distribute(current_solution)` explicitly; we
 * omit this step because this will happen at the end of the call to
 * `set_boundary_values()` below, and it is not necessary to do it twice.
 * 
 * @code
 *     hanging_node_constraints.clear();
 * 
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 * @endcode
 * 
 * Once we have the interpolated solution and all information about
 * hanging nodes, we have to make sure that the $u^n$ we now have
 * actually has the correct boundary values. As explained at the end of
 * the introduction, this is not automatically the case even if the
 * solution before refinement had the correct boundary values, and so we
 * have to explicitly make sure that it now has:
 * 
 * @code
 *     set_boundary_values();
 * 
 * @endcode
 * 
 * We end the function by updating all the remaining data structures,
 * indicating to <code>setup_dofs()</code> that this is not the first
 * go-around and that it needs to preserve the content of the solution
 * vector:
 * 
 * @code
 *     setup_system(false);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemset_boundary_values"></a> 
 * <h4>MinimalSurfaceProblem::set_boundary_values</h4>
 * 

 * 
 * The next function ensures that the solution vector's entries respect the
 * boundary values for our problem.  Having refined the mesh (or just
 * started computations), there might be new nodal points on the
 * boundary. These have values that are simply interpolated from the
 * previous mesh in `refine_mesh()`, instead of the correct boundary
 * values. This is fixed up by setting all boundary nodes of the current
 * solution vector explicit to the right value.
 *   

 * 
 * There is one issue we have to pay attention to, though: If we have
 * a hanging node right next to a new boundary node, then its value
 * must also be adjusted to make sure that the finite element field
 * remains continuous. This is what the call in the last line of this
 * function does.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::set_boundary_values()
 *   {
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              BoundaryValues<dim>(),
 *                                              boundary_values);
 *     for (auto &boundary_value : boundary_values)
 *       current_solution(boundary_value.first) = boundary_value.second;
 * 
 *     hanging_node_constraints.distribute(current_solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemcompute_residual"></a> 
 * <h4>MinimalSurfaceProblem::compute_residual</h4>
 * 

 * 
 * In order to monitor convergence, we need a way to compute the norm of the
 * (discrete) residual, i.e., the norm of the vector
 * $\left<F(u^n),\varphi_i\right>$ with $F(u)=-\nabla \cdot \left(
 * \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)$ as discussed in the
 * introduction. It turns out that (although we don't use this feature in
 * the current version of the program) one needs to compute the residual
 * $\left<F(u^n+\alpha^n\;\delta u^n),\varphi_i\right>$ when determining
 * optimal step lengths, and so this is what we implement here: the function
 * takes the step length $\alpha^n$ as an argument. The original
 * functionality is of course obtained by passing a zero as argument.
 *   

 * 
 * In the function below, we first set up a vector for the residual, and
 * then a vector for the evaluation point $u^n+\alpha^n\;\delta u^n$. This
 * is followed by the same boilerplate code we use for all integration
 * operations:
 * 
 * @code
 *   template <int dim>
 *   double MinimalSurfaceProblem<dim>::compute_residual(const double alpha) const
 *   {
 *     Vector<double> residual(dof_handler.n_dofs());
 * 
 *     Vector<double> evaluation_point(dof_handler.n_dofs());
 *     evaluation_point = current_solution;
 *     evaluation_point.add(alpha, newton_update);
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_gradients | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double>              cell_residual(dofs_per_cell);
 *     std::vector<Tensor<1, dim>> gradients(n_q_points);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_residual = 0;
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * The actual computation is much as in
 * <code>assemble_system()</code>. We first evaluate the gradients of
 * $u^n+\alpha^n\,\delta u^n$ at the quadrature points, then compute
 * the coefficient $a_n$, and then plug it all into the formula for
 * the residual:
 * 
 * @code
 *         fe_values.get_function_gradients(evaluation_point, gradients);
 * 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double coeff =
 *               1. / std::sqrt(1 + gradients[q] * gradients[q]);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_residual(i) -= (fe_values.shape_grad(i, q) // \nabla \phi_i
 *                                    * coeff                    // * a_n
 *                                    * gradients[q]             // * u_n
 *                                    * fe_values.JxW(q));       // * dx
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           residual(local_dof_indices[i]) += cell_residual(i);
 *       }
 * 
 * @endcode
 * 
 * At the end of this function we also have to deal with the hanging node
 * constraints and with the issue of boundary values. With regard to the
 * latter, we have to set to zero the elements of the residual vector for
 * all entries that correspond to degrees of freedom that sit at the
 * boundary. The reason is that because the value of the solution there is
 * fixed, they are of course no "real" degrees of freedom and so, strictly
 * speaking, we shouldn't have assembled entries in the residual vector
 * for them. However, as we always do, we want to do exactly the same
 * thing on every cell and so we didn't not want to deal with the question
 * of whether a particular degree of freedom sits at the boundary in the
 * integration above. Rather, we will simply set to zero these entries
 * after the fact. To this end, we need to determine which degrees
 * of freedom do in fact belong to the boundary and then loop over all of
 * those and set the residual entry to zero. This happens in the following
 * lines which we have already seen used in step-11, using the appropriate
 * function from namespace DoFTools:
 * 
 * @code
 *     hanging_node_constraints.condense(residual);
 * 
 *     for (types::global_dof_index i :
 *          DoFTools::extract_boundary_dofs(dof_handler))
 *       residual(i) = 0;
 * 
 * @endcode
 * 
 * At the end of the function, we return the norm of the residual:
 * 
 * @code
 *     return residual.l2_norm();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemdetermine_step_length"></a> 
 * <h4>MinimalSurfaceProblem::determine_step_length</h4>
 * 

 * 
 * As discussed in the introduction, Newton's method frequently does not
 * converge if we always take full steps, i.e., compute $u^{n+1}=u^n+\delta
 * u^n$. Rather, one needs a damping parameter (step length) $\alpha^n$ and
 * set $u^{n+1}=u^n+\alpha^n\delta u^n$. This function is the one called
 * to compute $\alpha^n$.
 *   

 * 
 * Here, we simply always return 0.1. This is of course a sub-optimal
 * choice: ideally, what one wants is that the step size goes to one as we
 * get closer to the solution, so that we get to enjoy the rapid quadratic
 * convergence of Newton's method. We will discuss better strategies below
 * in the results section, and step-77 also covers this aspect.
 * 
 * @code
 *   template <int dim>
 *   double MinimalSurfaceProblem<dim>::determine_step_length() const
 *   {
 *     return 0.1;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemoutput_results"></a> 
 * <h4>MinimalSurfaceProblem::output_results</h4>
 * 

 * 
 * This last function to be called from `run()` outputs the current solution
 * (and the Newton update) in graphical form as a VTU file. It is entirely the
 * same as what has been used in previous tutorials.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::output_results(
 *     const unsigned int refinement_cycle) const
 *   {
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(current_solution, "solution");
 *     data_out.add_data_vector(newton_update, "update");
 *     data_out.build_patches();
 * 
 *     const std::string filename =
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu";
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemrun"></a> 
 * <h4>MinimalSurfaceProblem::run</h4>
 * 

 * 
 * In the run function, we build the first grid and then have the top-level
 * logic for the Newton iteration.
 *   

 * 
 * As described in the introduction, the domain is the unit disk around
 * the origin, created in the same way as shown in step-6. The mesh is
 * globally refined twice followed later on by several adaptive cycles.
 *   

 * 
 * Before starting the Newton loop, we also need to do a bit of
 * setup work: We need to create the basic data structures and
 * ensure that the first Newton iterate already has the correct
 * boundary values, as discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_ball(triangulation);
 *     triangulation.refine_global(2);
 * 
 *     setup_system(/*first time=*/true);
 *     set_boundary_values();
 * 
 * @endcode
 * 
 * The Newton iteration starts next. We iterate until the (norm of the)
 * residual computed at the end of the previous iteration is less than
 * $10^{-3}$, as checked at the end of the `do { ... } while` loop that
 * starts here. Because we don't have a reasonable value to initialize
 * the variable, we just use the largest value that can be represented
 * as a `double`.
 * 
 * @code
 *     double       last_residual_norm = std::numeric_limits<double>::max();
 *     unsigned int refinement_cycle   = 0;
 *     do
 *       {
 *         std::cout << "Mesh refinement step " << refinement_cycle << std::endl;
 * 
 *         if (refinement_cycle != 0)
 *           refine_mesh();
 * 
 * @endcode
 * 
 * On every mesh we do exactly five Newton steps. We print the initial
 * residual here and then start the iterations on this mesh.
 *         

 * 
 * In every Newton step the system matrix and the right hand side have
 * to be computed first, after which we store the norm of the right
 * hand side as the residual to check against when deciding whether to
 * stop the iterations. We then solve the linear system (the function
 * also updates $u^{n+1}=u^n+\alpha^n\;\delta u^n$) and output the
 * norm of the residual at the end of this Newton step.
 *         

 * 
 * After the end of this loop, we then also output the solution on the
 * current mesh in graphical form and increment the counter for the
 * mesh refinement cycle.
 * 
 * @code
 *         std::cout << "  Initial residual: " << compute_residual(0) << std::endl;
 * 
 *         for (unsigned int inner_iteration = 0; inner_iteration < 5;
 *              ++inner_iteration)
 *           {
 *             assemble_system();
 *             last_residual_norm = system_rhs.l2_norm();
 * 
 *             solve();
 * 
 *             std::cout << "  Residual: " << compute_residual(0) << std::endl;
 *           }
 * 
 *         output_results(refinement_cycle);
 * 
 *         ++refinement_cycle;
 *         std::cout << std::endl;
 *       }
 *     while (last_residual_norm > 1e-3);
 *   }
 * } // namespace Step15
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h4>The main function</h4>
 * 

 * 
 * Finally the main function. This follows the scheme of all other main
 * functions:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step15;
 * 
 *       MinimalSurfaceProblem<2> laplace_problem_2d;
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
 *   return 0;
 * }
 * @endcode
examples/step-15/doc/results.dox



<a name="Results"></a><h1>Results</h1>



程序的输出看起来如下。

@code
Mesh refinement step 0
  Initial residual: 1.53143
  Residual: 1.08746
  Residual: 0.966748
  Residual: 0.859602
  Residual: 0.766462
  Residual: 0.685475


Mesh refinement step 1
  Initial residual: 0.868959
  Residual: 0.762125
  Residual: 0.677792
  Residual: 0.605762
  Residual: 0.542748
  Residual: 0.48704


Mesh refinement step 2
  Initial residual: 0.426445
  Residual: 0.382731
  Residual: 0.343865
  Residual: 0.30918
  Residual: 0.278147
  Residual: 0.250327


Mesh refinement step 3
  Initial residual: 0.282026
  Residual: 0.253146
  Residual: 0.227414
  Residual: 0.20441
  Residual: 0.183803
  Residual: 0.165319


Mesh refinement step 4
  Initial residual: 0.154404
  Residual: 0.138723
  Residual: 0.124694
  Residual: 0.112124
  Residual: 0.100847
  Residual: 0.0907222


....
@endcode



很明显，该方案会收敛，即使不是非常快。我们将在下面讨论加速该方法的策略。

我们可以在每一组五次牛顿迭代之后，即在我们近似解决方案的每一个网格上，直观地看到解决方案。这就产生了以下一组图像。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_15_solution_1.png" alt="带等高线的零周期后的解决方案。" width="230" height="273"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_15_solution_2.png" alt="带等高线的一个周期后的解决方案。" width="230" height="273"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_15_solution_3.png" alt="带轮廓线的两个周期后的解决方案。" width="230" height="273"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_15_solution_4.png" alt="带轮廓线的三个周期后的解决方案。" width="230" height="273"> </div> <div>

可以清楚地看到，每次细化后的解决方案都能使表面最小化。解决方案收敛于人们想象中的肥皂泡，它位于一个像边界一样弯曲的线环内。同样可以看出，每次细化后，边界是如何被平滑化的。在粗略的网格上，边界看起来并不像正弦，而网格越细越像。

网格主要是在边界附近被细化，在那里解的增加或减少很强烈，而在域的内部则被粗化，在那里没有什么有趣的事情发生，因为解没有什么变化。这里显示的是第九个解和网格。

<div class="onecolumn" style="width: 60%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_15_solution_9.png" alt="第九个周期的网格和解决方案与等高线。" width="507" height="507"> </div> </div>




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


该程序显示了一个非线性、静止问题的求解器的基本结构。然而，它的收敛速度不是特别快，这是有原因的。

- 该程序总是采取0.1的步长。这就排除了牛顿方法通常选择的快速、二次收敛。

- 它没有将非线性迭代与网格细化迭代联系起来。

很明显，一个更好的方案必须解决这两点。我们将在下文中讨论它们。




<a name="Steplengthcontrol"></a><h4> Step length control </h4>


牛顿方法有两个众所周知的特性。

- 它可能不会从任意选择的起点收敛。相反，一个起点必须足够接近解决方案以保证收敛。然而，我们可以通过使用 <i>step length</i> 0<  $\alpha^n\le
  1$  的阻尼迭代来扩大牛顿方法的收敛区域。

- 如果(i)步长选择为 $\alpha^n=1$ ，并且(ii)事实上在选择步长的情况下，它表现出快速收敛的二次方阶。

这两个观察的结果是，一个成功的策略是为初始迭代选择 $\alpha^n<1$ ，直到迭代已经足够接近，允许以全步长收敛，这时我们要切换到 $\alpha^n=1$  。问题是如何以自动方式选择 $\alpha^n$ ，以满足这些标准。

我们不想在这里回顾关于这个主题的文献，只是简单地提到有两种基本的方法来解决这个问题：回溯线搜索和信任区域方法。前者更广泛地用于偏微分方程，基本上是这样做的。

- 计算一个搜索方向

- 看看 $u^n + \alpha^n\;\delta u^n$ 与 $\alpha^n=1$ 产生的残差是否比 $u^n$ 单独产生的残差 "大大减少"。

- 如果是这样，那么就采取  $\alpha^n=1$  。

- 如果不是，用  $\alpha^n=2/3$  试试残差是否 "大大缩小"。

- 如果是这样，则取 $\alpha^n=2/3$  。

- 如果不是，用 $\alpha^n=(2/3)^2$ 试试残差是否 "大大缩小"。

- 等等。当然，除了上面选择的 $2/3,
(2/3)^2, \ldots$ ，我们还可以选择其他因素 $r, r^2, \ldots$ ，用于 $0<r<1$  。很明显，"回溯 "一词的来源是：我们尝试一个长的步骤，但如果不成功，我们就尝试一个更短的步骤，越来越短的步骤，等等。函数 <code>determine_step_length()</code> 的编写方式正是为了支持这种用例。

我们是否接受一个特定的步长 $\alpha^n$ 取决于我们如何定义 "大大小于"。有很多方法，但不详细介绍，我们只说最常见的是使用沃尔夫和阿米约-戈尔德斯坦条件。对于这些，人们可以证明如下。

- 总有一个步长 $\alpha^n$ 可以满足条件，也就是说，只要问题是凸的，迭代就不会卡住。

- 如果我们足够接近解决方案，那么条件允许 $\alpha^n=1$  ，从而实现二次收敛。

我们在此不再赘述，而是将这种算法的实现作为一个练习。然而，我们注意到，如果实施得当，大多数合理的非线性问题可以在5到15次牛顿迭代中得到解决，达到工程精度&mdash；比我们目前版本的程序所需要的次数要少得多，这是一个普遍现象。

关于包括回溯在内的全局化方法的更多细节，例如可以在  @cite GNS08  和  @cite NW99  找到。

然而，非常值得一提的是，在实践中，高效非线性求解器的实现与高效有限元方法的实现一样复杂。我们不应该试图通过自己实现所有的必要步骤来重新发明车轮。在 LineMinimization::line_search() 函数中已经有了大量的拼图，可以用来实现这一目的。但是，相反，就像在deal.II等库上构建有限元求解器一样，人们应该在[SUNDIALS](https://computing.llnl.gov/projects/sundials)等库上构建非线性求解器。事实上，deal.II有与SUNDIALS的接口，特别是通过 SUNDIALS::KINSOL 类与它的非线性求解器子包KINSOL的接口。将目前的问题建立在该接口上并不十分困难--事实上，这正是step-77所做的。




<a name="Integratingmeshrefinementandnonlinearandlinearsolvers"></a><h4> Integrating mesh refinement and nonlinear and linear solvers </h4>


我们目前在每个网格上正好做了5次迭代。但这是最优的吗？人们可以提出以下问题。

- 也许在初始网格上做更多的迭代是值得的，因为那里的计算很便宜。

- 另一方面，我们不希望在每个网格上做太多的迭代：是的，我们可以在每个网格上将残差驱动到零，但这只意味着非线性迭代误差远远小于离散化误差。

- 我们应该用更高还是更低的精度来解决每个牛顿步骤中的线性系统？

最终，这归结为我们需要将当前网格上的离散化误差与我们希望在特定网格上通过牛顿迭代实现的非线性残差，以及我们希望在每个牛顿迭代中通过CG方法实现的线性迭代结合起来。

如何做到这一点，同样不是完全微不足道的，我们再次将其作为未来的练习。




<a name="UsingautomaticdifferentiationtocomputetheJacobianmatrix"></a><h4> Using automatic differentiation to compute the Jacobian matrix </h4>


正如介绍中所概述的，当解决一个形式为@f[
    F(u) \dealcoloneq


    -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right)
    = 0
  @f]的非线性问题时

我们使用牛顿迭代，要求我们反复解决线性偏微分方程@f{align*}
    F'(u^{n},\delta u^{n}) &=- F(u^{n})
  @f}。

这样，我们就可以计算出更新@f{align*}
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}。

与牛顿步骤的解 $\delta u^{n}$ 。对于这里的问题，我们可以用手计算导数 $F'(u,\delta u)$ ，得到@f[
  F'(u,\delta u)
  =


  - \nabla \cdot \left( \frac{1}{\left(1+|\nabla u|^{2}\right)^{\frac{1}{2}}}\nabla
  \delta u \right) +
  \nabla \cdot \left( \frac{\nabla u \cdot
  \nabla \delta u}{\left(1+|\nabla u|^{2}\right)^{\frac{3}{2}}} \nabla u
  \right).
  @f] 。

但这已经是一个相当大的表达方式了，无论是推导还是实现都很麻烦。在某种意义上，它也是重复的。如果我们在代码的某个地方实现了 $F(u)$ 是什么，那么 $F'(u,\delta u)$ 就不是一个独立的信息，而是至少在原则上计算机应该能够自己推断出来的东西。如果这真的能发生，那不是很好吗？也就是说，如果我们真的只需要实现 $F(u)$ ，而 $F'(u,\delta u)$ 是以某种方式隐含完成的，那不是很好吗？这实际上是可能的，并以 "自动微分 "的名义运行。步骤-71讨论了这个概念的一般术语，步骤-72说明了如何在实践中应用于我们在这里考虑的问题。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-15.cc"
*/
