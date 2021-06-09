/**
@page step_52 The step-52 tutorial program
This tutorial depends on step-26.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemstatement">Problem statement</a>
        <li><a href="#RungeKuttamethods">Runge-Kutta methods</a>
      <ul>
        <li><a href="#ExplicitRungeKuttamethods">Explicit Runge-Kutta methods</a>
        <li><a href="#EmbeddedRungeKuttamethods">Embedded Runge-Kutta methods</a>
        <li><a href="#ImplicitRungeKuttamethods">Implicit Runge-Kutta methods</a>
      </ul>
        <li><a href="#Spatiallydiscreteformulation">Spatially discrete formulation</a>
        <li><a href="#Notesonthetestcase">Notes on the testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeDiffusioncodeclass">The <code>Diffusion</code> class</a>
      <ul>
        <li><a href="#codeDiffusionsetup_systemcode"><code>Diffusion::setup_system</code></a>
        <li><a href="#codeDiffusionassemble_systemcode"><code>Diffusion::assemble_system</code></a>
        <li><a href="#codeDiffusionget_sourcecode"><code>Diffusion::get_source</code></a>
        <li><a href="#codeDiffusionevaluate_diffusioncode"><code>Diffusion::evaluate_diffusion</code></a>
        <li><a href="#codeDiffusionid_minus_tau_J_inversecode"><code>Diffusion::id_minus_tau_J_inverse</code></a>
        <li><a href="#codeDiffusionoutput_resultscode"><code>Diffusion::output_results</code></a>
        <li><a href="#codeDiffusionexplicit_methodcode"><code>Diffusion::explicit_method</code></a>
        <li><a href="#codeDiffusionimplicit_methodcode"><code>Diffusion::implicit_method</code></a>
        <li><a href="#codeDiffusionembedded_explicit_methodcode"><code>Diffusion::embedded_explicit_method</code></a>
        <li><a href="#codeDiffusionruncode"><code>Diffusion::run</code></a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-52/doc/intro.dox

 <br> 

<i>This program was contributed by Bruno Turcksin and Damien Lebrun-Grandie.</i>

 @note  为了运行这个程序，deal.II必须被配置为使用UMFPACK稀疏直接解算器。请参考<a
href="../../readme.html#umfpack">ReadMe</a>中的说明如何做到这一点。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序展示了如何使用Runge-Kutta方法来解决一个随时间变化的问题。它解决了首先在步骤26中讨论的热方程的一个小变化，但是由于这个程序的目的只是演示使用更高级的方法与deal.II的时间步进算法对接，所以只解决了一个均匀细化网格上的简单问题。




<a name="Problemstatement"></a><h3>Problem statement</h3>


在这个例子中，我们求解中子输运方程的单组时间依赖性扩散近似（关于时间依赖性多组扩散，见步骤28）。这是一个关于中子如何在高散射介质中移动的模型，因此它是时间依赖性扩散方程的一个变体--它只是步骤26中讨论的热方程的一个不同名称，加上一些额外的条款。我们假设介质是不可逆的，因此，中子通量满足以下方程。

@f{eqnarray*}
\frac{1}{v}\frac{\partial \phi(x,t)}{\partial t} = \nabla \cdot D(x) \nabla \phi(x,t)


- \Sigma_a(x) \phi(x,t) + S(x,t)


@f}

通过适当的边界条件增强。这里， $v$ 是中子的速度（为简单起见，我们假设它等于1，这可以通过简单地缩放时间变量来实现）， $D$ 是扩散系数， $\Sigma_a$ 是吸收截面， $S$ 是一个源。因为我们只对时间依赖性感兴趣，我们假设 $D$ 和 $\Sigma_a$ 是常数。

由于这个程序只打算演示如何使用先进的时间步进算法，我们将只寻找相对简单问题的解。具体来说，我们要在一个正方形域 $[0,b]\times[0,b]$ 上寻找一个解，其形式为

@f{eqnarray*}
\phi(x,t) = A\sin(\omega t)(bx-x^2).


@f}

通过使用二次有限元，我们可以在任何特定时间精确地表示这个函数，所有的误差都是由于时间离散化造成的。我们这样做是因为这样就很容易观察到我们将要考虑的各种时间步进方案的收敛顺序，而不需要将空间和时间误差分开。

我们施加以下边界条件：对 $x=0$ 和 $x=b$ 施加同质的迪里希特条件，对 $y=0$ 和 $y=b$ 施加同质的纽曼条件。我们选择源项，以便相应的解决方案实际上是上述的形式。

@f{eqnarray*}
S=A\left(\frac{1}{v}\omega \cos(\omega t)(bx -x^2) + \sin(\omega t)
\left(\Sigma_a (bx-x^2)+2D\right) \right).


@f}

因为解是时间上的正弦，我们知道精确解满足 $\phi\left(x,\frac{\pi}{\omega}\right) = 0$  。因此，时间 $t=\frac{\pi}{\omega}$ 的误差只是数值解的规范，即 $\|e(\cdot,t=\frac{\pi}{\omega})\|_{L_2} = \|\phi_h(\cdot,t=\frac{\pi}{\omega})\|_{L_2}$ ，而且特别容易评估。在代码中，我们评估 $l_2$ 的节点值的规范，而不是相关空间函数的 $L_2$ 规范，因为前者的计算更简单；然而，在均匀网格上，这两者只是由一个常数相关，因此我们可以用其中一个观察时间收敛顺序。




<a name="RungeKuttamethods"></a><h3>Runge-Kutta methods</h3>


在deal.II中实现的Runge-Kutta方法假定要解决的方程可以写成。

@f{eqnarray*}
\frac{dy}{dt} = g(t,y).


@f}

另一方面，当使用有限元时，离散化的时间导数总是导致左手边存在一个质量矩阵。这可以很容易地看出，如果上述方程中的解向量 $y(t)$ 实际上是节点系数的向量 $U(t)$ ，其形式为变量

@f{eqnarray*}
  u_h(x,t) = \sum_j U_j(t) \varphi_j(x)


@f}

用空间形状函数 $\varphi_j(x)$ ，然后乘以一个形式的方程

@f{eqnarray*}
  \frac{\partial u(x,t)}{\partial t} = q(t,u(x,t))


@f}

通过测试函数，对 $\Omega$ 进行积分，代入 $u\rightarrow u_h$ 并将测试函数限制在上面的 $\varphi_i(x)$ ，那么这个空间离散方程的形式为

@f{eqnarray*}
M\frac{dU}{dt} = f(t,U),


@f}

其中 $M$ 是质量矩阵， $f(t,U)$ 是 $q(t,u(x,t))$ 的空间离散版本（其中 $q$ 通常是出现空间导数的地方，但鉴于我们只考虑时间导数，这一点目前并不太关心）。换句话说，这种形式符合上面的一般方案，如果我们写成

@f{eqnarray*}
\frac{dy}{dt} = g(t,y) = M^{-1}f(t,y).


@f}



Runk-Kutta方法是一种时间步进方案，通过特定的一步法对 $y(t_n)\approx
y_{n}$ 进行近似。它们通常被写成以下形式

@f{eqnarray*}
y_{n+1} = y_n + \sum_{i=1}^s b_i k_i


@f}

其中对于上面的右手边的形式

@f{eqnarray*}
k_i = h M^{-1} f\left(t_n+c_ih,y_n+\sum_{j=1}^sa_{ij}k_j\right).


@f}

这里 $a_{ij}$ ,  $b_i$ , 和 $c_i$ 是已知的系数，确定你要使用的特定Runge-Kutta方案， $h=t_{n+1}-t_n$ 是使用的时间步长。Runge-Kutta类的不同时间步长方法在级数 $s$ 和系数 $a_{ij}$ 、 $b_i$ 和 $c_i$ 上有所不同，但由于可以查找这些系数的表格值，所以很容易实施。这些表格通常被称为Butcher tableaus）。

在编写本教程时，deal.II中实现的方法可分为三类。<ol>  <li>  显式Runge-Kutta；为了使一个方法成为显式，必须在上述定义 $k_i$ 的公式中， $k_i$ 不出现在右侧。换句话说，这些方法必须满足  $a_{ii}=0, i=1,\ldots,s$  。   <li>  嵌入式（或自适应）Runge-Kutta；我们将在下面讨论其特性。   <li>  隐式Runge-Kutta；这类方法需要解决可能是非线性系统的上述阶段 $k_i$ ，即它们至少有 $a_{ii}\neq 0$  个阶段 $i=1,\ldots,s$  。   </ol>  许多众所周知的时间步进方案，人们通常不会将其与Runge或Kutta的名字联系起来，事实上，它们也可以用这些类别来表达。它们往往代表这些系列的最低阶成员。




<a name="ExplicitRungeKuttamethods"></a><h4>Explicit Runge-Kutta methods</h4>


这些方法，只需要一个函数来评估 $M^{-1}f(t,y)$ ，但不需要（作为隐式方法）来解决涉及 $f(t,y)$ 的 $y$ 的方程。与所有显式时间步长方法一样，当选择的时间步长过大时，它们会变得不稳定。

这一类众所周知的方法包括正向欧拉、三阶Runge-Kutta和四阶Runge-Kutta（通常缩写为RK4）。




<a name="EmbeddedRungeKuttamethods"></a><h4>Embedded Runge-Kutta methods</h4>


这些方法同时使用低阶和高阶方法来估计误差，并决定是否需要缩短时间步长或可以增加。术语 "嵌入 "是指低阶方法不需要对函数 $M^{-1}f(\cdot,\cdot)$ 进行额外的评估，而是重复使用那些必须为高阶方法计算的数据。换句话说，它基本上是免费的，而我们得到的误差估计是使用高阶方法的副产品。

这类方法包括Heun-Euler、Bogacki-Shampine、Dormand-Prince（在Matlab中为ode45，通常缩写为RK45，表示这里使用的低阶和高阶方法分别为4阶和5阶Runge-Kutta方法），Fehlberg和Cash-Karp。

在撰写本文时，只有嵌入式的显式方法得到了实现。




<a name="ImplicitRungeKuttamethods"></a><h4>Implicit Runge-Kutta methods</h4>


隐式方法要求在每个（子）时间步中解决 $\alpha y = f(t,y)$ 形式的 $y$ 的（可能是非线性）系统。在内部，这是用牛顿式方法完成的，因此，它们要求用户提供能够评估 $M^{-1}f(t,y)$ 和 $\left(I-\tau M^{-1} \frac{\partial f}{\partial y}\right)^{-1}$ 或等价的 $\left(M - \tau \frac{\partial f}{\partial y}\right)^{-1} M$ 的函数。

这个算子的特殊形式来自于这样一个事实，即每一个牛顿步骤都需要解一个形式的方程

@f{align*}
  \left(M - \tau \frac{\partial f}{\partial y}\right) \Delta y
  = -M h(t,y)


@f}

对于一些（给定的） $h(t,y)$  。无论时间步长如何，隐式方法总是稳定的，但过大的时间步长当然会影响到解的<i>accuracy</i>，即使数值解仍然稳定且有界。

这类方法包括后退欧拉法、隐式中点法、Crank-Nicolson法和两阶段SDIRK法（"单对角隐式Runge-Kutta "的简称，这个术语是用来表示定义时间步进方法的对角线元素 $a_{ii}$ 都是相等的；这个特性使得牛顿矩阵 $I-\tau M^{-1}\frac{\partial f}{\partial y}$ 可以在各阶段之间重复使用，因为 $\tau$ 每次都是相同的）。




<a name="Spatiallydiscreteformulation"></a><h3>Spatially discrete formulation</h3>


通过扩大我们的模型问题的解决方案，一如既往地使用形状函数 $\psi_j$ 并写出

@f{eqnarray*}
\phi_h(x,t) = \sum_j U_j(t) \psi_j(x),


@f}

我们立即得到扩散方程的空间离散化版本为

@f{eqnarray*}
  M \frac{dU(t)}{dt}
  = -{\cal D} U(t) - {\cal A} U(t) + {\cal S}(t)


@f}

其中

@f{eqnarray*}
  M_{ij}  &=& (\psi_i,\psi_j), \\
  {\cal D}_{ij}  &=& (D\nabla\psi_i,\nabla\psi_j)_\Omega, \\
  {\cal A}_{ij}  &=& (\Sigma_a\psi_i,\psi_j)_\Omega, \\
  {\cal S}_{i}(t)  &=& (\psi_i,S(x,t))_\Omega.


@f}

参见第24步和第26步以了解我们如何到达这里。由于当前问题所选择的边界条件，边界项是没有必要的。为了使用Runge-Kutta方法，我们将其改写如下。

@f{eqnarray*}
f(y) = -{\cal D}y - {\cal A}y + {\cal S}.


@f}

在代码中，我们将需要能够评估这个函数 $f(U)$ 以及它的导数。

@f{eqnarray*}
\frac{\partial f}{\partial y} = -{\cal D} - {\cal A}.


@f}






<a name="Notesonthetestcase"></a><h3>Notes on the testcase</h3>


为了简化问题，域是二维的，网格是均匀细化的（不需要调整网格，因为我们使用的是二次有限元，而且精确解是二次的）。从二维域到三维域并不是很有挑战性。然而，如果你打算解决更复杂的问题，必须对网格进行调整（例如在步骤26中），那么就必须记住以下问题。

<ol>  <li>  在改变网格时，你需要将解投影到新的网格上。当然，从每个时间步长的开始到结束，所使用的网格应该是相同的，这个问题的出现是因为Runge-Kutta方法在每个时间步长内使用了多次方程求值。   <li>  每次改变网格时，你都需要更新质量矩阵和它的逆值。   </ol>  这些步骤的技术可以通过查看步骤26轻易获得。


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
 * The first task as usual is to include the functionality of these well-known
 * deal.II library files and some C++ header files.
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_out.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparse_direct.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * #include <cmath>
 * #include <map>
 * 
 * @endcode
 * 
 * This is the only include file that is new: It includes all the Runge-Kutta
 * methods.
 * 
 * @code
 * #include <deal.II/base/time_stepping.h>
 * 
 * 
 * @endcode
 * 
 * The next step is like in all previous tutorial programs: We put everything
 * into a namespace of its own and then import the deal.II classes and functions
 * into it.
 * 
 * @code
 * namespace Step52
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeDiffusioncodeclass"></a> 
 * <h3>The <code>Diffusion</code> class</h3>
 * 

 * 
 * The next piece is the declaration of the main class. Most of the
 * functions in this class are not new and have been explained in previous
 * tutorials. The only interesting functions are
 * <code>evaluate_diffusion()</code> and
 * <code>id_minus_tau_J_inverse()</code>. <code>evaluate_diffusion()</code>
 * evaluates the diffusion equation, $M^{-1}(f(t,y))$, at a given time and a
 * given $y$. <code>id_minus_tau_J_inverse()</code> evaluates $\left(I-\tau
 * M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ or equivalently
 * $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$ at a given
 * time, for a given $\tau$ and $y$. This function is needed when an
 * implicit method is used.
 * 
 * @code
 *   class Diffusion
 *   {
 *   public:
 *     Diffusion();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 * 
 *     void assemble_system();
 * 
 *     double get_source(const double time, const Point<2> &point) const;
 * 
 *     Vector<double> evaluate_diffusion(const double          time,
 *                                       const Vector<double> &y) const;
 * 
 *     Vector<double> id_minus_tau_J_inverse(const double          time,
 *                                           const double          tau,
 *                                           const Vector<double> &y);
 * 
 *     void output_results(const double                     time,
 *                         const unsigned int               time_step,
 *                         TimeStepping::runge_kutta_method method) const;
 * 
 * @endcode
 * 
 * The next three functions are the drivers for the explicit methods, the
 * implicit methods, and the embedded explicit methods respectively. The
 * driver function for embedded explicit methods returns the number of
 * steps executed given that it only takes the number of time steps passed
 * as an argument as a hint, but internally computed the optimal time step
 * itself.
 * 
 * @code
 *     void explicit_method(const TimeStepping::runge_kutta_method method,
 *                          const unsigned int                     n_time_steps,
 *                          const double                           initial_time,
 *                          const double                           final_time);
 * 
 *     void implicit_method(const TimeStepping::runge_kutta_method method,
 *                          const unsigned int                     n_time_steps,
 *                          const double                           initial_time,
 *                          const double                           final_time);
 * 
 *     unsigned int
 *     embedded_explicit_method(const TimeStepping::runge_kutta_method method,
 *                              const unsigned int n_time_steps,
 *                              const double       initial_time,
 *                              const double       final_time);
 * 
 * 
 *     const unsigned int fe_degree;
 * 
 *     const double diffusion_coefficient;
 *     const double absorption_cross_section;
 * 
 *     Triangulation<2> triangulation;
 * 
 *     const FE_Q<2> fe;
 * 
 *     DoFHandler<2> dof_handler;
 * 
 *     AffineConstraints<double> constraint_matrix;
 * 
 *     SparsityPattern sparsity_pattern;
 * 
 *     SparseMatrix<double> system_matrix;
 *     SparseMatrix<double> mass_matrix;
 *     SparseMatrix<double> mass_minus_tau_Jacobian;
 * 
 *     SparseDirectUMFPACK inverse_mass_matrix;
 * 
 *     Vector<double> solution;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * We choose quadratic finite elements and we initialize the parameters.
 * 
 * @code
 *   Diffusion::Diffusion()
 *     : fe_degree(2)
 *     , diffusion_coefficient(1. / 30.)
 *     , absorption_cross_section(1.)
 *     , fe(fe_degree)
 *     , dof_handler(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionsetup_systemcode"></a> 
 * <h4><code>Diffusion::setup_system</code></h4>
 * Now, we create the constraint matrix and the sparsity pattern. Then, we
 * initialize the matrices and the solution vector.
 * 
 * @code
 *   void Diffusion::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              1,
 *                                              Functions::ZeroFunction<2>(),
 *                                              constraint_matrix);
 *     constraint_matrix.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraint_matrix);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *     mass_matrix.reinit(sparsity_pattern);
 *     mass_minus_tau_Jacobian.reinit(sparsity_pattern);
 *     solution.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionassemble_systemcode"></a> 
 * <h4><code>Diffusion::assemble_system</code></h4>
 * In this function, we compute $-\int D \nabla b_i \cdot \nabla b_j
 * d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$ and the mass
 * matrix $\int b_i b_j d\boldsymbol{r}$. The mass matrix is then
 * inverted using a direct solver; the <code>inverse_mass_matrix</code>
 * variable will then store the inverse of the mass matrix so that
 * $M^{-1}$ can be applied to a vector using the <code>vmult()</code>
 * function of that object. (Internally, UMFPACK does not really store
 * the inverse of the matrix, but its LU factors; applying the inverse
 * matrix is then equivalent to doing one forward and one backward solves
 * with these two factors, which has the same complexity as applying an
 * explicit inverse of the matrix).
 * 
 * @code
 *   void Diffusion::assemble_system()
 *   {
 *     system_matrix = 0.;
 *     mass_matrix   = 0.;
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1);
 * 
 *     FEValues<2> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients | update_JxW_values);
 * 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix      = 0.;
 *         cell_mass_matrix = 0.;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               {
 *                 cell_matrix(i, j) +=
 *                   ((-diffusion_coefficient *                // (-D
 *                       fe_values.shape_grad(i, q_point) *    //  * grad phi_i
 *                       fe_values.shape_grad(j, q_point)      //  * grad phi_j
 *                     - absorption_cross_section *            //  -Sigma
 *                         fe_values.shape_value(i, q_point) * //  * phi_i
 *                         fe_values.shape_value(j, q_point))  //  * phi_j)
 *                    * fe_values.JxW(q_point));               // * dx
 *                 cell_mass_matrix(i, j) += fe_values.shape_value(i, q_point) *
 *                                           fe_values.shape_value(j, q_point) *
 *                                           fe_values.JxW(q_point);
 *               }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraint_matrix.distribute_local_to_global(cell_matrix,
 *                                                      local_dof_indices,
 *                                                      system_matrix);
 *         constraint_matrix.distribute_local_to_global(cell_mass_matrix,
 *                                                      local_dof_indices,
 *                                                      mass_matrix);
 *       }
 * 
 *     inverse_mass_matrix.initialize(mass_matrix);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionget_sourcecode"></a> 
 * <h4><code>Diffusion::get_source</code></h4>
 *   

 * 
 * In this function, the source term of the equation for a given time and a
 * given point is computed.
 * 
 * @code
 *   double Diffusion::get_source(const double time, const Point<2> &point) const
 *   {
 *     const double intensity = 10.;
 *     const double frequency = numbers::PI / 10.;
 *     const double b         = 5.;
 *     const double x         = point(0);
 * 
 *     return intensity *
 *            (frequency * std::cos(frequency * time) * (b * x - x * x) +
 *             std::sin(frequency * time) *
 *               (absorption_cross_section * (b * x - x * x) +
 *                2. * diffusion_coefficient));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionevaluate_diffusioncode"></a> 
 * <h4><code>Diffusion::evaluate_diffusion</code></h4>
 *   

 * 
 * Next, we evaluate the weak form of the diffusion equation at a given time
 * $t$ and for a given vector $y$. In other words, as outlined in the
 * introduction, we evaluate $M^{-1}(-{\cal D}y - {\cal A}y + {\cal
 * S})$. For this, we have to apply the matrix $-{\cal D} - {\cal A}$
 * (previously computed and stored in the variable
 * <code>system_matrix</code>) to $y$ and then add the source term which we
 * integrate as we usually do. (Integrating up the solution could be done
 * using VectorTools::create_right_hand_side() if you wanted to save a few
 * lines of code, or wanted to take advantage of doing the integration in
 * parallel.) The result is then multiplied by $M^{-1}$.
 * 
 * @code
 *   Vector<double> Diffusion::evaluate_diffusion(const double          time,
 *                                                const Vector<double> &y) const
 *   {
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     tmp = 0.;
 *     system_matrix.vmult(tmp, y);
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1);
 * 
 *     FEValues<2> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_quadrature_points |
 *                             update_JxW_values);
 * 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double> cell_source(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_source = 0.;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           {
 *             const double source =
 *               get_source(time, fe_values.quadrature_point(q_point));
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_source(i) += fe_values.shape_value(i, q_point) * // phi_i(x)
 *                                 source *                            // * S(x)
 *                                 fe_values.JxW(q_point);             // * dx
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraint_matrix.distribute_local_to_global(cell_source,
 *                                                      local_dof_indices,
 *                                                      tmp);
 *       }
 * 
 *     Vector<double> value(dof_handler.n_dofs());
 *     inverse_mass_matrix.vmult(value, tmp);
 * 
 *     return value;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionid_minus_tau_J_inversecode"></a> 
 * <h4><code>Diffusion::id_minus_tau_J_inverse</code></h4>
 *   

 * 
 * We compute $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$. This
 * is done in several steps:
 * - compute $M-\tau \frac{\partial f}{\partial y}$
 * - invert the matrix to get $\left(M-\tau \frac{\partial f}
 * {\partial y}\right)^{-1}$
 * - compute $tmp=My$
 * - compute $z=\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} tmp =
 * \left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} My$
 * - return z.
 * 
 * @code
 *   Vector<double> Diffusion::id_minus_tau_J_inverse(const double /*time*/,
 *                                                    const double          tau,
 *                                                    const Vector<double> &y)
 *   {
 *     SparseDirectUMFPACK inverse_mass_minus_tau_Jacobian;
 * 
 *     mass_minus_tau_Jacobian.copy_from(mass_matrix);
 *     mass_minus_tau_Jacobian.add(-tau, system_matrix);
 * 
 *     inverse_mass_minus_tau_Jacobian.initialize(mass_minus_tau_Jacobian);
 * 
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     mass_matrix.vmult(tmp, y);
 * 
 *     Vector<double> result(y);
 *     inverse_mass_minus_tau_Jacobian.vmult(result, tmp);
 * 
 *     return result;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionoutput_resultscode"></a> 
 * <h4><code>Diffusion::output_results</code></h4>
 *   

 * 
 * The following function then outputs the solution in vtu files indexed by
 * the number of the time step and the name of the time stepping method. Of
 * course, the (exact) result should really be the same for all time
 * stepping method, but the output here at least allows us to compare them.
 * 
 * @code
 *   void Diffusion::output_results(const double                     time,
 *                                  const unsigned int               time_step,
 *                                  TimeStepping::runge_kutta_method method) const
 *   {
 *     std::string method_name;
 * 
 *     switch (method)
 *       {
 *         case TimeStepping::FORWARD_EULER:
 *           {
 *             method_name = "forward_euler";
 *             break;
 *           }
 *         case TimeStepping::RK_THIRD_ORDER:
 *           {
 *             method_name = "rk3";
 *             break;
 *           }
 *         case TimeStepping::RK_CLASSIC_FOURTH_ORDER:
 *           {
 *             method_name = "rk4";
 *             break;
 *           }
 *         case TimeStepping::BACKWARD_EULER:
 *           {
 *             method_name = "backward_euler";
 *             break;
 *           }
 *         case TimeStepping::IMPLICIT_MIDPOINT:
 *           {
 *             method_name = "implicit_midpoint";
 *             break;
 *           }
 *         case TimeStepping::SDIRK_TWO_STAGES:
 *           {
 *             method_name = "sdirk";
 *             break;
 *           }
 *         case TimeStepping::HEUN_EULER:
 *           {
 *             method_name = "heun_euler";
 *             break;
 *           }
 *         case TimeStepping::BOGACKI_SHAMPINE:
 *           {
 *             method_name = "bocacki_shampine";
 *             break;
 *           }
 *         case TimeStepping::DOPRI:
 *           {
 *             method_name = "dopri";
 *             break;
 *           }
 *         case TimeStepping::FEHLBERG:
 *           {
 *             method_name = "fehlberg";
 *             break;
 *           }
 *         case TimeStepping::CASH_KARP:
 *           {
 *             method_name = "cash_karp";
 *             break;
 *           }
 *         default:
 *           {
 *             break;
 *           }
 *       }
 * 
 *     DataOut<2> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 * 
 *     data_out.build_patches();
 * 
 *     data_out.set_flags(DataOutBase::VtkFlags(time, time_step));
 * 
 *     const std::string filename = "solution_" + method_name + "-" +
 *                                  Utilities::int_to_string(time_step, 3) +
 *                                  ".vtu";
 *     std::ofstream output(filename);
 *     data_out.write_vtu(output);
 * 
 *     static std::vector<std::pair<double, std::string>> times_and_names;
 * 
 *     static std::string method_name_prev = "";
 *     static std::string pvd_filename;
 *     if (method_name_prev != method_name)
 *       {
 *         times_and_names.clear();
 *         method_name_prev = method_name;
 *         pvd_filename     = "solution_" + method_name + ".pvd";
 *       }
 *     times_and_names.emplace_back(time, filename);
 *     std::ofstream pvd_output(pvd_filename);
 *     DataOutBase::write_pvd_record(pvd_output, times_and_names);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionexplicit_methodcode"></a> 
 * <h4><code>Diffusion::explicit_method</code></h4>
 *   

 * 
 * This function is the driver for all the explicit methods. At the
 * top it initializes the time stepping and the solution (by setting
 * it to zero and then ensuring that boundary value and hanging node
 * constraints are respected; of course, with the mesh we use here,
 * hanging node constraints are not in fact an issue). It then calls
 * <code>evolve_one_time_step</code> which performs one time step.
 * Time is stored and incremented through a DiscreteTime object.
 *   

 * 
 * For explicit methods, <code>evolve_one_time_step</code> needs to
 * evaluate $M^{-1}(f(t,y))$, i.e, it needs
 * <code>evaluate_diffusion</code>. Because
 * <code>evaluate_diffusion</code> is a member function, it needs to
 * be bound to <code>this</code>. After each evolution step, we
 * again apply the correct boundary values and hanging node
 * constraints.
 *   

 * 
 * Finally, the solution is output
 * every 10 time steps.
 * 
 * @code
 *   void Diffusion::explicit_method(const TimeStepping::runge_kutta_method method,
 *                                   const unsigned int n_time_steps,
 *                                   const double       initial_time,
 *                                   const double       final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::ExplicitRungeKutta<Vector<double>> explicit_runge_kutta(
 *       method);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         explicit_runge_kutta.evolve_one_time_step(
 *           [this](const double time, const Vector<double> &y) {
 *             return this->evaluate_diffusion(time, y);
 *           },
 *           time.get_current_time(),
 *           time.get_next_step_size(),
 *           solution);
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionimplicit_methodcode"></a> 
 * <h4><code>Diffusion::implicit_method</code></h4>
 * This function is equivalent to <code>explicit_method</code> but for
 * implicit methods. When using implicit methods, we need to evaluate
 * $M^{-1}(f(t,y))$ and $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial
 * y}\right)^{-1}$ for which we use the two member functions previously
 * introduced.
 * 
 * @code
 *   void Diffusion::implicit_method(const TimeStepping::runge_kutta_method method,
 *                                   const unsigned int n_time_steps,
 *                                   const double       initial_time,
 *                                   const double       final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::ImplicitRungeKutta<Vector<double>> implicit_runge_kutta(
 *       method);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         implicit_runge_kutta.evolve_one_time_step(
 *           [this](const double time, const Vector<double> &y) {
 *             return this->evaluate_diffusion(time, y);
 *           },
 *           [this](const double time, const double tau, const Vector<double> &y) {
 *             return this->id_minus_tau_J_inverse(time, tau, y);
 *           },
 *           time.get_current_time(),
 *           time.get_next_step_size(),
 *           solution);
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionembedded_explicit_methodcode"></a> 
 * <h4><code>Diffusion::embedded_explicit_method</code></h4>
 * This function is the driver for the embedded explicit methods. It requires
 * more parameters:
 * - coarsen_param: factor multiplying the current time step when the error
 * is below the threshold.
 * - refine_param: factor multiplying the current time step when the error
 * is above the threshold.
 * - min_delta: smallest time step acceptable.
 * - max_delta: largest time step acceptable.
 * - refine_tol: threshold above which the time step is refined.
 * - coarsen_tol: threshold below which the time step is coarsen.
 *   

 * 
 * Embedded methods use a guessed time step. If the error using this time step
 * is too large, the time step will be reduced. If the error is below the
 * threshold, a larger time step will be tried for the next time step.
 * <code>delta_t_guess</code> is the guessed time step produced by the
 * embedded method. In summary, time step size is potentially modified in
 * three ways:
 * - Reducing or increasing time step size within
 * TimeStepping::EmbeddedExplicitRungeKutta::evolve_one_time_step().
 * - Using the calculated <code>delta_t_guess</code>.
 * - Automatically adjusting the step size of the last time step to ensure
 * simulation ends precisely at <code>final_time</code>. This adjustment
 * is handled inside the DiscreteTime instance.
 * 
 * @code
 *   unsigned int Diffusion::embedded_explicit_method(
 *     const TimeStepping::runge_kutta_method method,
 *     const unsigned int                     n_time_steps,
 *     const double                           initial_time,
 *     const double                           final_time)
 *   {
 *     const double time_step =
 *       (final_time - initial_time) / static_cast<double>(n_time_steps);
 *     const double coarsen_param = 1.2;
 *     const double refine_param  = 0.8;
 *     const double min_delta     = 1e-8;
 *     const double max_delta     = 10 * time_step;
 *     const double refine_tol    = 1e-1;
 *     const double coarsen_tol   = 1e-5;
 * 
 *     solution = 0.;
 *     constraint_matrix.distribute(solution);
 * 
 *     TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>>
 *       embedded_explicit_runge_kutta(method,
 *                                     coarsen_param,
 *                                     refine_param,
 *                                     min_delta,
 *                                     max_delta,
 *                                     refine_tol,
 *                                     coarsen_tol);
 *     output_results(initial_time, 0, method);
 *     DiscreteTime time(initial_time, final_time, time_step);
 *     while (time.is_at_end() == false)
 *       {
 *         const double new_time =
 *           embedded_explicit_runge_kutta.evolve_one_time_step(
 *             [this](const double time, const Vector<double> &y) {
 *               return this->evaluate_diffusion(time, y);
 *             },
 *             time.get_current_time(),
 *             time.get_next_step_size(),
 *             solution);
 *         time.set_next_step_size(new_time - time.get_current_time());
 *         time.advance_time();
 * 
 *         constraint_matrix.distribute(solution);
 * 
 *         if (time.get_step_number() % 10 == 0)
 *           output_results(time.get_current_time(),
 *                          time.get_step_number(),
 *                          method);
 * 
 *         time.set_desired_next_step_size(
 *           embedded_explicit_runge_kutta.get_status().delta_t_guess);
 *       }
 * 
 *     return time.get_step_number();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionruncode"></a> 
 * <h4><code>Diffusion::run</code></h4>
 *   

 * 
 * The following is the main function of the program. At the top, we create
 * the grid (a [0,5]x[0,5] square) and refine it four times to get a mesh
 * that has 16 by 16 cells, for a total of 256.  We then set the boundary
 * indicator to 1 for those parts of the boundary where $x=0$ and $x=5$.
 * 
 * @code
 *   void Diffusion::run()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0., 5.);
 *     triangulation.refine_global(4);
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             if ((face->center()[0] == 0.) || (face->center()[0] == 5.))
 *               face->set_boundary_id(1);
 *             else
 *               face->set_boundary_id(0);
 *           }
 * 
 * @endcode
 * 
 * Next, we set up the linear systems and fill them with content so that
 * they can be used throughout the time stepping process:
 * 
 * @code
 *     setup_system();
 * 
 *     assemble_system();
 * 
 * @endcode
 * 
 * Finally, we solve the diffusion problem using several of the
 * Runge-Kutta methods implemented in namespace TimeStepping, each time
 * outputting the error at the end time. (As explained in the
 * introduction, since the exact solution is zero at the final time, the
 * error equals the numerical solution and can be computed by just taking
 * the $l_2$ norm of the solution vector.)
 * 
 * @code
 *     unsigned int       n_steps      = 0;
 *     const unsigned int n_time_steps = 200;
 *     const double       initial_time = 0.;
 *     const double       final_time   = 10.;
 * 
 *     std::cout << "Explicit methods:" << std::endl;
 *     explicit_method(TimeStepping::FORWARD_EULER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Forward Euler:            error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     explicit_method(TimeStepping::RK_THIRD_ORDER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Third order Runge-Kutta:  error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     explicit_method(TimeStepping::RK_CLASSIC_FOURTH_ORDER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Fourth order Runge-Kutta: error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << std::endl;
 * 
 * 
 *     std::cout << "Implicit methods:" << std::endl;
 *     implicit_method(TimeStepping::BACKWARD_EULER,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Backward Euler:           error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::IMPLICIT_MIDPOINT,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Implicit Midpoint:        error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::CRANK_NICOLSON,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   Crank-Nicolson:           error=" << solution.l2_norm()
 *               << std::endl;
 * 
 *     implicit_method(TimeStepping::SDIRK_TWO_STAGES,
 *                     n_time_steps,
 *                     initial_time,
 *                     final_time);
 *     std::cout << "   SDIRK:                    error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << std::endl;
 * 
 * 
 *     std::cout << "Embedded explicit methods:" << std::endl;
 *     n_steps = embedded_explicit_method(TimeStepping::HEUN_EULER,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Heun-Euler:               error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::BOGACKI_SHAMPINE,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Bogacki-Shampine:         error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::DOPRI,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Dopri:                    error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::FEHLBERG,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Fehlberg:                 error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::CASH_KARP,
 *                                        n_time_steps,
 *                                        initial_time,
 *                                        final_time);
 *     std::cout << "   Cash-Karp:                error=" << solution.l2_norm()
 *               << std::endl;
 *     std::cout << "                   steps performed=" << n_steps << std::endl;
 *   }
 * } // namespace Step52
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * The following <code>main</code> function is similar to previous examples
 * and need not be commented on.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step52::Diffusion diffusion;
 *       diffusion.run();
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
 *     };
 * 
 *   return 0;
 * }
 * @endcode
examples/step-52/doc/results.dox



<a name="Results"></a><h1>Results</h1>


这个程序的重点不在于显示特定的结果，而在于显示它是如何做到的。这一点我们已经通过讨论上面的代码证明过了。因此，该程序的输出相对较少，只包括控制台输出和用于可视化的VTU格式的解决方案。

控制台输出既包含错误，也包含对某些方法所执行的步骤数量。

@code
Explicit methods:
   Forward Euler:            error=1.00883
   Third order Runge-Kutta:  error=0.000227982
   Fourth order Runge-Kutta: error=1.90541e-06


Implicit methods:
   Backward Euler:           error=1.03428
   Implicit Midpoint:        error=0.00862702
   Crank-Nicolson:           error=0.00862675
   SDIRK:                    error=0.0042349


Embedded explicit methods:
   Heun-Euler:               error=0.0073012
                   steps performed=284
   Bogacki-Shampine:         error=0.000408407
                   steps performed=181
   Dopri:                    error=0.000836695
                   steps performed=120
   Fehlberg:                 error=0.00248922
                   steps performed=106
   Cash-Karp:                error=0.0787735
                   steps performed=106
@endcode



正如预期的那样，高阶方法给出了（更）准确的解决方案。我们还看到，（相当不准确的）Heun-Euler方法增加了时间步数，以满足公差要求。另一方面，其他嵌入式方法使用的时间步数比规定的要少得多。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-52.cc"
*/
