/**
@page step_72 The step-72 tutorial program
This tutorial depends on step-71, step-15.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Motivation">Motivation</a>
        <li><a href="#ComputingtheJacobianfromtheresidual"> Computing the Jacobian from the residual </a>
        <li><a href="#ComputingtheJacobianandtheresidualfromtheenergyfunctional"> Computing the Jacobian and the residual from the energy functional </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMinimalSurfaceProblemParameterscodeclass">The <code>MinimalSurfaceProblemParameters</code> class</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a>
        <li><a href="#Boundarycondition">Boundary condition</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#MinimalSurfaceProblemMinimalSurfaceProblem">MinimalSurfaceProblem::MinimalSurfaceProblem</a>
        <li><a href="#MinimalSurfaceProblemsetup_system">MinimalSurfaceProblem::setup_system</a>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
      <ul>
        <li><a href="#Manualassembly">Manual assembly</a>
        <li><a href="#Assemblyviadifferentiationoftheresidualvector">Assembly via differentiation of the residual vector</a>
        <li><a href="#Assemblyviadifferentiationoftheenergyfunctional">Assembly via differentiation of the energy functional</a>
      </ul>
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
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-72/doc/intro.dox

 <br> 

<i>This program was contributed by Jean-Paul Pelteret and Wolfgang Bangerth.


Wolfgang Bangerth's work is partially supported by National Science
Foundation grants OCI-1148116, OAC-1835673, DMS-1821210, and EAR-1925595;
and by the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-1550901 and The University of California-Davis.
</i>




<a name="Introduction"></a><h1>Introduction</h1>


<a name="Motivation"></a><h3>Motivation</h3>


这个程序解决的问题与步骤15相同，即求解[最小表面方程](https://en.wikipedia.org/wiki/Minimal_surface) @f{align*}
    F(u) \dealcoloneq -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.
  @f}



我们在那里发现的问题（见<a href="step_15#extensions">Possibilities for extensions</a>部分）是，当想要使用牛顿迭代时，我们需要计算方程残差对解的导数 $u$ （这里，因为右手边是零，残差只是左手边）。对于我们这里的方程来说，这很麻烦，但并非不可能 -- 但我们很容易想象出更复杂的方程，仅仅正确实现残差本身就是一个挑战，更不用说为计算雅各布矩阵所需的导数而这样做。我们将在这个程序中解决这个问题。使用在步骤-71中详细讨论的自动微分技术，我们将想出一个办法，我们只需要实现残差，就可以免费得到雅各布矩阵。

事实上，我们甚至可以更进一步。虽然在第15步中，我们只是把方程作为一个给定值，但最小表面方程实际上是最小化一个能量的产物。具体来说，最小曲面方程是对应于最小化能量的欧拉-拉格朗日方程@f[
    E(u) = \int_\Omega \Psi \left( u \right)
  @f]

其中*能量密度*由@f[
    \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}}.
  @f]给出。

这等于说，我们寻求找到能量函数变化的静止点@f[
    \min\limits_{u} E \left( u \right)
      \quad \rightarrow \quad
      \delta E \left( u, \varphi \right) \dealcoloneq
      \left(\varphi, F(u)\right) = 0
      \qquad
      \forall \varphi,
  @f] 。

因为这是边界值问题的平衡解所在。

那么关键的一点是，也许，我们甚至不需要实现残差，但实现更简单的能量密度 $\Psi(u)$ 可能实际上已经足够了。

那么我们的目标是这样的。当使用牛顿迭代时，我们需要反复解决线性偏微分方程@f{align*}
    F'(u^{n},\delta u^{n}) &=- F(u^{n})
  @f}。

这样我们就可以计算出更新@f{align*}
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}。

与牛顿步骤的解 $\delta u^{n}$ 。正如步骤15所讨论的，我们可以用手计算导数 $F'(u,\delta u)$ ，得到@f[
  F'(u,\delta u)
  =


  - \nabla \cdot \left( \frac{1}{\left(1+|\nabla u|^{2}\right)^{\frac{1}{2}}}\nabla
  \delta u \right) +
  \nabla \cdot \left( \frac{\nabla u \cdot
  \nabla \delta u}{\left(1+|\nabla u|^{2}\right)^{\frac{3}{2}}} \nabla u
  \right).
  @f]。



那么，这里就是这个计划的内容。它是关于可以帮助我们计算 $F'(u,\delta u)$ 的技术，而不必明确地实现它，要么提供 $F(u)$ 的实现，要么提供 $E(u)$  的实现。更确切地说，我们将实现三种不同的方法，并在运行时间方面进行比较，但同时--也许更重要的是--实现这些方法需要多少人力。

- 第15步中使用的方法，形成雅各布矩阵。

- 从残差 $F(u)$ 的实现中计算雅各布矩阵，使用自动微分法。

- 从能量函数 $E(u)$ 的实现中计算残差和雅各布矩阵，也使用自动微分法。

对于这些方法中的第一个，与步骤15相比，没有任何概念上的变化。




<a name="ComputingtheJacobianfromtheresidual"></a><h3> Computing the Jacobian from the residual </h3>


对于第二种方法，让我们概述一下我们将如何利用自动微分来计算残差向量的线性化。为此，让我们暂时改变一下符号，用 $F(U)$ 表示的不是微分方程的残差，而实际上是*残差向量*，即*离散残差。我们这样做是因为当我们在给定的网格上对问题进行离散时，这就是我们*实际*做的事情。我们解决 $F(U)=0$ 问题，其中 $U$ 是未知数的矢量。

更准确地说，残差的 $i$ th分量由以下公式给出

@f[
  F(U)_i \dealcoloneq
  \int\limits_{\Omega}\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,


@f]

其中 $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$  。鉴于此，单元格 $K$ 的贡献是

@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,


@f]

它的一阶泰勒展开为

@f[
  F(U + \delta U)_i^K
  \approx F(U)_i^K
  + \sum_{j}^{n_{\textrm{dofs}}} \left[ \frac{\partial F(U)_i^K}{\partial
  U_j} \delta U_j \right],


@f]

因此我们可以计算出 $K$ 单元格对雅各布矩阵 $J$ 的贡献为 $J(U)_{ij}^K = \frac{\partial F(U)_i^K}{\partial U_j}$  。这里重要的一点是，在单元格 $K$ 上，我们可以表示为

@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+\left|
  \sum_{j'}^{n_\textrm{dofs}} U_{j'} \nabla \varphi_{j'}\right|^{2}}}
  \left(\sum_{j''}^{n_\textrm{dofs}} U_{j''} \nabla \varphi_{j''}\right)\right] \, dV.


@f]

为了清楚起见，我们用 $j'$ 和 $j''$ 作为计数索引，以明确它们彼此之间以及与上述 $j$ 的区别。因为在这个公式中， $F(U)$ 只取决于系数 $U_j$ ，我们可以通过自动微分 $F(U)_i^K$ 来计算导数 $J(U)_{ij}^K$ 作为一个矩阵。通过我们一直使用的相同论证，很明显 $F(U)^K$ 实际上并不依赖于*所有*未知数 $U_j$ ，而只是依赖于 $j$ 是住在单元格 $K$ 的形状函数的那些未知数。] ，因此在实践中，我们将 $F(U)^K$ 和 $J(U)^K$ 限制为矢量和矩阵中对应于*本地*DoF指数的部分，然后从本地单元 $K$ 分布到全球对象。

使用所有这些实现，然后的方法将是在程序中实现 $F(U)^K$ ，并让自动微分机械从中计算导数 $J(U)^K$ 。




<a name="ComputingtheJacobianandtheresidualfromtheenergyfunctional"></a><h3> Computing the Jacobian and the residual from the energy functional </h3>


对于装配过程的最终实现，我们将比残差高一个层次：我们的整个线性系统将直接由支配这个边界值问题的物理学的能量函数决定。我们可以利用这样一个事实：我们可以直接从局部贡献中计算出域中的总能量，即。

@f[
  E \left( U \right) \dealcoloneq \int\limits_{\Omega} \Psi \left( u
  \right) \, dV .


@f]

在离散设置中，这意味着在每个有限元上我们有

@f[
   E \left( U \right)^K
    \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
    \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
    \mathbf{x}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times W_{q}}_{\text{JxW(q)}} .


@f]

如果我们实现细胞能量，它取决于场解，我们可以计算它的第一个（离散）变化

@f[
  F(U)^K_i
    = \frac{\partial E(U)^K}{\partial U_i}


@f]

此后，它的第二个（离散）变化

@f[
  J(U)^K_{ij}
    = \frac{\partial^{2}  E(U)^K}{\partial U_i \partial U_j}.


@f]

因此，从单元格对总能量函数的贡献来看，只要我们能够提供局部能量的实现，我们就可以期望为我们生成近似的残差和正切贡献  $E(U)^K$  。同样，由于本教程中使用的自动微分变量的设计，在实践中，这些对残差向量和正切矩阵贡献的近似值实际上是精确到机器精度的。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The majority of this tutorial is an exact replica of step-15. So, in the
 * interest of brevity and maintaining a focus on the changes implemented here,
 * we will only document what's new and simply indicate which sections of
 * code are a repetition of what has come before.
 * 

 * 
 * 

 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * There are a few new header files that have been included in this tutorial.
 * The first is the one that provides the declaration of the ParameterAcceptor
 * class.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/parameter_acceptor.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/utilities.h>
 * 
 * @endcode
 * 
 * This is the second, which is an all-inclusive header that will allow us
 * to incorporate the automatic differentiation (AD) functionality within this
 * code.
 * 
 * @code
 * #include <deal.II/differentiation/ad.h>
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
 * #include <deal.II/fe/fe_values_extractors.h>
 * #include <deal.II/fe/fe_q.h>
 * 
 * @endcode
 * 
 * And the next three provide some multi-threading capability using the generic
 * MeshWorker::mesh_loop() framework.
 * 
 * @code
 * #include <deal.II/meshworker/copy_data.h>
 * #include <deal.II/meshworker/mesh_loop.h>
 * #include <deal.II/meshworker/scratch_data.h>
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
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * @endcode
 * 
 * We then open a namespace for this program and import everything from the
 * dealii namespace into it, as in previous programs:
 * 
 * @code
 * namespace Step72
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemParameterscodeclass"></a> 
 * <h3>The <code>MinimalSurfaceProblemParameters</code> class</h3>
 * 

 * 
 * In this tutorial we will implement three different approaches for
 * assembling the linear system. One mirrors the hand implementation
 * originally provided in step-15, while the other two use the Sacado
 * automatic differentiation library that is provided as a part of the
 * Trilinos framework.
 *   

 * 
 * To facilitate switching between the three implementations, we have
 * this really basic parameters class that has only two options that are
 * configurable.
 * 
 * @code
 *   class MinimalSurfaceProblemParameters : public ParameterAcceptor
 *   {
 *   public:
 *     MinimalSurfaceProblemParameters();
 * 
 * @endcode
 * 
 * Selection for the formulation and corresponding AD framework to be used:
 * -  formulation = 0 : Unassisted implementation (full hand linearization).
 * -  formulation = 1 : Automated linearization of the finite element
 * residual.
 * -  formulation = 2 : Automated computation of finite element
 * residual and linearization using a
 * variational formulation.
 * 
 * @code
 *     unsigned int formulation = 0;
 * 
 * @endcode
 * 
 * The maximum acceptable tolerance for the linear system residual.
 * We will see that the assembly time becomes appreciable once we use
 * the AD framework, so we have increased the tolerance selected in
 * step-15 by one order of magnitude. This way, the computations do
 * not take too long to complete.
 * 
 * @code
 *     double tolerance = 1e-2;
 *   };
 * 
 * 
 *   MinimalSurfaceProblemParameters::MinimalSurfaceProblemParameters()
 *     : ParameterAcceptor("Minimal Surface Problem/")
 *   {
 *     add_parameter(
 *       "Formulation", formulation, "", this->prm, Patterns::Integer(0, 2));
 *     add_parameter("Tolerance", tolerance, "", this->prm, Patterns::Double(0.0));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class template</h3>
 * 

 * 
 * The class template is essentially the same as in step-15.
 * The only functional changes to the class are that:
 * - the run() function now takes in two arguments: one to choose which
 * assembly approach is to be adopted, and one for the tolerance for
 * the permissible final residual is, and
 * - there are now three different assembly functions that implement the
 * three methods of assembling the linear system. We'll provide details
 * on these later on.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class MinimalSurfaceProblem
 *   {
 *   public:
 *     MinimalSurfaceProblem();
 * 
 *     void run(const int formulation, const double tolerance);
 * 
 *   private:
 *     void   setup_system(const bool initial_step);
 *     void   assemble_system_unassisted();
 *     void   assemble_system_with_residual_linearization();
 *     void   assemble_system_using_energy_functional();
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
 *     QGauss<dim>     quadrature_formula;
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
 * There are no changes to the boundary conditions applied to the problem.
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
 * There have been no changes made to the class constructor.
 * 
 * @code
 *   template <int dim>
 *   MinimalSurfaceProblem<dim>::MinimalSurfaceProblem()
 *     : dof_handler(triangulation)
 *     , fe(2)
 *     , quadrature_formula(fe.degree + 1)
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
 * There have been no changes made to the function that sets up the class
 * data structures, namely the DoFHandler, the hanging node constraints
 * applied to the problem, and the linear system.
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
 * <a name="Assemblingthelinearsystem"></a> 
 * <h4>Assembling the linear system</h4>
 * 

 * 
 * 
 * <a name="Manualassembly"></a> 
 * <h5>Manual assembly</h5>
 * 

 * 
 * The assembly functions are the interesting contributions to this tutorial.
 * The assemble_system_unassisted() method implements exactly the same
 * assembly function as is detailed in step-15, but in this instance we
 * use the MeshWorker::mesh_loop() function to multithread the assembly
 * process. The reason for doing this is quite simple: When using
 * automatic differentiation, we know that there is to be some additional
 * computational overhead incurred. In order to mitigate this performance
 * loss, we'd like to take advantage of as many (easily available)
 * computational resources as possible. The MeshWorker::mesh_loop() concept
 * makes this a relatively straightforward task. At the same time, for the
 * purposes of fair comparison, we need to do the same to the implementation
 * that uses no assistance when computing the residual or its linearization.
 * (The MeshWorker::mesh_loop() function is first discussed in step-12 and
 * step-16, if you'd like to read up on it.)
 *   

 * 
 * The steps required to implement the multithreading are the same between the
 * three functions, so we'll use the assemble_system_unassisted() function
 * as an opportunity to focus on the multithreading itself.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::assemble_system_unassisted()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 * @endcode
 * 
 * The MeshWorker::mesh_loop() expects that we provide two exemplar data
 * structures. The first, `ScratchData`, is to store all large data that
 * is to be reused between threads. The `CopyData` will hold the
 * contributions to the linear system that come from each cell. These
 * independent matrix-vector pairs must be accumulated into the
 * global linear system sequentially. Since we don't need anything
 * on top of what the MeshWorker::ScratchData and MeshWorker::CopyData
 * classes already provide, we use these exact class definitions for
 * our problem. Note that we only require a single instance of a local
 * matrix, local right-hand side vector, and cell degree of freedom index
 * vector -- the MeshWorker::CopyData therefore has `1` for all three
 * of its template arguments.
 * 
 * @code
 *     using ScratchData = MeshWorker::ScratchData<dim>;
 *     using CopyData    = MeshWorker::CopyData<1, 1, 1>;
 * 
 * @endcode
 * 
 * We also need to know what type of iterator we'll be working with
 * during assembly. For simplicity, we just ask the compiler to work
 * this out for us using the decltype() specifier, knowing that we'll
 * be iterating over active cells owned by the @p dof_handler.
 * 
 * @code
 *     using CellIteratorType = decltype(dof_handler.begin_active());
 * 
 * @endcode
 * 
 * Here we initialize the exemplar data structures. Since we know that
 * we need to compute the shape function gradients, weighted Jacobian,
 * and the position of the quadrate points in real space, we pass these
 * flags into the class constructor.
 * 
 * @code
 *     const ScratchData sample_scratch_data(fe,
 *                                           quadrature_formula,
 *                                           update_gradients |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 *     const CopyData    sample_copy_data(dofs_per_cell);
 * 
 * @endcode
 * 
 * Now we define a lambda function that will perform the assembly on
 * a single cell. The three arguments are those that will be expected by
 * MeshWorker::mesh_loop(), due to the arguments that we'll pass to that
 * final call. We also capture the @p this pointer, which means that we'll
 * have access to "this" (i.e., the current `MinimalSurfaceProblem<dim>`)
 * class instance, and its private member data (since the lambda function is
 * defined within a MinimalSurfaceProblem<dim> method).
 *     

 * 
 * At the top of the function, we initialize the data structures
 * that are dependent on the cell for which the work is being
 * performed. Observe that the reinitialization call actually
 * returns an instance to an FEValues object that is initialized
 * and stored within (and, therefore, reused by) the
 * `scratch_data` object.
 *     

 * 
 * Similarly, we get aliases to the local matrix, local RHS
 * vector, and local cell DoF indices from the `copy_data`
 * instance that MeshWorker::mesh_loop() provides. We then
 * initialize the cell DoF indices, knowing that the local matrix
 * and vector are already correctly sized.
 * 
 * @code
 *     const auto cell_worker = [this](const CellIteratorType &cell,
 *                                     ScratchData &           scratch_data,
 *                                     CopyData &              copy_data) {
 *       const auto &fe_values = scratch_data.reinit(cell);
 * 
 *       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
 *       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
 *       std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 *       cell->get_dof_indices(local_dof_indices);
 * 
 * @endcode
 * 
 * For Newton's method, we require the gradient of the solution at the
 * point about which the problem is being linearized.
 *       

 * 
 * Once we have that, we can perform assembly for this cell in
 * the usual way.  One minor difference to step-15 is that we've
 * used the (rather convenient) range-based loops to iterate
 * over all quadrature points and degrees-of-freedom.
 * 
 * @code
 *       std::vector<Tensor<1, dim>> old_solution_gradients(
 *         fe_values.n_quadrature_points);
 *       fe_values.get_function_gradients(current_solution,
 *                                        old_solution_gradients);
 * 
 *       for (const unsigned int q : fe_values.quadrature_point_indices())
 *         {
 *           const double coeff =
 *             1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
 *                                     old_solution_gradients[q]);
 * 
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               for (const unsigned int j : fe_values.dof_indices())
 *                 cell_matrix(i, j) +=
 *                   (((fe_values.shape_grad(i, q)      // ((\nabla \phi_i
 *                      * coeff                         //   * a_n
 *                      * fe_values.shape_grad(j, q))   //   * \nabla \phi_j)
 *                     -                                //  -
 *                     (fe_values.shape_grad(i, q)      //  (\nabla \phi_i
 *                      * coeff * coeff * coeff         //   * a_n^3
 *                      * (fe_values.shape_grad(j, q)   //   * (\nabla \phi_j
 *                         * old_solution_gradients[q]) //      * \nabla u_n)
 *                      * old_solution_gradients[q]))   //   * \nabla u_n)))
 *                    * fe_values.JxW(q));              // * dx
 * 
 *               cell_rhs(i) -= (fe_values.shape_grad(i, q)  // \nabla \phi_i
 *                               * coeff                     // * a_n
 *                               * old_solution_gradients[q] // * u_n
 *                               * fe_values.JxW(q));        // * dx
 *             }
 *         }
 *     };
 * 
 * @endcode
 * 
 * The second lambda function that MeshWorker::mesh_loop() requires is
 * one that performs the task of accumulating the local contributions
 * in the global linear system. That is precisely what this one does,
 * and the details of the implementation have been seen before. The
 * primary point to recognize is that the local contributions are stored
 * in the `copy_data` instance that is passed into this function. This
 * `copy_data` has been filled with data during @a some call to the
 * `cell_worker`.
 * 
 * @code
 *     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
 *       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
 *       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
 *       const std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 * 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               cell_matrix(i, j));
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     };
 * 
 * @endcode
 * 
 * We have all of the required functions definitions in place, so
 * now we call the MeshWorker::mesh_loop() to perform the actual
 * assembly.  We pass a flag as the last parameter which states
 * that we only want to perform the assembly on the
 * cells. Internally, MeshWorker::mesh_loop() then distributes the
 * available work to different threads, making efficient use of
 * the multiple cores almost all of today's processors have to
 * offer.
 * 
 * @code
 *     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
 *                           cell_worker,
 *                           copier,
 *                           sample_scratch_data,
 *                           sample_copy_data,
 *                           MeshWorker::assemble_own_cells);
 * 
 * @endcode
 * 
 * And finally, as is done in step-15, we remove hanging nodes from the
 * system and apply zero boundary values to the linear system that defines
 * the Newton updates $\delta u^n$.
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
 * @endcode
 * 
 * 
 * <a name="Assemblyviadifferentiationoftheresidualvector"></a> 
 * <h5>Assembly via differentiation of the residual vector</h5>
 * 

 * 
 * As outlined in the introduction, what we need to do for this
 * second approach is implement the local contributions $F(U)^K$
 * from cell $K$ to the residual vector, and then let the
 * AD machinery deal with how to compute the
 * derivatives $J(U)_{ij}^K=\frac{\partial F(U)^K_i}{\partial U_j}$
 * from it.
 *   

 * 
 * For the following, recall that
 * @f[
 * F(U)_i^K \dealcoloneq
 * \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
 * u|^{2}}} \nabla u \right] \, dV ,
 * @f]
 * where $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$.
 *   

 * 
 * Let us see how this is implemented in practice:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::assemble_system_with_residual_linearization()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *     using ScratchData      = MeshWorker::ScratchData<dim>;
 *     using CopyData         = MeshWorker::CopyData<1, 1, 1>;
 *     using CellIteratorType = decltype(dof_handler.begin_active());
 * 
 *     const ScratchData sample_scratch_data(fe,
 *                                           quadrature_formula,
 *                                           update_gradients |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 *     const CopyData    sample_copy_data(dofs_per_cell);
 * 
 * @endcode
 * 
 * We'll define up front the AD data structures that we'll be using,
 * utilizing the techniques shown in step-71.
 * In this case, we choose the helper class that will automatically compute
 * the linearization of the finite element residual using Sacado forward
 * automatic differentiation types. These number types can be used to
 * compute first derivatives only. This is exactly what we want, because we
 * know that we'll only be linearizing the residual, which means that we
 * only need to compute first-order derivatives. The return values from the
 * calculations are to be of type `double`.
 *     

 * 
 * We also need an extractor to retrieve some data related to the field
 * solution to the problem.
 * 
 * @code
 *     using ADHelper = Differentiation::AD::ResidualLinearization<
 *       Differentiation::AD::NumberTypes::sacado_dfad,
 *       double>;
 *     using ADNumberType = typename ADHelper::ad_type;
 * 
 *     const FEValuesExtractors::Scalar u_fe(0);
 * 
 * @endcode
 * 
 * With this, let us define the lambda function that will be used
 * to compute the cell contributions to the Jacobian matrix and
 * the right hand side:
 * 
 * @code
 *     const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
 *                                            ScratchData &           scratch_data,
 *                                            CopyData &              copy_data) {
 *       const auto &       fe_values     = scratch_data.reinit(cell);
 *       const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell();
 * 
 *       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
 *       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
 *       std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 *       cell->get_dof_indices(local_dof_indices);
 * 
 * @endcode
 * 
 * We'll now create and initialize an instance of the AD helper class.
 * To do this, we need to specify how many independent variables and
 * dependent variables there are. The independent variables will be the
 * number of local degrees of freedom that our solution vector has,
 * i.e., the number $j$ in the per-element representation of the
 * discretized solution vector
 * $u (\mathbf{x})|_K = \sum\limits_{j} U^K_i \varphi_j(\mathbf{x})$
 * that indicates how many solution coefficients are associated with
 * each finite element. In deal.II, this equals
 * FiniteElement::dofs_per_cell. The number of dependent variables will be
 * the number of entries in the local residual vector that we will be
 * forming. In this particular problem (like many others that employ the
 * [standard Galerkin
 * method](https://en.wikipedia.org/wiki/Galerkin_method)) the number of
 * local solution coefficients matches the number of local residual
 * equations.
 * 
 * @code
 *       const unsigned int n_independent_variables = local_dof_indices.size();
 *       const unsigned int n_dependent_variables   = dofs_per_cell;
 *       ADHelper ad_helper(n_independent_variables, n_dependent_variables);
 * 
 * @endcode
 * 
 * Next we inform the helper of the values of the solution, i.e., the
 * actual values for $U_j$ about which we
 * wish to linearize. As this is done on each element individually, we
 * have to extract the solution coefficients from the global solution
 * vector. In other words, we define all of those coefficients $U_j$
 * where $j$ is a local degree of freedom as the independent variables
 * that enter the computation of the vector $F(U)^{K}$ (the dependent
 * function).
 *       

 * 
 * Then we get the complete set of degree of freedom values as
 * represented by auto-differentiable numbers. The operations
 * performed with these variables are tracked by the AD library
 * from this point until the object goes out of scope. So it is
 * <em>precisely these variables</em> with respect to which we will
 * compute derivatives of the residual entries.
 * 
 * @code
 *       ad_helper.register_dof_values(current_solution, local_dof_indices);
 * 
 *       const std::vector<ADNumberType> &dof_values_ad =
 *         ad_helper.get_sensitive_dof_values();
 * 
 * @endcode
 * 
 * Then we do some problem specific tasks, the first being to
 * compute all values, (spatial) gradients, and the like based on
 * "sensitive" AD degree of freedom values. In this instance we are
 * retrieving the solution gradients at each quadrature point. Observe
 * that the solution gradients are now sensitive
 * to the values of the degrees of freedom as they use the @p ADNumberType
 * as the scalar type and the @p dof_values_ad vector provides the local
 * DoF values.
 * 
 * @code
 *       std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
 *         fe_values.n_quadrature_points);
 *       fe_values[u_fe].get_function_gradients_from_local_dof_values(
 *         dof_values_ad, old_solution_gradients);
 * 
 * @endcode
 * 
 * The next variable that we declare will store the cell residual vector
 * contributions. This is rather self-explanatory, save for one
 * <b>very important</b> detail:
 * Note that each entry in the vector is hand-initialized with a value
 * of zero. This is a <em>highly recommended</em> practice, as some AD
 * libraries appear not to safely initialize the internal data
 * structures of these number types. Not doing so could lead to some
 * very hard to understand or detect bugs (appreciate that the author
 * of this program mentions this out of, generally bad, experience). So
 * out of an abundance of caution it's worthwhile zeroing the initial
 * value explicitly. After that, apart from a sign change the residual
 * assembly looks much the same as we saw for the cell RHS vector before:
 * We loop over all quadrature points, ensure that the coefficient now
 * encodes its dependence on the (sensitive) finite element DoF values by
 * using the correct `ADNumberType`, and finally we assemble the
 * components of the residual vector. For complete clarity, the finite
 * element shape functions (and their gradients, etc.) as well as the
 * "JxW" values remain scalar
 * valued, but the @p coeff and the  @p old_solution_gradients at each
 * quadrature point are computed in terms of the independent
 * variables.
 * 
 * @code
 *       std::vector<ADNumberType> residual_ad(n_dependent_variables,
 *                                             ADNumberType(0.0));
 *       for (const unsigned int q : fe_values.quadrature_point_indices())
 *         {
 *           const ADNumberType coeff =
 *             1.0 / std::sqrt(1.0 + old_solution_gradients[q] *
 *                                     old_solution_gradients[q]);
 * 
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               residual_ad[i] += (fe_values.shape_grad(i, q)   // \nabla \phi_i
 *                                  * coeff                      // * a_n
 *                                  * old_solution_gradients[q]) // * u_n
 *                                 * fe_values.JxW(q);           // * dx
 *             }
 *         }
 * 
 * @endcode
 * 
 * Once we have the full cell residual vector computed, we can register
 * it with the helper class.
 *       

 * 
 * Thereafter, we compute the residual values (basically,
 * extracting the real values from what we already computed) and
 * their Jacobian (the linearization of each residual component
 * with respect to all cell DoFs) at the evaluation point. For
 * the purposes of assembly into the global linear system, we
 * have to respect the sign difference between the residual and
 * the RHS contribution: For Newton's method, the right hand
 * side vector needs to be equal to the *negative* residual
 * vector.
 * 
 * @code
 *       ad_helper.register_residual_vector(residual_ad);
 * 
 *       ad_helper.compute_residual(cell_rhs);
 *       cell_rhs *= -1.0;
 * 
 *       ad_helper.compute_linearization(cell_matrix);
 *     };
 * 
 * @endcode
 * 
 * The remainder of the function equals what we had previously:
 * 
 * @code
 *     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
 *       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
 *       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
 *       const std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 * 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               cell_matrix(i, j));
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     };
 * 
 *     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
 *                           cell_worker,
 *                           copier,
 *                           sample_scratch_data,
 *                           sample_copy_data,
 *                           MeshWorker::assemble_own_cells);
 * 
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
 * @endcode
 * 
 * 
 * <a name="Assemblyviadifferentiationoftheenergyfunctional"></a> 
 * <h5>Assembly via differentiation of the energy functional</h5>
 * 

 * 
 * In this third approach, we compute residual and Jacobian as first
 * and second derivatives of the local energy functional
 * @f[
 * E\left( U \right)^K
 * \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
 * \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
 * \mathbf{X}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times
 * W_{q}}_{\text{JxW(q)}}
 * @f]
 * with the energy density given by
 * @f[
 * \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}} .
 * @f]
 *   

 * 
 * Let us again see how this is done:
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::assemble_system_using_energy_functional()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *     using ScratchData      = MeshWorker::ScratchData<dim>;
 *     using CopyData         = MeshWorker::CopyData<1, 1, 1>;
 *     using CellIteratorType = decltype(dof_handler.begin_active());
 * 
 *     const ScratchData sample_scratch_data(fe,
 *                                           quadrature_formula,
 *                                           update_gradients |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 *     const CopyData    sample_copy_data(dofs_per_cell);
 * 
 * @endcode
 * 
 * In this implementation of the assembly process, we choose the helper
 * class that will automatically compute both the residual and its
 * linearization from the cell contribution to an energy functional using
 * nested Sacado forward automatic differentiation types.
 * The selected number types can be used to compute both first and
 * second derivatives. We require this, as the residual defined as the
 * sensitivity of the potential energy with respect to the DoF values (i.e.
 * its gradient). We'll then need to linearize the residual, implying that
 * second derivatives of the potential energy must be computed. You might
 * want to compare this with the definition of `ADHelper` used int
 * previous function, where we used
 * `Differentiation::AD::ResidualLinearization<Differentiation::AD::NumberTypes::sacado_dfad,double>`.
 * 
 * @code
 *     using ADHelper = Differentiation::AD::EnergyFunctional<
 *       Differentiation::AD::NumberTypes::sacado_dfad_dfad,
 *       double>;
 *     using ADNumberType = typename ADHelper::ad_type;
 * 
 *     const FEValuesExtractors::Scalar u_fe(0);
 * 
 * @endcode
 * 
 * Let us then again define the lambda function that does the integration on
 * a cell.
 *     

 * 
 * To initialize an instance of the helper class, we now only require
 * that the number of independent variables (that is, the number
 * of degrees of freedom associated with the element solution vector)
 * are known up front. This is because the second-derivative matrix that
 * results from an energy functional is necessarily square (and also,
 * incidentally, symmetric).
 * 
 * @code
 *     const auto cell_worker = [&u_fe, this](const CellIteratorType &cell,
 *                                            ScratchData &           scratch_data,
 *                                            CopyData &              copy_data) {
 *       const auto &fe_values = scratch_data.reinit(cell);
 * 
 *       FullMatrix<double> &                  cell_matrix = copy_data.matrices[0];
 *       Vector<double> &                      cell_rhs    = copy_data.vectors[0];
 *       std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 *       cell->get_dof_indices(local_dof_indices);
 * 
 *       const unsigned int n_independent_variables = local_dof_indices.size();
 *       ADHelper           ad_helper(n_independent_variables);
 * 
 * @endcode
 * 
 * Once more, we register all cell DoFs values with the helper, followed
 * by extracting the "sensitive" variant of these values that are to be
 * used in subsequent operations that must be differentiated -- one of
 * those being the calculation of the solution gradients.
 * 
 * @code
 *       ad_helper.register_dof_values(current_solution, local_dof_indices);
 * 
 *       const std::vector<ADNumberType> &dof_values_ad =
 *         ad_helper.get_sensitive_dof_values();
 * 
 *       std::vector<Tensor<1, dim, ADNumberType>> old_solution_gradients(
 *         fe_values.n_quadrature_points);
 *       fe_values[u_fe].get_function_gradients_from_local_dof_values(
 *         dof_values_ad, old_solution_gradients);
 * 
 * @endcode
 * 
 * We next create a variable that stores the cell total energy.
 * Once more we emphasize that we explicitly zero-initialize this value,
 * thereby ensuring the integrity of the data for this starting value.
 *       

 * 
 * The aim for our approach is then to compute the cell total
 * energy, which is the sum of the internal (due to right hand
 * side functions, typically linear in $U$) and external
 * energies. In this particular case, we have no external
 * energies (e.g., from source terms or Neumann boundary
 * conditions), so we'll focus on the internal energy part.
 *       

 * 
 * In fact, computing $E(U)^K$ is almost trivial, requiring only
 * the following lines:
 * 
 * @code
 *       ADNumberType energy_ad = ADNumberType(0.0);
 *       for (const unsigned int q : fe_values.quadrature_point_indices())
 *         {
 *           const ADNumberType psi = std::sqrt(1.0 + old_solution_gradients[q] *
 *                                                      old_solution_gradients[q]);
 * 
 *           energy_ad += psi * fe_values.JxW(q);
 *         }
 * 
 * @endcode
 * 
 * After we've computed the total energy on this cell, we'll
 * register it with the helper.  Based on that, we may now
 * compute the desired quantities, namely the residual values
 * and their Jacobian at the evaluation point. As before, the
 * Newton right hand side needs to be the negative of the
 * residual:
 * 
 * @code
 *       ad_helper.register_energy_functional(energy_ad);
 * 
 *       ad_helper.compute_residual(cell_rhs);
 *       cell_rhs *= -1.0;
 * 
 *       ad_helper.compute_linearization(cell_matrix);
 *     };
 * 
 * @endcode
 * 
 * As in the previous two functions, the remainder of the function is as
 * before:
 * 
 * @code
 *     const auto copier = [dofs_per_cell, this](const CopyData &copy_data) {
 *       const FullMatrix<double> &cell_matrix = copy_data.matrices[0];
 *       const Vector<double> &    cell_rhs    = copy_data.vectors[0];
 *       const std::vector<types::global_dof_index> &local_dof_indices =
 *         copy_data.local_dof_indices[0];
 * 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               cell_matrix(i, j));
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     };
 * 
 *     MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
 *                           cell_worker,
 *                           copier,
 *                           sample_scratch_data,
 *                           sample_copy_data,
 *                           MeshWorker::assemble_own_cells);
 * 
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
 * @endcode
 * 
 * 
 * <a name="MinimalSurfaceProblemsolve"></a> 
 * <h4>MinimalSurfaceProblem::solve</h4>
 * 

 * 
 * The solve function is the same as is used in step-15.
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
 * Nothing has changed since step-15 with respect to the mesh refinement
 * procedure and transfer of the solution between adapted meshes.
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
 *     triangulation.prepare_coarsening_and_refinement();
 *     SolutionTransfer<dim> solution_transfer(dof_handler);
 *     solution_transfer.prepare_for_coarsening_and_refinement(current_solution);
 *     triangulation.execute_coarsening_and_refinement();
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     Vector<double> tmp(dof_handler.n_dofs());
 *     solution_transfer.interpolate(current_solution, tmp);
 *     current_solution = tmp;
 * 
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 *     set_boundary_values();
 * 
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
 * The choice of boundary conditions remains identical to step-15...
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
 * ... as does the function used to compute the residual during the
 * solution iteration procedure. One could replace this by
 * differentiation of the energy functional if one really wanted,
 * but for simplicity we here simply copy what we already had in
 * step-15.
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
 *         fe_values.get_function_gradients(evaluation_point, gradients);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double coeff =
 *               1.0 / std::sqrt(1.0 + gradients[q] * gradients[q]);
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
 *     hanging_node_constraints.condense(residual);
 * 
 *     for (types::global_dof_index i :
 *          DoFTools::extract_boundary_dofs(dof_handler))
 *       residual(i) = 0;
 * 
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
 * The choice of step length (or, under-relaxation factor) for the nonlinear
 * iterations procedure remains fixed at the value chosen and discussed in
 * step-15.
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
 * In the run function, most remains the same as was first implemented
 * in step-15. The only observable changes are that we can now choose (via
 * the parameter file) what the final acceptable tolerance for the system
 * residual is, and that we can choose which method of assembly we wish to
 * utilize. To make the second choice clear, we output to the console some
 * message which indicates the selection. Since we're interested in comparing
 * the time taken to assemble for each of the three methods, we've also
 * added a timer that keeps a track of how much time is spent during assembly.
 * We also track the time taken to solve the linear system, so that we can
 * contrast those numbers to the part of the code which would normally take
 * the longest time to execute.
 * 
 * @code
 *   template <int dim>
 *   void MinimalSurfaceProblem<dim>::run(const int    formulation,
 *                                        const double tolerance)
 *   {
 *     std::cout << "******** Assembly approach ********" << std::endl;
 *     const std::array<std::string, 3> method_descriptions = {
 *       {"Unassisted implementation (full hand linearization).",
 *        "Automated linearization of the finite element residual.",
 *        "Automated computation of finite element residual and linearization using a variational formulation."}};
 *     AssertIndexRange(formulation, method_descriptions.size());
 *     std::cout << method_descriptions[formulation] << std::endl << std::endl;
 * 
 * 
 *     TimerOutput timer(std::cout, TimerOutput::summary, TimerOutput::wall_times);
 * 
 *     GridGenerator::hyper_ball(triangulation);
 *     triangulation.refine_global(2);
 * 
 *     setup_system(/*first time=*/true);
 *     set_boundary_values();
 * 
 *     double       last_residual_norm = std::numeric_limits<double>::max();
 *     unsigned int refinement_cycle   = 0;
 *     do
 *       {
 *         std::cout << "Mesh refinement step " << refinement_cycle << std::endl;
 * 
 *         if (refinement_cycle != 0)
 *           refine_mesh();
 * 
 *         std::cout << "  Initial residual: " << compute_residual(0) << std::endl;
 * 
 *         for (unsigned int inner_iteration = 0; inner_iteration < 5;
 *              ++inner_iteration)
 *           {
 *             {
 *               TimerOutput::Scope t(timer, "Assemble");
 * 
 *               if (formulation == 0)
 *                 assemble_system_unassisted();
 *               else if (formulation == 1)
 *                 assemble_system_with_residual_linearization();
 *               else if (formulation == 2)
 *                 assemble_system_using_energy_functional();
 *               else
 *                 AssertThrow(false, ExcNotImplemented());
 *             }
 * 
 *             last_residual_norm = system_rhs.l2_norm();
 * 
 *             {
 *               TimerOutput::Scope t(timer, "Solve");
 *               solve();
 *             }
 * 
 * 
 *             std::cout << "  Residual: " << compute_residual(0) << std::endl;
 *           }
 * 
 *         output_results(refinement_cycle);
 * 
 *         ++refinement_cycle;
 *         std::cout << std::endl;
 *       }
 *     while (last_residual_norm > tolerance);
 *   }
 * } // namespace Step72
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h4>The main function</h4>
 * 

 * 
 * Finally the main function. This follows the scheme of most other main
 * functions, with two obvious exceptions:
 * - We call Utilities::MPI::MPI_InitFinalize in order to set up (via a hidden
 * default parameter) the number of threads using the execution of
 * multithreaded tasks.
 * - We also have a few lines dedicates to reading in or initializing the
 * user-defined parameters that will be considered during the execution of the
 * program.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step72;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
 * 
 *       std::string prm_file;
 *       if (argc > 1)
 *         prm_file = argv[1];
 *       else
 *         prm_file = "parameters.prm";
 * 
 *       const MinimalSurfaceProblemParameters parameters;
 *       ParameterAcceptor::initialize(prm_file);
 * 
 *       MinimalSurfaceProblem<2> minimal_surface_problem_2d;
 *       minimal_surface_problem_2d.run(parameters.formulation,
 *                                      parameters.tolerance);
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
examples/step-72/doc/results.dox



<a name="Results"></a><h1>Results</h1>


由于在步骤15中首先分析的问题的物理学没有变化，所以没有什么可报告的。它们之间唯一外显的区别是，在默认情况下，这个程序只运行9个网格细化步骤（相对于第15步，执行11个细化）。这可以从模拟状态中观察到，该状态出现在打印出正在使用的装配方法的标题文本和最终的时间。下面报告的所有时间都是在发布模式下获得的）。

@code
Mesh refinement step 0
  Initial residual: 1.53143
  Residual: 1.08746
  Residual: 0.966748
  Residual: 0.859602
  Residual: 0.766462
  Residual: 0.685475


...


Mesh refinement step 9
  Initial residual: 0.00924594
  Residual: 0.00831928
  Residual: 0.0074859
  Residual: 0.0067363
  Residual: 0.00606197
  Residual: 0.00545529
@endcode



因此，我们感兴趣的是比较三种不同实现方式的装配过程需要多长时间，并把它放到更大的背景中。下面是手部线性化的输出结果（在2012年左右的四核八线程笔记本电脑上计算的结果--但我们真正感兴趣的只是不同实现方式之间的相对时间）。

@code
******** Assembly approach ********
Unassisted implementation (full hand linearization).


...


+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      35.1s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |      1.56s |       4.5% |
| Solve                           |        50 |      30.8s |        88% |
+---------------------------------+-----------+------------+------------+
@endcode

而对于使用萨卡多动态正向AD数字类型，以自动方式将残差线性化的实施。

@code
******** Assembly approach ********
Automated linearization of the finite element residual.


...


+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      40.1s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |       8.8s |        22% |
| Solve                           |        50 |      28.6s |        71% |
+---------------------------------+-----------+------------+------------+
@endcode

最后，对于直接从能量函数（使用嵌套的Sacado动态前向AD数）计算残差和其线性化的实现。

@code
******** Assembly approach ********
Automated computation of finite element residual and linearization using a variational formulation.


...


+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      48.8s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |        50 |      16.7s |        34% |
| Solve                           |        50 |      29.3s |        60% |
+---------------------------------+-----------+------------+------------+
@endcode



很明显，交给自动分化框架执行的工作越多，在装配过程中花费的时间就越多。在所有细化步骤中累积起来，与无辅助装配相比，使用一级自动微分导致在装配阶段花费了 $5.65 \times$ 的计算时间，而直接从能量函数推导时，装配离散线性系统花费了 $10.7 \times$ 的时间。不足为奇的是，解决线性系统的总体时间保持不变。这意味着，随着在有限元水平上进行自动微分的次数的增加，花在求解阶段的时间与装配阶段的时间比例发生了明显的转变。对许多人来说，这可能意味着在生产代码中利用高阶微分（在有限元水平）会导致不可接受的开销，但在原型设计阶段，它可能仍然有用。因此，两者之间的一个很好的折衷办法是有限元残差的自动线性化，它以可衡量的、但也许不是不可接受的成本提供了很多便利。另外，我们可以考虑不在每一步中重新建立牛顿矩阵--这个主题在步骤77中有大量的深入探讨。

当然，在实践中，实际的开销在很大程度上取决于被评估的问题（例如，解决方案中有多少成分，被微分的函数的性质是什么，等等）。因此，这里提出的确切结果应该仅在这个标量问题的背景下进行解释，当涉及到其他问题时，用户肯定需要进行一些初步调查。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


与步骤-71一样，有几个与自动区分有关的项目可以进一步评估。

- 应调查其他AD框架的使用情况，并展望其他实施方式可能提供性能优势。

- 除了本教程中硬编码的数字类型外，还值得对AD数字类型进行评估。关于在有限元水平上采用的两次微分类型，混合微分模式（"RAD-FAD"）原则上应该比这里采用的单一模式（"FAD-FAD"）类型的计算效率更高。RAD-FAD类型没有被默认选择的原因是，在撰写本文时，在Sacado库中，它的实现仍然存在一些错误，导致内存泄漏。   这在 @ref auto_symb_diff 模块中有所记载。

- 也许使用低精度类型（即 "浮动"）作为AD数字的标量类型可以减少装配时的计算费用。使用 "float "作为矩阵和残差的数据类型并不是不合理的，因为牛顿更新只是为了让我们更接近解决方案，而不是实际*到解决方案；因此，考虑使用降低精度的数据类型来计算这些更新，然后将这些更新累积到使用全 "双 "精度的解决方案向量中，是有意义的。

- 在装配过程中可能减少资源的另一个方法是将AD的实现作为一个构成模型。这类似于步骤71中采用的方法，并将自动微分的起点推到了计算链的上一级。这反过来意味着AD库跟踪的操作更少，从而降低了微分的成本（尽管我们会在每个单元的正交点进行微分）。

- 第77步是第15步的另一个变化，解决了问题的一个非常不同的部分：直线搜索以及是否有必要在每次非线性迭代中重新建立牛顿矩阵。鉴于上述结果表明，使用自动微分是有代价的，第77步的技术有可能在一定程度上抵消这些代价。因此，将目前的程序与第77步中的技术结合起来是相当有趣的。对于生产代码来说，这肯定是个好办法。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-72.cc"
*/
