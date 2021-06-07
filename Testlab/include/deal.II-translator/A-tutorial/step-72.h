//include/deal.II-translator/A-tutorial/step-72_0.txt
/**
  @page step_72 The step-72 tutorial program  
* 本教程依赖于  step-71  ,  step-15  。
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
        <li><a href="#Manualassembly">Manual assembly</a> ]<a href="#Manualassembly">Manual assembly</a>
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
*  <br>  
* <i>This program was contributed by Jean-Paul Pelteret and Wolfgang Bangerth.
* 
Wolfgang Bangerth's work is partially supported by National Science
Foundation grants OCI-1148116, OAC-1835673, DMS-1821210, and EAR-1925595;
and by the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-1550901 and The University of California-Davis.
</i>
* 

*<a name="Introduction"></a><h1>Introduction</h1>
* 

*<a name="Motivation"></a><h3>Motivation</h3>
 

* 这个程序解决的问题与 step-15 相同，即解决[最小表面方程](https://en.wikipedia.org/wiki/Minimal_surface) 
* @f{align*}
    F(u) \dealcoloneq
* 
-\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.
  @f}
* 
* 我们在那里发现的问题（见<a href="step_15#extensions">Possibilities for extensions</a>部分）是，当想要使用牛顿迭代时，我们需要计算方程残差对解的导数 $u$ （这里，因为右手边是零，残差只是左手边的）。对于我们这里的方程，这是很麻烦的，但不是不可能的
* 
* 但是我们可以很容易地想象出更复杂的方程，仅仅正确地实现残差本身就是一个挑战，更不用说为计算雅各布矩阵所需的导数而这样做了。我们将在这个程序中解决这个问题。使用在 step-71 中详细讨论过的自动微分技术，我们将想出一个办法，我们只需要实现残差就可以免费得到雅各布矩阵。
*事实上，我们甚至可以更进一步。在 step-15 中，我们只是把方程当作一个给定的东西，而最小表面方程实际上是最小化一个能量的产物。具体来说，最小表面方程是对应于能量最小化的欧拉-拉格朗日方程@f[
    E(u) = \int_\Omega \Psi \left( u \right)
  @f]，其中能量密度*由@f[
    \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}}.
  @f]给出，这就等于说我们要找到能量函数@f[
    \min\limits_{u} E \left( u \right)
      \quad \rightarrow \quad
      \delta E \left( u, \varphi \right) \dealcoloneq
      \left(\varphi, F(u)\right) = 0
      \qquad
      \forall \varphi,
  @f]的变化的静止点，因为这就是边界值问题的平衡解所在。
* 那么关键的一点是，也许，我们甚至不需要实现常数，但实现更简单的能量密度 $\Psi(u)$ 实际上就足够了。
* 那么我们的目标是这样的。当使用牛顿迭代时，我们需要反复解决线性偏微分方程 
* @f{align*}
    F'(u^{n},\delta u^{n}) &=- F(u^{n})
  @f}
* 这样我们就可以计算出更新 
* @f{align*}
    u^{n+1}&=u^{n}+\alpha^n \delta u^{n}
  @f}
*与牛顿步骤的解决方案 $\delta u^{n}$ 。正如在 step-15 中所讨论的，我们可以用手计算导数 $F'(u,\delta u)$ ，并获得@f[
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
  @f]。
* 因此，这就是这个程序的意义所在。它是关于可以帮助我们计算 $F'(u,\delta u)$ 的技术，而不必明确地实现它，或者提供 $F(u)$ 的实现或 $E(u)$ 的实现。更确切地说，我们将实现三种不同的方法，并在运行时间方面对它们进行比较，但也要注意
* 
* - 也许更重要的是
* 
* - 实现它们需要多少人力。
* 
* - 在 step-15 中使用的方法来形成雅各布矩阵。
* 
* - 从残差的实现中计算雅各布矩阵  $F(u)$  ，使用自动微分。
* 
* - 从能量函数 $E(u)$ 的实现中计算残差和雅各布矩阵，同样使用自动微分。
* 对于这些方法中的第一个，与 step-15 相比，没有概念上的变化。
* 

* <a name="ComputingtheJacobianfromtheresidual"></a><h3> Computing the Jacobian from the residual </h3> 。
* 

* 对于第二种方法，让我们概述一下我们将如何利用自动微分来计算残差向量的线性化来处理这个问题。为此，让我们暂时改变一下符号，用 $F(U)$ 表示的不是微分方程的残差，而是事实上的残差向量*。
* 
*--即离散的残差*。我们这样做是因为这就是我们
*当我们在一个给定的网格上对问题进行离散化时，我们实际上是这样做的。我们解决 $F(U)=0$ 问题，其中 $U$ 是未知数的向量。
* 更确切地说，残差的 $i$ 分量由@f[
  F(U)_i \dealcoloneq
  \int\limits_{\Omega}\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,
@f]给出，其中 $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$  。鉴于此，单元格 $K$ 的贡献为@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
  u|^{2}}} \nabla u \right] \, dV ,
@f]，其一阶泰勒展开为@f[
  F(U + \delta U)_i^K
  \approx F(U)_i^K
  + \sum_{j}^{n_{\textrm{dofs}}} \left[ \frac{\partial F(U)_i^K}{\partial
  U_j} \delta U_j \right],
@f]，因此我们可以计算单元格 $K$ 对雅各布矩阵 $J$ 的贡献为 $J(U)_{ij}^K = \frac{\partial F(U)_i^K}{\partial U_j}$  。这里重要的一点是，在单元格 $K$ 上，我们可以表示@f[
  F(U)_i^K \dealcoloneq
  \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+\left|
  \sum_{j'}^{n_\textrm{dofs}} U_{j'} \nabla \varphi_{j'}\right|^{2}}}
  \left(\sum_{j''}^{n_\textrm{dofs}} U_{j''} \nabla \varphi_{j''}\right)\right] \, dV.
@f]为了清楚起见，我们用 $j'$ 和 $j''$ 作为计数索引，以明确它们彼此之间以及与 $j$ 不同。 因为在这个公式中， $F(U)$ 只取决于系数 $U_j$ ，我们可以通过 $F(U)_i^K$ 的自动微分，将导数 $J(U)_{ij}^K$ 计算为一个矩阵。通过我们一直使用的相同论证，很明显 $F(U)^K$ 实际上并不取决于
*所有*未知数 $U_j$ ，而只取决于 $j$ 是住在单元格 $K$ 上的形状函数的那些未知数，因此在实践中，我们将 $F(U)^K$ 和 $J(U)^K$ 限制在矢量和矩阵的那一部分，对应于本地*DoF指数，然后从本地单元格 $K$ 分配到全球对象。
* 使用所有这些实现，然后的方法是在程序中实现 $F(U)^K$ ，并让自动微分机械从中计算导数 $J(U)^K$ 。
* 

*<a name="ComputingtheJacobianandtheresidualfromtheenergyfunctional"></a><h3> Computing the Jacobian and the residual from the energy functional </h3>。
* 

* 对于装配过程的最终实现，我们将比残差高一个层次：我们的整个线性系统将直接由支配这个边界值问题的物理学的能量函数决定。我们可以利用这样一个事实，即我们可以直接从局部贡献中计算域中的总能量，即@f[
  E \left( U \right) \dealcoloneq \int\limits_{\Omega} \Psi \left( u
  \right) \, dV .
@f]在离散设置中，这意味着在每个有限元上我们有@f[
   E \left( U \right)^K
    \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
    \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
    \mathbf{x}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times W_{q}}_{\text{JxW(q)}} .
@f]如果我们实现单元能量，这取决于场解，我们可以计算其第一个（离散）变化@f[
  F(U)^K_i
    = \frac{\partial E(U)^K}{\partial U_i}
@f] 因此，从单元对总能量函数的贡献来看，只要我们能够提供局部能量的实现，我们就可以期待为我们生成近似的残差和切线贡献 $E(U)^K$ 。同样，由于本教程中使用的自动微分变量的设计，在实践中，这些对残差向量和正切矩阵的贡献的近似值实际上是精确到机器精度的。
* 

* <a name="CommProg"></a> <h1> The commented program</h1>。
* 本教程的大部分内容是对 step-15 的精确复制。因此，为了简洁起见，并保持对这里实现的变化的关注，我们将只记录新的内容，并简单地指出哪些部分的代码是对以前内容的重复。
* 

* 
*  
*  
* <a name="Includefiles"></a> <h3>Include files</h3>.
 

* 
* 本教程中包含了几个新的头文件。第一个是提供ParameterAcceptor类的声明。
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
* 这是第二个，它是一个包罗万象的头，可以让我们在这段代码中加入自动区分（AD）功能。
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
* 而接下来的三个则使用通用的 MeshWorker::mesh_loop() 框架提供一些多线程能力。
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
* 然后我们为这个程序打开一个命名空间，并将dealii命名空间中的所有东西导入其中，就像以前的程序一样。
* 

* 
* @code
 namespace Step72
 {
   using namespace dealii;
* 
 @endcode
* 
* <a name="ThecodeMinimalSurfaceProblemParameterscodeclass"></a> <h3>The <code>MinimalSurfaceProblemParameters</code> class</h3>。
 

* 
* 在本教程中，我们将实现三种不同的方法来组装线性系统。其中一个反映了最初在 step-15 中提供的手工实现，而另外两个则使用Sacado自动微分库，该库作为Trilinos框架的一部分提供。  
* 为了方便在这三种实现之间进行切换，我们有这个非常基本的参数类，它只有两个可配置的选项。
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
* 选择要使用的配方和相应的AD框架。
* 

* 
* 
* - 配方 = 0 : 无援助的实施（全手工线性化）。
* 

* 
* 
* - 配方 = 1 : 有限元残差的自动线性化。
* 

* 
* 
* - 配方 = 2 : 使用变量配方自动计算有限元残差和线性化。
* 

* 
* @code
     unsigned int formulation = 0;
* 
 @endcode
* 
* 线性系统残差的最大可接受公差。我们将看到，一旦我们使用AD框架，装配时间就会变得很明显，所以我们将 step-15 中选择的公差提高了一个数量级。这样，计算就不会花费太长时间来完成。
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
* <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a> <h3>The <code>MinimalSurfaceProblem</code> class template</h3>.
 

* 
* 该类模板与  step-15  中的内容基本相同。该类唯一的功能变化是：。
* 

* 
* 
* - run()函数现在需要两个参数：一个是选择采用哪种装配方法，一个是允许的最终残差的公差，以及
* 

* 
* 
* - 现在有三个不同的装配函数来实现三种装配线性系统的方法。我们将在后面提供关于这些的细节。
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
 
* <a name="Boundarycondition"></a> <h3>Boundary condition</h3>。
 

* 
*应用于该问题的边界条件没有变化。
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
* <a name="ThecodeMinimalSurfaceProblemcodeclassimplementation"></a> <h3>The <code>MinimalSurfaceProblem</code> class implementation</h3>.
 

* 
* <a name="MinimalSurfaceProblemMinimalSurfaceProblem"></a> <h4>MinimalSurfaceProblem::MinimalSurfaceProblem</h4>.
 

* 
* 对类的构造函数没有做任何修改。
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
 
* <a name="MinimalSurfaceProblemsetup_system"></a> <h4>MinimalSurfaceProblem::setup_system</h4>.
 

* 
* 设置类数据结构的函数没有任何变化，即DoFHandler、应用于问题的悬挂节点约束以及线性系统。
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
* <a name="Assemblingthelinearsystem"></a> <h4>Assembling the linear system</h4>
 

 
* <a name="Manualassembly"></a> <h5>Manual assembly</h5>.
 

* 
* 汇编函数是本教程的有趣贡献。assemble_system_unassisted()方法实现了与 step-15 中详述的完全相同的装配函数，但在这个例子中我们使用 MeshWorker::mesh_loop() 函数来多线程装配过程。这样做的原因很简单。当使用自动分化时，我们知道会有一些额外的计算开销产生。为了减轻这种性能损失，我们希望尽可能多地利用（容易获得的）计算资源。 MeshWorker::mesh_loop() 的概念使这成为一个相对简单的任务。同时，为了公平比较，我们需要对在计算残差或其线性化时不使用任何援助的实现做同样的事情。( MeshWorker::mesh_loop() 函数首先在 step-12 和 step-16 中讨论，如果你想阅读它的话)。  
* 实现多线程所需的步骤在这三个函数之间是相同的，所以我们将利用assemble_system_unassisted()函数作为一个机会来关注多线程本身。
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
 
*  MeshWorker::mesh_loop() 希望我们提供两个示范性的数据结构。第一个，`ScratchData`，是用来存储所有要在线程之间重复使用的大数据。`CopyData`将保存来自每个单元的对线性系统的贡献。这些独立的矩阵-向量对必须按顺序累积到全局线性系统中。由于我们不需要 MeshWorker::ScratchData 和 MeshWorker::CopyData 类已经提供的任何东西，我们使用这些确切的类定义来解决我们的问题。请注意，我们只需要一个局部矩阵、局部右手向量和单元自由度索引向量的单一实例。
* 
* - 因此 MeshWorker::CopyData 的三个模板参数都是`1`。
* 

* 
* @code
     using ScratchData = MeshWorker::ScratchData<dim>;
     using CopyData    = MeshWorker::CopyData<1, 1, 1>;
* 
 @endcode
 
* 我们还需要知道在装配过程中我们将使用什么类型的迭代器。为了简单起见，我们只要求编译器使用decltype()指定器为我们解决这个问题，知道我们将在由 @p dof_handler. 拥有的活动单元上进行迭代。
* 

* 
* @code
     using CellIteratorType = decltype(dof_handler.begin_active());
* 
 @endcode
* 
* 这里我们初始化示例数据结构。因为我们知道我们需要计算形状函数梯度、加权雅各布和实空间中的四分位点的位置，所以我们将这些标志传给类的构造函数。
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
* 现在我们定义一个lambda函数，它将在单个单元上执行装配。三个参数是由于我们将传递给那个最终调用的参数，将被 MeshWorker::mesh_loop(), 所期望的参数。我们还捕获了 @p this 指针，这意味着我们将访问 "this"（即当前的`MinimalSurfaceProblem<dim>`）类实例，以及它的私有成员数据（因为lambda函数被定义在MinimalSurfaceProblem<dim>方法内）。    
* 在函数的顶部，我们初始化依赖于正在进行工作的单元的数据结构。请注意，重新初始化的调用实际上返回了一个FEValues对象的实例，该对象被初始化并存储在 "scratch_data "对象中（因此，被重用）。    
* 同样，我们从 MeshWorker::mesh_loop() 提供的`copy_data`实例中获得本地矩阵、本地RHS向量和本地单元格DoF指数的别名。然后我们初始化单元格的DoF指数，因为我们知道本地矩阵和向量的大小已经正确。
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
* 对于牛顿方法，我们需要在问题被线性化的那一点上获得解的梯度。      
* 一旦我们有了这个梯度，我们就可以用通常的方法对这个单元进行装配。 与 step-15 的一个小区别是，我们使用了（相当方便的）基于范围的循环来迭代所有正交点和自由度。
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
*  MeshWorker::mesh_loop() 要求的第二个lambda函数是一个执行累积全局线性系统中的局部贡献的任务。这正是这个函数所做的，其实现的细节在之前已经看到过。需要认识的主要一点是，局部贡献被存储在传入该函数的`copy_data`实例中。这个`copy_data`在 @a 对`cell_worker`的一些调用中已经被填充了数据。
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
* 我们已经有了所有需要的函数定义，所以现在我们调用 MeshWorker::mesh_loop() 来执行实际的装配。 我们传递一个标志作为最后的参数，说明我们只想对单元格进行装配。在内部， MeshWorker::mesh_loop() 将可用的工作分配给不同的线程，有效地利用当今几乎所有的处理器所提供的多核。
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
* 最后，正如在  step-15  中所做的那样，我们从系统中移除悬空的节点，并对定义牛顿更新的线性系统应用零边界值  $\delta u^n$  。
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
 
* <a name="Assemblyviadifferentiationoftheresidualvector"></a> <h5>Assembly via differentiation of the residual vector</h5> 。
 

* 
* 正如在介绍中所概述的，对于这第二种方法，我们需要做的是实现从单元 $K$ 到残差向量的局部贡献 $F(U)^K$ ，然后让AD机器处理如何计算来自它的导数 $J(U)_{ij}^K=\frac{\partial F(U)^K_i}{\partial U_j}$ 。  
* 对于以下情况，请记住，@f[
 F(U)_i^K \dealcoloneq
 \int\limits_K\nabla \varphi_i \cdot \left[ \frac{1}{\sqrt{1+|\nabla
 u|^{2}}} \nabla u \right] \, dV ,
 @f]即 $u(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$ 。  
* 让我们看看这在实践中是如何实现的。
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
* 我们将利用  step-71  中所示的技术，预先定义我们要使用的AD数据结构。在这种情况下，我们选择辅助类，它将使用Sacado向前自动微分类型自动计算有限元残差的线性化。这些数字类型可以只用来计算一阶导数。这正是我们想要的，因为我们知道我们将只对残差进行线性化，这意味着我们只需要计算一阶导数。计算的返回值将是`double`类型。    
* 我们还需要一个提取器来检索一些与问题的现场解决方案有关的数据。
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
* 有了这个，让我们定义lambda函数，它将被用来计算单元格对雅各布矩阵和右手边的贡献。
 

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
* 我们现在要创建并初始化一个AD辅助类的实例。要做到这一点，我们需要指定有多少个自变量和因变量。自变量将是我们的解向量所具有的局部自由度的数量，即离散化解向量 $u (\mathbf{x})|_K = \sum\limits_{j} U^K_i \varphi_j(\mathbf{x})$ 的每元素表示中的数字 $j$ ，它表示每个有限元素有多少个解系数。在deal.II中，这等于 FiniteElement::dofs_per_cell. ，自变量的数量将是我们将要形成的局部残差向量中的条目数。在这个特定的问题中（就像许多其他采用[标准Galerkin方法](https://en.wikipedia.org/wiki/Galerkin_method)的问题一样），局部求解系数的数量与局部残差方程的数量相符。
* 

* 
* @code
       const unsigned int n_independent_variables = local_dof_indices.size();
       const unsigned int n_dependent_variables   = dofs_per_cell;
       ADHelper ad_helper(n_independent_variables, n_dependent_variables);
* 
 @endcode
 
* 接下来，我们将解的数值告知帮助器，即我们希望线性化的 $U_j$ 的实际数值。由于这是在每个元素上单独进行的，我们必须从全局解决方案向量中提取解决方案的系数。换句话说，我们将所有这些系数 $U_j$ （其中 $j$ 是一个局部自由度）定义为进入向量 $F(U)^{K}$ （因果函数）计算的自变量。      
*然后我们得到了由自动可分的数字表示的自由度值的完整集合。对这些变量进行的操作从这一点开始被AD库跟踪，直到对象超出范围。所以正是这些变量 <em>  </em> ，我们将对其计算残差分录的导数。
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
* 然后我们做一些特定问题的任务，首先是计算所有的值，（空间）梯度，以及基于 "敏感的 "AD自由度值的类似东西。在这个例子中，我们要检索每个正交点的解梯度。请注意，现在解梯度对自由度值很敏感，因为它们使用 @p ADNumberType 作为标量类型， @p dof_values_ad 矢量提供局部自由度值。
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
* 我们声明的下一个变量将存储单元格残余向量贡献。这是相当不言自明的，除了一个<b>very important</b>细节。请注意，向量中的每个条目都是手工初始化的，数值为0。这是一个 <em> 强烈推荐的 </em> 做法，因为一些AD库似乎没有安全地初始化这些数字类型的内部数据结构。不这样做可能会导致一些非常难以理解或检测的bug（感谢这个程序的作者出于一般的坏经验而提到这一点）。因此，出于谨慎考虑，值得明确地将初始值归零。在这之后，除了符号的改变，残差集看起来和我们之前看到的单元格RHS向量差不多。我们在所有正交点上循环，确保系数现在通过使用正确的`ADNumberType'来编码它对（敏感的）有限元DoF值的依赖性，最后我们组装残差向量的组件。为了完全清楚，有限元形状函数（及其梯度等）以及 "JxW "值仍然是标量值，但每个正交点的 @p coeff 和 @p old_solution_gradients 是以独立变量计算的。
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
* 一旦我们计算出了全单元的残差向量，我们就可以把它注册到帮助类中。      
* 此后，我们在评估点计算残差值（基本上是从我们已经计算出来的东西中提取出真实的值）和它们的雅各布（每个残差分量相对于所有单元DoF的线性化）。为了组装成全局线性系统，我们必须尊重残差和RHS贡献之间的符号差异。对于牛顿方法，右手边的矢量需要等于负的*残差矢量。
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
* 该函数的余数等于我们之前的结果。
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
* <a name="Assemblyviadifferentiationoftheenergyfunctional"></a> <h5>Assembly via differentiation of the energy functional</h5>。
 

* 
* 在这第三种方法中，我们将残差和雅各布作为局部能量函数@f[
 E\left( U \right)^K
 \dealcoloneq \int\limits_{K} \Psi \left( u \right) \, dV
 \approx \sum\limits_{q}^{n_{\textrm{q-points}}} \Psi \left( u \left(
 \mathbf{X}_{q} \right) \right) \underbrace{\vert J_{q} \vert \times
 W_{q}}_{\text{JxW(q)}}
 @f]的第一和第二导数来计算，能量密度由@f[
 \Psi \left( u \right) = \sqrt{1+|\nabla u|^{2}} .
 @f]给出。  
* 让我们再次看看这是如何做到的。
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
 
* 在这个装配过程的实现中，我们选择了帮助类，它将使用嵌套的Sacado前向自动微分类型自动计算残差及其从单元贡献到能量函数的线性化。所选的数字类型可以用来计算第一和第二导数。我们需要这样做，因为残差定义为势能对DoF值的敏感性（即其梯度）。然后我们需要将残差线性化，这意味着必须计算势能的二阶导数。你可能想把这与之前函数中使用的 "ADHelper "的定义进行比较，在那里我们使用了 `Differentiation::AD::ResidualLinearization<Differentiation::AD::NumberTypes::sacado_dfad,double>`.  。
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
* 然后让我们再次定义lambda函数，对单元格进行整合。    
* 为了初始化辅助类的实例，我们现在只需要预先知道自变量的数量（也就是与元素解向量相关的自由度数量）。这是因为由能量函数产生的二阶导数矩阵必然是平方的（而且，顺便说一下，也是对称的）。
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
 
* 再一次，我们用帮助器注册所有的单元格DoFs值，然后提取这些值的 "敏感 "变体，这些变体将被用于必须进行区分的后续操作中。
* 
* - 其中之一是计算解决方案的梯度。
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
 
* 我们接下来创建一个变量来存储单元格的总能量。我们再次强调，我们明确地对这个值进行了零初始化，从而确保了这个起始值的数据的完整性。      
* 我们方法的目的是计算细胞总能量，这是内部（由于右手函数，通常是 $U$ 中的线性）和外部能量的总和。在这种特殊情况下，我们没有外部能量（例如，来自源项或诺伊曼边界条件），所以我们将关注内部能量部分。      
*事实上，计算 $E(U)^K$ 几乎是微不足道的，只需要以下几行。
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
* 在我们计算出这个单元上的总能量后，我们将在帮助器上注册它。 在此基础上，我们现在可以计算出所需的数量，即残差值和它们在评估点的雅各布系数。和以前一样，牛顿的右手边需要是残差的负数。
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
* 和前面两个函数一样，函数的剩余部分也是如此。
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
* <a name="MinimalSurfaceProblemsolve"></a> <h4>MinimalSurfaceProblem::solve</h4>。
 

* 
* 求解函数与  step-15  中使用的相同。
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
* <a name="MinimalSurfaceProblemrefine_mesh"></a> <h4>MinimalSurfaceProblem::refine_mesh</h4>.
 

* 
* 自 step-15 以来，在网格细化程序和适应性网格之间的解决方案的转移方面没有任何变化。
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
* <a name="MinimalSurfaceProblemset_boundary_values"></a> <h4>MinimalSurfaceProblem::set_boundary_values</h4>.
 

* 
* 边界条件的选择仍然与 step-15 相同 ...
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
* <a name="MinimalSurfaceProblemcompute_residual"></a> <h4>MinimalSurfaceProblem::compute_residual</h4>。
 

* 
*......就像在求解迭代过程中用来计算残差的函数一样。如果真的需要，我们可以用能量函数的微分来代替它，但为了简单起见，我们在这里简单地复制我们在 step-15 中已经有的东西。
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
* <a name="MinimalSurfaceProblemdetermine_step_length"></a> <h4>MinimalSurfaceProblem::determine_step_length</h4>.
 

* 
* 非线性迭代程序的步长（或欠松系数）的选择仍然固定在 step-15 中选择和讨论的值。
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
* <a name="MinimalSurfaceProblemoutput_results"></a> <h4>MinimalSurfaceProblem::output_results</h4>
 

* 
* 从`run()`调用的最后一个函数以图形形式输出当前的解决方案（和牛顿更新），作为VTU文件。它与以前的教程中所使用的完全相同。
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
* <a name="MinimalSurfaceProblemrun"></a> <h4>MinimalSurfaceProblem::run</h4>。
 

* 
* 在运行函数中，大部分仍然与首次实现的  step-15  相同。唯一可以观察到的变化是，我们现在可以（通过参数文件）选择系统残差的最终可接受的公差是什么，并且我们可以选择我们希望利用的装配方法。为了使第二个选择明确，我们向控制台输出一些信息，表明选择。由于我们对比较三种方法中每一种的装配时间感兴趣，我们还添加了一个计时器，跟踪装配过程中所花费的时间。我们还跟踪了解决线性系统所需的时间，这样我们就可以将这些数字与代码中通常需要最长时间执行的部分进行对比。
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
* <a name="Themainfunction"></a> <h4>The main function</h4> <h4>The main function</h4>.
 

* 
* 最后是主函数。它遵循大多数其他主函数的方案，但有两个明显的例外。
 

* 
* 
* - 我们调用 Utilities::MPI::MPI_InitFinalize ，以便（通过一个隐藏的默认参数）设置使用多线程任务执行的线程数。
* 

* 
* 
* - 我们还有几行专门用于读入或初始化用户定义的参数，这些参数将在程序执行过程中被考虑。
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
* <a name="Results"></a><h1>Results</h1>.
 

* 由于在 step-15 中首次分析的问题的物理学没有变化，所以没有什么可报告的。它们之间唯一明显的区别是，在默认情况下，这个程序只运行9个网格细化步骤（而 step-15 则执行11个细化步骤）。下面报告的所有时间都是在发布模式下获得的）。
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
* 以下是手部线性化的输出（在2012年左右的四核八线程笔记本上计算）。
* 
以下是手部线性化的输出结果（在2012年左右的四核八线程笔记本电脑上计算的 * - 但我们真正感兴趣的是不同实现方式之间的相对时间）。)
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
* 而对于使用Sacado动态正向AD数字类型，以自动方式线性化残差的实现。
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
* 最后，对于直接从能量函数（使用嵌套的萨卡多动态正向AD数）计算残差和其线性化的实现。
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
* 很明显，交给自动区分框架执行的工作越多，在装配过程中花费的时间就越多。累计所有的细化步骤，使用一级自动微分导致在装配阶段花费了 $5.65 \times$ 的计算时间，与无辅助装配相比，而直接从能量函数推导出离散线性系统时，装配时间要长 $10.7 \times$ .毫不奇怪，解决线性系统的总体时间保持不变。对许多人来说，这可能意味着在生产代码中利用高阶微分（在有限元水平）会导致不可接受的开销，但在原型设计阶段可能仍然有用。因此，两者之间的一个很好的折衷办法可能是有限元残差的自动线性化，它以可衡量但也许不是不可接受的成本提供了很多便利。另外，我们可以考虑不在每一步中重新建立牛顿矩阵。
* 
* --这个话题在 step-77 中进行了大量的深入探讨。
* 当然，在实践中，实际的开销在很大程度上取决于被评估的问题（例如，解中有多少成分，被微分的函数的性质是什么，等等）。因此，这里提出的确切结果应该仅在这个标量问题的背景下进行解释，当涉及到其他问题时，用户肯定需要进行一些初步调查。
* 

*<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>
* 

* 与 step-71 一样，有几个与自动区分有关的项目可以进一步评估。
* 
* - 应调查其他AD框架的使用情况，展望其他实现方式可能提供性能优势。
* 
* - 除了本教程中硬编码的数字类型外，还值得评估AD数字类型。关于在有限元水平上采用的两次可微分类型，混合微分模式（"RAD-FAD"）原则上应该比这里采用的单模式（"FAD-FAD"）类型更有计算效率。RAD-FAD类型没有被默认选择的原因是，在撰写本文时，在Sacado库中，它的实现仍然存在一些错误，导致内存泄漏。 这在 @ref auto_symb_diff 模块中有所记载。
* 
* - 可能的情况是，使用降低精度的类型（即`float`）作为AD数字的标量类型，可以在装配时减少计算费用。使用 "float "作为矩阵和残差的数据类型并不是不合理的，因为牛顿更新只是为了让我们更接近解决方案，而不是真正的*解决方案；因此，考虑使用降低精度的数据类型来计算这些更新，然后将这些更新累积到使用完全 "双 "精度的解决方案向量中，是有意义的。
* 
* - 在装配过程中可能减少资源的另一种方法是将AD的实现作为一个构成性模型。这类似于 step-71 中采用的方法，并将自动分化的起点推到计算链的上一级。这反过来又意味着AD库跟踪的操作更少，从而降低了微分的成本（尽管我们会在每个单元正交点进行微分）。
* 
* -  step-77 是 step-15 的另一个变体，解决了问题的一个非常不同的部分：直线搜索以及是否有必要在每次非线性迭代中重新建立牛顿矩阵。鉴于上述结果表明，使用自动微分是有代价的， step-77 中的技术有可能在一定程度上抵消这些代价。因此，将目前的程序与 step-77 中的技术结合起来是相当有趣的。对于生产代码来说，这肯定是个好办法。
* 

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-72.cc"  。
* */


