//include/deal.II-translator/A-tutorial/step-63_0.txt
/**
 * @page step_63 The step-63 tutorial program
 * 本教程依赖于  step-16  。 @htmlonly <table class="tutorial"
 * width="50%"> <tr><th colspan="2"><b><small>Table of
 * contents</small></b><b><small>Table of contents</small></b></th></tr>
 * <tr><td width="50%" valign="top">
 * <ol>
 * <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro"
 * class=bold>Introduction</a>
 * <ul>
 * <li><a href="#Equation">Equation</a><a href="#Equation">Equation</a>
 * <ul>
 * <li><a href="#Streamlinediffusion">Streamline diffusion</a><a
 * href="#Streamlinediffusion">Streamline diffusion</a>
 * </ul>
 * <li><a href="#Smoothers">Smoothers</a><a href="#Smoothers">Smoothers</a>
 * <li><a href="#Testproblem">Test problem</a><a href="#Testproblem">Test
 * problem</a>
 * </ul>
 * <li> <a href="#CommProg" class=bold>The commented program</a><a
 * href="#CommProg" class=bold>The commented program</a>
 * <ul>
 * <li><a href="#Includefiles">Include files</a><a
 * href="#Includefiles">Include files</a>
 * <li><a href="#MeshWorkerdata">MeshWorker data</a><a
 * href="#MeshWorkerdata">MeshWorker data</a>
 * <li><a href="#Problemparameters">Problem parameters</a><a
 * href="#Problemparameters">Problem parameters</a>
 * <li><a href="#Cellpermutations">Cell permutations</a><a
 * href="#Cellpermutations">Cell permutations</a>
 * <li><a href="#Righthandsideandboundaryvalues">Right-hand side and boundary
 * values</a><a href="#Righthandsideandboundaryvalues">Right-hand side and
 * boundary values</a>
 * <li><a href="#Streamlinediffusionimplementation">Streamline diffusion
 * implementation</a><a href="#Streamlinediffusionimplementation">Streamline
 * diffusion implementation</a>
 * <li><a href="#codeAdvectionProlemcodeclass"><code>AdvectionProlem</code>
 * class</a><a
 * href="#codeAdvectionProlemcodeclass"><code>AdvectionProlem</code> class</a>
 * <ul>
 * <li><a
 * href="#codeAdvectionProblemsetup_systemcode"><code>AdvectionProblem::setup_system()</code></a><a
 * href="#codeAdvectionProblemsetup_systemcode"><code>AdvectionProblem::setup_system()</code></a>
 * <li><a
 * href="#codeAdvectionProblemassemble_cellcode"><code>AdvectionProblem::assemble_cell()</code></a><a
 * href="#codeAdvectionProblemassemble_cellcode"><code>AdvectionProblem::assemble_cell()</code></a>
 * <li><a
 * href="#codeAdvectionProblemassemble_system_and_multigridcode"><code>AdvectionProblem::assemble_system_and_multigrid()</code></a><a
 * href="#codeAdvectionProblemassemble_system_and_multigridcode"><code>AdvectionProblem::assemble_system_and_multigrid()</code></a>
 * <li><a
 * href="#codeAdvectionProblemsetup_smoothercode"><code>AdvectionProblem::setup_smoother()</code></a><a
 * href="#codeAdvectionProblemsetup_smoothercode"><code>AdvectionProblem::setup_smoother()</code></a>
 * <li><a
 * href="#codeAdvectionProblemsolvecode"><code>AdvectionProblem::solve()</code></a><a
 * href="#codeAdvectionProblemsolvecode"><code>AdvectionProblem::solve()</code></a>
 * <li><a
 * href="#codeAdvectionProblemoutput_resultscode"><code>AdvectionProblem::output_results()</code></a><a
 * href="#codeAdvectionProblemoutput_resultscode"><code>AdvectionProblem::output_results()</code></a>
 * <li><a
 * href="#codeAdvectionProblemruncode"><code>AdvectionProblem::run()</code></a><a
 * href="#codeAdvectionProblemruncode"><code>AdvectionProblem::run()</code></a>
 * </ul>
 * <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a><a
 * href="#Thecodemaincodefunction">The <code>main</code> function</a>
 * </ul>
 * </ol></td><td width="50%" valign="top"><ol>
 * <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results"
 * class=bold>Results</a>
 * <ul>
 * <li><a href="#GMRESIterationNumbers"> GMRES Iteration Numbers </a><a
 * href="#GMRESIterationNumbers"> GMRES Iteration Numbers </a>
 * <ul>
 * <li><a href="#DoFCellRenumbering"> DoF/Cell Renumbering </a><a
 * href="#DoFCellRenumbering"> DoF/Cell Renumbering </a>
 * <li><a href="#Pointvsblocksmoothers"> Point vs. block smoothers </a><a
 * href="#Pointvsblocksmoothers"> Point vs. block smoothers </a>
 * </ul>
 * <li><a href="#Cost"> Cost </a><a href="#Cost"> Cost </a>
 * <li><a href="#Additionalpoints"> Additional points </a><a
 * href="#Additionalpoints"> Additional points </a>
 * <li><a href="#Possibilitiesforextensions"> Possibilities for extensions
 * </a><a href="#Possibilitiesforextensions"> Possibilities for extensions
 * </a>
 * <ul>
 * <li><a href="#ConstantiterationsforQsub5sub"> Constant iterations for
 * Q<sub>5</sub> </a><a href="#ConstantiterationsforQsub5sub"> Constant
 * iterations for Q<sub>5</sub> </a>
 * <li><a href="#Effectivenessofrenumberingforchangingepsilon"> Effectiveness
 * of renumbering for changing epsilon </a><a
 * href="#Effectivenessofrenumberingforchangingepsilon"> Effectiveness of
 * renumbering for changing epsilon </a>
 * <li><a href="#Meshadaptivity"> Mesh adaptivity </a><a
 * href="#Meshadaptivity"> Mesh adaptivity </a>
 * </ul>
 * </ul>
 * <li> <a href="#PlainProg" class=bold>The plain program</a><a
 * href="#PlainProg" class=bold>The plain program</a>
 * </ol> </td> </tr> </table>
 * @endhtmlonly
 * <br>
 * <i>This program was contributed by Thomas C. Clevenger and Timo Heister.
 *
 * The creation of this tutorial was partially supported by NSF Award
 * DMS-1522191, DMS-1901529, OAC-1835452, by the Computational Infrastructure
 * in Geodynamics initiative (CIG), through the NSF under Award EAR-0949446
 * and EAR-1550901 and The University of California
 *
 * - Davis. </i>
 * @dealiiTutorialDOI{10.5281/zenodo.3382899,https://zenodo.org/badge/DOI/10.5281/zenodo.3382899.svg}
 * <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1> 。
 * 该程序使用几何多网格(GMG)预处理程序解决了一个平流-扩散问题。这个预处理程序的基本原理在
 * step-16
 * 中讨论过；这里我们讨论非对称PDE所需的必要变化。此外，我们还介绍了块平滑的概念（与
 * step-16
 * 中的点平滑相比），并研究了加法和乘法平滑器的DoF重编号的效果。
 * <a name="Equation"></a><h3>Equation</h3>
 * 平流-扩散方程由以下公式给出
 * @f{align*}{
 *
 * -\varepsilon \Delta u + \boldsymbol{\beta}\cdot \nabla u & = f &
 * \text{ in } \Omega\\
 * u &= g & \text{ on } \partial\Omega
 * @f}
 * 其中 $\varepsilon>0$  ， $\boldsymbol{\beta}$ 是<i>advection
 * direction</i>，而 $f$ 是一个源。有几个注意事项。 1.如果
 * $\boldsymbol{\beta}=\boldsymbol{0}$ ，这是在 step-16
 * 中解决的拉普拉斯方程（以及其他许多地方）。 2.如果
 * $\varepsilon=0$ ，这是在 step-9 中解决的固定平流方程。
 * 3.人们可以为这个问题定义一个无维数，称为<i>Peclet
 * number</i>。   $\mathcal{P} \dealcoloneq \frac{\|\boldsymbol{\beta}\|
 * L}{\varepsilon}$  ，其中 $L$
 * 是域的长度尺度。它描述了我们所考虑的那种方程。如果
 * $\mathcal{P}>1$
 * ，我们说问题是<i>advection-dominated</i>，否则如果
 * $\mathcal{P}<1$ ，我们会说问题是<i>diffusion-dominated</i>。
 * 在本教程的讨论中，我们将关注以对流为主的流动。这是很复杂的情况。我们知道，对于扩散主导的问题，标准的Galerkin方法可以很好地工作，我们也知道简单的多网格方法，如
 * step-16
 * 中定义的方法，非常有效。另一方面，对于以对流为主的问题，标准Galerkin方法会导致震荡和不稳定的离散，而且简单的求解器往往不是很有效。因此，本教程的目的是解决这两个问题。
 *
 *  <a name="Streamlinediffusion"></a><h4>Streamline diffusion</h4>
 * 使用标准的Galerkin有限元方法，对于合适的测试函数  $v_h$
 * ，PDE的离散弱形式将为
 * @f{align*}{
 * a(u_h,v_h) = F(v_h)
 * @f}
 * 其中
 * @f{align*}{
 * a(u_h,v_h) &= (\varepsilon \nabla v_h,\, \nabla u_h) +
 * (v_h,\,\boldsymbol{\beta}\cdot \nabla u_h),\\
 * F(v_h) &= (v_h,\,f).
 * @f}
 *
 * 不幸的是，用这种方法通常会得到震荡解。事实上，对于这种形式，可以显示出以下误差估计。
 * @f{align*}{
 * \|\nabla (u-u_h)\| \leq (1+\mathcal{P}) \inf_{v_h} \|\nabla (u-v_h)\|.
 * @f}
 * 如果精确解足够平滑，右边的下限可以被估计为如下。
 * @f{align*}{
 * \inf_{v_h} \|\nabla (u-v_h)\|.
 * \le
 * \|\nabla (u-I_h u)\|
 * \le
 * h^k
 * C
 * \|\nabla^k u)\|
 * @f}
 * 其中 $k$
 * 是所用有限元的多项式程度。因此，我们得到的估计是
 * @f{align*}{
 * \|\nabla (u-u_h)\|
 * \leq (1+\mathcal{P}) C h^k
 * \|\nabla^k u)\|.
 * @f}
 * 换句话说，数值解会收敛。另一方面，鉴于上述
 * $\mathcal{P}$ 的定义，我们不得不期待在 $\varepsilon \ll
 * \|\boldsymbol{\beta}\| L$
 * 时，即在问题只有少量扩散的情况下，有较大误差的数值解。
 * 为了解决这个问题，我们将考虑新的弱形式
 * @f{align*}{
 * a(u_h,\,v_h) + \sum_K (-\varepsilon \Delta u_h +
 * \boldsymbol{\beta}\cdot \nabla u_h-f,\,\delta_K
 * \boldsymbol{\beta}\cdot \nabla v_h)_K = F(v_h)
 * @f}
 * 其中，总和是在所有单元 $K$
 * 上进行的，每个单元都有内积， $\delta_K$ 是定义在 @cite
 * john2006discontinuity 中的一个单元常数稳定参数。
 * 本质上，加入离散强形式残差会增强双线性形式
 * $a(\cdot,\cdot)$
 * 的矫捷性，从而提高离散解的稳定性。这种方法通常被称为<i>streamline
 * diffusion</i>或<i>SUPG</i>（流线上风/Petrov-Galerkin）。
 *
 *  <a name="Smoothers"></a><h3>Smoothers</h3>
 * 本教程的目标之一是扩展使用一个简单的（point-wise）高斯-赛德尔（SOR）平滑器，该平滑器在
 * step-16
 * 中使用。对于标量问题，这意味着使用一个求解器，每次更新线性系统的一个未知数，同时保持所有其他未知数的固定；然后对问题中的所有未知数进行迭代，一旦完成，从第一个未知数重新开始，直到这些
 * "扫频
 * "收敛。Jacobi、Gauss-Seidel和SOR迭代都可以用这种方式来解释。在多重网格的背景下，人们不认为这些方法是
 * "求解器"，而是
 * "平滑器"。因此，人们对实际解决线性系统不感兴趣。为了使多网格方法发挥作用，只需去除残差中的高频部分即可，因为这样可以将解限制在更粗的网格中。
 * 因此，我们只需对所有未知数进行少数固定数量的
 * "扫描"。在本教程的代码中，这是由 "平滑步骤
 * "参数控制的。
 * 但众所周知，这些方法在用于解析器时收敛得相当慢。虽然作为多网格平滑器，它们出乎意料地好，但它们也可以被改进。特别是，我们在这里也考虑
 * "基于单元
 * "的平滑器。这些方法一次解决一个单元的所有未知数，保持所有其他未知数的固定；然后转到下一个单元，如此反复。我们可以把它们看作是雅可比（Jacobi）、高斯-赛德尔（Gauss-Seidel）或SOR的
 * "区块
 * "版本，但由于自由度是由多个单元共享的，这些区块是重叠的，这些方法实际上最好是在加法和乘法施瓦茨方法的框架内解释。
 * 与 step-16
 * 相比，我们的测试问题包含一个平流项。特别是在扩散常数较小的情况下
 * $\varepsilon$
 * ，信息是沿着给定的平流方向的流线进行传输的。这意味着，如果平滑器在一次平滑应用中允许信息沿下游方向传播，那么平滑器可能会更有效。如果我们想按照这些未知数（或未知数块）被列举的顺序每次解决一个未知数（或未知数块），那么这种信息传播特性要求我们相应地重新排列自由度或单元（对于基于单元的平滑器），以便更上游的自由度被提前处理（指数较低），更下游的自由度被推迟处理（指数较大）。排序的影响将在结果部分显现。
 * 现在让我们简单地定义一下本教程中使用的平滑器。关于更详细的介绍，我们参考
 * @cite KanschatNotesIterative 和书籍 @cite smith2004domain 和 @cite
 * toselli2006domain 。 一个Schwarzpreconditioner需要一个分解
 * @f{align*}{
 * V = \sum_{j=1}^J V_j
 * @f}
 * 我们的有限元空间  $V$  。每个子问题  $V_j$
 * 也有一个基于双线性形式  $a(\cdot,\cdot)$  的里兹投影  $P_j:
 * V \rightarrow V_j$  。这个投影为每个子问题 $V_j$
 * 诱导了一个局部算子 $A_j$ 。如果 $\Pi_j:V\rightarrow V_j$ 是对
 * $V_j$ 的正交投影，可以证明 $A_jP_j=\Pi_j^TA$  。
 * 由此，我们可以为运算符 $A$ 定义一个<i>additive Schwarz
 * preconditioner</i>为
 * @f{align*}{
 * B^{-1} = \sum_{j=1}^J P_j A^{-1} = \sum_{j=1}^J A_j^{-1} \Pi_j^T.
 * @f}
 * 换句话说，我们将我们的解决方案投射到每个子问题中，应用子问题
 * $A_j$ 的逆，并将所有 $j$ 的贡献相加。
 * 请注意，我们可以通过为每个自由度定义一个子问题 $V_j$
 * ，将逐点（一次一个未知数）的雅可比方法解释为加性施瓦茨方法。然后，
 * $A_j^{-1}$ 成为与 $A$ 的对角线条目的逆数相乘的结果。
 * 对于本教程中使用的 "Block Jacobi
 * "方法，我们为当前层次上的每个网格单元定义一个子问题
 * $V_j$
 * 。注意，我们使用的是连续有限元，所以这些块是重叠的，因为两个单元之间的界面上的自由度属于两个子问题。在子问题上操作施瓦兹算子的逻辑（在deal.II中它们被称为
 * "块"）是在RelaxationBlock类中实现的。BlockJacobi
 * "方法在RelaxationBlockJacobi类中实现。该类的许多方面（例如如何定义块以及如何反转局部子问题
 * $A_j$ ）可以在平滑数据中配置，详见
 * RelaxationBlock::AdditionalData 和 DoFTools::make_cell_patches() 。
 * 到目前为止，我们讨论的是加法平滑器，其中的更新可以独立应用，并且在单一平滑器应用中没有信息流动。一个<i>multiplicative
 * Schwarz preconditioner</i>解决了这个问题，其定义为
 * @f{align*}{
 * B^{-1} = \left( I- \prod_{j=1}^J \left(I-P_j\right) \right) A^{-1}.
 * @f}
 * 与上述不同的是，对子问题 $V_j$
 * 的更新是依次应用的。这意味着在反转子问题 $A_j$
 * 时得到的更新立即用于 $A_{j+1}$
 * 。这在写出项目时就可以看到。
 * @f{align*}{
 * B^{-1}
 * =
 * \left(
 * I
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
 * -
 * \left(I-P_1\right)\left(I-P_2\right)\cdots\left(I-P_J\right)
 * \right)
 * A^{-1}
 * =
 * A^{-1}
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
 * -
 * \left[ \left(I-P_1\right)
 * \left[ \left(I-P_2\right)\cdots
 *   \left[\left(I-P_J\right) A^{-1}\right] \cdots \right] \right]
 * @f}
 *
 * 当把子空间 $V_j$
 * 定义为整个自由度块时，这个方法在RelaxationBlockSOR类中实现，当你在本教程中选择
 * "Block SOR
 * "时使用。RelaxationBlockSOR类也是从RelaxationBlock派生出来的。因此，加法和乘法施瓦茨方法都在一个统一的框架内实现。
 * 最后，让我们注意到，标准的高斯-赛德尔（或SOR）方法可以被看作是一个乘法施瓦兹方法，每个DoF都有一个子问题。
 *
 *  <a name="Testproblem"></a><h3>Test problem</h3>
 * 我们将考虑以下测试问题： $\Omega =
 * [-1,\,1]\times[-1,\,1]\backslash B_{0.3}(0)$
 * ，即一个以原点为中心的半径为0.3的正方形被移除。此外，我们使用
 * $\varepsilon=0.005$ 、 $\boldsymbol{\beta} = [-\sin(\pi/6),\,\cos(\pi/6)]$
 * 、 $f=0$ 和迪里切特边界值。
 * @f{align*}{
 * g = \left\{\begin{array}{ll} 1 & \text{if } x=-1 \text{ or } y=-1,\,x\geq 0.5 \\
 * 0 & \text{otherwise} \end{array}\right.
 * @f}
 *
 * 以下数字描述了有（左）和无（右）流线扩散的解决方案。在没有流线扩散的情况下，我们看到边界层周围有很大的振荡，这表明标准Galerkin有限元方法对这个问题的不稳定性。
 * <table width="60%" align="center"> <tr> <td align="center"> <img
 * src="https://www.dealii.org/images/steps/developer/step-63-solution.png"
 * alt=""> </td> <td align="center"> <img
 * src="https://www.dealii.org/images/steps/developer/step-63-solution-no-sd.png"
 * alt=""> </td> </tr> </table>
 *
 *  <a name="CommProg"></a> <h1> The commented program</h1>。 <a
 * name="Includefiles"></a> <h3>Include files</h3>。
 *
 * 标准交易所需的典型文件.II:
 *
 * @code
 * #include <deal.II/base/tensor_function.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/parameter_handler.h>
 *
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/relaxation_block.h>
 *
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_out.h>
 *
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 *
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/fe/fe_values.h>
 *
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 *
 * @endcode
 *
 * 包括所有相关的多层次文件。
 *
 * @code
 * #include <deal.II/multigrid/mg_constrained_dofs.h>
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 *
 * @endcode
 *
 * C++:
 *
 * @code
 * #include <algorithm>
 * #include <fstream>
 * #include <iostream>
 * #include <random>
 *
 * @endcode
 *  我们将使用 MeshWorker::mesh_loop 功能来组装矩阵。
 *
 * @code
 * #include <deal.II/meshworker/mesh_loop.h>
 *
 *
 * @endcode
 *
 * <a name="MeshWorkerdata"></a> <h3>MeshWorker data</h3>。
 *
 * 像往常一样，我们将把与这个程序有关的所有东西放入一个自己的命名空间。
 *
 *
 * 由于我们将使用MeshWorker框架，第一步是定义以下结构，这些结构是由
 * MeshWorker::mesh_loop(): 使用的assemble_cell()函数所需要的
 * `ScratchData`包含一个FEValues对象，这是组装一个单元的局部贡献所需要的，而`CopyData`包含一个单元的局部贡献的输出和复制到全局系统的必要信息。它们的目的在WorkStream类的文档中也有解释）。
 *
 * @code
 * namespace Step63
 * {
 * using namespace dealii;
 *
 * template <int dim>
 * struct ScratchData
 * {
 *   ScratchData(const FiniteElement<dim> &fe,
 *               const unsigned int        quadrature_degree)
 *     : fe_values(fe,
 *                 QGauss<dim>(quadrature_degree),
 *                 update_values | update_gradients | update_hessians |
 *                   update_quadrature_points | update_JxW_values)
 *   {}
 *
 *   ScratchData(const ScratchData<dim> &scratch_data)
 *     : fe_values(scratch_data.fe_values.get_fe(),
 *                 scratch_data.fe_values.get_quadrature(),
 *                 update_values | update_gradients | update_hessians |
 *                   update_quadrature_points | update_JxW_values)
 *   {}
 *
 *   FEValues<dim> fe_values;
 * };
 *
 *
 *
 * struct CopyData
 * {
 *   CopyData() = default;
 *
 *   unsigned int level;
 *   unsigned int dofs_per_cell;
 *
 *   FullMatrix<double>                   cell_matrix;
 *   Vector<double>                       cell_rhs;
 *   std::vector<types::global_dof_index> local_dof_indices;
 * };
 *
 *
 *
 * @endcode
 *
 * <a name="Problemparameters"></a> <h3>Problem parameters</h3>。
 * 第二步是定义处理要从输入文件中读取的运行时参数的类。
 * 我们将使用ParameterHandler在运行时传入参数。结构`Settings`解析并存储整个程序要查询的参数。
 *
 * @code
 * struct Settings
 * {
 *   enum DoFRenumberingStrategy
 *   {
 *     none,
 *     downstream,
 *     upstream,
 *     random
 *   };
 *
 *   void get_parameters(const std::string &prm_filename);
 *
 *   double                 epsilon;
 *   unsigned int           fe_degree;
 *   std::string            smoother_type;
 *   unsigned int           smoothing_steps;
 *   DoFRenumberingStrategy dof_renumbering;
 *   bool                   with_streamline_diffusion;
 *   bool                   output;
 * };
 *
 *
 *
 * void Settings::get_parameters(const std::string &prm_filename)
 * {
 *   // First declare the parameters...
 *   ParameterHandler prm;
 *
 *   prm.declare_entry("Epsilon",
 *                     "0.005",
 *                     Patterns::Double(0),
 *                     "Diffusion parameter");
 *
 *   prm.declare_entry("Fe degree",
 *                     "1",
 *                     Patterns::Integer(1),
 *                     "Finite Element degree");
 *   prm.declare_entry("Smoother type",
 *                     "block SOR",
 *                     Patterns::Selection("SOR|Jacobi|block SOR|block Jacobi"),
 *                     "Select smoother: SOR|Jacobi|block SOR|block Jacobi");
 *   prm.declare_entry("Smoothing steps",
 *                     "2",
 *                     Patterns::Integer(1),
 *                     "Number of smoothing steps");
 *   prm.declare_entry(
 *     "DoF renumbering",
 *     "downstream",
 *     Patterns::Selection("none|downstream|upstream|random"),
 *     "Select DoF renumbering: none|downstream|upstream|random");
 *   prm.declare_entry("With streamline diffusion",
 *                     "true",
 *                     Patterns::Bool(),
 *                     "Enable streamline diffusion stabilization: true|false");
 *   prm.declare_entry("Output",
 *                     "true",
 *                     Patterns::Bool(),
 *                     "Generate graphical output: true|false");
 *
 *   // ...and then try to read their values from the input file:
 *   if (prm_filename.empty())
 *     {
 *       prm.print_parameters(std::cout, ParameterHandler::Text);
 *       AssertThrow(
 *         false, ExcMessage("Please pass a .prm file as the first argument!"));
 *     }
 *
 *   prm.parse_input(prm_filename);
 *
 *   epsilon         = prm.get_double("Epsilon");
 *   fe_degree       = prm.get_integer("Fe degree");
 *   smoother_type   = prm.get("Smoother type");
 *   smoothing_steps = prm.get_integer("Smoothing steps");
 *
 *   const std::string renumbering = prm.get("DoF renumbering");
 *   if (renumbering == "none")
 *     dof_renumbering = DoFRenumberingStrategy::none;
 *   else if (renumbering == "downstream")
 *     dof_renumbering = DoFRenumberingStrategy::downstream;
 *   else if (renumbering == "upstream")
 *     dof_renumbering = DoFRenumberingStrategy::upstream;
 *   else if (renumbering == "random")
 *     dof_renumbering = DoFRenumberingStrategy::random;
 *   else
 *     AssertThrow(false,
 *                 ExcMessage("The <DoF renumbering> parameter has "
 *                            "an invalid value."));
 *
 *   with_streamline_diffusion = prm.get_bool("With streamline diffusion");
 *   output                    = prm.get_bool("Output");
 * }
 *
 *
 * @endcode
 *
 * <a name="Cellpermutations"></a> <h3>Cell permutations</h3>.
 * 遍历单元和自由度的顺序将对乘法的收敛速度起作用。这里我们定义了一些函数，这些函数返回单元格的特定排序，供块平滑器使用。
 * 对于每种类型的单元格排序，我们定义了一个用于活动网格的函数和一个用于水平网格的函数（即用于多网格层次结构中的一个水平的单元格）。虽然求解系统所需的唯一重新排序是在水平网格上进行的，但我们在output_results()中包含了活动网格的重新排序，以达到可视化的目的。
 * 对于两个下游排序函数，我们首先创建一个包含所有相关单元的数组，然后使用一个
 * "比较器
 * "对象在下游方向进行排序。然后，函数的输出是一个简单的数组，包含了刚刚计算出的排序中的单元格的索引。
 *
 * @code
 * template <int dim>
 * std::vector<unsigned int>
 * create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                                 const Tensor<1, dim>   direction,
 *                                 const unsigned int     level)
 * {
 *   std::vector<typename DoFHandler<dim>::level_cell_iterator> ordered_cells;
 *   ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level));
 *   for (const auto &cell : dof_handler.cell_iterators_on_level(level))
 *     ordered_cells.push_back(cell);
 *
 *   const DoFRenumbering::
 *     CompareDownstream<typename DoFHandler<dim>::level_cell_iterator, dim>
 *       comparator(direction);
 *   std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);
 *
 *   std::vector<unsigned> ordered_indices;
 *   ordered_indices.reserve(dof_handler.get_triangulation().n_cells(level));
 *
 *   for (const auto &cell : ordered_cells)
 *     ordered_indices.push_back(cell->index());
 *
 *   return ordered_indices;
 * }
 *
 *
 *
 * template <int dim>
 * std::vector<unsigned int>
 * create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                                 const Tensor<1, dim>   direction)
 * {
 *   std::vector<typename DoFHandler<dim>::active_cell_iterator> ordered_cells;
 *   ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells());
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     ordered_cells.push_back(cell);
 *
 *   const DoFRenumbering::
 *     CompareDownstream<typename DoFHandler<dim>::active_cell_iterator, dim>
 *       comparator(direction);
 *   std::sort(ordered_cells.begin(), ordered_cells.end(), comparator);
 *
 *   std::vector<unsigned int> ordered_indices;
 *   ordered_indices.reserve(dof_handler.get_triangulation().n_active_cells());
 *
 *   for (const auto &cell : ordered_cells)
 *     ordered_indices.push_back(cell->index());
 *
 *   return ordered_indices;
 * }
 *
 *
 * @endcode
 *
 * 产生随机排序的函数在精神上是相似的，它们首先将所有单元的信息放入一个数组。但是，它们不是对它们进行排序，而是利用C++提供的生成随机数的设施对元素进行随机洗牌。这样做的方式是在数组的所有元素上进行迭代，为之前的另一个元素抽取一个随机数，然后交换这些元素。其结果是对数组中的元素进行随机洗牌。
 *
 * @code
 * template <int dim>
 * std::vector<unsigned int>
 * create_random_cell_ordering(const DoFHandler<dim> &dof_handler,
 *                             const unsigned int     level)
 * {
 *   std::vector<unsigned int> ordered_cells;
 *   ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level));
 *   for (const auto &cell : dof_handler.cell_iterators_on_level(level))
 *     ordered_cells.push_back(cell->index());
 *
 *   std::mt19937 random_number_generator;
 *   std::shuffle(ordered_cells.begin(),
 *                ordered_cells.end(),
 *                random_number_generator);
 *
 *   return ordered_cells;
 * }
 *
 *
 *
 * template <int dim>
 * std::vector<unsigned int>
 * create_random_cell_ordering(const DoFHandler<dim> &dof_handler)
 * {
 *   std::vector<unsigned int> ordered_cells;
 *   ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells());
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     ordered_cells.push_back(cell->index());
 *
 *   std::mt19937 random_number_generator;
 *   std::shuffle(ordered_cells.begin(),
 *                ordered_cells.end(),
 *                random_number_generator);
 *
 *   return ordered_cells;
 * }
 *
 *
 * @endcode
 *
 * <a name="Righthandsideandboundaryvalues"></a> <h3>Right-hand side and
 * boundary values</h3>。
 *
 * 本教程中所解决的问题是对<a
 * href="https://global.oup.com/academic/product/finite-elements-and-fast-iterative-solvers-9780199678808">
 * Finite Elements and Fast Iterative Solvers: with Applications in
 * Incompressible Fluid Dynamics by Elman, Silvester, and
 * Wathen</a>第118页上的例3.1.3的改编。主要的区别是我们在域的中心增加了一个洞，其边界条件为零的Dirichlet。
 * 为了完整的描述，我们需要首先实现零右手边的类（当然，我们可以直接使用
 * Functions::ZeroFunction): ）。
 *
 * @code
 * template <int dim>
 * class RightHandSide : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 *
 *   virtual void value_list(const std::vector<Point<dim>> &points,
 *                           std::vector<double> &          values,
 *                           const unsigned int component = 0) const override;
 * };
 *
 *
 *
 * template <int dim>
 * double RightHandSide<dim>::value(const Point<dim> &,
 *                                  const unsigned int component) const
 * {
 *   Assert(component == 0, ExcIndexRange(component, 0, 1));
 *   (void)component;
 *
 *   return 0.0;
 * }
 *
 *
 *
 * template <int dim>
 * void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                     std::vector<double> &          values,
 *                                     const unsigned int component) const
 * {
 *   Assert(values.size() == points.size(),
 *          ExcDimensionMismatch(values.size(), points.size()));
 *
 *   for (unsigned int i = 0; i < points.size(); ++i)
 *     values[i] = RightHandSide<dim>::value(points[i], component);
 * }
 *
 *
 * @endcode
 * 我们也有迪里希特的边界条件。在外部正方形边界的连接部分，我们将数值设置为1，其他地方（包括内部圆形边界）的数值设置为0。
 *
 * @code
 * template <int dim>
 * class BoundaryValues : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 *
 *   virtual void value_list(const std::vector<Point<dim>> &points,
 *                           std::vector<double> &          values,
 *                           const unsigned int component = 0) const override;
 * };
 *
 *
 *
 * template <int dim>
 * double BoundaryValues<dim>::value(const Point<dim> & p,
 *                                   const unsigned int component) const
 * {
 *   Assert(component == 0, ExcIndexRange(component, 0, 1));
 *   (void)component;
 *
 * @endcode
 *
 * 如果  $x=1$  ，或如果  $x>0.5$  和  $y=-1$  ，则将边界设置为
 * 1 。
 *
 * @code
 *   if (std::fabs(p[0]
 *
 * - 1) < 1e-8 ||
 *       (std::fabs(p[1] + 1) < 1e-8 && p[0] >= 0.5))
 *     {
 *       return 1.0;
 *     }
 *   else
 *     {
 *       return 0.0;
 *     }
 * }
 *
 *
 *
 * template <int dim>
 * void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                      std::vector<double> &          values,
 *                                      const unsigned int component) const
 * {
 *   Assert(values.size() == points.size(),
 *          ExcDimensionMismatch(values.size(), points.size()));
 *
 *   for (unsigned int i = 0; i < points.size(); ++i)
 *     values[i] = BoundaryValues<dim>::value(points[i], component);
 * }
 *
 *
 *
 * @endcode
 *
 * <a name="Streamlinediffusionimplementation"></a> <h3>Streamline diffusion
 * implementation</h3>。
 *
 * 流水线扩散方法有一个稳定常数，我们需要能够计算出来。这个参数的计算方式的选择取自于<a
 * href="https://link.springer.com/chapter/10.1007/978-3-540-34288-5_27">On
 * Discontinuity-Capturing Methods for Convection-Diffusion Equations by
 * Volker John and Petr Knobloch</a>。
 *
 * @code
 * template <int dim>
 * double compute_stabilization_delta(const double         hk,
 *                                    const double         eps,
 *                                    const Tensor<1, dim> dir,
 *                                    const double         pk)
 * {
 *   const double Peclet = dir.norm() hk / (2.0 eps pk);
 *   const double coth =
 *     (1.0 + std::exp(-2.0 Peclet)) / (1.0
 *
 * - std::exp(-2.0 Peclet));
 *
 *   return hk / (2.0 dir.norm() pk) (coth
 *
 * - 1.0 / Peclet);
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProlemcodeclass"></a>
 * <h3><code>AdvectionProlem</code> class</h3>。
 *
 * 这是程序的主类，看起来应该与  step-16
 * 非常相似。主要的区别是，由于我们是在运行时定义我们的多网格平滑器，我们选择定义一个函数`create_smoother()`和一个类对象`mg_smoother`，这是一个
 * `std::unique_ptr`
 * 派生于MGSmoother的平滑器。请注意，对于从RelaxationBlock派生的平滑器，我们必须为每个级别包括一个`smoother_data`对象。这将包含关于单元排序和单元矩阵倒置方法的信息。
 *
 *
 * @code
 * template <int dim>
 * class AdvectionProblem
 * {
 * public:
 *   AdvectionProblem(const Settings &settings);
 *   void run();
 *
 * private:
 *   void setup_system();
 *
 *   template <class IteratorType>
 *   void assemble_cell(const IteratorType &cell,
 *                      ScratchData<dim> &  scratch_data,
 *                      CopyData &          copy_data);
 *   void assemble_system_and_multigrid();
 *
 *   void setup_smoother();
 *
 *   void solve();
 *   void refine_grid();
 *   void output_results(const unsigned int cycle) const;
 *
 *   Triangulation<dim> triangulation;
 *   DoFHandler<dim>    dof_handler;
 *
 *   const FE_Q<dim>     fe;
 *   const MappingQ<dim> mapping;
 *
 *   AffineConstraints<double> constraints;
 *
 *   SparsityPattern      sparsity_pattern;
 *   SparseMatrix<double> system_matrix;
 *
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 *
 *   MGLevelObject<SparsityPattern> mg_sparsity_patterns;
 *   MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns;
 *
 *   MGLevelObject<SparseMatrix<double>> mg_matrices;
 *   MGLevelObject<SparseMatrix<double>> mg_interface_in;
 *   MGLevelObject<SparseMatrix<double>> mg_interface_out;
 *
 *   mg::Matrix<Vector<double>> mg_matrix;
 *   mg::Matrix<Vector<double>> mg_interface_matrix_in;
 *   mg::Matrix<Vector<double>> mg_interface_matrix_out;
 *
 *   std::unique_ptr<MGSmoother<Vector<double>>> mg_smoother;
 *
 *   using SmootherType =
 *     RelaxationBlock<SparseMatrix<double>, double, Vector<double>>;
 *   using SmootherAdditionalDataType = SmootherType::AdditionalData;
 *   MGLevelObject<SmootherAdditionalDataType> smoother_data;
 *
 *   MGConstrainedDoFs mg_constrained_dofs;
 *
 *   Tensor<1, dim> advection_direction;
 *
 *   const Settings settings;
 * };
 *
 *
 *
 * template <int dim>
 * AdvectionProblem<dim>::AdvectionProblem(const Settings &settings)
 *   : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *   , dof_handler(triangulation)
 *   , fe(settings.fe_degree)
 *   , mapping(settings.fe_degree)
 *   , settings(settings)
 * {
 *   advection_direction[0] =
 *
 * -std::sin(numbers::PI / 6.0);
 *   if (dim >= 2)
 *     advection_direction[1] = std::cos(numbers::PI / 6.0);
 *   if (dim >= 3)
 *     AssertThrow(false, ExcNotImplemented());
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemsetup_systemcode"></a>
 * <h4><code>AdvectionProblem::setup_system()</code></h4>。
 *
 * 在这里我们首先为活动和多网格级别的网格设置DoFHandler、AffineConstraints和SparsityPattern对象。
 * 我们可以用DoFRenumbering类对活动的DoFs进行重新编号，但是平滑器只作用于多网格层，因此，这对计算并不重要。相反，我们将对每个多网格层次上的DoFs进行重新编号。
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::setup_system()
 * {
 *   const unsigned int n_levels = triangulation.n_levels();
 *
 *   dof_handler.distribute_dofs(fe);
 *
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 *
 *   constraints.clear();
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *
 *   VectorTools::interpolate_boundary_values(
 *     mapping, dof_handler, 0, BoundaryValues<dim>(), constraints);
 *   VectorTools::interpolate_boundary_values(
 *     mapping, dof_handler, 1, BoundaryValues<dim>(), constraints);
 *   constraints.close();
 *
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler,
 *                                   dsp,
 *                                   constraints,
 *                                   //keep_constrained_dofs ！改 =  false);
 *
 *   sparsity_pattern.copy_from(dsp);
 *   system_matrix.reinit(sparsity_pattern);
 *
 *   dof_handler.distribute_mg_dofs();
 *
 * @endcode
 *
 * 在列举了全局自由度以及（上面最后一行）水平自由度之后，让我们对水平自由度进行重新编号，以获得一个更好的平滑器，正如介绍中所解释的那样。
 * 如果需要的话，下面的第一个区块会对下游或上游方向的每个层次的自由度进行重新编号。这只对点平滑器（SOR和Jacobi）是必要的，因为块平滑器在单元上操作（见`create_smoother()`）。下面的块也实现了随机编号。
 *
 * @code
 *   if (settings.smoother_type == "SOR" || settings.smoother_type == "Jacobi")
 *     {
 *       if (settings.dof_renumbering ==
 *             Settings::DoFRenumberingStrategy::downstream ||
 *           settings.dof_renumbering ==
 *             Settings::DoFRenumberingStrategy::upstream)
 *         {
 *           const Tensor<1, dim> direction =
 *             (settings.dof_renumbering ==
 *                  Settings::DoFRenumberingStrategy::upstream ?
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
 * -1.0 :
 *                1.0)
 *             advection_direction;
 *
 *           for (unsigned int level = 0; level < n_levels; ++level)
 *             DoFRenumbering::downstream(dof_handler,
 *                                        level,
 *                                        direction,
 *                                        //dof_wise_renumbering = !改 true);
 *         }
 *       else if (settings.dof_renumbering ==
 *                Settings::DoFRenumberingStrategy::random)
 *         {
 *           for (unsigned int level = 0; level < n_levels; ++level)
 *             DoFRenumbering::random(dof_handler, level);
 *         }
 *       else
 *         Assert(false, ExcNotImplemented());
 *     }
 *
 * @endcode
 *
 * 该函数的其余部分只是设置了数据结构。下面代码的最后几行与其他GMG教程不同，因为它同时设置了接口输入和输出矩阵。我们需要这样做，因为我们的问题是非对称性的。
 *
 * @code
 *   mg_constrained_dofs.clear();
 *   mg_constrained_dofs.initialize(dof_handler);
 *
 *   mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, {0, 1});
 *
 *   mg_matrices.resize(0, n_levels
 *
 * - 1);
 *   mg_matrices.clear_elements();
 *   mg_interface_in.resize(0, n_levels
 *
 * - 1);
 *   mg_interface_in.clear_elements();
 *   mg_interface_out.resize(0, n_levels
 *
 * - 1);
 *   mg_interface_out.clear_elements();
 *   mg_sparsity_patterns.resize(0, n_levels
 *
 * - 1);
 *   mg_interface_sparsity_patterns.resize(0, n_levels
 *
 * - 1);
 *
 *   for (unsigned int level = 0; level < n_levels; ++level)
 *     {
 *       {
 *         DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
 *                                    dof_handler.n_dofs(level));
 *         MGTools::make_sparsity_pattern(dof_handler, dsp, level);
 *         mg_sparsity_patterns[level].copy_from(dsp);
 *         mg_matrices[level].reinit(mg_sparsity_patterns[level]);
 *       }
 *       {
 *         DynamicSparsityPattern dsp(dof_handler.n_dofs(level),
 *                                    dof_handler.n_dofs(level));
 *         MGTools::make_interface_sparsity_pattern(dof_handler,
 *                                                  mg_constrained_dofs,
 *                                                  dsp,
 *                                                  level);
 *         mg_interface_sparsity_patterns[level].copy_from(dsp);
 *
 *         mg_interface_in[level].reinit(mg_interface_sparsity_patterns[level]);
 *         mg_interface_out[level].reinit(mg_interface_sparsity_patterns[level]);
 *       }
 *     }
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemassemble_cellcode"></a>
 * <h4><code>AdvectionProblem::assemble_cell()</code></h4>。
 *
 * 这里我们定义了每个单元上的线性系统的装配，以便被下面的Mesh_loop()函数使用。这个函数为活动单元或水平单元(不管它的第一个参数是什么)装配单元矩阵，并且只有在调用活动单元时才装配右手边。
 *
 *
 * @code
 * template <int dim>
 * template <class IteratorType>
 * void AdvectionProblem<dim>::assemble_cell(const IteratorType &cell,
 *                                           ScratchData<dim> &  scratch_data,
 *                                           CopyData &          copy_data)
 * {
 *   copy_data.level = cell->level();
 *
 *   const unsigned int dofs_per_cell =
 *     scratch_data.fe_values.get_fe().n_dofs_per_cell();
 *   copy_data.dofs_per_cell = dofs_per_cell;
 *   copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 *
 *   const unsigned int n_q_points =
 *     scratch_data.fe_values.get_quadrature().size();
 *
 *   if (cell->is_level_cell() == false)
 *     copy_data.cell_rhs.reinit(dofs_per_cell);
 *
 *   copy_data.local_dof_indices.resize(dofs_per_cell);
 *   cell->get_active_or_mg_dof_indices(copy_data.local_dof_indices);
 *
 *   scratch_data.fe_values.reinit(cell);
 *
 *   RightHandSide<dim>  right_hand_side;
 *   std::vector<double> rhs_values(n_q_points);
 *
 *   right_hand_side.value_list(scratch_data.fe_values.get_quadrature_points(),
 *                              rhs_values);
 *
 * @endcode
 *
 * 如果我们使用流线扩散，我们必须将其贡献加入到单元矩阵和单元的右手边。如果我们不使用流线扩散，设置
 * $\delta=0$
 * 就可以否定下面的贡献，我们就可以使用标准的Galerkin有限元组合。
 *
 * @code
 *   const double delta = (settings.with_streamline_diffusion ?
 *                           compute_stabilization_delta(cell->diameter(),
 *                                                       settings.epsilon,
 *                                                       advection_direction,
 *                                                       settings.fe_degree) :
 *                           0.0);
 *
 *   for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *       {
 *         for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *           {
 * @endcode
 *
 * 局部矩阵的组装有两个部分。首先是Galerkin贡献。
 *
 * @code
 *             copy_data.cell_matrix(i, j) +=
 *               (settings.epsilon
 *                scratch_data.fe_values.shape_grad(i, q_point)
 *                scratch_data.fe_values.shape_grad(j, q_point)
 *                scratch_data.fe_values.JxW(q_point)) +
 *               (scratch_data.fe_values.shape_value(i, q_point)
 *                (advection_direction
 *                 scratch_data.fe_values.shape_grad(j, q_point))
 *                scratch_data.fe_values.JxW(q_point))
 * @endcode
 *  然后是流线扩散贡献。
 *
 * @code
 *               + delta
 *                   (advection_direction
 *                    scratch_data.fe_values.shape_grad(j, q_point))
 *                   (advection_direction
 *                    scratch_data.fe_values.shape_grad(i, q_point))
 *                   scratch_data.fe_values.JxW(q_point)
 *
 * -
 *               delta settings.epsilon
 *                 trace(scratch_data.fe_values.shape_hessian(j, q_point))
 *                 (advection_direction
 *                  scratch_data.fe_values.shape_grad(i, q_point))
 *                 scratch_data.fe_values.JxW(q_point);
 *           }
 *         if (cell->is_level_cell() == false)
 *           {
 * @endcode
 *  同样适用于右手边。首先是Galerkin贡献。
 *
 * @code
 *             copy_data.cell_rhs(i) +=
 *               scratch_data.fe_values.shape_value(i, q_point)
 *                 rhs_values[q_point] scratch_data.fe_values.JxW(q_point)
 * @endcode
 *  然后是流线扩散贡献。
 *
 * @code
 *               + delta rhs_values[q_point] advection_direction
 *                   scratch_data.fe_values.shape_grad(i, q_point)
 *                   scratch_data.fe_values.JxW(q_point);
 *           }
 *       }
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemassemble_system_and_multigridcode"></a>
 * <h4><code>AdvectionProblem::assemble_system_and_multigrid()</code></h4>。
 *
 * 这里我们采用 MeshWorker::mesh_loop()
 * 来翻阅单元格，为我们组装system_matrix、system_rhs和所有mg_matrices。
 *
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::assemble_system_and_multigrid()
 * {
 *   const auto cell_worker_active =
 *     [&](const decltype(dof_handler.begin_active()) &cell,
 *         ScratchData<dim> &                          scratch_data,
 *         CopyData &                                  copy_data) {
 *       this->assemble_cell(cell, scratch_data, copy_data);
 *     };
 *
 *   const auto copier_active = [&](const CopyData &copy_data) {
 *     constraints.distribute_local_to_global(copy_data.cell_matrix,
 *                                            copy_data.cell_rhs,
 *                                            copy_data.local_dof_indices,
 *                                            system_matrix,
 *                                            system_rhs);
 *   };
 *
 *
 *   MeshWorker::mesh_loop(dof_handler.begin_active(),
 *                         dof_handler.end(),
 *                         cell_worker_active,
 *                         copier_active,
 *                         ScratchData<dim>(fe, fe.degree + 1),
 *                         CopyData(),
 *                         MeshWorker::assemble_own_cells);
 *
 * @endcode
 *
 * 与活动层的约束不同，我们选择为每个多网格层在这个函数的本地创建约束对象，因为它们在程序的其他地方从来不需要。
 *
 * @code
 *   std::vector<AffineConstraints<double>> boundary_constraints(
 *     triangulation.n_global_levels());
 *   for (unsigned int level = 0; level < triangulation.n_global_levels();
 *        ++level)
 *     {
 *       IndexSet locally_owned_level_dof_indices;
 *       DoFTools::extract_locally_relevant_level_dofs(
 *         dof_handler, level, locally_owned_level_dof_indices);
 *       boundary_constraints[level].reinit(locally_owned_level_dof_indices);
 *       boundary_constraints[level].add_lines(
 *         mg_constrained_dofs.get_refinement_edge_indices(level));
 *       boundary_constraints[level].add_lines(
 *         mg_constrained_dofs.get_boundary_indices(level));
 *       boundary_constraints[level].close();
 *     }
 *
 *   const auto cell_worker_mg =
 *     [&](const decltype(dof_handler.begin_mg()) &cell,
 *         ScratchData<dim> &                      scratch_data,
 *         CopyData &                              copy_data) {
 *       this->assemble_cell(cell, scratch_data, copy_data);
 *     };
 *
 *   const auto copier_mg = [&](const CopyData &copy_data) {
 *     boundary_constraints[copy_data.level].distribute_local_to_global(
 *       copy_data.cell_matrix,
 *       copy_data.local_dof_indices,
 *       mg_matrices[copy_data.level]);
 *
 * @endcode
 *
 * 如果 $(i,j)$ 是一个 "interface_out "关系对，那么 $(j,i)$
 * 是一个 "interface_in "关系对。注意：对于
 * "interface_in"，我们加载接口条目的转置，即，dof对 $(j,i)$
 * 的条目被存储在
 * "interface_in(i,j)"。这是对对称情况的优化，允许在solve()中设置边缘矩阵时只使用一个矩阵。然而，在这里，由于我们的问题是非对称的，我们必须同时存储`interface_in`和`interface_out`矩阵。
 *
 * @code
 *     for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i)
 *       for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j)
 *         if (mg_constrained_dofs.is_interface_matrix_entry(
 *               copy_data.level,
 *               copy_data.local_dof_indices[i],
 *               copy_data.local_dof_indices[j]))
 *           {
 *             mg_interface_out[copy_data.level].add(
 *               copy_data.local_dof_indices[i],
 *               copy_data.local_dof_indices[j],
 *               copy_data.cell_matrix(i, j));
 *             mg_interface_in[copy_data.level].add(
 *               copy_data.local_dof_indices[i],
 *               copy_data.local_dof_indices[j],
 *               copy_data.cell_matrix(j, i));
 *           }
 *   };
 *
 *   MeshWorker::mesh_loop(dof_handler.begin_mg(),
 *                         dof_handler.end_mg(),
 *                         cell_worker_mg,
 *                         copier_mg,
 *                         ScratchData<dim>(fe, fe.degree + 1),
 *                         CopyData(),
 *                         MeshWorker::assemble_own_cells);
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemsetup_smoothercode"></a>
 * <h4><code>AdvectionProblem::setup_smoother()</code></h4>。
 *
 * 接下来，我们根据`.prm`文件中的设置来设置平滑器。两个重要的选项是多网格v周期每一级的平滑前和平滑后步骤的数量以及松弛参数。
 *
 *
 * 因为乘法往往比加法更强大，所以需要更少的平滑步骤来看到与网格大小无关的收敛。这一点对于块平滑器和点平滑器同样适用。这反映在对每一种平滑器的平滑步数的选择上，如下所示。
 *
 *
 * 点平滑器的松弛参数是根据试验和错误选择的，反映了在我们细化网格时保持GMRES求解的迭代次数不变（或尽可能接近）的必要值。在`.prm`文件中给
 * "Jacobi "和 "SOR
 * "的两个值是针对1度和3度有限元的。如果用户想改成其他度数，他们可能需要调整这些数字。对于块平滑器，这个参数有一个更直接的解释，即对于二维的加法，一个DoF可以有多达4个单元的重复贡献，因此我们必须将这些方法放松0.25来补偿。对于乘法来说，这不是一个问题，因为每个单元的逆向应用都会给其所有的DoF带来新的信息。
 *
 *
 * 最后，如上所述，点平滑器只对DoF进行操作，而块平滑器对单元进行操作，所以只有块平滑器需要被赋予有关单元排序的信息。点平滑器的DoF排序已经在`setup_system()`中得到了处理。
 *
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::setup_smoother()
 * {
 *   if (settings.smoother_type == "SOR")
 *     {
 *       using Smoother = PreconditionSOR<SparseMatrix<double>>;
 *
 *       auto smoother =
 *         std::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
 *                                                 Smoother,
 *                                                 Vector<double>>>();
 *       smoother->initialize(mg_matrices,
 *                            Smoother::AdditionalData(fe.degree == 1 ? 1.0 :
 *                                                                      0.62));
 *       smoother->set_steps(settings.smoothing_steps);
 *       mg_smoother = std::move(smoother);
 *     }
 *   else if (settings.smoother_type == "Jacobi")
 *     {
 *       using Smoother = PreconditionJacobi<SparseMatrix<double>>;
 *       auto smoother =
 *         std::make_unique<MGSmootherPrecondition<SparseMatrix<double>,
 *                                                 Smoother,
 *                                                 Vector<double>>>();
 *       smoother->initialize(mg_matrices,
 *                            Smoother::AdditionalData(fe.degree == 1 ? 0.6667 :
 *                                                                      0.47));
 *       smoother->set_steps(settings.smoothing_steps);
 *       mg_smoother = std::move(smoother);
 *     }
 *   else if (settings.smoother_type == "block SOR" ||
 *            settings.smoother_type == "block Jacobi")
 *     {
 *       smoother_data.resize(0, triangulation.n_levels()
 *
 * - 1);
 *
 *       for (unsigned int level = 0; level < triangulation.n_levels(); ++level)
 *         {
 *           DoFTools::make_cell_patches(smoother_data[level].block_list,
 *                                       dof_handler,
 *                                       level);
 *
 *           smoother_data[level].relaxation =
 *             (settings.smoother_type == "block SOR" ? 1.0 : 0.25);
 *           smoother_data[level].inversion = PreconditionBlockBase<double>::svd;
 *
 *           std::vector<unsigned int> ordered_indices;
 *           switch (settings.dof_renumbering)
 *             {
 *               case Settings::DoFRenumberingStrategy::downstream:
 *                 ordered_indices =
 *                   create_downstream_cell_ordering(dof_handler,
 *                                                   advection_direction,
 *                                                   level);
 *                 break;
 *
 *               case Settings::DoFRenumberingStrategy::upstream:
 *                 ordered_indices =
 *                   create_downstream_cell_ordering(dof_handler,
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
 * -1.0 advection_direction,
 *                                                   level);
 *                 break;
 *
 *               case Settings::DoFRenumberingStrategy::random:
 *                 ordered_indices =
 *                   create_random_cell_ordering(dof_handler, level);
 *                 break;
 *
 *               case Settings::DoFRenumberingStrategy::none:
 *                 break;
 *
 *               default:
 *                 AssertThrow(false, ExcNotImplemented());
 *                 break;
 *             }
 *
 *           smoother_data[level].order =
 *             std::vector<std::vector<unsigned int>>(1, ordered_indices);
 *         }
 *
 *       if (settings.smoother_type == "block SOR")
 *         {
 *           auto smoother = std::make_unique<MGSmootherPrecondition<
 *             SparseMatrix<double>,
 *             RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>,
 *             Vector<double>>>();
 *           smoother->initialize(mg_matrices, smoother_data);
 *           smoother->set_steps(settings.smoothing_steps);
 *           mg_smoother = std::move(smoother);
 *         }
 *       else if (settings.smoother_type == "block Jacobi")
 *         {
 *           auto smoother = std::make_unique<
 *             MGSmootherPrecondition<SparseMatrix<double>,
 *                                    RelaxationBlockJacobi<SparseMatrix<double>,
 *                                                          double,
 *                                                          Vector<double>>,
 *                                    Vector<double>>>();
 *           smoother->initialize(mg_matrices, smoother_data);
 *           smoother->set_steps(settings.smoothing_steps);
 *           mg_smoother = std::move(smoother);
 *         }
 *     }
 *   else
 *     AssertThrow(false, ExcNotImplemented());
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemsolvecode"></a>
 * <h4><code>AdvectionProblem::solve()</code></h4>。
 *
 * 在我们解决这个系统之前，我们必须首先设置多网格预处理程序。这需要设置各级之间的转换，粗略的矩阵求解器和平滑器。这个设置几乎与
 * Step-16
 * 相同，主要区别在于上面定义的各种平滑器，以及由于我们的问题是非对称的，我们需要不同的界面边缘矩阵。实际上，在本教程中，这些接口矩阵是空的，因为我们只使用全局细化，因此没有细化边。然而，我们在这里仍然包括了这两个矩阵，因为如果我们简单地切换到自适应细化方法，程序仍然可以正常运行）。)
 *
 *
 * 最后要注意的是，由于我们的问题是非对称性的，我们必须使用适当的Krylov子空间方法。我们在这里选择使用GMRES，因为它能保证在每次迭代中减少残差。GMRES的主要缺点是，每次迭代，存储的临时向量的数量都会增加一个，而且还需要计算与之前存储的所有向量的标量积。这是很昂贵的。通过使用重启的GMRES方法可以放松这一要求，该方法对我们在任何时候需要存储的向量数量设置了上限（这里我们在50个临时向量后重启，即48次迭代）。这样做的缺点是我们失去了在整个迭代过程中收集的信息，因此我们可以看到收敛速度较慢。因此，在哪里重启是一个平衡内存消耗、CPU工作量和收敛速度的问题。然而，本教程的目标是通过使用强大的GMG预处理程序来实现非常低的迭代次数，所以我们选择了重启长度，使下面显示的所有结果在重启发生之前就能收敛，因此我们有一个标准的GMRES方法。如果用户有兴趣，deal.II中提供的另一种合适的方法是BiCGStab。
 *
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::solve()
 * {
 *   const unsigned int max_iters       = 200;
 *   const double       solve_tolerance = 1e-8 system_rhs.l2_norm();
 *   SolverControl      solver_control(max_iters, solve_tolerance, true, true);
 *   solver_control.enable_history_data();
 *
 *   using Transfer = MGTransferPrebuilt<Vector<double>>;
 *   Transfer mg_transfer(mg_constrained_dofs);
 *   mg_transfer.build(dof_handler);
 *
 *   FullMatrix<double> coarse_matrix;
 *   coarse_matrix.copy_from(mg_matrices[0]);
 *   MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver;
 *   coarse_grid_solver.initialize(coarse_matrix);
 *
 *   setup_smoother();
 *
 *   mg_matrix.initialize(mg_matrices);
 *   mg_interface_matrix_in.initialize(mg_interface_in);
 *   mg_interface_matrix_out.initialize(mg_interface_out);
 *
 *   Multigrid<Vector<double>> mg(
 *     mg_matrix, coarse_grid_solver, mg_transfer,mg_smoother,mg_smoother);
 *   mg.set_edge_matrices(mg_interface_matrix_out, mg_interface_matrix_in);
 *
 *   PreconditionMG<dim, Vector<double>, Transfer> preconditioner(dof_handler,
 *                                                                mg,
 *                                                                mg_transfer);
 *
 *   std::cout << "     Solving with GMRES to tol " << solve_tolerance << "..."
 *             << std::endl;
 *   SolverGMRES<Vector<double>> solver(
 *     solver_control, SolverGMRES<Vector<double>>::AdditionalData(50, true));
 *
 *   Timer time;
 *   time.start();
 *   solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *   time.stop();
 *
 *   std::cout << "          converged in " << solver_control.last_step()
 *             << " iterations"
 *             << " in " << time.last_wall_time() << " seconds " << std::endl;
 *
 *   constraints.distribute(solution);
 *
 *   mg_smoother.release();
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemoutput_resultscode"></a>
 * <h4><code>AdvectionProblem::output_results()</code></h4>。
 *
 * 最后感兴趣的函数产生图形输出。这里我们以.vtu格式输出解决方案和单元格排序。
 *
 *
 * 在函数的顶部，我们为每个单元格生成一个索引，以可视化平滑器所使用的排序。请注意，我们只对活动单元而不是平滑器实际使用的层次做这个工作。对于点平滑器，我们对DoFs而不是单元进行重新编号，所以这只是对现实中发生的情况的一种近似。最后，这个随机排序并不是我们实际使用的随机排序（见`create_smoother()`）。
 * 单元的（整数）排序然后被复制到一个（浮点）矢量，用于图形输出。
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
 * {
 *   const unsigned int n_active_cells = triangulation.n_active_cells();
 *   Vector<double>     cell_indices(n_active_cells);
 *   {
 *     std::vector<unsigned int> ordered_indices;
 *     switch (settings.dof_renumbering)
 *       {
 *         case Settings::DoFRenumberingStrategy::downstream:
 *           ordered_indices =
 *             create_downstream_cell_ordering(dof_handler, advection_direction);
 *           break;
 *
 *         case Settings::DoFRenumberingStrategy::upstream:
 *           ordered_indices =
 *             create_downstream_cell_ordering(dof_handler,
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
 * -1.0 advection_direction);
 *           break;
 *
 *         case Settings::DoFRenumberingStrategy::random:
 *           ordered_indices = create_random_cell_ordering(dof_handler);
 *           break;
 *
 *         case Settings::DoFRenumberingStrategy::none:
 *           ordered_indices.resize(n_active_cells);
 *           for (unsigned int i = 0; i < n_active_cells; ++i)
 *             ordered_indices[i] = i;
 *           break;
 *
 *         default:
 *           AssertThrow(false, ExcNotImplemented());
 *           break;
 *       }
 *
 *     for (unsigned int i = 0; i < n_active_cells; ++i)
 *       cell_indices(ordered_indices[i]) = static_cast<double>(i);
 *   }
 *
 * @endcode
 *
 * 然后，考虑到以前的教程程序，该函数的其余部分就很简单了。
 *
 * @code
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 *   data_out.add_data_vector(cell_indices, "cell_index");
 *   data_out.build_patches();
 *
 *   const std::string filename =
 *     "solution-" + Utilities::int_to_string(cycle) + ".vtu";
 *   std::ofstream output(filename.c_str());
 *   data_out.write_vtu(output);
 * }
 *
 *
 * @endcode
 *
 * <a name="codeAdvectionProblemruncode"></a>
 * <h4><code>AdvectionProblem::run()</code></h4>。
 *
 * 和大多数教程一样，这个函数创建/细化网格并调用上面定义的各种函数来设置、装配、求解和输出结果。
 *
 *
 * 在第0个周期中，我们在正方形 <code>[-1,1]^dim</code>
 * 上生成网格，以原点为中心有一个半径为3/10单位的洞。对于`manifold_id`等于1的对象（即与洞相邻的面），我们分配一个球形流形。
 *
 *
 * @code
 * template <int dim>
 * void AdvectionProblem<dim>::run()
 * {
 *   for (unsigned int cycle = 0; cycle < (settings.fe_degree == 1 ? 7 : 5);
 *        ++cycle)
 *     {
 *       std::cout << "  Cycle " << cycle << ':' << std::endl;
 *
 *       if (cycle == 0)
 *         {
 *           GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,
 *                                                           0.3,
 *                                                           1.0);
 *
 *           const SphericalManifold<dim> manifold_description(Point<dim>(0, 0));
 *           triangulation.set_manifold(1, manifold_description);
 *         }
 *
 *       triangulation.refine_global();
 *
 *       setup_system();
 *
 *       std::cout << "     Number of active cells:       "
 *                 << triangulation.n_active_cells() << " ("
 *                 << triangulation.n_levels() << " levels)" << std::endl;
 *       std::cout << "     Number of degrees of freedom: "
 *                 << dof_handler.n_dofs() << std::endl;
 *
 *       assemble_system_and_multigrid();
 *
 *       solve();
 *
 *       if (settings.output)
 *         output_results(cycle);
 *
 *       std::cout << std::endl;
 *     }
 * }
 * } // namespace Step63
 *
 *
 * @endcode
 *  <a name="Thecodemaincodefunction"></a> <h3>The <code>main</code>
 * function</h3>。
 * 最后，主函数和大多数教程一样。唯一有趣的一点是，我们要求用户传递一个`.prm`文件作为唯一的命令行参数。如果没有给出参数文件，程序将在屏幕上输出一个带有所有默认值的样本参数文件的内容，然后用户可以复制并粘贴到自己的`.prm`文件中。
 *
 *
 * @code
 * int main(int argc, charargv[])
 * {
 * try
 *   {
 *     Step63::Settings settings;
 *     settings.get_parameters((argc > 1) ? (argv[1]) : "");
 *
 *     Step63::AdvectionProblem<2> advection_problem_2d(settings);
 *     advection_problem_2d.run();
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
 * <a name="Results"></a><h1>Results</h1> 。   <a
 * name="GMRESIterationNumbers"></a><h3> GMRES Iteration Numbers </h3> 。
 * GMG的主要优势在于它是一种 $\mathcal{O}(n)$
 * 方法，也就是说，问题的复杂性随着问题的大小而线性增加。为了证明本教程中介绍的线性求解器实际上是
 * $\mathcal{O}(n)$
 * ，我们需要做的就是证明GMRES求解的迭代次数随着网格的细化而保持大致不变。
 * 以下各表给出了GMRES的迭代次数，以减少初始残差的系数
 * $10^8$
 * 。我们选择了足够数量的平滑步骤(基于该方法)来获得与网格大小无关的迭代数。从下面的表格可以看出，该方法确实是
 * $\mathcal{O}(n)$  。 <a name="DoFCellRenumbering"></a><h4> DoF/Cell
 * Renumbering </h4> 。
 *
 *  逐点平滑器（"Jacobi "和
 * "SOR"）是按照DoF在每一层的编号顺序应用的。我们可以使用DoFRenumbering命名空间来影响这一点。块平滑器的应用是基于我们在`setup_smoother()`中设置的顺序。我们可以直观地看到这个编号。下面的图片显示了在下游、随机和上游编号中的活动单元的编号（从左到右）。
 * <img
 * src="https://www.dealii.org/images/steps/developer/step-63-cell-order.png"
 * alt="">
 * 让我们从加法平滑器开始。下表显示了从GMRES获得收敛所需的迭代次数。
 * <table align="center" class="doxtable"> <tr> <th></th> <th></th> <th
 * colspan="1">$Q_1$</th> <th colspan="7">Smoother (smoothing steps)</th>
 * </tr> <tr> <th></th> <th></th> <th></th> <th colspan="3">Jacobi (6)</th>
 * <th></th> <th colspan="3">Block Jacobi (3)</th> </tr> <tr> <th></th>
 * <th></th> <th></th> <th colspan="3">Renumbering Strategy</th> <th></th> <th
 * colspan="3">Renumbering Strategy</th> </tr> <tr> <th>Cells</th> <th></th>
 * <th>DoFs</th> <th>Downstream</th> <th>Random</th> <th>Upstream</th>
 * <th></th> <th>Downstream</th> <th>Random</th> <th>Upstream</th> </tr> <tr>
 * <th>32</th> <th></th> <th>48</th> <td>3</th> <td>3</th> <td>3</th>
 * <th></th> <td>3</th> <td>3</th> <td>3</th> </tr> <tr> <th>128</th>
 * <th></th> <th>160</th> <td>6</th> <td>6</th> <td>6</th> <th></th>
 * <td>6</th> <td>6</th> <td>6</th> </tr> <tr> <th>512</th> <th></th>
 * <th>576</th> <td>11</th> <td>11</th> <td>11</th> <th></th> <td>9</th>
 * <td>9</th> <td>9</th> </tr> <tr> <th>2048</th> <th></th> <th>2176</th>
 * <td>15</th> <td>15</th> <td>15</th> <th></th> <td>13</th> <td>13</th>
 * <td>13</th> </tr> <tr> <th>8192</th> <th></th> <th>8448</th> <td>18</th>
 * <td>18</th> <td>18</th> <th></th> <td>15</th> <td>15</th> <td>15</th> </tr>
 * <tr> <th>32768</th> <th></th> <th>33280</th> <td>20</th> <td>20</th>
 * <td>20</th> <th></th> <td>16</th> <td>16</th> <td>16</th> </tr> <tr>
 * <th>131072</th> <th></th> <th>132096</th> <td>20</th> <td>20</th>
 * <td>20</th> <th></th> <td>16</th> <td>16</th> <td>16</th> </tr> </table>
 * 我们看到，对DoFs/单元进行重新编号对收敛速度没有影响。这是因为这些平滑器对每个DoF(点平滑器)或单元(块平滑器)进行独立计算，并将结果相加。由于我们可以将这些平滑器定义为矩阵之和的应用，而矩阵加法是交换性的，所以我们对不同组件的加法顺序不会影响最终结果。
 * 另一方面，乘法平滑器的情况则不同。 <table align="center"
 * class="doxtable"> <tr> <th></th> <th></th> <th colspan="1">$Q_1$</th> <th
 * colspan="7">Smoother (smoothing steps)</th> </tr> <tr> <th></th> <th></th>
 * <th></th> <th colspan="3">SOR (3)</th> <th></th> <th colspan="3">Block SOR
 * (1)</th> </tr> <tr> <th></th> <th></th> <th></th> <th
 * colspan="3">Renumbering Strategy</th> <th></th> <th colspan="3">Renumbering
 * Strategy</th> </tr> <tr> <th>Cells</th> <th></th> <th>DoFs</th>
 * <th>Downstream</th> <th>Random</th> <th>Upstream</th> <th></th>
 * <th>Downstream</th> <th>Random</th> <th>Upstream</th> </tr> <tr>
 * <th>32</th> <th></th> <th>48</th> <td>2</th> <td>2</th> <td>3</th>
 * <th></th> <td>2</th> <td>2</th> <td>3</th> </tr> <tr> <th>128</th>
 * <th></th> <th>160</th> <td>5</th> <td>5</th> <td>7</th> <th></th>
 * <td>5</th> <td>5</th> <td>7</th> </tr> <tr> <th>512</th> <th></th>
 * <th>576</th> <td>7</th> <td>9</th> <td>11</th> <th></th> <td>7</th>
 * <td>7</th> <td>12</th> </tr> <tr> <th>2048</th> <th></th> <th>2176</th>
 * <td>10</th> <td>12</th> <td>15</th> <th></th> <td>8</th> <td>10</th>
 * <td>17</th> </tr> <tr> <th>8192</th> <th></th> <th>8448</th> <td>11</th>
 * <td>15</th> <td>19</th> <th></th> <td>10</th> <td>11</th> <td>20</th> </tr>
 * <tr> <th>32768</th> <th></th> <th>33280</th> <td>12</th> <td>16</th>
 * <td>20</th> <th></th> <td>10</th> <td>12</th> <td>21</th> </tr> <tr>
 * <th>131072</th> <th></th> <th>132096</th> <td>12</th> <td>16</th>
 * <td>19</th> <th></th> <td>11</th> <td>12</th> <td>21</th> </tr> </table>
 * 在这里，我们可以通过在平流方向上对DoFs/单元进行重新编号来加快收敛，同样，如果我们在相反的方向上进行重新编号，也会减慢收敛。这是因为平流主导的问题有一个信息流的方向性（平流方向），如果对DoFs/单元进行正确的重新编号，乘法就能够抓住这个方向。
 * 然而，乘法的这一特点取决于 $\varepsilon$
 * 的值。当我们增加 $\varepsilon$
 * ，问题变得更加以扩散为主时，我们在网格上的信息传播更加均匀，在平流方向上的编号优势就会减弱。相反，在
 * $\varepsilon=0$
 * 的极端情况下（仅有平流），我们有一个一阶PDE和乘法，正确的重编号方法成为有效的求解器。正确的下游编号可能会导致只需要一次迭代的方法，因为信息可以从流入边界向下游传播，而没有相反方向的信息传输。然而，请注意，在
 * $\varepsilon=0$
 * 的情况下，必须对这种情况下的边界条件给予特别关注）。
 *
 *  <a name="Pointvsblocksmoothers"></a><h4> %Point vs. block smoothers
 * </h4>。
 * 我们将把结果限制在使用下游编号的运行上。下面是对
 * $Q_1$ 和 $Q_3$ 两个元素的所有四个平滑器的交叉比较。
 * <table align="center" class="doxtable"> <tr> <th></th> <td></th> <th
 * colspan="1">$Q_1$</th> <th colspan="4">Smoother (smoothing steps)</th>
 * <th></th> <th colspan="1">$Q_3$</th> <th colspan="4">Smoother (smoothing
 * steps)</th> </tr> <tr> <th colspan="1">Cells</th> <td></th> <th
 * colspan="1">DoFs</th> <th colspan="1">Jacobi (6)</th> <th colspan="1">Block
 * Jacobi (3)</th> <th colspan="1">SOR (3)</th> <th colspan="1">Block SOR
 * (1)</th> <th></th> <th colspan="1">DoFs</th> <th colspan="1">Jacobi
 * (6)</th> <th colspan="1">Block Jacobi (3)</th> <th colspan="1">SOR (3)</th>
 * <th colspan="1">Block SOR (1)</th> </tr> <tr> <th>32</th> <td></th>
 * <th>48</th> <td>3</th> <td>3</th> <td>2</th> <td>2</th> <td></th>
 * <th>336</th> <td>15</th> <td>14</th> <td>15</th> <td>6</th> </tr> <tr>
 * <th>128</th> <td></th> <th>160</th> <td>6</th> <td>6</th> <td>5</th>
 * <td>5</th> <td></th> <th>1248</th> <td>23</th> <td>18</th> <td>21</th>
 * <td>9</th> </tr> <tr> <th>512</th> <td></th> <th>576</th> <td>11</th>
 * <td>9</th> <td>7</th> <td>7</th> <td></th> <th>4800</th> <td>29</th>
 * <td>21</th> <td>28</th> <td>9</th> </tr> <tr> <th>2048</th> <td></th>
 * <th>2176</th> <td>15</th> <td>13</th> <td>10</th> <td>8</th> <td></th>
 * <th>18816</th> <td>33</th> <td>22</th> <td>32</th> <td>9</th> </tr> <tr>
 * <th>8192</th> <td></th> <th>8448</th> <td>18</th> <td>15</th> <td>11</th>
 * <td>10</th> <td></th> <th>74496</th> <td>35</th> <td>22</th> <td>34</th>
 * <td>10</th> </tr> <tr> <th>32768</th> <td></th> <th>33280</th> <td>20</th>
 * <td>16</th> <td>12</th> <td>10</th> <td></th> </tr> <tr> <th>131072</th>
 * <td></th> <th>132096</th> <td>20</th> <td>16</th> <td>12</th> <td>11</th>
 * <td></th> </tr> </table> 我们看到，对于 $Q_1$
 * ，两个乘法平滑器所需要的平滑步骤和迭代次数的组合都比两个加法平滑器小。然而，当我们将度数增加到
 * $Q_3$
 * 元素时，在平滑步骤和迭代次数方面，块平滑器有明显的优势。具体来说，块状SOR平滑器在度数上给出了恒定的迭代次数，而块状Jacobi平滑器的迭代次数只增加了约38%，而Jacobi和SOR则分别增加了75%和183%。
 * <a name="Cost"></a><h3> Cost </h3> 。
 *
 * 迭代次数并不能完全说明一个单一平滑器对另一个平滑器的最优性。很明显，我们必须检查迭代的成本。块状平滑器在这里处于不利地位，因为它们必须为每个单元构建和反转一个单元矩阵。下面是一个具有74,496个DoF的
 * $Q_3$ 元素的求解时间的比较。 <table align="center"
 * class="doxtable"> <tr> <th colspan="1">$Q_3$</th> <th colspan="4">Smoother
 * (smoothing steps)</th> </tr> <tr> <th colspan="1">DoFs</th> <th
 * colspan="1">Jacobi (6)</th> <th colspan="1">Block Jacobi (3)</th> <th
 * colspan="1">SOR (3)</th> <th colspan="1">Block SOR (1)</th> </tr> <tr>
 * <th>74496</th> <td>0.68s</th> <td>5.82s</th> <td>1.18s</th> <td>1.02s</th>
 * </tr> </table>
 * 需要最多迭代的平滑器（Jacobi）实际上花费了最短的时间（大约是下一个快速方法的2/3）。这是因为应用Jacobismo的步骤所需要的只是乘以一个对角线矩阵，这是非常快的。另一方面，虽然SOR需要比块SOR多3倍的迭代（每个迭代有3倍的平滑步骤），但时间却完全相等，这意味着块SOR的平滑步骤比SOR的平滑步骤慢9倍。最后，Jacobi块的成本比SOR块高6倍，这从每个方法的1个步骤具有相同的成本（反转单元格矩阵并将它们相加或相乘），而Jacobi块每次迭代的平滑步骤数是2倍的事实中可以直观地理解。
 *
 *  <a name="Additionalpoints"></a><h3> Additional points </h3> 。
 * 还有几个重要的点需要提及。
 * <ol>   <li>  对于平行分布的网格，乘法不能在整个领域内执行。这是因为它们一次只操作一个单元，而下游的单元只有在上游的单元完成后才能被处理。这在单个处理器上是没有问题的。处理器只是一个接一个地浏览单元的列表。然而，在并行的情况下，这意味着一些处理器处于空闲状态，因为上游处理器还没有完成对当前处理器所拥有的上游单元的工作。一旦上游处理器完成工作，下游处理器就可以开始工作，但那时上游处理器已经没有工作了。换句话说，在这些平稳的步骤中，大多数时间，大多数处理器实际上是空闲的。这不是获得良好的并行可扩展性的方法!
 * 我们可以使用一种混合方法，在每个子域上应用乘法平滑器，但是当你增加子域的数量时，该方法接近于加法的行为。这是这些方法的一个主要缺点。
 * </li>
 * <li>
 * 目前对块平滑器的研究表明，很快我们就可以计算单元矩阵的逆值，比目前在交易里面做的便宜得多。这项研究是基于快速对角线化方法（可追溯到20世纪60年代），在光谱界已经使用了20年左右（见，例如，<a
 * href="https://doi.org/10.1007/s10915-004-4787-3"> Hybrid Multigrid/Schwarz
 * Algorithms for the Spectral Element Method by Lottes
 * and Fischer</a>）。目前，人们正在努力将这些方法推广到DG，并使其更加强大。此外，似乎应该能够利用无矩阵的实现，以及在域的内部，单元格矩阵往往看起来非常相似的事实，允许更少的矩阵逆计算。 </li>   </ol>  *结合1.
 * 结合1.和2.给出了一个很好的理由，期望类似块状雅可比的方法在未来变得非常强大，尽管目前对这些例子来说它是相当慢的。
 *
 *  <a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions
 * </h3> 。   <a name="ConstantiterationsforQsub5sub"></a><h4> Constant
 * iterations for Q<sub>5</sub> </h4> 。
 * 改变平滑步骤的数量和平滑器松弛参数（在
 * <code>Smoother::AdditionalData()</code> 中设置在
 * <code>create_smoother()</code>
 * 内，只对点平滑器有必要），以便我们对一个 $Q_5$
 * 元素保持恒定的迭代次数。 <a
 * name="Effectivenessofrenumberingforchangingepsilon"></a><h4> Effectiveness
 * of renumbering for changing epsilon </h4>。
 *
 *  增加/减少乘法的`.prm`文件中的参数
 * "Epsilon"，观察哪些数值的重新编号不再影响收敛速度。 <a
 * name="Meshadaptivity"></a><h4> Mesh adaptivity </h4>
 *
 * 该代码被设置为可以在自适应细化的网格中正常工作（界面矩阵被创建和设置）。设计一个合适的细化临界点或者尝试KellyErrorEstimator类。
 *
*<a name="PlainProg"></a><h1> The plain program</h1>  @include "step-63.cc"
 *
 */


