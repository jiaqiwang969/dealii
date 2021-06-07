//include/deal.II-translator/A-tutorial/step-75_0.txt
/**
 * @page step_75 The step-75 tutorial program
 * 本教程依赖于  step-27  ,  step-37  ,  step-40  。 @htmlonly <table
 * class="tutorial" width="50%"> <tr><th colspan="2"><b><small>Table of
 * contents</small></b><b><small>Table of contents</small></b></th></tr>
 * <tr><td width="50%" valign="top">
 * <ol>
 * <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro"
 * class=bold>Introduction</a>
 * <ul>
 * <li><a href="#Loadbalancing">Load balancing</a><a
 * href="#Loadbalancing">Load balancing</a>
 * <li><a href="#hpdecisionindicators">hp-decision indicators</a><a
 * href="#hpdecisionindicators">hp-decision indicators</a>
 * <li><a href="#Hybridgeometricmultigrid">Hybrid geometric multigrid</a><a
 * href="#Hybridgeometricmultigrid">Hybrid geometric multigrid</a>
 * <li><a href="#Thetestcase">The test case</a><a href="#Thetestcase">The test
 * case</a>
 * </ul>
 * <li> <a href="#CommProg" class=bold>The commented program</a><a
 * href="#CommProg" class=bold>The commented program</a>
 * <ul>
 * <li><a href="#Includefiles">Include files</a><a
 * href="#Includefiles">Include files</a>
 * <li><a href="#ThecodeSolutioncodeclasstemplate">The <code>Solution</code>
 * class template</a><a href="#ThecodeSolutioncodeclasstemplate">The
 * <code>Solution</code> class template</a>
 * <li><a href="#Parameters">Parameters</a><a
 * href="#Parameters">Parameters</a>
 * <li><a href="#MatrixfreeLaplaceoperator">Matrix-free Laplace operator</a><a
 * href="#MatrixfreeLaplaceoperator">Matrix-free Laplace operator</a>
 * <li><a href="#Solverandpreconditioner">Solver and preconditioner</a><a
 * href="#Solverandpreconditioner">Solver and preconditioner</a>
 * <ul>
 * <li><a
 * href="#Conjugategradientsolverwithmultigridpreconditioner">Conjugate-gradient
 * solver with multigrid preconditioner</a><a
 * href="#Conjugategradientsolverwithmultigridpreconditioner">Conjugate-gradient
 * solver with multigrid preconditioner</a>
 * <li><a
 * href="#Hybridpolynomialgeometricglobalcoarseningmultigridpreconditioner">Hybrid
 * polynomial/geometric-global-coarsening multigrid preconditioner</a><a
 * href="#Hybridpolynomialgeometricglobalcoarseningmultigridpreconditioner">Hybrid
 * polynomial/geometric-global-coarsening multigrid preconditioner</a>
 * </ul>
 * <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The
 * <code>LaplaceProblem</code> class template</a><a
 * href="#ThecodeLaplaceProblemcodeclasstemplate">The
 * <code>LaplaceProblem</code> class template</a>
 * <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The
 * <code>LaplaceProblem</code> class implementation</a><a
 * href="#ThecodeLaplaceProblemcodeclassimplementation">The
 * <code>LaplaceProblem</code> class implementation</a>
 * <ul>
 * <li><a href="#Constructor">Constructor</a> ]<a
 * href="#Constructor">Constructor</a>
 * <li><a
 * href="#LaplaceProbleminitialize_grid">LaplaceProblem::initialize_grid</a><a
 * href="#LaplaceProbleminitialize_grid">LaplaceProblem::initialize_grid</a>
 * <li><a
 * href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a><a
 * href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
 * <li><a
 * href="#LaplaceProblemprint_diagnostics">LaplaceProblem::print_diagnostics</a><a
 * href="#LaplaceProblemprint_diagnostics">LaplaceProblem::print_diagnostics</a>
 * <li><a
 * href="#LaplaceProblemsolve_system">LaplaceProblem::solve_system</a><a
 * href="#LaplaceProblemsolve_system">LaplaceProblem::solve_system</a>
 * <li><a
 * href="#LaplaceProblemcompute_indicators">LaplaceProblem::compute_indicators</a><a
 * href="#LaplaceProblemcompute_indicators">LaplaceProblem::compute_indicators</a>
 * <li><a
 * href="#LaplaceProblemadapt_resolution">LaplaceProblem::adapt_resolution</a><a
 * href="#LaplaceProblemadapt_resolution">LaplaceProblem::adapt_resolution</a>
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
 * <ul>
 * <li><a href="#Differenthpdecisionstrategies">Different hp-decision
 * strategies</a><a href="#Differenthpdecisionstrategies">Different
 * hp-decision strategies</a>
 * <li><a href="#Solvewithmatrixbasedmethods">Solve with matrix-based
 * methods</a><a href="#Solvewithmatrixbasedmethods">Solve with matrix-based
 * methods</a>
 * <li><a href="#Multigridvariants">Multigrid variants</a><a
 * href="#Multigridvariants">Multigrid variants</a>
 * </ul>
 * </ul>
 * <li> <a href="#PlainProg" class=bold>The plain program</a><a
 * href="#PlainProg" class=bold>The plain program</a>
 * </ol> </td> </tr> </table>
 * @endhtmlonly
 * <br>
 * <i>This program was contributed by Marc Fehling, Peter Munch and Wolfgang
 * Bangerth.
 * <br>
 * This material is based upon work partly supported by the National Science
 * Foundation under Award No. DMS-1821210, EAR-1550901, and OAC-1835673. Any
 * opinions, findings, and conclusions or recommendations expressed in this
 * publication are those of the authors and do not necessarily reflect the
 * views of the National Science Foundation.
 * <br>
 * Peter Munch would like to thank Timo Heister, Martin Kronbichler, and Laura
 * Prieto Saavedra for many very interesting discussions. </i>
 *
 * @note
 * 作为这个程序的先决条件，你需要安装p4estlibrary和Trilinos库。在<a
 * href="../../readme.html"
 * target="body">README</a>文件中描述了deal.IItogether与这些附加库的安装。
 *
 *
 * <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1> 。
 *
 * 在有限元背景下，更多的自由度通常会产生更精确的解决方案，但也需要更多的计算工作。
 * 在以前的整个教程中，我们找到了通过将网格分辨率与解的复杂性进行局部调整来有效分布自由度的方法（自适应网格细化，
 * step-6
 * ）。如果我们不仅调整网格，而且还局部调整每个单元上的相关无限元素的多项式程度，这种方法就特别有效（hp-adaptation,
 * step-27  ）。
 * 此外，分配更多的进程来同时运行你的程序，有助于在更短的时间内解决计算的工作量。根据你的机器的硬件架构，你的程序必须为所有进程都能访问相同的内存（共享内存，
 * step-18
 * ），或者进程被托管在几个独立的节点上（分布式内存，
 * step-40 ）的情况做好准备。
 * 在高性能计算领域，内存访问被证明是当前超级计算机的瓶颈。我们可以通过MatrixFree方法（
 * step-37
 * ）在飞行中计算矩阵向量乘积的效果来避免完全存储矩阵。它们可用于几何多网格方法（
 * step-50
 * ），也可用于多项式多网格方法，以极大地加快方程组的求解速度。
 * 本教程结合所有这些特点，以最先进的方式介绍了如何解决一个简单的拉普拉斯问题：在具有分布式内存的机器上利用hp-adaptation和无矩阵混合多栅方法。
 *
 *  <a name="Loadbalancing"></a><h3>Load balancing</h3> 。
 *
 * 对于有限元的并行应用，我们将网格划分为子域（又称域分解），并将其分配给进程。正如在
 * step-40
 * 中所展示的，这种划分发生在deal.II的活动单元上。在那里，每个单元都有相同的有限元和相同的自由度分配，以及大约相同的工作负荷。为了平衡所有进程的工作量，我们必须平衡所有参与进程中的单元数量。
 * 有了hp自适应方法，情况就不一样了：有限元类型可能因单元而异，因此自由度的数量也不同。匹配单元的数量并不能产生平衡的工作量。在无矩阵的情况下，工作量可以假定为与每个过程的自由度数量成正比，因为在最好的情况下，只有源向量和目的向量需要被加载。
 * 我们可以通过给每个单元分配权重来平衡工作量，这些权重与自由度的数量成正比，并在所有进程之间平衡所有权重的总和。给每个单元分配单独的权重可以用我们以后要用的
 * parallel::CellWeights 类来实现。
 *
 *  <a name="hpdecisionindicators"></a><h3>hp-decision indicators</h3>。
 *
 * 使用hp-adaptive方法，我们不仅要决定哪些单元要细化或粗化，而且还要选择如何做：要么调整网格分辨率，要么调整有限元的多项式程度。
 * 我们将再次根据当前解决方案的（后验）计算误差估计值来决定哪些单元需要调整，例如，使用KellyErrorEstimator。我们将同样决定如何适应（事后）计算的平滑度估计：大的多项式度数在解决方案的平滑部分效果最好，而细的网格分辨率在不规则部分是有利的。在
 * step-27
 * 中，我们提出了基于傅里叶系数的衰减来计算光滑度估计值。让我们利用这个机会，提出另一种遵循同样思路的方法，但采用Legendre系数。
 * 我们将简要介绍这种新技术的思路，但为了简单起见，将其描述限于一维。假设 $u_\text{hp}(x)$ 是一个定义在单元 $K$ 上的有限元函数@f[
 * u_\text{hp}(x) = \sum c_i \varphi_i(x)
 * @f]，其中每个 $\varphi_i(x)$ 是一个形状函数。我们可以在Legendre多项式 $P_k$ 的基础上等效地表示 $u_\text{hp}(x)$ 为@f[
 * u_\text{hp}(x) = \sum l_k P_k(x).
 * @f]。我们的目标是在有限元系数 $c_i$ 和Legendre系数 $l_k$ 之间获得一个映射。我们将通过把问题写成 $L^2$ 对 $u_\text{hp}(x)$ 在Legendre基础上的投影来实现这一目标。每个系数 $l_k$ 都可以通过@f[
 * l_k = \int_K u_\text{hp}(x) P_k(x) dx.
 * @f]计算。根据结构，Legendre多项式在 $L^2$ 的内积下是正交的。此外，我们假设它们已经被归一化，所以它们的内积可以写成@f[
 * \int_K P_i(x) P_j(x) dx = \det(J_K) \, \delta_{ij} @f]，其中
 * $\delta_{ij}$ 是Kronecker delta， $J_K$ 是 $\hat{K}$ 到 $K$
 * 的映射的Jacobian，（在本教程中）假设它是常数（即，该映射必须是仿射）。
 * 因此，结合所有这些假设，在Legendre基础上表达
 * $u_\text{hp}(x)$ 的投影矩阵只是 $\det(J_K) \, \mathbb{I}$  。
 *
 * *也就是说， $\det(J_K)$ 乘以身份矩阵。让 $F_K$ 成为从 $K$ 到其参考单元 $\hat{K}$ 的映射。因此，投影系统中右侧的条目为：@f[
 * \int_K u_\text{hp}(x) P_k(x) dx
 * = \det(J_K) \int_\hat{K} u_\text{hp}(F_K(\hat{x})) P_k(F_K(\hat{x})) d\hat{x}.
 * @f]回顾 $u_\text{hp}(x)$ 的形状函数表示，我们可以将其写成 $\det(J_K) \, \mathbf{C} \, \mathbf{c}$ ，其中 $\mathbf{C}$  ]是具有条目@f[
 * \int_K P_i(x) \varphi_j(x) dx = \det(J_K) \int_{\hat{K}} P_i(F_K(\hat{x}))
 * \varphi_j(F_K(\hat{x})) d\hat{x} = \det(J_K) \int_{\hat{K}}
 * \hat{P}_i(\hat{x}) \hat{\varphi}_j(\hat{x}) d\hat{x} \dealcoloneq \det(J_K)
 * \, C_{ij}
 * @f]的基数变化矩阵，因此 $\mathbf{C}$ 的值可以写成 <em> 独立于 </em> 的 $K$ ，在转换为参考坐标后，将 $\det(J_K)$ 前置因数。因此，把这一切放在一起，投影问题可以写成@f[
 * \det(J_K) \, \mathbb{I} \, \mathbf{l} = \det(J_K) \, \mathbf{C} \,
 * \mathbf{c}
 * @f]，可以简单地改写成@f[
 * \mathbf{l} = \mathbf{C} \, \mathbf{c}. @f]。
 * 在这一点上，我们需要强调的是，大多数有限元应用都使用非结构化的网格，对于这些网格的映射几乎总是非线性的。换句话说：
 * $J_K$
 * 在整个单元中是常数的假设对于一般的网格来说是不正确的。因此，对
 * $l_k$
 * 的正确计算不仅要求我们为每一个单元计算相应的变换矩阵
 * $\mathbf{C}$ ，而且还要求我们在单元 $K$
 * 上定义一组类似Legendre的正交函数，该单元可能有一个任意的和非常复杂的几何形状。尤其是第二部分，计算成本非常高。目前FESeries变换类的实现是依靠有一个恒定的Jacobian所带来的简化来提高性能，因此只对仿射映射产生正确的结果。除此之外，这种情况对本教程来说并不构成问题，因为我们只使用方形的单元。
 * Eibner和Melenk @cite eibner2007hp 认为，当且仅当Legendre系数的绝对值随指数增加而呈指数衰减时，一个函数是解析的，即可以用幂级数表示 $k$ ：@f[
 * \exists C,\sigma > 0 : \quad \forall k \in \mathbb{N}_0 : \quad |l_k|
 * \leq C \exp\left(
 *
 * - \sigma k \right) .
 * @f]衰减率 $\sigma$ 可以解释为衡量该函数的平稳性。我们可以把它理解为对转换系数进行线性回归拟合的斜率：@f[
 * \ln(|l_k|) \sim \ln(C)
 *
 * - \sigma k . @f] 。 我们将对每个单元进行这种拟合 $K$
 * ，以获得有限元近似的平滑度的局部估计。衰减率
 * $\sigma_K$ 将作为hp适应的决策指标。对于具有多项式程度
 * $p$ 的单元 $K$ 上的有限元，计算 $k \leq (p+1)$
 * 的系数被证明是估计平滑度的合理选择。你可以在 @cite
 * fehling2020 中找到一个更详细的、与维度无关的描述。
 * 以上所有内容已经在 FESeries::Legendre 类和
 * SmoothnessEstimator::Legendre
 * 命名空间中实现。有了误差估计和平滑度指标，我们就剩下对单元格进行实际细化和粗化的标记。来自
 * parallel::distributed::GridRefinement 和 hp::Refinement
 * 命名空间的一些函数将在后面帮助我们完成这些工作。
 *
 *  <a name="Hybridgeometricmultigrid"></a><h3>Hybrid geometric multigrid</h3>
 * 。
 *
 * 有限元矩阵通常是非常稀疏的。此外，hp-adaptive方法对应于每行非零项数量变化很大的矩阵。一些最先进的预处理程序，如
 * step-40
 * 中使用的代数多重网格（AMG），在这些情况下表现不佳。
 * 因此，我们将依靠一个无矩阵的混合多网格预处理器。
 * step-50
 * 已经证明了几何多网格方法与MatrixFree框架结合时的优越性。在hp-adaptive
 * FEM上的应用需要一些额外的工作，因为一个单元的子女可能有不同的多项式度数。作为补救措施，我们首先对线性元素进行p松弛（类似于Mitchell
 * @cite mitchell2010hpmg
 * ），然后以常规方式进行h松弛。在最粗的层次上，我们应用代数多网格解算器。p-multigrid、h-multigrid和AMG的结合使得该求解器成为一个混合多网格求解器。
 * 我们将通过使用MGTransferGlobalCoarsening，在现有的全局粗化基础设施的帮助下，创建一个具有上述特殊要求的自定义混合多网格预处理器。
 *
 *  <a name="Thetestcase"></a><h3>The test case</h3>
 *
 *  对于椭圆方程来说，每个再入角通常会调用asingularity
 * @cite brenner2008
 * 。我们可以利用这种情况对我们的hp决策算法进行测试：在所有需要调整的单元上，我们倾向于在奇点附近采用精细的网格，而在其他地方则采用高的多项式程度。
 * 作为在这些条件下要解决的最简单的椭圆问题，我们选择了L型域中的拉普拉斯方程，其重角在坐标系的原点。
 * 为了能够确定实际误差，我们制造了一个有已知解的边界值问题。在上述领域中，拉普拉斯方程的解是，在极坐标中， $(r, \varphi)$  : @f[
 * u_\text{sol} = r^{2/3} \sin(2/3 \varphi).
 * @f] 。
 * 另见 @cite brenner2008 或 @cite mitchell2014hp
 * 。解决方案看起来如下。 <div style="text-align:center;"> <img
 * src="https://www.dealii.org/images/steps/developer/step-75.solution.svg"
 * alt="Analytic solution."> </div>
 *
 * The singularity becomes obvious by investigating the solution's gradient in
 * the vicinity of the reentrant corner, i.e., the origin
 * @f[
 * \left\| \nabla u_\text{sol} \right\|_{2} = 2/3 r^{-1/3} , \quad
 * \lim\limits_{r \rightarrow 0} \left\| \nabla u_\text{sol} \right\|_{2} =
 * \infty .
 * @f]
 *
 * As we know where the singularity will be located, we expect that our
 * hp-decision algorithm decides for a fine grid resolution in this particular
 * region, and high polynomial degree anywhere else.
 *
 * So let's see if that is actually the case, and how hp-adaptation performs
 * compared to pure h-adaptation. But first let us have a detailed look at the
 * actual code.
 *
 *
 * <a name="CommProg"></a> <h1> The commented program</h1>
 *
 *  <a name="Includefiles"></a> <h3>Include files</h3>
 *
 *
 *
 * The following include files have been used and discussed in previous
 * tutorial programs, especially in step-27 and step-40.
 *
 *
 *
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/base/mpi.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/timer.h>
 *
 * #include <deal.II/distributed/grid_refinement.h>
 * #include <deal.II/distributed/tria.h>
 *
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 *
 * #include <deal.II/grid/grid_generator.h>
 *
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_series.h>
 *
 * #include <deal.II/hp/fe_collection.h>
 * #include <deal.II/hp/refinement.h>
 *
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 *
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/smoothness_estimator.h>
 * #include <deal.II/numerics/vector_tools.h>
 *
 * #include <algorithm>
 * #include <fstream>
 * #include <iostream>
 *
 * @endcode
 *
 *
 * For load balancing we will assign individual weights on cells, and for that
 * we will use the class parallel::CellWeights.
 *
 *
 *
 * @code
 * #include <deal.II/distributed/cell_weights.h>
 *
 * @endcode
 *
 *
 * The solution function requires a transformation from Cartesian to polar
 * coordinates. The GeometricUtilities::Coordinates namespace provides the
 * necessary tools.
 *
 *
 *
 * @code
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/geometric_utilities.h>
 *
 * @endcode
 *
 *
 * The following include files will enable the MatrixFree functionality.
 *
 *
 *
 * @code
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * #include <deal.II/matrix_free/tools.h>
 *
 * @endcode
 *
 *
 * We will use LinearAlgebra::distributed::Vector for linear algebra
 * operations.
 *
 *
 *
 * @code
 * #include <deal.II/lac/la_parallel_vector.h>
 *
 * @endcode
 *
 *
 * We are left to include the files needed by the multigrid solver.
 *
 *
 *
 * @code
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_constrained_dofs.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_transfer_global_coarsening.h>
 * #include <deal.II/multigrid/multigrid.h>
 *
 * namespace Step75
 * {
 * using namespace dealii;
 *
 * @endcode
 *
 *
 *  <a name="ThecodeSolutioncodeclasstemplate"></a> <h3>The
 * <code>Solution</code> class template</h3>
 *
 *
 *
 * We have an analytic solution to the scenario at our disposal. We will use
 * this solution to impose boundary conditions for the numerical solution of
 * the problem. The formulation of the solution requires a transformation to
 * polar coordinates. To transform from Cartesian to spherical coordinates, we
 * will use a helper function from the GeometricUtilities::Coordinates
 * namespace. The first two coordinates of this transformation correspond to
 * polar coordinates in the x-y-plane.
 *
 *
 *
 * @code
 * template <int dim>
 * class Solution : public Function<dim>
 * {
 * public:
 *   Solution()
 *     : Function<dim>()
 *   {}
 *
 *   virtual double value(const Point<dim> &p,
 *                        const unsigned int  [4.x.0] ) const override
 *   {
 *     const std::array<double, dim> p_sphere =
 *       GeometricUtilities::Coordinates::to_spherical(p);
 *
 *     constexpr const double alpha = 2. / 3.;
 *     return std::pow(p_sphere[0], alpha) std::sin(alpha p_sphere[1]);
 *   }
 * };
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="Parameters"></a> <h3>Parameters</h3>
 *
 *
 *
 * For this tutorial, we will use a simplified set of parameters. It is also
 * possible to use a ParameterHandler class here, but to keep this tutorial
 * short we decided on using simple structs. The actual intention of all these
 * parameters will be described in the upcoming classes at their respective
 * location where they are used.
 *
 *  The following parameter set controls the coarse-grid solver, the
 * smoothers, and the inter-grid transfer scheme of the multigrid mechanism.
 * We populate it with default parameters.
 *
 *
 *
 * @code
 * struct MultigridParameters
 * {
 *   struct
 *   {
 *     std::string  type            = "cg_with_amg";
 *     unsigned int maxiter         = 10000;
 *     double       abstol          = 1e-20;
 *     double       reltol          = 1e-4;
 *     unsigned int smoother_sweeps = 1;
 *     unsigned int n_cycles        = 1;
 *     std::string  smoother_type   = "ILU";
 *   } coarse_solver;
 *
 *   struct
 *   {
 *     std::string  type                = "chebyshev";
 *     double       smoothing_range     = 20;
 *     unsigned int degree              = 5;
 *     unsigned int eig_cg_n_iterations = 20;
 *   } smoother;
 *
 *   struct
 *   {
 *     MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
 *       p_sequence = MGTransferGlobalCoarseningTools::
 *         PolynomialCoarseningSequenceType::decrease_by_one;
 *     bool perform_h_transfer = true;
 *   } transfer;
 * };
 *
 *
 *
 * @endcode
 *
 *
 * This is the general parameter struct for the problem class. You will find
 * this struct divided into several categories, including general runtime
 * parameters, level limits, refine and coarsen fractions, as well as
 * parameters for cell weighting. It also contains an instance of the above
 * struct for multigrid parameters which will be passed to the multigrid
 * algorithm.
 *
 *
 *
 * @code
 * struct Parameters
 * {
 *   unsigned int n_cycles         = 8;
 *   double       tolerance_factor = 1e-12;
 *
 *   MultigridParameters mg_data;
 *
 *   unsigned int min_h_level            = 5;
 *   unsigned int max_h_level            = 12;
 *   unsigned int min_p_degree           = 2;
 *   unsigned int max_p_degree           = 6;
 *   unsigned int max_p_level_difference = 1;
 *
 *   double refine_fraction    = 0.3;
 *   double coarsen_fraction   = 0.03;
 *   double p_refine_fraction  = 0.9;
 *   double p_coarsen_fraction = 0.9;
 *
 *   double weighting_factor   = 1e6;
 *   double weighting_exponent = 1.;
 * };
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="MatrixfreeLaplaceoperator"></a> <h3>Matrix-free Laplace
 * operator</h3>
 *
 *
 *
 * This is a matrix-free implementation of the Laplace operator that will
 * basically take over the part of the `assemble_system()` function from other
 * tutorials. The meaning of all member functions will be explained at their
 * definition later.
 *
 *  We will use the FEEvaluation class to evaluate the solution vector at the
 * quadrature points and to perform the integration. In contrast to other
 * tutorials, the template arguments `degree` is set to $-1$ and `number of
 * quadrature in 1D` to $0$. In this case, FEEvaluation selects dynamically
 * the correct polynomial degree and number of quadrature points. Here, we
 * introduce an alias to FEEvaluation with the correct template parameters so
 * that we do not have to worry about them later on.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * class LaplaceOperator : public Subscriptor
 * {
 * public:
 *   using VectorType = LinearAlgebra::distributed::Vector<number>;
 *
 *   using FECellIntegrator = FEEvaluation<dim,
 *
 * -1, 0, 1, number>;
 *
 *   LaplaceOperator() = default;
 *
 *   LaplaceOperator(const hp::MappingCollection<dim> &mapping,
 *                   const DoFHandler<dim> &           dof_handler,
 *                   const hp::QCollection<dim> &      quad,
 *                   const AffineConstraints<number> & constraints,
 *                   VectorType &                      system_rhs);
 *
 *   void reinit(const hp::MappingCollection<dim> &mapping,
 *               const DoFHandler<dim> &           dof_handler,
 *               const hp::QCollection<dim> &      quad,
 *               const AffineConstraints<number> & constraints,
 *               VectorType &                      system_rhs);
 *
 *   types::global_dof_index m() const;
 *
 *   number el(unsigned int, unsigned int) const;
 *
 *   void initialize_dof_vector(VectorType &vec) const;
 *
 *   void vmult(VectorType &dst, const VectorType &src) const;
 *
 *   void Tvmult(VectorType &dst, const VectorType &src) const;
 *
 *   const TrilinosWrappers::SparseMatrix &get_system_matrix() const;
 *
 *   void compute_inverse_diagonal(VectorType &diagonal) const;
 *
 * private:
 *   void do_cell_integral_local(FECellIntegrator &integrator) const;
 *
 *   void do_cell_integral_global(FECellIntegrator &integrator,
 *                                VectorType &      dst,
 *                                const VectorType &src) const;
 *
 *
 *   void do_cell_integral_range(
 *     const MatrixFree<dim, number> &              matrix_free,
 *     VectorType &                                 dst,
 *     const VectorType &                           src,
 *     const std::pair<unsigned int, unsigned int> &range) const;
 *
 *   MatrixFree<dim, number> matrix_free;
 *
 * @endcode
 *
 *
 * To solve the equation system on the coarsest level with an AMG
 * preconditioner, we need an actual system matrix on the coarsest level. For
 * this purpose, we provide a mechanism that optionally computes a matrix from
 * the matrix-free formulation, for which we introduce a dedicated
 * SparseMatrix object. In the default case, this matrix stays empty. Once
 * `get_system_matrix()` is called, this matrix is filled (lazy allocation).
 * Since this is a `const` function, we need the "mutable" keyword here. We
 * also need a the constraints object to build the matrix.
 *
 *
 *
 * @code
 *   AffineConstraints<number>              constraints;
 *   mutable TrilinosWrappers::SparseMatrix system_matrix;
 * };
 *
 *
 *
 * @endcode
 *
 *
 * The following section contains functions to initialize and reinitialize the
 * class. In particular, these functions initialize the internal MatrixFree
 * instance. For sake of simplicity, we also compute the system
 * right-hand-side vector.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * LaplaceOperator<dim, number>::LaplaceOperator(
 *   const hp::MappingCollection<dim> &mapping,
 *   const DoFHandler<dim> &           dof_handler,
 *   const hp::QCollection<dim> &      quad,
 *   const AffineConstraints<number> & constraints,
 *   VectorType &                      system_rhs)
 * {
 *   this->reinit(mapping, dof_handler, quad, constraints, system_rhs);
 * }
 *
 *
 *
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::reinit(
 *   const hp::MappingCollection<dim> &mapping,
 *   const DoFHandler<dim> &           dof_handler,
 *   const hp::QCollection<dim> &      quad,
 *   const AffineConstraints<number> & constraints,
 *   VectorType &                      system_rhs)
 * {
 * @endcode
 *
 *
 * Clear internal data structures (in the case that the operator is reused).
 *
 *
 *
 * @code
 *   this->system_matrix.clear();
 *
 * @endcode
 *
 *
 * Copy the constraints, since they might be needed for computation of the
 * system matrix later on.
 *
 *
 *
 * @code
 *   this->constraints.copy_from(constraints);
 *
 * @endcode
 *
 *
 * Set up MatrixFree. At the quadrature points, we only need to evaluate the
 * gradient of the solution and test with the gradient of the shape functions
 * so that we only need to set the flag `update_gradients`.
 *
 *
 *
 * @code
 *   typename MatrixFree<dim, number>::AdditionalData data;
 *   data.mapping_update_flags = update_gradients;
 *
 *   matrix_free.reinit(mapping, dof_handler, constraints, quad, data);
 *
 * @endcode
 *
 *
 * Compute the right-hand side vector. For this purpose, we set up a second
 * MatrixFree instance that uses a modified AffineConstraints not containing
 * the constraints due to Dirichlet-boundary conditions. This modified
 * operator is applied to a vector with only the Dirichlet values set. The
 * result is the negative right-hand-side vector.
 *
 *
 *
 * @code
 *   {
 *     AffineConstraints<number> constraints_without_dbc;
 *
 *     IndexSet locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(dof_handler,
 *                                             locally_relevant_dofs);
 *     constraints_without_dbc.reinit(locally_relevant_dofs);
 *
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             constraints_without_dbc);
 *     constraints_without_dbc.close();
 *
 *     VectorType b, x;
 *
 *     this->initialize_dof_vector(system_rhs);
 *
 *     MatrixFree<dim, number> matrix_free;
 *     matrix_free.reinit(
 *       mapping, dof_handler, constraints_without_dbc, quad, data);
 *
 *     matrix_free.initialize_dof_vector(b);
 *     matrix_free.initialize_dof_vector(x);
 *
 *     constraints.distribute(x);
 *
 *     matrix_free.cell_loop(&LaplaceOperator::do_cell_integral_range,
 *                           this,
 *                           b,
 *                           x);
 *
 *     constraints.set_zero(b);
 *
 *     system_rhs
 *
 * -= b;
 *   }
 * }
 *
 *
 *
 * @endcode
 *
 *
 * The following functions are implicitly needed by the multigrid algorithm,
 * including the smoothers.
 *
 *
 *
 * Since we do not have a matrix, query the DoFHandler for the number of
 * degrees of freedom.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * types::global_dof_index LaplaceOperator<dim, number>::m() const
 * {
 *   return matrix_free.get_dof_handler().n_dofs();
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Access a particular element in the matrix. This function is neither needed
 * nor implemented, however, is required to compile the program.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * number LaplaceOperator<dim, number>::el(unsigned int, unsigned int) const
 * {
 *   Assert(false, ExcNotImplemented());
 *   return 0;
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Initialize the given vector. We simply delegate the task to the MatrixFree
 * function with the same name.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void
 * LaplaceOperator<dim, number>::initialize_dof_vector(VectorType &vec) const
 * {
 *   matrix_free.initialize_dof_vector(vec);
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Perform an operator evaluation by looping with the help of MatrixFree over
 * all cells and evaluating the effect of the cell integrals (see also:
 * `do_cell_integral_local()` and `do_cell_integral_global()`).
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::vmult(VectorType &      dst,
 *                                          const VectorType &src) const
 * {
 *   this->matrix_free.cell_loop(
 *     &LaplaceOperator::do_cell_integral_range, this, dst, src, true);
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Perform the transposed operator evaluation. Since we are considering
 * symmetric "matrices", this function can simply delegate it task to vmult().
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::Tvmult(VectorType &      dst,
 *                                           const VectorType &src) const
 * {
 *   this->vmult(dst, src);
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Since we do not have a system matrix, we cannot loop over the the diagonal
 * entries of the matrix. Instead, we compute the diagonal by performing a
 * sequence of operator evaluations to unit basis vectors. For this purpose,
 * an optimized function from the MatrixFreeTools namespace is used. The
 * inversion is performed manually afterwards.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::compute_inverse_diagonal(
 *   VectorType &diagonal) const
 * {
 *   MatrixFreeTools::compute_diagonal(matrix_free,
 *                                     diagonal,
 *                                     &LaplaceOperator::do_cell_integral_local,
 *                                     this);
 *
 *   for (auto &i : diagonal)
 *     i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;
 * }
 *
 *
 *
 * @endcode
 *
 *
 * In the matrix-free context, no system matrix is set up during
 * initialization of this class. As a consequence, it has to be computed here
 * if it should be requested. Since the matrix is only computed in this
 * tutorial for linear elements (on the coarse grid), this is acceptable. The
 * matrix entries are obtained via sequence of operator evaluations. For this
 * purpose, the optimized function MatrixFreeTools::compute_matrix() is used.
 * The matrix will only be computed if it has not been set up yet (lazy
 * allocation).
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * const TrilinosWrappers::SparseMatrix &
 * LaplaceOperator<dim, number>::get_system_matrix() const
 * {
 *   if (system_matrix.m() == 0 && system_matrix.n() == 0)
 *     {
 *       const auto &dof_handler = this->matrix_free.get_dof_handler();
 *
 *       TrilinosWrappers::SparsityPattern dsp(
 *         dof_handler.locally_owned_dofs(),
 *         dof_handler.get_triangulation().get_communicator());
 *
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints);
 *
 *       dsp.compress();
 *       system_matrix.reinit(dsp);
 *
 *       MatrixFreeTools::compute_matrix(
 *         matrix_free,
 *         constraints,
 *         system_matrix,
 *         &LaplaceOperator::do_cell_integral_local,
 *         this);
 *     }
 *
 *   return this->system_matrix;
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Perform cell integral on a cell batch without gathering and scattering the
 * values. This function is needed for the MatrixFreeTools functions since
 * these functions operate directly on the buffers of FEEvaluation.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::do_cell_integral_local(
 *   FECellIntegrator &integrator) const
 * {
 *   integrator.evaluate(EvaluationFlags::gradients);
 *
 *   for (unsigned int q = 0; q < integrator.n_q_points; ++q)
 *     integrator.submit_gradient(integrator.get_gradient(q), q);
 *
 *   integrator.integrate(EvaluationFlags::gradients);
 * }
 *
 *
 *
 * @endcode
 *
 *
 * Same as above but with access to the global vectors.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::do_cell_integral_global(
 *   FECellIntegrator &integrator,
 *   VectorType &      dst,
 *   const VectorType &src) const
 * {
 *   integrator.gather_evaluate(src, EvaluationFlags::gradients);
 *
 *   for (unsigned int q = 0; q < integrator.n_q_points; ++q)
 *     integrator.submit_gradient(integrator.get_gradient(q), q);
 *
 *   integrator.integrate_scatter(EvaluationFlags::gradients, dst);
 * }
 *
 *
 *
 * @endcode
 *
 *
 * This function loops over all cell batches within a cell-batch range and
 * calls the above function.
 *
 *
 *
 * @code
 * template <int dim, typename number>
 * void LaplaceOperator<dim, number>::do_cell_integral_range(
 *   const MatrixFree<dim, number> &              matrix_free,
 *   VectorType &                                 dst,
 *   const VectorType &                           src,
 *   const std::pair<unsigned int, unsigned int> &range) const
 * {
 *   FECellIntegrator integrator(matrix_free, range);
 *
 *   for (unsigned cell = range.first; cell < range.second; ++cell)
 *     {
 *       integrator.reinit(cell);
 *
 *       do_cell_integral_global(integrator, dst, src);
 *     }
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="Solverandpreconditioner"></a> <h3>Solver and preconditioner</h3>
 *
 *
 *
 *  <a name="Conjugategradientsolverwithmultigridpreconditioner"></a>
 * <h4>Conjugate-gradient solver with multigrid preconditioner</h4>
 *
 *
 *
 * This function solves the equation system with a sequence of provided
 * multigrid objects. It is meant to be treated as general as possible, hence
 * the multitude of template parameters.
 *
 *
 *
 * @code
 * template <typename VectorType,
 *           int dim,
 *           typename SystemMatrixType,
 *           typename LevelMatrixType,
 *           typename MGTransferType>
 * static void
 * mg_solve(SolverControl &            solver_control,
 *          VectorType &               dst,
 *          const VectorType &         src,
 *          const MultigridParameters &mg_data,
 *          const DoFHandler<dim> &    dof,
 *          const SystemMatrixType &   fine_matrix,
 *          const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices,
 *          const MGTransferType &                                 mg_transfer)
 * {
 *   AssertThrow(mg_data.coarse_solver.type == "cg_with_amg",
 *               ExcNotImplemented());
 *   AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented());
 *
 *   const unsigned int min_level = mg_matrices.min_level();
 *   const unsigned int max_level = mg_matrices.max_level();
 *
 *   using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
 *   using SmootherType               = PreconditionChebyshev<LevelMatrixType,
 *                                              VectorType,
 *                                              SmootherPreconditionerType>;
 *   using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;
 *
 * @endcode
 *
 *
 * We initialize level operators and Chebyshev smoothers here.
 *
 *
 *
 * @code
 *   mg::Matrix<VectorType> mg_matrix(mg_matrices);
 *
 *   MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
 *     min_level, max_level);
 *
 *   for (unsigned int level = min_level; level <= max_level; level++)
 *     {
 *       smoother_data[level].preconditioner =
 *         std::make_shared<SmootherPreconditionerType>();
 *       mg_matrices[level]->compute_inverse_diagonal(
 *         smoother_data[level].preconditioner->get_vector());
 *       smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range;
 *       smoother_data[level].degree          = mg_data.smoother.degree;
 *       smoother_data[level].eig_cg_n_iterations =
 *         mg_data.smoother.eig_cg_n_iterations;
 *     }
 *
 *   MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType>
 *     mg_smoother;
 *   mg_smoother.initialize(mg_matrices, smoother_data);
 *
 * @endcode
 *
 *
 * Next, we initialize the coarse-grid solver. We use conjugate-gradient
 * method with AMG as preconditioner.
 *
 *
 *
 * @code
 *   ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
 *                                               mg_data.coarse_solver.abstol,
 *                                               mg_data.coarse_solver.reltol,
 *                                               false,
 *                                               false);
 *   SolverCG<VectorType> coarse_grid_solver(coarse_grid_solver_control);
 *
 *   std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;
 *
 *   TrilinosWrappers::PreconditionAMG                 precondition_amg;
 *   TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
 *   amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
 *   amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
 *   amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();
 *
 *   precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(),
 *                               amg_data);
 *
 *   mg_coarse =
 *     std::make_unique<MGCoarseGridIterativeSolver<VectorType,
 *                                                  SolverCG<VectorType>,
 *                                                  LevelMatrixType,
 *                                                  decltype(precondition_amg)>>(
 *       coarse_grid_solver,mg_matrices[min_level], precondition_amg);
 *
 * @endcode
 *
 *
 * Finally, we create the Multigrid object, convert it to a preconditioner,
 * and use it inside of a conjugate-gradient solver to solve the linear system
 * of equations.
 *
 *
 *
 * @code
 *   Multigrid<VectorType> mg(
 *     mg_matrix,mg_coarse, mg_transfer, mg_smoother, mg_smoother);
 *
 *   PreconditionerType preconditioner(dof, mg, mg_transfer);
 *
 *   SolverCG<VectorType>(solver_control)
 *     .solve(fine_matrix, dst, src, preconditioner);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a
 * name="Hybridpolynomialgeometricglobalcoarseningmultigridpreconditioner"></a>
 * <h4>Hybrid polynomial/geometric-global-coarsening multigrid
 * preconditioner</h4>
 *
 *
 *
 * The above function deals with the actual solution for a given sequence of
 * multigrid objects. This functions creates the actual multigrid levels, in
 * particular the operators, and the transfer operator as a
 * MGTransferGlobalCoarsening object.
 *
 *
 *
 * @code
 * template <typename VectorType, typename OperatorType, int dim>
 * void solve_with_gmg(SolverControl &                  solver_control,
 *                     const OperatorType &             system_matrix,
 *                     VectorType &                     dst,
 *                     const VectorType &               src,
 *                     const MultigridParameters &      mg_data,
 *                     const hp::MappingCollection<dim> mapping_collection,
 *                     const DoFHandler<dim> &          dof_handler,
 *                     const hp::QCollection<dim> &     quadrature_collection)
 * {
 * @endcode
 *
 *
 * Create a DoFHandler and operator for each multigrid level, as well as,
 * create transfer operators. To be able to set up the operators, we need a
 * set of DoFHandler that we create via global coarsening of p or h. For
 * latter, we need also a sequence of Triangulation objects that are obtained
 * by Triangulation::coarsen_global().
 *
 *  In case no h-transfer is requested, we provide an empty deleter for the
 * `emplace_back()` function, since the Triangulation of our DoFHandler is an
 * external field and its destructor is called somewhere else.
 *
 *
 *
 * @code
 *   MGLevelObject<DoFHandler<dim>>                     dof_handlers;
 *   MGLevelObject<std::unique_ptr<OperatorType>>       operators;
 *   MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;
 *
 *   std::vector<std::shared_ptr<const Triangulation<dim>>>
 *     coarse_grid_triangulations;
 *   if (mg_data.transfer.perform_h_transfer)
 *     coarse_grid_triangulations =
 *       MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
 *         dof_handler.get_triangulation());
 *   else
 *     coarse_grid_triangulations.emplace_back(
 *       const_cast<Triangulation<dim>>(&(dof_handler.get_triangulation())),
 *       [](auto &) {});
 *
 * @endcode
 *
 *
 * Determine the total number of levels for the multigrid operation and
 * allocate sufficient memory for all levels.
 *
 *
 *
 * @code
 *   const unsigned int n_h_levels = coarse_grid_triangulations.size()
 *
 * - 1;
 *
 *   const auto get_max_active_fe_degree = [&](const auto &dof_handler) {
 *     unsigned int max = 0;
 *
 *     for (auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         max =
 *           std::max(max, dof_handler.get_fe(cell->active_fe_index()).degree);
 *
 *     return Utilities::MPI::max(max, MPI_COMM_WORLD);
 *   };
 *
 *   const unsigned int n_p_levels =
 *     MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
 *       get_max_active_fe_degree(dof_handler), mg_data.transfer.p_sequence)
 *       .size();
 *
 *   std::map<unsigned int, unsigned int> fe_index_for_degree;
 *   for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
 *     {
 *       const unsigned int degree = dof_handler.get_fe(i).degree;
 *       Assert(fe_index_for_degree.find(degree) == fe_index_for_degree.end(),
 *              ExcMessage("FECollection does not contain unique degrees."));
 *       fe_index_for_degree[degree] = i;
 *     }
 *
 *   unsigned int minlevel   = 0;
 *   unsigned int minlevel_p = n_h_levels;
 *   unsigned int maxlevel   = n_h_levels + n_p_levels
 *
 * - 1;
 *
 *   dof_handlers.resize(minlevel, maxlevel);
 *   operators.resize(minlevel, maxlevel);
 *   transfers.resize(minlevel, maxlevel);
 *
 * @endcode
 *
 *
 * Loop from the minimum (coarsest) to the maximum (finest) level and set up
 * DoFHandler accordingly. We start with the h-levels, where we distribute on
 * increasingly finer meshes linear elements.
 *
 *
 *
 * @code
 *   for (unsigned int l = 0; l < n_h_levels; ++l)
 *     {
 *       dof_handlers[l].reinit(*coarse_grid_triangulations[l]);
 *       dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
 *     }
 *
 * @endcode
 *
 *
 * After we reached the finest mesh, we will adjust the polynomial degrees on
 * each level. We reverse iterate over our data structure and start at the
 * finest mesh that contains all information about the active FE indices. We
 * then lower the polynomial degree of each cell level by level.
 *
 *
 *
 * @code
 *   for (unsigned int i = 0, l = maxlevel; i < n_p_levels; ++i,
 *
 * --l)
 *     {
 *       dof_handlers[l].reinit(dof_handler.get_triangulation());
 *
 *       if (l == maxlevel) // finest level
 *         {
 *           auto &dof_handler_mg = dof_handlers[l];
 *
 *           auto cell_other = dof_handler.begin_active();
 *           for (auto &cell : dof_handler_mg.active_cell_iterators())
 *             {
 *               if (cell->is_locally_owned())
 *                 cell->set_active_fe_index(cell_other->active_fe_index());
 *               cell_other++;
 *             }
 *         }
 *       else // coarse level
 *         {
 *           auto &dof_handler_fine   = dof_handlers[l + 1];
 *           auto &dof_handler_coarse = dof_handlers[l + 0];
 *
 *           auto cell_other = dof_handler_fine.begin_active();
 *           for (auto &cell : dof_handler_coarse.active_cell_iterators())
 *             {
 *               if (cell->is_locally_owned())
 *                 {
 *                   const unsigned int next_degree =
 *                     MGTransferGlobalCoarseningTools::
 *                       create_next_polynomial_coarsening_degree(
 *                         cell_other->get_fe().degree,
 *                         mg_data.transfer.p_sequence);
 *                   Assert(fe_index_for_degree.find(next_degree) !=
 *                            fe_index_for_degree.end(),
 *                          ExcMessage("Next polynomial degree in sequence "
 *                                     "does not exist in FECollection."));
 *
 *                   cell->set_active_fe_index(fe_index_for_degree[next_degree]);
 *                 }
 *               cell_other++;
 *             }
 *         }
 *
 *       dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
 *     }
 *
 * @endcode
 *
 *
 * Next, we will create all data structures additionally needed on each
 * multigrid level. This involves determining constraints with homogeneous
 * Dirichlet boundary conditions, and building the operator just like on the
 * active level.
 *
 *
 *
 * @code
 *   MGLevelObject<AffineConstraints<typename VectorType::value_type>>
 *     constraints(minlevel, maxlevel);
 *
 *   for (unsigned int level = minlevel; level <= maxlevel; ++level)
 *     {
 *       const auto &dof_handler = dof_handlers[level];
 *       auto &      constraint  = constraints[level];
 *
 *       IndexSet locally_relevant_dofs;
 *       DoFTools::extract_locally_relevant_dofs(dof_handler,
 *                                               locally_relevant_dofs);
 *       constraint.reinit(locally_relevant_dofs);
 *
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraint);
 *       VectorTools::interpolate_boundary_values(mapping_collection,
 *                                                dof_handler,
 *                                                0,
 *                                                Functions::ZeroFunction<dim>(),
 *                                                constraint);
 *       constraint.close();
 *
 *       VectorType dummy;
 *
 *       operators[level] = std::make_unique<OperatorType>(mapping_collection,
 *                                                         dof_handler,
 *                                                         quadrature_collection,
 *                                                         constraint,
 *                                                         dummy);
 *     }
 *
 * @endcode
 *
 *
 * Set up intergrid operators and collect transfer operators within a single
 * operator as needed by the Multigrid solver class.
 *
 *
 *
 * @code
 *   for (unsigned int level = minlevel; level < minlevel_p; ++level)
 *     transfers[level + 1].reinit_geometric_transfer(dof_handlers[level + 1],
 *                                                    dof_handlers[level],
 *                                                    constraints[level + 1],
 *                                                    constraints[level]);
 *
 *   for (unsigned int level = minlevel_p; level < maxlevel; ++level)
 *     transfers[level + 1].reinit_polynomial_transfer(dof_handlers[level + 1],
 *                                                     dof_handlers[level],
 *                                                     constraints[level + 1],
 *                                                     constraints[level]);
 *
 *   MGTransferGlobalCoarsening<dim, VectorType> transfer(
 *     transfers, [&](const auto l, auto &vec) {
 *       operators[l]->initialize_dof_vector(vec);
 *     });
 *
 * @endcode
 *
 *
 * Finally, proceed to solve the problem with multigrid.
 *
 *
 *
 * @code
 *   mg_solve(solver_control,
 *            dst,
 *            src,
 *            mg_data,
 *            dof_handler,
 *            system_matrix,
 *            operators,
 *            transfer);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> <h3>The
 * <code>LaplaceProblem</code> class template</h3>
 *
 *
 *
 * Now we will finally declare the main class of this program, which solves
 * the Laplace equation on subsequently refined function spaces. Its structure
 * will look familiar as it is similar to the main classes of step-27 and
 * step-40. There are basically just two additions:
 *
 *
 *
 * - The SparseMatrix object that would hold the system matrix has been
 * replaced by an object of the LaplaceOperator class for the MatrixFree
 * formulation.
 *
 *
 *
 * - An object of parallel::CellWeights, which will help us with load
 * balancing, has been added.
 *
 *
 *
 * @code
 * template <int dim>
 * class LaplaceProblem
 * {
 * public:
 *   LaplaceProblem(const Parameters &parameters);
 *
 *   void run();
 *
 * private:
 *   void initialize_grid();
 *   void setup_system();
 *   void print_diagnostics();
 *   void solve_system();
 *   void compute_indicators();
 *   void adapt_resolution();
 *   void output_results(const unsigned int cycle);
 *
 *   MPI_Comm mpi_communicator;
 *
 *   const Parameters prm;
 *
 *   parallel::distributed::Triangulation<dim> triangulation;
 *   DoFHandler<dim>                           dof_handler;
 *
 *   hp::MappingCollection<dim> mapping_collection;
 *   hp::FECollection<dim>      fe_collection;
 *   hp::QCollection<dim>       quadrature_collection;
 *   hp::QCollection<dim
 *
 * - 1>   face_quadrature_collection;
 *
 *   IndexSet locally_owned_dofs;
 *   IndexSet locally_relevant_dofs;
 *
 *   AffineConstraints<double> constraints;
 *
 *   LaplaceOperator<dim, double>               laplace_operator;
 *   LinearAlgebra::distributed::Vector<double> locally_relevant_solution;
 *   LinearAlgebra::distributed::Vector<double> system_rhs;
 *
 *   std::unique_ptr<FESeries::Legendre<dim>>    legendre;
 *   std::unique_ptr<parallel::CellWeights<dim>> cell_weights;
 *
 *   Vector<float> estimated_error_per_cell;
 *   Vector<float> hp_decision_indicators;
 *
 *   ConditionalOStream pcout;
 *   TimerOutput        computing_timer;
 * };
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="ThecodeLaplaceProblemcodeclassimplementation"></a> <h3>The
 * <code>LaplaceProblem</code> class implementation</h3>
 *
 *
 *
 *  <a name="Constructor"></a> <h4>Constructor</h4>
 *
 *
 *
 * The constructor starts with an initializer list that looks similar to the
 * one of step-40. We again prepare the ConditionalOStream object to allow
 * only the first process to output anything over the console, and initialize
 * the computing timer properly.
 *
 *
 *
 * @code
 * template <int dim>
 * LaplaceProblem<dim>::LaplaceProblem(const Parameters &parameters)
 *   : mpi_communicator(MPI_COMM_WORLD)
 *   , prm(parameters)
 *   , triangulation(mpi_communicator)
 *   , dof_handler(triangulation)
 *   , pcout(std::cout,
 *           (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *   , computing_timer(mpi_communicator,
 *                     pcout,
 *                     TimerOutput::summary,
 *                     TimerOutput::wall_times)
 * {
 *   Assert(prm.min_h_level <= prm.max_h_level,
 *          ExcMessage(
 *            "Triangulation level limits have been incorrectly set up."));
 *   Assert(prm.min_p_degree <= prm.max_p_degree,
 *          ExcMessage("FECollection degrees have been incorrectly set up."));
 *
 * @endcode
 *
 *
 * We need to prepare the data structures for the hp-functionality in the
 * actual body of the constructor, and create corresponding objects for every
 * degree in the specified range from the parameter struct. As we are only
 * dealing with non-distorted rectangular cells, a linear mapping object is
 * sufficient in this context.
 *
 *  In the Parameters struct, we provide ranges for levels on which the
 * function space is operating with a reasonable resolution. The multigrid
 * algorithm requires linear elements on the coarsest possible level. So we
 * start with the lowest polynomial degree and fill the collection with
 * consecutively higher degrees until the user-specified maximum is reached.
 *
 *
 *
 * @code
 *   mapping_collection.push_back(MappingQ1<dim>());
 *
 *   for (unsigned int degree = 1; degree <= prm.max_p_degree; ++degree)
 *     {
 *       fe_collection.push_back(FE_Q<dim>(degree));
 *       quadrature_collection.push_back(QGauss<dim>(degree + 1));
 *       face_quadrature_collection.push_back(QGauss<dim
 *
 * - 1>(degree + 1));
 *     }
 *
 * @endcode
 *
 *
 * As our FECollection contains more finite elements than we want to use for
 * the finite element approximation of our solution, we would like to limit
 * the range on which active FE indices can operate on. For this, the
 * FECollection class allows to register a hierarchy that determines the
 * succeeding and preceding finite element in case of of p-refinement and
 * p-coarsening, respectively. All functions in the hp::Refinement namespace
 * consult this hierarchy to determine future FE indices. We will register
 * such a hierarchy that only works on finite elements with polynomial degrees
 * in the proposed range <code>[min_p_degree, max_p_degree]</code>.
 *
 *
 *
 * @code
 *   const unsigned int min_fe_index = prm.min_p_degree
 *
 * - 1;
 *   fe_collection.set_hierarchy(
 *      [4.x.1]
 *     [](const typename hp::FECollection<dim> &fe_collection,
 *        const unsigned int                    fe_index)
 *
 * -> unsigned int {
 *       return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 :
 *                                                        fe_index;
 *     },
 *      [4.x.2]
 *     [min_fe_index](const typename hp::FECollection<dim> &,
 *                    const unsigned int fe_index)
 *
 * -> unsigned int {
 *       Assert(fe_index >= min_fe_index,
 *              ExcMessage("Finite element is not part of hierarchy!"));
 *       return (fe_index > min_fe_index) ? fe_index
 *
 * - 1 : fe_index;
 *     });
 *
 * @endcode
 *
 *
 * We initialize the FESeries::Legendre object in the default configuration
 * for smoothness estimation.
 *
 *
 *
 * @code
 *   legendre = std::make_unique<FESeries::Legendre<dim>>(
 *     SmoothnessEstimator::Legendre::default_fe_series(fe_collection));
 *
 * @endcode
 *
 *
 * The next part is going to be tricky. During execution of refinement, a few
 * hp-algorithms need to interfere with the actual refinement process on the
 * Triangulation object. We do this by connecting several functions to
 * Triangulation::Signals: signals will be called at different stages during
 * the actual refinement process and trigger all connected functions. We
 * require this functionality for load balancing and to limit the polynomial
 * degrees of neighboring cells.
 *
 *  For the former, we would like to assign a weight to every cell that is
 * proportional to the number of degrees of freedom of its future finite
 * element. The library offers a class parallel::CellWeights that allows to
 * easily attach individual weights at the right place during the refinement
 * process, i.e., after all refine and coarsen flags have been set correctly
 * for hp-adaptation and right before repartitioning for load balancing is
 * about to happen. Functions can be registered that will attach weights in
 * the form that $a (n_\text{dofs})^b$ with a provided pair of parameters
 * $(a,b)$. We register such a function in the following. Every cell will be
 * charged with a constant weight at creation, which is a value of 1000 (see
 * Triangulation::Signals::cell_weight).
 *
 *  For load balancing, efficient solvers like the one we use should scale
 * linearly with the number of degrees of freedom owned. Further, to increase
 * the impact of the weights we would like to attach, make sure that the
 * individual weight will exceed this base weight by orders of magnitude. We
 * set the parameters for cell weighting correspondingly: A large weighting
 * factor of $10^6$ and an exponent of $1$.
 *
 *
 *
 * @code
 *   cell_weights = std::make_unique<parallel::CellWeights<dim>>(
 *     dof_handler,
 *     parallel::CellWeights<dim>::ndofs_weighting(
 *       {prm.weighting_factor, prm.weighting_exponent}));
 *
 * @endcode
 *
 *
 * In h-adaptive applications, we ensure a 2:1 mesh balance by limiting the
 * difference of refinement levels of neighboring cells to one. With the
 * second call in the following code snippet, we will ensure the same for
 * p-levels on neighboring cells: levels of future finite elements are not
 * allowed to differ by more than a specified difference. The function
 * hp::Refinement::limit_p_level_difference takes care of this, but needs to
 * be connected to a very specific signal in the parallel context. The issue
 * is that we need to know how the mesh will be actually refined to set future
 * FE indices accordingly. As we ask the p4est oracle to perform refinement,
 * we need to ensure that the Triangulation has been updated with the
 * adaptation flags of the oracle first. An instantiation of
 * parallel::distributed::TemporarilyMatchRefineFlags does exactly that for
 * the duration of its life. Thus, we will create an object of this class
 * right before limiting the p-level difference, and connect the corresponding
 * lambda function to the signal
 * Triangulation::Signals::post_p4est_refinement, which will be triggered
 * after the oracle got refined, but before the Triangulation is refined.
 * Furthermore, we specify that this function will be connected to the front
 * of the signal, to ensure that the modification is performed before any
 * other function connected to the same signal.
 *
 *
 *
 * @code
 *   triangulation.signals.post_p4est_refinement.connect(
 *     [&, min_fe_index]() {
 *       const parallel::distributed::TemporarilyMatchRefineFlags<dim>
 *         refine_modifier(triangulation);
 *       hp::Refinement::limit_p_level_difference(dof_handler,
 *                                                prm.max_p_level_difference,
 *                                                 [4.x.3] min_fe_index);
 *     },
 *     boost::signals2::at_front);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProbleminitialize_grid"></a>
 * <h4>LaplaceProblem::initialize_grid</h4>
 *
 *
 *
 * For a L-shaped domain, we could use the function GridGenerator::hyper_L()
 * as demonstrated in step-50. However in the 2D case, that particular
 * function removes the first quadrant, while we need the fourth quadrant
 * removed in our scenario. Thus, we will use a different function
 * GridGenerator::subdivided_hyper_L() which gives us more options to create
 * the mesh. Furthermore, we formulate that function in a way that it also
 * generates a 3D mesh: the 2D L-shaped domain will basically elongated by 1
 * in the positive z-direction.
 *
 *  We first pretend to build a GridGenerator::subdivided_hyper_rectangle().
 * The parameters that we need to provide are Point objects for the lower left
 * and top right corners, as well as the number of repetitions that the base
 * mesh will have in each direction. We provide them for the first two
 * dimensions and treat the higher third dimension separately.
 *
 *  To create a L-shaped domain, we need to remove the excess cells. For this,
 * we specify the <code>cells_to_remove</code> accordingly. We would like to
 * remove one cell in every cell from the negative direction, but remove one
 * from the positive x-direction.
 *
 *  In the end, we supply the number of initial refinements that corresponds
 * to the supplied minimal grid refinement level. Further, we set the initial
 * active FE indices accordingly.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::initialize_grid()
 * {
 *   TimerOutput::Scope t(computing_timer, "initialize grid");
 *
 *   std::vector<unsigned int> repetitions(dim);
 *   Point<dim>                bottom_left, top_right;
 *   for (unsigned int d = 0; d < dim; ++d)
 *     if (d < 2)
 *       {
 *         repetitions[d] = 2;
 *         bottom_left[d] =
 *
 * -1.;
 *         top_right[d]   = 1.;
 *       }
 *     else
 *       {
 *         repetitions[d] = 1;
 *         bottom_left[d] = 0.;
 *         top_right[d]   = 1.;
 *       }
 *
 *   std::vector<int> cells_to_remove(dim, 1);
 *   cells_to_remove[0] =
 *
 * -1;
 *
 *   GridGenerator::subdivided_hyper_L(
 *     triangulation, repetitions, bottom_left, top_right, cells_to_remove);
 *
 *   triangulation.refine_global(prm.min_h_level);
 *
 *   const unsigned int min_fe_index = prm.min_p_degree
 *
 * - 1;
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     if (cell->is_locally_owned())
 *       cell->set_active_fe_index(min_fe_index);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemsetup_system"></a>
 * <h4>LaplaceProblem::setup_system</h4>
 *
 *
 *
 * This function looks exactly the same to the one of step-40, but you will
 * notice the absence of the system matrix as well as the scaffold that
 * surrounds it. Instead, we will initialize the MatrixFree formulation of the
 * <code>laplace_operator</code> here. For boundary conditions, we will use
 * the Solution class introduced earlier in this tutorial.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::setup_system()
 * {
 *   TimerOutput::Scope t(computing_timer, "setup system");
 *
 *   dof_handler.distribute_dofs(fe_collection);
 *
 *   locally_owned_dofs = dof_handler.locally_owned_dofs();
 *   DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *
 *   locally_relevant_solution.reinit(locally_owned_dofs,
 *                                    locally_relevant_dofs,
 *                                    mpi_communicator);
 *   system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *
 *   constraints.clear();
 *   constraints.reinit(locally_relevant_dofs);
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *   VectorTools::interpolate_boundary_values(
 *     mapping_collection, dof_handler, 0, Solution<dim>(), constraints);
 *   constraints.close();
 *
 *   laplace_operator.reinit(mapping_collection,
 *                           dof_handler,
 *                           quadrature_collection,
 *                           constraints,
 *                           system_rhs);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemprint_diagnostics"></a>
 * <h4>LaplaceProblem::print_diagnostics</h4>
 *
 *
 *
 * This is a function that prints additional diagnostics about the equation
 * system and its partitioning. In addition to the usual global number of
 * active cells and degrees of freedom, we also output their local
 * equivalents. For a regulated output, we will communicate the local
 * quantities with a Utilities::MPI::gather operation to the first process
 * which will then output all information. Output of local quantities is
 * limited to the first 8 processes to avoid cluttering the terminal.
 *
 *  Furthermore, we would like to print the frequencies of the polynomial
 * degrees in the numerical discretization. Since this information is only
 * stored locally, we will count the finite elements on locally owned cells
 * and later communicate them via Utilities::MPI::sum.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::print_diagnostics()
 * {
 *   const unsigned int first_n_processes =
 *     std::min<unsigned int>(8,
 *                            Utilities::MPI::n_mpi_processes(mpi_communicator));
 *   const bool output_cropped =
 *     first_n_processes < Utilities::MPI::n_mpi_processes(mpi_communicator);
 *
 *   {
 *     pcout << "   Number of active cells:       "
 *           << triangulation.n_global_active_cells() << std::endl
 *           << "     by partition:              ";
 *
 *     std::vector<unsigned int> n_active_cells_per_subdomain =
 *       Utilities::MPI::gather(mpi_communicator,
 *                              triangulation.n_locally_owned_active_cells());
 *     for (unsigned int i = 0; i < first_n_processes; ++i)
 *       pcout << ' ' << n_active_cells_per_subdomain[i];
 *     if (output_cropped)
 *       pcout << " ...";
 *     pcout << std::endl;
 *   }
 *
 *   {
 *     pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl
 *           << "     by partition:              ";
 *
 *     std::vector<types::global_dof_index> n_dofs_per_subdomain =
 *       Utilities::MPI::gather(mpi_communicator,
 *                              dof_handler.n_locally_owned_dofs());
 *     for (unsigned int i = 0; i < first_n_processes; ++i)
 *       pcout << ' ' << n_dofs_per_subdomain[i];
 *     if (output_cropped)
 *       pcout << " ...";
 *     pcout << std::endl;
 *   }
 *
 *   {
 *     std::vector<types::global_dof_index> n_constraints_per_subdomain =
 *       Utilities::MPI::gather(mpi_communicator, constraints.n_constraints());
 *
 *     pcout << "   Number of constraints:        "
 *           << std::accumulate(n_constraints_per_subdomain.begin(),
 *                              n_constraints_per_subdomain.end(),
 *                              0)
 *           << std::endl
 *           << "     by partition:              ";
 *     for (unsigned int i = 0; i < first_n_processes; ++i)
 *       pcout << ' ' << n_constraints_per_subdomain[i];
 *     if (output_cropped)
 *       pcout << " ...";
 *     pcout << std::endl;
 *   }
 *
 *   {
 *     std::vector<unsigned int> n_fe_indices(fe_collection.size(), 0);
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         n_fe_indices[cell->active_fe_index()]++;
 *
 *     Utilities::MPI::sum(n_fe_indices, mpi_communicator, n_fe_indices);
 *
 *     pcout << "   Frequencies of poly. degrees:";
 *     for (unsigned int i = 0; i < fe_collection.size(); ++i)
 *       if (n_fe_indices[i] > 0)
 *         pcout << ' ' << fe_collection[i].degree << ":" << n_fe_indices[i];
 *     pcout << std::endl;
 *   }
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemsolve_system"></a>
 * <h4>LaplaceProblem::solve_system</h4>
 *
 *
 *
 * The scaffold around the solution is similar to the one of step-40. We
 * prepare a vector that matches the requirements of MatrixFree and collect
 * the locally-relevant degrees of freedoms we solved the equation system. The
 * solution happens with the function introduced earlier.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::solve_system()
 * {
 *   TimerOutput::Scope t(computing_timer, "solve system");
 *
 *   LinearAlgebra::distributed::Vector<double> completely_distributed_solution;
 *   laplace_operator.initialize_dof_vector(completely_distributed_solution);
 *
 *   SolverControl solver_control(system_rhs.size(),
 *                                prm.tolerance_factor system_rhs.l2_norm());
 *
 *   solve_with_gmg(solver_control,
 *                  laplace_operator,
 *                  completely_distributed_solution,
 *                  system_rhs,
 *                  prm.mg_data,
 *                  mapping_collection,
 *                  dof_handler,
 *                  quadrature_collection);
 *
 *   pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *         << std::endl;
 *
 *   constraints.distribute(completely_distributed_solution);
 *
 *   locally_relevant_solution.copy_locally_owned_data_from(
 *     completely_distributed_solution);
 *   locally_relevant_solution.update_ghost_values();
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemcompute_indicators"></a>
 * <h4>LaplaceProblem::compute_indicators</h4>
 *
 *
 *
 * This function contains only a part of the typical <code>refine_grid</code>
 * function from other tutorials and is new in that sense. Here, we will only
 * calculate all indicators for adaptation with actually refining the grid. We
 * do this for the purpose of writing all indicators to the file system, so we
 * store them for later.
 *
 *  Since we are dealing the an elliptic problem, we will make use of the
 * KellyErrorEstimator again, but with a slight difference. Modifying the
 * scaling factor of the underlying face integrals to be dependent on the
 * actual polynomial degree of the neighboring elements is favorable in
 * hp-adaptive applications @cite davydov2017hp. We can do this by specifying
 * the very last parameter from the additional ones you notices. The others
 * are actually just the defaults.
 *
 *  For the purpose of hp-adaptation, we will calculate smoothness estimates
 * with the strategy presented in the tutorial introduction and use the
 * implementation in SmoothnessEstimator::Legendre. In the Parameters struct,
 * we set the minimal polynomial degree to 2 as it seems that the smoothness
 * estimation algorithms have trouble with linear elements.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::compute_indicators()
 * {
 *   TimerOutput::Scope t(computing_timer, "compute indicators");
 *
 *   estimated_error_per_cell.grow_or_shrink(triangulation.n_active_cells());
 *   KellyErrorEstimator<dim>::estimate(
 *     dof_handler,
 *     face_quadrature_collection,
 *     std::map<types::boundary_id, const Function<dim>>(),
 *     locally_relevant_solution,
 *     estimated_error_per_cell,
 *      [4.x.4] ComponentMask(),
 *      [4.x.5] nullptr,
 *      [4.x.6] numbers::invalid_unsigned_int,
 *      [4.x.7] numbers::invalid_subdomain_id,
 *      [4.x.8] numbers::invalid_material_id,
 *      [4.x.9]
 *     KellyErrorEstimator<dim>::Strategy::face_diameter_over_twice_max_degree);
 *
 *   hp_decision_indicators.grow_or_shrink(triangulation.n_active_cells());
 *   SmoothnessEstimator::Legendre::coefficient_decay(*legendre,
 *                                                    dof_handler,
 *                                                    locally_relevant_solution,
 *                                                    hp_decision_indicators);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemadapt_resolution"></a>
 * <h4>LaplaceProblem::adapt_resolution</h4>
 *
 *
 *
 * With the previously calculated indicators, we will finally flag all cells
 * for adaptation and also execute refinement in this function. As in previous
 * tutorials, we will use the "fixed number" strategy, but now for
 * hp-adaptation.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::adapt_resolution()
 * {
 *   TimerOutput::Scope t(computing_timer, "adapt resolution");
 *
 * @endcode
 *
 *
 * First, we will set refine and coarsen flags based on the error estimates on
 * each cell. There is nothing new here.
 *
 *  We will use general refine and coarsen fractions that have been elaborated
 * in the other deal.II tutorials: using the fixed number strategy, we will
 * flag 30% of all cells for refinement and 3% for coarsening, as provided in
 * the Parameters struct.
 *
 *
 *
 * @code
 *   parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
 *     triangulation,
 *     estimated_error_per_cell,
 *     prm.refine_fraction,
 *     prm.coarsen_fraction);
 *
 * @endcode
 *
 *
 * Next, we will make all adjustments for hp-adaptation. We want to refine and
 * coarsen those cells flagged in the previous step, but need to decide if we
 * would like to do it by adjusting the grid resolution or the polynomial
 * degree.
 *
 *  The next function call sets future FE indices according to the previously
 * calculated smoothness indicators as p-adaptation indicators. These indices
 * will only be set on those cells that have refine or coarsen flags assigned.
 *
 *  For the p-adaptation fractions, we will take an educated guess. Since we
 * only expect a single singularity in our scenario, i.e., in the origin of
 * the domain, and a smooth solution anywhere else, we would like to strongly
 * prefer to use p-adaptation over h-adaptation. This reflects in our choice
 * of a fraction of 90% for both p-refinement and p-coarsening.
 *
 *
 *
 * @code
 *   hp::Refinement::p_adaptivity_fixed_number(dof_handler,
 *                                             hp_decision_indicators,
 *                                             prm.p_refine_fraction,
 *                                             prm.p_coarsen_fraction);
 *
 * @endcode
 *
 *
 * At this stage, we have both the future FE indices and the classic refine
 * and coarsen flags set, from which the latter will be interpreted by
 * Triangulation::execute_coarsening_and_refinement() for h-adaptation. We
 * would like to only impose one type of adaptation on cells, which is what
 * the next function will sort out for us. In short, on cells which have both
 * types of indicators assigned, we will favor the p-adaptation one and remove
 * the h-adaptation one.
 *
 *
 *
 * @code
 *   hp::Refinement::choose_p_over_h(dof_handler);
 *
 * @endcode
 *
 *
 * After setting all indicators, we will remove those that exceed the
 * specified limits of the provided level ranges in the Parameters struct.
 * This limitation naturally arises for p-adaptation as the number of supplied
 * finite elements is limited. In addition, we registered a custom hierarchy
 * for p-adaptation in the constructor. Now, we need to do this manually in
 * the h-adaptive context like in step-31.
 *
 *  We will iterate over all cells on the designated min and max levels and
 * remove the corresponding flags. As an alternative, we could also flag these
 * cells for p-adaptation by setting future FE indices accordingly instead of
 * simply clearing the refine and coarsen flags.
 *
 *
 *
 * @code
 *   Assert(triangulation.n_levels() >= prm.min_h_level + 1 &&
 *            triangulation.n_levels() <= prm.max_h_level + 1,
 *          ExcInternalError());
 *
 *   if (triangulation.n_levels() > prm.max_h_level)
 *     for (const auto &cell :
 *          triangulation.active_cell_iterators_on_level(prm.max_h_level))
 *       cell->clear_refine_flag();
 *
 *   for (const auto &cell :
 *        triangulation.active_cell_iterators_on_level(prm.min_h_level))
 *     cell->clear_coarsen_flag();
 *
 * @endcode
 *
 *
 * In the end, we are left to execute coarsening and refinement. Here, not
 * only the grid will be updated, but also all previous future FE indices will
 * become active.
 *
 *  Remember that we have attached functions to triangulation signals in the
 * constructor, will be triggered in this function call. So there is even more
 * happening: weighted repartitioning will be performed to ensure load
 * balancing, as well as we will limit the difference of p-levels between
 * neighboring cells.
 *
 *
 *
 * @code
 *   triangulation.execute_coarsening_and_refinement();
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemoutput_results"></a>
 * <h4>LaplaceProblem::output_results</h4>
 *
 *
 *
 * Writing results to the file system in parallel applications works exactly
 * like in step-40. In addition to the data containers that we prepared
 * throughout the tutorial, we would also like to write out the polynomial
 * degree of each finite element on the grid as well as the subdomain each
 * cell belongs to. We prepare necessary containers for this in the scope of
 * this function.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::output_results(const unsigned int cycle)
 * {
 *   TimerOutput::Scope t(computing_timer, "output results");
 *
 *   Vector<float> fe_degrees(triangulation.n_active_cells());
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     if (cell->is_locally_owned())
 *       fe_degrees(cell->active_cell_index()) = cell->get_fe().degree;
 *
 *   Vector<float> subdomain(triangulation.n_active_cells());
 *   for (auto &subd : subdomain)
 *     subd = triangulation.locally_owned_subdomain();
 *
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(locally_relevant_solution, "solution");
 *   data_out.add_data_vector(fe_degrees, "fe_degree");
 *   data_out.add_data_vector(subdomain, "subdomain");
 *   data_out.add_data_vector(estimated_error_per_cell, "error");
 *   data_out.add_data_vector(hp_decision_indicators, "hp_indicator");
 *   data_out.build_patches(mapping_collection);
 *
 *   data_out.write_vtu_with_pvtu_record(
 *     "./", "solution", cycle, mpi_communicator, 2, 1);
 * }
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="LaplaceProblemrun"></a> <h4>LaplaceProblem::run</h4>
 *
 *
 *
 * The actual run function again looks very familiar to step-40. The only
 * addition is the bracketed section that precedes the actual cycle loop.
 * Here, we will pre-calculate the Legendre transformation matrices. In
 * general, these will be calculated on the fly via lazy allocation whenever a
 * certain matrix is needed. For timing purposes however, we would like to
 * calculate them all at once before the actual time measurement begins. We
 * will thus designate their calculation to their own scope.
 *
 *
 *
 * @code
 * template <int dim>
 * void LaplaceProblem<dim>::run()
 * {
 *   pcout << "Running with Trilinos on "
 *         << Utilities::MPI::n_mpi_processes(mpi_communicator)
 *         << " MPI rank(s)..." << std::endl;
 *
 *   {
 *     pcout << "Calculating transformation matrices..." << std::endl;
 *     TimerOutput::Scope t(computing_timer, "calculate transformation");
 *     legendre->precalculate_all_transformation_matrices();
 *   }
 *
 *   for (unsigned int cycle = 0; cycle < prm.n_cycles; ++cycle)
 *     {
 *       pcout << "Cycle " << cycle << ':' << std::endl;
 *
 *       if (cycle == 0)
 *         initialize_grid();
 *       else
 *         adapt_resolution();
 *
 *       setup_system();
 *
 *       print_diagnostics();
 *
 *       solve_system();
 *
 *       compute_indicators();
 *
 *       if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
 *         output_results(cycle);
 *
 *       computing_timer.print_summary();
 *       computing_timer.reset();
 *
 *       pcout << std::endl;
 *     }
 * }
 * } // namespace Step75
 *
 *
 *
 * @endcode
 *
 *
 *  <a name="main"></a> <h4>main()</h4>
 *
 *
 *
 * The final function is the <code>main</code> function that will ultimately
 * create and run a LaplaceOperator instantiation. Its structure is similar to
 * most other tutorial programs.
 *
 *
 *
 * @code
 * int main(int argc, charargv[])
 * {
 * try
 *   {
 *     using namespace dealii;
 *     using namespace Step75;
 *
 *     Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 *
 *     Parameters        prm;
 *     LaplaceProblem<2> laplace_problem(prm);
 *     laplace_problem.run();
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
 *
 * <a name="Results"></a><h1>Results</h1>
 *
 *
 * When you run the program with the given parameters on four processes in
 * release mode, your terminal output should look like this:
 *
 * @code
 * Running with Trilinos on 4 MPI rank(s)...
 * Calculating transformation matrices...
 * Cycle 0:
 * Number of active cells:       3072
 *   by partition:               768 768 768 768
 * Number of degrees of freedom: 12545
 *   by partition:               3201 3104 3136 3104
 * Number of constraints:        542
 *   by partition:               165 74 138 165
 * Frequencies of poly. degrees: 2:3072
 * Solved in 7 iterations.
 *
 *
 *
 *
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |     0.598s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | calculate transformation        |         1 |    0.0533s |       8.9% |
 * | compute indicators              |         1 |    0.0177s |         3% |
 * | initialize grid                 |         1 |    0.0397s |       6.6% |
 * | output results                  |         1 |    0.0844s |        14% |
 * | setup system                    |         1 |    0.0351s |       5.9% |
 * | solve system                    |         1 |     0.362s |        61% |
 * +---------------------------------+-----------+------------+------------+
 *
 *
 *
 *
 * Cycle 1:
 * Number of active cells:       3351
 *   by partition:               875 761 843 872
 * Number of degrees of freedom: 18223
 *   by partition:               4535 4735 4543 4410
 * Number of constraints:        1202
 *   by partition:               303 290 326 283
 * Frequencies of poly. degrees: 2:2523 3:828
 * Solved in 7 iterations.
 *
 *
 *
 *
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |     0.442s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | adapt resolution                |         1 |    0.0189s |       4.3% |
 * | compute indicators              |         1 |    0.0135s |         3% |
 * | output results                  |         1 |     0.064s |        14% |
 * | setup system                    |         1 |    0.0232s |       5.2% |
 * | solve system                    |         1 |     0.322s |        73% |
 * +---------------------------------+-----------+------------+------------+
 *
 *
 *
 *
 * ...
 *
 *
 *
 *
 * Cycle 7:
 * Number of active cells:       5610
 *   by partition:               1324 1483 1482 1321
 * Number of degrees of freedom: 82062
 *   by partition:               21116 19951 20113 20882
 * Number of constraints:        14383
 *   by partition:               3825 3225 3557 3776
 * Frequencies of poly. degrees: 2:1130 3:1283 4:2727 5:465 6:5
 * Solved in 7 iterations.
 *
 *
 *
 *
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |     0.932s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | adapt resolution                |         1 |    0.0182s |       1.9% |
 * | compute indicators              |         1 |    0.0173s |       1.9% |
 * | output results                  |         1 |    0.0572s |       6.1% |
 * | setup system                    |         1 |    0.0252s |       2.7% |
 * | solve system                    |         1 |     0.813s |        87% |
 * +---------------------------------+-----------+------------+------------+
 * @endcode
 *
 *
 * When running the code with more processes, you will notice slight
 * differences in the number of active cells and degrees of freedom. This is
 * due to the fact that solver and preconditioner depend on the partitioning
 * of the problem, which might yield to slight differences of the solution in
 * the last digits and ultimately yields to different adaptation behavior.
 *
 * Furthermore, the number of iterations for the solver stays about the same
 * in all cycles despite hp-adaptation, indicating the robustness of the
 * proposed algorithms and promising good scalability for even larger problem
 * sizes and on more processes.
 *
 * Let us have a look at the graphical output of the program. After all
 * refinement cycles in the given parameter configuration, the actual
 * discretized function space looks like the following with its partitioning
 * on twelve processes on the left and the polynomial degrees of finite
 * elements on the right. In the left picture, each color represents a unique
 * subdomain. In the right picture, the lightest color corresponds to the
 * polynomial degree two and the darkest one corresponds to degree six:
 *
 * <div class="twocolumn" style="width: 80%; text-align: center;"> <div> <img
 * src="https://www.dealii.org/images/steps/developer/step-75.subdomains-07.svg"
 * alt="Partitioning after seven refinements."> </div> <div> <img
 * src="https://www.dealii.org/images/steps/developer/step-75.fedegrees-07.svg"
 * alt="Local approximation degrees after seven refinements."> </div> </div>
 *
 *   <a name="extensions"></a><a
 * name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>
 * 。
 *
 *  <a name="Differenthpdecisionstrategies"></a><h4>Different hp-decision
 * strategies</h4>。
 *
 * deal.II库提供了多种策略来决定对单元的适应类型：要么调整网格分辨率，要么改变多项式程度。我们在本教程中只介绍了<i>Legendre
 * coefficient decay</i>策略，而 step-27
 * 则展示了相同想法的<i>Fourier</i>等价物。
 * 关于这些策略的概述，请参见 step-27 中的 "扩展的可能性
 * "一节，或相应的文档中的详细描述。
 * 在那里，提到了另一个到目前为止还没有在任何教程中展示过的策略：基于<i>refinement
 * history</i>的策略。这个方法在并行分布式应用中的使用比其他方法更棘手，所以我们将强调与之相伴的挑战。我们需要有关细化标志的最终状态的信息，并且我们需要在细化的网格之间转移解决方案。对于前者，我们需要将
 * hp::Refinement::predict_error() 函数附加到
 * Triangulation::Signals::post_p4est_refinement
 * 信号上，使其被称为<i>after</i>的
 * hp::Refinement::limit_p_level_difference()
 * 函数。在这个阶段，所有的精炼标志和未来的FE指数都被最终设定，并且可以对错误进行可靠的预测。预测的误差需要借助于
 * parallel::distributed::CellDataTransfer.
 * 在精炼网格之间进行转移。
 * 试着在本教程中实施这些策略之一，并观察结果的微妙变化。你会注意到，所有的策略都能识别出重心角附近的奇异点，并在这些区域进行
 * $h$ 精化，而在体域中则倾向于 $p$
 * 精化。这些策略的详细比较见 @cite fehling2020 。
 *
 *  <a name="Solvewithmatrixbasedmethods"></a><h4>Solve with matrix-based
 * methods</h4> 。
 *
 * 本教程只关注无矩阵策略。然而，所有的hp-adaptivealgorithms在并行分布式环境中也可以使用基于矩阵的方法。
 * 要创建一个系统矩阵，你可以使用
 * LaplaceOperator::get_system_matrix() 函数，或者使用类似于 step-27
 * 的 <code>assemble_system()</code>
 * 函数，然后你可以像往常一样将系统矩阵传递给求解器。
 * 你可以对基于矩阵和无矩阵实现的结果进行计时，量化加速，并说服自己哪种变量更快。
 *
 *  <a name="Multigridvariants"></a><h4>Multigrid variants</h4>
 *
 * 为了简单起见，我们将自己限制在单一类型的粗网格求解器（带AMG的CG）、平滑器（带point
 * Jacobi预处理的Chebyshev平滑器），以及多网格算法中的几何粗化方案（globalcoarsening）。请自由尝试其他方法，研究它们的性能和稳健性。
 *
* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-75.cc"
 *
 */


