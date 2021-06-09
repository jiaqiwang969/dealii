/**
@page step_59 The step-59 tutorial program
This tutorial depends on step-37.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#ThesymmetricinteriorpenaltyformulationfortheLaplacian">The symmetric interior penalty formulation for the Laplacian</a>
        <li><a href="#FaceintegrationsupportinMatrixFreeandFEFaceEvaluation">Face integration support in MatrixFree and FEFaceEvaluation</a>
        <li><a href="#TheFE_DGQHermiteelement">The FE_DGQHermite element</a>
        <li><a href="#AnapproximateblockJacobismootherusingthefastdiagonalizationmethod">An approximate block-Jacobi smoother using the fast diagonalization method</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Matrixfreeimplementation">Matrix-free implementation</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Comparisonofefficiencyatdifferentpolynomialdegrees">Comparison of efficiency at different polynomial degrees</a>
        <li><a href="#Evaluationofefficiencyofingredients">Evaluation of efficiency of ingredients</a>
        <li><a href="#Possibilitiesforextension">Possibilities for extension</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-59/doc/intro.dox

 <br> 

<i>
This program was contributed by Katharina Kormann and Martin
Kronbichler.


This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA). </i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


无矩阵算子评估使高阶多项式基的离散化得到了非常有效的实现，这是由于一种叫做和因子化的方法。这个概念已经在step-37和step-48教程程序中介绍过。在这个教程程序中，我们将这些概念扩展到包括面积分的非连续Galerkin（DG）方案，这是一类高阶特别普遍的方法。

无矩阵评估的基本思想与连续元素的评估是一样的。迭代求解器或多网格平滑器中出现的矩阵-向量乘积不是由经典的稀疏矩阵核来实现的，而是通过对基础积分的实时评估来隐含地应用。对于用张量积正交规则集成的张量积形状函数，通过使用和-因子化技术，这种评估特别有效，它将最初涉及 $(k+1)^{2d}$ 的插值操作与相关形状函数在 $k$ 度的 $d$ 维的张量积正交点分解为 $d$ 的一维操作，每个成本 $(k+1)^{d+1}$ 。在三维中，这将复杂度降低了 $k$ 的两个次方。当以每个自由度的复杂度来衡量时，复杂度是 $\mathcal O(k)$ 的多项式程度。由于DG中存在面积分，也由于对正交点的操作涉及更多的内存传输，这两者的规模为 $\mathcal O(1)$ ，观察到的复杂度通常是中等 $k\leq 10$ 的常数。这意味着，以每秒自由度计算，高阶方法的评估吞吐量与低阶方法相同。

关于算法的更多信息可参见Martin Kronbichler和Katharina Kormann的预印本 <br>  <a href="https://arxiv.org/abs/1711.03590">Fast matrix-free evaluation of
discontinuous Galerkin finite element operators</a>，arXiv:1711.03590。

<a name="ThesymmetricinteriorpenaltyformulationfortheLaplacian"></a><h3>The symmetric interior penalty formulation for the Laplacian</h3>


在本教程中，我们以无矩阵的DG框架为例，对拉普拉斯进行内部惩罚离散化，即与Step-39教程中使用的方案相同。拉普拉斯的离散化是由以下弱形式给出的

@f{align*}
&\sum_{K\in\text{cells}} \left(\nabla v_h, \nabla u_h\right)_{K}+\\
&\sum_{F\in\text{faces}}\Big(-\left<\jump{v_h}, \average{\nabla u_h}\right>_{F} - \left<\average{\nabla v_h}, \jump{u_h}\right>_{F} + \left<\jump{v_h}, \sigma \jump{u_h}\right>_{F}\Big) \\
&= \sum_{K\in\text{cells}}\left(v_h, f\right)_{K},


@f}

其中 $\jump{v} = v^- \mathbf{n}^- + v^+ \mathbf{n}^+ = \mathbf n^{-}
\left(v^- - v^+\right)$ 表示数量 $v$ 从两个相关单元 $K^-$ 和 $K^+$ 的定向跳跃， $\average{v}=\frac{v^- + v^+}{2}$ 是两边的平均值。

方程中的项代表分项积分后的单元积分，由于分项积分和插入平均通量而在元素界面产生的原始一致性项，为恢复基础矩阵的对称性而添加的邻接一致性项，以及系数为 $\sigma$ 的惩罚项，其大小等于单元在面的法线方向的长度乘以 $k(k+1)$  ，见步骤39。惩罚项的选择是为了使逆向估计成立，并且最终的弱形式是强制性的，即在离散设置中是正定的。邻接一致性项和惩罚项涉及元素界面上的跳跃 $\jump{u_h}$ ，这在解析解中消失了 $u$  。因此，这些项与原始PDE一致，保证了该方法可以保留最佳收敛顺序。

在下面的实现中，我们通过将法向量 $\mathbf{n}^-$ 从跳跃项移到导数中来实现上面的弱形式，形成<i>normal</i>形式的导数。这使得在正交点上的实现稍微有效一些，因为我们只需要处理标量项而不是张量，而且在数学上是等价的。

对于边界条件，我们使用所谓的镜像原理，通过从内部解 $u^-$ 结合给定的边界数据进行外推来定义<i>artificial</i>外部值 $u^+$ ，在迪里希特边界设置 $u^+ = -u^- + 2
g_\text{D}$ 和 $\mathbf{n}^-\cdot \nabla u^+ = \mathbf{n}^-\cdot \nabla u^-$ ，在纽曼边界设置 $u^+=u^-$ 和 $\mathbf{n}^-\cdot \nabla u^+ =


-\mathbf{n}^-\cdot \nabla u^- + 2 g_\text{N}$ ，对于给定的迪里希特值 $g_\text{D}$ 和纽曼值 $g_\text{N}$  。然后将这些表达式插入上述弱形式中。涉及已知量 $g_\text{D}$ 和 $g_\text{N}$ 的贡献最终被移到右手边，而未知值 $\varepsilon(\Delta \mathbf{u}^n(\mathbf{x}_q))$ 被保留在左手边，对矩阵项的贡献与内部面类似。经过这些操作，得到了与步骤39中相同的弱形式。

<a name="FaceintegrationsupportinMatrixFreeandFEFaceEvaluation"></a><h3>Face integration support in MatrixFree and FEFaceEvaluation</h3>


deal.II的无矩阵框架为实现上述离散化方程的作用提供了必要的基础设施。相对于我们在步骤37和步骤48中使用的 MatrixFree::cell_loop() ，我们现在用 MatrixFree::loop() 建立一个代码，它需要三个函数指针，一个用于单元积分，一个用于内面积分，一个用于边界面积分（与步骤39教程程序中使用的MeshWorker的设计相类似）。在这三个函数中的每一个，我们都会在正交点上实现各自的条款。对于向量条目和正交点上的值和梯度之间的插值，我们使用FEEvaluation类来处理单元贡献，FEFaceEvaluation类来处理面的贡献。这些函数的基本用法已经在step-37教程程序中得到了广泛的讨论。

在 MatrixFree::loop(), 中，所有内部面都正好被访问一次，所以必须确保计算来自测试函数 $v_h^-$ 和 $v_h^+$ 的贡献。考虑到两边的测试函数确实是独立的，上面的弱形式实际上意味着我们在用测试函数的法向导数进行测试时，向名为`phi_inner`和`phi_outer`的FEFaceEvaluation对象提交相同的贡献，而在用测试函数的值进行测试时，则提交相反符号的值，因为后者由于跳跃项而涉及相反符号。对于不同细化程度的单元之间的面，从细化的一侧进行整合，FEFaceEvaluation自动执行内插到粗略一侧的子面。因此，一个悬空的节点永远不会明确地出现在用户实现的弱形式中。

每个面被精确访问一次的事实也适用于用MPI并行化时不同处理器之间的子域边界的那些面，其中一个单元属于一个处理器，一个属于另一个。 MatrixFree::reinit() 中的设置将面分成了两边，最终只报告了 MatrixFree::n_inner_face_batches() 和 MatrixFree::n_boundary_face_batches(), 中分别实际处理的面。注意，与步骤37中讨论的单元积分相类似，deal.II在几个面上应用矢量化，以使用SIMD，在我们称之为<i>batch of faces</i>的东西上用一条指令工作。面的批次与单元的批次是独立的，尽管处理面的积分的时间与处理各自单元的积分的时间保持一致，以增加数据定位。

这个程序中的另一个新东西是，我们不再像 FEEvaluation::read_dof_values() 或 FEEvaluation::distribute_local_to_global() 那样将向量访问从求值和积分步骤中分割出来，而是分别调用组合函数 FEEvaluation::gather_evaluate() 和 FEEvaluation::integrate_scatter(), 。这对面积分很有用，因为根据面的评估内容，并非所有单元的向量项都必须首先被触及。例如，想想节点元素FE_DGQ的情况，节点点在元素表面。如果我们对面的形状函数值感兴趣，只有 $(k+ 1)^{d-1}$ 个自由度以非显著的方式对它们做出贡献（用更专业的方式来说，只有 $(k+1)^{d-1}$ 个形状函数在面有非零支持，并对 FiniteElement::has_support_on_face()). 个自由度返回真值。 当与单元的 $(k+1)^d$ 个自由度相比，这要少一个幂。

现在，我们当然不只对函数值感兴趣，而且对单元格上的导数也感兴趣。幸运的是，在deal.II中，有一个元素将这种减少访问的属性也扩展到面的导数上，即FE_DGQHermite元素。

<a name="TheFE_DGQHermiteelement"></a><h3>The FE_DGQHermite element</h3>


FE_DGQHermite元素属于FE_DGQ元素家族，即其形状函数是一维多项式的张量乘积，该元素是完全不连续的。与通常的FE_DGQ元素中的节点特征相反，FE_DGQHermite元素是一个基于Hermite-like概念的节点贡献和导数贡献的混合物。基本的多项式类是 Polynomials::HermiteLikeInterpolation ，可以总结为以下几点。对于三次多项式，我们用两个多项式来表示单位区间左端的函数值和第一次导数， $x=0$  ，用两个多项式来表示函数值和第一次导数以及单位区间的右端， $x=1$  。在相反的两端，形状函数的值和一阶导数都是零，确保四个基础函数中只有两个对各自端部的值和导数有贡献。然而，我们偏离了经典的Hermite内插法，没有严格地为值和一导数指定一个自由度，而是允许一导数是第一和第二形状函数的线性组合。这样做是为了改善插值的调节。另外，当度数超过三时，我们以类似拉格朗日的方式在元素内部增加节点点，在 $x=0$ 和 $x=1$ 这两个点上结合双零。这些额外节点的位置是由一些雅可比多项式的零点决定的，在类的描述中解释过  Polynomials::HermiteLikeInterpolation.  。

使用这个元素，我们只需要访问 $2(k+1)^{d-1}$ 自由度来计算一个面的值和导数。检查是否满足Hermite属性是在 FEFaceEvaluation::gather_evaluate() 和 FEFaceEvaluation::integrate_scatter() 中透明地完成的，它们检查基础的类型，并尽可能地减少对数据的访问。显然，如果我们将 FEFaceEvaluation::read_dof_values() 与 FEFaceEvaluation::evaluate(), 分开，这将是不可能的，因为我们需要读取的条目数量取决于导数的类型（只有值，一阶导数，等等），因此必须交给`read_dof_values()`。

这种优化不仅对计算面的积分有用，而且对MPI的鬼层交换也有用。在一个天真的交换中，如果一个单元的所有自由度由另一个处理器负责计算面的贡献，我们就需要把这个单元的所有自由度发送给另一个处理器。由于我们知道在用FEFaceEvaluation进行的评估中只有部分自由度被触及，所以自然只交换相关的自由度。 MatrixFree::loop() 函数在与 LinearAlgebra::distributed::Vector. 结合时支持选定的数据交换。为了实现这一点，我们需要告诉循环我们要对脸部做什么样的评价，使用 MatrixFree::DataAccessOnFaces, 类型的参数，正如在下面 `LaplaceOperator::vmult()` 的实现中可以看到。在这种情况下，数据交换的方式是如下的。矢量中的幽灵层数据仍然假装代表所有的自由度，这样，FEFaceEvaluation可以继续读取数值，就像单元格是本地拥有的一样。数据交换例程负责将数据打包和解包成这种格式的任务。虽然这听起来很复杂，但我们将在下面的结果部分显示，通过与不指定面孔上的数据访问的基线代码进行性能比较，这确实得到了回报。

<a name="AnapproximateblockJacobismootherusingthefastdiagonalizationmethod"></a><h3>An approximate block-Jacobi smoother using the fast diagonalization method</h3>


按照step-37程序的传统，我们再次用共轭梯度求解器内的几何多网格预处理器解决泊松问题。在这个教程程序中，我们没有计算对角线并使用基本的PreconditionChebyshev作为平滑器，而是选择了一种不同的策略。我们实现了一个块状Jacobi预处理程序，其中块状指的是一个单元上的所有自由度。我们没有在预处理程序中建立完整的单元格矩阵并应用它的LU因子化（或逆），这种操作将受到严重的内存带宽限制，因此相当缓慢；我们通过一种称为快速对角化法的特殊技术对块的逆进行近似处理。

该方法的思路是利用单元矩阵的结构。对于在笛卡尔网格上离散的恒定系数的拉普拉斯，单元矩阵 $L$ 可以写为

@f{align*}{
L &= A_1 \otimes M_0 + M_1 \otimes A_0


@f}

在二维和

@f{align*}{
L &= A_2 \otimes M_1 \otimes M_0 + M_2 \otimes A_1 \otimes M_0 + M_2 \otimes M_1 \otimes A_0


@f}

在三维中。矩阵 $A_0$ 和 $A_1$ 表示一维拉普拉斯矩阵（包括与当前单元值 $u^-_h$ 和 $v^-_h$ 相关的单元和面项）， $M_0$ 和 $M_1$ 是质量矩阵。请注意，一旦单元上有非恒定系数或几何形状不再恒定，这种简单的张量乘积结构就会消失。我们提到，类似的设置也可以用来用这个最终的张量积形式的矩阵来替代计算的积分，这将把算子评估的操作减少到一半以下。然而，考虑到这只适用于直角坐标单元和恒定系数的情况，这是一个相当狭窄的情况，我们避免继续探讨这个想法。

有趣的是，由于1964年<a
href="http://dl.acm.org/citation.cfm?id=2716130">R. E. Lynch, J. R. Rice,
D. H. Thomas, Direct solution of partial difference equations by tensor
product methods, Numerische Mathematik 6, 185-199</a>引入的方法，矩阵 $L$ 的精确逆值可以通过张量积找到。

@f{align*}{
L^{-1} &= S_1 \otimes S_0 (\Lambda_1 \otimes I + I \otimes \Lambda_0)^{-1}
S_1^\mathrm T \otimes S_0^\mathrm T,


@f}

其中  $S_d$  是给定张量方向上广义特征值问题的特征向量矩阵  $d$  。

@f{align*}{
A_d s  &= \lambda M_d s, \quad d = 0, \ldots,\mathrm{dim-1},


@f}

和 $\Lambda_d$ 是代表广义特征值 $\lambda$ 的对角矩阵。请注意，向量 $s$ 是这样的：它们同时对角化 $A_d$ 和 $M_d$ ，即 $S_d^{\mathrm T} A_d S_d =
\Lambda_d$ 和 $S_d^{\mathrm T} M_d S_d = I$  。

deal.II库使用这个概念实现了一个类，叫做TensorProductMatrixSymmetricSum。

为了这个程序，我们坚持使用常数系数和直角坐标系网格，尽管基于张量积的近似版本对于更一般的网格仍然是可能的，而且算子评估本身当然是通用的。另外，我们也不关心自适应网格，因为多网格算法需要获得不同细化边缘的通量矩阵，如步骤39所解释的。然而，我们做的一件事是仍然将我们的块状Jacobi预处理包在PreconditionChebyshev里面。这个类使我们不必寻找适当的松弛参数（对于块-雅各比平滑器来说，二维的松弛参数约为0.7，三维的松弛参数约为0.5），而且通常比普通的雅各比平滑器提高了一些平滑效率，因为当设置切比雪夫多项式的度数为1或2时，它可以降低解的时间。

请注意，块状Jacobi平滑器有一个额外的好处：快速对角线化方法也可以解释为从FE_DGQHermite的Hermite-like多项式转变为一个单元拉普拉斯是对角线的基础。因此，它抵消了基础的影响，无论我们使用FE_DGQHermite还是FE_DGQ，都会得到相同的迭代次数。这与使用只有对角线的PreconditionChebyshev类（点Jacobi方案）相比，FE_DGQ和FE_DGQHermite确实表现不同，FE_DGQ需要的迭代次数比FE_DGQHermite少2-5次，尽管对类似Hermite的形状函数做了修改以确保良好的调节。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The include files are essentially the same as in step-37, with the
 * exception of the finite element class FE_DGQHermite instead of FE_Q. All
 * functionality for matrix-free computations on face integrals is already
 * contained in `fe_evaluation.h`.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/tensor_product_matrix.h>
 * 
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_tools.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer_matrix_free.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * 
 * namespace Step59
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * As in step-37, we collect the dimension and polynomial degree as
 * constants here at the top of the program for simplicity. As opposed to
 * step-37, we choose a really high order method this time with degree 8
 * where any implementation not using sum factorization would become
 * prohibitively slow compared to the implementation with MatrixFree which
 * provides an efficiency that is essentially the same as at degrees two or
 * three. Furthermore, all classes in this tutorial program are templated,
 * so it would be easy to select the degree at run time from an input file
 * or a command-line argument by adding instantiations of the appropriate
 * degrees in the `main()` function.
 * 

 * 
 * 
 * @code
 *   const unsigned int degree_finite_element = 8;
 *   const unsigned int dimension             = 3;
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * In analogy to step-7, we define an analytic solution that we try to
 * reproduce with our discretization. Since the aim of this tutorial is to
 * show matrix-free methods, we choose one of the simplest possibilities,
 * namely a cosine function whose derivatives are simple enough for us to
 * compute analytically. Further down, the wave number 2.4 we select here
 * will be matched with the domain extent in $x$-direction that is 2.5, such
 * that we obtain a periodic solution at $x = 2.5$ including $6pi$ or three
 * full wave revolutions in the cosine. The first function defines the
 * solution and its gradient for expressing the analytic solution for the
 * Dirichlet and Neumann boundary conditions, respectively. Furthermore, a
 * class representing the negative Laplacian of the solution is used to
 * represent the right hand side (forcing) function that we use to match the
 * given analytic solution in the discretized version (manufactured
 * solution).
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int = 0) const override final
 *     {
 *       double val = 1.;
 *       for (unsigned int d = 0; d < dim; ++d)
 *         val *= std::cos(numbers::PI * 2.4 * p[d]);
 *       return val;
 *     }
 * 
 *     virtual Tensor<1, dim> gradient(const Point<dim> &p,
 *                                     const unsigned int = 0) const override final
 *     {
 *       const double   arg = numbers::PI * 2.4;
 *       Tensor<1, dim> grad;
 *       for (unsigned int d = 0; d < dim; ++d)
 *         {
 *           grad[d] = 1.;
 *           for (unsigned int e = 0; e < dim; ++e)
 *             if (d == e)
 *               grad[d] *= -arg * std::sin(arg * p[e]);
 *             else
 *               grad[d] *= std::cos(arg * p[e]);
 *         }
 *       return grad;
 *     }
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int = 0) const override final
 *     {
 *       const double arg = numbers::PI * 2.4;
 *       double       val = 1.;
 *       for (unsigned int d = 0; d < dim; ++d)
 *         val *= std::cos(arg * p[d]);
 *       return dim * arg * arg * val;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Matrixfreeimplementation"></a> 
 * <h3>Matrix-free implementation</h3>
 * 

 * 
 * The `LaplaceOperator` class is similar to the respective class in
 * step-37. A significant difference is that we do not derive the class from
 * MatrixFreeOperators::Base because we want to present some additional
 * features of MatrixFree::loop() that are not available in the
 * general-purpose class MatrixFreeOperators::Base. We derive the class from
 * the Subscriptor class to be able to use the operator within the Chebyshev
 * preconditioner because that preconditioner stores the underlying matrix
 * via a SmartPointer.
 *   

 * 
 * Given that we implement a complete matrix interface by hand, we need to
 * add an `initialize()` function, an `m()` function, a `vmult()` function,
 * and a `Tvmult()` function that were previously provided by
 * MatrixFreeOperators::Base. Our LaplaceOperator also contains a member
 * function `get_penalty_factor()` that centralizes the selection of the
 * penalty parameter in the symmetric interior penalty method according to
 * step-39.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   class LaplaceOperator : public Subscriptor
 *   {
 *   public:
 *     using value_type = number;
 * 
 *     LaplaceOperator() = default;
 * 
 *     void initialize(std::shared_ptr<const MatrixFree<dim, number>> data);
 * 
 *     void clear();
 * 
 *     types::global_dof_index m() const;
 * 
 *     void initialize_dof_vector(
 *       LinearAlgebra::distributed::Vector<number> &vec) const;
 * 
 *     std::shared_ptr<const MatrixFree<dim, number>> get_matrix_free() const;
 * 
 *     void vmult(LinearAlgebra::distributed::Vector<number> &      dst,
 *                const LinearAlgebra::distributed::Vector<number> &src) const;
 * 
 *     void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst,
 *                 const LinearAlgebra::distributed::Vector<number> &src) const;
 * 
 *     number get_penalty_factor() const
 *     {
 *       return 1.0 * fe_degree * (fe_degree + 1);
 *     }
 * 
 *   private:
 *     void
 *     apply_cell(const MatrixFree<dim, number> &                   data,
 *                LinearAlgebra::distributed::Vector<number> &      dst,
 *                const LinearAlgebra::distributed::Vector<number> &src,
 *                const std::pair<unsigned int, unsigned int> &cell_range) const;
 * 
 *     void
 *     apply_face(const MatrixFree<dim, number> &                   data,
 *                LinearAlgebra::distributed::Vector<number> &      dst,
 *                const LinearAlgebra::distributed::Vector<number> &src,
 *                const std::pair<unsigned int, unsigned int> &face_range) const;
 * 
 *     void apply_boundary(
 *       const MatrixFree<dim, number> &                   data,
 *       LinearAlgebra::distributed::Vector<number> &      dst,
 *       const LinearAlgebra::distributed::Vector<number> &src,
 *       const std::pair<unsigned int, unsigned int> &     face_range) const;
 * 
 *     std::shared_ptr<const MatrixFree<dim, number>> data;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * The `%PreconditionBlockJacobi` class defines our custom preconditioner for
 * this problem. As opposed to step-37 which was based on the matrix
 * diagonal, we here compute an approximate inversion of the diagonal blocks
 * in the discontinuous Galerkin method by using the so-called fast
 * diagonalization method discussed in the introduction.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   class PreconditionBlockJacobi
 *   {
 *   public:
 *     using value_type = number;
 * 
 *     void clear()
 *     {
 *       cell_matrices.clear();
 *     }
 * 
 *     void initialize(const LaplaceOperator<dim, fe_degree, number> &op);
 * 
 *     void vmult(LinearAlgebra::distributed::Vector<number> &      dst,
 *                const LinearAlgebra::distributed::Vector<number> &src) const;
 * 
 *     void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst,
 *                 const LinearAlgebra::distributed::Vector<number> &src) const
 *     {
 *       vmult(dst, src);
 *     }
 * 
 *   private:
 *     std::shared_ptr<const MatrixFree<dim, number>> data;
 *     std::vector<TensorProductMatrixSymmetricSum<dim,
 *                                                 VectorizedArray<number>,
 *                                                 fe_degree + 1>>
 *       cell_matrices;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * This free-standing function is used in both the `LaplaceOperator` and
 * `%PreconditionBlockJacobi` classes to adjust the ghost range. This function
 * is necessary because some of the vectors that the `vmult()` functions are
 * supplied with are not initialized properly with
 * `LaplaceOperator::initialize_dof_vector` that includes the correct layout
 * of ghost entries, but instead comes from the MGTransferMatrixFree class
 * that has no notion on the ghost selection of the matrix-free classes. To
 * avoid index confusion, we must adjust the ghost range before actually
 * doing something with these vectors. Since the vectors are kept around in
 * the multigrid smoother and transfer classes, a vector whose ghost range
 * has once been adjusted will remain in this state throughout the lifetime
 * of the object, so we can use a shortcut at the start of the function to
 * see whether the partitioner object of the distributed vector, which is
 * stored as a shared pointer, is the same as the layout expected by
 * MatrixFree, which is stored in a data structure accessed by
 * MatrixFree::get_dof_info(0), where the 0 indicates the DoFHandler number
 * from which this was extracted; we only use a single DoFHandler in
 * MatrixFree, so the only valid number is 0 here.
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number>
 *   void adjust_ghost_range_if_necessary(
 *     const MatrixFree<dim, number> &                   data,
 *     const LinearAlgebra::distributed::Vector<number> &vec)
 *   {
 *     if (vec.get_partitioner().get() ==
 *         data.get_dof_info(0).vector_partitioner.get())
 *       return;
 * 
 *     LinearAlgebra::distributed::Vector<number> copy_vec(vec);
 *     const_cast<LinearAlgebra::distributed::Vector<number> &>(vec).reinit(
 *       data.get_dof_info(0).vector_partitioner);
 *     const_cast<LinearAlgebra::distributed::Vector<number> &>(vec)
 *       .copy_locally_owned_data_from(copy_vec);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next five functions to clear and initialize the `LaplaceOperator`
 * class, to return the shared pointer holding the MatrixFree data
 * container, as well as the correct initialization of the vector and
 * operator sizes are the same as in step-37 or rather
 * MatrixFreeOperators::Base.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::clear()
 *   {
 *     data.reset();
 *   }
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::initialize(
 *     std::shared_ptr<const MatrixFree<dim, number>> data)
 *   {
 *     this->data = data;
 *   }
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   std::shared_ptr<const MatrixFree<dim, number>>
 *   LaplaceOperator<dim, fe_degree, number>::get_matrix_free() const
 *   {
 *     return data;
 *   }
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::initialize_dof_vector(
 *     LinearAlgebra::distributed::Vector<number> &vec) const
 *   {
 *     data->initialize_dof_vector(vec);
 *   }
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   types::global_dof_index LaplaceOperator<dim, fe_degree, number>::m() const
 *   {
 *     Assert(data.get() != nullptr, ExcNotInitialized());
 *     return data->get_dof_handler().n_dofs();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This function implements the action of the LaplaceOperator on a vector
 * `src` and stores the result in the vector `dst`. When compared to
 * step-37, there are four new features present in this call.
 *   

 * 
 * The first new feature is the `adjust_ghost_range_if_necessary` function
 * mentioned above that is needed to fit the vectors to the layout expected
 * by FEEvaluation and FEFaceEvaluation in the cell and face functions.
 *   

 * 
 * The second new feature is the fact that we do not implement a
 * `vmult_add()` function as we did in step-37 (through the virtual function
 * MatrixFreeOperators::Base::vmult_add()), but directly implement a
 * `vmult()` functionality. Since both cell and face integrals will sum into
 * the destination vector, we must of course zero the vector somewhere. For
 * DG elements, we are given two options &ndash; one is to use
 * FEEvaluation::set_dof_values() instead of
 * FEEvaluation::distribute_local_to_global() in the `apply_cell` function
 * below. This works because the loop layout in MatrixFree is such that cell
 * integrals always touch a given vector entry before the face
 * integrals. However, this really only works for fully discontinuous bases
 * where every cell has its own degrees of freedom, without any sharing with
 * neighboring results. An alternative setup, the one chosen here, is to let
 * the MatrixFree::loop() take care of zeroing the vector. This can be
 * thought of as simply calling `dst = 0;` somewhere in the code. The
 * implementation is more involved for supported vectors such as
 * `LinearAlgebra::distributed::Vector`, because we aim to not zero the
 * whole vector at once. Doing the zero operation on a small enough pieces
 * of a few thousands of vector entries has the advantage that the vector
 * entries that get zeroed remain in caches before they are accessed again
 * in FEEvaluation::distribute_local_to_global() and
 * FEFaceEvaluation::distribute_local_to_global(). Since matrix-free
 * operator evaluation is really fast, just zeroing a large vector can
 * amount to up to a 25% of the operator evaluation time, and we obviously
 * want to avoid this cost. This option of zeroing the vector is also
 * available for MatrixFree::cell_loop and for continuous bases, even though
 * it was not used in the step-37 or step-48 tutorial programs.
 *   

 * 
 * The third new feature is the way we provide the functions to compute on
 * cells, inner faces, and boundary faces: The class MatrixFree has a
 * function called `loop` that takes three function pointers to the three
 * cases, allowing to separate the implementations of different things. As
 * explained in step-37, these function pointers can be `std::function`
 * objects or member functions of a class. In this case, we use pointers to
 * member functions.
 *   

 * 
 * The final new feature are the last two arguments of type
 * MatrixFree::DataAccessOnFaces that can be given to
 * MatrixFree::loop(). This class passes the type of data access for face
 * integrals to the MPI data exchange routines
 * LinearAlgebra::distributed::Vector::update_ghost_values() and
 * LinearAlgebra::distributed::Vector::compress() of the parallel
 * vectors. The purpose is to not send all degrees of freedom of a
 * neighboring element, but to reduce the amount of data to what is really
 * needed for the computations at hand. The data exchange is a real
 * bottleneck in particular for high-degree DG methods, therefore a more
 * restrictive way of exchange is often beneficial. The enum field
 * MatrixFree::DataAccessOnFaces can take the value `none`, which means that
 * no face integrals at all are done, which would be analogous to
 * MatrixFree::cell_loop(), the value `values` meaning that only shape
 * function values (but no derivatives) are used on faces, and the value
 * `gradients` when also first derivatives on faces are accessed besides the
 * values. A value `unspecified` means that all degrees of freedom will be
 * exchanged for the faces that are located at the processor boundaries and
 * designated to be worked on at the local processor.
 *   

 * 
 * To see how the data can be reduced, think of the case of the nodal
 * element FE_DGQ with node points on the element surface, where only
 * $(k+1)^{d-1}$ degrees of freedom contribute to the values on a face for
 * polynomial degree $k$ in $d$ space dimensions, out of the $(k+1)^d$
 * degrees of freedom of a cell. A similar reduction is also possible for
 * the interior penalty method that evaluates values and first derivatives
 * on the faces. When using a Hermite-like basis in 1D, only up to two basis
 * functions contribute to the value and derivative. The class FE_DGQHermite
 * implements a tensor product of this concept, as discussed in the
 * introduction. Thus, only $2(k+1)^{d-1}$ degrees of freedom must be
 * exchanged for each face, which is a clear win once $k$ gets larger than
 * four or five. Note that this reduced exchange of FE_DGQHermite is valid
 * also on meshes with curved boundaries, as the derivatives are taken on
 * the reference element, whereas the geometry only mixes them on the
 * inside. Thus, this is different from the attempt to obtain $C^1$
 * continuity with continuous Hermite-type shape functions where the
 * non-Cartesian case changes the picture significantly. Obviously, on
 * non-Cartesian meshes the derivatives also include tangential derivatives
 * of shape functions beyond the normal derivative, but those only need the
 * function values on the element surface, too. Should the element not
 * provide any compression, the loop automatically exchanges all entries for
 * the affected cells.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::vmult(
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src) const
 *   {
 *     adjust_ghost_range_if_necessary(*data, dst);
 *     adjust_ghost_range_if_necessary(*data, src);
 *     data->loop(&LaplaceOperator::apply_cell,
 *                &LaplaceOperator::apply_face,
 *                &LaplaceOperator::apply_boundary,
 *                this,
 *                dst,
 *                src,
 *                /*zero_dst =*/true,
 *                MatrixFree<dim, number>::DataAccessOnFaces::gradients,
 *                MatrixFree<dim, number>::DataAccessOnFaces::gradients);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Since the Laplacian is symmetric, the `Tvmult()` (needed by the multigrid
 * smoother interfaces) operation is simply forwarded to the `vmult()` case.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::Tvmult(
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src) const
 *   {
 *     vmult(dst, src);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The cell operation is very similar to step-37. We do not use a
 * coefficient here, though. The second difference is that we replaced the
 * two steps of FEEvaluation::read_dof_values() followed by
 * FEEvaluation::evaluate() by a single function call
 * FEEvaluation::gather_evaluate() which internally calls the sequence of
 * the two individual methods. Likewise, FEEvaluation::integrate_scatter()
 * implements the sequence of FEEvaluation::integrate() followed by
 * FEEvaluation::distribute_local_to_global(). In this case, these new
 * functions merely save two lines of code. However, we use them for the
 * analogy with FEFaceEvaluation where they are more important as
 * explained below.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::apply_cell(
 *     const MatrixFree<dim, number> &                   data,
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(src, EvaluationFlags::gradients);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_gradient(phi.get_gradient(q), q);
 *         phi.integrate_scatter(EvaluationFlags::gradients, dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The face operation implements the terms of the interior penalty method in
 * analogy to step-39, as explained in the introduction. We need two
 * evaluator objects for this task, one for handling the solution that comes
 * from the cell on one of the two sides of an interior face, and one for
 * handling the solution from the other side. The evaluators for face
 * integrals are called FEFaceEvaluation and take a boolean argument in the
 * second slot of the constructor to indicate which of the two sides the
 * evaluator should belong two. In FEFaceEvaluation and MatrixFree, we call
 * one of the two sides the `interior` one and the other the `exterior`
 * one. The name `exterior` refers to the fact that the evaluator from both
 * sides will return the same normal vector. For the `interior` side, the
 * normal vector points outwards, whereas it points inwards on the other
 * side, and is opposed to the outer normal vector of that cell. Apart from
 * the new class name, we again get a range of items to work with in
 * analogy to what was discussed in step-37, but for the interior faces in
 * this case. Note that the data structure of MatrixFree forms batches of
 * faces that are analogous to the batches of cells for the cell
 * integrals. All faces within a batch involve different cell numbers but
 * have the face number within the reference cell, have the same refinement
 * configuration (no refinement or the same subface), and the same
 * orientation, to keep SIMD operations simple and efficient.
 *   

 * 
 * Note that there is no implied meaning in interior versus exterior except
 * the logic decision of the orientation of the normal, which is pretty
 * random internally. One can in no way rely on a certain pattern of
 * assigning interior versus exterior flags, as the decision is made for the
 * sake of access regularity and uniformity in the MatrixFree setup
 * routines. Since most sane DG methods are conservative, i.e., fluxes look
 * the same from both sides of an interface, the mathematics are unaltered
 * if the interior/exterior flags are switched and normal vectors get the
 * opposite sign.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::apply_face(
 *     const MatrixFree<dim, number> &                   data,
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data,
 *                                                                          true);
 *     FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_outer(data,
 *                                                                          false);
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 * @endcode
 * 
 * On a given batch of faces, we first update the pointers to the
 * current face and then access the vector. As mentioned above, we
 * combine the vector access with the evaluation. In the case of face
 * integrals, the data access into the vector can be reduced for the
 * special case of an FE_DGQHermite basis as explained for the data
 * exchange above: Since only $2(k+1)^{d-1}$ out of the $(k+1)^d$ cell
 * degrees of freedom get multiplied by a non-zero value or derivative
 * of a shape function, this structure can be utilized for the
 * evaluation, significantly reducing the data access. The reduction
 * of the data access is not only beneficial because it reduces the
 * data in flight and thus helps caching, but also because the data
 * access to faces is often more irregular than for cell integrals when
 * gathering values from cells that are farther apart in the index
 * list of cells.
 * 
 * @code
 *         phi_inner.reinit(face);
 *         phi_inner.gather_evaluate(src,
 *                                   EvaluationFlags::values |
 *                                     EvaluationFlags::gradients);
 *         phi_outer.reinit(face);
 *         phi_outer.gather_evaluate(src,
 *                                   EvaluationFlags::values |
 *                                     EvaluationFlags::gradients);
 * 
 * @endcode
 * 
 * The next two statements compute the penalty parameter for the
 * interior penalty method. As explained in the introduction, we would
 * like to have a scaling like $\frac{1}{h_\text{i}}$ of the length
 * $h_\text{i}$ normal to the face. For a general non-Cartesian mesh,
 * this length must be computed by the product of the inverse Jacobian
 * times the normal vector in real coordinates. From this vector of
 * `dim` components, we must finally pick the component that is
 * oriented normal to the reference cell. In the geometry data stored
 * in MatrixFree, a permutation of the components in the Jacobian is
 * applied such that this latter direction is always the last
 * component `dim-1` (this is beneficial because reference-cell
 * derivative sorting can be made agnostic of the direction of the
 * face). This means that we can simply access the last entry `dim-1`
 * and must not look up the local face number in
 * `data.get_face_info(face).interior_face_no` and
 * `data.get_face_info(face).exterior_face_no`. Finally, we must also
 * take the absolute value of these factors as the normal could point
 * into either positive or negative direction.
 * 
 * @code
 *         const VectorizedArray<number> inverse_length_normal_to_face =
 *           0.5 * (std::abs((phi_inner.get_normal_vector(0) *
 *                            phi_inner.inverse_jacobian(0))[dim - 1]) +
 *                  std::abs((phi_outer.get_normal_vector(0) *
 *                            phi_outer.inverse_jacobian(0))[dim - 1]));
 *         const VectorizedArray<number> sigma =
 *           inverse_length_normal_to_face * get_penalty_factor();
 * 
 * @endcode
 * 
 * In the loop over the quadrature points, we eventually compute all
 * contributions to the interior penalty scheme. According to the
 * formulas in the introduction, the value of the test function gets
 * multiplied by the difference of the jump in the solution times the
 * penalty parameter and the average of the normal derivative in real
 * space. Since the two evaluators for interior and exterior sides get
 * different signs due to the jump, we pass the result with a
 * different sign here. The normal derivative of the test function
 * gets multiplied by the negative jump in the solution between the
 * interior and exterior side. This term, coined adjoint consistency
 * term, must also include the factor of $\frac{1}{2}$ in the code in
 * accordance with its relation to the primal consistency term that
 * gets the factor of one half due to the average in the test function
 * slot.
 * 
 * @code
 *         for (unsigned int q = 0; q < phi_inner.n_q_points; ++q)
 *           {
 *             const VectorizedArray<number> solution_jump =
 *               (phi_inner.get_value(q) - phi_outer.get_value(q));
 *             const VectorizedArray<number> average_normal_derivative =
 *               (phi_inner.get_normal_derivative(q) +
 *                phi_outer.get_normal_derivative(q)) *
 *               number(0.5);
 *             const VectorizedArray<number> test_by_value =
 *               solution_jump * sigma - average_normal_derivative;
 * 
 *             phi_inner.submit_value(test_by_value, q);
 *             phi_outer.submit_value(-test_by_value, q);
 * 
 *             phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q);
 *             phi_outer.submit_normal_derivative(-solution_jump * number(0.5), q);
 *           }
 * 
 * @endcode
 * 
 * Once we are done with the loop over quadrature points, we can do
 * the sum factorization operations for the integration loops on faces
 * and sum the results into the result vector, using the
 * `integrate_scatter` function. The name `scatter` reflects the
 * distribution of the vector data into scattered positions in the
 * vector using the same pattern as in `gather_evaluate`. Like before,
 * the combined integrate + write operation allows us to reduce the
 * data access.
 * 
 * @code
 *         phi_inner.integrate_scatter(EvaluationFlags::values |
 *                                       EvaluationFlags::gradients,
 *                                     dst);
 *         phi_outer.integrate_scatter(EvaluationFlags::values |
 *                                       EvaluationFlags::gradients,
 *                                     dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The boundary face function follows by and large the interior face
 * function. The only difference is the fact that we do not have a separate
 * FEFaceEvaluation object that provides us with exterior values $u^+$, but
 * we must define them from the boundary conditions and interior values
 * $u^-$. As explained in the introduction, we use $u^+ = -u^- + 2
 * g_\text{D}$ and $\mathbf{n}^-\cdot \nabla u^+ = \mathbf{n}^-\cdot \nabla
 * u^-$ on Dirichlet boundaries and $u^+=u^-$ and $\mathbf{n}^-\cdot \nabla
 * u^+ = -\mathbf{n}^-\cdot \nabla u^- + 2 g_\text{N}$ on Neumann
 * boundaries. Since this operation implements the homogeneous part, i.e.,
 * the matrix-vector product, we must neglect the boundary functions
 * $g_\text{D}$ and $g_\text{N}$ here, and added them to the right hand side
 * in `LaplaceProblem::compute_rhs()`. Note that due to extension of the
 * solution $u^-$ to the exterior via $u^+$, we can keep all factors $0.5$
 * the same as in the inner face function, see also the discussion in
 * step-39.
 *   

 * 
 * There is one catch at this point: The implementation below uses a boolean
 * variable `is_dirichlet` to switch between the Dirichlet and the Neumann
 * cases. However, we solve a problem where we also want to impose periodic
 * boundary conditions on some boundaries, namely along those in the $x$
 * direction. One might wonder how those conditions should be handled
 * here. The answer is that MatrixFree automatically treats periodic
 * boundaries as what they are technically, namely an inner face where the
 * solution values of two adjacent cells meet and must be treated by proper
 * numerical fluxes. Thus, all the faces on the periodic boundaries will
 * appear in the `apply_face()` function and not in this one.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::apply_boundary(
 *     const MatrixFree<dim, number> &                   data,
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data,
 *                                                                          true);
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 *         phi_inner.reinit(face);
 *         phi_inner.gather_evaluate(src,
 *                                   EvaluationFlags::values |
 *                                     EvaluationFlags::gradients);
 * 
 *         const VectorizedArray<number> inverse_length_normal_to_face =
 *           std::abs((phi_inner.get_normal_vector(0) *
 *                     phi_inner.inverse_jacobian(0))[dim - 1]);
 *         const VectorizedArray<number> sigma =
 *           inverse_length_normal_to_face * get_penalty_factor();
 * 
 *         const bool is_dirichlet = (data.get_boundary_id(face) == 0);
 * 
 *         for (unsigned int q = 0; q < phi_inner.n_q_points; ++q)
 *           {
 *             const VectorizedArray<number> u_inner = phi_inner.get_value(q);
 *             const VectorizedArray<number> u_outer =
 *               is_dirichlet ? -u_inner : u_inner;
 *             const VectorizedArray<number> normal_derivative_inner =
 *               phi_inner.get_normal_derivative(q);
 *             const VectorizedArray<number> normal_derivative_outer =
 *               is_dirichlet ? normal_derivative_inner : -normal_derivative_inner;
 *             const VectorizedArray<number> solution_jump = (u_inner - u_outer);
 *             const VectorizedArray<number> average_normal_derivative =
 *               (normal_derivative_inner + normal_derivative_outer) * number(0.5);
 *             const VectorizedArray<number> test_by_value =
 *               solution_jump * sigma - average_normal_derivative;
 *             phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q);
 *             phi_inner.submit_value(test_by_value, q);
 *           }
 *         phi_inner.integrate_scatter(EvaluationFlags::values |
 *                                       EvaluationFlags::gradients,
 *                                     dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Next we turn to the preconditioner initialization. As explained in the
 * introduction, we want to construct an (approximate) inverse of the cell
 * matrices from a product of 1D mass and Laplace matrices. Our first task
 * is to compute the 1D matrices, which we do by first creating a 1D finite
 * element. Instead of anticipating FE_DGQHermite<1> here, we get the finite
 * element's name from DoFHandler, replace the @p dim argument (2 or 3) by 1
 * to create a 1D name, and construct the 1D element by using FETools.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void PreconditionBlockJacobi<dim, fe_degree, number>::initialize(
 *     const LaplaceOperator<dim, fe_degree, number> &op)
 *   {
 *     data = op.get_matrix_free();
 * 
 *     std::string name = data->get_dof_handler().get_fe().get_name();
 *     name.replace(name.find('<') + 1, 1, "1");
 *     std::unique_ptr<FiniteElement<1>> fe_1d = FETools::get_fe_by_name<1>(name);
 * 
 * @endcode
 * 
 * As for computing the 1D matrices on the unit element, we simply write
 * down what a typical assembly procedure over rows and columns of the
 * matrix as well as the quadrature points would do. We select the same
 * Laplace matrices once and for all using the coefficients 0.5 for
 * interior faces (but possibly scaled differently in different directions
 * as a result of the mesh). Thus, we make a slight mistake at the
 * Dirichlet boundary (where the correct factor would be 1 for the
 * derivative terms and 2 for the penalty term, see step-39) or at the
 * Neumann boundary where the factor should be zero. Since we only use
 * this class as a smoother inside a multigrid scheme, this error is not
 * going to have any significant effect and merely affects smoothing
 * quality.
 * 
 * @code
 *     const unsigned int                                 N = fe_degree + 1;
 *     FullMatrix<double>                                 laplace_unscaled(N, N);
 *     std::array<Table<2, VectorizedArray<number>>, dim> mass_matrices;
 *     std::array<Table<2, VectorizedArray<number>>, dim> laplace_matrices;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         mass_matrices[d].reinit(N, N);
 *         laplace_matrices[d].reinit(N, N);
 *       }
 * 
 *     QGauss<1> quadrature(N);
 *     for (unsigned int i = 0; i < N; ++i)
 *       for (unsigned int j = 0; j < N; ++j)
 *         {
 *           double sum_mass = 0, sum_laplace = 0;
 *           for (unsigned int q = 0; q < quadrature.size(); ++q)
 *             {
 *               sum_mass += (fe_1d->shape_value(i, quadrature.point(q)) *
 *                            fe_1d->shape_value(j, quadrature.point(q))) *
 *                           quadrature.weight(q);
 *               sum_laplace += (fe_1d->shape_grad(i, quadrature.point(q))[0] *
 *                               fe_1d->shape_grad(j, quadrature.point(q))[0]) *
 *                              quadrature.weight(q);
 *             }
 *           for (unsigned int d = 0; d < dim; ++d)
 *             mass_matrices[d](i, j) = sum_mass;
 * 
 * @endcode
 * 
 * The left and right boundary terms assembled by the next two
 * statements appear to have somewhat arbitrary signs, but those are
 * correct as can be verified by looking at step-39 and inserting
 * the value -1 and 1 for the normal vector in the 1D case.
 * 
 * @code
 *           sum_laplace +=
 *             (1. * fe_1d->shape_value(i, Point<1>()) *
 *                fe_1d->shape_value(j, Point<1>()) * op.get_penalty_factor() +
 *              0.5 * fe_1d->shape_grad(i, Point<1>())[0] *
 *                fe_1d->shape_value(j, Point<1>()) +
 *              0.5 * fe_1d->shape_grad(j, Point<1>())[0] *
 *                fe_1d->shape_value(i, Point<1>()));
 * 
 *           sum_laplace +=
 *             (1. * fe_1d->shape_value(i, Point<1>(1.0)) *
 *                fe_1d->shape_value(j, Point<1>(1.0)) * op.get_penalty_factor() -
 *              0.5 * fe_1d->shape_grad(i, Point<1>(1.0))[0] *
 *                fe_1d->shape_value(j, Point<1>(1.0)) -
 *              0.5 * fe_1d->shape_grad(j, Point<1>(1.0))[0] *
 *                fe_1d->shape_value(i, Point<1>(1.0)));
 * 
 *           laplace_unscaled(i, j) = sum_laplace;
 *         }
 * 
 * @endcode
 * 
 * Next, we go through the cells and pass the scaled matrices to
 * TensorProductMatrixSymmetricSum to actually compute the generalized
 * eigenvalue problem for representing the inverse: Since the matrix
 * approximation is constructed as $A\otimes M + M\otimes A$ and the
 * weights are constant for each element, we can apply all weights on the
 * Laplace matrix and simply keep the mass matrices unscaled. In the loop
 * over cells, we want to make use of the geometry compression provided by
 * the MatrixFree class and check if the current geometry is the same as
 * on the last cell batch, in which case there is nothing to do. This
 * compression can be accessed by
 * FEEvaluation::get_mapping_data_index_offset() once `reinit()` has been
 * called.
 *     

 * 
 * Once we have accessed the inverse Jacobian through the FEEvaluation
 * access function (we take the one for the zeroth quadrature point as
 * they should be the same on all quadrature points for a Cartesian cell),
 * we check that it is diagonal and then extract the determinant of the
 * original Jacobian, i.e., the inverse of the determinant of the inverse
 * Jacobian, and set the weight as $\text{det}(J) / h_d^2$ according to
 * the 1D Laplacian times $d-1$ copies of the mass matrix.
 * 
 * @code
 *     cell_matrices.clear();
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data);
 *     unsigned int old_mapping_data_index = numbers::invalid_unsigned_int;
 *     for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 * 
 *         if (phi.get_mapping_data_index_offset() == old_mapping_data_index)
 *           continue;
 * 
 *         Tensor<2, dim, VectorizedArray<number>> inverse_jacobian =
 *           phi.inverse_jacobian(0);
 * 
 *         for (unsigned int d = 0; d < dim; ++d)
 *           for (unsigned int e = 0; e < dim; ++e)
 *             if (d != e)
 *               for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
 *                 AssertThrow(inverse_jacobian[d][e][v] == 0.,
 *                             ExcNotImplemented());
 * 
 *         VectorizedArray<number> jacobian_determinant = inverse_jacobian[0][0];
 *         for (unsigned int e = 1; e < dim; ++e)
 *           jacobian_determinant *= inverse_jacobian[e][e];
 *         jacobian_determinant = 1. / jacobian_determinant;
 * 
 *         for (unsigned int d = 0; d < dim; ++d)
 *           {
 *             const VectorizedArray<number> scaling_factor =
 *               inverse_jacobian[d][d] * inverse_jacobian[d][d] *
 *               jacobian_determinant;
 * 
 * @endcode
 * 
 * Once we know the factor by which we should scale the Laplace
 * matrix, we apply this weight to the unscaled DG Laplace matrix
 * and send the array to the class TensorProductMatrixSymmetricSum
 * for computing the generalized eigenvalue problem mentioned in
 * the introduction.
 * 

 * 
 * 
 * @code
 *             for (unsigned int i = 0; i < N; ++i)
 *               for (unsigned int j = 0; j < N; ++j)
 *                 laplace_matrices[d](i, j) =
 *                   scaling_factor * laplace_unscaled(i, j);
 *           }
 *         if (cell_matrices.size() <= phi.get_mapping_data_index_offset())
 *           cell_matrices.resize(phi.get_mapping_data_index_offset() + 1);
 *         cell_matrices[phi.get_mapping_data_index_offset()].reinit(
 *           mass_matrices, laplace_matrices);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The vmult function for the approximate block-Jacobi preconditioner is
 * very simple in the DG context: We simply need to read the values of the
 * current cell batch, apply the inverse for the given entry in the array of
 * tensor product matrix, and write the result back. In this loop, we
 * overwrite the content in `dst` rather than first setting the entries to
 * zero. This is legitimate for a DG method because every cell has
 * independent degrees of freedom. Furthermore, we manually write out the
 * loop over all cell batches, rather than going through
 * MatrixFree::cell_loop(). We do this because we know that we are not going
 * to need data exchange over the MPI network here as all computations are
 * done on the cells held locally on each processor.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void PreconditionBlockJacobi<dim, fe_degree, number>::vmult(
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src) const
 *   {
 *     adjust_ghost_range_if_necessary(*data, dst);
 *     adjust_ghost_range_if_necessary(*data, src);
 * 
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data);
 *     for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.read_dof_values(src);
 *         cell_matrices[phi.get_mapping_data_index_offset()].apply_inverse(
 *           ArrayView<VectorizedArray<number>>(phi.begin_dof_values(),
 *                                              phi.dofs_per_cell),
 *           ArrayView<const VectorizedArray<number>>(phi.begin_dof_values(),
 *                                                    phi.dofs_per_cell));
 *         phi.set_dof_values(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The definition of the LaplaceProblem class is very similar to
 * step-37. One difference is the fact that we add the element degree as a
 * template argument to the class, which would allow us to more easily
 * include more than one degree in the same program by creating different
 * instances in the `main()` function. The second difference is the
 * selection of the element, FE_DGQHermite, which is specialized for this
 * kind of equations.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void compute_rhs();
 *     void solve();
 *     void analyze_results() const;
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 * 
 *     FE_DGQHermite<dim> fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     MappingQ1<dim> mapping;
 * 
 *     using SystemMatrixType = LaplaceOperator<dim, fe_degree, double>;
 *     SystemMatrixType system_matrix;
 * 
 *     using LevelMatrixType = LaplaceOperator<dim, fe_degree, float>;
 *     MGLevelObject<LevelMatrixType> mg_matrices;
 * 
 *     LinearAlgebra::distributed::Vector<double> solution;
 *     LinearAlgebra::distributed::Vector<double> system_rhs;
 * 
 *     double             setup_time;
 *     ConditionalOStream pcout;
 *     ConditionalOStream time_details;
 *   };
 * 
 * 
 * 
 *   template <int dim, int fe_degree>
 *   LaplaceProblem<dim, fe_degree>::LaplaceProblem()
 *     :
 * #ifdef DEAL_II_WITH_P4EST
 *     triangulation(
 *       MPI_COMM_WORLD,
 *       Triangulation<dim>::limit_level_difference_at_vertices,
 *       parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
 *     ,
 * #else
 *     triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     ,
 * #endif
 *     fe(fe_degree)
 *     , dof_handler(triangulation)
 *     , setup_time(0.)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *     , time_details(std::cout,
 *                    false &&
 *                      Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * The setup function differs in two aspects from step-37. The first is that
 * we do not need to interpolate any constraints for the discontinuous
 * ansatz space, and simply pass a dummy AffineConstraints object into
 * Matrixfree::reinit(). The second change arises because we need to tell
 * MatrixFree to also initialize the data structures for faces. We do this
 * by setting update flags for the inner and boundary faces,
 * respectively. On the boundary faces, we need both the function values,
 * their gradients, JxW values (for integration), the normal vectors, and
 * quadrature points (for the evaluation of the boundary conditions),
 * whereas we only need shape function values, gradients, JxW values, and
 * normal vectors for interior faces. The face data structures in MatrixFree
 * are always built as soon as one of `mapping_update_flags_inner_faces` or
 * `mapping_update_flags_boundary_faces` are different from the default
 * value `update_default` of UpdateFlags.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void LaplaceProblem<dim, fe_degree>::setup_system()
 *   {
 *     Timer time;
 *     setup_time = 0;
 * 
 *     system_matrix.clear();
 *     mg_matrices.clear_elements();
 * 
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler.distribute_mg_dofs();
 * 
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl;
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Distribute DoFs               " << time.wall_time() << " s"
 *                  << std::endl;
 *     time.restart();
 * 
 *     AffineConstraints<double> dummy;
 *     dummy.close();
 * 
 *     {
 *       typename MatrixFree<dim, double>::AdditionalData additional_data;
 *       additional_data.tasks_parallel_scheme =
 *         MatrixFree<dim, double>::AdditionalData::none;
 *       additional_data.mapping_update_flags =
 *         (update_gradients | update_JxW_values | update_quadrature_points);
 *       additional_data.mapping_update_flags_inner_faces =
 *         (update_gradients | update_JxW_values | update_normal_vectors);
 *       additional_data.mapping_update_flags_boundary_faces =
 *         (update_gradients | update_JxW_values | update_normal_vectors |
 *          update_quadrature_points);
 *       const auto system_mf_storage =
 *         std::make_shared<MatrixFree<dim, double>>();
 *       system_mf_storage->reinit(
 *         mapping, dof_handler, dummy, QGauss<1>(fe.degree + 1), additional_data);
 *       system_matrix.initialize(system_mf_storage);
 *     }
 * 
 *     system_matrix.initialize_dof_vector(solution);
 *     system_matrix.initialize_dof_vector(system_rhs);
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free system      " << time.wall_time() << " s"
 *                  << std::endl;
 *     time.restart();
 * 
 *     const unsigned int nlevels = triangulation.n_global_levels();
 *     mg_matrices.resize(0, nlevels - 1);
 * 
 *     for (unsigned int level = 0; level < nlevels; ++level)
 *       {
 *         typename MatrixFree<dim, float>::AdditionalData additional_data;
 *         additional_data.tasks_parallel_scheme =
 *           MatrixFree<dim, float>::AdditionalData::none;
 *         additional_data.mapping_update_flags =
 *           (update_gradients | update_JxW_values);
 *         additional_data.mapping_update_flags_inner_faces =
 *           (update_gradients | update_JxW_values);
 *         additional_data.mapping_update_flags_boundary_faces =
 *           (update_gradients | update_JxW_values);
 *         additional_data.mg_level = level;
 *         const auto mg_mf_storage_level =
 *           std::make_shared<MatrixFree<dim, float>>();
 *         mg_mf_storage_level->reinit(mapping,
 *                                     dof_handler,
 *                                     dummy,
 *                                     QGauss<1>(fe.degree + 1),
 *                                     additional_data);
 * 
 *         mg_matrices[level].initialize(mg_mf_storage_level);
 *       }
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free levels      " << time.wall_time() << " s"
 *                  << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The computation of the right hand side is a bit more complicated than in
 * step-37. The cell term now consists of the negative Laplacian of the
 * analytical solution, `RightHandSide`, for which we need to first split up
 * the Point of VectorizedArray fields, i.e., a batch of points, into a
 * single point by evaluating all lanes in the VectorizedArray
 * separately. Remember that the number of lanes depends on the hardware; it
 * could be 1 for systems that do not offer vectorization (or where deal.II
 * does not have intrinsics), but it could also be 8 or 16 on AVX-512 of
 * recent Intel architectures.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void LaplaceProblem<dim, fe_degree>::compute_rhs()
 *   {
 *     Timer time;
 *     system_rhs                          = 0;
 *     const MatrixFree<dim, double> &data = *system_matrix.get_matrix_free();
 *     FEEvaluation<dim, fe_degree>   phi(data);
 *     RightHandSide<dim>             rhs_func;
 *     Solution<dim>                  exact_solution;
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             VectorizedArray<double> rhs_val = VectorizedArray<double>();
 *             Point<dim, VectorizedArray<double>> point_batch =
 *               phi.quadrature_point(q);
 *             for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
 *               {
 *                 Point<dim> single_point;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   single_point[d] = point_batch[d][v];
 *                 rhs_val[v] = rhs_func.value(single_point);
 *               }
 *             phi.submit_value(rhs_val, q);
 *           }
 *         phi.integrate_scatter(EvaluationFlags::values, system_rhs);
 *       }
 * 
 * @endcode
 * 
 * Secondly, we also need to apply the Dirichlet and Neumann boundary
 * conditions. This function is the missing part of to the function
 * `LaplaceOperator::apply_boundary()` function once the exterior solution
 * values $u^+ = -u^- + 2 g_\text{D}$ and $\mathbf{n}^-\cdot \nabla u^+ =
 * \mathbf{n}^-\cdot \nabla u^-$ on Dirichlet boundaries and $u^+=u^-$ and
 * $\mathbf{n}^-\cdot \nabla u^+ = -\mathbf{n}^-\cdot \nabla u^- + 2
 * g_\text{N}$ on Neumann boundaries are inserted and expanded in terms of
 * the boundary functions $g_\text{D}$ and $g_\text{N}$. One thing to
 * remember is that we move the boundary conditions to the right hand
 * side, so the sign is the opposite from what we imposed on the solution
 * part.
 *     

 * 
 * We could have issued both the cell and the boundary part through a
 * MatrixFree::loop part, but we choose to manually write the full loop
 * over all faces to learn how the index layout of face indices is set up
 * in MatrixFree: Both the inner faces and the boundary faces share the
 * index range, and all batches of inner faces have lower numbers than the
 * batches of boundary cells. A single index for both variants allows us
 * to easily use the same data structure FEFaceEvaluation for both cases
 * that attaches to the same data field, just at different positions. The
 * number of inner face batches (where a batch is due to the combination
 * of several faces into one for vectorization) is given by
 * MatrixFree::n_inner_face_batches(), whereas the number of boundary face
 * batches is given by MatrixFree::n_boundary_face_batches().
 * 
 * @code
 *     FEFaceEvaluation<dim, fe_degree> phi_face(data, true);
 *     for (unsigned int face = data.n_inner_face_batches();
 *          face < data.n_inner_face_batches() + data.n_boundary_face_batches();
 *          ++face)
 *       {
 *         phi_face.reinit(face);
 * 
 *         const VectorizedArray<double> inverse_length_normal_to_face =
 *           std::abs((phi_face.get_normal_vector(0) *
 *                     phi_face.inverse_jacobian(0))[dim - 1]);
 *         const VectorizedArray<double> sigma =
 *           inverse_length_normal_to_face * system_matrix.get_penalty_factor();
 * 
 *         for (unsigned int q = 0; q < phi_face.n_q_points; ++q)
 *           {
 *             VectorizedArray<double> test_value = VectorizedArray<double>(),
 *                                     test_normal_derivative =
 *                                       VectorizedArray<double>();
 *             Point<dim, VectorizedArray<double>> point_batch =
 *               phi_face.quadrature_point(q);
 * 
 *             for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
 *               {
 *                 Point<dim> single_point;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   single_point[d] = point_batch[d][v];
 * 
 * @endcode
 * 
 * The MatrixFree class lets us query the boundary_id of the
 * current face batch. Remember that MatrixFree sets up the
 * batches for vectorization such that all faces within a
 * batch have the same properties, which includes their
 * `boundary_id`. Thus, we can query that id here for the
 * current face index `face` and either impose the Dirichlet
 * case (where we add something to the function value) or the
 * Neumann case (where we add something to the normal
 * derivative).
 * 
 * @code
 *                 if (data.get_boundary_id(face) == 0)
 *                   test_value[v] = 2.0 * exact_solution.value(single_point);
 *                 else
 *                   {
 *                     Tensor<1, dim> normal;
 *                     for (unsigned int d = 0; d < dim; ++d)
 *                       normal[d] = phi_face.get_normal_vector(q)[d][v];
 *                     test_normal_derivative[v] =
 *                       -normal * exact_solution.gradient(single_point);
 *                   }
 *               }
 *             phi_face.submit_value(test_value * sigma - test_normal_derivative,
 *                                   q);
 *             phi_face.submit_normal_derivative(-0.5 * test_value, q);
 *           }
 *         phi_face.integrate_scatter(EvaluationFlags::values |
 *                                      EvaluationFlags::gradients,
 *                                    system_rhs);
 *       }
 * 
 * @endcode
 * 
 * Since we have manually run the loop over cells rather than using
 * MatrixFree::loop(), we must not forget to perform the data exchange
 * with MPI - or actually, we would not need that for DG elements here
 * because each cell carries its own degrees of freedom and cell and
 * boundary integrals only evaluate quantities on the locally owned
 * cells. The coupling to neighboring subdomain only comes in by the inner
 * face integrals, which we have not done here. That said, it does not
 * hurt to call this function here, so we do it as a reminder of what
 * happens inside MatrixFree::loop().
 * 
 * @code
 *     system_rhs.compress(VectorOperation::add);
 *     setup_time += time.wall_time();
 *     time_details << "Compute right hand side       " << time.wall_time()
 *                  << " s\n";
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The `solve()` function is copied almost verbatim from step-37. We set up
 * the same multigrid ingredients, namely the level transfer, a smoother,
 * and a coarse grid solver. The only difference is the fact that we do not
 * use the diagonal of the Laplacian for the preconditioner of the Chebyshev
 * iteration used for smoothing, but instead our newly resolved class
 * `%PreconditionBlockJacobi`. The mechanisms are the same, though.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void LaplaceProblem<dim, fe_degree>::solve()
 *   {
 *     Timer                            time;
 *     MGTransferMatrixFree<dim, float> mg_transfer;
 *     mg_transfer.build(dof_handler);
 *     setup_time += time.wall_time();
 *     time_details << "MG build transfer time        " << time.wall_time()
 *                  << " s\n";
 *     time.restart();
 * 
 *     using SmootherType =
 *       PreconditionChebyshev<LevelMatrixType,
 *                             LinearAlgebra::distributed::Vector<float>,
 *                             PreconditionBlockJacobi<dim, fe_degree, float>>;
 *     mg::SmootherRelaxation<SmootherType,
 *                            LinearAlgebra::distributed::Vector<float>>
 *                                                          mg_smoother;
 *     MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
 *     smoother_data.resize(0, triangulation.n_global_levels() - 1);
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       {
 *         if (level > 0)
 *           {
 *             smoother_data[level].smoothing_range     = 15.;
 *             smoother_data[level].degree              = 3;
 *             smoother_data[level].eig_cg_n_iterations = 10;
 *           }
 *         else
 *           {
 *             smoother_data[0].smoothing_range = 2e-2;
 *             smoother_data[0].degree          = numbers::invalid_unsigned_int;
 *             smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
 *           }
 *         smoother_data[level].preconditioner =
 *           std::make_shared<PreconditionBlockJacobi<dim, fe_degree, float>>();
 *         smoother_data[level].preconditioner->initialize(mg_matrices[level]);
 *       }
 *     mg_smoother.initialize(mg_matrices, smoother_data);
 * 
 *     MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
 *       mg_coarse;
 *     mg_coarse.initialize(mg_smoother);
 * 
 *     mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
 *       mg_matrices);
 * 
 *     Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
 *       mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
 * 
 *     PreconditionMG<dim,
 *                    LinearAlgebra::distributed::Vector<float>,
 *                    MGTransferMatrixFree<dim, float>>
 *       preconditioner(dof_handler, mg, mg_transfer);
 * 
 *     SolverControl solver_control(10000, 1e-12 * system_rhs.l2_norm());
 *     SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);
 *     setup_time += time.wall_time();
 *     time_details << "MG build smoother time        " << time.wall_time()
 *                  << "s\n";
 *     pcout << "Total setup time              " << setup_time << " s\n";
 * 
 *     time.reset();
 *     time.start();
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     pcout << "Time solve (" << solver_control.last_step() << " iterations)    "
 *           << time.wall_time() << " s" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Since we have solved a problem with analytic solution, we want to verify
 * the correctness of our implementation by computing the L2 error of the
 * numerical result against the analytic solution.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void LaplaceProblem<dim, fe_degree>::analyze_results() const
 *   {
 *     Vector<float> error_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Solution<dim>(),
 *                                       error_per_cell,
 *                                       QGauss<dim>(fe.degree + 2),
 *                                       VectorTools::L2_norm);
 *     pcout << "Verification via L2 error:    "
 *           << std::sqrt(
 *                Utilities::MPI::sum(error_per_cell.norm_sqr(), MPI_COMM_WORLD))
 *           << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The `run()` function sets up the initial grid and then runs the multigrid
 * program in the usual way. As a domain, we choose a rectangle with
 * periodic boundary conditions in the $x$-direction, a Dirichlet condition
 * on the front face in $y$ direction (i.e., the face with index number 2,
 * with boundary id equal to 0), and Neumann conditions on the back face as
 * well as the two faces in $z$ direction for the 3D case (with boundary id
 * equal to 1). The extent of the domain is a bit different in the $x$
 * direction (where we want to achieve a periodic solution given the
 * definition of `Solution`) as compared to the $y$ and $z$ directions.
 * 

 * 
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void LaplaceProblem<dim, fe_degree>::run()
 *   {
 *     const unsigned int n_ranks =
 *       Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
 *     pcout << "Running with " << n_ranks << " MPI process"
 *           << (n_ranks > 1 ? "es" : "") << ", element " << fe.get_name()
 *           << std::endl
 *           << std::endl;
 *     for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             Point<dim> upper_right;
 *             upper_right[0] = 2.5;
 *             for (unsigned int d = 1; d < dim; ++d)
 *               upper_right[d] = 2.8;
 *             GridGenerator::hyper_rectangle(triangulation,
 *                                            Point<dim>(),
 *                                            upper_right);
 *             triangulation.begin_active()->face(0)->set_boundary_id(10);
 *             triangulation.begin_active()->face(1)->set_boundary_id(11);
 *             triangulation.begin_active()->face(2)->set_boundary_id(0);
 *             for (unsigned int f = 3;
 *                  f < triangulation.begin_active()->n_faces();
 *                  ++f)
 *               triangulation.begin_active()->face(f)->set_boundary_id(1);
 * 
 *             std::vector<GridTools::PeriodicFacePair<
 *               typename Triangulation<dim>::cell_iterator>>
 *               periodic_faces;
 *             GridTools::collect_periodic_faces(
 *               triangulation, 10, 11, 0, periodic_faces);
 *             triangulation.add_periodicity(periodic_faces);
 * 
 *             triangulation.refine_global(6 - 2 * dim);
 *           }
 *         triangulation.refine_global(1);
 *         setup_system();
 *         compute_rhs();
 *         solve();
 *         analyze_results();
 *         pcout << std::endl;
 *       };
 *   }
 * } // namespace Step59
 * 
 * 
 * 
 * @endcode
 * 
 * There is nothing unexpected in the `main()` function. We call `MPI_Init()`
 * through the `MPI_InitFinalize` class, pass on the two parameters on the
 * dimension and the degree set at the top of the file, and run the Laplace
 * problem.
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step59;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
 * 
 *       LaplaceProblem<dimension, degree_finite_element> laplace_problem;
 *       laplace_problem.run();
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
 *     }
 * 
 *   return 0;
 * }
 * @endcode
examples/step-59/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


与第37步一样，我们在运行时间方面评估多网格求解器。  在两个空间维度的8度元素中，一个可能的输出可能如下。

@code
Running with 12 MPI processes, element FE_DGQHermite<2>(8)


Cycle 0
Number of degrees of freedom: 5184
Total setup time              0.0282445 s
Time solve (14 iterations)    0.0110712 s
Verification via L2 error:    1.66232e-07


Cycle 1
Number of degrees of freedom: 20736
Total setup time              0.0126282 s
Time solve (14 iterations)    0.0157021 s
Verification via L2 error:    2.91505e-10


Cycle 2
Number of degrees of freedom: 82944
Total setup time              0.0227573 s
Time solve (14 iterations)    0.026568 s
Verification via L2 error:    6.64514e-13


Cycle 3
Number of degrees of freedom: 331776
Total setup time              0.0604685 s
Time solve (14 iterations)    0.0628356 s
Verification via L2 error:    5.57513e-13


Cycle 4
Number of degrees of freedom: 1327104
Total setup time              0.154359 s
Time solve (13 iterations)    0.219555 s
Verification via L2 error:    3.08139e-12


Cycle 5
Number of degrees of freedom: 5308416
Total setup time              0.467764 s
Time solve (13 iterations)    1.1821 s
Verification via L2 error:    3.90334e-12


Cycle 6
Number of degrees of freedom: 21233664
Total setup time              1.73263 s
Time solve (13 iterations)    5.21054 s
Verification via L2 error:    4.94543e-12
@endcode



与第37步一样，随着问题大小的增加，CG的迭代次数保持不变。迭代次数要高一些，这是因为我们使用的切比雪夫多项式的度数较低（步骤37中为2比5），而且内部惩罚离散化的特征值分布也比较大。尽管如此，13次迭代将残差减少了12个数量级，或者说每次迭代几乎是9个系数，这表明总体上是一种非常有效的方法。特别是，当使用12个核心时，我们可以在5秒内解决一个具有2100万自由度的系统，这是一个非常好的效率。当然，在二维中，我们很好地进入了8级多项式的四舍五入体系；事实上，大约83000个自由度或0.025秒就足以完全收敛这个（简单的）分析解。

如果我们在三个空间维度上运行程序，并没有太大的变化，只是我们现在用做更高的多项式度数和不断增加的网格大小来做一些更有用的事情，因为舍入误差只在最细的网格上获得。尽管如此，令人瞩目的是，我们可以在一台12核的机器上非常容易地解决一个具有三个周期的波浪的三维拉普拉斯问题，达到四舍五入的精度--对于24m DoFs的第二至最大的情况，总共使用约3.5GB的内存，花费不超过8秒。最大的案例使用了30GB的内存，有1.91亿个DoFs。

@code
Running with 12 MPI processes, element FE_DGQHermite<3>(8)


Cycle 0
Number of degrees of freedom: 5832
Total setup time              0.0210681 s
Time solve (15 iterations)    0.0956945 s
Verification via L2 error:    0.0297194


Cycle 1
Number of degrees of freedom: 46656
Total setup time              0.0452428 s
Time solve (15 iterations)    0.113827 s
Verification via L2 error:    9.55733e-05


Cycle 2
Number of degrees of freedom: 373248
Total setup time              0.190423 s
Time solve (15 iterations)    0.218309 s
Verification via L2 error:    2.6868e-07


Cycle 3
Number of degrees of freedom: 2985984
Total setup time              0.627914 s
Time solve (15 iterations)    1.0595 s
Verification via L2 error:    4.6918e-10


Cycle 4
Number of degrees of freedom: 23887872
Total setup time              2.85215 s
Time solve (15 iterations)    8.30576 s
Verification via L2 error:    9.38583e-13


Cycle 5
Number of degrees of freedom: 191102976
Total setup time              16.1324 s
Time solve (15 iterations)    65.57 s
Verification via L2 error:    3.17875e-13
@endcode



<a name="Comparisonofefficiencyatdifferentpolynomialdegrees"></a><h3>Comparison of efficiency at different polynomial degrees</h3>


在介绍和代码中的评论中，多次提到用FEEvaluation和FEFaceEvaluation评估器可以非常有效地处理高阶。现在，我们想通过观察三维多网格求解器在不同多项式程度下的吞吐量来证实这些说法。我们收集的时间如下。我们首先在问题大小接近一千万的情况下运行求解器，如表前四行所示，并记录时间。然后，我们通过记录每秒解决的百万自由度数（MDoFs/s）来规范吞吐量，以便能够比较不同程度的效率，计算方法是自由度数除以求解器时间。

 <table align="center" class="doxtable">
  <tr>
   <th>degree</th>
   <th>1</th>
   <th>2</th>
   <th>3</th>
   <th>4</th>
   <th>5</th>
   <th>6</th>
   <th>7</th>
   <th>8</th>
   <th>9</th>
   <th>10</th>
   <th>11</th>
   <th>12</th>
  </tr>
  <tr>
   <th>Number of DoFs</th>
   <td>2097152</td>
   <td>7077888</td>
   <td>16777216</td>
   <td>32768000</td>
   <td>7077888</td>
   <td>11239424</td>
   <td>16777216</td>
   <td>23887872</td>
   <td>32768000</td>
   <td>43614208</td>
   <td>7077888</td>
   <td>8998912</td>
  </tr>
  <tr>
   <th>Number of iterations</th>
   <td>13</td>
   <td>12</td>
   <td>12</td>
   <td>12</td>
   <td>13</td>
   <td>13</td>
   <td>15</td>
   <td>15</td>
   <td>17</td>
   <td>19</td>
   <td>18</td>
   <td>18</td>
  </tr>
  <tr>
   <th>Solver time [s]</th>
   <td>0.713</td>
   <td>2.150</td>
   <td>4.638</td>
   <td>8.803</td>
   <td>2.041</td>
   <td>3.295</td>
   <td>5.723</td>
   <td>8.306</td>
   <td>12.75</td>
   <td>19.25</td>
   <td>3.530</td>
   <td>4.814</td>
  </tr>
  <tr>
   <th>MDoFs/s</th>
   <td>2.94</td>
   <td>3.29</td>
   <td>3.62</td>
   <td>3.72</td>
   <td>3.47</td>
   <td>3.41</td>
   <td>2.93</td>
   <td>2.88</td>
   <td>2.57</td>
   <td>2.27</td>
   <td>2.01</td>
   <td>1.87</td>
  </tr>
</table> 

我们清楚地看到每个DoF的效率最初是如何提高的，直到它达到多项式度数的最大值  $k=4$  。这种效果是令人惊讶的，不仅是因为较高的多项式度数通常会产生一个好得多的解决方案，而且特别是当考虑到基于矩阵的方案时，在较高的度数下更密集的耦合会导致单调的吞吐量下降（在3D中是一个急剧的下降， $k=4$ 比 $k=1$ 慢十倍以上！）。对于更高的度数，吞吐量减少了一些，这既是由于迭代次数的增加（从 $k=2,3,4$ 的12次增加到 $k=10$ 的19次），也是由于运算符评估的 $\mathcal O(k)$ 复杂性。尽管如此，对于更高的多项式度数，作为求解时间的效率仍然会更好，因为它们有更好的收敛率（至少对于像这个问题一样简单的问题）。对于 $k=12$ ，我们在100万个DoFs的情况下已经达到了舍入精度（求解时间不到一秒），而对于 $k=8$ ，我们需要2400万个DoFs和8秒钟。对于 $k=5$ ，误差约为 $10^{-9}$ 的5700万个DoFs，因此，尽管花了16秒，但离舍入还很远。

请注意，上述数字有点悲观，因为它们包括切比雪夫平滑器计算特征值估计的时间，这大约是求解器时间的10%。如果系统被多次求解（例如在流体力学中很常见），这个特征值的成本只需支付一次，更快的时间就可以得到。

<a name="Evaluationofefficiencyofingredients"></a><h3>Evaluation of efficiency of ingredients</h3>


最后，我们来看一下本教程程序中提出的一些特殊成分，即特别是FE_DGQHermite基础和 MatrixFree::DataAccessOnFaces. 的规格。 在下面的表格中，第三行显示了上面的优化求解器，第四行显示了只将 MatrixFree::DataAccessOnFaces 设置为 "未指定 "而不是最优的 "网格"，最后一行是用基本的FE_DGQ元素代替FE_DGQHermite，其中MPI交换更加昂贵，由 FEFaceEvaluation::gather_evaluate() 和 FEFaceEvaluation::integrate_scatter(). 完成的操作

 <table align="center" class="doxtable">
  <tr>
   <th>degree</th>
   <th>1</th>
   <th>2</th>
   <th>3</th>
   <th>4</th>
   <th>5</th>
   <th>6</th>
   <th>7</th>
   <th>8</th>
   <th>9</th>
   <th>10</th>
   <th>11</th>
   <th>12</th>
  </tr>
  <tr>
   <th>Number of DoFs</th>
   <td>2097152</td>
   <td>7077888</td>
   <td>16777216</td>
   <td>32768000</td>
   <td>7077888</td>
   <td>11239424</td>
   <td>16777216</td>
   <td>23887872</td>
   <td>32768000</td>
   <td>43614208</td>
   <td>7077888</td>
   <td>8998912</td>
  </tr>
  <tr>
   <th>Solver time optimized as in tutorial [s]</th>
   <td>0.713</td>
   <td>2.150</td>
   <td>4.638</td>
   <td>8.803</td>
   <td>2.041</td>
   <td>3.295</td>
   <td>5.723</td>
   <td>8.306</td>
   <td>12.75</td>
   <td>19.25</td>
   <td>3.530</td>
   <td>4.814</td>
  </tr>
  <tr>
   <th>Solver time MatrixFree::DataAccessOnFaces::unspecified [s]</th>
   <td>0.711</td>
   <td>2.151</td>
   <td>4.675</td>
   <td>8.968</td>
   <td>2.243</td>
   <td>3.655</td>
   <td>6.277</td>
   <td>9.082</td>
   <td>13.50</td>
   <td>20.05</td>
   <td>3.817</td>
   <td>5.178</td>
  </tr>
  <tr>
   <th>Solver time FE_DGQ [s]</th>
   <td>0.712</td>
   <td>2.041</td>
   <td>5.066</td>
   <td>9.335</td>
   <td>2.379</td>
   <td>3.802</td>
   <td>6.564</td>
   <td>9.714</td>
   <td>14.54</td>
   <td>22.76</td>
   <td>4.148</td>
   <td>5.857</td>
  </tr>
</table> 

表中的数据显示，不使用 MatrixFree::DataAccessOnFaces 对于较高的多项式度数来说，成本增加了10%左右。对于较低的度数，差异显然没有那么明显，因为体积与表面的比例更有利，需要交换的数据更少。如果只看矩阵-向量乘积，而不是这里显示的完整的多网格求解器，差异会更大，仅仅因为MPI通信，时间就会差20%左右。

对于 $k=1$ 和 $k=2$ ，类似Hermite的基函数显然没有真正得到回报（事实上，对于 $k=1$ ，多项式与FE_DGQ完全相同），其结果与FE_DGQ基函数相似。然而，对于从3开始的度数，我们看到FE_DGQHermite的优势越来越大，显示了这些基函数的有效性。

<a name="Possibilitiesforextension"></a><h3>Possibilities for extension</h3>


正如介绍中提到的，快速对角线化方法与具有恒定系数的直角坐标网相联系。如果我们想解决可变系数的问题，我们就需要在平滑参数的设计上投入更多的时间，选择适当的泛函（例如，在最近的箱形元素上近似反演）。

另一种扩展程序的方式是包括对自适应网格的支持，对于这种支持，在不同细化水平的边缘进行接口操作是必要的，如步骤39所讨论的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-59.cc"
*/
