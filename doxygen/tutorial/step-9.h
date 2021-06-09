/**
@page step_9 The step-9 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Discretizingtheadvectionequation">Discretizing the advection equation</a>
        <li><a href="#Whyisthismethodcalledstreamlinediffusion">Why is this method called "streamline diffusion"?</a>
        <li><a href="#WhyisthismethodalsocalledPetrovGalerkin">Why is this method also called "Petrov-Galerkin"?</a>
        <li><a href="#Whyisthismethodalsocalledstreamlineupwind">Why is this method also called "streamline-upwind"?</a>
        <li><a href="#Solvingthelinearsystemthatcorrespondstotheadvectionequation">Solving the linear system that corresponds to the advection equation</a>
        <li><a href="#Thetestcase">The test case</a>
        <li><a href="#Asimplerefinementcriterion">A simple refinement criterion</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdatadeclaration">Equation data declaration</a>
        <li><a href="#AdvectionProblemclassdeclaration">AdvectionProblem class declaration</a>
        <li><a href="#GradientEstimationclassdeclaration">GradientEstimation class declaration</a>
        <li><a href="#AdvectionProblemclassimplementation">AdvectionProblem class implementation</a>
        <li><a href="#GradientEstimationclassimplementation">GradientEstimation class implementation</a>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-9/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



在这个例子中，我们的目的如下。<ol>  <li>  解决平流方程  $\beta \cdot \nabla u = f$  ；  <li>  显示如果我们有一台多处理器机器，我们如何使用多线程来快速获得结果；  <li>  开发一个简单的细化准则。   </ol> 虽然第二个目的在不参考代码的情况下很难进行一般性描述，但我们将在下文中讨论其他两个目的。然后将在程序中的相关地方详细介绍多线程的使用。然而，我们将遵循 @ref threads "多处理器访问共享内存的并行计算 "文件模块中详述的关于WorkStream方法的一般讨论。




<a name="Discretizingtheadvectionequation"></a><h3>Discretizing the advection equation</h3>


在本例程序中，我们要对平流方程的解进行数值近似计算

@f[
  \beta \cdot \nabla u = f,


@f]

其中 $\beta$ 是描述平流方向和速度的矢量场（如果 $\beta=\beta(\mathbf x)$ ，它可能取决于空间变量）， $f$ 是一个源函数， $u$ 是解。该方程描述的物理过程是一个给定的流场 $\beta$ ，另一种物质随其流动，其密度或浓度由 $u$ 给出。该方程不包含这第二种物质在其载体物质内的扩散，但有源项。

很明显，在流入地，上述方程需要用边界条件来增加。

@f[
  u = g \qquad\qquad \mathrm{on}\ \partial\Omega_-,


@f]

其中 $\partial\Omega_-$ 描述了边界的流入部分，正式定义为

@f[
  \partial\Omega_-
  =
  \{{\mathbf x}\in \partial\Omega: \beta\cdot{\mathbf n}({\mathbf x}) < 0\},


@f]

和 ${\mathbf n}({\mathbf x})$ 是点 ${\mathbf x}\in\partial\Omega$ 处的域的向外法线。这个定义非常直观，因为由于 ${\mathbf n}$ 指向外侧，如果传输方向 $\beta$ 指向内侧，即在流入边界，那么与 $\beta$ 的标量乘积只能是负数。数学理论规定，我们不能在边界的流出部分提出任何边界条件。

不幸的是，上述方程不能以稳定的方式用标准的有限元方法进行求解。问题是，这个方程的解在垂直于传输方向上具有不充分的规则性：虽然它们沿 "风场" $\beta$ 定义的流线是平滑的，但它们在垂直于这个方向上可能是不连续的。这很容易理解：方程 $\beta \cdot
\nabla u = f$ 的意思实质上就是<i>rate of change of $u$ in
direction $\beta$ equals $f$</i>。但该方程对垂直方向的导数没有影响，因此，如果 $u$ 在流入边界的某一点上不连续，那么这个不连续将简单地沿着从这个边界点开始的风场流线传输。这些不连续会导致数值不稳定，使标准的连续有限元离散化不可能获得稳定的解。

解决这一困难的标准方法是  <em>  "流线-上风Petrov-Galerkin"  </em>  (SUPG)方法，有时也称为流线扩散法。对该方法的良好解释可以在  @cite elman2005  中找到。从形式上看，这种方法取代了我们从强形式推导出微分方程的弱形式的步骤。我们不是用测试函数 $v$ 乘以方程并在域上积分，而是乘以 $v + \delta \beta\cdot\nabla v$  ，其中 $\delta$ 是在（局部）网格宽度范围内选择的参数 $h$ ；通过设置 $\delta=0.1h$ 通常可以获得良好的结果。 为什么这被称为 "流线扩散 "将在下面解释；目前，让我们简单地认为这是我们如何推导出稳定的离散公式。这里 $\delta$ 的值小到足以使我们不引入过度的扩散，但大到足以使所产生的问题得到良好的解决。

使用上面定义的测试函数，问题的初始弱形式将要求找到一个函数 $u_h$ ，以便对于所有测试函数 $v_h$ ，我们有

@f[
  (\beta \cdot \nabla u_h, v_h + \delta \beta\cdot\nabla v_h)_\Omega
  =
  (f, v_h + \delta \beta\cdot\nabla v_h)_\Omega.


@f]

然而，我们希望将流入的边界条件 $u=g$ 弱化到这个问题中，这可以通过要求除了上述方程之外，我们还有

@f[
  (u_h, w_h)_{\partial\Omega_-}
  =
  (g, w_h)_{\partial\Omega_-}


@f]

为所有住在边界上的测试函数 $w_h$ ，它们来自一个合适的测试空间。事实证明，一个合适的检验函数空间恰好是 $\beta\cdot {\mathbf n}$ 乘以我们已经用于域内微分方程的检验空间中的函数 $v_h$ 的踪迹。因此，我们要求对于所有测试函数 $v_h$ ，我们有

@f[
  (u_h, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}
  =
  (g, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}.


@f]

在不试图进行论证的情况下（请再次参阅关于一般的有限元方法，特别是流线扩散方法的文献），我们可以将微分方程和边界值的方程结合在我们的稳定化问题的以下弱表述中：找到一个离散函数 $u_h$ ，使得对于所有离散测试函数 $v_h$ ，存在着

@f[
  (\beta \cdot \nabla u_h, v_h + \delta \beta\cdot\nabla v_h)_\Omega


  -
  (u_h, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}
  =
  (f, v_h + \delta \beta\cdot\nabla v_h)_\Omega


  -
  (g, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}.


@f]




人们会认为，这将导致一个系统矩阵被倒置，其形式为

@f[
  a_{ij} =
  (\beta \cdot \nabla \varphi_i,
   \varphi_j + \delta \beta\cdot\nabla \varphi_j)_\Omega


  -
  (\varphi_i, \beta\cdot {\mathbf n} \varphi_j)_{\partial\Omega_-},


@f]

与基函数  $\varphi_i,\varphi_j$  。  然而，这是每一个数值分析师至少会遇到一次的陷阱（包括作者）：我们在这里扩大了解决方案  $u_h = \sum_i U_i \varphi_i$  ，但如果我们这样做，我们将不得不解决这个问题

@f[
  U^T A = F^T,


@f]

其中 $U$ 是膨胀系数的向量，也就是说，我们必须解决我们可能天真的预期的转置问题。

这是我们在步骤3的介绍中提出的一个观点。在那里，我们认为为了避免这种问题，应该养成总是与检验函数<i>from the left</i>相乘的习惯，而不是从右边得到正确的矩阵。为了得到我们所需要的线性系统的形式，最好将弱式改写为

@f[
  (v_h + \delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega


  -
  (\beta\cdot {\mathbf n} v_h, u_h)_{\partial\Omega_-}
  =
  (v_h + \delta \beta\cdot\nabla v_h, f)_\Omega


  -
  (\beta\cdot {\mathbf n} v_h, g)_{\partial\Omega_-}


@f]

然后得到

@f[
  a_{ij} =
  (\varphi_i + \delta \beta \cdot \nabla \varphi_i,
   \beta\cdot\nabla \varphi_j)_\Omega


  -
  (\beta\cdot {\mathbf n} \varphi_i, \varphi_j)_{\partial\Omega_-},


@f]

作为系统矩阵。我们将在程序中组装这个矩阵。




<a name="Whyisthismethodcalledstreamlinediffusion"></a><h3>Why is this method called "streamline diffusion"?</h3>


看一下上面提到的双线性形式，我们看到离散解必须满足一个方程，其中弱形式的左手边有一个域项，即

@f[
  (v_h + \delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega,


@f]

或如果我们将其拆分，则形式为

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega
  +
  (\delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega.


@f]

如果我们想看看这将对应于什么强的方程形式，我们需要对第二项进行积分。这就产生了下面的公式，为了简单起见，我们暂时不考虑边界项。

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega


  -
  \left(v_h, \delta \nabla \cdot \left[\beta \left(\beta \cdot \nabla
  u_h\right)\right]\right)_\Omega
  +
  \text{boundary terms}.


@f]

让我们暂时假设风场 $\beta$ 是无发散的，即 $\nabla \cdot \beta = 0$  。然后将乘积法则应用于右边方括号内项的导数，并利用发散-绿色，我们将得到以下结果。

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega


  -
  \left(v_h, \delta \left[\beta \cdot \nabla\right] \left[\beta \cdot \nabla
  \right]u_h\right)_\Omega
  +
  \text{boundary terms}.


@f]

这意味着，方程的强势形式将是这样的

@f[
  \beta \cdot \nabla u_h


  -
  \delta
  \left[\beta \cdot \nabla\right] \left[\beta \cdot \nabla
  \right] u_h.


@f]

现在要认识到的是， $\beta\cdot\nabla$ 是 <em> 方向的导数 $\beta$   </em>  。因此，如果我们用 $\beta\cdot\nabla=\frac{\partial}{\partial \beta}$ 来表示（就像我们经常用 $\mathbf n\cdot\nabla=\frac{\partial}{\partial n}$ 来表示边界处法线方向的导数一样），那么方程的强形式是

@f[
  \beta \cdot \nabla u_h


  -
  \delta
  \frac{\partial^2}{\partial\beta^2} u_h.


@f]

换句话说，测试函数的不寻常选择相当于在强形式中增加了一个项，它对应于风场方向的二阶（即扩散）微分算子 $\beta$ ，即 "流线方向"。更全面的说明还必须探讨测试函数对边界值的影响，以及为什么有必要对右手边也使用相同的测试函数，但上面的讨论可能会使人明白该方法的 "流线扩散 "名称的由来。




<a name="WhyisthismethodalsocalledPetrovGalerkin"></a><h3>Why is this method also called "Petrov-Galerkin"?</h3>


"Galerkin方法 "是指通过将方程乘以测试函数 $v$ （然后在 $\Omega$ 上进行积分）来获得弱表述，其中函数 $v$ 与解 $u$ 来自同一空间（尽管可能具有不同的边界值）。但这并不是严格意义上的必要条件。我们也可以想象从不同的函数集中选择测试函数，只要这个不同的函数集具有与原始函数集 "同样多的维度"，这样我们最终就会有与自由度同样多的独立方程（在无限维的情况下，所有这些都需要适当地定义）。利用这种可能性的方法（即以不同的方式选择测试函数集和解决方案集）被称为 "Petrov-Galerkin "方法。在目前的情况下，测试函数的形式都是 $v+\beta\cdot\nabla v$ ，其中 $v$ 是来自解集。




<a name="Whyisthismethodalsocalledstreamlineupwind"></a><h3>Why is this method also called "streamline-upwind"?</h3>


[上风方法](https://en.wikipedia.org/wiki/Upwind_scheme)在推导平流方程的稳定方案方面有很长的历史。一般来说，这个想法是，我们不是在 "这里 "看一个函数，而是在 "上游 "或 "上风 "的一小段距离看它，也就是说，"这里 "的信息最初来自哪里。这可能意味着不考虑 $u(\mathbf x)$ ，而是像 $u(\mathbf x - \delta \beta)$ 这样的。或者，等于是在整合后，我们可以评估 $u(\mathbf x)$ ，而考虑 $v$ 的下游部分。   $v(\mathbf x+\delta \beta)$  .由于各种原因，这将是很麻烦的。首先，如果 $\mathbf x + \delta \beta$ 恰好在 $\Omega$ 之外，我们将不得不定义 $v$ 应该是什么；其次，数值计算积分将更加困难，因为我们不再在同一正交点评估 $u$ 和 $v$ 。但由于我们假设 $\delta$ 很小，我们可以做一个泰勒扩展。

@f[
  v(\mathbf x + \delta \beta)
  \approx
  v(\mathbf x) + \delta \beta \cdot \nabla v(\mathbf x).


@f]

这个测试函数的形式现在看起来应该很熟悉。




<a name="Solvingthelinearsystemthatcorrespondstotheadvectionequation"></a><h3>Solving the linear system that corresponds to the advection equation</h3>


由于产生的矩阵不再是对称正定的，我们不能使用通常的共轭梯度方法（在SolverCG类中实现）来解决这个系统。相反，我们使用GMRES（Generalized Minimum RESidual）方法（在SolverGMRES中实现），它适用于我们这里的那种问题。




<a name="Thetestcase"></a><h3>The test case</h3>


对于我们将在本教程程序中解决的问题，我们使用以下领域和函数（在 $d=2$ 空间维度）。

@f{eqnarray*}
  \Omega &=& [-1,1]^d \\
  \beta({\mathbf x})
  &=&
  \left(
    \begin{array}{c}2 \\ 1+\frac 45 \sin(8\pi x)\end{array}
  \right),
  \\
  s
  &=&
  0.1,
  \\
  f({\mathbf x})
  &=&
  \left\{
    \begin{array}{ll}
        \frac 1{10 s^d} &
        \mathrm{for}\ |{\mathbf x}-{\mathbf x}_0|<s, \\
        0 & \mathrm{else},
    \end{array}
  \right.
  \qquad\qquad
  {\mathbf x}_0
  =
  \left(
    \begin{array}{c} -\frac 34 \\ -\frac 34\end{array}
  \right),
  \\
  g
  &=&
  e^{5 (1 - |\mathbf x|^2)} \sin(16\pi|\mathbf x|^2).


@f}

对于 $d>2$ ，我们扩展了 $\beta$ 和 ${\mathbf x}_0$ ，只是将上面显示的最后一个组件再重复一次。

说了这么多，下面的评论是有必要的。<ol>  <li>  平流场 $\beta$ 大致以对角线方向从左下角向右上角传输解决方案，但叠加了一个摆动结构。   <li>  右手边在流入边界条件产生的场中加入了左下角的一个圆球，然后沿着这个圆球传输。   <li>  流入边界条件施加了一个加权的正弦结构，该结构与流场一起被传输。由于 $|{\mathbf x}|\ge 1$ 在边界上，加权项从未变得非常大。   </ol> 




<a name="Asimplerefinementcriterion"></a><h3>A simple refinement criterion</h3>


在以前所有的自适应细化的例子中，我们都使用了Kelly等人首先开发的误差估计器，它为每个单元 $K$ 分配了以下指标。

@f[
  \eta_K =
  \left(
    \frac {h_K}{24}
    \int_{\partial K}
      [\partial_n u_h]^2 \; d\sigma
  \right)^{1/2},


@f]

其中 $[\partial n u_h]$ 表示跨越单元格 $K$ 的一个面的法向导数的跳变。可以证明，这个误差指标使用的是二阶导数的离散类似物，由单元大小的一个幂加权，这个幂是根据这里假定使用的线性元素调整的。

@f[
  \eta_K \approx
  C h \| \nabla^2 u \|_K,


@f]

其本身与能量准则中的误差大小有关。

在目前的情况下，这个误差指标的问题是，它假定精确的解拥有二阶导数。在某些情况下，这对于拉普拉斯问题的解来说已经是个问题了，尽管那里大多数问题允许在 $H^2$ 中求解。如果解只在 $H^1$ 中，那么二阶导数在域的某些部分（低维）是奇异的，在网格细化的情况下，误差指标不会减少。因此，该算法将不断细化这些部分周围的单元，即细化为点或线（在2d中）。

然而，对于目前的情况，解通常不在 $H^1$ 中（而且这种缺失的规律性并不是像拉普拉斯方程那样的特殊情况），所以上述的误差指标并不真正适用。因此，我们将开发一个基于梯度的离散近似的指标。尽管梯度经常不存在，但这是我们唯一可用的标准，至少在我们使用连续元素时是如此。首先，我们注意到，给定两个单元 $K$ ， $K'$ ，其中心由矢量 ${\mathbf y}_{KK'}$ 连接，我们可以对一个函数 $u$ 的方向导数进行近似，如下所示。

@f[
  \frac{{\mathbf y}_{KK'}^T}{|{\mathbf y}_{KK'}|} \nabla u
  \approx
  \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|},


@f]

其中 $u(K)$ 和 $u(K')$ 表示 $u$ 在各自单元格中心的评价。现在我们将上述近似值乘以 ${\mathbf y}_{KK'}/|{\mathbf y}_{KK'}|$ ，并对 $K$ 的所有邻居 $K'$ 求和。

@f[
  \underbrace{
    \left(\sum_{K'} \frac{{\mathbf y}_{KK'} {\mathbf y}_{KK'}^T}
                         {|{\mathbf y}_{KK'}|^2}\right)}_{=:Y}
  \nabla u
  \approx
  \sum_{K'}
  \frac{{\mathbf y}_{KK'}}{|{\mathbf y}_{KK'}|}
  \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|}.


@f]

如果连接 ${\mathbf y}_{KK'}$ 与邻居的向量 $K$ 横跨整个空间（即大致为： $K$ 在所有方向都有邻居），那么左侧表达式中括号内的项形成一个正则矩阵，我们可以将其反转，得到 $u$ 对 $K$ 的梯度的近似。

@f[
  \nabla u
  \approx
  Y^{-1}
  \left(
    \sum_{K'}
    \frac{{\mathbf y}_{KK'}}{|{\mathbf y}_{KK'}|}
    \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|}
  \right).


@f]

我们将用 $\nabla_h u(K)$ 表示右手边的近似值，我们将使用以下数量作为细化标准。

@f[
  \eta_K = h^{1+d/2} |\nabla_h u_h(K)|,


@f]

这是受以下（不严谨的）论证的启发。

@f{eqnarray*}
  \|u-u_h\|^2_{L_2}
  &\le&
  C h^2 \|\nabla u\|^2_{L_2}
\\
  &\approx&
  C
  \sum_K
  h_K^2 \|\nabla u\|^2_{L_2(K)}
\\
  &\le&
  C
  \sum_K
  h_K^2 h_K^d \|\nabla u\|^2_{L_\infty(K)}
\\
  &\approx&
  C
  \sum_K
  h_K^{2+d} |\nabla_h u_h(K)|^2


@f}




 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * Just as in previous examples, we have to include several files of which the
 * meaning has already been discussed:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/grid/grid_out.h>
 * 
 * @endcode
 * 
 * The following two files provide classes and information for multithreaded
 * programs. In the first one, the classes and functions are declared which we
 * need to do assembly in parallel (i.e. the
 * <code>WorkStream</code> namespace). The
 * second file has a class MultithreadInfo which can be used to query the
 * number of processors in your system, which is often useful when deciding
 * how many threads to start in parallel.
 * 
 * @code
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/base/multithread_info.h>
 * 
 * @endcode
 * 
 * The next new include file declares a base class <code>TensorFunction</code>
 * not unlike the <code>Function</code> class, but with the difference that
 * TensorFunction::value returns a Tensor instead of a scalar.
 * 
 * @code
 * #include <deal.II/base/tensor_function.h>
 * 
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * This is C++, as we want to write some output to disk:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * The last step is as in previous programs:
 * 
 * @code
 * namespace Step9
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdatadeclaration"></a> 
 * <h3>Equation data declaration</h3>
 * 

 * 
 * Next we declare a class that describes the advection field. This, of
 * course, is a vector field with as many components as there are space
 * dimensions. One could now use a class derived from the
 * <code>Function</code> base class, as we have done for boundary values and
 * coefficients in previous examples, but there is another possibility in
 * the library, namely a base class that describes tensor valued
 * functions. This is more convenient than overriding Function::value() with
 * a method that knows about multiple function components: at the end of the
 * day we need a Tensor, so we may as well just use a class that returns a
 * Tensor.
 * 
 * @code
 *   template <int dim>
 *   class AdvectionField : public TensorFunction<1, dim>
 *   {
 *   public:
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override;
 * 
 * @endcode
 * 
 * In previous examples, we have used assertions that throw exceptions in
 * several places. However, we have never seen how such exceptions are
 * declared. This can be done as follows:
 * 
 * @code
 *     DeclException2(ExcDimensionMismatch,
 *                    unsigned int,
 *                    unsigned int,
 *                    << "The vector has size " << arg1 << " but should have "
 *                    << arg2 << " elements.");
 * @endcode
 * 
 * The syntax may look a little strange, but is reasonable. The format is
 * basically as follows: use the name of one of the macros
 * <code>DeclExceptionN</code>, where <code>N</code> denotes the number of
 * additional parameters which the exception object shall take. In this
 * case, as we want to throw the exception when the sizes of two vectors
 * differ, we need two arguments, so we use
 * <code>DeclException2</code>. The first parameter then describes the
 * name of the exception, while the following declare the data types of
 * the parameters. The last argument is a sequence of output directives
 * that will be piped into the <code>std::cerr</code> object, thus the
 * strange format with the leading <code>@<@<</code> operator and the
 * like. Note that we can access the parameters which are passed to the
 * exception upon construction (i.e. within the <code>Assert</code> call)
 * by using the names <code>arg1</code> through <code>argN</code>, where
 * <code>N</code> is the number of arguments as defined by the use of the
 * respective macro <code>DeclExceptionN</code>.
 *     

 * 
 * To learn how the preprocessor expands this macro into actual code,
 * please refer to the documentation of the exception classes. In brief,
 * this macro call declares and defines a class
 * <code>ExcDimensionMismatch</code> inheriting from ExceptionBase which
 * implements all necessary error output functions.
 * 
 * @code
 *   };
 * 
 * @endcode
 * 
 * The following two functions implement the interface described above. The
 * first simply implements the function as described in the introduction,
 * while the second uses the same trick to avoid calling a virtual function
 * as has already been introduced in the previous example program. Note the
 * check for the right sizes of the arguments in the second function, which
 * should always be present in such functions; it is our experience that
 * many if not most programming errors result from incorrectly initialized
 * arrays, incompatible parameters to functions and the like; using
 * assertion as in this case can eliminate many of these problems.
 * 
 * @code
 *   template <int dim>
 *   Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const
 *   {
 *     Tensor<1, dim> value;
 *     value[0] = 2;
 *     for (unsigned int i = 1; i < dim; ++i)
 *       value[i] = 1 + 0.8 * std::sin(8. * numbers::PI * p[0]);
 * 
 *     return value;
 *   }
 * 
 * @endcode
 * 
 * Besides the advection field, we need two functions describing the source
 * terms (<code>right hand side</code>) and the boundary values. As
 * described in the introduction, the source is a constant function in the
 * vicinity of a source point, which we denote by the constant static
 * variable <code>center_point</code>. We set the values of this center
 * using the same template tricks as we have shown in the step-7 example
 * program. The rest is simple and has been shown previously.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 *     static const Point<dim> center_point;
 *   };
 * 
 * 
 *   template <>
 *   const Point<1> RightHandSide<1>::center_point = Point<1>(-0.75);
 * 
 *   template <>
 *   const Point<2> RightHandSide<2>::center_point = Point<2>(-0.75, -0.75);
 * 
 *   template <>
 *   const Point<3> RightHandSide<3>::center_point = Point<3>(-0.75, -0.75, -0.75);
 * 
 * 
 * 
 * @endcode
 * 
 * The only new thing here is that we check for the value of the
 * <code>component</code> parameter. As this is a scalar function, it is
 * obvious that it only makes sense if the desired component has the index
 * zero, so we assert that this is indeed the
 * case. <code>ExcIndexRange</code> is a global predefined exception
 * (probably the one most often used, we therefore made it global instead of
 * local to some class), that takes three parameters: the index that is
 * outside the allowed range, the first element of the valid range and the
 * one past the last (i.e. again the half-open interval so often used in the
 * C++ standard library):
 * 
 * @code
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 *     (void)component;
 *     Assert(component == 0, ExcIndexRange(component, 0, 1));
 *     const double diameter = 0.1;
 *     return ((p - center_point).norm_square() < diameter * diameter ?
 *               0.1 / std::pow(diameter, dim) :
 *               0.0);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Finally for the boundary values, which is just another class derived from
 * the <code>Function</code> base class:
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
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> & p,
 *                                     const unsigned int component) const
 *   {
 *     (void)component;
 *     Assert(component == 0, ExcIndexRange(component, 0, 1));
 * 
 *     const double sine_term = std::sin(16. * numbers::PI * p.norm_square());
 *     const double weight    = std::exp(5. * (1. - p.norm_square()));
 *     return weight * sine_term;
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="AdvectionProblemclassdeclaration"></a> 
 * <h3>AdvectionProblem class declaration</h3>
 * 

 * 
 * Here comes the main class of this program. It is very much like the main
 * classes of previous examples, so we again only comment on the
 * differences.
 * 
 * @code
 *   template <int dim>
 *   class AdvectionProblem
 *   {
 *   public:
 *     AdvectionProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 * 
 * @endcode
 * 
 * The next set of functions will be used to assemble the
 * matrix. However, unlike in the previous examples, the
 * <code>assemble_system()</code> function will not do the work
 * itself, but rather will delegate the actual assembly to helper
 * functions <code>assemble_local_system()</code> and
 * <code>copy_local_to_global()</code>. The rationale is that
 * matrix assembly can be parallelized quite well, as the
 * computation of the local contributions on each cell is entirely
 * independent of other cells, and we only have to synchronize
 * when we add the contribution of a cell to the global
 * matrix.
 *     

 * 
 * The strategy for parallelization we choose here is one of the
 * possibilities mentioned in detail in the @ref threads module in
 * the documentation. Specifically, we will use the WorkStream
 * approach discussed there. Since there is so much documentation
 * in this module, we will not repeat the rationale for the design
 * choices here (for example, if you read through the module
 * mentioned above, you will understand what the purpose of the
 * <code>AssemblyScratchData</code> and
 * <code>AssemblyCopyData</code> structures is). Rather, we will
 * only discuss the specific implementation.
 *     

 * 
 * If you read the page mentioned above, you will find that in
 * order to parallelize assembly, we need two data structures --
 * one that corresponds to data that we need during local
 * integration ("scratch data", i.e., things we only need as
 * temporary storage), and one that carries information from the
 * local integration to the function that then adds the local
 * contributions to the corresponding elements of the global
 * matrix. The former of these typically contains the FEValues and
 * FEFaceValues objects, whereas the latter has the local matrix,
 * local right hand side, and information about which degrees of
 * freedom live on the cell for which we are assembling a local
 * contribution. With this information, the following should be
 * relatively self-explanatory:
 * 
 * @code
 *     struct AssemblyScratchData
 *     {
 *       AssemblyScratchData(const FiniteElement<dim> &fe);
 *       AssemblyScratchData(const AssemblyScratchData &scratch_data);
 * 
 * @endcode
 * 
 * FEValues and FEFaceValues are expensive objects to set up, so we
 * include them in the scratch object so that as much data is reused
 * between cells as possible.
 * 
 * @code
 *       FEValues<dim>     fe_values;
 *       FEFaceValues<dim> fe_face_values;
 * 
 * @endcode
 * 
 * We also store a few vectors that we will populate with values on each
 * cell. Setting these objects up is, in the usual case, cheap; however,
 * they require memory allocations, which can be expensive in
 * multithreaded applications. Hence we keep them here so that
 * computations on a cell do not require new allocations.
 * 
 * @code
 *       std::vector<double>         rhs_values;
 *       std::vector<Tensor<1, dim>> advection_directions;
 *       std::vector<double>         face_boundary_values;
 *       std::vector<Tensor<1, dim>> face_advection_directions;
 * 
 * @endcode
 * 
 * Finally, we need objects that describe the problem's data:
 * 
 * @code
 *       AdvectionField<dim> advection_field;
 *       RightHandSide<dim>  right_hand_side;
 *       BoundaryValues<dim> boundary_values;
 *     };
 * 
 *     struct AssemblyCopyData
 *     {
 *       FullMatrix<double>                   cell_matrix;
 *       Vector<double>                       cell_rhs;
 *       std::vector<types::global_dof_index> local_dof_indices;
 *     };
 * 
 *     void assemble_system();
 *     void local_assemble_system(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       AssemblyScratchData &                                 scratch,
 *       AssemblyCopyData &                                    copy_data);
 *     void copy_local_to_global(const AssemblyCopyData &copy_data);
 * 
 * 
 * @endcode
 * 
 * The following functions again are the same as they were in previous
 * examples, as are the subsequent variables:
 * 
 * @code
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     FE_Q<dim> fe;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="GradientEstimationclassdeclaration"></a> 
 * <h3>GradientEstimation class declaration</h3>
 * 

 * 
 * Now, finally, here comes the class that will compute the difference
 * approximation of the gradient on each cell and weighs that with a power
 * of the mesh size, as described in the introduction. This class is a
 * simple version of the <code>DerivativeApproximation</code> class in the
 * library, that uses similar techniques to obtain finite difference
 * approximations of the gradient of a finite element field, or of higher
 * derivatives.
 *   

 * 
 * The class has one public static function <code>estimate</code> that is
 * called to compute a vector of error indicators, and a few private functions
 * that do the actual work on all active cells. As in other parts of the
 * library, we follow an informal convention to use vectors of floats for
 * error indicators rather than the common vectors of doubles, as the
 * additional accuracy is not necessary for estimated values.
 *   

 * 
 * In addition to these two functions, the class declares two exceptions
 * which are raised when a cell has no neighbors in each of the space
 * directions (in which case the matrix described in the introduction would
 * be singular and can't be inverted), while the other one is used in the
 * more common case of invalid parameters to a function, namely a vector of
 * wrong size.
 *   

 * 
 * Two other comments: first, the class has no non-static member functions
 * or variables, so this is not really a class, but rather serves the
 * purpose of a <code>namespace</code> in C++. The reason that we chose a
 * class over a namespace is that this way we can declare functions that are
 * private. This can be done with namespaces as well, if one declares some
 * functions in header files in the namespace and implements these and other
 * functions in the implementation file. The functions not declared in the
 * header file are still in the namespace but are not callable from
 * outside. However, as we have only one file here, it is not possible to
 * hide functions in the present case.
 *   

 * 
 * The second comment is that the dimension template parameter is attached
 * to the function rather than to the class itself. This way, you don't have
 * to specify the template parameter yourself as in most other cases, but
 * the compiler can figure its value out itself from the dimension of the
 * DoFHandler object that one passes as first argument.
 *   

 * 
 * Before jumping into the fray with the implementation, let us also comment
 * on the parallelization strategy. We have already introduced the necessary
 * framework for using the WorkStream concept in the declaration of the main
 * class of this program above. We will use it again here. In the current
 * context, this means that we have to define
 * <ol>
 * <li>classes for scratch and copy objects,</li>
 * <li>a function that does the local computation on one cell, and</li>
 * <li>a function that copies the local result into a global object.</li>
 * </ol>
 * Given this general framework, we will, however, deviate from it a
 * bit. In particular, WorkStream was generally invented for cases where
 * each local computation on a cell <i>adds</i> to a global object -- for
 * example, when assembling linear systems where we add local contributions
 * into a global matrix and right hand side. WorkStream is designed to handle
 * the potential conflict of multiple threads trying to do this addition at
 * the same time, and consequently has to provide for some way to ensure that
 * only one thread gets to do this at a time. Here, however, the situation is
 * slightly different: we compute contributions from every cell
 * individually, but then all we need to do is put them into an element of
 * an output vector that is unique to each cell. Consequently, there is no
 * risk that the write operations from two cells might conflict, and the
 * elaborate machinery of WorkStream to avoid conflicting writes is not
 * necessary. Consequently, what we will do is this: We still need a scratch
 * object that holds, for example, the FEValues object. However, we only
 * create a fake, empty copy data structure. Likewise, we do need the
 * function that computes local contributions, but since it can already put
 * the result into its final location, we do not need a copy-local-to-global
 * function and will instead give the WorkStream::run() function an empty
 * function object -- the equivalent to a NULL function pointer.
 * 
 * @code
 *   class GradientEstimation
 *   {
 *   public:
 *     template <int dim>
 *     static void estimate(const DoFHandler<dim> &dof,
 *                          const Vector<double> & solution,
 *                          Vector<float> &        error_per_cell);
 * 
 *     DeclException2(ExcInvalidVectorLength,
 *                    int,
 *                    int,
 *                    << "Vector has length " << arg1 << ", but should have "
 *                    << arg2);
 *     DeclException0(ExcInsufficientDirections);
 * 
 *   private:
 *     template <int dim>
 *     struct EstimateScratchData
 *     {
 *       EstimateScratchData(const FiniteElement<dim> &fe,
 *                           const Vector<double> &    solution,
 *                           Vector<float> &           error_per_cell);
 *       EstimateScratchData(const EstimateScratchData &data);
 * 
 *       FEValues<dim> fe_midpoint_value;
 *       std::vector<typename DoFHandler<dim>::active_cell_iterator>
 *         active_neighbors;
 * 
 *       const Vector<double> &solution;
 *       Vector<float> &       error_per_cell;
 * 
 *       std::vector<double> cell_midpoint_value;
 *       std::vector<double> neighbor_midpoint_value;
 *     };
 * 
 *     struct EstimateCopyData
 *     {};
 * 
 *     template <int dim>
 *     static void
 *     estimate_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                   EstimateScratchData<dim> &scratch_data,
 *                   const EstimateCopyData &  copy_data);
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="AdvectionProblemclassimplementation"></a> 
 * <h3>AdvectionProblem class implementation</h3>
 * 

 * 
 * 

 * 
 * Now for the implementation of the main class. Constructor, destructor and
 * the function <code>setup_system</code> follow the same pattern that was
 * used previously, so we need not comment on these three function:
 * 
 * @code
 *   template <int dim>
 *   AdvectionProblem<dim>::AdvectionProblem()
 *     : dof_handler(triangulation)
 *     , fe(5)
 *   {}
 * 
 * 
 * 
 *   template <int dim>
 *   void AdvectionProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     hanging_node_constraints,
 *                                     /*keep_constrained_dofs =*/false);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In the following function, the matrix and right hand side are
 * assembled. As stated in the documentation of the main class above, it
 * does not do this itself, but rather delegates to the function following
 * next, utilizing the WorkStream concept discussed in @ref threads .
 *   

 * 
 * If you have looked through the @ref threads module, you will have
 * seen that assembling in parallel does not take an incredible
 * amount of extra code as long as you diligently describe what the
 * scratch and copy data objects are, and if you define suitable
 * functions for the local assembly and the copy operation from local
 * contributions to global objects. This done, the following will do
 * all the heavy lifting to get these operations done on multiple
 * threads on as many cores as you have in your system:
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::assemble_system()
 *   {
 *     WorkStream::run(dof_handler.begin_active(),
 *                     dof_handler.end(),
 *                     *this,
 *                     &AdvectionProblem::local_assemble_system,
 *                     &AdvectionProblem::copy_local_to_global,
 *                     AssemblyScratchData(fe),
 *                     AssemblyCopyData());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * As already mentioned above, we need to have scratch objects for
 * the parallel computation of local contributions. These objects
 * contain FEValues and FEFaceValues objects (as well as some arrays), and so
 * we will need to have constructors and copy constructors that allow us to
 * create them. For the cell terms we need the values
 * and gradients of the shape functions, the quadrature points in
 * order to determine the source density and the advection field at
 * a given point, and the weights of the quadrature points times the
 * determinant of the Jacobian at these points. In contrast, for the
 * boundary integrals, we don't need the gradients, but rather the
 * normal vectors to the cells. This determines which update flags
 * we will have to pass to the constructors of the members of the
 * class:
 * 
 * @code
 *   template <int dim>
 *   AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
 *     const FiniteElement<dim> &fe)
 *     : fe_values(fe,
 *                 QGauss<dim>(fe.degree + 1),
 *                 update_values | update_gradients | update_quadrature_points |
 *                   update_JxW_values)
 *     , fe_face_values(fe,
 *                      QGauss<dim - 1>(fe.degree + 1),
 *                      update_values | update_quadrature_points |
 *                        update_JxW_values | update_normal_vectors)
 *     , rhs_values(fe_values.get_quadrature().size())
 *     , advection_directions(fe_values.get_quadrature().size())
 *     , face_boundary_values(fe_face_values.get_quadrature().size())
 *     , face_advection_directions(fe_face_values.get_quadrature().size())
 *   {}
 * 
 * 
 * 
 *   template <int dim>
 *   AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData(
 *     const AssemblyScratchData &scratch_data)
 *     : fe_values(scratch_data.fe_values.get_fe(),
 *                 scratch_data.fe_values.get_quadrature(),
 *                 update_values | update_gradients | update_quadrature_points |
 *                   update_JxW_values)
 *     , fe_face_values(scratch_data.fe_face_values.get_fe(),
 *                      scratch_data.fe_face_values.get_quadrature(),
 *                      update_values | update_quadrature_points |
 *                        update_JxW_values | update_normal_vectors)
 *     , rhs_values(scratch_data.rhs_values.size())
 *     , advection_directions(scratch_data.advection_directions.size())
 *     , face_boundary_values(scratch_data.face_boundary_values.size())
 *     , face_advection_directions(scratch_data.face_advection_directions.size())
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * Now, this is the function that does the actual work. It is not very
 * different from the <code>assemble_system</code> functions of previous
 * example programs, so we will again only comment on the differences. The
 * mathematical stuff closely follows what we have said in the introduction.
 *   

 * 
 * There are a number of points worth mentioning here, though. The
 * first one is that we have moved the FEValues and FEFaceValues
 * objects into the ScratchData object. We have done so because the
 * alternative would have been to simply create one every time we
 * get into this function -- i.e., on every cell. It now turns out
 * that the FEValues classes were written with the explicit goal of
 * moving everything that remains the same from cell to cell into
 * the construction of the object, and only do as little work as
 * possible in FEValues::reinit() whenever we move to a new
 * cell. What this means is that it would be very expensive to
 * create a new object of this kind in this function as we would
 * have to do it for every cell -- exactly the thing we wanted to
 * avoid with the FEValues class. Instead, what we do is create it
 * only once (or a small number of times) in the scratch objects and
 * then re-use it as often as we can.
 *   

 * 
 * This begs the question of whether there are other objects we
 * create in this function whose creation is expensive compared to
 * its use. Indeed, at the top of the function, we declare all sorts
 * of objects. The <code>AdvectionField</code>,
 * <code>RightHandSide</code> and <code>BoundaryValues</code> do not
 * cost much to create, so there is no harm here. However,
 * allocating memory in creating the <code>rhs_values</code> and
 * similar variables below typically costs a significant amount of
 * time, compared to just accessing the (temporary) values we store
 * in them. Consequently, these would be candidates for moving into
 * the <code>AssemblyScratchData</code> class. We will leave this as
 * an exercise.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::local_assemble_system(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     AssemblyScratchData &                                 scratch_data,
 *     AssemblyCopyData &                                    copy_data)
 *   {
 * @endcode
 * 
 * We define some abbreviations to avoid unnecessarily long lines:
 * 
 * @code
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points =
 *       scratch_data.fe_values.get_quadrature().size();
 *     const unsigned int n_face_q_points =
 *       scratch_data.fe_face_values.get_quadrature().size();
 * 
 * @endcode
 * 
 * We declare cell matrix and cell right hand side...
 * 
 * @code
 *     copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 *     copy_data.cell_rhs.reinit(dofs_per_cell);
 * 
 * @endcode
 * 
 * ... an array to hold the global indices of the degrees of freedom of
 * the cell on which we are presently working...
 * 
 * @code
 *     copy_data.local_dof_indices.resize(dofs_per_cell);
 * 
 * @endcode
 * 
 * ... then initialize the <code>FEValues</code> object...
 * 
 * @code
 *     scratch_data.fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * ... obtain the values of right hand side and advection directions
 * at the quadrature points...
 * 
 * @code
 *     scratch_data.advection_field.value_list(
 *       scratch_data.fe_values.get_quadrature_points(),
 *       scratch_data.advection_directions);
 *     scratch_data.right_hand_side.value_list(
 *       scratch_data.fe_values.get_quadrature_points(), scratch_data.rhs_values);
 * 
 * @endcode
 * 
 * ... set the value of the streamline diffusion parameter as
 * described in the introduction...
 * 
 * @code
 *     const double delta = 0.1 * cell->diameter();
 * 
 * @endcode
 * 
 * ... and assemble the local contributions to the system matrix and
 * right hand side as also discussed above:
 * 
 * @code
 *     for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *         {
 * @endcode
 * 
 * Alias the AssemblyScratchData object to keep the lines from
 * getting too long:
 * 
 * @code
 *           const auto &sd = scratch_data;
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             copy_data.cell_matrix(i, j) +=
 *               ((sd.fe_values.shape_value(i, q_point) +           // (phi_i +
 *                 delta * (sd.advection_directions[q_point] *      // delta beta
 *                          sd.fe_values.shape_grad(i, q_point))) * // grad phi_i)
 *                sd.advection_directions[q_point] *                // beta
 *                sd.fe_values.shape_grad(j, q_point)) *            // grad phi_j
 *               sd.fe_values.JxW(q_point);                         // dx
 * 
 *           copy_data.cell_rhs(i) +=
 *             (sd.fe_values.shape_value(i, q_point) +           // (phi_i +
 *              delta * (sd.advection_directions[q_point] *      // delta beta
 *                       sd.fe_values.shape_grad(i, q_point))) * // grad phi_i)
 *             sd.rhs_values[q_point] *                          // f
 *             sd.fe_values.JxW(q_point);                        // dx
 *         }
 * 
 * @endcode
 * 
 * Besides the cell terms which we have built up now, the bilinear
 * form of the present problem also contains terms on the boundary of
 * the domain. Therefore, we have to check whether any of the faces of
 * this cell are on the boundary of the domain, and if so assemble the
 * contributions of this face as well. Of course, the bilinear form
 * only contains contributions from the <code>inflow</code> part of
 * the boundary, but to find out whether a certain part of a face of
 * the present cell is part of the inflow boundary, we have to have
 * information on the exact location of the quadrature points and on
 * the direction of flow at this point; we obtain this information
 * using the FEFaceValues object and only decide within the main loop
 * whether a quadrature point is on the inflow boundary.
 * 
 * @code
 *     for (const auto &face : cell->face_iterators())
 *       if (face->at_boundary())
 *         {
 * @endcode
 * 
 * Ok, this face of the present cell is on the boundary of the
 * domain. Just as for the usual FEValues object which we have
 * used in previous examples and also above, we have to
 * reinitialize the FEFaceValues object for the present face:
 * 
 * @code
 *           scratch_data.fe_face_values.reinit(cell, face);
 * 
 * @endcode
 * 
 * For the quadrature points at hand, we ask for the values of
 * the inflow function and for the direction of flow:
 * 
 * @code
 *           scratch_data.boundary_values.value_list(
 *             scratch_data.fe_face_values.get_quadrature_points(),
 *             scratch_data.face_boundary_values);
 *           scratch_data.advection_field.value_list(
 *             scratch_data.fe_face_values.get_quadrature_points(),
 *             scratch_data.face_advection_directions);
 * 
 * @endcode
 * 
 * Now loop over all quadrature points and see whether this face is on
 * the inflow or outflow part of the boundary. The normal
 * vector points out of the cell: since the face is at
 * the boundary, the normal vector points out of the domain,
 * so if the advection direction points into the domain, its
 * scalar product with the normal vector must be negative (to see why
 * this is true, consider the scalar product definition that uses a
 * cosine):
 * 
 * @code
 *           for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
 *             if (scratch_data.fe_face_values.normal_vector(q_point) *
 *                   scratch_data.face_advection_directions[q_point] <
 *                 0.)
 * @endcode
 * 
 * If the face is part of the inflow boundary, then compute the
 * contributions of this face to the global matrix and right
 * hand side, using the values obtained from the
 * FEFaceValues object and the formulae discussed in the
 * introduction:
 * 
 * @code
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     copy_data.cell_matrix(i, j) -=
 *                       (scratch_data.face_advection_directions[q_point] *
 *                        scratch_data.fe_face_values.normal_vector(q_point) *
 *                        scratch_data.fe_face_values.shape_value(i, q_point) *
 *                        scratch_data.fe_face_values.shape_value(j, q_point) *
 *                        scratch_data.fe_face_values.JxW(q_point));
 * 
 *                   copy_data.cell_rhs(i) -=
 *                     (scratch_data.face_advection_directions[q_point] *
 *                      scratch_data.fe_face_values.normal_vector(q_point) *
 *                      scratch_data.face_boundary_values[q_point] *
 *                      scratch_data.fe_face_values.shape_value(i, q_point) *
 *                      scratch_data.fe_face_values.JxW(q_point));
 *                 }
 *         }
 * 
 * @endcode
 * 
 * The final piece of information the copy routine needs is the global
 * indices of the degrees of freedom on this cell, so we end by writing
 * them to the local array:
 * 
 * @code
 *     cell->get_dof_indices(copy_data.local_dof_indices);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The second function we needed to write was the one that copies
 * the local contributions the previous function computed (and
 * put into the AssemblyCopyData object) into the global matrix and right
 * hand side vector objects. This is essentially what we always had
 * as the last block of code when assembling something on every
 * cell. The following should therefore be pretty obvious:
 * 
 * @code
 *   template <int dim>
 *   void
 *   AdvectionProblem<dim>::copy_local_to_global(const AssemblyCopyData &copy_data)
 *   {
 *     hanging_node_constraints.distribute_local_to_global(
 *       copy_data.cell_matrix,
 *       copy_data.cell_rhs,
 *       copy_data.local_dof_indices,
 *       system_matrix,
 *       system_rhs);
 *   }
 * 
 * @endcode
 * 
 * Here comes the linear solver routine. As the system is no longer
 * symmetric positive definite as in all the previous examples, we cannot
 * use the Conjugate Gradient method anymore. Rather, we use a solver that
 * is more general and does not rely on any special properties of the
 * matrix: the GMRES method. GMRES, like the conjugate gradient method,
 * requires a decent preconditioner: we use a Jacobi preconditioner here,
 * which works well enough for this problem.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::solve()
 *   {
 *     SolverControl               solver_control(std::max<std::size_t>(1000,
 *                                                        system_rhs.size() / 10),
 *                                  1e-10 * system_rhs.l2_norm());
 *     SolverGMRES<Vector<double>> solver(solver_control);
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.0);
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     Vector<double> residual(dof_handler.n_dofs());
 * 
 *     system_matrix.vmult(residual, solution);
 *     residual -= system_rhs;
 *     std::cout << "   Iterations required for convergence: "
 *               << solver_control.last_step() << '\n'
 *               << "   Max norm of residual:                "
 *               << residual.linfty_norm() << '\n';
 * 
 *     hanging_node_constraints.distribute(solution);
 *   }
 * 
 * @endcode
 * 
 * The following function refines the grid according to the quantity
 * described in the introduction. The respective computations are made in
 * the class <code>GradientEstimation</code>.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::refine_grid()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     GradientEstimation::estimate(dof_handler,
 *                                  solution,
 *                                  estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * @endcode
 * 
 * This function is similar to the one in step 6, but since we use a higher
 * degree finite element we save the solution in a different
 * way. Visualization programs like VisIt and Paraview typically only
 * understand data that is associated with nodes: they cannot plot
 * fifth-degree basis functions, which results in a very inaccurate picture
 * of the solution we computed. To get around this we save multiple
 * <em>patches</em> per cell: in 2D we save 64 bilinear `cells' to the VTU
 * file for each cell, and in 3D we save 512. The end result is that the
 * visualization program will use a piecewise linear interpolation of the
 * cubic basis functions: this captures the solution detail and, with most
 * screen resolutions, looks smooth. We save the grid in a separate step
 * with no extra patches so that we have a visual representation of the cell
 * faces.
 *   

 * 
 * Version 9.1 of deal.II gained the ability to write higher degree
 * polynomials (i.e., write piecewise bicubic visualization data for our
 * piecewise bicubic solution) VTK and VTU output: however, not all recent
 * versions of ParaView and VisIt (as of 2018) can read this format, so we
 * use the older, more general (but less efficient) approach here.
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     {
 *       GridOut       grid_out;
 *       std::ofstream output("grid-" + std::to_string(cycle) + ".vtu");
 *       grid_out.write_vtu(triangulation, output);
 *     }
 * 
 *     {
 *       DataOut<dim> data_out;
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "solution");
 *       data_out.build_patches(8);
 * 
 * @endcode
 * 
 * VTU output can be expensive, both to compute and to write to
 * disk. Here we ask ZLib, a compression library, to compress the data
 * in a way that maximizes throughput.
 * 
 * @code
 *       DataOutBase::VtkFlags vtk_flags;
 *       vtk_flags.compression_level =
 *         DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
 *       data_out.set_flags(vtk_flags);
 * 
 *       std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
 *       data_out.write_vtu(output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * ... as is the main loop (setup -- solve -- refine), aside from the number
 * of cycles and the initial grid:
 * 
 * @code
 *   template <int dim>
 *   void AdvectionProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 10; ++cycle)
 *       {
 *         std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, -1, 1);
 *             triangulation.refine_global(3);
 *           }
 *         else
 *           {
 *             refine_grid();
 *           }
 * 
 * 
 *         std::cout << "   Number of active cells:              "
 *                   << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         std::cout << "   Number of degrees of freedom:        "
 *                   << dof_handler.n_dofs() << std::endl;
 * 
 *         assemble_system();
 *         solve();
 *         output_results(cycle);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="GradientEstimationclassimplementation"></a> 
 * <h3>GradientEstimation class implementation</h3>
 * 

 * 
 * Now for the implementation of the <code>GradientEstimation</code> class.
 * Let us start by defining constructors for the
 * <code>EstimateScratchData</code> class used by the
 * <code>estimate_cell()</code> function:
 * 
 * @code
 *   template <int dim>
 *   GradientEstimation::EstimateScratchData<dim>::EstimateScratchData(
 *     const FiniteElement<dim> &fe,
 *     const Vector<double> &    solution,
 *     Vector<float> &           error_per_cell)
 *     : fe_midpoint_value(fe,
 *                         QMidpoint<dim>(),
 *                         update_values | update_quadrature_points)
 *     , solution(solution)
 *     , error_per_cell(error_per_cell)
 *     , cell_midpoint_value(1)
 *     , neighbor_midpoint_value(1)
 *   {
 * @endcode
 * 
 * We allocate a vector to hold iterators to all active neighbors of
 * a cell. We reserve the maximal number of active neighbors in order to
 * avoid later reallocations. Note how this maximal number of active
 * neighbors is computed here.
 * 
 * @code
 *     active_neighbors.reserve(GeometryInfo<dim>::faces_per_cell *
 *                              GeometryInfo<dim>::max_children_per_face);
 *   }
 * 
 * 
 *   template <int dim>
 *   GradientEstimation::EstimateScratchData<dim>::EstimateScratchData(
 *     const EstimateScratchData &scratch_data)
 *     : fe_midpoint_value(scratch_data.fe_midpoint_value.get_fe(),
 *                         scratch_data.fe_midpoint_value.get_quadrature(),
 *                         update_values | update_quadrature_points)
 *     , solution(scratch_data.solution)
 *     , error_per_cell(scratch_data.error_per_cell)
 *     , cell_midpoint_value(1)
 *     , neighbor_midpoint_value(1)
 *   {}
 * 
 * 
 * @endcode
 * 
 * Next comes the implementation of the <code>GradientEstimation</code>
 * class. The first function does not much except for delegating work to the
 * other function, but there is a bit of setup at the top.
 *   

 * 
 * Before starting with the work, we check that the vector into
 * which the results are written has the right size. Programming
 * mistakes in which one forgets to size arguments correctly at the
 * calling site are quite common. Because the resulting damage from
 * not catching such errors is often subtle (e.g., corruption of
 * data somewhere in memory, or non-reproducible results), it is
 * well worth the effort to check for such things.
 * 
 * @code
 *   template <int dim>
 *   void GradientEstimation::estimate(const DoFHandler<dim> &dof_handler,
 *                                     const Vector<double> & solution,
 *                                     Vector<float> &        error_per_cell)
 *   {
 *     Assert(
 *       error_per_cell.size() == dof_handler.get_triangulation().n_active_cells(),
 *       ExcInvalidVectorLength(error_per_cell.size(),
 *                              dof_handler.get_triangulation().n_active_cells()));
 * 
 *     WorkStream::run(dof_handler.begin_active(),
 *                     dof_handler.end(),
 *                     &GradientEstimation::template estimate_cell<dim>,
 *                     std::function<void(const EstimateCopyData &)>(),
 *                     EstimateScratchData<dim>(dof_handler.get_fe(),
 *                                              solution,
 *                                              error_per_cell),
 *                     EstimateCopyData());
 *   }
 * 
 * 
 * @endcode
 * 
 * Here comes the function that estimates the local error by computing the
 * finite difference approximation of the gradient. The function first
 * computes the list of active neighbors of the present cell and then
 * computes the quantities described in the introduction for each of
 * the neighbors. The reason for this order is that it is not a one-liner
 * to find a given neighbor with locally refined meshes. In principle, an
 * optimized implementation would find neighbors and the quantities
 * depending on them in one step, rather than first building a list of
 * neighbors and in a second step their contributions but we will gladly
 * leave this as an exercise. As discussed before, the worker function
 * passed to WorkStream::run works on "scratch" objects that keep all
 * temporary objects. This way, we do not need to create and initialize
 * objects that are expensive to initialize within the function that does
 * the work every time it is called for a given cell. Such an argument is
 * passed as the second argument. The third argument would be a "copy-data"
 * object (see @ref threads for more information) but we do not actually use
 * any of these here. Since WorkStream::run() insists on passing three
 * arguments, we declare this function with three arguments, but simply
 * ignore the last one.
 *   

 * 
 * (This is unsatisfactory from an aesthetic perspective. It can be avoided
 * by using an anonymous (lambda) function. If you allow, let us here show
 * how. First, assume that we had declared this function to only take two
 * arguments by omitting the unused last one. Now, WorkStream::run still
 * wants to call this function with three arguments, so we need to find a
 * way to "forget" the third argument in the call. Simply passing
 * WorkStream::run the pointer to the function as we do above will not do
 * this -- the compiler will complain that a function declared to have two
 * arguments is called with three arguments. However, we can do this by
 * passing the following as the third argument to WorkStream::run():
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * [](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *    EstimateScratchData<dim> &                            scratch_data,
 *    EstimateCopyData &)
 * {
 *   GradientEstimation::estimate_cell<dim>(cell, scratch_data);
 * }
 * @endcode
 * </div>
 * This is not much better than the solution implemented below: either the
 * routine itself must take three arguments or it must be wrapped by
 * something that takes three arguments. We don't use this since adding the
 * unused argument at the beginning is simpler.
 *   

 * 
 * Now for the details:
 * 
 * @code
 *   template <int dim>
 *   void GradientEstimation::estimate_cell(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     EstimateScratchData<dim> &                            scratch_data,
 *     const EstimateCopyData &)
 *   {
 * @endcode
 * 
 * We need space for the tensor <code>Y</code>, which is the sum of
 * outer products of the y-vectors.
 * 
 * @code
 *     Tensor<2, dim> Y;
 * 
 * @endcode
 * 
 * First initialize the <code>FEValues</code> object, as well as the
 * <code>Y</code> tensor:
 * 
 * @code
 *     scratch_data.fe_midpoint_value.reinit(cell);
 * 
 * @endcode
 * 
 * Now, before we go on, we first compute a list of all active neighbors
 * of the present cell. We do so by first looping over all faces and see
 * whether the neighbor there is active, which would be the case if it
 * is on the same level as the present cell or one level coarser (note
 * that a neighbor can only be once coarser than the present cell, as
 * we only allow a maximal difference of one refinement over a face in
 * deal.II). Alternatively, the neighbor could be on the same level
 * and be further refined; then we have to find which of its children
 * are next to the present cell and select these (note that if a child
 * of a neighbor of an active cell that is next to this active cell,
 * needs necessarily be active itself, due to the one-refinement rule
 * cited above).
 *     

 * 
 * Things are slightly different in one space dimension, as there the
 * one-refinement rule does not exist: neighboring active cells may
 * differ in as many refinement levels as they like. In this case, the
 * computation becomes a little more difficult, but we will explain
 * this below.
 *     

 * 
 * Before starting the loop over all neighbors of the present cell, we
 * have to clear the array storing the iterators to the active
 * neighbors, of course.
 * 
 * @code
 *     scratch_data.active_neighbors.clear();
 *     for (const auto face_n : cell->face_indices())
 *       if (!cell->at_boundary(face_n))
 *         {
 * @endcode
 * 
 * First define an abbreviation for the iterator to the face and
 * the neighbor
 * 
 * @code
 *           const auto face     = cell->face(face_n);
 *           const auto neighbor = cell->neighbor(face_n);
 * 
 * @endcode
 * 
 * Then check whether the neighbor is active. If it is, then it
 * is on the same level or one level coarser (if we are not in
 * 1D), and we are interested in it in any case.
 * 
 * @code
 *           if (neighbor->is_active())
 *             scratch_data.active_neighbors.push_back(neighbor);
 *           else
 *             {
 * @endcode
 * 
 * If the neighbor is not active, then check its children.
 * 
 * @code
 *               if (dim == 1)
 *                 {
 * @endcode
 * 
 * To find the child of the neighbor which bounds to the
 * present cell, successively go to its right child if
 * we are left of the present cell (n==0), or go to the
 * left child if we are on the right (n==1), until we
 * find an active cell.
 * 
 * @code
 *                   auto neighbor_child = neighbor;
 *                   while (neighbor_child->has_children())
 *                     neighbor_child = neighbor_child->child(face_n == 0 ? 1 : 0);
 * 
 * @endcode
 * 
 * As this used some non-trivial geometrical intuition,
 * we might want to check whether we did it right,
 * i.e., check whether the neighbor of the cell we found
 * is indeed the cell we are presently working
 * on. Checks like this are often useful and have
 * frequently uncovered errors both in algorithms like
 * the line above (where it is simple to involuntarily
 * exchange <code>n==1</code> for <code>n==0</code> or
 * the like) and in the library (the assumptions
 * underlying the algorithm above could either be wrong,
 * wrongly documented, or are violated due to an error
 * in the library). One could in principle remove such
 * checks after the program works for some time, but it
 * might be a good things to leave it in anyway to check
 * for changes in the library or in the algorithm above.
 *                   

 * 
 * Note that if this check fails, then this is certainly
 * an error that is irrecoverable and probably qualifies
 * as an internal error. We therefore use a predefined
 * exception class to throw here.
 * 
 * @code
 *                   Assert(neighbor_child->neighbor(face_n == 0 ? 1 : 0) == cell,
 *                          ExcInternalError());
 * 
 * @endcode
 * 
 * If the check succeeded, we push the active neighbor
 * we just found to the stack we keep:
 * 
 * @code
 *                   scratch_data.active_neighbors.push_back(neighbor_child);
 *                 }
 *               else
 * @endcode
 * 
 * If we are not in 1d, we collect all neighbor children
 * `behind' the subfaces of the current face and move on:
 * 
 * @code
 *                 for (unsigned int subface_n = 0; subface_n < face->n_children();
 *                      ++subface_n)
 *                   scratch_data.active_neighbors.push_back(
 *                     cell->neighbor_child_on_subface(face_n, subface_n));
 *             }
 *         }
 * 
 * @endcode
 * 
 * OK, now that we have all the neighbors, lets start the computation
 * on each of them. First we do some preliminaries: find out about the
 * center of the present cell and the solution at this point. The
 * latter is obtained as a vector of function values at the quadrature
 * points, of which there are only one, of course. Likewise, the
 * position of the center is the position of the first (and only)
 * quadrature point in real space.
 * 
 * @code
 *     const Point<dim> this_center =
 *       scratch_data.fe_midpoint_value.quadrature_point(0);
 * 
 *     scratch_data.fe_midpoint_value.get_function_values(
 *       scratch_data.solution, scratch_data.cell_midpoint_value);
 * 
 * @endcode
 * 
 * Now loop over all active neighbors and collect the data we
 * need.
 * 
 * @code
 *     Tensor<1, dim> projected_gradient;
 *     for (const auto &neighbor : scratch_data.active_neighbors)
 *       {
 * @endcode
 * 
 * Then get the center of the neighbor cell and the value of the
 * finite element function at that point. Note that for this
 * information we have to reinitialize the <code>FEValues</code>
 * object for the neighbor cell.
 * 
 * @code
 *         scratch_data.fe_midpoint_value.reinit(neighbor);
 *         const Point<dim> neighbor_center =
 *           scratch_data.fe_midpoint_value.quadrature_point(0);
 * 
 *         scratch_data.fe_midpoint_value.get_function_values(
 *           scratch_data.solution, scratch_data.neighbor_midpoint_value);
 * 
 * @endcode
 * 
 * Compute the vector <code>y</code> connecting the centers of the
 * two cells. Note that as opposed to the introduction, we denote
 * by <code>y</code> the normalized difference vector, as this is
 * the quantity used everywhere in the computations.
 * 
 * @code
 *         Tensor<1, dim> y        = neighbor_center - this_center;
 *         const double   distance = y.norm();
 *         y /= distance;
 * 
 * @endcode
 * 
 * Then add up the contribution of this cell to the Y matrix...
 * 
 * @code
 *         for (unsigned int i = 0; i < dim; ++i)
 *           for (unsigned int j = 0; j < dim; ++j)
 *             Y[i][j] += y[i] * y[j];
 * 
 * @endcode
 * 
 * ... and update the sum of difference quotients:
 * 
 * @code
 *         projected_gradient += (scratch_data.neighbor_midpoint_value[0] -
 *                                scratch_data.cell_midpoint_value[0]) /
 *                               distance * y;
 *       }
 * 
 * @endcode
 * 
 * If now, after collecting all the information from the neighbors, we
 * can determine an approximation of the gradient for the present
 * cell, then we need to have passed over vectors <code>y</code> which
 * span the whole space, otherwise we would not have all components of
 * the gradient. This is indicated by the invertibility of the matrix.
 *     

 * 
 * If the matrix is not invertible, then the present
 * cell had an insufficient number of active neighbors. In contrast to
 * all previous cases (where we raised exceptions) this is, however,
 * not a programming error: it is a runtime error that can happen in
 * optimized mode even if it ran well in debug mode, so it is
 * reasonable to try to catch this error also in optimized mode. For
 * this case, there is the <code>AssertThrow</code> macro: it checks
 * the condition like the <code>Assert</code> macro, but not only in
 * debug mode; it then outputs an error message, but instead of
 * aborting the program as in the case of the <code>Assert</code>
 * macro, the exception is thrown using the <code>throw</code> command
 * of C++. This way, one has the possibility to catch this error and
 * take reasonable counter actions. One such measure would be to
 * refine the grid globally, as the case of insufficient directions
 * can not occur if every cell of the initial grid has been refined at
 * least once.
 * 
 * @code
 *     AssertThrow(determinant(Y) != 0, ExcInsufficientDirections());
 * 
 * @endcode
 * 
 * If, on the other hand, the matrix is invertible, then invert it,
 * multiply the other quantity with it, and compute the estimated error
 * using this quantity and the correct powers of the mesh width:
 * 
 * @code
 *     const Tensor<2, dim> Y_inverse = invert(Y);
 * 
 *     const Tensor<1, dim> gradient = Y_inverse * projected_gradient;
 * 
 * @endcode
 * 
 * The last part of this function is the one where we write into
 * the element of the output vector what we have just
 * computed. The address of this vector has been stored in the
 * scratch data object, and all we have to do is know how to get
 * at the correct element inside this vector -- but we can ask the
 * cell we're on the how-manyth active cell it is for this:
 * 
 * @code
 *     scratch_data.error_per_cell(cell->active_cell_index()) =
 *       (std::pow(cell->diameter(), 1 + 1.0 * dim / 2) * gradient.norm());
 *   }
 * } // namespace Step9
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Mainfunction"></a> 
 * <h3>Main function</h3>
 * 

 * 
 * The <code>main</code> function is similar to the previous examples. The
 * primary difference is that we use MultithreadInfo to set the maximum
 * number of threads (see the documentation module @ref threads
 * "Parallel computing with multiple processors accessing shared memory"
 * for more information). The number of threads used is the minimum of the
 * environment variable DEAL_II_NUM_THREADS and the parameter of
 * <code>set_thread_limit</code>. If no value is given to
 * <code>set_thread_limit</code>, the default value from the Intel Threading
 * Building Blocks (TBB) library is used. If the call to
 * <code>set_thread_limit</code> is omitted, the number of threads will be
 * chosen by TBB independently of DEAL_II_NUM_THREADS.
 * 
 * @code
 * int main()
 * {
 *   using namespace dealii;
 *   try
 *     {
 *       MultithreadInfo::set_thread_limit();
 * 
 *       Step9::AdvectionProblem<2> advection_problem_2d;
 *       advection_problem_2d.run();
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
examples/step-9/doc/results.dox



<a name="Results"></a><h1>Results</h1>



这个程序的结果并不特别引人注目。它们由控制台输出、一些网格文件和每个网格的解决方案组成。首先是控制台的输出。

@code
Cycle 0:
   Number of active cells:              64
   Number of degrees of freedom:        1681
   Iterations required for convergence: 298
   Max norm of residual:                3.60316e-12
Cycle 1:
   Number of active cells:              124
   Number of degrees of freedom:        3537
   Iterations required for convergence: 415
   Max norm of residual:                3.70682e-12
Cycle 2:
   Number of active cells:              247
   Number of degrees of freedom:        6734
   Iterations required for convergence: 543
   Max norm of residual:                7.19716e-13
Cycle 3:
   Number of active cells:              502
   Number of degrees of freedom:        14105
   Iterations required for convergence: 666
   Max norm of residual:                3.45628e-13
Cycle 4:
   Number of active cells:              1003
   Number of degrees of freedom:        27462
   Iterations required for convergence: 1064
   Max norm of residual:                1.86495e-13
Cycle 5:
   Number of active cells:              1993
   Number of degrees of freedom:        55044
   Iterations required for convergence: 1251
   Max norm of residual:                1.28765e-13
Cycle 6:
   Number of active cells:              3985
   Number of degrees of freedom:        108492
   Iterations required for convergence: 2035
   Max norm of residual:                6.78085e-14
Cycle 7:
   Number of active cells:              7747
   Number of degrees of freedom:        210612
   Iterations required for convergence: 2187
   Max norm of residual:                2.61457e-14
Cycle 8:
   Number of active cells:              15067
   Number of degrees of freedom:        406907
   Iterations required for convergence: 3079
   Max norm of residual:                2.9932e-14
Cycle 9:
   Number of active cells:              29341
   Number of degrees of freedom:        780591
   Iterations required for convergence: 3913
   Max norm of residual:                8.15689e-15
@endcode



相当多的单元格被用在最精细的层面上，以解决解决方案的特征。下面是第四和第十个网格。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-grid-3.png" alt="细化周期中的第四个网格，显示对特征的一些适应性。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-grid-9.png" alt="细化周期中的第十个网格，显示完全捕捉到波。" width="400" height="400"> </div> <div> 以及第四和第十的解决方案。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-3.png" alt="第四个解决方案，显示我们解决了大多数特征，但有些仍然没有解决，显得很模糊。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-9.png" alt="第十个解决方案，显示完全解决的流动。" width="400" height="400"> </div> </div> 以及网格和解决方案都被放大。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-3-zoom.png" alt="第四个解决方案的细节，显示我们解决了大多数特征，但有些仍然没有解决，显得模糊不清。特别是，较大的单元格需要细化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-9-zoom.png" alt="第十个解决方案的细节，显示我们需要比第四个解决方案中存在的更多单元格。" width="400" height="400"> </div> </div>

解决方案是由那部分沿摆动的平流场从左、下边界传送到右上方的部分，以及由左下角的源所产生的部分，其结果也是沿传送的。上面显示的网格很好地适应了解决这些特征。图中的比较表明，即使我们使用的是高阶近似，我们仍然需要自适应的网格细化来完全解决摆动。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-9.cc"
*/
