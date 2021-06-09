/**
@page step_21 The step-21 tutorial program
This tutorial depends on step-20.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetwophaseflowproblem">The two phase flow problem</a>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Spacediscretization">Space discretization</a>
        <li><a href="#Linearsolvers">Linear solvers</a>
        <li><a href="#Choosingatimestep">Choosing a time step</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeTwoPhaseFlowProblemcodeclass">The <code>TwoPhaseFlowProblem</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
      <ul>
        <li><a href="#Pressurerighthandside">Pressure right hand side</a>
        <li><a href="#Pressureboundaryvalues">Pressure boundary values</a>
        <li><a href="#Saturationboundaryvalues">Saturation boundary values</a>
        <li><a href="#Initialdata">Initial data</a>
      </ul>
        <li><a href="#Theinversepermeabilitytensor">The inverse permeability tensor</a>
      <ul>
        <li><a href="#Singlecurvingcrackpermeability">Single curving crack permeability</a>
        <li><a href="#Randommediumpermeability">Random medium permeability</a>
      </ul>
        <li><a href="#Theinversemobilityandsaturationfunctions">The inverse mobility and saturation functions</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
        <li><a href="#codeTwoPhaseFlowProblemcodeclassimplementation"><code>TwoPhaseFlowProblem</code> class implementation</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemTwoPhaseFlowProblem">TwoPhaseFlowProblem::TwoPhaseFlowProblem</a>
        <li><a href="#TwoPhaseFlowProblemmake_grid_and_dofs">TwoPhaseFlowProblem::make_grid_and_dofs</a>
        <li><a href="#TwoPhaseFlowProblemassemble_system">TwoPhaseFlowProblem::assemble_system</a>
        <li><a href="#TwoPhaseFlowProblemassemble_rhs_S">TwoPhaseFlowProblem::assemble_rhs_S</a>
        <li><a href="#TwoPhaseFlowProblemsolve">TwoPhaseFlowProblem::solve</a>
        <li><a href="#TwoPhaseFlowProblemoutput_results">TwoPhaseFlowProblem::output_results</a>
        <li><a href="#TwoPhaseFlowProblemproject_back_saturation">TwoPhaseFlowProblem::project_back_saturation</a>
        <li><a href="#TwoPhaseFlowProblemget_maximal_velocity">TwoPhaseFlowProblem::get_maximal_velocity</a>
        <li><a href="#TwoPhaseFlowProblemrun">TwoPhaseFlowProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Solvers">Solvers</a>
        <li><a href="#Timestepping">Time stepping</a>
        <li><a href="#Adaptivity">Adaptivity</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-21/doc/intro.dox

<a name="Introduction"></a>［<a name="Intro"></a>］ ［<h1>Introduction</h1>］


这个项目是由德克萨斯A&amp;M大学的李艳的一个学生项目发展而来的。这个项目的大部分工作都是由她完成的。

在这个项目中，我们提出了一个针对多孔介质中两相流动问题的数值模拟。这个问题包括一个椭圆方程和一个非线性的、随时间变化的传输方程。因此，这也是第一个时间相关的教程程序（除了有点奇怪的时间相关的 @ref step_18 "step-18"）。

这里涉及的方程是步骤20中已经涉及的材料的延伸。特别是，它们属于矢量值问题的范畴。这个主题的顶层概述可以在 @ref vector_valued 模块中找到。




<a name="Thetwophaseflowproblem"></a><h3>The two phase flow problem</h3>


多孔介质中两相流动的建模对于环境修复和石油及地下水库的管理都很重要。涉及两相流动的实际情况包括非水相液体在含水层中的分散，或混合液体（如油和水）在储层中的联合运动。仿真模型如果要提供真实的预测，就必须准确地考虑到这些影响。

为了推导出管理方程，考虑储层中的两相流动 $\Omega$ ，假设流体的运动由粘性效应主导；即我们忽略了重力、压缩性和毛细压力的影响。孔隙率将被认为是恒定的。我们将使用下标 $w$ 和 $o$ 来表示两相中的任何一个变量，即水和油的简称。然而，方程的推导对其他流体对也是适用的。

两相中每一相的分子移动的速度由达西定律决定，该定律指出速度与压力梯度成正比。

@f{eqnarray*}
  \mathbf{u}_{j}
  =


  -\frac{k_{rj}(S)}{\mu_{j}} \mathbf{K} \cdot \nabla p


@f}

其中 $\mathbf{u}_{j}$ 是相 $j=o,w$ 的速度， $K$ 是渗透率张量， $k_{rj}$ 是相 $j$ 的相对渗透率， $p$ 是压力， $\mu_{j}$ 是相 $j$ 的粘性。最后， $S$ 是饱和度（体积分数），即一个数值在0和1之间的函数，表示流体混合物的组成。一般来说，系数 $K, k_{rj}, \mu$ 可能是空间上的变量，在下文中我们将始终把它们作为非常数函数。

我们将达西定律与各相的质量守恒声明结合起来。

@f[
  \textrm{div}\ \mathbf{u}_{j} = q_j,


@f]

每个相都有一个源项。通过对两相求和，我们可以用所谓的压力方程来表达治理方程。

@f{eqnarray*}


- \nabla \cdot (\mathbf{K}\lambda(S) \nabla p)= q.


@f}

这里， $q$ 是和源项，而

@f[
  \lambda(S) = \frac{k_{rw}(S)}{\mu_{w}}+\frac{k_{ro}(S)}{\mu_{o}}


@f]

是总的流动性。

到目前为止，这看起来是一个普通的静止的、类似泊松的方程，我们可以用前几个教程的技术马上解决（例如，看一下步骤6，非常类似的东西）。然而，我们还没有说到饱和度，这当然会随着流体的移动而改变。

方程的第二部分是对饱和度的动态描述，即两种流体的相对浓度如何随时间变化。置换流体（水）的饱和度方程由以下守恒定律给出。

@f{eqnarray*}
  S_{t} + \nabla \cdot (F(S) \mathbf{u}) = q_{w},


@f}

这可以通过使用前一个方程中发散算子的乘积规则来重写。

@f{eqnarray*}
  S_{t} + F(S) \left[\nabla \cdot \mathbf{u}\right]
        + \mathbf{u} \cdot \left[ \nabla F(S)\right]
  = S_{t} + F(S) q + \mathbf{u} \cdot \nabla F(S) = q_{w}.


@f}

这里， $q=\nabla\cdot \mathbf{u}$ 是上面介绍的总流入量， $q_{w}$ 是置换流体（水）的流速。这两者与分流量 $F(S)$ 的关系如下。

@f[
  q_{w} = F(S) q,


@f]

其中分数流通常通过（启发式）表达式进行参数化

@f[
  F(S)
  =
  \frac{k_{rw}(S)/\mu_{w}}{k_{rw}(S)/\mu_{w} + k_{ro}(S)/\mu_{o}}.


@f]

将所有这些放在一起，可以得到以下形式的饱和度方程，即平流式。

@f{eqnarray*}
  S_{t} + \mathbf{u} \cdot \nabla F(S) = 0,


@f}

其中 $\mathbf u$ 是总速度

@f[
  \mathbf{u} =
  \mathbf{u}_{o} + \mathbf{u}_{w} = -\lambda(S) \mathbf{K}\cdot\nabla p.


@f]

注意，平流方程包含术语 $\mathbf{u} \cdot \nabla
F(S)$ 而不是 $\mathbf{u} \cdot \nabla S$ ，以表明饱和度不是简单地沿途传送；相反，由于两相以不同的速度移动，即使在平流坐标系中，饱和度实际上也可以改变。为了看到这一点，重写 $\mathbf{u} \cdot \nabla F(S)
= \mathbf{u} F'(S) \cdot \nabla S$ ，观察到具有饱和度 $S$ 的相的<i>actual</i>传输速度是 $\mathbf u F'(S)$ ，而另一相的传输速度是 $\mathbf u (1-F'(S))$  。  因此， $F(S)$ 通常被称为<i>fractional flow</i>。

综上所述，我们得到的是以下两个方程式。

@f{eqnarray*}


  - \nabla \cdot (\mathbf{K}\lambda(S) \nabla p) &=& q
  \qquad \textrm{in}\ \Omega\times[0,T],
  \\
  S_{t} + \mathbf{u} \cdot \nabla F(S) &=& 0
  \qquad \textrm{in}\ \Omega\times[0,T].


@f}

这里， $p=p(\mathbf x, t), S=S(\mathbf x, t)$ 现在是随时间变化的函数：虽然在每个时间瞬间，流场与压力处于平衡状态（即我们忽略了动态加速），但饱和度随着流动而运输，因此随时间变化，反过来又通过第一个方程对 $S$ 的依赖性影响流场。

这组方程有一个奇特的特点：两个方程中的一个有时间导数，另一个没有。这与压力和速度通过瞬时约束耦合的特点相对应，而饱和度在有限的时间尺度上演变。

这样的方程组被称为微分代数方程（DAE），因为其中一个方程是微分方程，另一个不是（至少不是相对于时间变量），因此是一个 "代数 "方程。这个符号来自常微分方程领域，在这个领域中，所有没有关于时间变量的导数的东西都必然是一个代数方程）。这类方程包含相当知名的情况：例如，与时间相关的斯托克斯和纳维-斯托克斯方程（其中代数约束是流场的发散， $\textrm{div}\ \mathbf u$ ，必须为零）以及与时间相关的麦克斯韦方程（这里，代数约束是电位移场的发散等于电荷密度， $\textrm{div}\ \mathbf D = \rho$  ，磁通密度的发散为零。   $\textrm{div}\ \mathbf
B = 0$ ）；即使是Step-18的准静态模型也属于这个类别。我们将看到，这两个方程的不同特征将告知我们这两个方程的离散化策略。




<a name="Timediscretization"></a><h3>Time discretization</h3>


在储层模拟界，通常是通过回到一阶混合公式来解决上面得出的方程。为此，我们重新引入总速度 $\mathbf u$ ，并将方程写成以下形式。

@f{eqnarray*}
  \mathbf{u}+\mathbf{K}\lambda(S) \nabla p&=&0 \\
  \nabla \cdot\mathbf{u} &=& q \\
  S_{t} + \mathbf{u} \cdot \nabla F(S) &=& 0.


@f}

这种表述方式还有一个好处，即我们不必将出现在传输方程中的总速度 $\mathbf u$ 表示为压力的函数，而是可以将其作为主变量。鉴于前两个方程的鞍点结构以及它们与我们在第20步中介绍的混合拉普拉斯公式的相似性，我们将再次使用混合离散化，这并不奇怪。

但是，让我们先把这个问题推迟一下。我们处理这些方程的第一件事是考虑时间离散化。在储层模拟中，有一个相当标准的算法，我们将在这里使用。它首先使用隐式方程解决压力问题，然后使用显式时间步进方案解决饱和问题。该算法被称为IMplicit Pressure Explicit Saturation（隐式压力显式饱和），很早以前就被提出：Sheldon等人在1959年提出，Stone和Gardner在1961年提出（J.W.Sheldon, B. Zondek and W. T. Cardwell:<i>One-dimensional, incompressible, non-capillary, two-phase
fluid flow in a porous medium</i>, Trans.SPE AIME, 216 (1959), pp. 290-296; H. L. Stone and A. O. Gardner Jr: <i>Analysis of gas-cap or dissolved-gas
reservoirs</i>, Trans.SPE AIME, 222 (1961), pp. 92-104)。在一个稍加修改的形式中，这个算法可以写成如下：对于每一个时间步长，解决

@f{eqnarray*}
  \mathbf{u}^{n+1}+\mathbf{K}\lambda(S^n) \nabla p^{n+1}&=&0 \\
  \nabla \cdot\mathbf{u}^{n+1} &=& q^{n+1} \\
  \frac {S^{n+1}-S^n}{\triangle t} + \mathbf{u}^{n+1} \cdot \nabla F(S^n) &=& 0,


@f}

其中 $\triangle t$ 是一个时间步长。请注意我们是如何解决隐式压力-速度系统的，它只取决于先前计算的饱和度 $S^n$ ，然后对 $S^{n+1}$ 做一个显式时间步长，它只取决于先前已知的 $S^n$ 和刚刚计算的 $\mathbf{u}^{n+1}$ 。这样一来，我们就不必像使用全隐式方法那样，对系统的非线性进行迭代。从更现代的角度来看，这应该被看作是一种 "算子分割 "方法。

然后我们可以将问题以弱的形式陈述如下，用测试函数 $\mathbf v$ 、 $\phi$ 和 $\sigma$ 乘以每个方程，并通过部分整合条款。

@f{eqnarray*}
  \left((\mathbf{K}\lambda(S^n))^{-1} \mathbf{u}^{n+1},\mathbf v\right)_\Omega -
  (p^{n+1}, \nabla\cdot\mathbf v)_\Omega &=&


  - (p^{n+1}, \mathbf v)_{\partial\Omega}
  \\
  (\nabla \cdot\mathbf{u}^{n+1}, \phi)_\Omega &=& (q^{n+1},\phi)_\Omega


@f}

注意，在第一项中，我们必须规定边界 $p^{n+1}$ 上的压力 $\partial\Omega$ 作为我们问题的边界值。   $\mathbf n$ 表示对 $\partial K$ 的单位外向法向量，如常。

对于饱和度方程，我们通过部分积分后得到

@f{eqnarray*}
  (S^{n+1}, \sigma)_\Omega


  -
  \triangle t
  \sum_K
  \left\{
  \left(F(S^n), \nabla \cdot (\mathbf{u}^{n+1} \sigma)\right)_K


  -
  \left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}, \sigma\right)_{\partial K}
  \right\}
  &=&
  (S^n,\sigma)_\Omega.


@f}

利用 $\nabla \cdot \mathbf{u}^{n+1}=q^{n+1}$ 这一事实，我们可以重写细胞项，得到如下方程。

@f{eqnarray*}
  (S^{n+1}, \sigma)_\Omega


  -
  \triangle t
  \sum_K
  \left\{
  \left(F(S^n) \mathbf{u}^{n+1}, \nabla \sigma\right)_K


  -
  \left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}), \sigma\right)_{\partial K}
  \right\}
  &=&
  (S^n,\sigma)_\Omega +
  \triangle t \sum_K  \left(F(S^n) q^{n+1}, \sigma\right)_K.


@f}

我们引入了一个DiscreteTime类型的对象，以便在代码中保持对时间和时间步长的当前值的跟踪。这个类封装了许多关于调整时间步长和在指定的最终时间停止的复杂情况。




<a name="Spacediscretization"></a><h3>Space discretization</h3>


在每个时间步长中，我们再对速度和压力应用 @ref step_20 "step-20 "的混合有限方法。为了得到良好的解决，我们对 $\mathbf{u}$ 选择Raviart-Thomas空间 $RT_{k}$ ，对 $p$ 选择 $DGQ_{k}$ 类的不连续元素。对于饱和度，我们也将选择 $DGQ_{k}$ 空间。

由于我们有不连续的空间，我们必须考虑如何评估细胞之间界面上的项，因为不连续的函数在那里没有真正定义。特别是，我们必须给饱和度方程左边的最后一个项赋予一个意义。为此，让我们定义，我们要在以下意义上评估它。

@f{eqnarray*}
  &&\left(F(S^n) (\mathbf n \cdot \mathbf{u}^{n+1}), \sigma\right)_{\partial K}
  \\
  &&\qquad =
  \left(F(S^n_+) (\mathbf n \cdot \mathbf{u}^{n+1}_+), \sigma\right)_{\partial K_+}
  +
  \left(F(S^n_-) (\mathbf n \cdot \mathbf{u}^{n+1}_-), \sigma\right)_{\partial K_-},


@f}

其中 $\partial K_{-} \dealcoloneq \{x\in \partial K, \mathbf{u}(x) \cdot \mathbf{n}<0\}$ 表示流入边界， $\partial K_{+} \dealcoloneq \{\partial K \setminus
\partial K_{-}\}$ 是边界的流出部分。数量 $S_+,\mathbf{u}_+$ 对应于当前单元上的这些变量值，而 $S_-,\mathbf{u}_-$ （需要在 $K$ 边界的流入部分）是取自邻近单元的数量。关于非连续元素技术和通量评估的更多背景，也可以在步骤12和步骤12b中找到。




<a name="Linearsolvers"></a><h3>Linear solvers</h3>


这个程序中使用的线性求解器是对步骤20中使用的线性求解器的直接扩展（但没有LinearOperator）。从本质上讲，我们只需将一切从两个解元扩展到三个解元。如果我们使用上面提到的离散空间，并将形状函数放入双线性形式中，我们得出以下线性系统，以解决时间步长 $n+1$  。

@f[
\left(
\begin{array}{ccc}
M^u(S^{n}) & B^{T}& 0\\
B &    0 & 0\\
\triangle t\; H &    0& M^S
\end{array}
\right)
\left(
\begin{array}{c}
\mathbf{U}^{n+1} \\ P^{n+1} \\ S^{n+1}
\end{array}
\right)
=
\left(
\begin{array}{c}
0 \\ F_2 \\ F_3
\end{array}
\right)


@f]

其中各个矩阵和向量的定义如下：使用形状函数 $\mathbf v_i$ （类型为Raviart Thomas $RT_k$ ）定义速度，使用 $\phi_i$ （类型为 $DGQ_k$  ）定义压力和饱和度。

@f{eqnarray*}
M^u(S^n)_{ij} &=&
\left((\mathbf{K}\lambda(S^n))^{-1} \mathbf{v}_i,\mathbf
v_j\right)_\Omega,
\\
B_{ij} &=&


-(\nabla \cdot \mathbf v_j, \phi_i)_\Omega,
\\
H_{ij} &=&


  -
  \sum_K
  \left\{
  \left(F(S^n) \mathbf v_i, \nabla \phi_j)\right)_K


  -
  \left(F(S^n_+) (\mathbf n \cdot (\mathbf v_i)_+), \phi_j\right)_{\partial K_+}


  -
  \left(F(S^n_-) (\mathbf n \cdot (\mathbf v_i)_-), \phi_j\right)_{\partial K_-},
  \right\}
\\
M^S_{ij} &=&
(\phi_i, \phi_j)_\Omega,
\\
(F_2)_i &=&


-(q^{n+1},\phi_i)_\Omega,
\\
(F_3)_i &=&
(S^n,\phi_i)_\Omega +\triangle t \sum_K  \left(F(S^n) q^{n+1}, \phi_i\right)_K.


@f}



 @note  由于历史原因，与第20步相比，矩阵 $B$ 和 $B^T$ 的作用在本程序中被还原了。换句话说，这里 $B$ 指的是发散， $B^T$ 指的是梯度算子，而在第20步中则是相反。

上面的系统出现了一个复杂的问题。由于矩阵 $H_{ij}$ 隐含地依赖于 $\mathbf u^{n+1}$ （需要速度来确定细胞边界 $\partial K$ 的哪些部分是流入或流出的部分），我们只能在解决了速度问题之后才能组装这个矩阵。

然后，求解方案包括以下步骤。<ol>  <li>  使用步骤20中介绍的Schur补足技术求解压力 $p^{n+1}$ 。

    <li>  求解速度 $\mathbf u^{n+1}$ ，也是在步骤20中讨论的。

    <li>  计算项 $F_3-\triangle t\; H \mathbf u^{n+1}$  ，使用刚刚计算的速度。

    <li>  求解饱和度  $S^{n+1}$  。   </ol> 

在这个方案中，我们实际上从未建立过矩阵 $H$ ，而是在我们准备好后生成第三个方程的右手边。

在程序中，我们使用一个变量 <code>solution</code> 来存储当前时间步骤的解决方案。在每一步结束时，我们将其内容，即其所有的三个块状成分，复制到变量 <code>old_solution</code> 中，以便在下一个时间步骤中使用。




<a name="Choosingatimestep"></a><h3>Choosing a time step</h3>


在双曲输运方程中，像我们要解决的饱和方程的一般经验法则是，如果我们使用显式时间步长方案，那么我们应该使用一个时间步长，使粒子在一个时间步长内所能走的距离不大于一个细胞的直径。换句话说，在这里，我们应该选择

@f[
  \triangle t_{n+1} \le \frac h{|\mathbf{u}^{n+1}(\mathbf{x})|}.


@f]

幸运的是，我们处在一个可以做到这一点的位置：我们只需要当我们想集合饱和方程的右边时的时间步长，也就是在我们已经解出 $\mathbf{u}^{n+1}$ 之后。因此，在求解速度之后，我们要做的就是在域中的所有正交点上循环，确定速度的最大幅度。然后我们可以将饱和方程的时间步长设定为

@f[
  \triangle t_{n+1} = \frac {\min_K h_K}{\max_{\mathbf{x}}|\mathbf{u}^{n+1}(\mathbf{x})|}.


@f]



为什么要这样做呢？如果我们不这样做，那么我们就会发现很多地方的饱和度大于1或小于0，这一点很容易得到验证。请记住，饱和度对应于流体混合物中的水比例，因此在物理上必须在0和1之间）。另一方面，如果我们根据上面列出的标准选择时间步长，这种情况只会非常非常少地发生，事实上在整个程序运行中只有一次。然而，为了安全起见，我们在每个时间步长结束时运行一个函数 <code>project_back_saturation</code> ，如果饱和度已经超出了物理范围，则简单地将其投射回区间 $[0,1]$ 。这很有用，因为函数 $\lambda(S)$ 和 $F(S)$ 并不代表这个范围之外的任何物理现象，而且一旦我们有负的饱和度或大于1的饱和度，我们不应该期望程序做任何有用的事情。

请注意，我们在第23步和第24步中也会对时间步长有类似的限制，在这两步中我们要解决与时间有关的波浪方程，即另一个双曲问题。我们还将在下面的<a href="#extensions">possible
extensions to this program</a>一节中再来讨论时间步长的选择问题。




<a name="Thetestcase"></a><h3>The test case</h3>


为了简单起见，本程序假定没有源头，  $q=0$  ，并且异质多孔介质是各向同性的  $\mathbf{K}(\mathbf{x}) =
k(\mathbf{x}) \mathbf{I}$  。其中第一个假设在油藏中是一个现实的假设：除了注水井和生产井之外，通常没有液体突然出现或消失的机制。第二个假设更难证明：在微观层面上，大多数岩石是各向同性的，因为它们是由相互连接的孔隙网络组成的。然而，这种微观尺度超出了今天计算机模拟的范围，我们不得不满足于模拟米级的东西。然而，在这个尺度上，流体运输通常是通过岩石中的裂缝网络，而不是通过孔隙发生的。然而，裂缝通常是由岩层中的外部应力场造成的（例如由构造断层造成的），因此裂缝是大致排列的。这就导致了这样一种情况：在平行于裂缝的方向上，渗透率往往比垂直于裂缝的方向大几个数量级。然而，在储层模拟中通常面临的一个问题是，建模者不知道裂缝的方向，因为油藏不容易被检查到。在这种情况下，唯一的解决办法是假设有效的、各向同性的渗透率。

无论怎样，这两个限制，即无源和各向同性，只要在程序中写上几行代码就能轻松解除。

接下来，为了简单起见，我们的数值模拟将在 $\Omega = [0,1]\times [0,1]$ 的单元格上进行，即 $t\in [0,T]$ 。我们的初始条件是 $S(\mathbf{x},0)=0$ ；在油藏图片中， $S$ 将表示水的饱和度，这意味着油藏一开始就含有纯油。请注意，我们不需要任何压力或速度的初始条件，因为这些方程不包含这些变量的时间导数。最后，我们施加以下压力边界条件。

@f[
  p(\mathbf{x},t)=1-x_1 \qquad \textrm{on}\ \partial\Omega.


@f]

由于压力和速度求解的是混合形式的泊松方程，所以施加的压力导致了速度的流场。另一方面，这个流场决定了边界的某一部分是流入还是流出，这很重要，因为我们必须在边界的流入部分施加饱和度的边界条件。

@f[
  \Gamma_{in}(t) = \{\mathbf{x}\in\partial\Omega:
                     \mathbf{n} \cdot \mathbf{u}(\mathbf{x},t) < 0\}.


@f]

在这个流入的边界上，我们施加以下的饱和值。

@f{eqnarray}
  S(\mathbf{x},t) = 1 & \textrm{on}\ \Gamma_{in}\cap\{x_1=0\},
  \\
  S(\mathbf{x},t) = 0 & \textrm{on}\ \Gamma_{in}\backslash \{x_1=0\}.


@f}

换句话说，我们有纯水在左边进入储层，而边界的其他部分与储层的未受干扰部分接触，只要这些边界上发生流入，纯油就会进入。

在我们的模拟中，我们选择总流动性为

@f[
  \lambda (S) = \frac{1.0}{\mu} S^2 +(1-S)^2


@f]

其中我们用 $\mu=0.2$ 表示粘度。此外，水的分流量由以下公式给出

@f[
  F(S)=\frac{S^2}{S^2+\mu (1-S)^2}


@f]



 @note  几年后在step-43中再来看这个测试案例，发现这个测试案例的设置有一个奇怪之处。为此，考虑我们可以将饱和度的平流方程改写为  $S_{t} + (\mathbf{u}
F'(S)) \cdot \nabla S = 0$  。现在，在初始时间，我们有 $S=0$ ，在给定的函数 $F(S)$ 的选择下，我们正好有 $F'(0)=0$ 。换句话说，在 $t=0$ 处，方程对所有 $\mathbf x$ 都还原为 $S_t=0$ ，所以饱和度在任何地方都是零，而且在任何地方都会保持零！这就是为什么在 $\mathbf x$ 处的饱和度为零。尽管 $\mathbf u$ 不一定是零：组合流体在移动，但我们选择的部分通量 $F(S)$ 是这样的：无穷小量的润湿流体也只以无穷小的速度移动（也就是说，它们粘附在介质上的程度比它们所嵌入的非润湿相要大）。也就是说，我们如何将这一点与润湿性液体从左边侵入，导致<a href="#Results">results section</a>中看到的流动模式的知识联系起来？这就是我们进入数学的地方。像我们在这里考虑的传输方程有无限多的解决方案，但其中只有一个是物理的：由所谓的粘性极限产生的解决方案，称为<a
href="http://en.wikipedia.org/wiki/Viscosity_solution">viscosity
solution</a>。事情是这样的，用不连续的元素，我们到达了这个粘性极限，因为使用数值通量在数值方案中引入了有限量的人工粘性。另一方面，在step-43中，我们在每个单元上使用与 $\|\mathbf u F'(S)\|$ 成比例的人工粘度，在初始时间是零。因此，那里的饱和度为零，并保持为零；然后我们得到的解是<i>one</i>的平流方程解，但如果不进一步改变，该方法不会收敛到粘性解。因此，我们将在该程序中使用一个不同的初始条件。


最后，回到测试案例的描述，我们将展示用 @ref step_20  "step-20 "的结果部分末尾介绍的两个渗透率函数计算的结果。   <ul>   <li>  一个函数，模拟一个蜿蜒穿过领域的单一裂缝。与step-20相类似，但考虑到我们这里的几何形状略有不同，我们用以下函数来描述。   @f[
    k(\mathbf x)
    =
    \max \left\{ e^{-\left(\frac{x_2-\frac 12 - 0.1\sin(10x_1)}{0.1}\right)^2}, 0.01 \right\}.
  @f]

  取最大值是必要的，以确保最大和最小渗透率之间的比率保持有界。如果我们不这样做，渗透率将跨越许多数量级。另一方面，最大和最小渗透率之间的比率是舒尔补矩阵的条件数的一个因素，如果太大，会导致我们的线性求解器不再正常收敛的问题。

    <li>  一个模拟某种随机介质的函数。在这里，我们选择@f{eqnarray*}
    k(\mathbf x)
    &=&
    \min \left\{ \max \left\{ \sum_{i=1}^N \sigma_i(\mathbf{x}), 0.01 \right\}, 4\right\},
    \\
    \sigma_i(\mathbf x)
    &=&
    e^{-\left(\frac{|\mathbf{x}-\mathbf{x}_i|}{0.05}\right)^2},
  @f}。

  其中中心 $\mathbf{x}_i$ 是域内 $N$ 随机选择的位置。这个函数模拟了一个领域，其中有 $N$ 个渗透率较高的中心（例如，岩石已经开裂）嵌入到一个更原始的、未受干扰的背景岩石矩阵中。请注意，在这里我们切断了上方和下方的渗透率函数，以确保有界的条件数。   </ul> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * This program is an adaptation of step-20 and includes some technique of DG
 * methods from step-12. A good part of the program is therefore very similar
 * to step-20 and we will not comment again on these parts. Only the new stuff
 * will be discussed in more detail.
 * 

 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * All of these include files have been used before:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_raviart_thomas.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * @endcode
 * 
 * In this program, we use a tensor-valued coefficient. Since it may have a
 * spatial dependence, we consider it a tensor-valued function. The following
 * include file provides the <code>TensorFunction</code> class that offers
 * such functionality:
 * 
 * @code
 * #include <deal.II/base/tensor_function.h>
 * 
 * @endcode
 * 
 * Additionally, we use the class <code>DiscreteTime</code> to perform
 * operations related to time incrementation.
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step21
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeTwoPhaseFlowProblemcodeclass"></a> 
 * <h3>The <code>TwoPhaseFlowProblem</code> class</h3>
 * 

 * 
 * This is the main class of the program. It is close to the one of step-20,
 * but with a few additional functions:
 *   

 * 
 * <ul> <li><code>assemble_rhs_S</code> assembles the right hand side of the
 * saturation equation. As explained in the introduction, this can't be
 * integrated into <code>assemble_rhs</code> since it depends on the
 * velocity that is computed in the first part of the time step.
 *   

 * 
 * <li><code>get_maximal_velocity</code> does as its name suggests. This
 * function is used in the computation of the time step size.
 *   

 * 
 * <li><code>project_back_saturation</code> resets all saturation degrees
 * of freedom with values less than zero to zero, and all those with
 * saturations greater than one to one.  </ul>
 *   

 * 
 * The rest of the class should be pretty much obvious. The
 * <code>viscosity</code> variable stores the viscosity $\mu$ that enters
 * several of the formulas in the nonlinear equations. The variable
 * <code>time</code> keeps track of the time information within the
 * simulation.
 * 
 * @code
 *   template <int dim>
 *   class TwoPhaseFlowProblem
 *   {
 *   public:
 *     TwoPhaseFlowProblem(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void   make_grid_and_dofs();
 *     void   assemble_system();
 *     void   assemble_rhs_S();
 *     double get_maximal_velocity() const;
 *     void   solve();
 *     void   project_back_saturation();
 *     void   output_results() const;
 * 
 *     const unsigned int degree;
 * 
 *     Triangulation<dim> triangulation;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     BlockSparsityPattern      sparsity_pattern;
 *     BlockSparseMatrix<double> system_matrix;
 * 
 *     const unsigned int n_refinement_steps;
 * 
 *     DiscreteTime time;
 *     double       viscosity;
 * 
 *     BlockVector<double> solution;
 *     BlockVector<double> old_solution;
 *     BlockVector<double> system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 
 * <a name="Pressurerighthandside"></a> 
 * <h4>Pressure right hand side</h4>
 * 

 * 
 * At present, the right hand side of the pressure equation is simply the
 * zero function. However, the rest of the program is fully equipped to deal
 * with anything else, if this is desired:
 * 
 * @code
 *   template <int dim>
 *   class PressureRightHandSide : public Function<dim>
 *   {
 *   public:
 *     PressureRightHandSide()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> & /*p*/,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Pressureboundaryvalues"></a> 
 * <h4>Pressure boundary values</h4>
 * 

 * 
 * The next are pressure boundary values. As mentioned in the introduction,
 * we choose a linear pressure field:
 * 
 * @code
 *   template <int dim>
 *   class PressureBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     PressureBoundaryValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       return 1 - p[0];
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Saturationboundaryvalues"></a> 
 * <h4>Saturation boundary values</h4>
 * 

 * 
 * Then we also need boundary values on the inflow portions of the
 * boundary. The question whether something is an inflow part is decided
 * when assembling the right hand side, we only have to provide a functional
 * description of the boundary values. This is as explained in the
 * introduction:
 * 
 * @code
 *   template <int dim>
 *   class SaturationBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     SaturationBoundaryValues()
 *       : Function<dim>(1)
 *     {}
 * 
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       if (p[0] == 0)
 *         return 1;
 *       else
 *         return 0;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Initialdata"></a> 
 * <h4>Initial data</h4>
 * 

 * 
 * Finally, we need initial data. In reality, we only need initial data for
 * the saturation, but we are lazy, so we will later, before the first time
 * step, simply interpolate the entire solution for the previous time step
 * from a function that contains all vector components.
 *   

 * 
 * We therefore simply create a function that returns zero in all
 * components. We do that by simply forward every function to the
 * Functions::ZeroFunction class. Why not use that right away in the places of
 * this program where we presently use the <code>InitialValues</code> class?
 * Because this way it is simpler to later go back and choose a different
 * function for initial values.
 * 
 * @code
 *   template <int dim>
 *   class InitialValues : public Function<dim>
 *   {
 *   public:
 *     InitialValues()
 *       : Function<dim>(dim + 2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override
 *     {
 *       return Functions::ZeroFunction<dim>(dim + 2).value(p, component);
 *     }
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override
 *     {
 *       Functions::ZeroFunction<dim>(dim + 2).vector_value(p, values);
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversepermeabilitytensor"></a> 
 * <h3>The inverse permeability tensor</h3>
 * 

 * 
 * As announced in the introduction, we implement two different permeability
 * tensor fields. Each of them we put into a namespace of its own, so that
 * it will be easy later to replace use of one by the other in the code.
 * 

 * 
 * 
 * <a name="Singlecurvingcrackpermeability"></a> 
 * <h4>Single curving crack permeability</h4>
 * 

 * 
 * The first function for the permeability was the one that models a single
 * curving crack. It was already used at the end of step-20, and its
 * functional form is given in the introduction of the present tutorial
 * program. As in some previous programs, we have to declare a (seemingly
 * unnecessary) default constructor of the KInverse class to avoid warnings
 * from some compilers:
 * 
 * @code
 *   namespace SingleCurvingCrack
 *   {
 *     template <int dim>
 *     class KInverse : public TensorFunction<2, dim>
 *     {
 *     public:
 *       KInverse()
 *         : TensorFunction<2, dim>()
 *       {}
 * 
 *       virtual void
 *       value_list(const std::vector<Point<dim>> &points,
 *                  std::vector<Tensor<2, dim>> &  values) const override
 *       {
 *         Assert(points.size() == values.size(),
 *                ExcDimensionMismatch(points.size(), values.size()));
 * 
 *         for (unsigned int p = 0; p < points.size(); ++p)
 *           {
 *             values[p].clear();
 * 
 *             const double distance_to_flowline =
 *               std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0]));
 * 
 *             const double permeability =
 *               std::max(std::exp(-(distance_to_flowline * distance_to_flowline) /
 *                                 (0.1 * 0.1)),
 *                        0.01);
 * 
 *             for (unsigned int d = 0; d < dim; ++d)
 *               values[p][d][d] = 1. / permeability;
 *           }
 *       }
 *     };
 *   } // namespace SingleCurvingCrack
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Randommediumpermeability"></a> 
 * <h4>Random medium permeability</h4>
 * 

 * 
 * This function does as announced in the introduction, i.e. it creates an
 * overlay of exponentials at random places. There is one thing worth
 * considering for this class. The issue centers around the problem that the
 * class creates the centers of the exponentials using a random function. If
 * we therefore created the centers each time we create an object of the
 * present type, we would get a different list of centers each time. That's
 * not what we expect from classes of this type: they should reliably
 * represent the same function.
 *   

 * 
 * The solution to this problem is to make the list of centers a static
 * member variable of this class, i.e. there exists exactly one such
 * variable for the entire program, rather than for each object of this
 * type. That's exactly what we are going to do.
 *   

 * 
 * The next problem, however, is that we need a way to initialize this
 * variable. Since this variable is initialized at the beginning of the
 * program, we can't use a regular member function for that since there may
 * not be an object of this type around at the time. The C++ standard
 * therefore says that only non-member and static member functions can be
 * used to initialize a static variable. We use the latter possibility by
 * defining a function <code>get_centers</code> that computes the list of
 * center points when called.
 *   

 * 
 * Note that this class works just fine in both 2d and 3d, with the only
 * difference being that we use more points in 3d: by experimenting we find
 * that we need more exponentials in 3d than in 2d (we have more ground to
 * cover, after all, if we want to keep the distance between centers roughly
 * equal), so we choose 40 in 2d and 100 in 3d. For any other dimension, the
 * function does presently not know what to do so simply throws an exception
 * indicating exactly this.
 * 
 * @code
 *   namespace RandomMedium
 *   {
 *     template <int dim>
 *     class KInverse : public TensorFunction<2, dim>
 *     {
 *     public:
 *       KInverse()
 *         : TensorFunction<2, dim>()
 *       {}
 * 
 *       virtual void
 *       value_list(const std::vector<Point<dim>> &points,
 *                  std::vector<Tensor<2, dim>> &  values) const override
 *       {
 *         Assert(points.size() == values.size(),
 *                ExcDimensionMismatch(points.size(), values.size()));
 * 
 *         for (unsigned int p = 0; p < points.size(); ++p)
 *           {
 *             values[p].clear();
 * 
 *             double permeability = 0;
 *             for (unsigned int i = 0; i < centers.size(); ++i)
 *               permeability += std::exp(-(points[p] - centers[i]).norm_square() /
 *                                        (0.05 * 0.05));
 * 
 *             const double normalized_permeability =
 *               std::min(std::max(permeability, 0.01), 4.);
 * 
 *             for (unsigned int d = 0; d < dim; ++d)
 *               values[p][d][d] = 1. / normalized_permeability;
 *           }
 *       }
 * 
 *     private:
 *       static std::vector<Point<dim>> centers;
 * 
 *       static std::vector<Point<dim>> get_centers()
 *       {
 *         const unsigned int N =
 *           (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented()));
 * 
 *         std::vector<Point<dim>> centers_list(N);
 *         for (unsigned int i = 0; i < N; ++i)
 *           for (unsigned int d = 0; d < dim; ++d)
 *             centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX;
 * 
 *         return centers_list;
 *       }
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     std::vector<Point<dim>>
 *       KInverse<dim>::centers = KInverse<dim>::get_centers();
 *   } // namespace RandomMedium
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversemobilityandsaturationfunctions"></a> 
 * <h3>The inverse mobility and saturation functions</h3>
 * 

 * 
 * There are two more pieces of data that we need to describe, namely the
 * inverse mobility function and the saturation curve. Their form is also
 * given in the introduction:
 * 
 * @code
 *   double mobility_inverse(const double S, const double viscosity)
 *   {
 *     return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S));
 *   }
 * 
 *   double fractional_flow(const double S, const double viscosity)
 *   {
 *     return S * S / (S * S + viscosity * (1 - S) * (1 - S));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * The linear solvers we use are also completely analogous to the ones used
 * in step-20. The following classes are therefore copied verbatim from
 * there. Note that the classes here are not only copied from
 * step-20, but also duplicate classes in deal.II. In a future version of this
 * example, they should be replaced by an efficient method, though. There is a
 * single change: if the size of a linear system is small, i.e. when the mesh
 * is very coarse, then it is sometimes not sufficient to set a maximum of
 * <code>src.size()</code> CG iterations before the solver in the
 * <code>vmult()</code> function converges. (This is, of course, a result of
 * numerical round-off, since we know that on paper, the CG method converges
 * in at most <code>src.size()</code> steps.) As a consequence, we set the
 * maximum number of iterations equal to the maximum of the size of the linear
 * system and 200.
 * 
 * @code
 *   template <class MatrixType>
 *   class InverseMatrix : public Subscriptor
 *   {
 *   public:
 *     InverseMatrix(const MatrixType &m)
 *       : matrix(&m)
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       SolverControl solver_control(std::max<unsigned int>(src.size(), 200),
 *                                    1e-8 * src.l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       dst = 0;
 * 
 *       cg.solve(*matrix, dst, src, PreconditionIdentity());
 *     }
 * 
 *   private:
 *     const SmartPointer<const MatrixType> matrix;
 *   };
 * 
 * 
 * 
 *   class SchurComplement : public Subscriptor
 *   {
 *   public:
 *     SchurComplement(const BlockSparseMatrix<double> &          A,
 *                     const InverseMatrix<SparseMatrix<double>> &Minv)
 *       : system_matrix(&A)
 *       , m_inverse(&Minv)
 *       , tmp1(A.block(0, 0).m())
 *       , tmp2(A.block(0, 0).m())
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       system_matrix->block(0, 1).vmult(tmp1, src);
 *       m_inverse->vmult(tmp2, tmp1);
 *       system_matrix->block(1, 0).vmult(dst, tmp2);
 *     }
 * 
 *   private:
 *     const SmartPointer<const BlockSparseMatrix<double>>           system_matrix;
 *     const SmartPointer<const InverseMatrix<SparseMatrix<double>>> m_inverse;
 * 
 *     mutable Vector<double> tmp1, tmp2;
 *   };
 * 
 * 
 * 
 *   class ApproximateSchurComplement : public Subscriptor
 *   {
 *   public:
 *     ApproximateSchurComplement(const BlockSparseMatrix<double> &A)
 *       : system_matrix(&A)
 *       , tmp1(A.block(0, 0).m())
 *       , tmp2(A.block(0, 0).m())
 *     {}
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const
 *     {
 *       system_matrix->block(0, 1).vmult(tmp1, src);
 *       system_matrix->block(0, 0).precondition_Jacobi(tmp2, tmp1);
 *       system_matrix->block(1, 0).vmult(dst, tmp2);
 *     }
 * 
 *   private:
 *     const SmartPointer<const BlockSparseMatrix<double>> system_matrix;
 * 
 *     mutable Vector<double> tmp1, tmp2;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="codeTwoPhaseFlowProblemcodeclassimplementation"></a> 
 * <h3><code>TwoPhaseFlowProblem</code> class implementation</h3>
 * 

 * 
 * Here now the implementation of the main class. Much of it is actually
 * copied from step-20, so we won't comment on it in much detail. You should
 * try to get familiar with that program first, then most of what is
 * happening here should be mostly clear.
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemTwoPhaseFlowProblem"></a> 
 * <h4>TwoPhaseFlowProblem::TwoPhaseFlowProblem</h4>
 * 

 * 
 * First for the constructor. We use $RT_k \times DQ_k \times DQ_k$
 * spaces. For initializing the DiscreteTime object, we don't set the time
 * step size in the constructor because we don't have its value yet.
 * The time step size is initially set to zero, but it will be computed
 * before it is needed to increment time, as described in a subsection of
 * the introduction. The time object internally prevents itself from being
 * incremented when $dt = 0$, forcing us to set a non-zero desired size for
 * $dt$ before advancing time.
 * 
 * @code
 *   template <int dim>
 *   TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree)
 *     : degree(degree)
 *     , fe(FE_RaviartThomas<dim>(degree),
 *          1,
 *          FE_DGQ<dim>(degree),
 *          1,
 *          FE_DGQ<dim>(degree),
 *          1)
 *     , dof_handler(triangulation)
 *     , n_refinement_steps(5)
 *     , time(/*start time*/ 0., /*end time*/ 1.)
 *     , viscosity(0.2)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemmake_grid_and_dofs"></a> 
 * <h4>TwoPhaseFlowProblem::make_grid_and_dofs</h4>
 * 

 * 
 * This next function starts out with well-known functions calls that create
 * and refine a mesh, and then associate degrees of freedom with it. It does
 * all the same things as in step-20, just now for three components instead
 * of two.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0, 1);
 *     triangulation.refine_global(n_refinement_steps);
 * 
 *     dof_handler.distribute_dofs(fe);
 *     DoFRenumbering::component_wise(dof_handler);
 * 
 *     const std::vector<types::global_dof_index> dofs_per_component =
 *       DoFTools::count_dofs_per_fe_component(dof_handler);
 *     const unsigned int n_u = dofs_per_component[0],
 *                        n_p = dofs_per_component[dim],
 *                        n_s = dofs_per_component[dim + 1];
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (" << n_u << '+' << n_p << '+' << n_s << ')' << std::endl
 *               << std::endl;
 * 
 *     const unsigned int n_couplings = dof_handler.max_couplings_between_dofs();
 * 
 *     sparsity_pattern.reinit(3, 3);
 *     sparsity_pattern.block(0, 0).reinit(n_u, n_u, n_couplings);
 *     sparsity_pattern.block(1, 0).reinit(n_p, n_u, n_couplings);
 *     sparsity_pattern.block(2, 0).reinit(n_s, n_u, n_couplings);
 *     sparsity_pattern.block(0, 1).reinit(n_u, n_p, n_couplings);
 *     sparsity_pattern.block(1, 1).reinit(n_p, n_p, n_couplings);
 *     sparsity_pattern.block(2, 1).reinit(n_s, n_p, n_couplings);
 *     sparsity_pattern.block(0, 2).reinit(n_u, n_s, n_couplings);
 *     sparsity_pattern.block(1, 2).reinit(n_p, n_s, n_couplings);
 *     sparsity_pattern.block(2, 2).reinit(n_s, n_s, n_couplings);
 * 
 *     sparsity_pattern.collect_sizes();
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
 *     sparsity_pattern.compress();
 * 
 * 
 *     system_matrix.reinit(sparsity_pattern);
 * 
 * 
 *     solution.reinit(3);
 *     solution.block(0).reinit(n_u);
 *     solution.block(1).reinit(n_p);
 *     solution.block(2).reinit(n_s);
 *     solution.collect_sizes();
 * 
 *     old_solution.reinit(3);
 *     old_solution.block(0).reinit(n_u);
 *     old_solution.block(1).reinit(n_p);
 *     old_solution.block(2).reinit(n_s);
 *     old_solution.collect_sizes();
 * 
 *     system_rhs.reinit(3);
 *     system_rhs.block(0).reinit(n_u);
 *     system_rhs.block(1).reinit(n_p);
 *     system_rhs.block(2).reinit(n_s);
 *     system_rhs.collect_sizes();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemassemble_system"></a> 
 * <h4>TwoPhaseFlowProblem::assemble_system</h4>
 * 

 * 
 * This is the function that assembles the linear system, or at least
 * everything except the (1,3) block that depends on the still-unknown
 * velocity computed during this time step (we deal with this in
 * <code>assemble_rhs_S</code>). Much of it is again as in step-20, but we
 * have to deal with some nonlinearity this time.  However, the top of the
 * function is pretty much as usual (note that we set matrix and right hand
 * side to zero at the beginning &mdash; something we didn't have to do for
 * stationary problems since there we use each matrix object only once and
 * it is empty at the beginning anyway).
 *   

 * 
 * Note that in its present form, the function uses the permeability
 * implemented in the RandomMedium::KInverse class. Switching to the single
 * curved crack permeability function is as simple as just changing the
 * namespace name.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     QGauss<dim>     quadrature_formula(degree + 2);
 *     QGauss<dim - 1> face_quadrature_formula(degree + 2);
 * 
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const PressureRightHandSide<dim>  pressure_right_hand_side;
 *     const PressureBoundaryValues<dim> pressure_boundary_values;
 *     const RandomMedium::KInverse<dim> k_inverse;
 * 
 *     std::vector<double>         pressure_rhs_values(n_q_points);
 *     std::vector<double>         boundary_values(n_face_q_points);
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);
 * 
 *     std::vector<Vector<double>>              old_solution_values(n_q_points,
 *                                                                  Vector<double>(dim + 2));
 *     std::vector<std::vector<Tensor<1, dim>>> old_solution_grads(
 *       n_q_points, std::vector<Tensor<1, dim>>(dim + 2));
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 *     const FEValuesExtractors::Scalar saturation(dim + 1);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 * @endcode
 * 
 * Here's the first significant difference: We have to get the values
 * of the saturation function of the previous time step at the
 * quadrature points. To this end, we can use the
 * FEValues::get_function_values (previously already used in step-9,
 * step-14 and step-15), a function that takes a solution vector and
 * returns a list of function values at the quadrature points of the
 * present cell. In fact, it returns the complete vector-valued
 * solution at each quadrature point, i.e. not only the saturation but
 * also the velocities and pressure:
 * 
 * @code
 *         fe_values.get_function_values(old_solution, old_solution_values);
 * 
 * @endcode
 * 
 * Then we also have to get the values of the pressure right hand side
 * and of the inverse permeability tensor at the quadrature points:
 * 
 * @code
 *         pressure_right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                             pressure_rhs_values);
 *         k_inverse.value_list(fe_values.get_quadrature_points(),
 *                              k_inverse_values);
 * 
 * @endcode
 * 
 * With all this, we can now loop over all the quadrature points and
 * shape functions on this cell and assemble those parts of the matrix
 * and right hand side that we deal with in this function. The
 * individual terms in the contributions should be self-explanatory
 * given the explicit form of the bilinear form stated in the
 * introduction:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const double old_s = old_solution_values[q](dim + 1);
 * 
 *               const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q);
 *               const double div_phi_i_u = fe_values[velocities].divergence(i, q);
 *               const double phi_i_p     = fe_values[pressure].value(i, q);
 *               const double phi_i_s     = fe_values[saturation].value(i, q);
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const Tensor<1, dim> phi_j_u =
 *                     fe_values[velocities].value(j, q);
 *                   const double div_phi_j_u =
 *                     fe_values[velocities].divergence(j, q);
 *                   const double phi_j_p = fe_values[pressure].value(j, q);
 *                   const double phi_j_s = fe_values[saturation].value(j, q);
 * 
 *                   local_matrix(i, j) +=
 *                     (phi_i_u * k_inverse_values[q] *
 *                        mobility_inverse(old_s, viscosity) * phi_j_u -
 *                      div_phi_i_u * phi_j_p - phi_i_p * div_phi_j_u +
 *                      phi_i_s * phi_j_s) *
 *                     fe_values.JxW(q);
 *                 }
 * 
 *               local_rhs(i) +=
 *                 (-phi_i_p * pressure_rhs_values[q]) * fe_values.JxW(q);
 *             }
 * 
 * 
 * @endcode
 * 
 * Next, we also have to deal with the pressure boundary values. This,
 * again is as in step-20:
 * 
 * @code
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary())
 *             {
 *               fe_face_values.reinit(cell, face);
 * 
 *               pressure_boundary_values.value_list(
 *                 fe_face_values.get_quadrature_points(), boundary_values);
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   {
 *                     const Tensor<1, dim> phi_i_u =
 *                       fe_face_values[velocities].value(i, q);
 * 
 *                     local_rhs(i) +=
 *                       -(phi_i_u * fe_face_values.normal_vector(q) *
 *                         boundary_values[q] * fe_face_values.JxW(q));
 *                   }
 *             }
 * 
 * @endcode
 * 
 * The final step in the loop over all cells is to transfer local
 * contributions into the global matrix and right hand side vector:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               local_matrix(i, j));
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           system_rhs(local_dof_indices[i]) += local_rhs(i);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * So much for assembly of matrix and right hand side. Note that we do not
 * have to interpolate and apply boundary values since they have all been
 * taken care of in the weak form already.
 * 

 * 
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemassemble_rhs_S"></a> 
 * <h4>TwoPhaseFlowProblem::assemble_rhs_S</h4>
 * 

 * 
 * As explained in the introduction, we can only evaluate the right hand
 * side of the saturation equation once the velocity has been computed. We
 * therefore have this separate function to this end.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::assemble_rhs_S()
 *   {
 *     QGauss<dim>       quadrature_formula(degree + 2);
 *     QGauss<dim - 1>   face_quadrature_formula(degree + 2);
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 *     FEFaceValues<dim> fe_face_values_neighbor(fe,
 *                                               face_quadrature_formula,
 *                                               update_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     Vector<double> local_rhs(dofs_per_cell);
 * 
 *     std::vector<Vector<double>> old_solution_values(n_q_points,
 *                                                     Vector<double>(dim + 2));
 *     std::vector<Vector<double>> old_solution_values_face(n_face_q_points,
 *                                                          Vector<double>(dim +
 *                                                                         2));
 *     std::vector<Vector<double>> old_solution_values_face_neighbor(
 *       n_face_q_points, Vector<double>(dim + 2));
 *     std::vector<Vector<double>> present_solution_values(n_q_points,
 *                                                         Vector<double>(dim +
 *                                                                        2));
 *     std::vector<Vector<double>> present_solution_values_face(
 *       n_face_q_points, Vector<double>(dim + 2));
 * 
 *     std::vector<double>                  neighbor_saturation(n_face_q_points);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     SaturationBoundaryValues<dim> saturation_boundary_values;
 * 
 *     const FEValuesExtractors::Scalar saturation(dim + 1);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         local_rhs = 0;
 *         fe_values.reinit(cell);
 * 
 *         fe_values.get_function_values(old_solution, old_solution_values);
 *         fe_values.get_function_values(solution, present_solution_values);
 * 
 * @endcode
 * 
 * First for the cell terms. These are, following the formulas in the
 * introduction, $(S^n,\sigma)-(F(S^n) \mathbf{v}^{n+1},\nabla
 * \sigma)$, where $\sigma$ is the saturation component of the test
 * function:
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const double   old_s = old_solution_values[q](dim + 1);
 *               Tensor<1, dim> present_u;
 *               for (unsigned int d = 0; d < dim; ++d)
 *                 present_u[d] = present_solution_values[q](d);
 * 
 *               const double         phi_i_s = fe_values[saturation].value(i, q);
 *               const Tensor<1, dim> grad_phi_i_s =
 *                 fe_values[saturation].gradient(i, q);
 * 
 *               local_rhs(i) +=
 *                 (time.get_next_step_size() * fractional_flow(old_s, viscosity) *
 *                    present_u * grad_phi_i_s +
 *                  old_s * phi_i_s) *
 *                 fe_values.JxW(q);
 *             }
 * 
 * @endcode
 * 
 * Secondly, we have to deal with the flux parts on the face
 * boundaries. This was a bit more involved because we first have to
 * determine which are the influx and outflux parts of the cell
 * boundary. If we have an influx boundary, we need to evaluate the
 * saturation on the other side of the face (or the boundary values,
 * if we are at the boundary of the domain).
 *         

 * 
 * All this is a bit tricky, but has been explained in some detail
 * already in step-9. Take a look there how this is supposed to work!
 * 
 * @code
 *         for (const auto face_no : cell->face_indices())
 *           {
 *             fe_face_values.reinit(cell, face_no);
 * 
 *             fe_face_values.get_function_values(old_solution,
 *                                                old_solution_values_face);
 *             fe_face_values.get_function_values(solution,
 *                                                present_solution_values_face);
 * 
 *             if (cell->at_boundary(face_no))
 *               saturation_boundary_values.value_list(
 *                 fe_face_values.get_quadrature_points(), neighbor_saturation);
 *             else
 *               {
 *                 const auto         neighbor = cell->neighbor(face_no);
 *                 const unsigned int neighbor_face =
 *                   cell->neighbor_of_neighbor(face_no);
 * 
 *                 fe_face_values_neighbor.reinit(neighbor, neighbor_face);
 * 
 *                 fe_face_values_neighbor.get_function_values(
 *                   old_solution, old_solution_values_face_neighbor);
 * 
 *                 for (unsigned int q = 0; q < n_face_q_points; ++q)
 *                   neighbor_saturation[q] =
 *                     old_solution_values_face_neighbor[q](dim + 1);
 *               }
 * 
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 Tensor<1, dim> present_u_face;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   present_u_face[d] = present_solution_values_face[q](d);
 * 
 *                 const double normal_flux =
 *                   present_u_face * fe_face_values.normal_vector(q);
 * 
 *                 const bool is_outflow_q_point = (normal_flux >= 0);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   local_rhs(i) -=
 *                     time.get_next_step_size() * normal_flux *
 *                     fractional_flow((is_outflow_q_point == true ?
 *                                        old_solution_values_face[q](dim + 1) :
 *                                        neighbor_saturation[q]),
 *                                     viscosity) *
 *                     fe_face_values[saturation].value(i, q) *
 *                     fe_face_values.JxW(q);
 *               }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           system_rhs(local_dof_indices[i]) += local_rhs(i);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemsolve"></a> 
 * <h4>TwoPhaseFlowProblem::solve</h4>
 * 

 * 
 * After all these preparations, we finally solve the linear system for
 * velocity and pressure in the same way as in step-20. After that, we have
 * to deal with the saturation equation (see below):
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::solve()
 *   {
 *     const InverseMatrix<SparseMatrix<double>> m_inverse(
 *       system_matrix.block(0, 0));
 *     Vector<double> tmp(solution.block(0).size());
 *     Vector<double> schur_rhs(solution.block(1).size());
 *     Vector<double> tmp2(solution.block(2).size());
 * 
 * 
 * @endcode
 * 
 * First the pressure, using the pressure Schur complement of the first
 * two equations:
 * 
 * @code
 *     {
 *       m_inverse.vmult(tmp, system_rhs.block(0));
 *       system_matrix.block(1, 0).vmult(schur_rhs, tmp);
 *       schur_rhs -= system_rhs.block(1);
 * 
 * 
 *       SchurComplement schur_complement(system_matrix, m_inverse);
 * 
 *       ApproximateSchurComplement approximate_schur_complement(system_matrix);
 * 
 *       InverseMatrix<ApproximateSchurComplement> preconditioner(
 *         approximate_schur_complement);
 * 
 * 
 *       SolverControl            solver_control(solution.block(1).size(),
 *                                    1e-12 * schur_rhs.l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       cg.solve(schur_complement, solution.block(1), schur_rhs, preconditioner);
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " CG Schur complement iterations for pressure." << std::endl;
 *     }
 * 
 * @endcode
 * 
 * Now the velocity:
 * 
 * @code
 *     {
 *       system_matrix.block(0, 1).vmult(tmp, solution.block(1));
 *       tmp *= -1;
 *       tmp += system_rhs.block(0);
 * 
 *       m_inverse.vmult(solution.block(0), tmp);
 *     }
 * 
 * @endcode
 * 
 * Finally, we have to take care of the saturation equation. The first
 * business we have here is to determine the time step using the formula
 * in the introduction. Knowing the shape of our domain and that we
 * created the mesh by regular subdivision of cells, we can compute the
 * diameter of each of our cells quite easily (in fact we use the linear
 * extensions in coordinate directions of the cells, not the
 * diameter). Note that we will learn a more general way to do this in
 * step-24, where we use the GridTools::minimal_cell_diameter function.
 *     

 * 
 * The maximal velocity we compute using a helper function to compute the
 * maximal velocity defined below, and with all this we can evaluate our
 * new time step length. We use the method
 * DiscreteTime::set_desired_next_time_step() to suggest the new
 * calculated value of the time step to the DiscreteTime object. In most
 * cases, the time object uses the exact provided value to increment time.
 * It some case, the step size may be modified further by the time object.
 * For example, if the calculated time increment overshoots the end time,
 * it is truncated accordingly.
 * 
 * @code
 *     time.set_desired_next_step_size(std::pow(0.5, double(n_refinement_steps)) /
 *                                     get_maximal_velocity());
 * 
 * @endcode
 * 
 * The next step is to assemble the right hand side, and then to pass
 * everything on for solution. At the end, we project back saturations
 * onto the physically reasonable range:
 * 
 * @code
 *     assemble_rhs_S();
 *     {
 *       SolverControl            solver_control(system_matrix.block(2, 2).m(),
 *                                    1e-8 * system_rhs.block(2).l2_norm());
 *       SolverCG<Vector<double>> cg(solver_control);
 *       cg.solve(system_matrix.block(2, 2),
 *                solution.block(2),
 *                system_rhs.block(2),
 *                PreconditionIdentity());
 * 
 *       project_back_saturation();
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " CG iterations for saturation." << std::endl;
 *     }
 * 
 * 
 *     old_solution = solution;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemoutput_results"></a> 
 * <h4>TwoPhaseFlowProblem::output_results</h4>
 * 

 * 
 * There is nothing surprising here. Since the program will do a lot of time
 * steps, we create an output file only every fifth time step and skip all
 * other time steps at the top of the file already.
 *   

 * 
 * When creating file names for output close to the bottom of the function,
 * we convert the number of the time step to a string representation that
 * is padded by leading zeros to four digits. We do this because this way
 * all output file names have the same length, and consequently sort well
 * when creating a directory listing.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::output_results() const
 *   {
 *     if (time.get_step_number() % 5 != 0)
 *       return;
 * 
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 2:
 *           solution_names = {"u", "v", "p", "S"};
 *           break;
 * 
 *         case 3:
 *           solution_names = {"u", "v", "w", "p", "S"};
 *           break;
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, solution_names);
 * 
 *     data_out.build_patches(degree + 1);
 * 
 *     std::ofstream output("solution-" +
 *                          Utilities::int_to_string(time.get_step_number(), 4) +
 *                          ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemproject_back_saturation"></a> 
 * <h4>TwoPhaseFlowProblem::project_back_saturation</h4>
 * 

 * 
 * In this function, we simply run over all saturation degrees of freedom
 * and make sure that if they should have left the physically reasonable
 * range, that they be reset to the interval $[0,1]$. To do this, we only
 * have to loop over all saturation components of the solution vector; these
 * are stored in the block 2 (block 0 are the velocities, block 1 are the
 * pressures).
 *   

 * 
 * It may be instructive to note that this function almost never triggers
 * when the time step is chosen as mentioned in the introduction. However,
 * if we choose the timestep only slightly larger, we get plenty of values
 * outside the proper range. Strictly speaking, the function is therefore
 * unnecessary if we choose the time step small enough. In a sense, the
 * function is therefore only a safety device to avoid situations where our
 * entire solution becomes unphysical because individual degrees of freedom
 * have become unphysical a few time steps earlier.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::project_back_saturation()
 *   {
 *     for (unsigned int i = 0; i < solution.block(2).size(); ++i)
 *       if (solution.block(2)(i) < 0)
 *         solution.block(2)(i) = 0;
 *       else if (solution.block(2)(i) > 1)
 *         solution.block(2)(i) = 1;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemget_maximal_velocity"></a> 
 * <h4>TwoPhaseFlowProblem::get_maximal_velocity</h4>
 * 

 * 
 * The following function is used in determining the maximal allowable time
 * step. What it does is to loop over all quadrature points in the domain
 * and find what the maximal magnitude of the velocity is.
 * 
 * @code
 *   template <int dim>
 *   double TwoPhaseFlowProblem<dim>::get_maximal_velocity() const
 *   {
 *     QGauss<dim>        quadrature_formula(degree + 2);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(fe, quadrature_formula, update_values);
 *     std::vector<Vector<double>> solution_values(n_q_points,
 *                                                 Vector<double>(dim + 2));
 *     double                      max_velocity = 0;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         fe_values.get_function_values(solution, solution_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             Tensor<1, dim> velocity;
 *             for (unsigned int i = 0; i < dim; ++i)
 *               velocity[i] = solution_values[q](i);
 * 
 *             max_velocity = std::max(max_velocity, velocity.norm());
 *           }
 *       }
 * 
 *     return max_velocity;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemrun"></a> 
 * <h4>TwoPhaseFlowProblem::run</h4>
 * 

 * 
 * This is the final function of our main class. Its brevity speaks for
 * itself. There are only two points worth noting: First, the function
 * projects the initial values onto the finite element space at the
 * beginning; the VectorTools::project function doing this requires an
 * argument indicating the hanging node constraints. We have none in this
 * program (we compute on a uniformly refined mesh), but the function
 * requires the argument anyway, of course. So we have to create a
 * constraint object. In its original state, constraint objects are
 * unsorted, and have to be sorted (using the AffineConstraints::close
 * function) before they can be used. This is what we do here, and which is
 * why we can't simply call the VectorTools::project function with an
 * anonymous temporary object <code>AffineConstraints<double>()</code> as the
 * second argument.
 *   

 * 
 * The second point worth mentioning is that we only compute the length of
 * the present time step in the middle of solving the linear system
 * corresponding to each time step. We can therefore output the present
 * time of a time step only at the end of the time step.
 * We increment time by calling the method DiscreteTime::advance_time()
 * inside the loop. Since we are reporting the time and dt after we
 * increment it, we have to call the method
 * DiscreteTime::get_previous_step_size() instead of
 * DiscreteTime::get_next_step_size(). After many steps, when the simulation
 * reaches the end time, the last dt is chosen by the DiscreteTime class in
 * such a way that the last step finishes exactly at the end time.
 * 
 * @code
 *   template <int dim>
 *   void TwoPhaseFlowProblem<dim>::run()
 *   {
 *     make_grid_and_dofs();
 * 
 *     {
 *       AffineConstraints<double> constraints;
 *       constraints.close();
 * 
 *       VectorTools::project(dof_handler,
 *                            constraints,
 *                            QGauss<dim>(degree + 2),
 *                            InitialValues<dim>(),
 *                            old_solution);
 *     }
 * 
 *     do
 *       {
 *         std::cout << "Timestep " << time.get_step_number() + 1 << std::endl;
 * 
 *         assemble_system();
 * 
 *         solve();
 * 
 *         output_results();
 * 
 *         time.advance_time();
 *         std::cout << "   Now at t=" << time.get_current_time()
 *                   << ", dt=" << time.get_previous_step_size() << '.'
 *                   << std::endl
 *                   << std::endl;
 *       }
 *     while (time.is_at_end() == false);
 *   }
 * } // namespace Step21
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * That's it. In the main function, we pass the degree of the finite element
 * space to the constructor of the TwoPhaseFlowProblem object.  Here, we use
 * zero-th degree elements, i.e. $RT_0\times DQ_0 \times DQ_0$. The rest is as
 * in all the other programs.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step21;
 * 
 *       TwoPhaseFlowProblem<2> two_phase_flow_problem(0);
 *       two_phase_flow_problem.run();
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
examples/step-21/doc/results.dox



<a name="Results"></a><h1>Results</h1>


这里介绍的代码并没有真正计算出网页上的结果。原因是，即使在一台好的电脑上，它也要运行一天以上。如果你想重现这些结果，在TwoPhaseFlowProblem的构造函数中把DiscreteTime对象的结束时间修改为`250`。

如果我们运行该程序，我们会得到以下这种输出。

@code
Number of active cells: 1024
Number of degrees of freedom: 4160 (2112+1024+1024)


Timestep 1
   22 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0326742, dt=0.0326742.


Timestep 2
   17 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0653816, dt=0.0327074.


Timestep 3
   17 CG Schur complement iterations for pressure.
   1 CG iterations for saturation.
   Now at t=0.0980651, dt=0.0326836.


...
@endcode

我们可以看到，时间步长从一开始就基本恒定，这表明域中的速度并不强烈依赖于饱和度的变化，尽管它们肯定是通过压力方程中的因子 $\lambda(S)$ 来决定的。

我们的第二个观察结果是，在第一和第二时间步之间，解决压力舒尔补足方程所需的CG迭代次数从22次下降到17次（事实上，在其余的计算中，它保持在17次左右）。原因其实很简单。在我们求解一个时间步长的压力之前，我们没有将 <code>solution</code> 变量重置为零。因此，在我们进入CG求解器时，压力（和其他变量）具有前一个时间步骤的值。由于速度和压力在计算过程中变化不大，前一个时间步骤的压力实际上是对这个时间步骤压力的一个很好的初始猜测。因此，一旦我们计算了一次压力，我们需要的迭代次数就会大大减少。

最后的观察是关于求解饱和度所需的迭代次数，也就是一次。这不应该让我们太惊讶：我们必须解决的矩阵是质量矩阵。然而，这是 $DGQ_0$ 元素的分片常数的质量矩阵，其中没有元素与相邻单元的自由度耦合。因此，该矩阵是一个对角线矩阵，很明显，我们应该能够在一次CG迭代中反转该矩阵。


说了这么多，这里有几个电影，显示了饱和度是如何随时间推移而发展的。首先，这是针对单一裂缝模型的，正如在 <code>SingleCurvingCrack::KInverse</code> 类中实现的那样。

 <img src="https://www.dealii.org/images/steps/developer/step-21.centerline.gif" alt=""> 

可以看出，富水流体主要是沿着域中间的高渗透区蜿蜒前行，而域的其他部分则大部分是不渗透的。这部电影和下一部电影是用 <code>n_refinement_steps=7</code> 生成的，导致 $128\times 128$ 的网格有大约16000个单元和大约66000个未知数。


第二部电影显示了 <code>RandomMedium::KInverse</code> 类的随机介质模型的饱和度，我们有随机分布的高渗透率中心，流体从这些区域中的一个跳到另一个。

 <img src="https://www.dealii.org/images/steps/developer/step-21.random2d.gif" alt=""> 


最后，这里是在三个空间维度上的相同情况，在一个具有 <code>n_refinement_steps=5</code> 的网格上，产生一个大约32000个单元和167000个自由度的网格。

 <img src="https://www.dealii.org/images/steps/developer/step-21.random3d.gif" alt=""> 

要重复这些计算，你所要做的就是改变行数

@code
      TwoPhaseFlowProblem<2> two_phase_flow_problem(0);
@endcode

在主函数中为

@code
      TwoPhaseFlowProblem<3> two_phase_flow_problem(0);
@endcode

可视化采用了云技术，每个单元的饱和度都由彩色但透明的云来表示。这样，人们也可以在一定程度上看到域的深处发生了什么。另一种可视化的方式是显示饱和度随时间变化的等值面。有一些技术可以透明地绘制等值面，这样就可以像洋葱的层次一样同时看到几个等值面。

那么，为什么我们不显示这样的等值面呢？问题在于等值面的计算方式：它们要求要可视化的场是连续的，所以等值面可以通过至少在单个细胞中遵循轮廓线来生成。然而，我们的饱和场是片状常数和不连续的。如果我们想为一个饱和度 $S=0.5$ 绘制一个等值面，那么在这个领域中就很可能没有一个点是真正达到饱和度的。如果我们必须在这种情况下定义等值面，我们将不得不采取细胞之间的界面，其中相邻的两个细胞之一的饱和度大于，另一个细胞的饱和度小于0.5。然而，大多数可视化程序似乎并不具备做这种转换的能力。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个项目有许多可以改进的地方。下面列出了其中的三个。事实上，所有这些问题都在构成当前程序的延续的辅导程序中得到了解决：Step-43。




<a name="Solvers"></a><h4>Solvers</h4>


目前，该程序并不是特别快：二维随机介质的计算在1000个左右的时间步长中花费了大约一天时间。相应的三维计算在800个时间步骤中几乎花了两天时间。没有比这更快的原因有两个方面。首先，我们在每个时间步骤中都要重建整个矩阵，尽管有些部分如 $B$ 、 $B^T$ 和 $M^S$ 块从未改变。

第二，我们可以在求解器和预处理器方面做得更好。目前，我们用CG方法解决Schur补数 $B^TM^u(S)^{-1}B$ ，使用 $[B^T (\textrm{diag}(M^u(S)))^{-1} B]^{-1}$ 作为预处理程序。应用这个预处理程序是很昂贵的，因为它涉及到每次解决一个线性系统。这可能适合于 @ref
step_20 的 "第20步"，在那里我们只需要解决整个问题一次。然而，在这里我们必须求解数百次，在这种情况下，值得考虑使用一种第一次设置起来比较昂贵，但以后应用起来比较便宜的预处理程序。

一种可能性是意识到我们用作预处理的矩阵， $B^T (\textrm{diag}(M^u(S)))^{-1} B$ 仍然是稀疏的，而且是对称的。如果看一下流场随时间的演变，我们还可以看到，虽然 $S$ 随时间变化很大，但压力几乎没有变化，因此 $B^T (\textrm{diag}(M^u(S)))^{-1} B \approx B^T (\textrm{diag}(M^u(S^0)))^{-1}
B$  。换句话说，第一个时间步骤的矩阵应该是一个很好的前提条件，也适用于所有后来的时间步骤。  通过一些反反复复的操作，实际上并不难得到一个SparseMatrix对象的表示。然后我们可以把它交给SparseMIC类，以形成一个稀疏的不完全Cholesky分解。形成这种分解是很昂贵的，但是我们只需要在第一个时间步骤中做一次，然后就可以在未来把它作为一个廉价的预处理程序。我们甚至可以通过使用SparseDirectUMFPACK类来做得更好，它不仅能产生一个不完整的，而且是一个完整的矩阵分解，这应该会产生一个更好的预处理程序。

最后，为什么使用近似值 $B^T (\textrm{diag}(M^u(S)))^{-1} B$ 来预设 $B^T M^u(S)^{-1} B$ ？后者的矩阵毕竟是压力空间上拉普拉斯算子的混合形式，我们对其使用线性元素。因此，我们可以在直接对应于拉普拉斯的非混合形式的一侧建立一个单独的矩阵 $A^p$ ，例如使用双线性形式 $(\mathbf{K}\lambda(S^n) \nabla
\varphi_i,\nabla\varphi_j)$  。然后我们可以形成这个非混合矩阵的不完全或完全分解，并将其作为混合形式的预处理。

使用这样的技术，可以合理地预期，求解过程将至少快一个数量级。




<a name="Timestepping"></a><h4>Time stepping</h4>


在介绍中，我们已经确定了时间步长的限制

@f[
  \triangle t_{n+1} \le \frac h{|\mathbf{u}^{n+1}(\mathbf{x})|}


@f]

必须是全局性的，即对所有的 $\mathbf x$ 。离散化后，我们通过选择以下方式来满足它

@f[
  \triangle t_{n+1} = \frac {\min_K h_K}{\max_{\mathbf{x}}|\mathbf{u}^{n+1}(\mathbf{x})|}.


@f]



这种对时间步长的限制有些烦人：我们把网格做得越细，时间步长就越小；换句话说，我们受到了两次惩罚：每个时间步长的求解成本更高，我们必须做更多的时间步长。

这一点特别令人讨厌，因为大部分的额外工作是用来解决方程的隐含部分，即压力-速度系统，而正是饱和度的双曲传输方程施加了时间步长的限制。

为了避免这个瓶颈，人们发明了一些方法。例如，他们可能每隔几步时间才重新计算压力-速度场（或者，如果你愿意，对压力/速度和饱和度方程使用不同的时间步长）。这就保持了对廉价显式部分的时间步长限制，而使隐式部分的求解不那么频繁。这个方向的实验当然是值得的；这种方法的一个起点是陈章新、桓冠仁和李宝岩的论文：<i>An improved IMPES method for
two-phase flow in porous media</i>，Transport in Porous Media，54（2004），第361&mdash；376页。当然也有很多其他关于这个主题的论文，但这篇论文前段时间刚好落在我们的桌上。




<a name="Adaptivity"></a><h4>Adaptivity</h4>


适应性显然也会有帮助。看一下电影，我们可以清楚地看到，大部分的行动都局限于领域的一个相对较小的部分（这对饱和度来说特别明显，但对速度和压力也是如此）。因此，自适应性可望保持必要的低自由度数量，或者增加精确度。

另一方面，对于时间相关问题的自适应性也不是一件小事：我们必须每隔几个时间步数改变网格，而且每次改变网格时，我们都必须将目前的解决方案传送到下一个网格（SolutionTransfer类可以帮助解决这个问题）。这些并不是无法克服的障碍，但它们确实需要一些额外的编码，而且比我们认为值得打包到这个教程程序中的更多。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-21.cc"
*/
