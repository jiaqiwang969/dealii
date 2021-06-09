/**
@page step_31 The step-31 tutorial program
This tutorial depends on step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#TheBoussinesqequations">The Boussinesq equations</a>
        <li><a href="#Boundaryandinitialconditions">Boundary and initial conditions</a>
        <li><a href="#Solutionapproach">Solution approach</a>
      <ul>
        <li><a href="#Timestepping">Time stepping</a>
        <li><a href="#WeakformandspacediscretizationfortheStokespart">Weak form and space discretization for the Stokes part</a>
        <li><a href="#Stabilizationweakformandspacediscretizationforthetemperatureequation">Stabilization, weak form and space discretization for the temperature equation</a>
        <li><a href="#Linearsolvers">Linear solvers</a>
      <ul>
        <li><a href="#LinearsolversfortheStokesproblem">Linear solvers for the Stokes problem</a>
        <li><a href="#Linearsolversforthetemperatureequation">Linear solvers for the temperature equation</a>
      </ul>
      </ul>
        <li><a href="#Implementationdetails">Implementation details</a>
      <ul>
        <li><a href="#UsingdifferentDoFHandlerobjects">Using different DoFHandler objects</a>
        <li><a href="#UsingTrilinos">Using Trilinos</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
      <ul>
        <li><a href="#ThecodeInverseMatrixcodeclasstemplate">The <code>InverseMatrix</code> class template</a>
        <li><a href="#Schurcomplementpreconditioner">Schur complement preconditioner</a>
      </ul>
        <li><a href="#ThecodeBoussinesqFlowProblemcodeclasstemplate">The <code>BoussinesqFlowProblem</code> class template</a>
        <li><a href="#BoussinesqFlowProblemclassimplementation">BoussinesqFlowProblem class implementation</a>
      <ul>
        <li><a href="#BoussinesqFlowProblemBoussinesqFlowProblem">BoussinesqFlowProblem::BoussinesqFlowProblem</a>
        <li><a href="#BoussinesqFlowProblemget_maximal_velocity">BoussinesqFlowProblem::get_maximal_velocity</a>
        <li><a href="#BoussinesqFlowProblemget_extrapolated_temperature_range">BoussinesqFlowProblem::get_extrapolated_temperature_range</a>
        <li><a href="#BoussinesqFlowProblemcompute_viscosity">BoussinesqFlowProblem::compute_viscosity</a>
        <li><a href="#BoussinesqFlowProblemsetup_dofs">BoussinesqFlowProblem::setup_dofs</a>
        <li><a href="#BoussinesqFlowProblemassemble_stokes_preconditioner">BoussinesqFlowProblem::assemble_stokes_preconditioner</a>
        <li><a href="#BoussinesqFlowProblembuild_stokes_preconditioner">BoussinesqFlowProblem::build_stokes_preconditioner</a>
        <li><a href="#BoussinesqFlowProblemassemble_stokes_system">BoussinesqFlowProblem::assemble_stokes_system</a>
        <li><a href="#BoussinesqFlowProblemassemble_temperature_matrix">BoussinesqFlowProblem::assemble_temperature_matrix</a>
        <li><a href="#BoussinesqFlowProblemassemble_temperature_system">BoussinesqFlowProblem::assemble_temperature_system</a>
        <li><a href="#BoussinesqFlowProblemsolve">BoussinesqFlowProblem::solve</a>
        <li><a href="#BoussinesqFlowProblemoutput_results">BoussinesqFlowProblem::output_results</a>
        <li><a href="#BoussinesqFlowProblemrefine_mesh">BoussinesqFlowProblem::refine_mesh</a>
        <li><a href="#BoussinesqFlowProblemrun">BoussinesqFlowProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Resultsin2d"> Results in 2d </a>
        <li><a href="#Resultsin3d"> Results in 3d </a>
        <li><a href="#Numericalexperimentstodetermineoptimalparameters"> Numerical experiments to determine optimal parameters </a>
      <ul>
        <li><a href="#Choosingicsubksubiicsubksubiandbeta"> Choosing <i>c<sub>k</sub></i><i>c<sub>k</sub></i> and beta </a>
      <ul>
        <li><a href="#ResultsforQsub1subelements">Results for Q<sub>1</sub> elements</a>
        <li><a href="#ResultsforQsub2subelements">Results for Q<sub>2</sub> elements</a>
        <li><a href="#Resultsfor3d">Results for 3d</a>
        <li><a href="#Conclusions">Conclusions</a>
      </ul>
      </ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-31/doc/intro.dox

 <br> 

<i>This program was contributed by Martin Kronbichler and Wolfgang
Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation or of The
California Institute of Technology.
</i>


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="TheBoussinesqequations"></a><h3>The Boussinesq equations</h3>


这个程序涉及一个有趣的物理问题：如果流体（即液体或气体）遇到由温度差异引起的浮力差异，它是如何表现的？很明显，流体中温度较高（因此较轻）的部分会上升，温度较低（密度较大）的部分会在重力作用下下沉。

在流体运动速度足够慢，以至于惯性效应可以被忽略的情况下，描述这种行为的方程是布西尼斯克方程，其内容如下。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&


  -\rho\; \beta \; T\; \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0,
  \\
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma.


@f}

这些方程属于矢量值问题的范畴（这个主题的顶层概述可以在 @ref vector_valued 模块中找到）。这里， $\mathbf u$ 是速度场， $p$ 是压力， $T$ 是流体的温度。   $\varepsilon ({\mathbf u}) = \frac 12
[(\nabla{\mathbf u}) + (\nabla {\mathbf u})^T]$ 是速度的对称梯度。可以看出，速度和压力解决了描述不可压缩流体运动的斯托克斯方程，这个方程我们以前在步骤22中考虑过；我们将广泛借鉴在该程序中获得的经验，特别是关于高效线性斯托克斯求解器的经验。

流体运动的强制项是流体的浮力，表示为密度 $\rho$ 、热膨胀系数 $\beta$ 、温度 $T$ 和指向下方的重力矢量 $\mathbf{g}$ 的积。在第32步的介绍中给出了为什么右手边看起来像它的推导）。前两个方程描述了流体如何通过移动对温差做出反应，第三个方程说明了流体运动如何影响温度场：它是一个平流扩散方程，即温度附着在流体颗粒上，并在流场中平流，还有一个额外的扩散（热传导）项。在许多应用中，扩散系数相当小，温度方程实际上是传输的，而不是扩散主导的，因此其特征是双曲而不是椭圆；我们在开发一个稳定的离散化时必须考虑到这一点。

在上述方程中，右侧的 $\gamma$ 项表示热源，可能是一个空间和时间上的变化函数。   $\eta$ 和 $\kappa$ 表示粘度和扩散系数，在本教程程序中我们假定这两个系数为常数。当 $\eta$ 取决于温度时，更普遍的情况是物理应用中的一个重要因素。大多数材料随着温度的升高而变得更加流动（即 $\eta$ 随着 $T$ 的降低而降低）；有时，如在温度接近熔点的岩石矿物的情况下， $\eta$ 可能在典型的温度范围内发生数量级的变化。

我们注意到，上述斯托克斯方程可以通过引入<a target="_top"
href="http://en.wikipedia.org/wiki/Rayleigh_number">Rayleigh
number</a>  $\mathrm{Ra}=\frac{\|\mathbf{g}\| \beta \rho}{\eta \kappa} \delta T L^3$ 来实现非维度化，使用的是典型长度尺度 $L$ 、典型温差 $\delta T$ 、密度 $\rho$ 、热扩散率 $\eta$ 和热导率 $\kappa$  。   $\mathrm{Ra}$ 是一个无尺寸的数字，它描述了由温差引起的浮力变化导致的热传输和热扩散导致的热传输的比率。一个小的瑞利数意味着浮力相对于粘度来说并不强，流体运动 $\mathbf{u}$ 足够慢，因此热扩散 $\kappa\nabla T$ 是主要的热传输项。另一方面，高瑞利数的流体将显示出主导热传导的强烈对流。

对于我们感兴趣的计算热对流的大多数流体，瑞利数是非常大的，通常是 $10^6$ 或更大。从方程的结构中，我们看到这将导致大的压力差和大的速度。因此， $T$ 的对流-扩散方程中的对流项也将非常大，这个方程的精确解将要求我们选择小的时间步长。因此，具有大雷利数的问题很难用数值来解决，其原因与<a
href="http://en.wikipedia.org/wiki/Navier-stokes_equations">Navier-Stokes
equations</a>大时难以解决<a
href="http://en.wikipedia.org/wiki/Reynolds_number">Reynolds number
$\mathrm{Re}$</a>的问题相似。

请注意，大的瑞利数不一定涉及大的绝对速度。例如，地幔中的瑞利数大于 $10^6$  。然而，速度却很小：该材料实际上是固体岩石，但它是如此之热，而且处于压力之下，它可以非常缓慢地流动，每年最多只有几厘米的速度。然而，这可以导致在数百万年的时间尺度上的混合，这个时间尺度比相同数量的热量通过热传导分布要短得多，而且这个时间尺度与影响地球内部和表面结构的演变有关。

 @note 如果你对使用该程序作为你自己实验的基础感兴趣，你也会想看看它在step-32中的延续。此外，step-32后来被发展成更大的开放源代码ASPECT（见https://aspect.geodynamics.org/），它可以解决现实的问题，在试图将step-31变形为可以解决任何你想解决的问题之前，你可能想研究一下它。




<a name="Boundaryandinitialconditions"></a><h3>Boundary and initial conditions</h3>


由于Boussinesq方程是在流体运动的惯性不起作用的假设下推导出来的，所以流场在每个时间段完全由该时间段的浮力差决定，而不是由以前的流场决定。这反映在上面的前两个方程是不包含时间导数的稳态斯托克斯方程的事实。因此，我们不需要速度或压力的初始条件。另一方面，温度场确实满足一个有时间导数的方程，所以我们需要初始条件 $T$  。

至于边界条件：如果 $\kappa>0$ ，那么温度满足一个二阶微分方程，需要边界周围所有时间的边界数据。这些数据可以是规定的边界温度 $T|_{\partial\Omega}=T_b$ （Dirichlet边界条件），也可以是规定的热通量 $\mathbf{n}\cdot\kappa\nabla
T|_{\partial\Omega}=\phi$ ；在这个程序中，我们将使用一个绝缘的边界条件，即规定没有热通量。   $\phi=0$  .

同样地，速度场要求我们提出边界条件。这些条件可以是 $\mathbf{u}=0$ 上的无滑移无通量条件 $\partial\Omega$ ，如果流体粘在边界上，或者无正常通量条件 $\mathbf n \cdot \mathbf
u = 0$ ，如果流体可以沿边界流动但不能穿过边界，或者任何数量的其他物理上合理的条件。在这个程序中，我们将使用无正常通量条件。




<a name="Solutionapproach"></a><h3>Solution approach</h3>


与步骤21中解决的方程一样，我们这里有一个微分代数方程（DAE）系统：就时间变量而言，只有温度方程是微分方程，而 $\mathbf{u}$ 和 $p$ 的斯托克斯系统没有时间导数，因此属于必须在每个时间瞬间保持的那种代数约束。与第21步的主要区别是，那里的代数约束是一个混合拉普拉斯系统，其形式为

@f{eqnarray*}
  \mathbf u + {\mathbf K}\lambda \nabla p &=& 0, \\
  \nabla\cdot \mathbf u &=& f,


@f}

现在我们有一个斯托克斯系统

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=& f, \\
  \nabla\cdot \mathbf u &=& 0,


@f}

其中 $\nabla \cdot \eta \varepsilon (\cdot)$ 是一个类似于拉普拉斯 $\Delta$ 的算子，适用于一个矢量场。

鉴于与我们在步骤21中所做的相似，我们选择类似的方法可能并不令人惊讶，尽管我们将不得不对微分算子左上角的算子变化进行调整。




<a name="Timestepping"></a><h4>Time stepping</h4>


作为DAE的问题结构允许我们使用与我们在步骤21中已经使用的相同的策略，即我们使用一个时间滞后方案：我们首先解决温度方程（使用外推的速度场），然后将新的温度解插入速度方程的右侧。不过，我们在代码中实现这一方案的方式是从一个稍微不同的角度来看问题。我们首先使用前一个时间步长的温度场来求解速度和压力的斯托克斯方程，这意味着我们得到前一个时间步长的速度。换句话说，我们首先求解时间步长 $n - 1$ 的斯托克斯系统，即

@f{eqnarray*}


  -\nabla \cdot (2\eta \varepsilon ({\mathbf u}^{n-1})) + \nabla p^{n-1} &=&


  -\rho\; \beta \; T^{n-1} \mathbf{g},
  \\
  \nabla \cdot {\mathbf u}^{n-1} &=& 0,


@f}

然后用外推速度场的温度方程到时间  $n$  。

与第21步相比，我们在这里将使用一个高阶时间步进方案，即用（单边）差分商 $\frac{\frac 32 T^{n}-2T^{n-1}+\frac 12 T^{n-2}}{k}$ 取代时间导数 $\frac{\partial T}{\partial t}$ ， $k$ 为时间步长。这就得到了离散化的时间温度方程

@f{eqnarray*}
  \frac 32 T^n


  -
  k\nabla \cdot \kappa \nabla T^n
  &=&
  2 T^{n-1}


  -
  \frac 12 T^{n-2}


  -
  k(2{\mathbf u}^{n-1} - {\mathbf u}^{n-2} ) \cdot \nabla (2T^{n-1}-T^{n-2})
  +
  k\gamma.


@f}

请注意温度方程是如何被半显式解决的：扩散被隐式处理，而平流被显式处理，使用温度和速度的外推法（或前推法），包括刚刚计算的速度  ${\mathbf u}^{n-1}$  。对当前时间水平的正向投影  $n$  是由泰勒扩展得出的，  $T^n
\approx T^{n-1} + k_n \frac{\partial T}{\partial t} \approx T^{n-1} + k_n
\frac{T^{n-1}-T^{n-2}}{k_n} = 2T^{n-1}-T^{n-2}$  。我们需要这个投影来保持BDF-2方案的精度。换句话说，我们在显式右手边使用的温度场是当前温度场的二阶近似值&mdash；不完全是显式时间步进方案，但从特征上看也不会太远。

温度外推的引入将时间步长限制在<a href="http://en.wikipedia.org/wiki/Courant–Friedrichs–Lewy_condition">
Courant-Friedrichs-Lewy (CFL) condition</a>，就像在 @ref step_21 "步骤-21 "中一样。(如果我们隐含地处理平流项，我们就不会有这个稳定条件，因为BDF-2方案是A级稳定的，代价是我们需要在每个时间步长建立一个新的温度矩阵。)我们将在<a href="#Results">results
section</a>中讨论时间步长的确切选择，但目前重要的是，这个CFL条件意味着时间步长 $k$ 可能在不同的时间步长中发生变化，我们必须稍微修改上述公式。如果 $k_n,k_{n-1}$ 是当前和前一个时间步长的时间步长，那么我们使用近似值

@f{align*}{
\frac{\partial T}{\partial t} \approx
 \frac 1{k_n}
 \left(
       \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^{n}


       -
       \frac{k_n+k_{n-1}}{k_{n-1}}T^{n-1}
       +
       \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}
 \right)
 @f}

和

@f{align*}{
T^n \approx
   T^{n-1} + k_n \frac{\partial T}{\partial t}
   \approx
   T^{n-1} + k_n
   \frac{T^{n-1}-T^{n-2}}{k_{n-1}}
   =
   \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2},


@f}

并将上述方程概括如下。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^n


  -
  k_n\nabla \cdot \kappa \nabla T^n
  &=&
  \frac{k_n+k_{n-1}}{k_{n-1}} T^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}


  -
  k_n{\mathbf u}^{*,n} \cdot \nabla T^{*,n}
  +
  k_n\gamma,


@f}



其中 ${(\cdot)}^{*,n} = \left(1+\frac{k_n}{k_{n-1}}\right)(\cdot)^{n-1} -
\frac{k_n}{k_{n-1}}(\cdot)^{n-2}$ 表示速度 $\mathbf u$ 和温度 $T$ 外推到时间级别 $n$ ，使用前两个时间步骤的数值。这不是一个容易读懂的方程，但会为我们提供所需的高阶精度。作为一致性检查，很容易验证，如果 $k_n=k_{n-1}$  ，它可以还原成与上面相同的方程。

最后我们注意到，选择高阶时间步进方案当然会迫使我们在内存中保留更多的时间步进；特别是，我们在这里需要保留 $T^{n-2}$ ，这是一个我们以前可以抛弃的向量。这似乎是一个麻烦，我们以前可以通过使用一阶时间步进方案来避免，但是正如我们在下面讨论稳定化问题时看到的那样，我们无论如何都需要这个向量，因此在时间离散化中保留它基本上是免费的，并给我们提供了使用高阶方案的机会。




<a name="WeakformandspacediscretizationfortheStokespart"></a><h4>Weak form and space discretization for the Stokes part</h4>


像解决混合拉普拉斯方程一样，解决斯托克斯方程需要我们为速度和压力变量选择特定的有限元对。因为这在步骤22中已经讨论过了，所以我们只简单介绍一下这个话题。这里，我们使用稳定对 $Q_{p+1}^d \times Q_p, p\ge 1$  。这些都是连续元素，所以我们可以通过部分积分和用离散函数替代连续函数来形成斯托克斯方程的弱形式，没有问题。

@f{eqnarray*}
  (\nabla {\mathbf v}_h, 2\eta \varepsilon ({\mathbf u}^{n-1}_h))


  -
  (\nabla \cdot {\mathbf v}_h, p^{n-1}_h)
  &=&


  -({\mathbf v}_h, \rho\; \beta \; T^{n-1}_h \mathbf{g}),
  \\
  (q_h, \nabla \cdot {\mathbf u}^{n-1}_h) &=& 0,


@f}

为所有测试函数  $\mathbf v_h, q_h$  。第一个方程的第一项被认为是张量之间的内积，即 $(\nabla {\mathbf v}_h, \eta \varepsilon ({\mathbf u}^{n-1}_h))_\Omega
 = \int_\Omega \sum_{i,j=1}^d [\nabla {\mathbf v}_h]_{ij}
           \eta [\varepsilon ({\mathbf u}^{n-1}_h)]_{ij}\, dx$  。因为这个乘积中的第二个张量是对称的，所以 $\nabla {\mathbf v}_h$ 的反对称分量不起作用，如果我们用 $\mathbf v_h$ 的对称梯度代替，会导致完全一样的形式。因此，我们考虑并实施的表述是

@f{eqnarray*}
  (\varepsilon({\mathbf v}_h), 2\eta \varepsilon ({\mathbf u}^{n-1}_h))


  -
  (\nabla \cdot {\mathbf v}_h, p^{n-1}_h)
  &=&


  -({\mathbf v}_h, \rho\; \beta \; T^{n-1}_h \mathbf{g}),
  \\
  (q_h, \nabla \cdot {\mathbf u}^{n-1}_h) &=& 0.


@f}



这与我们在第22步中已经讨论过的完全一样，这里就不多说了。




<a name="Stabilizationweakformandspacediscretizationforthetemperatureequation"></a><h4>Stabilization, weak form and space discretization for the temperature equation</h4>


更有趣的问题是如何处理温度平流-扩散方程。默认情况下，并不是所有这个方程的离散化都是同样稳定的，除非我们要么做一些像上卷、稳定化，或者所有这些的事情。实现这一点的方法之一是使用不连续元素（即我们在步骤12中离散传输方程或在步骤20和步骤21中离散压力时使用的FE_DGQ类），并在单元间的界面上定义一个考虑到上卷的流量。如果我们有一个纯粹的平流问题，这可能是最简单的方法。然而，这里我们也有一些扩散，用不连续元素对拉普拉斯算子进行离散化是很麻烦的，因为有大量的附加项需要在单元间的每个面上进行积分。不连续元素还有一个缺点，即使用数值通量会带来额外的数值扩散，这种扩散无处不在，而我们真的希望将数值扩散的影响降到最低，只在需要稳定方案的地方应用它。

因此，一个更好的选择是在模型中加入一些非线性粘度。从本质上讲，这样做的目的是将温度方程的形式从

@f{eqnarray*}
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma


@f}

到类似于

@f{eqnarray*}
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot (\kappa+\nu(T)) \nabla T &=& \gamma,


@f}

其中 $\nu(T)$ 是一个额外的粘度（扩散）项，只在冲击和其他不连续点附近发挥作用。   $\nu(T)$ 的选择方式是，如果 $T$ 满足原始方程，则额外的粘性为零。

为了实现这一点，文献中包含了许多方法。我们在这里将遵循Guermond和Popov开发的一种方法，它建立在一个适当定义的残差和一个额外粘度的极限程序之上。为此，让我们定义一个残差 $R_\alpha(T)$ 如下。

@f{eqnarray*}
  R_\alpha(T)
  =
  \left(
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T - \gamma
  \right)
  T^{\alpha-1}


@f}

其中，我们以后将从 $[1,2]$ 范围内选择稳定指数 $\alpha$ 。请注意，如果 $T$ 满足温度方程， $R_\alpha(T)$ 将为零，因为此时括号内的项将为零。将条款相乘，我们得到以下完全等同的形式。

@f{eqnarray*}
  R_\alpha(T)
  =
  \frac 1\alpha
  \frac{\partial (T^\alpha)}{\partial t}
  +
  \frac 1\alpha
  {\mathbf u} \cdot \nabla (T^\alpha)


  -
  \frac 1\alpha
  \nabla \cdot \kappa \nabla (T^\alpha)
  +
  \kappa(\alpha-1)
  T^{\alpha-2} |\nabla T|^2


  -
  \gamma
  T^{\alpha-1}


@f}



有了这个残差，我们现在可以把人工黏度定义为一个片状常数函数，在直径为 $K$ 的每个单元上分别定义如下。

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
  \min\left\{
    h_K,
    h_K^\alpha
    \frac{\|R_\alpha(T)\|_{L^\infty(K)}}{c(\mathbf{u},T)}
  \right\}


@f}



这里， $\beta$ 是一个稳定常数（通过维度分析发现它是无单位的，因此与比例无关；我们将在<a href="#Results">results section</a>中讨论其选择）， $c(\mathbf{u},T)$ 是一个归一化常数，其单位必须是 $\frac{m^{\alpha-1}K^\alpha}{s}$  。我们将选择它作为 $c(\mathbf{u},T) =
 c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
 \ |\mathrm{diam}(\Omega)|^{\alpha-2}$  ，其中 $\mathrm{var}(T)=\max_\Omega T - \min_\Omega T$ 是目前温度值的范围（记住，浮力是由温度变化驱动的，而不是绝对温度）， $c_R$ 是一个无尺寸常数。为了理解这个方法为什么有效，请考虑这个问题。如果在一个特定的单元 $K$ 上，温度场是平滑的，那么我们希望那里的残差很小（事实上是在 ${\cal O}(h_K)$ 的数量级上），注入人工扩散的稳定项在那里的大小将是 $h_K^{\alpha+1}$ &mdash；也就是说，相当小，就像我们希望它在没有必要进行额外扩散时那样。另一方面，如果我们处于或接近温度场的不连续性，那么残差将很大； $\nu_\alpha(T)$ 定义中的最小操作将确保稳定项的大小为 $h_K$ &mdash；这是确保方案稳定的最佳人工粘性量。

这种方案是否真的有效是个好问题。Guermond和Popov的计算表明，这种形式的稳定方案实际上比其他大多数稳定方案（例如流线扩散，仅举最简单的一种）表现得更好。此外，对于 $\alpha\in
[1,2)$ ，他们甚至可以证明，对于线性传输方程，它比流线扩散产生更好的收敛阶数。对于 $\alpha=2$ ，目前还没有理论结果，但数值测试表明，其结果比 $\alpha=1$ 好得多。

一个更实际的问题是如何将这种人工扩散引入我们想要解决的方程。请注意，数值粘度 $\nu(T)$ 是随温度变化的，所以我们要解决的方程在 $T$ 中是非线性的&mdash；这不是人们对稳定方程的简单方法的期望，如果我们意识到 $\nu(T)$ 在 $T$ 中是不可分的，那就更不可能了。然而，我们没有理由绝望：我们仍然要在时间上进行离散，我们可以明确地处理这个术语。

在稳定参数的定义中，我们用  $\frac{\partial T}{\partial t} \approx
\frac{T^{n-1}-T^{n-2}}{k^{n-1}}$  对时间导数进行近似。这种近似只利用了可用的时间数据，这就是我们需要存储前两个时间步骤的数据的原因（这使我们能够使用BDF-2方案而不需要额外的存储成本）。我们现在可以简单地在 $t_{n-1}$ 处评估其余的项，但这样一来，离散残差无非是一个向后的欧拉近似，它只有一阶精度。因此，在平滑解的情况下，尽管外部BDF-2方案和空间FE离散化的时间精度为二阶，但残差仍为 $h$ 阶。这当然不是我们想要的（事实上，我们希望在解决方案表现良好的区域有较小的残差），所以需要更谨慎一些。这个问题的关键是观察我们构造的第一导数实际上是以 $t_{n-\frac{3}{2}}$ 为中心的。如果我们通过使用近似值 $\frac 12 T^{n-1}+\frac 12 T^{n-2}$ 来评估 $t_{n-\frac{3}{2}}$ 处的所有空间项，我们就可以得到所需的二阶精确残差计算，这意味着我们将非线性粘度计算为这个中间温度的函数， $\nu_\alpha =
\nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)$  。请注意，这种对残差的评估无非是一个Crank-Nicholson方案，所以我们可以肯定，现在一切正常了。人们可能会想，现在的数值粘度没有在时间 $n$ 进行评估（相对于方程的其余部分），这是否是一个问题。然而，这种偏移是不严谨的。对于平滑解， $\nu_\alpha$ 将连续变化，所以时间偏移的误差比非线性粘度本身要小 $k$ 倍，也就是说，它是被遗漏的一个小的高阶贡献。这很好，因为该项本身已经达到了光滑区域的离散化误差水平。

使用上面介绍的BDF-2方案，这就得到了更简单的大小为 $k$ 的均匀时间步长的情况。

@f{eqnarray*}
  \frac 32 T^n


  -
  k\nabla \cdot \kappa \nabla T^n
  &=&
  2 T^{n-1}


  -
  \frac 12 T^{n-2}
  \\
  &&
  +
  k\nabla \cdot
  \left[
    \nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)
    \ \nabla (2T^{n-1}-T^{n-2})
  \right]
  \\
  &&


  -
  k(2{\mathbf u}^{n-1}-{\mathbf u}^{n-2}) \cdot \nabla (2T^{n-1}-T^{n-2})
  \\
  &&
  +
  k\gamma.


@f}

在这个方程的左侧仍然是来自时间导数的项和我们隐含处理的原始（物理）扩散（这实际上是一个很好的项：从左侧产生的矩阵是质量矩阵和拉普拉斯矩阵的倍数&mdash；两者都是正定的，如果时间步长 $k$ 很小，和很容易反转）。在右侧，第一行的条款是时间导数的结果；第二行是时间 $t_{n-\frac
32}$ 的人工扩散；第三行包含平流条款，第四行是来源。请注意，人工扩散对当前时间的外推温度的作用，与我们在时间步进一节中讨论的平流作用相同。

我们在现实中必须使用的非均匀时间步长的形式要复杂一些（这就是为什么我们先展示了上面的简单形式），其内容为：。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^n


  -
  k_n\nabla \cdot \kappa \nabla T^n
  &=&
  \frac{k_n+k_{n-1}}{k_{n-1}} T^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}
  \\
  &&
  +
  k_n\nabla \cdot
  \left[
    \nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)
    \ \nabla  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \right]
  \\
  &&


  -
  k_n
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right){\mathbf u}^{n-1} -
    \frac{k_n}{k_{n-1}}{\mathbf u}^{n-2}
  \right]
  \cdot \nabla
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \\
  &&
  +
  k_n\gamma.


@f}



在解决了所有这些问题之后，弱形式自然而然地从最后一个方程中显示的强形式中产生，我们立即得出了离散化方程的弱形式。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} (\tau_h,T_h^n)
  +
  k_n (\nabla \tau_h, \kappa \nabla T_h^n)
  &=&
  \biggl(\tau_h,
  \frac{k_n+k_{n-1}}{k_{n-1}} T_h^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T_h^{n-2}
  \\
  &&\qquad


  -
  k_n
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right){\mathbf u}^{n-1} -
    \frac{k_n}{k_{n-1}}{\mathbf u}^{n-2}
  \right]
  \cdot \nabla
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  +
  k_n\gamma \biggr)
  \\
  &&


  -
  k_n \left(\nabla \tau_h,
    \nu_\alpha\left(\frac 12 T_h^{n-1}+\frac 12 T_h^{n-2}\right)
    \ \nabla \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \right)


@f}

为所有离散测试函数  $\tau_h$  。在这里，扩散项已经被部分整合，我们已经使用，我们将施加没有热通量，  $\mathbf{n}\cdot\kappa\nabla T|_{\partial\Omega}=0$  。

这就产生了一个矩阵方程，其形式为

@f{eqnarray*}
  \left( \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} M+k_n A_T\right) T_h^n
  = F(U_h^{n-1}, U_h^{n-2},T_h^{n-1},T_h^{n-2}),


@f}

考虑到左边的矩阵结构（两个正定矩阵之和），使用共轭梯度法很容易解决这个问题。




<a name="Linearsolvers"></a><h4>Linear solvers</h4>


如上所述，我们解决速度/压力和温度的联合系统的方法是使用算子分割，我们首先用旧的温度场解决速度和压力的斯托克斯系统，然后用刚刚计算的速度场解决新的温度场。关于算子分割方法的更广泛的讨论可以在步骤58中找到）。




<a name="LinearsolversfortheStokesproblem"></a><h5>Linear solvers for the Stokes problem</h5>


解决来自斯托克斯系统的线性方程已经在步骤22中进行了详细的讨论。特别是在该程序的结果部分，我们讨论了一些替代的线性求解器策略，结果发现这些策略比原来的方法更有效。在那里确定的最佳替代方案是使用一个由涉及舒尔补码的块状矩阵预处理的GMRES求解器。具体来说，斯托克斯算子导致了一个块状结构的矩阵

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)


@f}

正如那里所讨论的，一个好的预处理程序是

@f{eqnarray*}
  P
  =
  \left(\begin{array}{cc}
    A & 0 \\ B & -S
  \end{array}\right),
  \qquad
  \text{or equivalently}
  \qquad
  P^{-1}
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
  \end{array}\right)


@f}

其中 $S$ 是斯托克斯算子的舒尔补 $S=B^TA^{-1}B$  。当然，这个预处理程序是没有用的，因为我们不能形成矩阵的各种倒数，但我们可以用下面的方法作为预处理程序。

@f{eqnarray*}
  \tilde P^{-1}
  =
  \left(\begin{array}{cc}
    \tilde A^{-1} & 0 \\ \tilde S^{-1} B \tilde A^{-1} & -\tilde S^{-1}
  \end{array}\right)


@f}

其中 $\tilde A^{-1},\tilde S^{-1}$ 是反矩阵的近似值。特别是，事实证明 $S$ 在光谱上等同于质量矩阵，因此，用适用于压力空间上的质量矩阵的CG求解器取代 $\tilde
S^{-1}$ 是一个不错的选择。与步骤22稍有不同的是，我们在这里的动量方程中有一个系数 $\eta$ ，通过与那里相同的推导，我们应该得出结论，我们应该使用的是具有条目 $\tilde S_{ij}=(\eta^{-1}\varphi_i,\varphi_j)$ 的加权质量矩阵。

想出一个好的替代方案 $\tilde
A^{-1}$ 更为复杂，它对应于矢量值速度场的离散化对称拉普拉斯，即 $A_{ij} = (\varepsilon {\mathbf v}_i, 2\eta \varepsilon ({\mathbf
v}_j))$  。在步骤22中，我们用 $A$ 的稀疏LU分解（使用SparseDirectUMFPACK类）来代替 $\tilde A^{-1}$ &mdash; 完美的前置条件&mdash; 在2D中，但对于3D来说，内存和计算时间通常不足以实际计算这个分解；因此，我们在3D中只使用不完全LU分解（ILU，使用稀疏ILU类）。

对于这个项目，我们想走得更远一点。为此，请注意，矢量场上的对称化双线性形式 $(\varepsilon {\mathbf v}_i, 2 \eta \varepsilon ({\mathbf v}_j))$ 与非对称化版本 $(\nabla {\mathbf v}_i, \eta \nabla {\mathbf v}_j)
= \sum_{k,l=1}^d
  (\partial_k ({\mathbf v}_i)_l, \eta \partial_k ({\mathbf v}_j)_l)
$ 相差不大（请注意，在这个形式中因子2已经消失了）。然而，后者的优点是测试函数的 <code>dim</code> 矢量分量不是耦合的（好吧，几乎是，见下文），也就是说，得到的矩阵是块对角线的：每个矢量分量有一个块，这些块中的每个都等于这个矢量分量的拉普拉斯矩阵。因此，假设我们以这样的方式排列自由度，即首先对速度的所有 $x$ 分量进行编号，然后是 $y$ 分量，然后是 $z$ 分量，那么与这种稍有不同的双线性形式相关的矩阵 $\hat A$ 具有如下形式

@f{eqnarray*}
  \hat A =
  \left(\begin{array}{ccc}
    A_s & 0 & 0 \\ 0 & A_s & 0 \\ 0 & 0 & A_s
  \end{array}\right)


@f}

其中 $A_s$ 是一个拉普拉斯矩阵，其大小等于与矢量值速度的每个分量相关的形状函数数量。有了这个矩阵，我们就可以对速度矩阵 $A$ 的预处理进行如下定义。

@f{eqnarray*}
  \tilde A^{-1} =
  \left(\begin{array}{ccc}
    \tilde A_s^{-1} & 0 & 0 \\
    0 & \tilde A_s^{-1} & 0 \\
    0 & 0 & \tilde A_s^{-1}
  \end{array}\right),


@f}

其中 $\tilde A_s^{-1}$ 是拉普拉斯矩阵的预处理程序&mdash;我们非常清楚如何建立良好的预处理程序!

在现实中，故事并不那么简单。为了使矩阵 $\tilde A$ 确定，我们需要通过应用边界条件使各个块 $\tilde
A_s$ 确定。我们可以尝试通过在边界周围应用狄氏边界条件来做到这一点，然后，如果后者的矩阵是由斯托克斯问题产生的，我们在领域周围的速度分量上也有狄氏边界条件，即如果我们执行 $\mathbf{u} =
0$ ，那么如此定义的前置条件 $\tilde A^{-1}$ 就变成了 $A$ 的良好前置条件。

不幸的是，这个 "如果 "是 "如果且仅是如果"：在下面的程序中，我们将希望使用 $\mathbf u
\cdot \mathbf n = 0$ 形式的无流量边界条件（即允许与边界平行的流量%，但没有通过边界的流量）。在这种情况下，事实证明，上面定义的块状对角线矩阵不是一个好的预处理程序，因为它忽略了边界上的成分耦合。因此，更好的方法是如果我们将矩阵 $\hat A$ 建立为矢量拉普拉斯矩阵 $\hat A_{ij} = (\nabla {\mathbf v}_i,
\eta \nabla {\mathbf v}_j)$ ，然后应用与我们应用于 $A$ 相同的边界条件。如果这是一个围绕域的迪里希特边界条件， $\hat A$ 将像上面那样解耦为三个对角线块，如果边界条件是 $\mathbf u
\cdot \mathbf n = 0$ 的形式，那么这将在边界引入自由度的耦合，但只在那里。事实上，这被证明是一个比上面介绍的更好的预处理程序，而且几乎具有我们希望得到的所有好处。


总结这整个故事，我们可以看到。   <ul>   <li>  与我们在步骤22中从对称梯度产生的原始矩阵 $A$ 建立一个预处理程序相比，我们不得不期待基于拉普拉斯双线性形式的预处理程序表现得更差，因为它没有考虑到向量分量之间的耦合。

    <li> 另一方面，拉普拉斯矩阵的预处理程序通常比矢量问题的预处理程序更成熟，性能更好。例如，在写这篇文章的时候，代数%多重网格（AMG）算法对于标量问题已经非常成熟，但对于矢量问题却不是如此。

    <li> 在建立这个预处理程序时，我们将不得不建立矩阵 $\hat A$ 及其预处理程序。虽然这意味着我们必须存储一个之前不需要的额外矩阵，但与存储耦合矩阵 $A$ 的预处理程序相比，预处理程序 $\tilde A_s^{-1}$ 可能需要的内存要少得多。这是因为矩阵 $A_s$ 每行只有三分之一的条目对应于内部自由度，并且只在边界条件引入耦合的部分包含向量分量之间的耦合。因此，存储该矩阵是比较便宜的，我们可以预期，计算和存储预处理程序 $\tilde A_s$ 也将比为完全耦合的矩阵做这些事情便宜得多。   </ul> 




<a name="Linearsolversforthetemperatureequation"></a><h5>Linear solvers for the temperature equation</h5>


这是最容易的部分。温度方程的矩阵具有 $\alpha M + \beta A$ 的形式，其中 $M,A$ 是温度空间上的质量和刚度矩阵， $\alpha,\beta$ 是与时间步进方案以及当前和前一个时间步进有关的常数。这是一个对称正定和一个对称正半定矩阵之和，其结果也是对称正定的。此外， $\frac\beta\alpha$ 是一个与时间步长成正比的数字，因此只要网格很细就会变小，从而阻尼当时条件不好的刚度矩阵的影响。

因此，用共轭梯度算法反转这个矩阵，使用一个简单的预处理程序，与反转斯托克斯矩阵相比是微不足道和非常便宜的。




<a name="Implementationdetails"></a><h3>Implementation details</h3>


<a name="UsingdifferentDoFHandlerobjects"></a><h4>Using different DoFHandler objects</h4>


关于下面的程序，值得事先解释的一件事是使用了两个不同的DoFHandler对象。如果看一下上述方程的结构和它们的求解方案，就会发现几乎没有什么共同点能使斯托克斯部分和温度部分保持一致。在我们以前讨论 @ref
vector_valued "矢量值问题 "的所有教程程序中，我们总是只使用一个具有几个矢量分量的单一有限元，以及一个DoFHandler对象。有时，我们将得到的矩阵分解成若干块，以方便特定的求解器方案；例如，在目前程序所依据的斯托克斯方程的第22步程序中就是如此。

当然，我们在这里也可以这样做。我们将得到的线性系统看起来像这样。

@f{eqnarray*}
  \left(\begin{array}{ccc}
    A & B^T & 0 \\ B & 0 &0 \\ C & 0 & K
  \end{array}\right)
  \left(\begin{array}{ccc}
    U^{n-1} \\ P^{n-1} \\ T^n
  \end{array}\right)
  =
  \left(\begin{array}{ccc}
    F_U(T^{n-1}) \\ 0 \\ F_T(U^{n-1},U^{n-2},T^{n-1},T^{n-2})
  \end{array}\right).


@f}

这方面的问题是。我们从未同时使用整个矩阵。事实上，它从未真正同时存在。如上所述， $K$ 和 $F_T$ 依赖于已经计算出的解 $U^n$ ，在第一种情况下，通过时间步长（这依赖于 $U^n$ ，因为它必须满足CFL条件）。所以我们只有在已经解决了左上角 $2\times 2$ 块斯托克斯系统后才能组装它，一旦我们转向温度方程，我们就不再需要斯托克斯部分了；我们为一个在任何时候都不会以整体存在于内存中的矩阵建立一个对象，这导致我们在步骤21中跳了一些圈套，所以我们不要重复这类错误。此外，我们实际上并没有建立矩阵 $C$ ：因为当我们进入温度方程时，我们已经知道了 $U^n$ ，而且因为我们必须在这个时候组装右手边的 $F_T$ ，我们只是将项 $CU^n$ 移到右手边，并将其与所有其他项组装在一起。这意味着矩阵中不存在温度变量和斯托克斯变量耦合的部分，因此所有自由度的全局列举不再重要：如果我们有所有斯托克斯自由度的列举，以及所有温度自由度的独立列举就足够了。

从本质上讲，将<i>everything</i>放入一个块状矩阵中并没有什么用处（当然，对于 $2\times 2$ 斯托克斯部分，也有同样好的理由这样做），或者，就这一点而言，将所有东西放入同一个DoFHandler对象。

但这样做是否有<i>downsides</i>的好处？这些问题是存在的，尽管它们一开始可能并不明显。主要问题是，如果我们需要创建一个包含速度、压力和温度形状函数的全局有限元，并使用它来初始化DoFHandler。但是我们也用这个有限元对象来初始化我们使用的所有FEValues或FEFaceValues对象。这可能看起来不是什么大问题，但是想象一下，例如，当我们评估我们需要计算人工粘度 $
  R_\alpha(T)
  =
  \left(
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T - \gamma
  \right)
  T^{\alpha-1}
$ 的残差 $\nu_\alpha(T)|_K$ 时会发生什么。  为此，我们需要温度的拉普拉斯，我们使用形状函数的二阶导数（Hessians）张量来计算（为此我们必须给FEValues对象加上 <code>update_hessians</code> 标志）。现在，如果我们有一个包含速度、压力和温度的形状函数的有限性，这意味着我们必须计算<i>all</i>形状函数的Hessians，包括速度的许多高阶形状函数。这是很多我们不需要的计算，事实上，如果一个人要这样做（就像我们在程序的早期版本中那样），组装右手边需要大约四分之一的整体计算时间。

所以我们要做的是使用两个不同的有限元对象，一个用于斯托克斯成分，一个用于温度。这样就有两个不同的DoFHandlers，两个稀疏模式和两个用于斯托克斯和温度部分的矩阵，等等。每当我们要组装包含温度和斯托克斯形状函数的东西时（特别是斯托克斯和温度方程的右侧），我们就使用两个FEValues对象，用两个单元格迭代器进行初始化，通过与同一三角化对象相关的两个DoFHandler对象进行平行行走。对于这两个FEValues对象，我们当然使用相同的正交对象，这样我们就可以在同一组正交点上进行迭代，但是每个FEValues对象将只根据它实际需要计算的内容来获得更新标志。特别是，当我们像上面那样计算残差时，我们只要求得到斯托克斯形状函数的值，但也要求得到温度形状函数的Hessians &mdash；确实便宜得多，而且事实证明：组装温度方程的右手边现在是程序中几乎无法测量的一个组成部分。

有了这些变化，对程序进行计时，可以得出只有以下操作与整个运行时间有关。   <ul>   <li>  解决斯托克斯系统：72%的运行时间。     <li>  组装斯托克斯预处理程序，并使用Trilinos ML包计算代数多网格层次结构：占运行时间的11%。     <li>  函数  <code>BoussinesqFlowProblem::setup_dofs</code>  : 占整体运行时间的7%。     <li>  组装斯托克斯和温度右侧向量以及组装矩阵。7%.   </ul>  实质上这意味着除了代数多重网格之外，所有的瓶颈都已经被移除。




<a name="UsingTrilinos"></a><h4>Using Trilinos</h4>


与我们在第17步和第18步中使用PETSc来支持我们的线性代数需求一样，我们在这个程序中使用了<a
href="http://trilinos.org">Trilinos</a>库的接口（安装说明见deal.II README文件）。Trilinos是一个非常大的集合，包括与线性和非线性代数有关的所有东西，以及围绕这些东西的各种工具（看起来它在未来也会向许多其他方向发展）。

使用Trilinos的主要原因，类似于我们探索的PETSc，是它是一个非常强大的库，比deal.II自己的线性代数库提供了很多工具。这尤其包括在集群上以%parallel方式工作的能力，使用MPI，以及更多种类的前置条件器。在后一类中，最有趣的能力之一是Trilinos ML包的存在，它实现了代数多栅（AMG）方法。我们将使用这个预处理程序对动量方程的二阶算子部分进行预处理。在步骤32中，我们将使用与这里讨论的相同的问题，探索以%并行方式解决问题的能力。

我们在第17步和第18步中使用的PETSc无疑是一个强大的库，它提供了大量处理矩阵、向量、迭代求解器和预处理器的函数，还有很多其他的东西，其中大部分在%parallel中运行得相当好。然而，它比Trilinos早了几年，是用C语言编写的，而且一般来说不像其他一些库那样容易使用。因此，deal.II也获得了与Trilinos的接口，Trilinos与PETSc有很多相同的功能。然而，它是一个年轻了好几年的项目，是用C++编写的，其作者一般都非常重视软件设计。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们在这里要解决的情况如下：我们用 $\kappa=10^{-6}, \eta=1, \rho=1, \beta=10$ 来解决上述的Boussinesq方程，即一个相对缓慢运动的流体，它几乎没有热扩散传导性，主要通过对流来传输热量。在边界上，我们将要求速度（ $\mathrm{n}\cdot\mathrm{u}=0$ ）和温度（ $\mathrm{n}\cdot\nabla T=0$ ）没有正态流量。这是在步骤22的介绍中讨论的情况之一，它固定了速度的一个分量，同时允许流动与边界%平行。还有 <code>dim-1</code> 分量需要固定，即法向应力的切向分量；对于这些分量，我们选择同质条件，这意味着我们不需要任何特殊条件。初始条件只对温度场是必要的，我们选择它为恒定的零。

然后，问题的演变完全由温度方程的右手边 $\gamma(\mathrm{x},t)$ 驱动，即由热源和汇驱动。在这里，我们选择了一个在圣诞讲座前发明的设置：美国的教室里当然禁止使用真实的蜡烛，但允许使用虚拟的蜡烛。因此，我们选择了三个球形的热源，不等距地靠近领域的底部，模仿三个蜡烛的样子。位于这些热源处的流体，最初处于静止状态，然后被加热，随着温度的升高，获得浮力，上升；更多的流体被拖上来，穿过热源，导致三个热羽上升，直到它们被外面下沉的流体循环所捕获，取代了因加热而上升的空气。


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
 * The first step, as always, is to include the functionality of these
 * well-known deal.II library files and some C++ header files.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/block_sparsity_pattern.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * @endcode
 * 
 * Then we need to include some header files that provide vector, matrix, and
 * preconditioner classes that implement interfaces to the respective Trilinos
 * classes. In particular, we will need interfaces to the matrix and vector
 * classes based on Trilinos as well as Trilinos preconditioners:
 * 
 * @code
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_vector.h>
 * #include <deal.II/lac/trilinos_parallel_block_vector.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * 
 * @endcode
 * 
 * Finally, here are a few C++ headers that haven't been included yet by one of
 * the aforelisted header files:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * #include <memory>
 * #include <limits>
 * 
 * 
 * @endcode
 * 
 * At the end of this top-matter, we import all deal.II names into the global
 * namespace:
 * 
 * @code
 * namespace Step31
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * Again, the next stage in the program is the definition of the equation
 * data, that is, the various boundary conditions, the right hand sides and
 * the initial condition (remember that we're about to solve a
 * time-dependent system). The basic strategy for this definition is the
 * same as in step-22. Regarding the details, though, there are some
 * differences.
 * 

 * 
 * The first thing is that we don't set any inhomogeneous boundary
 * conditions on the velocity, since as is explained in the introduction we
 * will use no-flux conditions $\mathbf{n}\cdot\mathbf{u}=0$. So what is
 * left are <code>dim-1</code> conditions for the tangential part of the
 * normal component of the stress tensor, $\textbf{n} \cdot [p \textbf{1} -
 * \eta\varepsilon(\textbf{u})]$; we assume homogeneous values for these
 * components, i.e., a natural boundary condition that requires no specific
 * action (it appears as a zero term in the right hand side of the weak
 * form).
 *   

 * 
 * For the temperature $T$, we assume no thermal energy flux,
 * i.e., $\mathbf{n} \cdot \kappa \nabla T=0$. This, again, is a boundary
 * condition that does not require us to do anything in particular.
 *   

 * 
 * Secondly, we have to set initial conditions for the temperature (no
 * initial conditions are required for the velocity and pressure, since the
 * Stokes equations for the quasi-stationary case we consider here have no
 * time derivatives of the velocity or pressure). Here, we choose a very
 * simple test case, where the initial temperature is zero, and all dynamics
 * are driven by the temperature right hand side.
 *   

 * 
 * Thirdly, we need to define the right hand side of the temperature
 * equation. We choose it to be constant within three circles (or spheres in
 * 3d) somewhere at the bottom of the domain, as explained in the
 * introduction, and zero outside.
 *   

 * 
 * Finally, or maybe firstly, at the top of this namespace, we define the
 * various material constants we need ($\eta,\kappa$, density $\rho$ and the
 * thermal expansion coefficient $\beta$):
 * 
 * @code
 *   namespace EquationData
 *   {
 *     constexpr double eta     = 1;
 *     constexpr double kappa   = 1e-6;
 *     constexpr double beta    = 10;
 *     constexpr double density = 1;
 * 
 * 
 *     template <int dim>
 *     class TemperatureInitialValues : public Function<dim>
 *     {
 *     public:
 *       TemperatureInitialValues()
 *         : Function<dim>(1)
 *       {}
 * 
 *       virtual double value(const Point<dim> & /*p*/,
 *                            const unsigned int /*component*/ = 0) const override
 *       {
 *         return 0;
 *       }
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  value) const override
 *       {
 *         for (unsigned int c = 0; c < this->n_components; ++c)
 *           value(c) = TemperatureInitialValues<dim>::value(p, c);
 *       }
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     class TemperatureRightHandSide : public Function<dim>
 *     {
 *     public:
 *       TemperatureRightHandSide()
 *         : Function<dim>(1)
 *       {}
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override
 *       {
 *         (void)component;
 *         Assert(component == 0,
 *                ExcMessage("Invalid operation for a scalar function."));
 * 
 *         Assert((dim == 2) || (dim == 3), ExcNotImplemented());
 * 
 *         static const Point<dim> source_centers[3] = {
 *           (dim == 2 ? Point<dim>(.3, .1) : Point<dim>(.3, .5, .1)),
 *           (dim == 2 ? Point<dim>(.45, .1) : Point<dim>(.45, .5, .1)),
 *           (dim == 2 ? Point<dim>(.75, .1) : Point<dim>(.75, .5, .1))};
 *         static const double source_radius = (dim == 2 ? 1. / 32 : 1. / 8);
 * 
 *         return ((source_centers[0].distance(p) < source_radius) ||
 *                     (source_centers[1].distance(p) < source_radius) ||
 *                     (source_centers[2].distance(p) < source_radius) ?
 *                   1 :
 *                   0);
 *       }
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  value) const override
 *       {
 *         for (unsigned int c = 0; c < this->n_components; ++c)
 *           value(c) = TemperatureRightHandSide<dim>::value(p, c);
 *       }
 *     };
 *   } // namespace EquationData
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
 * This section introduces some objects that are used for the solution of
 * the linear equations of the Stokes system that we need to solve in each
 * time step. Many of the ideas used here are the same as in step-20, where
 * Schur complement based preconditioners and solvers have been introduced,
 * with the actual interface taken from step-22 (in particular the
 * discussion in the "Results" section of step-22, in which we introduce
 * alternatives to the direct Schur complement approach). Note, however,
 * that here we don't use the Schur complement to solve the Stokes
 * equations, though an approximate Schur complement (the mass matrix on the
 * pressure space) appears in the preconditioner.
 * 
 * @code
 *   namespace LinearSolvers
 *   {
 * @endcode
 * 
 * 
 * <a name="ThecodeInverseMatrixcodeclasstemplate"></a> 
 * <h4>The <code>InverseMatrix</code> class template</h4>
 * 

 * 
 * This class is an interface to calculate the action of an "inverted"
 * matrix on a vector (using the <code>vmult</code> operation) in the same
 * way as the corresponding class in step-22: when the product of an
 * object of this class is requested, we solve a linear equation system
 * with that matrix using the CG method, accelerated by a preconditioner
 * of (templated) class <code>PreconditionerType</code>.
 *     

 * 
 * In a minor deviation from the implementation of the same class in
 * step-22, we make the <code>vmult</code> function take any
 * kind of vector type (it will yield compiler errors, however, if the
 * matrix does not allow a matrix-vector product with this kind of
 * vector).
 *     

 * 
 * Secondly, we catch any exceptions that the solver may have thrown. The
 * reason is as follows: When debugging a program like this one
 * occasionally makes a mistake of passing an indefinite or nonsymmetric
 * matrix or preconditioner to the current class. The solver will, in that
 * case, not converge and throw a run-time exception. If not caught here
 * it will propagate up the call stack and may end up in
 * <code>main()</code> where we output an error message that will say that
 * the CG solver failed. The question then becomes: Which CG solver? The
 * one that inverted the mass matrix? The one that inverted the top left
 * block with the Laplace operator? Or a CG solver in one of the several
 * other nested places where we use linear solvers in the current code? No
 * indication about this is present in a run-time exception because it
 * doesn't store the stack of calls through which we got to the place
 * where the exception was generated.
 *     

 * 
 * So rather than letting the exception propagate freely up to
 * <code>main()</code> we realize that there is little that an outer
 * function can do if the inner solver fails and rather convert the
 * run-time exception into an assertion that fails and triggers a call to
 * <code>abort()</code>, allowing us to trace back in a debugger how we
 * got to the current place.
 * 
 * @code
 *     template <class MatrixType, class PreconditionerType>
 *     class InverseMatrix : public Subscriptor
 *     {
 *     public:
 *       InverseMatrix(const MatrixType &        m,
 *                     const PreconditionerType &preconditioner);
 * 
 * 
 *       template <typename VectorType>
 *       void vmult(VectorType &dst, const VectorType &src) const;
 * 
 *     private:
 *       const SmartPointer<const MatrixType> matrix;
 *       const PreconditionerType &           preconditioner;
 *     };
 * 
 * 
 *     template <class MatrixType, class PreconditionerType>
 *     InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(
 *       const MatrixType &        m,
 *       const PreconditionerType &preconditioner)
 *       : matrix(&m)
 *       , preconditioner(preconditioner)
 *     {}
 * 
 * 
 * 
 *     template <class MatrixType, class PreconditionerType>
 *     template <typename VectorType>
 *     void InverseMatrix<MatrixType, PreconditionerType>::vmult(
 *       VectorType &      dst,
 *       const VectorType &src) const
 *     {
 *       SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm());
 *       SolverCG<VectorType> cg(solver_control);
 * 
 *       dst = 0;
 * 
 *       try
 *         {
 *           cg.solve(*matrix, dst, src, preconditioner);
 *         }
 *       catch (std::exception &e)
 *         {
 *           Assert(false, ExcMessage(e.what()));
 *         }
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Schurcomplementpreconditioner"></a> 
 * <h4>Schur complement preconditioner</h4>
 * 

 * 
 * This is the implementation of the Schur complement preconditioner as
 * described in detail in the introduction. As opposed to step-20 and
 * step-22, we solve the block system all-at-once using GMRES, and use the
 * Schur complement of the block structured matrix to build a good
 * preconditioner instead.
 *     

 * 
 * Let's have a look at the ideal preconditioner matrix
 * $P=\left(\begin{array}{cc} A & 0 \\ B & -S \end{array}\right)$
 * described in the introduction. If we apply this matrix in the solution
 * of a linear system, convergence of an iterative GMRES solver will be
 * governed by the matrix @f{eqnarray*} P^{-1}\left(\begin{array}{cc} A &
 * B^T \\ B & 0 \end{array}\right) = \left(\begin{array}{cc} I & A^{-1}
 * B^T \\ 0 & I \end{array}\right), @f} which indeed is very simple. A
 * GMRES solver based on exact matrices would converge in one iteration,
 * since all eigenvalues are equal (any Krylov method takes at most as
 * many iterations as there are distinct eigenvalues). Such a
 * preconditioner for the blocked Stokes system has been proposed by
 * Silvester and Wathen ("Fast iterative solution of stabilised Stokes
 * systems part II.  Using general block preconditioners", SIAM
 * J. Numer. Anal., 31 (1994), pp. 1352-1367).
 *     

 * 
 * Replacing $P$ by $\tilde{P}$ keeps that spirit alive: the product
 * $P^{-1} A$ will still be close to a matrix with eigenvalues 1 with a
 * distribution that does not depend on the problem size. This lets us
 * hope to be able to get a number of GMRES iterations that is
 * problem-size independent.
 *     

 * 
 * The deal.II users who have already gone through the step-20 and step-22
 * tutorials can certainly imagine how we're going to implement this.  We
 * replace the exact inverse matrices in $P^{-1}$ by some approximate
 * inverses built from the InverseMatrix class, and the inverse Schur
 * complement will be approximated by the pressure mass matrix $M_p$
 * (weighted by $\eta^{-1}$ as mentioned in the introduction). As pointed
 * out in the results section of step-22, we can replace the exact inverse
 * of $A$ by just the application of a preconditioner, in this case
 * on a vector Laplace matrix as was explained in the introduction. This
 * does increase the number of (outer) GMRES iterations, but is still
 * significantly cheaper than an exact inverse, which would require
 * between 20 and 35 CG iterations for <em>each</em> outer solver step
 * (using the AMG preconditioner).
 *     

 * 
 * Having the above explanations in mind, we define a preconditioner class
 * with a <code>vmult</code> functionality, which is all we need for the
 * interaction with the usual solver functions further below in the
 * program code.
 *     

 * 
 * First the declarations. These are similar to the definition of the
 * Schur complement in step-20, with the difference that we need some more
 * preconditioners in the constructor and that the matrices we use here
 * are built upon Trilinos:
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     class BlockSchurPreconditioner : public Subscriptor
 *     {
 *     public:
 *       BlockSchurPreconditioner(
 *         const TrilinosWrappers::BlockSparseMatrix &S,
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                             PreconditionerTypeMp> &Mpinv,
 *         const PreconditionerTypeA &                Apreconditioner);
 * 
 *       void vmult(TrilinosWrappers::MPI::BlockVector &      dst,
 *                  const TrilinosWrappers::MPI::BlockVector &src) const;
 * 
 *     private:
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix>
 *         stokes_matrix;
 *       const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                                              PreconditionerTypeMp>>
 *                                  m_inverse;
 *       const PreconditionerTypeA &a_preconditioner;
 * 
 *       mutable TrilinosWrappers::MPI::Vector tmp;
 *     };
 * 
 * 
 * 
 * @endcode
 * 
 * When using a TrilinosWrappers::MPI::Vector or a
 * TrilinosWrappers::MPI::BlockVector, the Vector is initialized using an
 * IndexSet. IndexSet is used not only to resize the
 * TrilinosWrappers::MPI::Vector but it also associates an index in the
 * TrilinosWrappers::MPI::Vector with a degree of freedom (see step-40 for
 * a more detailed explanation). The function complete_index_set() creates
 * an IndexSet where every valid index is part of the set. Note that this
 * program can only be run sequentially and will throw an exception if used
 * in parallel.
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::
 *       BlockSchurPreconditioner(
 *         const TrilinosWrappers::BlockSparseMatrix &S,
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                             PreconditionerTypeMp> &Mpinv,
 *         const PreconditionerTypeA &                Apreconditioner)
 *       : stokes_matrix(&S)
 *       , m_inverse(&Mpinv)
 *       , a_preconditioner(Apreconditioner)
 *       , tmp(complete_index_set(stokes_matrix->block(1, 1).m()))
 *     {}
 * 
 * 
 * @endcode
 * 
 * Next is the <code>vmult</code> function. We implement the action of
 * $P^{-1}$ as described above in three successive steps.  In formulas, we
 * want to compute $Y=P^{-1}X$ where $X,Y$ are both vectors with two block
 * components.
 *     

 * 
 * The first step multiplies the velocity part of the vector by a
 * preconditioner of the matrix $A$, i.e., we compute $Y_0={\tilde
 * A}^{-1}X_0$.  The resulting velocity vector is then multiplied by $B$
 * and subtracted from the pressure, i.e., we want to compute $X_1-BY_0$.
 * This second step only acts on the pressure vector and is accomplished
 * by the residual function of our matrix classes, except that the sign is
 * wrong. Consequently, we change the sign in the temporary pressure
 * vector and finally multiply by the inverse pressure mass matrix to get
 * the final pressure vector, completing our work on the Stokes
 * preconditioner:
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     void
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult(
 *       TrilinosWrappers::MPI::BlockVector &      dst,
 *       const TrilinosWrappers::MPI::BlockVector &src) const
 *     {
 *       a_preconditioner.vmult(dst.block(0), src.block(0));
 *       stokes_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1));
 *       tmp *= -1;
 *       m_inverse->vmult(dst.block(1), tmp);
 *     }
 *   } // namespace LinearSolvers
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBoussinesqFlowProblemcodeclasstemplate"></a> 
 * <h3>The <code>BoussinesqFlowProblem</code> class template</h3>
 * 

 * 
 * The definition of the class that defines the top-level logic of solving
 * the time-dependent Boussinesq problem is mainly based on the step-22
 * tutorial program. The main differences are that now we also have to solve
 * for the temperature equation, which forces us to have a second DoFHandler
 * object for the temperature variable as well as matrices, right hand
 * sides, and solution vectors for the current and previous time steps. As
 * mentioned in the introduction, all linear algebra objects are going to
 * use wrappers of the corresponding Trilinos functionality.
 *   

 * 
 * The member functions of this class are reminiscent of step-21, where we
 * also used a staggered scheme that first solve the flow equations (here
 * the Stokes equations, in step-21 Darcy flow) and then update the advected
 * quantity (here the temperature, there the saturation). The functions that
 * are new are mainly concerned with determining the time step, as well as
 * the proper size of the artificial viscosity stabilization.
 *   

 * 
 * The last three variables indicate whether the various matrices or
 * preconditioners need to be rebuilt the next time the corresponding build
 * functions are called. This allows us to move the corresponding
 * <code>if</code> into the respective function and thereby keeping our main
 * <code>run()</code> function clean and easy to read.
 * 
 * @code
 *   template <int dim>
 *   class BoussinesqFlowProblem
 *   {
 *   public:
 *     BoussinesqFlowProblem();
 *     void run();
 * 
 *   private:
 *     void   setup_dofs();
 *     void   assemble_stokes_preconditioner();
 *     void   build_stokes_preconditioner();
 *     void   assemble_stokes_system();
 *     void   assemble_temperature_system(const double maximal_velocity);
 *     void   assemble_temperature_matrix();
 *     double get_maximal_velocity() const;
 *     std::pair<double, double> get_extrapolated_temperature_range() const;
 *     void                      solve();
 *     void                      output_results() const;
 *     void                      refine_mesh(const unsigned int max_grid_level);
 * 
 *     double compute_viscosity(
 *       const std::vector<double> &        old_temperature,
 *       const std::vector<double> &        old_old_temperature,
 *       const std::vector<Tensor<1, dim>> &old_temperature_grads,
 *       const std::vector<Tensor<1, dim>> &old_old_temperature_grads,
 *       const std::vector<double> &        old_temperature_laplacians,
 *       const std::vector<double> &        old_old_temperature_laplacians,
 *       const std::vector<Tensor<1, dim>> &old_velocity_values,
 *       const std::vector<Tensor<1, dim>> &old_old_velocity_values,
 *       const std::vector<double> &        gamma_values,
 *       const double                       global_u_infty,
 *       const double                       global_T_variation,
 *       const double                       cell_diameter) const;
 * 
 * 
 *     Triangulation<dim> triangulation;
 *     double             global_Omega_diameter;
 * 
 *     const unsigned int        stokes_degree;
 *     FESystem<dim>             stokes_fe;
 *     DoFHandler<dim>           stokes_dof_handler;
 *     AffineConstraints<double> stokes_constraints;
 * 
 *     std::vector<IndexSet>               stokes_partitioning;
 *     TrilinosWrappers::BlockSparseMatrix stokes_matrix;
 *     TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;
 * 
 *     TrilinosWrappers::MPI::BlockVector stokes_solution;
 *     TrilinosWrappers::MPI::BlockVector old_stokes_solution;
 *     TrilinosWrappers::MPI::BlockVector stokes_rhs;
 * 
 * 
 *     const unsigned int        temperature_degree;
 *     FE_Q<dim>                 temperature_fe;
 *     DoFHandler<dim>           temperature_dof_handler;
 *     AffineConstraints<double> temperature_constraints;
 * 
 *     TrilinosWrappers::SparseMatrix temperature_mass_matrix;
 *     TrilinosWrappers::SparseMatrix temperature_stiffness_matrix;
 *     TrilinosWrappers::SparseMatrix temperature_matrix;
 * 
 *     TrilinosWrappers::MPI::Vector temperature_solution;
 *     TrilinosWrappers::MPI::Vector old_temperature_solution;
 *     TrilinosWrappers::MPI::Vector old_old_temperature_solution;
 *     TrilinosWrappers::MPI::Vector temperature_rhs;
 * 
 * 
 *     double       time_step;
 *     double       old_time_step;
 *     unsigned int timestep_number;
 * 
 *     std::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner;
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC>  Mp_preconditioner;
 * 
 *     bool rebuild_stokes_matrix;
 *     bool rebuild_temperature_matrices;
 *     bool rebuild_stokes_preconditioner;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemclassimplementation"></a> 
 * <h3>BoussinesqFlowProblem class implementation</h3>
 * 

 * 
 * 
 * <a name="BoussinesqFlowProblemBoussinesqFlowProblem"></a> 
 * <h4>BoussinesqFlowProblem::BoussinesqFlowProblem</h4>
 *   

 * 
 * The constructor of this class is an extension of the constructor in
 * step-22. We need to add the various variables that concern the
 * temperature. As discussed in the introduction, we are going to use
 * $Q_2\times Q_1$ (Taylor-Hood) elements again for the Stokes part, and
 * $Q_2$ elements for the temperature. However, by using variables that
 * store the polynomial degree of the Stokes and temperature finite
 * elements, it is easy to consistently modify the degree of the elements as
 * well as all quadrature formulas used on them downstream. Moreover, we
 * initialize the time stepping as well as the options for matrix assembly
 * and preconditioning:
 * 
 * @code
 *   template <int dim>
 *   BoussinesqFlowProblem<dim>::BoussinesqFlowProblem()
 *     : triangulation(Triangulation<dim>::maximum_smoothing)
 *     , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN())
 *     , stokes_degree(1)
 *     , stokes_fe(FE_Q<dim>(stokes_degree + 1), dim, FE_Q<dim>(stokes_degree), 1)
 *     , stokes_dof_handler(triangulation)
 *     ,
 * 
 *     temperature_degree(2)
 *     , temperature_fe(temperature_degree)
 *     , temperature_dof_handler(triangulation)
 *     ,
 * 
 *     time_step(0)
 *     , old_time_step(0)
 *     , timestep_number(0)
 *     , rebuild_stokes_matrix(true)
 *     , rebuild_temperature_matrices(true)
 *     , rebuild_stokes_preconditioner(true)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_maximal_velocity"></a> 
 * <h4>BoussinesqFlowProblem::get_maximal_velocity</h4>
 * 

 * 
 * Starting the real functionality of this class is a helper function that
 * determines the maximum ($L_\infty$) velocity in the domain (at the
 * quadrature points, in fact). How it works should be relatively obvious to
 * all who have gotten to this point of the tutorial. Note that since we are
 * only interested in the velocity, rather than using
 * <code>stokes_fe_values.get_function_values</code> to get the values of
 * the entire Stokes solution (velocities and pressures) we use
 * <code>stokes_fe_values[velocities].get_function_values</code> to extract
 * only the velocities part. This has the additional benefit that we get it
 * as a Tensor<1,dim>, rather than some components in a Vector<double>,
 * allowing us to process it right away using the <code>norm()</code>
 * function to get the magnitude of the velocity.
 *   

 * 
 * The only point worth thinking about a bit is how to choose the quadrature
 * points we use here. Since the goal of this function is to find the
 * maximal velocity over a domain by looking at quadrature points on each
 * cell. So we should ask how we should best choose these quadrature points
 * on each cell. To this end, recall that if we had a single $Q_1$ field
 * (rather than the vector-valued field of higher order) then the maximum
 * would be attained at a vertex of the mesh. In other words, we should use
 * the QTrapezoid class that has quadrature points only at the vertices of
 * cells.
 *   

 * 
 * For higher order shape functions, the situation is more complicated: the
 * maxima and minima may be attained at points between the support points of
 * shape functions (for the usual $Q_p$ elements the support points are the
 * equidistant Lagrange interpolation points); furthermore, since we are
 * looking for the maximum magnitude of a vector-valued quantity, we can
 * even less say with certainty where the set of potential maximal points
 * are. Nevertheless, intuitively if not provably, the Lagrange
 * interpolation points appear to be a better choice than the Gauss points.
 *   

 * 
 * There are now different methods to produce a quadrature formula with
 * quadrature points equal to the interpolation points of the finite
 * element. One option would be to use the
 * FiniteElement::get_unit_support_points() function, reduce the output to a
 * unique set of points to avoid duplicate function evaluations, and create
 * a Quadrature object using these points. Another option, chosen here, is
 * to use the QTrapezoid class and combine it with the QIterated class that
 * repeats the QTrapezoid formula on a number of sub-cells in each coordinate
 * direction. To cover all support points, we need to iterate it
 * <code>stokes_degree+1</code> times since this is the polynomial degree of
 * the Stokes element in use:
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::get_maximal_velocity() const
 *   {
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(), stokes_degree + 1);
 *     const unsigned int   n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(stokes_fe, quadrature_formula, update_values);
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points);
 *     double                      max_velocity = 0;
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 *         fe_values[velocities].get_function_values(stokes_solution,
 *                                                   velocity_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           max_velocity = std::max(max_velocity, velocity_values[q].norm());
 *       }
 * 
 *     return max_velocity;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_extrapolated_temperature_range"></a> 
 * <h4>BoussinesqFlowProblem::get_extrapolated_temperature_range</h4>
 * 

 * 
 * Next a function that determines the minimum and maximum temperature at
 * quadrature points inside $\Omega$ when extrapolated from the two previous
 * time steps to the current one. We need this information in the
 * computation of the artificial viscosity parameter $\nu$ as discussed in
 * the introduction.
 *   

 * 
 * The formula for the extrapolated temperature is
 * $\left(1+\frac{k_n}{k_{n-1}} \right)T^{n-1} + \frac{k_n}{k_{n-1}}
 * T^{n-2}$. The way to compute it is to loop over all quadrature points and
 * update the maximum and minimum value if the current value is
 * bigger/smaller than the previous one. We initialize the variables that
 * store the max and min before the loop over all quadrature points by the
 * smallest and the largest number representable as a double. Then we know
 * for a fact that it is larger/smaller than the minimum/maximum and that
 * the loop over all quadrature points is ultimately going to update the
 * initial value with the correct one.
 *   

 * 
 * The only other complication worth mentioning here is that in the first
 * time step, $T^{k-2}$ is not yet available of course. In that case, we can
 * only use $T^{k-1}$ which we have from the initial temperature. As
 * quadrature points, we use the same choice as in the previous function
 * though with the difference that now the number of repetitions is
 * determined by the polynomial degree of the temperature field.
 * 
 * @code
 *   template <int dim>
 *   std::pair<double, double>
 *   BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range() const
 *   {
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(),
 *                                             temperature_degree);
 *     const unsigned int   n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim> fe_values(temperature_fe, quadrature_formula, update_values);
 *     std::vector<double> old_temperature_values(n_q_points);
 *     std::vector<double> old_old_temperature_values(n_q_points);
 * 
 *     if (timestep_number != 0)
 *       {
 *         double min_temperature = std::numeric_limits<double>::max(),
 *                max_temperature = -std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *           {
 *             fe_values.reinit(cell);
 *             fe_values.get_function_values(old_temperature_solution,
 *                                           old_temperature_values);
 *             fe_values.get_function_values(old_old_temperature_solution,
 *                                           old_old_temperature_values);
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q)
 *               {
 *                 const double temperature =
 *                   (1. + time_step / old_time_step) * old_temperature_values[q] -
 *                   time_step / old_time_step * old_old_temperature_values[q];
 * 
 *                 min_temperature = std::min(min_temperature, temperature);
 *                 max_temperature = std::max(max_temperature, temperature);
 *               }
 *           }
 * 
 *         return std::make_pair(min_temperature, max_temperature);
 *       }
 *     else
 *       {
 *         double min_temperature = std::numeric_limits<double>::max(),
 *                max_temperature = -std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *           {
 *             fe_values.reinit(cell);
 *             fe_values.get_function_values(old_temperature_solution,
 *                                           old_temperature_values);
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q)
 *               {
 *                 const double temperature = old_temperature_values[q];
 * 
 *                 min_temperature = std::min(min_temperature, temperature);
 *                 max_temperature = std::max(max_temperature, temperature);
 *               }
 *           }
 * 
 *         return std::make_pair(min_temperature, max_temperature);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemcompute_viscosity"></a> 
 * <h4>BoussinesqFlowProblem::compute_viscosity</h4>
 * 

 * 
 * The last of the tool functions computes the artificial viscosity
 * parameter $\nu|_K$ on a cell $K$ as a function of the extrapolated
 * temperature, its gradient and Hessian (second derivatives), the velocity,
 * the right hand side $\gamma$ all on the quadrature points of the current
 * cell, and various other parameters as described in detail in the
 * introduction.
 *   

 * 
 * There are some universal constants worth mentioning here. First, we need
 * to fix $\beta$; we choose $\beta=0.017\cdot dim$, a choice discussed in
 * detail in the results section of this tutorial program. The second is the
 * exponent $\alpha$; $\alpha=1$ appears to work fine for the current
 * program, even though some additional benefit might be expected from
 * choosing $\alpha = 2$. Finally, there is one thing that requires special
 * casing: In the first time step, the velocity equals zero, and the formula
 * for $\nu|_K$ is not defined. In that case, we return $\nu|_K=5\cdot 10^3
 * \cdot h_K$, a choice admittedly more motivated by heuristics than
 * anything else (it is in the same order of magnitude, however, as the
 * value returned for most cells on the second time step).
 *   

 * 
 * The rest of the function should be mostly obvious based on the material
 * discussed in the introduction:
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::compute_viscosity(
 *     const std::vector<double> &        old_temperature,
 *     const std::vector<double> &        old_old_temperature,
 *     const std::vector<Tensor<1, dim>> &old_temperature_grads,
 *     const std::vector<Tensor<1, dim>> &old_old_temperature_grads,
 *     const std::vector<double> &        old_temperature_laplacians,
 *     const std::vector<double> &        old_old_temperature_laplacians,
 *     const std::vector<Tensor<1, dim>> &old_velocity_values,
 *     const std::vector<Tensor<1, dim>> &old_old_velocity_values,
 *     const std::vector<double> &        gamma_values,
 *     const double                       global_u_infty,
 *     const double                       global_T_variation,
 *     const double                       cell_diameter) const
 *   {
 *     constexpr double beta  = 0.017 * dim;
 *     constexpr double alpha = 1.0;
 * 
 *     if (global_u_infty == 0)
 *       return 5e-3 * cell_diameter;
 * 
 *     const unsigned int n_q_points = old_temperature.size();
 * 
 *     double max_residual = 0;
 *     double max_velocity = 0;
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         const Tensor<1, dim> u =
 *           (old_velocity_values[q] + old_old_velocity_values[q]) / 2;
 * 
 *         const double dT_dt =
 *           (old_temperature[q] - old_old_temperature[q]) / old_time_step;
 *         const double u_grad_T =
 *           u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2;
 * 
 *         const double kappa_Delta_T =
 *           EquationData::kappa *
 *           (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) /
 *           2;
 * 
 *         const double residual =
 *           std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) *
 *                    std::pow((old_temperature[q] + old_old_temperature[q]) / 2,
 *                             alpha - 1.));
 * 
 *         max_residual = std::max(residual, max_residual);
 *         max_velocity = std::max(std::sqrt(u * u), max_velocity);
 *       }
 * 
 *     const double c_R            = std::pow(2., (4. - 2 * alpha) / dim);
 *     const double global_scaling = c_R * global_u_infty * global_T_variation *
 *                                   std::pow(global_Omega_diameter, alpha - 2.);
 * 
 *     return (
 *       beta * max_velocity *
 *       std::min(cell_diameter,
 *                std::pow(cell_diameter, alpha) * max_residual / global_scaling));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemsetup_dofs"></a> 
 * <h4>BoussinesqFlowProblem::setup_dofs</h4>
 *   

 * 
 * This is the function that sets up the DoFHandler objects we have here
 * (one for the Stokes part and one for the temperature part) as well as set
 * to the right sizes the various objects required for the linear algebra in
 * this program. Its basic operations are similar to what we do in step-22.
 *   

 * 
 * The body of the function first enumerates all degrees of freedom for the
 * Stokes and temperature systems. For the Stokes part, degrees of freedom
 * are then sorted to ensure that velocities precede pressure DoFs so that
 * we can partition the Stokes matrix into a $2\times 2$ matrix. As a
 * difference to step-22, we do not perform any additional DoF
 * renumbering. In that program, it paid off since our solver was heavily
 * dependent on ILU's, whereas we use AMG here which is not sensitive to the
 * DoF numbering. The IC preconditioner for the inversion of the pressure
 * mass matrix would of course take advantage of a Cuthill-McKee like
 * renumbering, but its costs are low compared to the velocity portion, so
 * the additional work does not pay off.
 *   

 * 
 * We then proceed with the generation of the hanging node constraints that
 * arise from adaptive grid refinement for both DoFHandler objects. For the
 * velocity, we impose no-flux boundary conditions $\mathbf{u}\cdot
 * \mathbf{n}=0$ by adding constraints to the object that already stores the
 * hanging node constraints matrix. The second parameter in the function
 * describes the first of the velocity components in the total dof vector,
 * which is zero here. The variable <code>no_normal_flux_boundaries</code>
 * denotes the boundary indicators for which to set the no flux boundary
 * conditions; here, this is boundary indicator zero.
 *   

 * 
 * After having done so, we count the number of degrees of freedom in the
 * various blocks:
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::setup_dofs()
 *   {
 *     std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
 *     stokes_sub_blocks[dim] = 1;
 * 
 *     {
 *       stokes_dof_handler.distribute_dofs(stokes_fe);
 *       DoFRenumbering::component_wise(stokes_dof_handler, stokes_sub_blocks);
 * 
 *       stokes_constraints.clear();
 *       DoFTools::make_hanging_node_constraints(stokes_dof_handler,
 *                                               stokes_constraints);
 *       std::set<types::boundary_id> no_normal_flux_boundaries;
 *       no_normal_flux_boundaries.insert(0);
 *       VectorTools::compute_no_normal_flux_constraints(stokes_dof_handler,
 *                                                       0,
 *                                                       no_normal_flux_boundaries,
 *                                                       stokes_constraints);
 *       stokes_constraints.close();
 *     }
 *     {
 *       temperature_dof_handler.distribute_dofs(temperature_fe);
 * 
 *       temperature_constraints.clear();
 *       DoFTools::make_hanging_node_constraints(temperature_dof_handler,
 *                                               temperature_constraints);
 *       temperature_constraints.close();
 *     }
 * 
 *     const std::vector<types::global_dof_index> stokes_dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(stokes_dof_handler, stokes_sub_blocks);
 * 
 *     const unsigned int n_u = stokes_dofs_per_block[0],
 *                        n_p = stokes_dofs_per_block[1],
 *                        n_T = temperature_dof_handler.n_dofs();
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *               << " (on " << triangulation.n_levels() << " levels)" << std::endl
 *               << "Number of degrees of freedom: " << n_u + n_p + n_T << " ("
 *               << n_u << '+' << n_p << '+' << n_T << ')' << std::endl
 *               << std::endl;
 * 
 * @endcode
 * 
 * The next step is to create the sparsity pattern for the Stokes and
 * temperature system matrices as well as the preconditioner matrix from
 * which we build the Stokes preconditioner. As in step-22, we choose to
 * create the pattern by
 * using the blocked version of DynamicSparsityPattern.
 *     

 * 
 * So, we first release the memory stored in the matrices, then set up an
 * object of type BlockDynamicSparsityPattern consisting of
 * $2\times 2$ blocks (for the Stokes system matrix and preconditioner) or
 * DynamicSparsityPattern (for the temperature part). We then
 * fill these objects with the nonzero pattern, taking into account that
 * for the Stokes system matrix, there are no entries in the
 * pressure-pressure block (but all velocity vector components couple with
 * each other and with the pressure). Similarly, in the Stokes
 * preconditioner matrix, only the diagonal blocks are nonzero, since we
 * use the vector Laplacian as discussed in the introduction. This
 * operator only couples each vector component of the Laplacian with
 * itself, but not with the other vector components. (Application of the
 * constraints resulting from the no-flux boundary conditions will couple
 * vector components at the boundary again, however.)
 *     

 * 
 * When generating the sparsity pattern, we directly apply the constraints
 * from hanging nodes and no-flux boundary conditions. This approach was
 * already used in step-27, but is different from the one in early
 * tutorial programs where we first built the original sparsity pattern
 * and only then added the entries resulting from constraints. The reason
 * for doing so is that later during assembly we are going to distribute
 * the constraints immediately when transferring local to global
 * dofs. Consequently, there will be no data written at positions of
 * constrained degrees of freedom, so we can let the
 * DoFTools::make_sparsity_pattern function omit these entries by setting
 * the last Boolean flag to <code>false</code>. Once the sparsity pattern
 * is ready, we can use it to initialize the Trilinos matrices. Since the
 * Trilinos matrices store the sparsity pattern internally, there is no
 * need to keep the sparsity pattern around after the initialization of
 * the matrix.
 * 
 * @code
 *     stokes_partitioning.resize(2);
 *     stokes_partitioning[0] = complete_index_set(n_u);
 *     stokes_partitioning[1] = complete_index_set(n_p);
 *     {
 *       stokes_matrix.clear();
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2);
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u);
 *       dsp.block(0, 1).reinit(n_u, n_p);
 *       dsp.block(1, 0).reinit(n_p, n_u);
 *       dsp.block(1, 1).reinit(n_p, n_p);
 * 
 *       dsp.collect_sizes();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 * 
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (!((c == dim) && (d == dim)))
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       DoFTools::make_sparsity_pattern(
 *         stokes_dof_handler, coupling, dsp, stokes_constraints, false);
 * 
 *       stokes_matrix.reinit(dsp);
 *     }
 * 
 *     {
 *       Amg_preconditioner.reset();
 *       Mp_preconditioner.reset();
 *       stokes_preconditioner_matrix.clear();
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2);
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u);
 *       dsp.block(0, 1).reinit(n_u, n_p);
 *       dsp.block(1, 0).reinit(n_p, n_u);
 *       dsp.block(1, 1).reinit(n_p, n_p);
 * 
 *       dsp.collect_sizes();
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *       for (unsigned int c = 0; c < dim + 1; ++c)
 *         for (unsigned int d = 0; d < dim + 1; ++d)
 *           if (c == d)
 *             coupling[c][d] = DoFTools::always;
 *           else
 *             coupling[c][d] = DoFTools::none;
 * 
 *       DoFTools::make_sparsity_pattern(
 *         stokes_dof_handler, coupling, dsp, stokes_constraints, false);
 * 
 *       stokes_preconditioner_matrix.reinit(dsp);
 *     }
 * 
 * @endcode
 * 
 * The creation of the temperature matrix (or, rather, matrices, since we
 * provide a temperature mass matrix and a temperature stiffness matrix,
 * that will be added together for time discretization) follows the
 * generation of the Stokes matrix &ndash; except that it is much easier
 * here since we do not need to take care of any blocks or coupling
 * between components. Note how we initialize the three temperature
 * matrices: We only use the sparsity pattern for reinitialization of the
 * first matrix, whereas we use the previously generated matrix for the
 * two remaining reinits. The reason for doing so is that reinitialization
 * from an already generated matrix allows Trilinos to reuse the sparsity
 * pattern instead of generating a new one for each copy. This saves both
 * some time and memory.
 * 
 * @code
 *     {
 *       temperature_mass_matrix.clear();
 *       temperature_stiffness_matrix.clear();
 *       temperature_matrix.clear();
 * 
 *       DynamicSparsityPattern dsp(n_T, n_T);
 *       DoFTools::make_sparsity_pattern(temperature_dof_handler,
 *                                       dsp,
 *                                       temperature_constraints,
 *                                       false);
 * 
 *       temperature_matrix.reinit(dsp);
 *       temperature_mass_matrix.reinit(temperature_matrix);
 *       temperature_stiffness_matrix.reinit(temperature_matrix);
 *     }
 * 
 * @endcode
 * 
 * Lastly, we set the vectors for the Stokes solutions $\mathbf u^{n-1}$
 * and $\mathbf u^{n-2}$, as well as for the temperatures $T^{n}$,
 * $T^{n-1}$ and $T^{n-2}$ (required for time stepping) and all the system
 * right hand sides to their correct sizes and block structure:
 * 
 * @code
 *     IndexSet temperature_partitioning = complete_index_set(n_T);
 *     stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD);
 *     old_stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD);
 *     stokes_rhs.reinit(stokes_partitioning, MPI_COMM_WORLD);
 * 
 *     temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD);
 *     old_temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD);
 *     old_old_temperature_solution.reinit(temperature_partitioning,
 *                                         MPI_COMM_WORLD);
 * 
 *     temperature_rhs.reinit(temperature_partitioning, MPI_COMM_WORLD);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_stokes_preconditioner"></a> 
 * <h4>BoussinesqFlowProblem::assemble_stokes_preconditioner</h4>
 *   

 * 
 * This function assembles the matrix we use for preconditioning the Stokes
 * system. What we need are a vector Laplace matrix on the velocity
 * components and a mass matrix weighted by $\eta^{-1}$ on the pressure
 * component. We start by generating a quadrature object of appropriate
 * order, the FEValues object that can give values and gradients at the
 * quadrature points (together with quadrature weights). Next we create data
 * structures for the cell matrix and the relation between local and global
 * DoFs. The vectors <code>grad_phi_u</code> and <code>phi_p</code> are
 * going to hold the values of the basis functions in order to faster build
 * up the local matrices, as was already done in step-22. Before we start
 * the loop over all active cells, we have to specify which components are
 * pressure and which are velocity.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner()
 *   {
 *     stokes_preconditioner_matrix = 0;
 * 
 *     const QGauss<dim> quadrature_formula(stokes_degree + 2);
 *     FEValues<dim>     stokes_fe_values(stokes_fe,
 *                                    quadrature_formula,
 *                                    update_JxW_values | update_values |
 *                                      update_gradients);
 * 
 *     const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
 *     std::vector<double>         phi_p(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators())
 *       {
 *         stokes_fe_values.reinit(cell);
 *         local_matrix = 0;
 * 
 * @endcode
 * 
 * The creation of the local matrix is rather simple. There are only a
 * Laplace term (on the velocity) and a mass matrix weighted by
 * $\eta^{-1}$ to be generated, so the creation of the local matrix is
 * done in two lines. Once the local matrix is ready (loop over rows
 * and columns in the local matrix on each quadrature point), we get
 * the local DoF indices and write the local information into the
 * global matrix. We do this as in step-27, i.e., we directly apply the
 * constraints from hanging nodes locally. By doing so, we don't have
 * to do that afterwards, and we don't also write into entries of the
 * matrix that will actually be set to zero again later when
 * eliminating constraints.
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 grad_phi_u[k] = stokes_fe_values[velocities].gradient(k, q);
 *                 phi_p[k]      = stokes_fe_values[pressure].value(k, q);
 *               }
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 local_matrix(i, j) +=
 *                   (EquationData::eta *
 *                      scalar_product(grad_phi_u[i], grad_phi_u[j]) +
 *                    (1. / EquationData::eta) * phi_p[i] * phi_p[j]) *
 *                   stokes_fe_values.JxW(q);
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         stokes_constraints.distribute_local_to_global(
 *           local_matrix, local_dof_indices, stokes_preconditioner_matrix);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblembuild_stokes_preconditioner"></a> 
 * <h4>BoussinesqFlowProblem::build_stokes_preconditioner</h4>
 *   

 * 
 * This function generates the inner preconditioners that are going to be
 * used for the Schur complement block preconditioner. Since the
 * preconditioners need only to be regenerated when the matrices change,
 * this function does not have to do anything in case the matrices have not
 * changed (i.e., the flag <code>rebuild_stokes_preconditioner</code> has
 * the value <code>false</code>). Otherwise its first task is to call
 * <code>assemble_stokes_preconditioner</code> to generate the
 * preconditioner matrices.
 *   

 * 
 * Next, we set up the preconditioner for the velocity-velocity matrix
 * $A$. As explained in the introduction, we are going to use an AMG
 * preconditioner based on a vector Laplace matrix $\hat{A}$ (which is
 * spectrally close to the Stokes matrix $A$). Usually, the
 * TrilinosWrappers::PreconditionAMG class can be seen as a good black-box
 * preconditioner which does not need any special knowledge. In this case,
 * however, we have to be careful: since we build an AMG for a vector
 * problem, we have to tell the preconditioner setup which dofs belong to
 * which vector component. We do this using the function
 * DoFTools::extract_constant_modes, a function that generates a set of
 * <code>dim</code> vectors, where each one has ones in the respective
 * component of the vector problem and zeros elsewhere. Hence, these are the
 * constant modes on each component, which explains the name of the
 * variable.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::build_stokes_preconditioner()
 *   {
 *     if (rebuild_stokes_preconditioner == false)
 *       return;
 * 
 *     std::cout << "   Rebuilding Stokes preconditioner..." << std::flush;
 * 
 *     assemble_stokes_preconditioner();
 * 
 *     Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
 * 
 *     std::vector<std::vector<bool>> constant_modes;
 *     FEValuesExtractors::Vector     velocity_components(0);
 *     DoFTools::extract_constant_modes(stokes_dof_handler,
 *                                      stokes_fe.component_mask(
 *                                        velocity_components),
 *                                      constant_modes);
 *     TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
 *     amg_data.constant_modes = constant_modes;
 * 
 * @endcode
 * 
 * Next, we set some more options of the AMG preconditioner. In
 * particular, we need to tell the AMG setup that we use quadratic basis
 * functions for the velocity matrix (this implies more nonzero elements
 * in the matrix, so that a more robust algorithm needs to be chosen
 * internally). Moreover, we want to be able to control how the coarsening
 * structure is build up. The way the Trilinos smoothed aggregation AMG
 * does this is to look which matrix entries are of similar size as the
 * diagonal entry in order to algebraically build a coarse-grid
 * structure. By setting the parameter <code>aggregation_threshold</code>
 * to 0.02, we specify that all entries that are more than two percent of
 * size of some diagonal pivots in that row should form one coarse grid
 * point. This parameter is rather ad hoc, and some fine-tuning of it can
 * influence the performance of the preconditioner. As a rule of thumb,
 * larger values of <code>aggregation_threshold</code> will decrease the
 * number of iterations, but increase the costs per iteration. A look at
 * the Trilinos documentation will provide more information on these
 * parameters. With this data set, we then initialize the preconditioner
 * with the matrix we want it to apply to.
 *     

 * 
 * Finally, we also initialize the preconditioner for the inversion of the
 * pressure mass matrix. This matrix is symmetric and well-behaved, so we
 * can chose a simple preconditioner. We stick with an incomplete Cholesky
 * (IC) factorization preconditioner, which is designed for symmetric
 * matrices. We could have also chosen an SSOR preconditioner with
 * relaxation factor around 1.2, but IC is cheaper for our example. We
 * wrap the preconditioners into a <code>std::shared_ptr</code>
 * pointer, which makes it easier to recreate the preconditioner next time
 * around since we do not have to care about destroying the previously
 * used object.
 * 
 * @code
 *     amg_data.elliptic              = true;
 *     amg_data.higher_order_elements = true;
 *     amg_data.smoother_sweeps       = 2;
 *     amg_data.aggregation_threshold = 0.02;
 *     Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0, 0),
 *                                    amg_data);
 * 
 *     Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>();
 *     Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1, 1));
 * 
 *     std::cout << std::endl;
 * 
 *     rebuild_stokes_preconditioner = false;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_stokes_system"></a> 
 * <h4>BoussinesqFlowProblem::assemble_stokes_system</h4>
 *   

 * 
 * The time lag scheme we use for advancing the coupled Stokes-temperature
 * system forces us to split up the assembly (and the solution of linear
 * systems) into two step. The first one is to create the Stokes system
 * matrix and right hand side, and the second is to create matrix and right
 * hand sides for the temperature dofs, which depends on the result of the
 * linear system for the velocity.
 *   

 * 
 * This function is called at the beginning of each time step. In the first
 * time step or if the mesh has changed, indicated by the
 * <code>rebuild_stokes_matrix</code>, we need to assemble the Stokes
 * matrix; on the other hand, if the mesh hasn't changed and the matrix is
 * already available, this is not necessary and all we need to do is
 * assemble the right hand side vector which changes in each time step.
 *   

 * 
 * Regarding the technical details of implementation, not much has changed
 * from step-22. We reset matrix and vector, create a quadrature formula on
 * the cells, and then create the respective FEValues object. For the update
 * flags, we require basis function derivatives only in case of a full
 * assembly, since they are not needed for the right hand side; as always,
 * choosing the minimal set of flags depending on what is currently needed
 * makes the call to FEValues::reinit further down in the program more
 * efficient.
 *   

 * 
 * There is one thing that needs to be commented &ndash; since we have a
 * separate finite element and DoFHandler for the temperature, we need to
 * generate a second FEValues object for the proper evaluation of the
 * temperature solution. This isn't too complicated to realize here: just
 * use the temperature structures and set an update flag for the basis
 * function values which we need for evaluation of the temperature
 * solution. The only important part to remember here is that the same
 * quadrature formula is used for both FEValues objects to ensure that we
 * get matching information when we loop over the quadrature points of the
 * two objects.
 *   

 * 
 * The declarations proceed with some shortcuts for array sizes, the
 * creation of the local matrix and right hand side as well as the vector
 * for the indices of the local dofs compared to the global system.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_system()
 *   {
 *     std::cout << "   Assembling..." << std::flush;
 * 
 *     if (rebuild_stokes_matrix == true)
 *       stokes_matrix = 0;
 * 
 *     stokes_rhs = 0;
 * 
 *     const QGauss<dim> quadrature_formula(stokes_degree + 2);
 *     FEValues<dim>     stokes_fe_values(
 *       stokes_fe,
 *       quadrature_formula,
 *       update_values | update_quadrature_points | update_JxW_values |
 *         (rebuild_stokes_matrix == true ? update_gradients : UpdateFlags(0)));
 * 
 *     FEValues<dim> temperature_fe_values(temperature_fe,
 *                                         quadrature_formula,
 *                                         update_values);
 * 
 *     const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Next we need a vector that will contain the values of the temperature
 * solution at the previous time level at the quadrature points to
 * assemble the source term in the right hand side of the momentum
 * equation. Let's call this vector <code>old_solution_values</code>.
 *     

 * 
 * The set of vectors we create next hold the evaluations of the basis
 * functions as well as their gradients and symmetrized gradients that
 * will be used for creating the matrices. Putting these into their own
 * arrays rather than asking the FEValues object for this information each
 * time it is needed is an optimization to accelerate the assembly
 * process, see step-22 for details.
 *     

 * 
 * The last two declarations are used to extract the individual blocks
 * (velocity, pressure, temperature) from the total FE system.
 * 
 * @code
 *     std::vector<double> old_temperature_values(n_q_points);
 * 
 *     std::vector<Tensor<1, dim>>          phi_u(dofs_per_cell);
 *     std::vector<SymmetricTensor<2, dim>> grads_phi_u(dofs_per_cell);
 *     std::vector<double>                  div_phi_u(dofs_per_cell);
 *     std::vector<double>                  phi_p(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 * @endcode
 * 
 * Now start the loop over all cells in the problem. We are working on two
 * different DoFHandlers for this assembly routine, so we must have two
 * different cell iterators for the two objects in use. This might seem a
 * bit peculiar, since both the Stokes system and the temperature system
 * use the same grid, but that's the only way to keep degrees of freedom
 * in sync. The first statements within the loop are again all very
 * familiar, doing the update of the finite element data as specified by
 * the update flags, zeroing out the local arrays and getting the values
 * of the old solution at the quadrature points. Then we are ready to loop
 * over the quadrature points on the cell.
 * 
 * @code
 *     auto       cell             = stokes_dof_handler.begin_active();
 *     const auto endc             = stokes_dof_handler.end();
 *     auto       temperature_cell = temperature_dof_handler.begin_active();
 * 
 *     for (; cell != endc; ++cell, ++temperature_cell)
 *       {
 *         stokes_fe_values.reinit(cell);
 *         temperature_fe_values.reinit(temperature_cell);
 * 
 *         local_matrix = 0;
 *         local_rhs    = 0;
 * 
 *         temperature_fe_values.get_function_values(old_temperature_solution,
 *                                                   old_temperature_values);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             const double old_temperature = old_temperature_values[q];
 * 
 * @endcode
 * 
 * Next we extract the values and gradients of basis functions
 * relevant to the terms in the inner products. As shown in
 * step-22 this helps accelerate assembly.
 *             

 * 
 * Once this is done, we start the loop over the rows and columns
 * of the local matrix and feed the matrix with the relevant
 * products. The right hand side is filled with the forcing term
 * driven by temperature in direction of gravity (which is
 * vertical in our example).  Note that the right hand side term
 * is always generated, whereas the matrix contributions are only
 * updated when it is requested by the
 * <code>rebuild_matrices</code> flag.
 * 
 * @code
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 phi_u[k] = stokes_fe_values[velocities].value(k, q);
 *                 if (rebuild_stokes_matrix)
 *                   {
 *                     grads_phi_u[k] =
 *                       stokes_fe_values[velocities].symmetric_gradient(k, q);
 *                     div_phi_u[k] =
 *                       stokes_fe_values[velocities].divergence(k, q);
 *                     phi_p[k] = stokes_fe_values[pressure].value(k, q);
 *                   }
 *               }
 * 
 *             if (rebuild_stokes_matrix)
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                   local_matrix(i, j) +=
 *                     (EquationData::eta * 2 * (grads_phi_u[i] * grads_phi_u[j]) -
 *                      div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
 *                     stokes_fe_values.JxW(q);
 * 
 *             const Point<dim> gravity =
 *               -((dim == 2) ? (Point<dim>(0, 1)) : (Point<dim>(0, 0, 1)));
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               local_rhs(i) += (-EquationData::density * EquationData::beta *
 *                                gravity * phi_u[i] * old_temperature) *
 *                               stokes_fe_values.JxW(q);
 *           }
 * 
 * @endcode
 * 
 * The last step in the loop over all cells is to enter the local
 * contributions into the global matrix and vector structures to the
 * positions specified in <code>local_dof_indices</code>.  Again, we
 * let the AffineConstraints class do the insertion of the cell
 * matrix elements to the global matrix, which already condenses the
 * hanging node constraints.
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         if (rebuild_stokes_matrix == true)
 *           stokes_constraints.distribute_local_to_global(local_matrix,
 *                                                         local_rhs,
 *                                                         local_dof_indices,
 *                                                         stokes_matrix,
 *                                                         stokes_rhs);
 *         else
 *           stokes_constraints.distribute_local_to_global(local_rhs,
 *                                                         local_dof_indices,
 *                                                         stokes_rhs);
 *       }
 * 
 *     rebuild_stokes_matrix = false;
 * 
 *     std::cout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_temperature_matrix"></a> 
 * <h4>BoussinesqFlowProblem::assemble_temperature_matrix</h4>
 *   

 * 
 * This function assembles the matrix in the temperature equation. The
 * temperature matrix consists of two parts, a mass matrix and the time step
 * size times a stiffness matrix given by a Laplace term times the amount of
 * diffusion. Since the matrix depends on the time step size (which varies
 * from one step to another), the temperature matrix needs to be updated
 * every time step. We could simply regenerate the matrices in every time
 * step, but this is not really efficient since mass and Laplace matrix do
 * only change when we change the mesh. Hence, we do this more efficiently
 * by generating two separate matrices in this function, one for the mass
 * matrix and one for the stiffness (diffusion) matrix. We will then sum up
 * the matrix plus the stiffness matrix times the time step size once we
 * know the actual time step.
 *   

 * 
 * So the details for this first step are very simple. In case we need to
 * rebuild the matrix (i.e., the mesh has changed), we zero the data
 * structures, get a quadrature formula and a FEValues object, and create
 * local matrices, local dof indices and evaluation structures for the basis
 * functions.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_temperature_matrix()
 *   {
 *     if (rebuild_temperature_matrices == false)
 *       return;
 * 
 *     temperature_mass_matrix      = 0;
 *     temperature_stiffness_matrix = 0;
 * 
 *     QGauss<dim>   quadrature_formula(temperature_degree + 2);
 *     FEValues<dim> temperature_fe_values(temperature_fe,
 *                                         quadrature_formula,
 *                                         update_values | update_gradients |
 *                                           update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);
 *     FullMatrix<double> local_stiffness_matrix(dofs_per_cell, dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     std::vector<double>         phi_T(dofs_per_cell);
 *     std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell);
 * 
 * @endcode
 * 
 * Now, let's start the loop over all cells in the triangulation. We need
 * to zero out the local matrices, update the finite element evaluations,
 * and then loop over the rows and columns of the matrices on each
 * quadrature point, where we then create the mass matrix and the
 * stiffness matrix (Laplace terms times the diffusion
 * <code>EquationData::kappa</code>. Finally, we let the constraints
 * object insert these values into the global matrix, and directly
 * condense the constraints into the matrix.
 * 
 * @code
 *     for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *       {
 *         local_mass_matrix      = 0;
 *         local_stiffness_matrix = 0;
 * 
 *         temperature_fe_values.reinit(cell);
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 grad_phi_T[k] = temperature_fe_values.shape_grad(k, q);
 *                 phi_T[k]      = temperature_fe_values.shape_value(k, q);
 *               }
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   local_mass_matrix(i, j) +=
 *                     (phi_T[i] * phi_T[j] * temperature_fe_values.JxW(q));
 *                   local_stiffness_matrix(i, j) +=
 *                     (EquationData::kappa * grad_phi_T[i] * grad_phi_T[j] *
 *                      temperature_fe_values.JxW(q));
 *                 }
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         temperature_constraints.distribute_local_to_global(
 *           local_mass_matrix, local_dof_indices, temperature_mass_matrix);
 *         temperature_constraints.distribute_local_to_global(
 *           local_stiffness_matrix,
 *           local_dof_indices,
 *           temperature_stiffness_matrix);
 *       }
 * 
 *     rebuild_temperature_matrices = false;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_temperature_system"></a> 
 * <h4>BoussinesqFlowProblem::assemble_temperature_system</h4>
 *   

 * 
 * This function does the second part of the assembly work on the
 * temperature matrix, the actual addition of pressure mass and stiffness
 * matrix (where the time step size comes into play), as well as the
 * creation of the velocity-dependent right hand side. The declarations for
 * the right hand side assembly in this function are pretty much the same as
 * the ones used in the other assembly routines, except that we restrict
 * ourselves to vectors this time. We are going to calculate residuals on
 * the temperature system, which means that we have to evaluate second
 * derivatives, specified by the update flag <code>update_hessians</code>.
 *   

 * 
 * The temperature equation is coupled to the Stokes system by means of the
 * fluid velocity. These two parts of the solution are associated with
 * different DoFHandlers, so we again need to create a second FEValues
 * object for the evaluation of the velocity at the quadrature points.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_temperature_system(
 *     const double maximal_velocity)
 *   {
 *     const bool use_bdf2_scheme = (timestep_number != 0);
 * 
 *     if (use_bdf2_scheme == true)
 *       {
 *         temperature_matrix.copy_from(temperature_mass_matrix);
 *         temperature_matrix *=
 *           (2 * time_step + old_time_step) / (time_step + old_time_step);
 *         temperature_matrix.add(time_step, temperature_stiffness_matrix);
 *       }
 *     else
 *       {
 *         temperature_matrix.copy_from(temperature_mass_matrix);
 *         temperature_matrix.add(time_step, temperature_stiffness_matrix);
 *       }
 * 
 *     temperature_rhs = 0;
 * 
 *     const QGauss<dim> quadrature_formula(temperature_degree + 2);
 *     FEValues<dim>     temperature_fe_values(temperature_fe,
 *                                         quadrature_formula,
 *                                         update_values | update_gradients |
 *                                           update_hessians |
 *                                           update_quadrature_points |
 *                                           update_JxW_values);
 *     FEValues<dim>     stokes_fe_values(stokes_fe,
 *                                    quadrature_formula,
 *                                    update_values);
 * 
 *     const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     Vector<double> local_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Next comes the declaration of vectors to hold the old and older
 * solution values (as a notation for time levels $n-1$ and
 * $n-2$, respectively) and gradients at quadrature points of the
 * current cell. We also declare an object to hold the temperature right
 * hand side values (<code>gamma_values</code>), and we again use
 * shortcuts for the temperature basis functions. Eventually, we need to
 * find the temperature extrema and the diameter of the computational
 * domain which will be used for the definition of the stabilization
 * parameter (we got the maximal velocity as an input to this function).
 * 
 * @code
 *     std::vector<Tensor<1, dim>> old_velocity_values(n_q_points);
 *     std::vector<Tensor<1, dim>> old_old_velocity_values(n_q_points);
 *     std::vector<double>         old_temperature_values(n_q_points);
 *     std::vector<double>         old_old_temperature_values(n_q_points);
 *     std::vector<Tensor<1, dim>> old_temperature_grads(n_q_points);
 *     std::vector<Tensor<1, dim>> old_old_temperature_grads(n_q_points);
 *     std::vector<double>         old_temperature_laplacians(n_q_points);
 *     std::vector<double>         old_old_temperature_laplacians(n_q_points);
 * 
 *     EquationData::TemperatureRightHandSide<dim> temperature_right_hand_side;
 *     std::vector<double>                         gamma_values(n_q_points);
 * 
 *     std::vector<double>         phi_T(dofs_per_cell);
 *     std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell);
 * 
 *     const std::pair<double, double> global_T_range =
 *       get_extrapolated_temperature_range();
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 * @endcode
 * 
 * Now, let's start the loop over all cells in the triangulation. Again,
 * we need two cell iterators that walk in parallel through the cells of
 * the two involved DoFHandler objects for the Stokes and temperature
 * part. Within the loop, we first set the local rhs to zero, and then get
 * the values and derivatives of the old solution functions at the
 * quadrature points, since they are going to be needed for the definition
 * of the stabilization parameters and as coefficients in the equation,
 * respectively. Note that since the temperature has its own DoFHandler
 * and FEValues object we get the entire solution at the quadrature point
 * (which is the scalar temperature field only anyway) whereas for the
 * Stokes part we restrict ourselves to extracting the velocity part (and
 * ignoring the pressure part) by using
 * <code>stokes_fe_values[velocities].get_function_values</code>.
 * 
 * @code
 *     auto       cell        = temperature_dof_handler.begin_active();
 *     const auto endc        = temperature_dof_handler.end();
 *     auto       stokes_cell = stokes_dof_handler.begin_active();
 * 
 *     for (; cell != endc; ++cell, ++stokes_cell)
 *       {
 *         local_rhs = 0;
 * 
 *         temperature_fe_values.reinit(cell);
 *         stokes_fe_values.reinit(stokes_cell);
 * 
 *         temperature_fe_values.get_function_values(old_temperature_solution,
 *                                                   old_temperature_values);
 *         temperature_fe_values.get_function_values(old_old_temperature_solution,
 *                                                   old_old_temperature_values);
 * 
 *         temperature_fe_values.get_function_gradients(old_temperature_solution,
 *                                                      old_temperature_grads);
 *         temperature_fe_values.get_function_gradients(
 *           old_old_temperature_solution, old_old_temperature_grads);
 * 
 *         temperature_fe_values.get_function_laplacians(
 *           old_temperature_solution, old_temperature_laplacians);
 *         temperature_fe_values.get_function_laplacians(
 *           old_old_temperature_solution, old_old_temperature_laplacians);
 * 
 *         temperature_right_hand_side.value_list(
 *           temperature_fe_values.get_quadrature_points(), gamma_values);
 * 
 *         stokes_fe_values[velocities].get_function_values(stokes_solution,
 *                                                          old_velocity_values);
 *         stokes_fe_values[velocities].get_function_values(
 *           old_stokes_solution, old_old_velocity_values);
 * 
 * @endcode
 * 
 * Next, we calculate the artificial viscosity for stabilization
 * according to the discussion in the introduction using the dedicated
 * function. With that at hand, we can get into the loop over
 * quadrature points and local rhs vector components. The terms here
 * are quite lengthy, but their definition follows the time-discrete
 * system developed in the introduction of this program. The BDF-2
 * scheme needs one more term from the old time step (and involves
 * more complicated factors) than the backward Euler scheme that is
 * used for the first time step. When all this is done, we distribute
 * the local vector into the global one (including hanging node
 * constraints).
 * 
 * @code
 *         const double nu =
 *           compute_viscosity(old_temperature_values,
 *                             old_old_temperature_values,
 *                             old_temperature_grads,
 *                             old_old_temperature_grads,
 *                             old_temperature_laplacians,
 *                             old_old_temperature_laplacians,
 *                             old_velocity_values,
 *                             old_old_velocity_values,
 *                             gamma_values,
 *                             maximal_velocity,
 *                             global_T_range.second - global_T_range.first,
 *                             cell->diameter());
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *               {
 *                 grad_phi_T[k] = temperature_fe_values.shape_grad(k, q);
 *                 phi_T[k]      = temperature_fe_values.shape_value(k, q);
 *               }
 * 
 *             const double T_term_for_rhs =
 *               (use_bdf2_scheme ?
 *                  (old_temperature_values[q] * (1 + time_step / old_time_step) -
 *                   old_old_temperature_values[q] * (time_step * time_step) /
 *                     (old_time_step * (time_step + old_time_step))) :
 *                  old_temperature_values[q]);
 * 
 *             const Tensor<1, dim> ext_grad_T =
 *               (use_bdf2_scheme ?
 *                  (old_temperature_grads[q] * (1 + time_step / old_time_step) -
 *                   old_old_temperature_grads[q] * time_step / old_time_step) :
 *                  old_temperature_grads[q]);
 * 
 *             const Tensor<1, dim> extrapolated_u =
 *               (use_bdf2_scheme ?
 *                  (old_velocity_values[q] * (1 + time_step / old_time_step) -
 *                   old_old_velocity_values[q] * time_step / old_time_step) :
 *                  old_velocity_values[q]);
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               local_rhs(i) +=
 *                 (T_term_for_rhs * phi_T[i] -
 *                  time_step * extrapolated_u * ext_grad_T * phi_T[i] -
 *                  time_step * nu * ext_grad_T * grad_phi_T[i] +
 *                  time_step * gamma_values[q] * phi_T[i]) *
 *                 temperature_fe_values.JxW(q);
 *           }
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         temperature_constraints.distribute_local_to_global(local_rhs,
 *                                                            local_dof_indices,
 *                                                            temperature_rhs);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemsolve"></a> 
 * <h4>BoussinesqFlowProblem::solve</h4>
 *   

 * 
 * This function solves the linear systems of equations. Following the
 * introduction, we start with the Stokes system, where we need to generate
 * our block Schur preconditioner. Since all the relevant actions are
 * implemented in the class <code>BlockSchurPreconditioner</code>, all we
 * have to do is to initialize the class appropriately. What we need to pass
 * down is an <code>InverseMatrix</code> object for the pressure mass
 * matrix, which we set up using the respective class together with the IC
 * preconditioner we already generated, and the AMG preconditioner for the
 * velocity-velocity matrix. Note that both <code>Mp_preconditioner</code>
 * and <code>Amg_preconditioner</code> are only pointers, so we use
 * <code>*</code> to pass down the actual preconditioner objects.
 *   

 * 
 * Once the preconditioner is ready, we create a GMRES solver for the block
 * system. Since we are working with Trilinos data structures, we have to
 * set the respective template argument in the solver. GMRES needs to
 * internally store temporary vectors for each iteration (see the discussion
 * in the results section of step-22) &ndash; the more vectors it can use,
 * the better it will generally perform. To keep memory demands in check, we
 * set the number of vectors to 100. This means that up to 100 solver
 * iterations, every temporary vector can be stored. If the solver needs to
 * iterate more often to get the specified tolerance, it will work on a
 * reduced set of vectors by restarting at every 100 iterations.
 *   

 * 
 * With this all set up, we solve the system and distribute the constraints
 * in the Stokes system, i.e., hanging nodes and no-flux boundary condition,
 * in order to have the appropriate solution values even at constrained
 * dofs. Finally, we write the number of iterations to the screen.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::solve()
 *   {
 *     std::cout << "   Solving..." << std::endl;
 * 
 *     {
 *       const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
 *                                          TrilinosWrappers::PreconditionIC>
 *         mp_inverse(stokes_preconditioner_matrix.block(1, 1),
 *                    *Mp_preconditioner);
 * 
 *       const LinearSolvers::BlockSchurPreconditioner<
 *         TrilinosWrappers::PreconditionAMG,
 *         TrilinosWrappers::PreconditionIC>
 *         preconditioner(stokes_matrix, mp_inverse, *Amg_preconditioner);
 * 
 *       SolverControl solver_control(stokes_matrix.m(),
 *                                    1e-6 * stokes_rhs.l2_norm());
 * 
 *       SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres(
 *         solver_control,
 *         SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(100));
 * 
 *       for (unsigned int i = 0; i < stokes_solution.size(); ++i)
 *         if (stokes_constraints.is_constrained(i))
 *           stokes_solution(i) = 0;
 * 
 *       gmres.solve(stokes_matrix, stokes_solution, stokes_rhs, preconditioner);
 * 
 *       stokes_constraints.distribute(stokes_solution);
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " GMRES iterations for Stokes subsystem." << std::endl;
 *     }
 * 
 * @endcode
 * 
 * Once we know the Stokes solution, we can determine the new time step
 * from the maximal velocity. We have to do this to satisfy the CFL
 * condition since convection terms are treated explicitly in the
 * temperature equation, as discussed in the introduction. The exact form
 * of the formula used here for the time step is discussed in the results
 * section of this program.
 *     

 * 
 * There is a snatch here. The formula contains a division by the maximum
 * value of the velocity. However, at the start of the computation, we
 * have a constant temperature field (we start with a constant
 * temperature, and it will be nonconstant only after the first time step
 * during which the source acts). Constant temperature means that no
 * buoyancy acts, and so the velocity is zero. Dividing by it will not
 * likely lead to anything good.
 *     

 * 
 * To avoid the resulting infinite time step, we ask whether the maximal
 * velocity is very small (in particular smaller than the values we
 * encounter during any of the following time steps) and if so rather than
 * dividing by zero we just divide by a small value, resulting in a large
 * but finite time step.
 * 
 * @code
 *     old_time_step                 = time_step;
 *     const double maximal_velocity = get_maximal_velocity();
 * 
 *     if (maximal_velocity >= 0.01)
 *       time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree *
 *                   GridTools::minimal_cell_diameter(triangulation) /
 *                   maximal_velocity;
 *     else
 *       time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree *
 *                   GridTools::minimal_cell_diameter(triangulation) / .01;
 * 
 *     std::cout << "   "
 *               << "Time step: " << time_step << std::endl;
 * 
 *     temperature_solution = old_temperature_solution;
 * 
 * @endcode
 * 
 * Next we set up the temperature system and the right hand side using the
 * function <code>assemble_temperature_system()</code>.  Knowing the
 * matrix and right hand side of the temperature equation, we set up a
 * preconditioner and a solver. The temperature matrix is a mass matrix
 * (with eigenvalues around one) plus a Laplace matrix (with eigenvalues
 * between zero and $ch^{-2}$) times a small number proportional to the
 * time step $k_n$. Hence, the resulting symmetric and positive definite
 * matrix has eigenvalues in the range $[1,1+k_nh^{-2}]$ (up to
 * constants). This matrix is only moderately ill conditioned even for
 * small mesh sizes and we get a reasonably good preconditioner by simple
 * means, for example with an incomplete Cholesky decomposition
 * preconditioner (IC) as we also use for preconditioning the pressure
 * mass matrix solver. As a solver, we choose the conjugate gradient
 * method CG. As before, we tell the solver to use Trilinos vectors via
 * the template argument <code>TrilinosWrappers::MPI::Vector</code>.
 * Finally, we solve, distribute the hanging node constraints and write out
 * the number of iterations.
 * 
 * @code
 *     assemble_temperature_system(maximal_velocity);
 *     {
 *       SolverControl solver_control(temperature_matrix.m(),
 *                                    1e-8 * temperature_rhs.l2_norm());
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);
 * 
 *       TrilinosWrappers::PreconditionIC preconditioner;
 *       preconditioner.initialize(temperature_matrix);
 * 
 *       cg.solve(temperature_matrix,
 *                temperature_solution,
 *                temperature_rhs,
 *                preconditioner);
 * 
 *       temperature_constraints.distribute(temperature_solution);
 * 
 *       std::cout << "   " << solver_control.last_step()
 *                 << " CG iterations for temperature." << std::endl;
 * 
 * @endcode
 * 
 * At the end of this function, we step through the vector and read out
 * the maximum and minimum temperature value, which we also want to
 * output. This will come in handy when determining the correct constant
 * in the choice of time step as discuss in the results section of this
 * program.
 * 
 * @code
 *       double min_temperature = temperature_solution(0),
 *              max_temperature = temperature_solution(0);
 *       for (unsigned int i = 0; i < temperature_solution.size(); ++i)
 *         {
 *           min_temperature =
 *             std::min<double>(min_temperature, temperature_solution(i));
 *           max_temperature =
 *             std::max<double>(max_temperature, temperature_solution(i));
 *         }
 * 
 *       std::cout << "   Temperature range: " << min_temperature << ' '
 *                 << max_temperature << std::endl;
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemoutput_results"></a> 
 * <h4>BoussinesqFlowProblem::output_results</h4>
 *   

 * 
 * This function writes the solution to a VTK output file for visualization,
 * which is done every tenth time step. This is usually quite a simple task,
 * since the deal.II library provides functions that do almost all the job
 * for us. There is one new function compared to previous examples: We want
 * to visualize both the Stokes solution and the temperature as one data
 * set, but we have done all the calculations based on two different
 * DoFHandler objects. Luckily, the DataOut class is prepared to deal with
 * it. All we have to do is to not attach one single DoFHandler at the
 * beginning and then use that for all added vector, but specify the
 * DoFHandler to each vector separately. The rest is done as in step-22. We
 * create solution names (that are going to appear in the visualization
 * program for the individual components). The first <code>dim</code>
 * components are the vector velocity, and then we have pressure for the
 * Stokes part, whereas temperature is scalar. This information is read out
 * using the DataComponentInterpretation helper class. Next, we actually
 * attach the data vectors with their DoFHandler objects, build patches
 * according to the degree of freedom, which are (sub-) elements that
 * describe the data for visualization programs. Finally, we open a file
 * (that includes the time step number) and write the vtk data into it.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::output_results() const
 *   {
 *     if (timestep_number % 10 != 0)
 *       return;
 * 
 *     std::vector<std::string> stokes_names(dim, "velocity");
 *     stokes_names.emplace_back("p");
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       stokes_component_interpretation(
 *         dim + 1, DataComponentInterpretation::component_is_scalar);
 *     for (unsigned int i = 0; i < dim; ++i)
 *       stokes_component_interpretation[i] =
 *         DataComponentInterpretation::component_is_part_of_vector;
 * 
 *     DataOut<dim> data_out;
 *     data_out.add_data_vector(stokes_dof_handler,
 *                              stokes_solution,
 *                              stokes_names,
 *                              stokes_component_interpretation);
 *     data_out.add_data_vector(temperature_dof_handler,
 *                              temperature_solution,
 *                              "T");
 *     data_out.build_patches(std::min(stokes_degree, temperature_degree));
 * 
 *     std::ofstream output("solution-" +
 *                          Utilities::int_to_string(timestep_number, 4) + ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemrefine_mesh"></a> 
 * <h4>BoussinesqFlowProblem::refine_mesh</h4>
 *   

 * 
 * This function takes care of the adaptive mesh refinement. The three tasks
 * this function performs is to first find out which cells to
 * refine/coarsen, then to actually do the refinement and eventually
 * transfer the solution vectors between the two different grids. The first
 * task is simply achieved by using the well-established Kelly error
 * estimator on the temperature (it is the temperature we're mainly
 * interested in for this program, and we need to be accurate in regions of
 * high temperature gradients, also to not have too much numerical
 * diffusion). The second task is to actually do the remeshing. That
 * involves only basic functions as well, such as the
 * <code>refine_and_coarsen_fixed_fraction</code> that refines those cells
 * with the largest estimated error that together make up 80 per cent of the
 * error, and coarsens those cells with the smallest error that make up for
 * a combined 10 per cent of the error.
 *   

 * 
 * If implemented like this, we would get a program that will not make much
 * progress: Remember that we expect temperature fields that are nearly
 * discontinuous (the diffusivity $\kappa$ is very small after all) and
 * consequently we can expect that a freely adapted mesh will refine further
 * and further into the areas of large gradients. This decrease in mesh size
 * will then be accompanied by a decrease in time step, requiring an
 * exceedingly large number of time steps to solve to a given final time. It
 * will also lead to meshes that are much better at resolving
 * discontinuities after several mesh refinement cycles than in the
 * beginning.
 *   

 * 
 * In particular to prevent the decrease in time step size and the
 * correspondingly large number of time steps, we limit the maximal
 * refinement depth of the mesh. To this end, after the refinement indicator
 * has been applied to the cells, we simply loop over all cells on the
 * finest level and unselect them from refinement if they would result in
 * too high a mesh level.
 * 
 * @code
 *   template <int dim>
 *   void
 *   BoussinesqFlowProblem<dim>::refine_mesh(const unsigned int max_grid_level)
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(temperature_dof_handler,
 *                                        QGauss<dim - 1>(temperature_degree + 1),
 *                                        {},
 *                                        temperature_solution,
 *                                        estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.8,
 *                                                       0.1);
 *     if (triangulation.n_levels() > max_grid_level)
 *       for (auto &cell :
 *            triangulation.active_cell_iterators_on_level(max_grid_level))
 *         cell->clear_refine_flag();
 * 
 * @endcode
 * 
 * As part of mesh refinement we need to transfer the solution vectors
 * from the old mesh to the new one. To this end we use the
 * SolutionTransfer class and we have to prepare the solution vectors that
 * should be transferred to the new grid (we will lose the old grid once
 * we have done the refinement so the transfer has to happen concurrently
 * with refinement). What we definitely need are the current and the old
 * temperature (BDF-2 time stepping requires two old solutions). Since the
 * SolutionTransfer objects only support to transfer one object per dof
 * handler, we need to collect the two temperature solutions in one data
 * structure. Moreover, we choose to transfer the Stokes solution, too,
 * since we need the velocity at two previous time steps, of which only
 * one is calculated on the fly.
 *     

 * 
 * Consequently, we initialize two SolutionTransfer objects for the Stokes
 * and temperature DoFHandler objects, by attaching them to the old dof
 * handlers. With this at place, we can prepare the triangulation and the
 * data vectors for refinement (in this order).
 * 
 * @code
 *     std::vector<TrilinosWrappers::MPI::Vector> x_temperature(2);
 *     x_temperature[0]                            = temperature_solution;
 *     x_temperature[1]                            = old_temperature_solution;
 *     TrilinosWrappers::MPI::BlockVector x_stokes = stokes_solution;
 * 
 *     SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> temperature_trans(
 *       temperature_dof_handler);
 *     SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> stokes_trans(
 *       stokes_dof_handler);
 * 
 *     triangulation.prepare_coarsening_and_refinement();
 *     temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
 *     stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);
 * 
 * @endcode
 * 
 * Now everything is ready, so do the refinement and recreate the dof
 * structure on the new grid, and initialize the matrix structures and the
 * new vectors in the <code>setup_dofs</code> function. Next, we actually
 * perform the interpolation of the solutions between the grids. We create
 * another copy of temporary vectors for temperature (now corresponding to
 * the new grid), and let the interpolate function do the job. Then, the
 * resulting array of vectors is written into the respective vector member
 * variables.
 *     

 * 
 * Remember that the set of constraints will be updated for the new
 * triangulation in the setup_dofs() call.
 * 
 * @code
 *     triangulation.execute_coarsening_and_refinement();
 *     setup_dofs();
 * 
 *     std::vector<TrilinosWrappers::MPI::Vector> tmp(2);
 *     tmp[0].reinit(temperature_solution);
 *     tmp[1].reinit(temperature_solution);
 *     temperature_trans.interpolate(x_temperature, tmp);
 * 
 *     temperature_solution     = tmp[0];
 *     old_temperature_solution = tmp[1];
 * 
 * @endcode
 * 
 * After the solution has been transferred we then enforce the constraints
 * on the transferred solution.
 * 
 * @code
 *     temperature_constraints.distribute(temperature_solution);
 *     temperature_constraints.distribute(old_temperature_solution);
 * 
 * @endcode
 * 
 * For the Stokes vector, everything is just the same &ndash; except that
 * we do not need another temporary vector since we just interpolate a
 * single vector. In the end, we have to tell the program that the matrices
 * and preconditioners need to be regenerated, since the mesh has changed.
 * 
 * @code
 *     stokes_trans.interpolate(x_stokes, stokes_solution);
 * 
 *     stokes_constraints.distribute(stokes_solution);
 * 
 *     rebuild_stokes_matrix         = true;
 *     rebuild_temperature_matrices  = true;
 *     rebuild_stokes_preconditioner = true;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemrun"></a> 
 * <h4>BoussinesqFlowProblem::run</h4>
 *   

 * 
 * This function performs all the essential steps in the Boussinesq
 * program. It starts by setting up a grid (depending on the spatial
 * dimension, we choose some different level of initial refinement and
 * additional adaptive refinement steps, and then create a cube in
 * <code>dim</code> dimensions and set up the dofs for the first time. Since
 * we want to start the time stepping already with an adaptively refined
 * grid, we perform some pre-refinement steps, consisting of all assembly,
 * solution and refinement, but without actually advancing in time. Rather,
 * we use the vilified <code>goto</code> statement to jump out of the time
 * loop right after mesh refinement to start all over again on the new mesh
 * beginning at the <code>start_time_iteration</code> label. (The use of the
 * <code>goto</code> is discussed in step-26.)
 *   

 * 
 * Before we start, we project the initial values to the grid and obtain the
 * first data for the <code>old_temperature_solution</code> vector. Then, we
 * initialize time step number and time step and start the time loop.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::run()
 *   {
 *     const unsigned int initial_refinement     = (dim == 2 ? 4 : 2);
 *     const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3);
 * 
 * 
 *     GridGenerator::hyper_cube(triangulation);
 *     global_Omega_diameter = GridTools::diameter(triangulation);
 * 
 *     triangulation.refine_global(initial_refinement);
 * 
 *     setup_dofs();
 * 
 *     unsigned int pre_refinement_step = 0;
 * 
 *   start_time_iteration:
 * 
 *     VectorTools::project(temperature_dof_handler,
 *                          temperature_constraints,
 *                          QGauss<dim>(temperature_degree + 2),
 *                          EquationData::TemperatureInitialValues<dim>(),
 *                          old_temperature_solution);
 * 
 *     timestep_number = 0;
 *     time_step = old_time_step = 0;
 * 
 *     double time = 0;
 * 
 *     do
 *       {
 *         std::cout << "Timestep " << timestep_number << ":  t=" << time
 *                   << std::endl;
 * 
 * @endcode
 * 
 * The first steps in the time loop are all obvious &ndash; we
 * assemble the Stokes system, the preconditioner, the temperature
 * matrix (matrices and preconditioner do actually only change in case
 * we've remeshed before), and then do the solve. Before going on with
 * the next time step, we have to check whether we should first finish
 * the pre-refinement steps or if we should remesh (every fifth time
 * step), refining up to a level that is consistent with initial
 * refinement and pre-refinement steps. Last in the loop is to advance
 * the solutions, i.e., to copy the solutions to the next "older" time
 * level.
 * 
 * @code
 *         assemble_stokes_system();
 *         build_stokes_preconditioner();
 *         assemble_temperature_matrix();
 * 
 *         solve();
 * 
 *         output_results();
 * 
 *         std::cout << std::endl;
 * 
 *         if ((timestep_number == 0) &&
 *             (pre_refinement_step < n_pre_refinement_steps))
 *           {
 *             refine_mesh(initial_refinement + n_pre_refinement_steps);
 *             ++pre_refinement_step;
 *             goto start_time_iteration;
 *           }
 *         else if ((timestep_number > 0) && (timestep_number % 5 == 0))
 *           refine_mesh(initial_refinement + n_pre_refinement_steps);
 * 
 *         time += time_step;
 *         ++timestep_number;
 * 
 *         old_stokes_solution          = stokes_solution;
 *         old_old_temperature_solution = old_temperature_solution;
 *         old_temperature_solution     = temperature_solution;
 *       }
 * @endcode
 * 
 * Do all the above until we arrive at time 100.
 * 
 * @code
 *     while (time <= 100);
 *   }
 * } // namespace Step31
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The main function looks almost the same as in all other programs.
 * 

 * 
 * There is one difference we have to be careful about. This program uses
 * Trilinos and, typically, Trilinos is configured so that it can run in
 * %parallel using MPI. This doesn't mean that it <i>has</i> to run in
 * %parallel, and in fact this program (unlike step-32) makes no attempt at
 * all to do anything in %parallel using MPI. Nevertheless, Trilinos wants the
 * MPI system to be initialized. We do that be creating an object of type
 * Utilities::MPI::MPI_InitFinalize that initializes MPI (if available) using
 * the arguments given to main() (i.e., <code>argc</code> and
 * <code>argv</code>) and de-initializes it again when the object goes out of
 * scope.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step31;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *         argc, argv, numbers::invalid_unsigned_int);
 * 
 * @endcode
 * 
 * This program can only be run in serial. Otherwise, throw an exception.
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
 *                   ExcMessage(
 *                     "This program can only be run in serial, use ./step-31"));
 * 
 *       BoussinesqFlowProblem<2> flow_problem;
 *       flow_problem.run();
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
examples/step-31/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Resultsin2d"></a><h3> Results in 2d </h3>


当你在2D中运行该程序时，输出将看起来像这样。<code> <pre> 活动单元的数量：256（在5层） 自由度的数量：3556（2178+289+1089)

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.919118温度的9次CG迭代。    温度范围：-0.16687 1.30011

活动单元的数量：280（在6层） 自由度的数量：4062（2490+327+1245）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.459559温度的9次CG迭代。    温度范围：-0.0982971 0.598503

活动单元的数量：520（在7个层面上） 自由度的数量：7432（4562+589+2281）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.229779 温度的9次CG迭代。    温度范围：-0.0551098 0.294493

活动单元的数量：1072（在8层） 自由度的数量：15294（9398+1197+4699）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.11489 温度的9次CG迭代。    温度范围：-0.0273524 0.156861

活动单元的数量：2116（在9层） 自由度的数量：30114（18518+2337+9259）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.0574449温度的9次CG迭代。    温度范围：-0.014993 0.0738328

时间步骤1：t=0.0574449 装配...    解决...    斯托克斯子系统的56次GMRES迭代。    时间步长：0.0574449 温度的9次CG迭代。    温度范围：-0.0273934 0.14488

...</pre> </code>

在开始的时候，我们自适应地细化了几次网格，并总是返回到时间步长为零的新细化的网格上重新开始。只有这样，我们才开始实际的时间迭代。

程序运行了一段时间。时间步数为0、500、1000、1500、2000、3000、4000和5000的温度字段看起来是这样的（注意温度使用的色标并不总是相同）。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.07.png" alt="">
    </td>
  </tr>
</table> 

这里显示的视觉效果是使用实例的一个版本生成的，该版本在传输网格后没有强制执行约束。

可以看出，我们有三个加热流体的热源，因此产生了一个浮力效应，让流体的热袋上升并旋转起来。通过烟囱效应，这三股气流被来自外部并想加入上升气流的流体压在一起。请注意，由于流体最初处于静止状态，那些最初在源头上的流体部分比后来被充分发展的流场拖到源头上的流体获得更长的加热时间。因此，它更热，这一事实可以从三个羽流的红色尖端看出。还要注意流场的相对精细的特征，这是我们选择的温度方程的复杂传输稳定的结果。

除了上面的图片外，下面的图片显示了自适应网格和同一时间步长的流场。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.07.png" alt="">
    </td>
  </tr>
</table> 




<a name="Resultsin3d"></a><h3> Results in 3d </h3>


当然，同样的事情也可以在3D中完成，将 <code>main()</code> 中的BoussinesqFlowProblem对象的模板参数从2改为3，这样，现在的输出看起来如下。

<code> <pre> 活动单元的数量：64（在3层） 自由度的数量：3041（2187+125+729）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：2.45098 温度的9次CG迭代。    温度范围：-0.675683 4.94725

活动单元的数量：288（在4层） 自由度的数量：12379（8943+455+2981）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：1.22549 温度的9次CG迭代。    温度范围：-0.527701 2.25764

活动单元的数量：1296（在5层） 自由度的数量：51497（37305+1757+12435）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.612745温度的10次CG迭代。    温度范围：-0.496942 0.847395

活动单元的数量：5048（在6层） 自由度的数量：192425（139569+6333+46523）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.306373 温度的10次CG迭代。    温度范围：-0.267683 0.497739

时间步数1：t=0.306373 正在组装...    解决...    斯托克斯子系统的27次GMRES迭代。    时间步长：0.306373 温度的10次CG迭代。    温度范围：-0.461787 0.958679

...</pre> </code>

在时间步数为0、50、100、150、200、300、400、500、600、700和800的情况下，将温度等值线可视化，得到以下图示。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.07.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.08.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.09.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.10.png" alt="">
    </td>
    <td>
    </td>
  </tr>
</table> 

第一张图片看起来像三只刺猬，这是因为我们的方案基本上是将源乘以第一时间步长投射到网格上，以获得第一时间步的温度场。由于源函数是不连续的，我们需要期待这个项目的过冲和欠冲。这就是事实上发生的情况（在2d中更容易检查），并导致等值面的皱缩外观。  这里显示的视觉效果是使用例子的一个版本生成的，该版本在传输网格后没有强制执行约束。




<a name="Numericalexperimentstodetermineoptimalparameters"></a><h3> Numerical experiments to determine optimal parameters </h3>


现在的程序有三个参数，我们在理论上并没有掌握如何以最佳方式进行选择。这三个参数是。   <ul>   <li>  时间步骤必须满足CFL条件  $k\le \min_K \frac{c_kh_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。这里， $c_k$ 是无量纲的，但什么是正确的值？     <li>  在计算人工黏度时。

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
  \min\left\{
    h_K,
    h_K^\alpha
    \frac{\|R_\alpha(T)\|_{L^\infty(K)}}{c(\mathbf{u},T)}
  \right\},


@f}

      与 $c(\mathbf{u},T) =
      c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
      \ |\mathrm{diam}(\Omega)|^{\alpha-2}$  。       这里，无量纲%数 $\beta,c_R$ 的选择是有意义的。   </ul>  在所有这些情况下，我们将不得不期望每个值的正确选择取决于其他值的正确选择，而且很可能也取决于用于温度的有限元的空间尺寸和多项式程度。下面我们将讨论一些数值实验来选择常数  $c_k$  和  $\beta$  。

下面，我们将不讨论 $c_R$ 的选择问题。在程序中，我们将其设定为 $c_R=2^{\frac{4-2\alpha}{d}}$  。这个值的原因有点复杂，与程序的历史而不是推理有关：虽然全局缩放参数 $c(\mathbf{u},T)$ 的正确公式如上所示，但程序（包括与deal.II 6.2一起出厂的版本）最初有一个错误，即我们计算的是 $c(\mathbf{u},T) =
      \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
      \ \frac{1}{|\mathrm{diam}(\Omega)|^{\alpha-2}}$ ，而在这里我们将缩放参数设置为1。由于我们只在 $\mathrm{diam}(\Omega)=2^{1/d}$ 的单位平方/立方体上进行计算，这完全等同于使用 $c_R=\left(2^{1/d}\right)^{4-2\alpha}=2^{\frac{4-2\alpha}{d}}$ 的正确公式。由于 $c_R$ 的这个值对于当前的程序来说似乎很好用，我们在程序中修正了公式，并将 $c_R$ 设置为一个值，正好再现了我们之前的结果。不过，我们将在第32步中再次审视这个问题。

然而，现在回到讨论 $c_k$ 和 $\beta$ 的什么值来选择。




<a name="Choosingicsubksubiicsubksubiandbeta"></a><h4> Choosing <i>c<sub>k</sub></i><i>c<sub>k</sub></i> and beta </h4> 。


这两个常数肯定在某种程度上有联系。原因很容易看出来。在纯平流问题的情况下， $\frac{\partial T}{\partial t} + \mathbf{u}\cdot\nabla T = \gamma$ ，任何显式方案都必须满足形式为 $k\le \min_K \frac{c_k^a h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 的CFL条件。另一方面，对于纯扩散问题， $\frac{\partial T}{\partial t} + \nu \Delta T = \gamma$ ，显式方案需要满足一个条件 $k\le \min_K \frac{c_k^d h_K^2}{\nu}$ 。因此，鉴于上述 $\nu$ 的形式，像我们这里要解决的平流扩散问题将导致一个 $
k\le \min_K \min \left\{
  \frac{c_k^a h_K}{\|\mathbf{u}\|_{L^\infty(K)}},
  \frac{c_k^d h_K^2}{\beta \|\mathbf{u}\|_{L^\infty(K)} h_K}\right\}
  =
  \min_K \left( \min \left\{
  c_k^a,
  \frac{c_k^d}{\beta}\right\}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}} \right)
$ 的条件。因此，我们必须面对这样一个事实：我们可能想选择 $\beta$ 大一些，以提高数值方案的稳定性（通过增加人工扩散量），但我们必须以更小的、因而更多的时间步骤为代价。因此，在实践中，人们希望尽可能地选择 $\beta$ ，以保持传输问题的充分稳定，同时尽量选择大的时间步长，以减少总体工作量。

要找到正确的平衡，唯一的办法是做一些计算实验。下面是我们的做法。我们稍微修改了程序，允许更少的网格细化（所以我们不一定要等那么久），并选择 $
  \nu(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)} h_K
$ 来消除常数 $c_R$ 的影响（我们知道通过使用这个版本的 $\nu(T)$ 作为人工粘度，解决方案是稳定的，但我们可以通过使用这个人工粘度的更复杂的公式来改善情况--即使解决方案更清晰）。然后我们对不同的值 $c_k,\beta$ 运行程序，观察域中的最大和最小温度。我们期望看到的情况是这样的。如果我们选择的时间步长过大（即选择一个比理论上允许的大的 $c_k$ ），那么我们将得到温度的指数式增长。如果我们选择 $\beta$ 太小，那么传输稳定变得不充分，解决方案将显示出明显的振荡，但不是指数级增长。




<a name="ResultsforQsub1subelements"></a><h5>Results for Q<sub>1</sub> elements</h5>


下面是我们对 $\beta=0.01, \beta=0.1$ ，和 $\beta=0.5$ ， $c_k$ 的不同选择，以及2d的双线性元素（ <code>temperature_degree=1</code> ）得到的结果。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.1.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.5.png" alt="">
    </td>
  </tr>
</table> 

解释这些图表的方法是这样的：对于 $\beta=0.01$ 和 $c_k=\frac 12,\frac 14$ ，我们看到指数增长或至少是大的变化，但如果我们选择 $k=\frac 18\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小，那么这个方案虽然有点摇摆不定，但还是稳定的。对于更多的人工扩散，我们可以选择 $k=\frac 14\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小的 $\beta=0.03$ ， $k=\frac 13\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小的 $\beta=0.1$ ，并再次需要 $k=\frac 1{15}\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 的 $\beta=0.5$ （这次是因为许多扩散需要一个小的时间步长）。

那么该如何选择呢？如果我们只是对大时间步长感兴趣，那么我们会选择 $\beta=0.1$ 和 $k=\frac 13\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。另一方面，我们也对准确性感兴趣，在这里，实际调查这些曲线所显示的内容可能会有兴趣。为此，请注意，我们从零温度开始，我们的来源是正的&mdash；所以我们会直观地期望温度永远不会降到零以下。但它确实如此，这是使用连续元素来近似不连续的解决方案时，吉布现象的结果。因此，我们可以看到，选择 $\beta$ 太小是不好的：太少的人工扩散会导致没有扩散掉的过冲和欠冲。另一方面，对于大的 $\beta$ ，最低温度在开始时下降到零以下，但随后迅速扩散回零。

另一方面，我们也来看看最高温度。观察溶液的电影，我们看到最初流体处于静止状态。源头不断加热相同体积的流体，其温度在开始时呈线性增长，直到其浮力能够使其向上移动。因此，流体中最热的部分被带离了溶液，取而代之的流体只被加热了很短的时间就被移出了源区，因此仍然比初始气泡要冷。如果 $\kappa=0$ （在程序中是非零的，但非常小），那么流体中最热的部分应该随着流动而平移，其温度不变。这就是我们在最小的 $\beta$ 图中可以看到的：一旦达到最高温度，它就几乎不再变化。另一方面，人工扩散越大，热点的扩散就越多。请注意，对于这个标准，时间步长的大小并不发挥重要作用。

因此，总结起来，可能最好的选择似乎是 $\beta=0.03$ 和 $k=\frac 14\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。曲线有点摇摆不定，但总的来说，图片看起来相当合理，除了由于吉布现象而在接近开始时间时出现一些过冲和欠冲的情况。




<a name="ResultsforQsub2subelements"></a><h5>Results for Q<sub>2</sub> elements</h5>


我们也可以对高阶元素重复同样的实验序列。这里是温度的双二次方形状函数（ <code>temperature_degree=2</code> ）的图形，同时我们保留了斯托克斯系统的 $Q_2/Q_1$ 稳定泰勒-胡德元素。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.1.png" alt="">
    </td>
  </tr>
</table> 

同样， $\beta$ 的小值会导致较少的扩散，但我们必须选择非常小的时间步长来保持事情的控制。太大的 $\beta$ 值会导致更多的扩散，但同样需要小的时间步骤。最佳值似乎是 $\beta=0.03$ ，和 $Q_1$ 元素一样，然后我们必须选择 $k=\frac 18\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ &mdash；正好是 $Q_1$ 元素的一半大小。]元素，如果我们把CFL条件说成是要求时间步长足够小，以便运输在每个时间步长中的移动距离不超过一个<i>grid point</i>距离（对于 $Q_1$ 元素是 $h_K$ ，但对于 $Q_2$ 元素是 $h_K/2$ ），这个事实可能并不令人惊讶。事实证明， $\beta$ 需要稍微大一点，以便在模拟后期获得稳定的结果，时间大于60，所以我们实际上在代码中选择它作为 $\beta = 0.034$ 。




<a name="Resultsfor3d"></a><h5>Results for 3d</h5>


我们可以在3D中重复这些实验，找到每个 $\beta$ 值的最佳时间步骤，并找到 $\beta$ 的最佳值。人们发现，对于2d中已经使用的相同的 $\beta$ ，时间步长需要小一点，大约是1.2倍左右。这很容易解释：时间步长的限制是 $k=\min_K \frac{ch_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ ，其中 $h_K$ 是单元的<i>diameter</i>。然而，真正需要的是网格点之间的距离，它是 $\frac{h_K}{\sqrt{d}}$  。所以更合适的形式是  $k=\min_K \frac{ch_K}{\|\mathbf{u}\|_{L^\infty(K)}\sqrt{d}}$  。

第二个发现是，需要把 $\beta$ 选得稍微大一点（大约 $\beta=0.05$ 左右）。这就再次减少了我们可以采取的时间步骤。







<a name="Conclusions"></a><h5>Conclusions</h5>


总之，从上面的简单计算来看， $\beta=0.034$ 似乎是2D中稳定参数的一个好选择，而 $\beta=0.05$ 则是3D中的稳定参数。以独立于维度的方式，我们可以将其建模为 $\beta=0.017d$  。如果在更细的网格上做更长时间的计算（几千个时间步长），就会意识到时间步长还不够小，为了稳定，就必须把上述数值再降低一些（大约是 $\frac 78$ 的一个系数）。

因此，调和2D、3D和可变多项式程度并考虑到所有因素的公式如下。

@f{eqnarray*}
  k =
  \frac 1{2 \cdot 1.7} \frac 1{\sqrt{d}}
  \frac 2d
  \frac 1{q_T}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}
  =
  \frac 1{1.7 d\sqrt{d}}
  \frac 1{q_T}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}.


@f}

在第一种形式中（方程中心）， $\frac
1{2 \cdot 1.7}$ 是一个通用常数， $\frac 1{\sqrt{d}}$ 是说明单元直径和网格点间距的因素， $\frac 2d$ 说明 $\beta$ 随着空间尺寸的增加而增加， $\frac 1{q_T}$ 说明高阶元素的网格点之间的距离， $\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 说明相对于单元尺寸的局部传输速度。这就是我们在程序中使用的公式。

至于对温度使用 $Q_1$ 或 $Q_2$ 元素的问题，以下考虑可能是有用的。首先，解决温度方程在整个方案中几乎不是一个因素，因为几乎所有的计算时间都用于解决每个时间步骤中的斯托克斯系统。因此，温度方程的高阶元素并不是一个重要的缺点。另一方面，如果比较一下由于不连续的源描述而产生的过冲和欠冲的大小，我们会注意到，对于上述 $\beta$ 和 $k$ 的选择， $Q_1$ 的解决方案下降到 $-0.47$ 左右，而 $Q_2$ 的解决方案只到 $-0.13$ （记住，精确解决方案根本不应该变成负数。这意味着 $Q_2$ 解明显更准确；因此程序使用这些高阶元素，尽管我们在较小的时间步长方面付出了代价。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


有各种方法来扩展当前的程序。当然，特别感兴趣的是使其更快和/或提高程序的分辨率，特别是在3D方面。这就是step-32教程程序的主题，它将实现在集群上以%并行方式解决这个问题的策略。它也是更大的开放源代码ASPECT（见https://aspect.geodynamics.org/）的基础，它可以解决现实问题，并构成step-32的进一步发展。

另一个方向是使流体流动更加真实。这个程序最初是为了模拟各种情况，模拟地幔中的物质对流，即外地核和固体地壳之间的区域：在那里，物质从下面被加热，从上面被冷却，导致热对流。然而，这种流体的物理学要比这个程序中显示的复杂得多。地幔材料的粘度与温度有很大的关系，即 $\eta=\eta(T)$ ，这种关系经常被模拟为粘度随温度升高而呈指数下降。其次，地幔的大部分动态是由化学反应决定的，主要是构成地幔的各种晶体的相变；然后，斯托克斯方程右边的浮力项不仅取决于温度，而且还取决于某个特定位置的化学成分，这些化学成分被流场平流，但也作为压力和温度的函数而变化。我们将在以后的教程程序中也研究其中的一些影响。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-31.cc"
*/
