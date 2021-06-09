/**
@page step_67 The step-67 tutorial program
This tutorial depends on step-33, step-48, step-59.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#TheEulerequations">The Euler equations</a>
        <li><a href="#HighorderdiscontinuousGalerkindiscretization">High-order discontinuous Galerkin discretization</a>
        <li><a href="#Explicittimeintegration">Explicit time integration</a>
        <li><a href="#Fastevaluationofintegralsbymatrixfreetechniques">Fast evaluation of integrals by matrix-free techniques</a>
        <li><a href="#Evaluationoftheinversemassmatrixwithmatrixfreetechniques">Evaluation of the inverse mass matrix with matrix-free techniques</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#LowstorageexplicitRungeKuttatimeintegrators">Low-storage explicit Runge&mdash;Kutta time integrators</a>
        <li><a href="#ImplementationofpointwiseoperationsoftheEulerequations">Implementation of point-wise operations of the Euler equations</a>
        <li><a href="#TheEulerOperationclass">The EulerOperation class</a>
      <ul>
        <li><a href="#Localevaluators">Local evaluators</a>
        <li><a href="#Theapplyandrelatedfunctions">The apply() and related functions</a>
      </ul>
        <li><a href="#TheEulerProblemclass">The EulerProblem class</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Convergenceratesfortheanalyticaltestcase">Convergence rates for the analytical test case</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin2D">Results for flow in channel around cylinder in 2D</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin3D">Results for flow in channel around cylinder in 3D</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Moreadvancednumericalfluxfunctionsandskewsymmetricformulations">More advanced numerical flux functions and skew-symmetric formulations</a>
        <li><a href="#Equippingthecodeforsupersoniccalculations">Equipping the code for supersonic calculations</a>
        <li><a href="#ExtensiontothelinearizedEulerequations">Extension to the linearized Euler equations</a>
        <li><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a>
        <li><a href="#Usingcellcentricloopsandsharedmemory">Using cell-centric loops and shared memory</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-67/doc/intro.dox



 <br> 

<i>
This program was contributed by Martin Kronbichler. Many ideas presented here
are the result of common code development with Niklas Fehn, Katharina Kormann,
Peter Munch, and Svenja Schoeder.


This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA).
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序使用显式时间积分器求解流体力学的欧拉方程，其无矩阵框架应用于空间的高阶非连续Galerkin离散化。关于欧拉系统的细节和另一种隐式方法，我们也参考了第33步教程程序。你可能还想看看第69步，看看解决这些方程的另一种方法。




<a name="TheEulerequations"></a><h3>The Euler equations</h3>


欧拉方程是一个守恒定律，描述了一个可压缩的无粘性气体的运动。

@f[
\frac{\partial \mathbf{w}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{w}) =
\mathbf{G}(\mathbf w),


@f]

其中解向量的 $d+2$ 分量为 $\mathbf{w}=(\rho, \rho
u_1,\ldots,\rho u_d,E)^{\mathrm T}$  。这里， $\rho$  表示流体密度， ${\mathbf u}=(u_1,\ldots, u_d)^\mathrm T$  表示流体速度， $E$  表示气体的能量密度。速度不直接求解，而是用变量 $\rho \mathbf{u}$ ，即线性动量（因为这是一个守恒量）。

欧拉通量函数是一个 $(d+2)\times d$ 矩阵，定义为

@f[
  \mathbf F(\mathbf w)
  =
  \begin{pmatrix}
  \rho \mathbf{u}\\
  \rho \mathbf{u} \otimes \mathbf{u} + \mathbb{I}p\\
  (E+p)\mathbf{u}
  \end{pmatrix}


@f]

其中 $\mathbb{I}$ 为 $d\times d$ 身份矩阵， $\otimes$ 为外积；其组成部分分别表示质量、动量和能量通量。右手边的强制力由以下公式给出

@f[
  \mathbf G(\mathbf w)
  =
  \begin{pmatrix}
  0\\
  \rho\mathbf{g}\\
  \rho \mathbf{u} \cdot \mathbf{g}
  \end{pmatrix},


@f]

其中矢量 $\mathbf g$ 表示重力的方向和大小。然而，它也可以表示作用于流体的任何其他单位质量的外力。例如，想想外部电场对带电粒子所施加的静电力）。

这三块方程，第二块涉及 $d$ 成分，描述了质量、动量和能量的守恒。压力不是一个解决方案的变量，但需要通过其他变量的 "闭合关系 "来表达；我们在此选择适合由两个原子组成的分子的气体的关系，在中等温度下，由 $p=(\gamma - 1) \left(E-\frac 12 \rho
\mathbf{u}\cdot \mathbf{u}\right)$ 和常数 $\gamma = 1.4$ 给出。




<a name="HighorderdiscontinuousGalerkindiscretization"></a><h3>High-order discontinuous Galerkin discretization</h3>


对于空间离散化，我们使用高阶非连续加勒金（DG）离散化，使用的解扩展形式为

@f[
\mathbf{w}_h(\mathbf{x}, t) =
\sum_{j=1}^{n_\mathbf{dofs}} \boldsymbol{\varphi}_j(\mathbf{x}) {w}_j(t).


@f]

这里， $\boldsymbol{\varphi}_j$ 表示第 $j$ 个基函数，以矢量形式写出不同成分的独立形状函数，让 $w_j(t)$ 分别通过密度、动量和能量变量。在这种形式下，空间依赖性包含在形状函数中，时间依赖性包含在未知系数中  $w_j$  。与连续有限元方法中一些形状函数跨越元素边界不同，在DG方法中，形状函数是单个元素的局部，从一个元素到下一个元素是不连续的。从一个单元到其相邻单元的解的连接是由下面规定的数值通量来实现的。这允许一些额外的灵活性，例如，在数值方法中引入方向性，例如，上卷。

DG方法是解决传输特性问题的流行方法，因为它们结合了低分散误差和勉强解决的尺度上的可控耗散。这使得它们在流体动力学领域的模拟中特别有吸引力，因为在这个领域中，需要代表广泛的活动尺度，不充分解决的特征很容易干扰重要的良好解决的特征。此外，高阶DG方法非常适用于现代硬件的正确实施。同时，DG方法也不是万能的。特别是当解出现不连续（冲击）时，就像欧拉方程在某些流态下的典型情况一样，高阶DG方法容易出现振荡解，就像所有不使用通量或坡度限制器的高阶方法一样。这是<a
href="https://en.wikipedia.org/wiki/Godunov%27s_theorem">Godunov's theorem</a>的结果，即任何线性的总变差（TVD）方案（如基本的DG离散化）最多只能达到一阶精度。换句话说，由于DG方法的目标是高阶精度，因此它们不可能对出现冲击的解进行TVD。尽管有些人声称DG方法中的数值通量可以控制耗散，但除非问题中的<b>all</b>冲击与单元边界对齐，否则这一点的价值有限。任何穿过单元内部的冲击都会因为高阶多项式而再次产生振荡分量。在有限元和DG界，存在许多不同的方法来处理冲击，例如在有问题的单元上引入人工扩散（使用基于解的模态分解等的有问题单元指标），在子网格上转换为耗散性低阶有限体积方法，或者增加一些限制程序。考虑到这种情况下的大量可能性，再加上相当大的实施努力，我们在这里不考虑带有明显冲击的欧拉方程系统，而是集中在带有波浪状现象的亚音速流动系统。对于一个能很好地处理冲击的方法（但每个未知数的成本较高），我们可以参考step-69教程程序。

对于DG公式的推导，我们将欧拉方程与测试函数 $\mathbf{v}$ 相乘，并对单个单元进行积分 $K$ ，从而得到

@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}
+ \left(\mathbf{v}, \nabla \cdot \mathbf{F}(\mathbf{w})\right)_{K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.


@f]



然后我们对第二项进行分项积分，将分歧从解槽移到测试函数槽，并产生一个元素边界上的积分。

@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}


- \left(\nabla \mathbf{v}, \mathbf{F}(\mathbf{w})\right)_{K}
+ \left<\mathbf{v}, \mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})
\right>_{\partial K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.


@f]

在表面积分中，我们用术语 $\widehat{\mathbf{F}}(\mathbf w)$ 代替了术语 $\mathbf{F}(\mathbf w)$ ，即数值通量。数字通量的作用是连接相邻元素上的解，并弱化解的连续性。这保证了PDE的全局耦合反映在离散化中，尽管单元上有独立的基函数。通过将数值通量定义为来自内部面两侧的解的函数 $\widehat{\mathbf{F}}(\mathbf w^-,
\mathbf w^+)$ 和 $\mathbf w^+$ ，包括与邻居的连接。我们要求的一个基本属性是，数值通量需要是<b>conservative</b>。也就是说，我们希望所有的信息（即质量、动量和能量）在一个面上离开一个单元时，都能完整地进入邻近的单元，反之亦然。这可以表示为 $\widehat{\mathbf{F}}(\mathbf w^-, \mathbf w^+) =
\widehat{\mathbf{F}}(\mathbf w^+, \mathbf w^-)$ ，也就是说，数值通量从任何一边都评估为相同的结果。结合数值通量与所考虑的面的单位外法向量相乘的事实，即从两边指向相反的方向，我们看到守恒被满足了。数值通量的另一个观点是作为一个单值的中间状态，从两边微弱地连接解决方案。

有大量的数值通量函数可用，也称为黎曼解算器。对于欧拉方程，存在所谓的精确黎曼求解器--意味着来自双方的状态以一种与欧拉方程沿线不连续的方式结合起来--以及近似黎曼求解器，它违反了一些物理特性，并依靠其他机制来使方案总体上准确。近似黎曼求解器的优点是计算起来比较便宜。大多数通量函数都起源于有限体积界，它们类似于单元（称为体积）内的多项式0度的DG方法。由于欧拉算子 $\mathbf{F}$ 的体积积分对于恒定解和检验函数会消失，所以数值通量必须完全代表物理算子，这也解释了为什么该界有大量的研究。对于DG方法，一致性是由单元内的高阶多项式保证的，这使得数值通量不再是一个问题，通常只影响收敛率，例如，对于度数为 $\mathcal O(h^p)$ 的多项式，解是否收敛为 $\mathcal O(h^{p+1/2})$ 或 $\mathcal
O(h^{p+1})$ 的准则。因此，数值通量可以被看作是一种机制，用于选择更有利的耗散/分散特性或关于离散化和线性化算子的极值特征，这影响到显式时间积分器中最大的可接受的时间步长。

在这个教程程序中，我们实现了两种通量的变体，可以通过程序中的开关来控制（当然，要使它们成为通过输入文件控制的运行时参数也很容易）。第一个通量是本地的Lax--Friedrichs通量

@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
   \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
   \mathbf{n^-}.


@f]



在Lax--Friedrichs通量的原始定义中，使用了一个系数 $\lambda =
\max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ （对应于信息在界面两边移动的最大速度），说明两个状态之间的差异， $[\![\mathbf{w}]\!]$ 被欧拉通量中的最大特征值惩罚，即 $\|\mathbf{u}\|+c$  ，其中 $c=\sqrt{\gamma p / \rho}$  是音速。在下面的实现中，我们对惩罚项进行了一些修改，因为无论如何惩罚都是近似的。我们使用

@f{align*}{
\lambda
&=
\frac{1}{2}\max\left(\sqrt{\|\mathbf{u^-}\|^2+(c^-)^2},
                     \sqrt{\|\mathbf{u}^+\|^2+(c^+)^2}\right)
\\
&=
\frac{1}{2}\sqrt{\max\left(\|\mathbf{u^-}\|^2+(c^-)^2,
                           \|\mathbf{u}^+\|^2+(c^+)^2\right)}.


@f}

额外的因子 $\frac 12$ 降低了惩罚强度（这导致特征值的负实部减少，从而增加了可接受的时间步长）。使用和内的平方允许我们减少昂贵的平方根操作的数量，对于原始的Lax--Friedrichs定义是4个，现在只需要一个。这种简化导致参数 $\lambda$ 的减少最多为2倍，因为 $\|\mathbf{u}\|^2+c^2 \leq
\|\mathbf{u}\|^2+2 c |\mathbf{u}\| + c^2 = \left(\|\mathbf{u}\|+c\right)^2
\leq 2 \left(\|\mathbf{u}\|^2+c^2\right)$ ，最后一个不等式来自杨氏不等式。

第二个数值通量是由Harten、Lax和van Leer提出的，称为HLL通量。它考虑到欧拉方程的不同传播方向，取决于声速。它利用一些中间状态  $\bar{\mathbf{u}}$  和  $\bar{c}$  来定义两个分支  $s^\mathrm{p} = \max\left(0, \bar{\mathbf{u}}\cdot \mathbf{n} +
\bar{c}\right)$  和  $s^\mathrm{n} = \min\left(0, \bar{\mathbf{u}}\cdot
\mathbf{n} - \bar{c}\right)$  。从这些分支中，人们再定义出通量

@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{s^\mathrm{p} \mathbf{F}(\mathbf{w}^-)-s^\mathrm{n} \mathbf{F}(\mathbf{w}^+)}
                   {s^\mathrm p - s^\mathrm{n} } +
\frac{s^\mathrm{p} s^\mathrm{n}}{s^\mathrm{p}-s^\mathrm{n}}
\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes \mathbf{n^-}.


@f]

关于中间状态的定义  $\bar{\mathbf{u}}$  和  $\bar{c}$  ，已经提出了几个变种。最初提出的变体使用密度平均的速度定义，  $\bar{\mathbf{u}}
= \frac{\sqrt{\rho^-} \mathbf{u}^- + \sqrt{\rho^+}\mathbf{u}^+}{\sqrt{\rho^-}
+ \sqrt{\rho^+}}$  。由于我们考虑的是没有冲击的欧拉方程，因此在本教程程序中，我们简单地使用算术平均值， $\bar{\mathbf{u}} = \frac{\mathbf{u}^- +
\mathbf{u}^+}{2}$ 和 $\bar{c} = \frac{c^- + c^+}{2}$ ，与 $c^{\pm} =
\sqrt{\gamma p^{\pm} / \rho^{\pm}}$ ，而把其他变体留给可能的扩展。我们还注意到，HLL通量在文献中被扩展为所谓的HLLC通量，其中C代表表示接触不连续的能力。

在没有邻接状态 $\mathbf{w}^+$ 的边界上，通常的做法是从边界条件中推导出合适的外部值（详见关于DG方法的一般文献）。在这个教程程序中，我们考虑三种类型的边界条件，即<b>inflow boundary conditions</b>，其中所有分量都是规定的。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho_\mathrm{D}(t)\\
(\rho \mathbf u)_{\mathrm D}(t) \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(Dirichlet)},


@f]

<b>subsonic outflow boundaries</b>，在这里我们不规定外部解，因为流场要离开域，而使用内部值；我们仍然需要规定能量，因为欧拉通量中还有一个传入特性。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(mixed Neumann/Dirichlet)},


@f]

和<b>wall boundary condition</b>，它们描述了无渗透配置。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- - 2 [(\rho \mathbf u)^-\cdot \mathbf n] \mathbf{n}
 \\ E^-\end{pmatrix}.


@f]



解的多项式展开最后被插入到弱形式，测试函数被基函数取代。这就得到了一个空间上离散、时间上连续的非线性系统，其未知系数的数量有限  $w_j$  ,  $j=1,\ldots,n_\text{dofs}$  。关于DG方法中多项式度数的选择，截至2019年，文献中并没有关于什么多项式度数最有效的共识，决定取决于问题。高阶多项式可以确保更好的收敛率，因此对于中等到高精确度要求的<b>smooth</b>解来说，高阶多项式更有优势。同时，自由度所在的体积与表面的比率，随着高阶的增加而增加，这使得数值通量的影响变弱，通常会减少耗散。然而，在大多数情况下，解决方案是不平滑的，至少与可以承受的分辨率相比是不平滑的。例如，在不可压缩流体力学、可压缩流体力学以及与之相关的波浪传播课题中都是如此。在这个前渐进制度中，误差大约与数值分辨率成正比，而其他因素，如分散误差或耗散行为变得更加重要。非常高阶的方法往往被排除在外，因为它们带有根据未知数衡量的更多限制性的CFL条件，而且当涉及到表示复杂几何形状时，它们也不那么灵活。因此，2到6的多项式度数在实践中是最受欢迎的，例如见 @cite FehnWallKronbichler2019 中的效率评估和其中引用的参考文献。

<a name="Explicittimeintegration"></a><h3>Explicit time integration</h3>


为了进行时间离散化，我们稍微重新排列了弱的形式，并在所有单元上求和。

@f[
\sum_{K \in \mathcal T_h} \left(\boldsymbol{\varphi}_i,
\frac{\partial \mathbf{w}}{\partial t}\right)_{K}
=
\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w})\right)_{K}


-\left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})\right>_{\partial K} +
\left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w)\right)_{K}
\right],


@f]

其中 $\boldsymbol{\varphi}_i$ 贯穿了从1到 $n_\text{dofs}$ 的所有基函数。

我们现在用 $\mathcal M$ 表示质量矩阵，其条目为 $\mathcal M_{ij} =
\sum_{K} \left(\boldsymbol{\varphi}_i,
\boldsymbol{\varphi}_j\right)_K$ ，并用

@f[
\mathcal L_h(t,\mathbf{w}_h) = \left[\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w}_h)\right)_{K}


- \left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w}_h)\right>_{\partial K}
+ \left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w_h)\right)_{K}
\right]\right]_{i=1,\ldots,n_\text{dofs}}.


@f]

给定一个与全局未知数矢量和使用中的有限元相关的函数 $\mathbf{w}_h$ ，对欧拉算子的右手边进行评估的算子。这个函数 $\mathcal L_h$ 是明确随时间变化的，因为在边界上评估的数值通量将涉及边界某些部分的随时间变化的数据 $\rho_\mathrm{D}$ 、 $(\rho \mathbf{u})_\mathrm{D}$ 和 $E_\mathbf{D}$ ，取决于边界条件的分配。有了这个符号，我们可以把空间上的离散、时间上的连续系统紧凑地写为

@f[
\mathcal M \frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal L_h(t, \mathbf{w}_h),


@f]

其中我们冒昧地用 $\mathbf{w}_h$ 表示全局解矢量（除了相应的有限元函数外）。等价地，上述系统的形式为

@f[
\frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal M^{-1} \mathcal L_h(t, \mathbf{w}_h).


@f]



对于用高阶非连续Galerkin方法离散的双曲系统，该系统的显式时间积分非常流行。这是由于质量矩阵 $\mathcal M$ 是块对角线的（每个块只对应于定义在同一单元上的同类变量），因此很容易倒置。在每个时间步长--或Runge-Kutta方案的阶段--我们只需要用给定的数据评估一次微分算子，然后应用质量矩阵的逆。另一方面，对于隐式时间步进，人们首先必须将方程线性化，然后迭代解决线性系统，这涉及到几个残差评估和至少十几个线性化算子的应用，正如在步骤33教程程序中所展示的那样。

当然，显式时间步长的简单性是有代价的，即由于所谓的Courant-Friedrichs-Lewy（CFL）条件而产生的条件稳定性。它指出，时间步长不能大于离散微分算子的最快信息传播速度。用更现代的术语来说，传播速度对应于离散算子的最大特征值，反过来又取决于网格大小、多项式程度 $p$ 和欧拉算子的物理学，即 $\mathbf F(\mathbf w)$ 相对于 $\mathbf{w}$ 的线性化的特征值。在这个程序中，我们设定的时间步长如下。

@f[
\Delta t = \frac{\mathrm{Cr}}{p^{1.5}}\left(\frac{1}
           {\max\left[\frac{\|\mathbf{u}\|}{h_u} + \frac{c}{h_c}\right]}\right),


@f]



在所有正交点和所有单元中取最大值。无量纲数 $\mathrm{Cr}$ 表示库朗数，可以选择最大稳定数 $\mathrm{Cr}_\text{max}$ ，其值取决于所选择的时间步进方法及其稳定性。用于多项式缩放的幂 $p^{1.5}$ 是启发式的，代表1到8之间的多项式度数最接近，例如，见 @cite SchoederKormann2018 。在更高的度数限制下， $p>10$ ， $p^2$ 的比例更准确，与通常用于内部惩罚方法的逆向估计有关。关于公式中使用的<i>effective</i>网格尺寸 $h_u$ 和 $h_c$ ，我们注意到对流传输是定向的。因此，一个合适的比例是使用速度方向的元素长度  $\mathbf u$  。下面的代码从参考单元到实际单元的雅各布系数的倒数得出这个比例，也就是说，我们近似于  $\frac{\|\mathbf{u}\|}{h_u} \approx \|J^{-1} \mathbf
u\|_{\infty}$  。相反，声波具有各向同性的特点，这就是为什么我们使用最小的特征尺寸，由 $J$ 的最小奇异值代表，用于声学缩放  $h_c$  。最后，我们需要增加对流和声学限制，因为欧拉方程可以以速度传输信息  $\|\mathbf{u}\|+c$  。

在这个教程程序中，我们使用<a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">explicit
Runge--Kutta methods</a>的一个特定变体，一般来说，它使用以下更新程序，从时间 $t^n$ 的状态 $\mathbf{w}_h^{n}$ 到新时间 $t^{n+1}$ 的 $\Delta t = t^{n+1}-t^n$  。

@f[
\begin{aligned}
\mathbf{k}_1 &= \mathcal M^{-1} \mathcal L_h\left(t^n, \mathbf{w}_h^n\right),
\\
\mathbf{k}_2 &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_2\Delta t,
                       \mathbf{w}_h^n + a_{21} \Delta t \mathbf{k}_1\right),
\\
&\vdots \\
\mathbf{k}_s &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_s\Delta t,
  \mathbf{w}_h^n + \sum_{j=1}^{s-1} a_{sj} \Delta t \mathbf{k}_j\right),
\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^n + \Delta t\left(b_1 \mathbf{k}_1 +
b_2 \mathbf{k}_2 + \ldots + b_s \mathbf{k}_s\right).
\end{aligned}


@f]

在 $\mathbf{k}_i$ 、 $i=1,\ldots,s$ 的阶段性方案中，向量 $s$ 是算子在某个中间状态下的评价，并通过某种线性组合用于定义阶段性结束值 $\mathbf{w}_h^{n+1}$ 。该方案中的标量系数 $c_i$ 、 $a_{ij}$ 和 $b_j$ 的定义，使得高阶方案满足某些条件，最基本的是 $c_i = \sum_{j=1}^{i-1}a_{ij}$  。参数通常以所谓的<a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge%E2%80%93Kutta_methods">Butcher
tableau</a>的形式收集，它收集了定义该方案的所有系数。对于一个五级方案，它看起来是这样的。

@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_{21} \\
c_3 & a_{31} & a_{32} \\
c_4 & a_{41} & a_{42} & a_{43} \\
c_5 & a_{51} & a_{52} & a_{53} & a_{54} \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}


@f]



在这个教程程序中，我们使用显式Runge--Kutta方法的一个子集，即所谓的低存储Runge--Kutta方法（LSRK），它假定了系数的额外结构。在参考文献 @cite KennedyCarpenterLewis2000 所使用的变体中，假设使用的是Butcher tableaus的形式

@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_1 \\
c_3 & b_1 & a_2 \\
c_4 & b_1 & b_2 & a_3 \\
c_5 & b_1 & b_2 & b_3 & a_4 \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}


@f]

有了这样的定义，对  $\mathbf{w}_h^n$  的更新与中间值  $\mathbf{k}_i$  的信息共享存储。从 $\mathbf{w}^{n+1}=\mathbf{w}^n$ 和 $\mathbf{r}_1 = \mathbf{w}^n$ 开始，每个 $s$ 阶段的更新都简化为

@f[
\begin{aligned}
\mathbf{k}_i &=
\mathcal M^{-1} \mathcal L_h\left(t^n+c_i\Delta t, \mathbf{r}_{i} \right),\\
\mathbf{r}_{i+1} &= \mathbf{w}_h^{n+1} + \Delta t \, a_i \mathbf{k}_i,\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^{n+1} + \Delta t \, b_i \mathbf{k}_i.
\end{aligned}


@f]

除了连续更新的向量 $\mathbf w_h^{n+1}$ ，这个方案只需要两个辅助向量，即保存微分算子的评估的向量 $\mathbf{k}_i$ ，以及保存微分算子应用的右手边的向量 $\mathbf{r}_i$ 。在后续阶段  $i$  ，值  $\mathbf{k}_i$  和  $\mathbf{r}_i$  可以使用相同的存储。

低存储量变体的主要优点是一方面减少了内存消耗（如果必须在内存中装入非常多的未知数，持有所有的 $\mathbf{k}_i$ 来计算随后的更新，对于 $s$ 来说已经是一个极限，在5到8之间--记得我们使用的是显式方案，所以我们不需要存储任何通常比几个向量大很多的矩阵），另一方面是减少内存访问。在这个程序中，我们对后一个方面特别感兴趣。由于运算符评估的成本只是简单地从内存中流转输入和输出向量的一小部分，我们必须考虑向量更新的成本，而低存储的变体可以提供传统显式Runge--Kutta方法两倍的吞吐量，原因就在于此，例如，见 @cite SchoederKormann2018 中的分析。

除了参考文献 @cite KennedyCarpenterLewis2000 中的三阶、四阶和五阶精度的三个变体外，我们还使用了一个四阶精度的七级变体，该变体是为声学设置而优化的 @cite TseliosSimos2007  。声学问题是欧拉方程的亚音速制度的有趣方面之一，其中可压缩性导致了声波的传播；通常，人们使用围绕背景状态的线性化欧拉方程的进一步简化，或围绕固定框架的声波方程。




<a name="Fastevaluationofintegralsbymatrixfreetechniques"></a><h3>Fast evaluation of integrals by matrix-free techniques</h3>


这个程序中使用的主要成分是我们用来评估算子  $\mathcal L_h$  和反质量矩阵  $\mathcal M$  的快速无矩阵技术。实际上，<i>matrix-free</i>这个术语有点名不副实，因为我们是在处理一个非线性算子，并没有将反过来可以用矩阵表示的算子线性化。然而，作为稀疏矩阵-向量乘积的替代品，积分的快速评估已经变得很流行，如步骤-37和步骤-59所示，为此我们在交易二中创造了这个基础设施<i>matrix-free functionality</i>。此外，反质量矩阵确实是以无矩阵的方式应用的，详见下文。

无矩阵基础设施使我们能够快速评估弱形式的积分。其成分是将解系数快速插值为正交点的值和导数，在正交点进行逐点运算（在这里我们实现了上述的微分算子），以及与所有测试函数相乘和对正交点求和。第一和第三部分利用了和因子化，并在步骤37的单元积分教程和步骤59的面积分教程中进行了广泛的讨论。唯一的区别是，我们现在处理的是一个 $d+2$ 分量的系统，而不是以前教程程序中的标量系统。在代码中，所有的变化是FEEvaluation和FEFaceEvaluation类的一个模板参数，即设置分量的数量。对向量的访问和以前一样，都由评价器透明地处理。我们还注意到，下面的代码中选择的带有单一评价器的变体并不是唯一的选择--我们也可以为单独的组件 $\rho$ 、 $\rho \mathbf u$ 和 $E$ 使用单独的评价器；鉴于我们对所有组件的处理是类似的（也反映在我们把方程作为一个矢量系统的方式），这里会更复杂。和以前一样，FEEvaluation类通过结合对几个单元（和面）的操作来提供显式的矢量化，涉及的数据类型称为VectorizedArray。由于这种类型的算术运算都是重载的，所以我们不必为它费心，除了通过函数接口对函数进行评估，我们需要同时为几个正交点的位置提供特殊的<i>vectorized</i>评估。

这个程序中更大的变化是在正交点的操作。在这里，多分量评估器为我们提供了之前没有讨论过的返回类型。 FEEvaluation::get_value() 将为第37步的拉普拉斯返回一个标量（更准确地说，由于跨单元的矢量化，是一个VectorizedArray类型），现在它返回的类型是`Tensor<1,dim+2,VectorizedArray<Number>'。同样，梯度类型现在是`张量<1,dim+2,张量<1,dim,矢量化数组<Number>>`，其中外部张量收集了欧拉系统的`dim+2'分量，内部张量是各个方向的偏导数。例如，欧拉系统的通量 $\mathbf{F}(\mathbf{w})$ 就属于这种类型。为了减少我们为拼出这些类型而写的代码量，我们尽可能使用C++的`自动'关键字。

从实施的角度来看，非线性并不是一个很大的困难。它是在我们表达欧拉弱形式的条款时自然引入的，例如以动量条款的形式  $\rho \mathbf{u}
\otimes \mathbf{u}$  。为了得到这个表达式，我们首先从动量变量  $\rho \mathbf{u}$  推导出速度  $\mathbf{u}$  。鉴于 $\rho
\mathbf{u}$ 和 $\rho$ 一样被表示为 $p$ 度的多项式，速度 $\mathbf{u}$ 是参考坐标 $\hat{\mathbf{x}}$ 的一个有理表达。当我们进行乘法 $(\rho
\mathbf{u})\otimes \mathbf{u}$ 时，我们得到一个表达式，它是两个多项式的比值，分子中的多项式程度 $2p$ 和分母中的多项式程度 $p$ 。结合测试函数的梯度，分子中的积分度为 $3p$ ，分母中的积分度为 $p$ ，对于仿生单元，即平行四边形/平行四边形，已经有了积分。对于弧形单元，当积分乘以映射的雅各布系数时，会出现额外的多项式和有理表达式。在这一点上，人们通常需要放弃坚持精确的积分，而采取高斯（更确切地说，高斯-勒格伦德）正交提供的任何精度。这时的情况与拉普拉斯方程的情况类似，积分项包含非affince单元上的有理表达式，也只能进行近似积分。由于这些公式只对多项式进行精确积分，我们不得不以积分错误的形式忍受<a
href="https://mathoverflow.net/questions/26018/what-are-variational-crimes-and-who-coined-the-term">variational
crime</a>的影响。

虽然对于椭圆问题来说，不精确的积分通常是可以容忍的，但对于双曲问题来说，不精确的积分会引起一些令人头痛的效应，这种效应称为<b>aliasing</b>。这个术语来自于信号处理，表达了不适当的、过于粗糙的采样情况。就正交而言，不适当的采样意味着我们使用的正交点与准确采样变系数积分所需的点相比太少。在DG文献中已经表明，别离误差会在<i>barely</i>解析模拟的数值解中引入非物理性的振荡。别名主要影响到粗略的分辨率--而采用相同方案的更细的网格则工作良好--这一事实并不令人惊讶，因为分辨率高的模拟往往在一个单元的长度尺度上是平滑的（即，它们在较高的多项式程度上有小的系数，由于正交点太少而被遗漏，而在较低的多项式程度上的主要解贡献仍然被很好地捕获--这只是泰勒定理的一个结果）。为了解决这个问题，DG文献中提出了各种方法。一种技术是过滤，它可以抑制与高次多项式度数有关的解成分。由于所选择的节点基不是分层的，这就意味着要从节点基转化为分层基（例如，基于Legendre多项式的模态基），其中单元内的贡献是按多项式程度划分的。在这个基础上，我们可以将与高度数相关的求解系数乘以一个小数，保持低度数不变（以避免破坏一致性），然后再转换回节点基础。然而，过滤器会降低该方法的准确性。另一个在某种意义上更简单的策略是使用更多的正交点来更准确地捕捉非线性项。每个坐标方向使用超过 $p+1$ 个正交点有时被称为过度积分或一致积分。后者在不可压缩的Navier-Stokes方程中最为常见，其中 $\mathbf{u}\otimes \mathbf{u}$ 非线性导致 $3p$ 度的多项式积分（当同时考虑测试函数时），只要元素的几何形状是仿生的，每个方向的 $\textrm{floor}\left(\frac{3p}{2}\right)+1$ 正交点就可以精确积分。在非多项式积分的欧拉方程的背景下，选择就不那么明确了。根据各种变量的变化， $\textrm{floor}\left(\frac{3p}{2}\right)+1$ 或 $2p+1$ 点（分别精确积分度为 $3p$ 或 $4p$ 的多项式）都很常见。

为了反映程序中正交选择的这种可变性，我们把正交点的数量作为一个变量来指定，就像多项式的度数一样，并注意到人们会根据流量配置做出不同的选择。默认选择是 $p+2$ 点--比最小可能的 $p+1$ 点多一点。FEEvaluation和FEFaceEvaluation类允许通过模板参数无缝地改变点的数量，这样程序就不会因此而变得更复杂。




<a name="Evaluationoftheinversemassmatrixwithmatrixfreetechniques"></a><h3>Evaluation of the inverse mass matrix with matrix-free techniques</h3>


最后一个要素是反质量矩阵的评估  $\mathcal
M^{-1}$  。在具有显式时间积分的DG方法中，质量矩阵是块状对角线，因此很容易反转--人们只需要反转对角线块。然而，考虑到无矩阵的积分评估在成本上更接近于只访问向量，即使应用块对角矩阵（例如通过LU因子数组）也会比评估 $\mathcal L_h$ 贵几倍，仅仅是因为对于高阶有限元来说，仅仅存储和加载大小为`dofs_per_cell`x`dofs_per_cell`的矩阵是昂贵的。由于这显然是不可取的，部分社区已经转移到质量矩阵是对角线的基础，例如<i>L<sub>2</sub></i>正交Legendre基础，使用分层多项式或高斯四分法点上的拉格朗日多项式（这只是利用Legendre信息的另一种方式）。虽然对角线属性对于变形元素来说是失效的，但通过采取对角线质量矩阵而忽略其余部分（质量包络的变种，尽管不是步骤-48中利用的具有额外积分误差的变种）所产生的误差已被证明不会改变离散化精度。高斯正交点中的拉格朗日基础有时也被称为同位设置，因为多项式的结点与正交点重合（="同位"），避免了一些内插操作。鉴于我们想在 $\mathcal L_h$ 中对非线性项使用更多的正交点，然而，拼合属性就失去了。(更确切地说，在改变基础后，它仍然用于FEEvaluation和FEFaceEvaluation，见无矩阵论文  @cite KronbichlerKormann2019  。)

在这个教程程序中，我们使用拼合思想来应用反质量矩阵，但有一个小的转折。与其在高斯四分法的点上通过拉格朗日多项式使用配位，我们更倾向于在高斯-洛巴托点上使用传统的拉格朗日基础，因为那些使面积分的评估变得便宜。这是因为对于高斯-洛巴托点来说，一些节点点位于单元格的面上，而且不难证明，在任何给定的面上，唯一具有非零值的形状函数正是其节点点实际上位于该面上的那些。当然，我们也可以像步骤48那样使用高斯-洛巴托正交（有一些额外的积分误差），但我们不想牺牲精度，因为这些正交公式通常比一般的高斯正交公式的阶数低。相反，我们使用参考文献 @cite KronbichlerSchoeder2016 中描述的一个想法，其中提出为了应用反质量矩阵而改变基础。让我们用 $S$ 表示在正交点评价的形状函数矩阵，形状函数在矩阵的行中，正交点在列中。那么，单元格 $K$ 上的质量矩阵由以下公式给出

@f[
\mathcal M^K = S J^K S^\mathrm T.


@f]

这里， $J^K$ 是以雅各布系数乘以正交权重（JxW）的行列式作为条目的对角矩阵。矩阵 $S$ 被构造为一维矩阵的克朗克积（张量积），例如，在三维中为

@f[
S = S_{\text{1D}}\otimes S_{\text{1D}}\otimes S_{\text{1D}},


@f]

这是基函数是一维形状函数的张量积，正交公式是一维正交公式的张量积的结果。对于多项式的数量等于正交点的数量的情况， $S J^K S^\mathrm T$ 中的所有矩阵都是方形的，同样，克朗克积中的 $S$ 的成分也是方形的。因此，人们可以对每个矩阵进行反转，形成整体的逆。

@f[
\left(\mathcal M^K\right)^{-1} = S_{\text{1D}}^{-\mathrm T}\otimes
S_{\text{1D}}^{-\mathrm T}\otimes S_{\text{1D}}^{-\mathrm T}
\left(J^K\right)^{-1}
S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}.


@f]

这个公式的结构与用和因子化技术对积分进行正向评价的步骤完全相同（即交易.II的FEEvaluation和MatrixFree框架）。因此，我们可以利用相同的代码路径，采用不同的插值矩阵， $S_{\mathrm{1D}}^{-\mathrm{T}}$ 而不是 $S_{\mathrm{1D}}$  。

类 MatrixFreeOperators::CellwiseInverseMassMatrix 实现了这个操作。它从有限元中包含的基（在这里是FE_DGQ）改变为高斯正交点中的拉格朗日基。在这里，可以评估对角线质量矩阵的逆值，这只是`JxW`因子的逆值（即正交权重乘以从参考坐标到实坐标的雅各布系数）。一旦这样做了，我们就可以变回标准的节点高斯-洛巴托基础。

这种应用反质量矩阵的特殊方式的优点是成本类似于质量矩阵的正向应用，这比用超积分和面积分评估空间算子 $\mathcal L_h$ 更便宜。(我们将在<a href="#Results">results section</a>中用详细的时间信息证明这一点)。事实上，它是如此便宜，以至于在大多数现代架构上，它被读取源向量、读取对角线和写入目的向量的带宽所限制。用于结果部分的硬件可以使计算的速度至少比从内存流向量的速度快一倍。




<a name="Thetestcase"></a><h3>The test case</h3>


在这个教程程序中，我们实现了两个测试案例。第一个案例是限于两个空间维度的收敛性测试。它运行一个所谓的等熵涡旋，它通过一个背景流场进行传输。第二个案例使用了一个更令人兴奋的设置。我们从一个浸在通道中的圆柱体开始，使用 GridGenerator::channel_with_cylinder() 函数。在这里，我们强加一个马赫数为 $\mathrm{Ma}=0.307$ 的亚音速初始场，在 $x$ 方向上速度不变。在顶壁和底壁以及圆柱体上，我们施加了一个无穿透（即切向流动）的条件。与初始条件相比，这种设置迫使气流重新定向，从而导致大的声波从圆柱体上传播出去。在上游方向，波的传播速度较慢（因为它必须逆着迎面而来的气体移动），包括密度和压力的不连续。在下游方向，由于声音的传播和流体的流动方向相同，传输速度较快，这在一定程度上抹去了不连续性。一旦声波碰到上下壁，声音就会被反射回来，形成一些漂亮的形状，如下图<a href="#Results">results section</a>所示。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The include files are similar to the previous matrix-free tutorial programs
 * step-37, step-48, and step-59
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/time_stepping.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/vectorization.h>
 * 
 * #include <deal.II/distributed/tria.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * 
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_system.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/tria.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * 
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * #include <deal.II/matrix_free/matrix_free.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * 
 * #include <fstream>
 * #include <iomanip>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The following file includes the CellwiseInverseMassMatrix data structure
 * that we will use for the mass matrix inversion, the only new include
 * file for this tutorial program:
 * 
 * @code
 * #include <deal.II/matrix_free/operators.h>
 * 
 * 
 * 
 * namespace Euler_DG
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Similarly to the other matrix-free tutorial programs, we collect all
 * parameters that control the execution of the program at the top of the
 * file. Besides the dimension and polynomial degree we want to run with, we
 * also specify a number of points in the Gaussian quadrature formula we
 * want to use for the nonlinear terms in the Euler equations. Furthermore,
 * we specify the time interval for the time-dependent problem, and
 * implement two different test cases. The first one is an analytical
 * solution in 2D, whereas the second is a channel flow around a cylinder as
 * described in the introduction. Depending on the test case, we also change
 * the final time up to which we run the simulation, and a variable
 * `output_tick` that specifies in which intervals we want to write output
 * (assuming that the tick is larger than the time step size).
 * 
 * @code
 *   constexpr unsigned int testcase             = 0;
 *   constexpr unsigned int dimension            = 2;
 *   constexpr unsigned int n_global_refinements = 3;
 *   constexpr unsigned int fe_degree            = 5;
 *   constexpr unsigned int n_q_points_1d        = fe_degree + 2;
 * 
 *   using Number = double;
 * 
 *   constexpr double gamma       = 1.4;
 *   constexpr double final_time  = testcase == 0 ? 10 : 2.0;
 *   constexpr double output_tick = testcase == 0 ? 1 : 0.05;
 * 
 * @endcode
 * 
 * Next off are some details of the time integrator, namely a Courant number
 * that scales the time step size in terms of the formula $\Delta t =
 * \text{Cr} n_\text{stages} \frac{h}{(p+1)^{1.5} (\|\mathbf{u} +
 * c)_\text{max}}$, as well as a selection of a few low-storage Runge--Kutta
 * methods. We specify the Courant number per stage of the Runge--Kutta
 * scheme, as this gives a more realistic expression of the numerical cost
 * for schemes of various numbers of stages.
 * 
 * @code
 *   const double courant_number = 0.15 / std::pow(fe_degree, 1.5);
 *   enum LowStorageRungeKuttaScheme
 *   {
 *     stage_3_order_3, /* Kennedy, Carpenter, Lewis, 2000 */
 *     stage_5_order_4, /* Kennedy, Carpenter, Lewis, 2000 */
 *     stage_7_order_4, /* Tselios, Simos, 2007 */
 *     stage_9_order_5, /* Kennedy, Carpenter, Lewis, 2000 */
 *   };
 *   constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4;
 * 
 * @endcode
 * 
 * Eventually, we select a detail of the spatial discretization, namely the
 * numerical flux (Riemann solver) at the faces between cells. For this
 * program, we have implemented a modified variant of the Lax--Friedrichs
 * flux and the Harten--Lax--van Leer (HLL) flux.
 * 
 * @code
 *   enum EulerNumericalFlux
 *   {
 *     lax_friedrichs_modified,
 *     harten_lax_vanleer,
 *   };
 *   constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified;
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * We now define a class with the exact solution for the test case 0 and one
 * with a background flow field for test case 1 of the channel. Given that
 * the Euler equations are a problem with $d+2$ equations in $d$ dimensions,
 * we need to tell the Function base class about the correct number of
 * components.
 * 
 * @code
 *   template <int dim>
 *   class ExactSolution : public Function<dim>
 *   {
 *   public:
 *     ExactSolution(const double time)
 *       : Function<dim>(dim + 2, time)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * As far as the actual function implemented is concerned, the analytical
 * test case is an isentropic vortex case (see e.g. the book by Hesthaven
 * and Warburton, Example 6.1 in Section 6.6 on page 209) which fulfills the
 * Euler equations with zero force term on the right hand side. Given that
 * definition, we return either the density, the momentum, or the energy
 * depending on which component is requested. Note that the original
 * definition of the density involves the $\frac{1}{\gamma -1}$-th power of
 * some expression. Since `std::pow()` has pretty slow implementations on
 * some systems, we replace it by logarithm followed by exponentiation (of
 * base 2), which is mathematically equivalent but usually much better
 * optimized. This formula might lose accuracy in the last digits
 * for very small numbers compared to `std::pow()`, but we are happy with
 * it anyway, since small numbers map to data close to 1.
 *   

 * 
 * For the channel test case, we simply select a density of 1, a velocity of
 * 0.4 in $x$ direction and zero in the other directions, and an energy that
 * corresponds to a speed of sound of 1.3 measured against the background
 * velocity field, computed from the relation $E = \frac{c^2}{\gamma (\gamma
 * -1)} + \frac 12 \rho \|u\|^2$.
 * 
 * @code
 *   template <int dim>
 *   double ExactSolution<dim>::value(const Point<dim> & x,
 *                                    const unsigned int component) const
 *   {
 *     const double t = this->get_time();
 * 
 *     switch (testcase)
 *       {
 *         case 0:
 *           {
 *             Assert(dim == 2, ExcNotImplemented());
 *             const double beta = 5;
 * 
 *             Point<dim> x0;
 *             x0[0] = 5.;
 *             const double radius_sqr =
 *               (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t;
 *             const double factor =
 *               beta / (numbers::PI * 2) * std::exp(1. - radius_sqr);
 *             const double density_log = std::log2(
 *               std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor));
 *             const double density = std::exp2(density_log * (1. / (gamma - 1.)));
 *             const double u       = 1. - factor * (x[1] - x0[1]);
 *             const double v       = factor * (x[0] - t - x0[0]);
 * 
 *             if (component == 0)
 *               return density;
 *             else if (component == 1)
 *               return density * u;
 *             else if (component == 2)
 *               return density * v;
 *             else
 *               {
 *                 const double pressure =
 *                   std::exp2(density_log * (gamma / (gamma - 1.)));
 *                 return pressure / (gamma - 1.) +
 *                        0.5 * (density * u * u + density * v * v);
 *               }
 *           }
 * 
 *         case 1:
 *           {
 *             if (component == 0)
 *               return 1.;
 *             else if (component == 1)
 *               return 0.4;
 *             else if (component == dim + 1)
 *               return 3.097857142857143;
 *             else
 *               return 0.;
 *           }
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *           return 0.;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LowstorageexplicitRungeKuttatimeintegrators"></a> 
 * <h3>Low-storage explicit Runge--Kutta time integrators</h3>
 * 

 * 
 * The next few lines implement a few low-storage variants of Runge--Kutta
 * methods. These methods have specific Butcher tableaux with coefficients
 * $b_i$ and $a_i$ as shown in the introduction. As usual in Runge--Kutta
 * method, we can deduce time steps, $c_i = \sum_{j=1}^{i-2} b_i + a_{i-1}$
 * from those coefficients. The main advantage of this kind of scheme is the
 * fact that only two vectors are needed per stage, namely the accumulated
 * part of the solution $\mathbf{w}$ (that will hold the solution
 * $\mathbf{w}^{n+1}$ at the new time $t^{n+1}$ after the last stage), the
 * update vector $\mathbf{r}_i$ that gets evaluated during the stages, plus
 * one vector $\mathbf{k}_i$ to hold the evaluation of the operator. Such a
 * Runge--Kutta setup reduces the memory storage and memory access. As the
 * memory bandwidth is often the performance-limiting factor on modern
 * hardware when the evaluation of the differential operator is
 * well-optimized, performance can be improved over standard time
 * integrators. This is true also when taking into account that a
 * conventional Runge--Kutta scheme might allow for slightly larger time
 * steps as more free parameters allow for better stability properties.
 *   

 * 
 * In this tutorial programs, we concentrate on a few variants of
 * low-storage schemes defined in the article by Kennedy, Carpenter, and
 * Lewis (2000), as well as one variant described by Tselios and Simos
 * (2007). There is a large series of other schemes available, which could
 * be addressed by additional sets of coefficients or slightly different
 * update formulas.
 *   

 * 
 * We define a single class for the four integrators, distinguished by the
 * enum described above. To each scheme, we then fill the vectors for the
 * $b_i$ and $a_i$ to the given variables in the class.
 * 
 * @code
 *   class LowStorageRungeKuttaIntegrator
 *   {
 *   public:
 *     LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme)
 *     {
 *       TimeStepping::runge_kutta_method lsrk;
 * @endcode
 * 
 * First comes the three-stage scheme of order three by Kennedy et al.
 * (2000). While its stability region is significantly smaller than for
 * the other schemes, it only involves three stages, so it is very
 * competitive in terms of the work per stage.
 * 
 * @code
 *       switch (scheme)
 *         {
 *           case stage_3_order_3:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The next scheme is a five-stage scheme of order four, again
 * defined in the paper by Kennedy et al. (2000).
 * 
 * @code
 *           case stage_5_order_4:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The following scheme of seven stages and order four has been
 * explicitly derived for acoustics problems. It is a balance of
 * accuracy for imaginary eigenvalues among fourth order schemes,
 * combined with a large stability region. Since DG schemes are
 * dissipative among the highest frequencies, this does not
 * necessarily translate to the highest possible time step per
 * stage. In the context of the present tutorial program, the
 * numerical flux plays a crucial role in the dissipation and thus
 * also the maximal stable time step size. For the modified
 * Lax--Friedrichs flux, this scheme is similar to the
 * `stage_5_order_4` scheme in terms of step size per stage if only
 * stability is considered, but somewhat less efficient for the HLL
 * flux.
 * 
 * @code
 *           case stage_7_order_4:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4;
 *               break;
 *             }
 * 
 * @endcode
 * 
 * The last scheme included here is the nine-stage scheme of order
 * five from Kennedy et al. (2000). It is the most accurate among
 * the schemes used here, but the higher order of accuracy
 * sacrifices some stability, so the step length normalized per
 * stage is less than for the fourth order schemes.
 * 
 * @code
 *           case stage_9_order_5:
 *             {
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5;
 *               break;
 *             }
 * 
 *           default:
 *             AssertThrow(false, ExcNotImplemented());
 *         }
 *       TimeStepping::LowStorageRungeKutta<
 *         LinearAlgebra::distributed::Vector<Number>>
 *         rk_integrator(lsrk);
 *       rk_integrator.get_coefficients(ai, bi, ci);
 *     }
 * 
 *     unsigned int n_stages() const
 *     {
 *       return bi.size();
 *     }
 * 
 * @endcode
 * 
 * The main function of the time integrator is to go through the stages,
 * evaluate the operator, prepare the $\mathbf{r}_i$ vector for the next
 * evaluation, and update the solution vector $\mathbf{w}$. We hand off
 * the work to the `pde_operator` involved in order to be able to merge
 * the vector operations of the Runge--Kutta setup with the evaluation of
 * the differential operator for better performance, so all we do here is
 * to delegate the vectors and coefficients.
 *     

 * 
 * We separately call the operator for the first stage because we need
 * slightly modified arguments there: We evaluate the solution from
 * the old solution $\mathbf{w}^n$ rather than a $\mathbf r_i$ vector, so
 * the first argument is `solution`. We here let the stage vector
 * $\mathbf{r}_i$ also hold the temporary result of the evaluation, as it
 * is not used otherwise. For all subsequent stages, we use the vector
 * `vec_ki` as the second vector argument to store the result of the
 * operator evaluation. Finally, when we are at the last stage, we must
 * skip the computation of the vector $\mathbf{r}_{s+1}$ as there is no
 * coefficient $a_s$ available (nor will it be used).
 * 
 * @code
 *     template <typename VectorType, typename Operator>
 *     void perform_time_step(const Operator &pde_operator,
 *                            const double    current_time,
 *                            const double    time_step,
 *                            VectorType &    solution,
 *                            VectorType &    vec_ri,
 *                            VectorType &    vec_ki) const
 *     {
 *       AssertDimension(ai.size() + 1, bi.size());
 * 
 *       pde_operator.perform_stage(current_time,
 *                                  bi[0] * time_step,
 *                                  ai[0] * time_step,
 *                                  solution,
 *                                  vec_ri,
 *                                  solution,
 *                                  vec_ri);
 * 
 *       for (unsigned int stage = 1; stage < bi.size(); ++stage)
 *         {
 *           const double c_i = ci[stage];
 *           pde_operator.perform_stage(current_time + c_i * time_step,
 *                                      bi[stage] * time_step,
 *                                      (stage == bi.size() - 1 ?
 *                                         0 :
 *                                         ai[stage] * time_step),
 *                                      vec_ri,
 *                                      vec_ki,
 *                                      solution,
 *                                      vec_ri);
 *         }
 *     }
 * 
 *   private:
 *     std::vector<double> bi;
 *     std::vector<double> ai;
 *     std::vector<double> ci;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofpointwiseoperationsoftheEulerequations"></a> 
 * <h3>Implementation of point-wise operations of the Euler equations</h3>
 * 

 * 
 * In the following functions, we implement the various problem-specific
 * operators pertaining to the Euler equations. Each function acts on the
 * vector of conserved variables $[\rho, \rho\mathbf{u}, E]$ that we hold in
 * the solution vectors, and computes various derived quantities.
 *   

 * 
 * First out is the computation of the velocity, that we derive from the
 * momentum variable $\rho \mathbf{u}$ by division by $\rho$. One thing to
 * note here is that we decorate all those functions with the keyword
 * `DEAL_II_ALWAYS_INLINE`. This is a special macro that maps to a
 * compiler-specific keyword that tells the compiler to never create a
 * function call for any of those functions, and instead move the
 * implementation <a
 * href="https://en.wikipedia.org/wiki/Inline_function">inline</a> to where
 * they are called. This is critical for performance because we call into some
 * of those functions millions or billions of times: For example, we both use
 * the velocity for the computation of the flux further down, but also for the
 * computation of the pressure, and both of these places are evaluated at
 * every quadrature point of every cell. Making sure these functions are
 * inlined ensures not only that the processor does not have to execute a jump
 * instruction into the function (and the corresponding return jump), but also
 * that the compiler can re-use intermediate information from one function's
 * context in code that comes after the place where the function was called.
 * (We note that compilers are generally quite good at figuring out which
 * functions to inline by themselves. Here is a place where compilers may or
 * may not have figured it out by themselves but where we know for sure that
 * inlining is a win.)
 *   

 * 
 * Another trick we apply is a separate variable for the inverse density
 * $\frac{1}{\rho}$. This enables the compiler to only perform a single
 * division for the flux, despite the division being used at several
 * places. As divisions are around ten to twenty times as expensive as
 * multiplications or additions, avoiding redundant divisions is crucial for
 * performance. We note that taking the inverse first and later multiplying
 * with it is not equivalent to a division in floating point arithmetic due
 * to roundoff effects, so the compiler is not allowed to exchange one way by
 * the other with standard optimization flags. However, it is also not
 * particularly difficult to write the code in the right way.
 *   

 * 
 * To summarize, the chosen strategy of always inlining and careful
 * definition of expensive arithmetic operations allows us to write compact
 * code without passing all intermediate results around, despite making sure
 * that the code maps to excellent machine code.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim, Number>
 *     euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Number inverse_density = Number(1.) / conserved_variables[0];
 * 
 *     Tensor<1, dim, Number> velocity;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       velocity[d] = conserved_variables[1 + d] * inverse_density;
 * 
 *     return velocity;
 *   }
 * 
 * @endcode
 * 
 * The next function computes the pressure from the vector of conserved
 * variables, using the formula $p = (\gamma - 1) \left(E - \frac 12 \rho
 * \mathbf{u}\cdot \mathbf{u}\right)$. As explained above, we use the
 * velocity from the `euler_velocity()` function. Note that we need to
 * specify the first template argument `dim` here because the compiler is
 * not able to deduce it from the arguments of the tensor, whereas the
 * second argument (number type) can be automatically deduced.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Number
 *     euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Tensor<1, dim, Number> velocity =
 *       euler_velocity<dim>(conserved_variables);
 * 
 *     Number rho_u_dot_u = conserved_variables[1] * velocity[0];
 *     for (unsigned int d = 1; d < dim; ++d)
 *       rho_u_dot_u += conserved_variables[1 + d] * velocity[d];
 * 
 *     return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u);
 *   }
 * 
 * @endcode
 * 
 * Here is the definition of the Euler flux function, i.e., the definition
 * of the actual equation. Given the velocity and pressure (that the
 * compiler optimization will make sure are done only once), this is
 * straight-forward given the equation stated in the introduction.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>>
 *     euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
 *   {
 *     const Tensor<1, dim, Number> velocity =
 *       euler_velocity<dim>(conserved_variables);
 *     const Number pressure = euler_pressure<dim>(conserved_variables);
 * 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       {
 *         flux[0][d] = conserved_variables[1 + d];
 *         for (unsigned int e = 0; e < dim; ++e)
 *           flux[e + 1][d] = conserved_variables[e + 1] * velocity[d];
 *         flux[d + 1][d] += pressure;
 *         flux[dim + 1][d] =
 *           velocity[d] * (conserved_variables[dim + 1] + pressure);
 *       }
 * 
 *     return flux;
 *   }
 * 
 * @endcode
 * 
 * This next function is a helper to simplify the implementation of the
 * numerical flux, implementing the action of a tensor of tensors (with
 * non-standard outer dimension of size `dim + 2`, so the standard overloads
 * provided by deal.II's tensor classes do not apply here) with another
 * tensor of the same inner dimension, i.e., a matrix-vector product.
 * 
 * @code
 *   template <int n_components, int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, n_components, Number>
 *     operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
 *               const Tensor<1, dim, Number> &                         vector)
 *   {
 *     Tensor<1, n_components, Number> result;
 *     for (unsigned int d = 0; d < n_components; ++d)
 *       result[d] = matrix[d] * vector;
 *     return result;
 *   }
 * 
 * @endcode
 * 
 * This function implements the numerical flux (Riemann solver). It gets the
 * state from the two sides of an interface and the normal vector, oriented
 * from the side of the solution $\mathbf{w}^-$ towards the solution
 * $\mathbf{w}^+$. In finite volume methods which rely on piece-wise
 * constant data, the numerical flux is the central ingredient as it is the
 * only place where the physical information is entered. In DG methods, the
 * numerical flux is less central due to the polynomials within the elements
 * and the physical flux used there. As a result of higher-degree
 * interpolation with consistent values from both sides in the limit of a
 * continuous solution, the numerical flux can be seen as a control of the
 * jump of the solution from both sides to weakly impose continuity. It is
 * important to realize that a numerical flux alone cannot stabilize a
 * high-order DG method in the presence of shocks, and thus any DG method
 * must be combined with further shock-capturing techniques to handle those
 * cases. In this tutorial, we focus on wave-like solutions of the Euler
 * equations in the subsonic regime without strong discontinuities where our
 * basic scheme is sufficient.
 *   

 * 
 * Nonetheless, the numerical flux is decisive in terms of the numerical
 * dissipation of the overall scheme and influences the admissible time step
 * size with explicit Runge--Kutta methods. We consider two choices, a
 * modified Lax--Friedrichs scheme and the widely used Harten--Lax--van Leer
 * (HLL) flux. For both variants, we first need to get the velocities and
 * pressures from both sides of the interface and evaluate the physical
 * Euler flux.
 *   

 * 
 * For the local Lax--Friedrichs flux, the definition is $\hat{\mathbf{F}}
 * =\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
 * \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
 * \mathbf{n^-}$, where the factor $\lambda =
 * \max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ gives the
 * maximal wave speed and $c = \sqrt{\gamma p / \rho}$ is the speed of
 * sound. Here, we choose two modifications of that expression for reasons
 * of computational efficiency, given the small impact of the flux on the
 * solution. For the above definition of the factor $\lambda$, we would need
 * to take four square roots, two for the two velocity norms and two for the
 * speed of sound on either side. The first modification is hence to rather
 * use $\sqrt{\|\mathbf{u}\|^2+c^2}$ as an estimate of the maximal speed
 * (which is at most a factor of 2 away from the actual maximum, as shown in
 * the introduction). This allows us to pull the square root out of the
 * maximum and get away with a single square root computation. The second
 * modification is to further relax on the parameter $\lambda$---the smaller
 * it is, the smaller the dissipation factor (which is multiplied by the
 * jump in $\mathbf{w}$, which might result in a smaller or bigger
 * dissipation in the end). This allows us to fit the spectrum into the
 * stability region of the explicit Runge--Kutta integrator with bigger time
 * steps. However, we cannot make dissipation too small because otherwise
 * imaginary eigenvalues grow larger. Finally, the current conservative
 * formulation is not energy-stable in the limit of $\lambda\to 0$ as it is
 * not skew-symmetric, and would need additional measures such as split-form
 * DG schemes in that case.
 *   

 * 
 * For the HLL flux, we follow the formula from literature, introducing an
 * additional weighting of the two states from Lax--Friedrichs by a
 * parameter $s$. It is derived from the physical transport directions of
 * the Euler equations in terms of the current direction of velocity and
 * sound speed. For the velocity, we here choose a simple arithmetic average
 * which is sufficient for DG scenarios and moderate jumps in material
 * parameters.
 *   

 * 
 * Since the numerical flux is multiplied by the normal vector in the weak
 * form, we multiply by the result by the normal vector for all terms in the
 * equation. In these multiplications, the `operator*` defined above enables
 * a compact notation similar to the mathematical definition.
 *   

 * 
 * In this and the following functions, we use variable suffixes `_m` and
 * `_p` to indicate quantities derived from $\mathbf{w}^-$ and $\mathbf{w}^+$,
 * i.e., values "here" and "there" relative to the current cell when looking
 * at a neighbor cell.
 * 
 * @code
 *   template <int dim, typename Number>
 *   inline DEAL_II_ALWAYS_INLINE 
 *     Tensor<1, dim + 2, Number>
 *     euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
 *                          const Tensor<1, dim + 2, Number> &u_p,
 *                          const Tensor<1, dim, Number> &    normal)
 *   {
 *     const auto velocity_m = euler_velocity<dim>(u_m);
 *     const auto velocity_p = euler_velocity<dim>(u_p);
 * 
 *     const auto pressure_m = euler_pressure<dim>(u_m);
 *     const auto pressure_p = euler_pressure<dim>(u_p);
 * 
 *     const auto flux_m = euler_flux<dim>(u_m);
 *     const auto flux_p = euler_flux<dim>(u_p);
 * 
 *     switch (numerical_flux_type)
 *       {
 *         case lax_friedrichs_modified:
 *           {
 *             const auto lambda =
 *               0.5 * std::sqrt(std::max(velocity_p.norm_square() +
 *                                          gamma * pressure_p * (1. / u_p[0]),
 *                                        velocity_m.norm_square() +
 *                                          gamma * pressure_m * (1. / u_m[0])));
 * 
 *             return 0.5 * (flux_m * normal + flux_p * normal) +
 *                    0.5 * lambda * (u_m - u_p);
 *           }
 * 
 *         case harten_lax_vanleer:
 *           {
 *             const auto avg_velocity_normal =
 *               0.5 * ((velocity_m + velocity_p) * normal);
 *             const auto   avg_c = std::sqrt(std::abs(
 *               0.5 * gamma *
 *               (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0]))));
 *             const Number s_pos =
 *               std::max(Number(), avg_velocity_normal + avg_c);
 *             const Number s_neg =
 *               std::min(Number(), avg_velocity_normal - avg_c);
 *             const Number inverse_s = Number(1.) / (s_pos - s_neg);
 * 
 *             return inverse_s *
 *                    ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) -
 *                     s_pos * s_neg * (u_m - u_p));
 *           }
 * 
 *         default:
 *           {
 *             Assert(false, ExcNotImplemented());
 *             return {};
 *           }
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This and the next function are helper functions to provide compact
 * evaluation calls as multiple points get batched together via a
 * VectorizedArray argument (see the step-37 tutorial for details). This
 * function is used for the subsonic outflow boundary conditions where we
 * need to set the energy component to a prescribed value. The next one
 * requests the solution on all components and is used for inflow boundaries
 * where all components of the solution are set.
 * 
 * @code
 *   template <int dim, typename Number>
 *   VectorizedArray<Number>
 *   evaluate_function(const Function<dim> &                      function,
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized,
 *                     const unsigned int                         component)
 *   {
 *     VectorizedArray<Number> result;
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *       {
 *         Point<dim> p;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           p[d] = p_vectorized[d][v];
 *         result[v] = function.value(p, component);
 *       }
 *     return result;
 *   }
 * 
 * 
 *   template <int dim, typename Number, int n_components = dim + 2>
 *   Tensor<1, n_components, VectorizedArray<Number>>
 *   evaluate_function(const Function<dim> &                      function,
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized)
 *   {
 *     AssertDimension(function.n_components, n_components);
 *     Tensor<1, n_components, VectorizedArray<Number>> result;
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *       {
 *         Point<dim> p;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           p[d] = p_vectorized[d][v];
 *         for (unsigned int d = 0; d < n_components; ++d)
 *           result[d][v] = function.value(p, d);
 *       }
 *     return result;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerOperationclass"></a> 
 * <h3>The EulerOperation class</h3>
 * 

 * 
 * This class implements the evaluators for the Euler problem, in analogy to
 * the `LaplaceOperator` class of step-37 or step-59. Since the present
 * operator is non-linear and does not require a matrix interface (to be
 * handed over to preconditioners), we skip the various `vmult` functions
 * otherwise present in matrix-free operators and only implement an `apply`
 * function as well as the combination of `apply` with the required vector
 * updates for the low-storage Runge--Kutta time integrator mentioned above
 * (called `perform_stage`). Furthermore, we have added three additional
 * functions involving matrix-free routines, namely one to compute an
 * estimate of the time step scaling (that is combined with the Courant
 * number for the actual time step size) based on the velocity and speed of
 * sound in the elements, one for the projection of solutions (specializing
 * VectorTools::project() for the DG case), and one to compute the errors
 * against a possible analytical solution or norms against some background
 * state.
 *   

 * 
 * The rest of the class is similar to other matrix-free tutorials. As
 * discussed in the introduction, we provide a few functions to allow a user
 * to pass in various forms of boundary conditions on different parts of the
 * domain boundary marked by types::boundary_id variables, as well as
 * possible body forces.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   class EulerOperator
 *   {
 *   public:
 *     static constexpr unsigned int n_quadrature_points_1d = n_points_1d;
 * 
 *     EulerOperator(TimerOutput &timer_output);
 * 
 *     void reinit(const Mapping<dim> &   mapping,
 *                 const DoFHandler<dim> &dof_handler);
 * 
 *     void set_inflow_boundary(const types::boundary_id       boundary_id,
 *                              std::unique_ptr<Function<dim>> inflow_function);
 * 
 *     void set_subsonic_outflow_boundary(
 *       const types::boundary_id       boundary_id,
 *       std::unique_ptr<Function<dim>> outflow_energy);
 * 
 *     void set_wall_boundary(const types::boundary_id boundary_id);
 * 
 *     void set_body_force(std::unique_ptr<Function<dim>> body_force);
 * 
 *     void apply(const double                                      current_time,
 *                const LinearAlgebra::distributed::Vector<Number> &src,
 *                LinearAlgebra::distributed::Vector<Number> &      dst) const;
 * 
 *     void
 *     perform_stage(const Number cur_time,
 *                   const Number factor_solution,
 *                   const Number factor_ai,
 *                   const LinearAlgebra::distributed::Vector<Number> &current_ri,
 *                   LinearAlgebra::distributed::Vector<Number> &      vec_ki,
 *                   LinearAlgebra::distributed::Vector<Number> &      solution,
 *                   LinearAlgebra::distributed::Vector<Number> &next_ri) const;
 * 
 *     void project(const Function<dim> &                       function,
 *                  LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     std::array<double, 3> compute_errors(
 *       const Function<dim> &                             function,
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     double compute_cell_transport_speed(
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const;
 * 
 *     void
 *     initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const;
 * 
 *   private:
 *     MatrixFree<dim, Number> data;
 * 
 *     TimerOutput &timer;
 * 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
 *       inflow_boundaries;
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
 *                                    subsonic_outflow_boundaries;
 *     std::set<types::boundary_id>   wall_boundaries;
 *     std::unique_ptr<Function<dim>> body_force;
 * 
 *     void local_apply_inverse_mass_matrix(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const;
 * 
 *     void local_apply_cell(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const;
 * 
 *     void local_apply_face(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     face_range) const;
 * 
 *     void local_apply_boundary_face(
 *       const MatrixFree<dim, Number> &                   data,
 *       LinearAlgebra::distributed::Vector<Number> &      dst,
 *       const LinearAlgebra::distributed::Vector<Number> &src,
 *       const std::pair<unsigned int, unsigned int> &     face_range) const;
 *   };
 * 
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer)
 *     : timer(timer)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * For the initialization of the Euler operator, we set up the MatrixFree
 * variable contained in the class. This can be done given a mapping to
 * describe possible curved boundaries as well as a DoFHandler object
 * describing the degrees of freedom. Since we use a discontinuous Galerkin
 * discretization in this tutorial program where no constraints are imposed
 * strongly on the solution field, we do not need to pass in an
 * AffineConstraints object and rather use a dummy for the
 * construction. With respect to quadrature, we want to select two different
 * ways of computing the underlying integrals: The first is a flexible one,
 * based on a template parameter `n_points_1d` (that will be assigned the
 * `n_q_points_1d` value specified at the top of this file). More accurate
 * integration is necessary to avoid the aliasing problem due to the
 * variable coefficients in the Euler operator. The second less accurate
 * quadrature formula is a tight one based on `fe_degree+1` and needed for
 * the inverse mass matrix. While that formula provides an exact inverse
 * only on affine element shapes and not on deformed elements, it enables
 * the fast inversion of the mass matrix by tensor product techniques,
 * necessary to ensure optimal computational efficiency overall.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::reinit(
 *     const Mapping<dim> &   mapping,
 *     const DoFHandler<dim> &dof_handler)
 *   {
 *     const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler};
 *     const AffineConstraints<double>            dummy;
 *     const std::vector<const AffineConstraints<double> *> constraints = {&dummy};
 *     const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
 *                                                     QGauss<1>(fe_degree + 1)};
 * 
 *     typename MatrixFree<dim, Number>::AdditionalData additional_data;
 *     additional_data.mapping_update_flags =
 *       (update_gradients | update_JxW_values | update_quadrature_points |
 *        update_values);
 *     additional_data.mapping_update_flags_inner_faces =
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors |
 *        update_values);
 *     additional_data.mapping_update_flags_boundary_faces =
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors |
 *        update_values);
 *     additional_data.tasks_parallel_scheme =
 *       MatrixFree<dim, Number>::AdditionalData::none;
 * 
 *     data.reinit(
 *       mapping, dof_handlers, constraints, quadratures, additional_data);
 *   }
 * 
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::initialize_vector(
 *     LinearAlgebra::distributed::Vector<Number> &vector) const
 *   {
 *     data.initialize_dof_vector(vector);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The subsequent four member functions are the ones that must be called from
 * outside to specify the various types of boundaries. For an inflow boundary,
 * we must specify all components in terms of density $\rho$, momentum $\rho
 * \mathbf{u}$ and energy $E$. Given this information, we then store the
 * function alongside the respective boundary id in a map member variable of
 * this class. Likewise, we proceed for the subsonic outflow boundaries (where
 * we request a function as well, which we use to retrieve the energy) and for
 * wall (no-penetration) boundaries where we impose zero normal velocity (no
 * function necessary, so we only request the boundary id). For the present
 * DG code where boundary conditions are solely applied as part of the weak
 * form (during time integration), the call to set the boundary conditions
 * can appear both before or after the `reinit()` call to this class. This
 * is different from continuous finite element codes where the boundary
 * conditions determine the content of the AffineConstraints object that is
 * sent into MatrixFree for initialization, thus requiring to be set before
 * the initialization of the matrix-free data structures.
 *   

 * 
 * The checks added in each of the four function are used to
 * ensure that boundary conditions are mutually exclusive on the various
 * parts of the boundary, i.e., that a user does not accidentally designate a
 * boundary as both an inflow and say a subsonic outflow boundary.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary(
 *     const types::boundary_id       boundary_id,
 *     std::unique_ptr<Function<dim>> inflow_function)
 *   {
 *     AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
 *                     subsonic_outflow_boundaries.end() &&
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as inflow"));
 *     AssertThrow(inflow_function->n_components == dim + 2,
 *                 ExcMessage("Expected function with dim+2 components"));
 * 
 *     inflow_boundaries[boundary_id] = std::move(inflow_function);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary(
 *     const types::boundary_id       boundary_id,
 *     std::unique_ptr<Function<dim>> outflow_function)
 *   {
 *     AssertThrow(inflow_boundaries.find(boundary_id) ==
 *                     inflow_boundaries.end() &&
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as subsonic outflow"));
 *     AssertThrow(outflow_function->n_components == dim + 2,
 *                 ExcMessage("Expected function with dim+2 components"));
 * 
 *     subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary(
 *     const types::boundary_id boundary_id)
 *   {
 *     AssertThrow(inflow_boundaries.find(boundary_id) ==
 *                     inflow_boundaries.end() &&
 *                   subsonic_outflow_boundaries.find(boundary_id) ==
 *                     subsonic_outflow_boundaries.end(),
 *                 ExcMessage("You already set the boundary with id " +
 *                            std::to_string(static_cast<int>(boundary_id)) +
 *                            " to another type of boundary before now setting " +
 *                            "it as wall boundary"));
 * 
 *     wall_boundaries.insert(boundary_id);
 *   }
 * 
 * 
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::set_body_force(
 *     std::unique_ptr<Function<dim>> body_force)
 *   {
 *     AssertDimension(body_force->n_components, dim);
 * 
 *     this->body_force = std::move(body_force);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Localevaluators"></a> 
 * <h4>Local evaluators</h4>
 * 

 * 
 * Now we proceed to the local evaluators for the Euler problem. The
 * evaluators are relatively simple and follow what has been presented in
 * step-37, step-48, or step-59. The first notable difference is the fact
 * that we use an FEEvaluation with a non-standard number of quadrature
 * points. Whereas we previously always set the number of quadrature points
 * to equal the polynomial degree plus one (ensuring exact integration on
 * affine element shapes), we now set the number quadrature points as a
 * separate variable (e.g. the polynomial degree plus two or three halves of
 * the polynomial degree) to more accurately handle nonlinear terms. Since
 * the evaluator is fed with the appropriate loop lengths via the template
 * argument and keeps the number of quadrature points in the whole cell in
 * the variable FEEvaluation::n_q_points, we now automatically operate on
 * the more accurate formula without further changes.
 *   

 * 
 * The second difference is due to the fact that we are now evaluating a
 * multi-component system, as opposed to the scalar systems considered
 * previously. The matrix-free framework provides several ways to handle the
 * multi-component case. The variant shown here utilizes an FEEvaluation
 * object with multiple components embedded into it, specified by the fourth
 * template argument `dim + 2` for the components in the Euler system. As a
 * consequence, the return type of FEEvaluation::get_value() is not a scalar
 * any more (that would return a VectorizedArray type, collecting data from
 * several elements), but a Tensor of `dim+2` components. The functionality
 * is otherwise similar to the scalar case; it is handled by a template
 * specialization of a base class, called FEEvaluationAccess. An alternative
 * variant would have been to use several FEEvaluation objects, a scalar one
 * for the density, a vector-valued one with `dim` components for the
 * momentum, and another scalar evaluator for the energy. To ensure that
 * those components point to the correct part of the solution, the
 * constructor of FEEvaluation takes three optional integer arguments after
 * the required MatrixFree field, namely the number of the DoFHandler for
 * multi-DoFHandler systems (taking the first by default), the number of the
 * quadrature point in case there are multiple Quadrature objects (see more
 * below), and as a third argument the component within a vector system. As
 * we have a single vector for all components, we would go with the third
 * argument, and set it to `0` for the density, `1` for the vector-valued
 * momentum, and `dim+1` for the energy slot. FEEvaluation then picks the
 * appropriate subrange of the solution vector during
 * FEEvaluationBase::read_dof_values() and
 * FEEvaluation::distributed_local_to_global() or the more compact
 * FEEvaluation::gather_evaluate() and FEEvaluation::integrate_scatter()
 * calls.
 *   

 * 
 * When it comes to the evaluation of the body force vector, we distinguish
 * between two cases for efficiency reasons: In case we have a constant
 * function (derived from Functions::ConstantFunction), we can precompute
 * the value outside the loop over quadrature points and simply use the
 * value everywhere. For a more general function, we instead need to call
 * the `evaluate_function()` method we provided above; this path is more
 * expensive because we need to access the memory associated with the
 * quadrature point data.
 *   

 * 
 * The rest follows the other tutorial programs. Since we have implemented
 * all physics for the Euler equations in the separate `euler_flux()`
 * function, all we have to do here is to call this function
 * given the current solution evaluated at quadrature points, returned by
 * `phi.get_value(q)`, and tell the FEEvaluation object to queue the flux
 * for testing it by the gradients of the shape functions (which is a Tensor
 * of outer `dim+2` components, each holding a tensor of `dim` components
 * for the $x,y,z$ component of the Euler flux). One final thing worth
 * mentioning is the order in which we queue the data for testing by the
 * value of the test function, `phi.submit_value()`, in case we are given an
 * external function: We must do this after calling `phi.get_value(q)`,
 * because `get_value()` (reading the solution) and `submit_value()`
 * (queuing the value for multiplication by the test function and summation
 * over quadrature points) access the same underlying data field. Here it
 * would be easy to achieve also without temporary variable `w_q` since
 * there is no mixing between values and gradients. For more complicated
 * setups, one has to first copy out e.g. both the value and gradient at a
 * quadrature point and then queue results again by
 * FEEvaluationBase::submit_value() and FEEvaluationBase::submit_gradient().
 *   

 * 
 * As a final note, we mention that we do not use the first MatrixFree
 * argument of this function, which is a call-back from MatrixFree::loop().
 * The interfaces imposes the present list of arguments, but since we are in
 * a member function where the MatrixFree object is already available as the
 * `data` variable, we stick with that to avoid confusion.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_cell(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data);
 * 
 *     Tensor<1, dim, VectorizedArray<Number>> constant_body_force;
 *     const Functions::ConstantFunction<dim> *constant_function =
 *       dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());
 * 
 *     if (constant_function)
 *       constant_body_force = evaluate_function<dim, Number, dim>(
 *         *constant_function, Point<dim, VectorizedArray<Number>>());
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto w_q = phi.get_value(q);
 *             phi.submit_gradient(euler_flux<dim>(w_q), q);
 *             if (body_force.get() != nullptr)
 *               {
 *                 const Tensor<1, dim, VectorizedArray<Number>> force =
 *                   constant_function ? constant_body_force :
 *                                       evaluate_function<dim, Number, dim>(
 *                                         *body_force, phi.quadrature_point(q));
 * 
 *                 Tensor<1, dim + 2, VectorizedArray<Number>> forcing;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   forcing[d + 1] = w_q[0] * force[d];
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   forcing[dim + 1] += force[d] * w_q[d + 1];
 * 
 *                 phi.submit_value(forcing, q);
 *               }
 *           }
 * 
 *         phi.integrate_scatter(((body_force.get() != nullptr) ?
 *                                  EvaluationFlags::values :
 *                                  EvaluationFlags::nothing) |
 *                                 EvaluationFlags::gradients,
 *                               dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function concerns the computation of integrals on interior
 * faces, where we need evaluators from both cells adjacent to the face. We
 * associate the variable `phi_m` with the solution component $\mathbf{w}^-$
 * and the variable `phi_p` with the solution component $\mathbf{w}^+$. We
 * distinguish the two sides in the constructor of FEFaceEvaluation by the
 * second argument, with `true` for the interior side and `false` for the
 * exterior side, with interior and exterior denoting the orientation with
 * respect to the normal vector.
 *   

 * 
 * Note that the calls FEFaceEvaluation::gather_evaluate() and
 * FEFaceEvaluation::integrate_scatter() combine the access to the vectors
 * and the sum factorization parts. This combined operation not only saves a
 * line of code, but also contains an important optimization: Given that we
 * use a nodal basis in terms of the Lagrange polynomials in the points of
 * the Gauss-Lobatto quadrature formula, only $(p+1)^{d-1}$ out of the
 * $(p+1)^d$ basis functions evaluate to non-zero on each face. Thus, the
 * evaluator only accesses the necessary data in the vector and skips the
 * parts which are multiplied by zero. If we had first read the vector, we
 * would have needed to load all data from the vector, as the call in
 * isolation would not know what data is required in subsequent
 * operations. If the subsequent FEFaceEvaluation::evaluate() call requests
 * values and derivatives, indeed all $(p+1)^d$ vector entries for each
 * component are needed, as the normal derivative is nonzero for all basis
 * functions.
 *   

 * 
 * The arguments to the evaluators as well as the procedure is similar to
 * the cell evaluation. We again use the more accurate (over-)integration
 * scheme due to the nonlinear terms, specified as the third template
 * argument in the list. At the quadrature points, we then go to our
 * free-standing function for the numerical flux. It receives the solution
 * evaluated at quadrature points from both sides (i.e., $\mathbf{w}^-$ and
 * $\mathbf{w}^+$), as well as the normal vector onto the minus side. As
 * explained above, the numerical flux is already multiplied by the normal
 * vector from the minus side. We need to switch the sign because the
 * boundary term comes with a minus sign in the weak form derived in the
 * introduction. The flux is then queued for testing both on the minus sign
 * and on the plus sign, with switched sign as the normal vector from the
 * plus side is exactly opposed to the one from the minus side.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_face(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data,
 *                                                                       true);
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data,
 *                                                                       false);
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 *         phi_p.reinit(face);
 *         phi_p.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         phi_m.reinit(face);
 *         phi_m.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
 *           {
 *             const auto numerical_flux =
 *               euler_numerical_flux<dim>(phi_m.get_value(q),
 *                                         phi_p.get_value(q),
 *                                         phi_m.get_normal_vector(q));
 *             phi_m.submit_value(-numerical_flux, q);
 *             phi_p.submit_value(numerical_flux, q);
 *           }
 * 
 *         phi_p.integrate_scatter(EvaluationFlags::values, dst);
 *         phi_m.integrate_scatter(EvaluationFlags::values, dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For faces located at the boundary, we need to impose the appropriate
 * boundary conditions. In this tutorial program, we implement four cases as
 * mentioned above. (A fifth case, for supersonic outflow conditions is
 * discussed in the "Results" section below.) The discontinuous Galerkin
 * method imposes boundary conditions not as constraints, but only
 * weakly. Thus, the various conditions are imposed by finding an appropriate
 * <i>exterior</i> quantity $\mathbf{w}^+$ that is then handed to the
 * numerical flux function also used for the interior faces. In essence,
 * we "pretend" a state on the outside of the domain in such a way that
 * if that were reality, the solution of the PDE would satisfy the boundary
 * conditions we want.
 *   

 * 
 * For wall boundaries, we need to impose a no-normal-flux condition on the
 * momentum variable, whereas we use a Neumann condition for the density and
 * energy with $\rho^+ = \rho^-$ and $E^+ = E^-$. To achieve the no-normal
 * flux condition, we set the exterior values to the interior values and
 * subtract two times the velocity in wall-normal direction, i.e., in the
 * direction of the normal vector.
 *   

 * 
 * For inflow boundaries, we simply set the given Dirichlet data
 * $\mathbf{w}_\mathrm{D}$ as a boundary value. An alternative would have been
 * to use $\mathbf{w}^+ = -\mathbf{w}^- + 2 \mathbf{w}_\mathrm{D}$, the
 * so-called mirror principle.
 *   

 * 
 * The imposition of outflow is essentially a Neumann condition, i.e.,
 * setting $\mathbf{w}^+ = \mathbf{w}^-$. For the case of subsonic outflow,
 * we still need to impose a value for the energy, which we derive from the
 * respective function. A special step is needed for the case of
 * <i>backflow</i>, i.e., the case where there is a momentum flux into the
 * domain on the Neumann portion. According to the literature (a fact that can
 * be derived by appropriate energy arguments), we must switch to another
 * variant of the flux on inflow parts, see Gravemeier, Comerford,
 * Yoshihara, Ismail, Wall, "A novel formulation for Neumann inflow
 * conditions in biomechanics", Int. J. Numer. Meth. Biomed. Eng., vol. 28
 * (2012). Here, the momentum term needs to be added once again, which
 * corresponds to removing the flux contribution on the momentum
 * variables. We do this in a post-processing step, and only for the case
 * when we both are at an outflow boundary and the dot product between the
 * normal vector and the momentum (or, equivalently, velocity) is
 * negative. As we work on data of several quadrature points at once for
 * SIMD vectorizations, we here need to explicitly loop over the array
 * entries of the SIMD array.
 *   

 * 
 * In the implementation below, we check for the various types
 * of boundaries at the level of quadrature points. Of course, we could also
 * have moved the decision out of the quadrature point loop and treat entire
 * faces as of the same kind, which avoids some map/set lookups in the inner
 * loop over quadrature points. However, the loss of efficiency is hardly
 * noticeable, so we opt for the simpler code here. Also note that the final
 * `else` clause will catch the case when some part of the boundary was not
 * assigned any boundary condition via `EulerOperator::set_..._boundary(...)`.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_boundary_face(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     face_range) const
 *   {
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true);
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face)
 *       {
 *         phi.reinit(face);
 *         phi.gather_evaluate(src, EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto w_m    = phi.get_value(q);
 *             const auto normal = phi.get_normal_vector(q);
 * 
 *             auto rho_u_dot_n = w_m[1] * normal[0];
 *             for (unsigned int d = 1; d < dim; ++d)
 *               rho_u_dot_n += w_m[1 + d] * normal[d];
 * 
 *             bool at_outflow = false;
 * 
 *             Tensor<1, dim + 2, VectorizedArray<Number>> w_p;
 *             const auto boundary_id = data.get_boundary_id(face);
 *             if (wall_boundaries.find(boundary_id) != wall_boundaries.end())
 *               {
 *                 w_p[0] = w_m[0];
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
 *                 w_p[dim + 1] = w_m[dim + 1];
 *               }
 *             else if (inflow_boundaries.find(boundary_id) !=
 *                      inflow_boundaries.end())
 *               w_p =
 *                 evaluate_function(*inflow_boundaries.find(boundary_id)->second,
 *                                   phi.quadrature_point(q));
 *             else if (subsonic_outflow_boundaries.find(boundary_id) !=
 *                      subsonic_outflow_boundaries.end())
 *               {
 *                 w_p          = w_m;
 *                 w_p[dim + 1] = evaluate_function(
 *                   *subsonic_outflow_boundaries.find(boundary_id)->second,
 *                   phi.quadrature_point(q),
 *                   dim + 1);
 *                 at_outflow = true;
 *               }
 *             else
 *               AssertThrow(false,
 *                           ExcMessage("Unknown boundary id, did "
 *                                      "you set a boundary condition for "
 *                                      "this part of the domain boundary?"));
 * 
 *             auto flux = euler_numerical_flux<dim>(w_m, w_p, normal);
 * 
 *             if (at_outflow)
 *               for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
 *                 {
 *                   if (rho_u_dot_n[v] < -1e-12)
 *                     for (unsigned int d = 0; d < dim; ++d)
 *                       flux[d + 1][v] = 0.;
 *                 }
 * 
 *             phi.submit_value(-flux, q);
 *           }
 * 
 *         phi.integrate_scatter(EvaluationFlags::values, dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function implements the inverse mass matrix operation. The
 * algorithms and rationale have been discussed extensively in the
 * introduction, so we here limit ourselves to the technicalities of the
 * MatrixFreeOperators::CellwiseInverseMassMatrix class. It does similar
 * operations as the forward evaluation of the mass matrix, except with a
 * different interpolation matrix, representing the inverse $S^{-1}$
 * factors. These represent a change of basis from the specified basis (in
 * this case, the Lagrange basis in the points of the Gauss--Lobatto
 * quadrature formula) to the Lagrange basis in the points of the Gauss
 * quadrature formula. In the latter basis, we can apply the inverse of the
 * point-wise `JxW` factor, i.e., the quadrature weight times the
 * determinant of the Jacobian of the mapping from reference to real
 * coordinates. Once this is done, the basis is changed back to the nodal
 * Gauss-Lobatto basis again. All of these operations are done by the
 * `apply()` function below. What we need to provide is the local fields to
 * operate on (which we extract from the global vector by an FEEvaluation
 * object) and write the results back to the destination vector of the mass
 * matrix operation.
 *   

 * 
 * One thing to note is that we added two integer arguments (that are
 * optional) to the constructor of FEEvaluation, the first being 0
 * (selecting among the DoFHandler in multi-DoFHandler systems; here, we
 * only have one) and the second being 1 to make the quadrature formula
 * selection. As we use the quadrature formula 0 for the over-integration of
 * nonlinear terms, we use the formula 1 with the default $p+1$ (or
 * `fe_degree+1` in terms of the variable name) points for the mass
 * matrix. This leads to square contributions to the mass matrix and ensures
 * exact integration, as explained in the introduction.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_inverse_mass_matrix(
 *     const MatrixFree<dim, Number> &,
 *     LinearAlgebra::distributed::Vector<Number> &      dst,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
 *       inverse(phi);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.read_dof_values(src);
 * 
 *         inverse.apply(phi.begin_dof_values(), phi.begin_dof_values());
 * 
 *         phi.set_dof_values(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theapplyandrelatedfunctions"></a> 
 * <h4>The apply() and related functions</h4>
 * 

 * 
 * We now come to the function which implements the evaluation of the Euler
 * operator as a whole, i.e., $\mathcal M^{-1} \mathcal L(t, \mathbf{w})$,
 * calling into the local evaluators presented above. The steps should be
 * clear from the previous code. One thing to note is that we need to adjust
 * the time in the functions we have associated with the various parts of
 * the boundary, in order to be consistent with the equation in case the
 * boundary data is time-dependent. Then, we call MatrixFree::loop() to
 * perform the cell and face integrals, including the necessary ghost data
 * exchange in the `src` vector. The seventh argument to the function,
 * `true`, specifies that we want to zero the `dst` vector as part of the
 * loop, before we start accumulating integrals into it. This variant is
 * preferred over explicitly calling `dst = 0.;` before the loop as the
 * zeroing operation is done on a subrange of the vector in parts that are
 * written by the integrals nearby. This enhances data locality and allows
 * for caching, saving one roundtrip of vector data to main memory and
 * enhancing performance. The last two arguments to the loop determine which
 * data is exchanged: Since we only access the values of the shape functions
 * one faces, typical of first-order hyperbolic problems, and since we have
 * a nodal basis with nodes at the reference element surface, we only need
 * to exchange those parts. This again saves precious memory bandwidth.
 *   

 * 
 * Once the spatial operator $\mathcal L$ is applied, we need to make a
 * second round and apply the inverse mass matrix. Here, we call
 * MatrixFree::cell_loop() since only cell integrals appear. The cell loop
 * is cheaper than the full loop as access only goes to the degrees of
 * freedom associated with the locally owned cells, which is simply the
 * locally owned degrees of freedom for DG discretizations. Thus, no ghost
 * exchange is needed here.
 *   

 * 
 * Around all these functions, we put timer scopes to record the
 * computational time for statistics about the contributions of the various
 * parts.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::apply(
 *     const double                                      current_time,
 *     const LinearAlgebra::distributed::Vector<Number> &src,
 *     LinearAlgebra::distributed::Vector<Number> &      dst) const
 *   {
 *     {
 *       TimerOutput::Scope t(timer, "apply - integrals");
 * 
 *       for (auto &i : inflow_boundaries)
 *         i.second->set_time(current_time);
 *       for (auto &i : subsonic_outflow_boundaries)
 *         i.second->set_time(current_time);
 * 
 *       data.loop(&EulerOperator::local_apply_cell,
 *                 &EulerOperator::local_apply_face,
 *                 &EulerOperator::local_apply_boundary_face,
 *                 this,
 *                 dst,
 *                 src,
 *                 true,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values);
 *     }
 * 
 *     {
 *       TimerOutput::Scope t(timer, "apply - inverse mass");
 * 
 *       data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix,
 *                      this,
 *                      dst,
 *                      dst);
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Let us move to the function that does an entire stage of a Runge--Kutta
 * update. It calls EulerOperator::apply() followed by some updates
 * to the vectors, namely `next_ri = solution + factor_ai * k_i` and
 * `solution += factor_solution * k_i`. Rather than performing these
 * steps through the vector interfaces, we here present an alternative
 * strategy that is faster on cache-based architectures. As the memory
 * consumed by the vectors is often much larger than what fits into caches,
 * the data has to effectively come from the slow RAM memory. The situation
 * can be improved by loop fusion, i.e., performing both the updates to
 * `next_ki` and `solution` within a single sweep. In that case, we would
 * read the two vectors `rhs` and `solution` and write into `next_ki` and
 * `solution`, compared to at least 4 reads and two writes in the baseline
 * case. Here, we go one step further and perform the loop immediately when
 * the mass matrix inversion has finished on a part of the
 * vector. MatrixFree::cell_loop() provides a mechanism to attach an
 * `std::function` both before the loop over cells first touches a vector
 * entry (which we do not use here, but is e.g. used for zeroing the vector)
 * and a second `std::function` to be called after the loop last touches
 * an entry. The callback is in form of a range over the given vector (in
 * terms of the local index numbering in the MPI universe) that can be
 * addressed by `local_element()` functions.
 *   

 * 
 * For this second callback, we create a lambda that works on a range and
 * write the respective update on this range. Ideally, we would add the
 * `DEAL_II_OPENMP_SIMD_PRAGMA` before the local loop to suggest to the
 * compiler to SIMD parallelize this loop (which means in practice that we
 * ensure that there is no overlap, also called aliasing, between the index
 * ranges of the pointers we use inside the loops). It turns out that at the
 * time of this writing, GCC 7.2 fails to compile an OpenMP pragma inside a
 * lambda function, so we comment this pragma out below. If your compiler is
 * newer, you should be able to uncomment these lines again.
 *   

 * 
 * Note that we select a different code path for the last
 * Runge--Kutta stage when we do not need to update the `next_ri`
 * vector. This strategy gives a considerable speedup. Whereas the inverse
 * mass matrix and vector updates take more than 60% of the computational
 * time with default vector updates on a 40-core machine, the percentage is
 * around 35% with the more optimized variant. In other words, this is a
 * speedup of around a third.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::perform_stage(
 *     const Number                                      current_time,
 *     const Number                                      factor_solution,
 *     const Number                                      factor_ai,
 *     const LinearAlgebra::distributed::Vector<Number> &current_ri,
 *     LinearAlgebra::distributed::Vector<Number> &      vec_ki,
 *     LinearAlgebra::distributed::Vector<Number> &      solution,
 *     LinearAlgebra::distributed::Vector<Number> &      next_ri) const
 *   {
 *     {
 *       TimerOutput::Scope t(timer, "rk_stage - integrals L_h");
 * 
 *       for (auto &i : inflow_boundaries)
 *         i.second->set_time(current_time);
 *       for (auto &i : subsonic_outflow_boundaries)
 *         i.second->set_time(current_time);
 * 
 *       data.loop(&EulerOperator::local_apply_cell,
 *                 &EulerOperator::local_apply_face,
 *                 &EulerOperator::local_apply_boundary_face,
 *                 this,
 *                 vec_ki,
 *                 current_ri,
 *                 true,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values,
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values);
 *     }
 * 
 * 
 *     {
 *       TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd");
 *       data.cell_loop(
 *         &EulerOperator::local_apply_inverse_mass_matrix,
 *         this,
 *         next_ri,
 *         vec_ki,
 *         std::function<void(const unsigned int, const unsigned int)>(),
 *         [&](const unsigned int start_range, const unsigned int end_range) {
 *           const Number ai = factor_ai;
 *           const Number bi = factor_solution;
 *           if (ai == Number())
 *             {
 *               /* DEAL_II_OPENMP_SIMD_PRAGMA */
 *               for (unsigned int i = start_range; i < end_range; ++i)
 *                 {
 *                   const Number k_i          = next_ri.local_element(i);
 *                   const Number sol_i        = solution.local_element(i);
 *                   solution.local_element(i) = sol_i + bi * k_i;
 *                 }
 *             }
 *           else
 *             {
 *               /* DEAL_II_OPENMP_SIMD_PRAGMA */
 *               for (unsigned int i = start_range; i < end_range; ++i)
 *                 {
 *                   const Number k_i          = next_ri.local_element(i);
 *                   const Number sol_i        = solution.local_element(i);
 *                   solution.local_element(i) = sol_i + bi * k_i;
 *                   next_ri.local_element(i)  = sol_i + ai * k_i;
 *                 }
 *             }
 *         });
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Having discussed the implementation of the functions that deal with
 * advancing the solution by one time step, let us now move to functions
 * that implement other, ancillary operations. Specifically, these are
 * functions that compute projections, evaluate errors, and compute the speed
 * of information transport on a cell.
 *   

 * 
 * The first of these functions is essentially equivalent to
 * VectorTools::project(), just much faster because it is specialized for DG
 * elements where there is no need to set up and solve a linear system, as
 * each element has independent basis functions. The reason why we show the
 * code here, besides a small speedup of this non-critical operation, is that
 * it shows additional functionality provided by
 * MatrixFreeOperators::CellwiseInverseMassMatrix.
 *   

 * 
 * The projection operation works as follows: If we denote the matrix of
 * shape functions evaluated at quadrature points by $S$, the projection on
 * cell $K$ is an operation of the form $\underbrace{S J^K S^\mathrm
 * T}_{\mathcal M^K} \mathbf{w}^K = S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$, where $J^K$ is the diagonal
 * matrix containing the determinant of the Jacobian times the quadrature
 * weight (JxW), $\mathcal M^K$ is the cell-wise mass matrix, and
 * $\tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ is the evaluation of the
 * field to be projected onto quadrature points. (In reality the matrix $S$
 * has additional structure through the tensor product, as explained in the
 * introduction.) This system can now equivalently be written as
 * $\mathbf{w}^K = \left(S J^K S^\mathrm T\right)^{-1} S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q} = S^{-\mathrm T}
 * \left(J^K\right)^{-1} S^{-1} S J^K
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. Now, the term $S^{-1} S$ and
 * then $\left(J^K\right)^{-1} J^K$ cancel, resulting in the final
 * expression $\mathbf{w}^K = S^{-\mathrm T}
 * \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. This operation is
 * implemented by
 * MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis().
 * The name is derived from the fact that this projection is simply
 * the multiplication by $S^{-\mathrm T}$, a basis change from the
 * nodal basis in the points of the Gaussian quadrature to the given finite
 * element basis. Note that we call FEEvaluation::set_dof_values() to write
 * the result into the vector, overwriting previous content, rather than
 * accumulating the results as typical in integration tasks -- we can do
 * this because every vector entry has contributions from only a single
 * cell for discontinuous Galerkin discretizations.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   void EulerOperator<dim, degree, n_points_1d>::project(
 *     const Function<dim> &                       function,
 *     LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
 *       inverse(phi);
 *     solution.zero_out_ghost_values();
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_dof_value(evaluate_function(function,
 *                                                  phi.quadrature_point(q)),
 *                                q);
 *         inverse.transform_from_q_points_to_basis(dim + 2,
 *                                                  phi.begin_dof_values(),
 *                                                  phi.begin_dof_values());
 *         phi.set_dof_values(solution);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function again repeats functionality also provided by the
 * deal.II library, namely VectorTools::integrate_difference(). We here show
 * the explicit code to highlight how the vectorization across several cells
 * works and how to accumulate results via that interface: Recall that each
 * <i>lane</i> of the vectorized array holds data from a different cell. By
 * the loop over all cell batches that are owned by the current MPI process,
 * we could then fill a VectorizedArray of results; to obtain a global sum,
 * we would need to further go on and sum across the entries in the SIMD
 * array. However, such a procedure is not stable as the SIMD array could in
 * fact not hold valid data for all its lanes. This happens when the number
 * of locally owned cells is not a multiple of the SIMD width. To avoid
 * invalid data, we must explicitly skip those invalid lanes when accessing
 * the data. While one could imagine that we could make it work by simply
 * setting the empty lanes to zero (and thus, not contribute to a sum), the
 * situation is more complicated than that: What if we were to compute a
 * velocity out of the momentum? Then, we would need to divide by the
 * density, which is zero -- the result would consequently be NaN and
 * contaminate the result. This trap is avoided by accumulating the results
 * from the valid SIMD range as we loop through the cell batches, using the
 * function MatrixFree::n_active_entries_per_cell_batch() to give us the
 * number of lanes with valid data. It equals VectorizedArray::size() on
 * most cells, but can be less on the last cell batch if the number of cells
 * has a remainder compared to the SIMD width.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors(
 *     const Function<dim> &                             function,
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     TimerOutput::Scope t(timer, "compute errors");
 *     double             errors_squared[3] = {};
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0);
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(solution, EvaluationFlags::values);
 *         VectorizedArray<Number> local_errors_squared[3] = {};
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto error =
 *               evaluate_function(function, phi.quadrature_point(q)) -
 *               phi.get_value(q);
 *             const auto JxW = phi.JxW(q);
 * 
 *             local_errors_squared[0] += error[0] * error[0] * JxW;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW;
 *             local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW;
 *           }
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
 *              ++v)
 *           for (unsigned int d = 0; d < 3; ++d)
 *             errors_squared[d] += local_errors_squared[d][v];
 *       }
 * 
 *     Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);
 * 
 *     std::array<double, 3> errors;
 *     for (unsigned int d = 0; d < 3; ++d)
 *       errors[d] = std::sqrt(errors_squared[d]);
 * 
 *     return errors;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This final function of the EulerOperator class is used to estimate the
 * transport speed, scaled by the mesh size, that is relevant for setting
 * the time step size in the explicit time integrator. In the Euler
 * equations, there are two speeds of transport, namely the convective
 * velocity $\mathbf{u}$ and the propagation of sound waves with sound
 * speed $c = \sqrt{\gamma p/\rho}$ relative to the medium moving at
 * velocity $\mathbf u$.
 *   

 * 
 * In the formula for the time step size, we are interested not by
 * these absolute speeds, but by the amount of time it takes for
 * information to cross a single cell. For information transported along with
 * the medium, $\mathbf u$ is scaled by the mesh size,
 * so an estimate of the maximal velocity can be obtained by computing
 * $\|J^{-\mathrm T} \mathbf{u}\|_\infty$, where $J$ is the Jacobian of the
 * transformation from real to the reference domain. Note that
 * FEEvaluationBase::inverse_jacobian() returns the inverse and transpose
 * Jacobian, representing the metric term from real to reference
 * coordinates, so we do not need to transpose it again. We store this limit
 * in the variable `convective_limit` in the code below.
 *   

 * 
 * The sound propagation is isotropic, so we need to take mesh sizes in any
 * direction into account. The appropriate mesh size scaling is then given
 * by the minimal singular value of $J$ or, equivalently, the maximal
 * singular value of $J^{-1}$. Note that one could approximate this quantity
 * by the minimal distance between vertices of a cell when ignoring curved
 * cells. To get the maximal singular value of the Jacobian, the general
 * strategy would be some LAPACK function. Since all we need here is an
 * estimate, we can avoid the hassle of decomposing a tensor of
 * VectorizedArray numbers into several matrices and go into an (expensive)
 * eigenvalue function without vectorization, and instead use a few
 * iterations (five in the code below) of the power method applied to
 * $J^{-1}J^{-\mathrm T}$. The speed of convergence of this method depends
 * on the ratio of the largest to the next largest eigenvalue and the
 * initial guess, which is the vector of all ones. This might suggest that
 * we get slow convergence on cells close to a cube shape where all
 * lengths are almost the same. However, this slow convergence means that
 * the result will sit between the two largest singular values, which both
 * are close to the maximal value anyway. In all other cases, convergence
 * will be quick. Thus, we can merely hardcode 5 iterations here and be
 * confident that the result is good.
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d>
 *   double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed(
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const
 *   {
 *     TimerOutput::Scope t(timer, "compute transport speed");
 *     Number             max_transport = 0;
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         phi.reinit(cell);
 *         phi.gather_evaluate(solution, EvaluationFlags::values);
 *         VectorizedArray<Number> local_max = 0.;
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           {
 *             const auto solution = phi.get_value(q);
 *             const auto velocity = euler_velocity<dim>(solution);
 *             const auto pressure = euler_pressure<dim>(solution);
 * 
 *             const auto inverse_jacobian = phi.inverse_jacobian(q);
 *             const auto convective_speed = inverse_jacobian * velocity;
 *             VectorizedArray<Number> convective_limit = 0.;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               convective_limit =
 *                 std::max(convective_limit, std::abs(convective_speed[d]));
 * 
 *             const auto speed_of_sound =
 *               std::sqrt(gamma * pressure * (1. / solution[0]));
 * 
 *             Tensor<1, dim, VectorizedArray<Number>> eigenvector;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               eigenvector[d] = 1.;
 *             for (unsigned int i = 0; i < 5; ++i)
 *               {
 *                 eigenvector = transpose(inverse_jacobian) *
 *                               (inverse_jacobian * eigenvector);
 *                 VectorizedArray<Number> eigenvector_norm = 0.;
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   eigenvector_norm =
 *                     std::max(eigenvector_norm, std::abs(eigenvector[d]));
 *                 eigenvector /= eigenvector_norm;
 *               }
 *             const auto jac_times_ev   = inverse_jacobian * eigenvector;
 *             const auto max_eigenvalue = std::sqrt(
 *               (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector));
 *             local_max =
 *               std::max(local_max,
 *                        max_eigenvalue * speed_of_sound + convective_limit);
 *           }
 * 
 * @endcode
 * 
 * Similarly to the previous function, we must make sure to accumulate
 * speed only on the valid cells of a cell batch.
 * 
 * @code
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
 *              ++v)
 *           for (unsigned int d = 0; d < 3; ++d)
 *             max_transport = std::max(max_transport, local_max[v]);
 *       }
 * 
 *     max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);
 * 
 *     return max_transport;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerProblemclass"></a> 
 * <h3>The EulerProblem class</h3>
 * 

 * 
 * This class combines the EulerOperator class with the time integrator and
 * the usual global data structures such as FiniteElement and DoFHandler, to
 * actually run the simulations of the Euler problem.
 *   

 * 
 * The member variables are a triangulation, a finite element, a mapping (to
 * create high-order curved surfaces, see e.g. step-10), and a DoFHandler to
 * describe the degrees of freedom. In addition, we keep an instance of the
 * EulerOperator described above around, which will do all heavy lifting in
 * terms of integrals, and some parameters for time integration like the
 * current time or the time step size.
 *   

 * 
 * Furthermore, we use a PostProcessor instance to write some additional
 * information to the output file, in similarity to what was done in
 * step-33. The interface of the DataPostprocessor class is intuitive,
 * requiring us to provide information about what needs to be evaluated
 * (typically only the values of the solution, except for the Schlieren plot
 * that we only enable in 2D where it makes sense), and the names of what
 * gets evaluated. Note that it would also be possible to extract most
 * information by calculator tools within visualization programs such as
 * ParaView, but it is so much more convenient to do it already when writing
 * the output.
 * 
 * @code
 *   template <int dim>
 *   class EulerProblem
 *   {
 *   public:
 *     EulerProblem();
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid_and_dofs();
 * 
 *     void output_results(const unsigned int result_number);
 * 
 *     LinearAlgebra::distributed::Vector<Number> solution;
 * 
 *     ConditionalOStream pcout;
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 * 
 *     FESystem<dim>        fe;
 *     MappingQGeneric<dim> mapping;
 *     DoFHandler<dim>      dof_handler;
 * 
 *     TimerOutput timer;
 * 
 *     EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator;
 * 
 *     double time, time_step;
 * 
 *     class Postprocessor : public DataPostprocessor<dim>
 *     {
 *     public:
 *       Postprocessor();
 * 
 *       virtual void evaluate_vector_field(
 *         const DataPostprocessorInputs::Vector<dim> &inputs,
 *         std::vector<Vector<double>> &computed_quantities) const override;
 * 
 *       virtual std::vector<std::string> get_names() const override;
 * 
 *       virtual std::vector<
 *         DataComponentInterpretation::DataComponentInterpretation>
 *       get_data_component_interpretation() const override;
 * 
 *       virtual UpdateFlags get_needed_update_flags() const override;
 * 
 *     private:
 *       const bool do_schlieren_plot;
 *     };
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   EulerProblem<dim>::Postprocessor::Postprocessor()
 *     : do_schlieren_plot(dim == 2)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * For the main evaluation of the field variables, we first check that the
 * lengths of the arrays equal the expected values (the lengths `2*dim+4` or
 * `2*dim+5` are derived from the sizes of the names we specify in the
 * get_names() function below). Then we loop over all evaluation points and
 * fill the respective information: First we fill the primal solution
 * variables of density $\rho$, momentum $\rho \mathbf{u}$ and energy $E$,
 * then we compute the derived velocity $\mathbf u$, the pressure $p$, the
 * speed of sound $c=\sqrt{\gamma p / \rho}$, as well as the Schlieren plot
 * showing $s = |\nabla \rho|^2$ in case it is enabled. (See step-69 for
 * another example where we create a Schlieren plot.)
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::Postprocessor::evaluate_vector_field(
 *     const DataPostprocessorInputs::Vector<dim> &inputs,
 *     std::vector<Vector<double>> &               computed_quantities) const
 *   {
 *     const unsigned int n_evaluation_points = inputs.solution_values.size();
 * 
 *     if (do_schlieren_plot == true)
 *       Assert(inputs.solution_gradients.size() == n_evaluation_points,
 *              ExcInternalError());
 * 
 *     Assert(computed_quantities.size() == n_evaluation_points,
 *            ExcInternalError());
 *     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
 *     Assert(computed_quantities[0].size() ==
 *              dim + 2 + (do_schlieren_plot == true ? 1 : 0),
 *            ExcInternalError());
 * 
 *     for (unsigned int q = 0; q < n_evaluation_points; ++q)
 *       {
 *         Tensor<1, dim + 2> solution;
 *         for (unsigned int d = 0; d < dim + 2; ++d)
 *           solution[d] = inputs.solution_values[q](d);
 * 
 *         const double         density  = solution[0];
 *         const Tensor<1, dim> velocity = euler_velocity<dim>(solution);
 *         const double         pressure = euler_pressure<dim>(solution);
 * 
 *         for (unsigned int d = 0; d < dim; ++d)
 *           computed_quantities[q](d) = velocity[d];
 *         computed_quantities[q](dim)     = pressure;
 *         computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density);
 * 
 *         if (do_schlieren_plot == true)
 *           computed_quantities[q](dim + 2) =
 *             inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0];
 *       }
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const
 *   {
 *     std::vector<std::string> names;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       names.emplace_back("velocity");
 *     names.emplace_back("pressure");
 *     names.emplace_back("speed_of_sound");
 * 
 *     if (do_schlieren_plot == true)
 *       names.emplace_back("schlieren_plot");
 * 
 *     return names;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For the interpretation of quantities, we have scalar density, energy,
 * pressure, speed of sound, and the Schlieren plot, and vectors for the
 * momentum and the velocity.
 * 
 * @code
 *   template <int dim>
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *   EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const
 *   {
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       interpretation;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       interpretation.push_back(
 *         DataComponentInterpretation::component_is_part_of_vector);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 * 
 *     if (do_schlieren_plot == true)
 *       interpretation.push_back(
 *         DataComponentInterpretation::component_is_scalar);
 * 
 *     return interpretation;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * With respect to the necessary update flags, we only need the values for
 * all quantities but the Schlieren plot, which is based on the density
 * gradient.
 * 
 * @code
 *   template <int dim>
 *   UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const
 *   {
 *     if (do_schlieren_plot == true)
 *       return update_values | update_gradients;
 *     else
 *       return update_values;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The constructor for this class is unsurprising: We set up a parallel
 * triangulation based on the `MPI_COMM_WORLD` communicator, a vector finite
 * element with `dim+2` components for density, momentum, and energy, a
 * high-order mapping of the same degree as the underlying finite element,
 * and initialize the time and time step to zero.
 * 
 * @code
 *   template <int dim>
 *   EulerProblem<dim>::EulerProblem()
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 * #ifdef DEAL_II_WITH_P4EST
 *     , triangulation(MPI_COMM_WORLD)
 * #endif
 *     , fe(FE_DGQ<dim>(fe_degree), dim + 2)
 *     , mapping(fe_degree)
 *     , dof_handler(triangulation)
 *     , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
 *     , euler_operator(timer)
 *     , time(0)
 *     , time_step(0)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * As a mesh, this tutorial program implements two options, depending on the
 * global variable `testcase`: For the analytical variant (`testcase==0`),
 * the domain is $(0, 10) \times (-5, 5)$, with Dirichlet boundary
 * conditions (inflow) all around the domain. For `testcase==1`, we set the
 * domain to a cylinder in a rectangular box, derived from the flow past
 * cylinder testcase for incompressible viscous flow by Sch&auml;fer and
 * Turek (1996). Here, we have a larger variety of boundaries. The inflow
 * part at the left of the channel is given the inflow type, for which we
 * choose a constant inflow profile, whereas we set a subsonic outflow at
 * the right. For the boundary around the cylinder (boundary id equal to 2)
 * as well as the channel walls (boundary id equal to 3) we use the wall
 * boundary type, which is no-normal flow. Furthermore, for the 3D cylinder
 * we also add a gravity force in vertical direction. Having the base mesh
 * in place (including the manifolds set by
 * GridGenerator::channel_with_cylinder()), we can then perform the
 * specified number of global refinements, create the unknown numbering from
 * the DoFHandler, and hand the DoFHandler and Mapping objects to the
 * initialization of the EulerOperator.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::make_grid_and_dofs()
 *   {
 *     switch (testcase)
 *       {
 *         case 0:
 *           {
 *             Point<dim> lower_left;
 *             for (unsigned int d = 1; d < dim; ++d)
 *               lower_left[d] = -5;
 * 
 *             Point<dim> upper_right;
 *             upper_right[0] = 10;
 *             for (unsigned int d = 1; d < dim; ++d)
 *               upper_right[d] = 5;
 * 
 *             GridGenerator::hyper_rectangle(triangulation,
 *                                            lower_left,
 *                                            upper_right);
 *             triangulation.refine_global(2);
 * 
 *             euler_operator.set_inflow_boundary(
 *               0, std::make_unique<ExactSolution<dim>>(0));
 * 
 *             break;
 *           }
 * 
 *         case 1:
 *           {
 *             GridGenerator::channel_with_cylinder(
 *               triangulation, 0.03, 1, 0, true);
 * 
 *             euler_operator.set_inflow_boundary(
 *               0, std::make_unique<ExactSolution<dim>>(0));
 *             euler_operator.set_subsonic_outflow_boundary(
 *               1, std::make_unique<ExactSolution<dim>>(0));
 * 
 *             euler_operator.set_wall_boundary(2);
 *             euler_operator.set_wall_boundary(3);
 * 
 *             if (dim == 3)
 *               euler_operator.set_body_force(
 *                 std::make_unique<Functions::ConstantFunction<dim>>(
 *                   std::vector<double>({0., 0., -0.2})));
 * 
 *             break;
 *           }
 * 
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     triangulation.refine_global(n_global_refinements);
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     euler_operator.reinit(mapping, dof_handler);
 *     euler_operator.initialize_vector(solution);
 * 
 * @endcode
 * 
 * In the following, we output some statistics about the problem. Because we
 * often end up with quite large numbers of cells or degrees of freedom, we
 * would like to print them with a comma to separate each set of three
 * digits. This can be done via "locales", although the way this works is
 * not particularly intuitive. step-32 explains this in slightly more
 * detail.
 * 
 * @code
 *     std::locale s = pcout.get_stream().getloc();
 *     pcout.get_stream().imbue(std::locale(""));
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << " ( = " << (dim + 2) << " [vars] x "
 *           << triangulation.n_global_active_cells() << " [cells] x "
 *           << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )"
 *           << std::endl;
 *     pcout.get_stream().imbue(s);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * For output, we first let the Euler operator compute the errors of the
 * numerical results. More precisely, we compute the error against the
 * analytical result for the analytical solution case, whereas we compute
 * the deviation against the background field with constant density and
 * energy and constant velocity in $x$ direction for the second test case.
 *   

 * 
 * The next step is to create output. This is similar to what is done in
 * step-33: We let the postprocessor defined above control most of the
 * output, except for the primal field that we write directly. For the
 * analytical solution test case, we also perform another projection of the
 * analytical solution and print the difference between that field and the
 * numerical solution. Once we have defined all quantities to be written, we
 * build the patches for output. Similarly to step-65, we create a
 * high-order VTK output by setting the appropriate flag, which enables us
 * to visualize fields of high polynomial degrees. Finally, we call the
 * `DataOutInterface::write_vtu_in_parallel()` function to write the result
 * to the given file name. This function uses special MPI parallel write
 * facilities, which are typically more optimized for parallel file systems
 * than the standard library's `std::ofstream` variants used in most other
 * tutorial programs. A particularly nice feature of the
 * `write_vtu_in_parallel()` function is the fact that it can combine output
 * from all MPI ranks into a single file, making it unnecessary to have a
 * central record of all such files (namely, the "pvtu" file).
 *   

 * 
 * For parallel programs, it is often instructive to look at the partitioning
 * of cells among processors. To this end, one can pass a vector of numbers
 * to DataOut::add_data_vector() that contains as many entries as the
 * current processor has active cells; these numbers should then be the
 * rank of the processor that owns each of these cells. Such a vector
 * could, for example, be obtained from
 * GridTools::get_subdomain_association(). On the other hand, on each MPI
 * process, DataOut will only read those entries that correspond to locally
 * owned cells, and these of course all have the same value: namely, the rank
 * of the current process. What is in the remaining entries of the vector
 * doesn't actually matter, and so we can just get away with a cheap trick: We
 * just fill *all* values of the vector we give to DataOut::add_data_vector()
 * with the rank of the current MPI process. The key is that on each process,
 * only the entries corresponding to the locally owned cells will be read,
 * ignoring the (wrong) values in other entries. The fact that every process
 * submits a vector in which the correct subset of entries is correct is all
 * that is necessary.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::output_results(const unsigned int result_number)
 *   {
 *     const std::array<double, 3> errors =
 *       euler_operator.compute_errors(ExactSolution<dim>(time), solution);
 *     const std::string quantity_name = testcase == 0 ? "error" : "norm";
 * 
 *     pcout << "Time:" << std::setw(8) << std::setprecision(3) << time
 *           << ", dt: " << std::setw(8) << std::setprecision(2) << time_step
 *           << ", " << quantity_name << " rho: " << std::setprecision(4)
 *           << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4)
 *           << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4)
 *           << std::setw(10) << errors[2] << std::endl;
 * 
 *     {
 *       TimerOutput::Scope t(timer, "output");
 * 
 *       Postprocessor postprocessor;
 *       DataOut<dim>  data_out;
 * 
 *       DataOutBase::VtkFlags flags;
 *       flags.write_higher_order_cells = true;
 *       data_out.set_flags(flags);
 * 
 *       data_out.attach_dof_handler(dof_handler);
 *       {
 *         std::vector<std::string> names;
 *         names.emplace_back("density");
 *         for (unsigned int d = 0; d < dim; ++d)
 *           names.emplace_back("momentum");
 *         names.emplace_back("energy");
 * 
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *           interpretation;
 *         interpretation.push_back(
 *           DataComponentInterpretation::component_is_scalar);
 *         for (unsigned int d = 0; d < dim; ++d)
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_part_of_vector);
 *         interpretation.push_back(
 *           DataComponentInterpretation::component_is_scalar);
 * 
 *         data_out.add_data_vector(dof_handler, solution, names, interpretation);
 *       }
 *       data_out.add_data_vector(solution, postprocessor);
 * 
 *       LinearAlgebra::distributed::Vector<Number> reference;
 *       if (testcase == 0 && dim == 2)
 *         {
 *           reference.reinit(solution);
 *           euler_operator.project(ExactSolution<dim>(time), reference);
 *           reference.sadd(-1., 1, solution);
 *           std::vector<std::string> names;
 *           names.emplace_back("error_density");
 *           for (unsigned int d = 0; d < dim; ++d)
 *             names.emplace_back("error_momentum");
 *           names.emplace_back("error_energy");
 * 
 *           std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *             interpretation;
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_scalar);
 *           for (unsigned int d = 0; d < dim; ++d)
 *             interpretation.push_back(
 *               DataComponentInterpretation::component_is_part_of_vector);
 *           interpretation.push_back(
 *             DataComponentInterpretation::component_is_scalar);
 * 
 *           data_out.add_data_vector(dof_handler,
 *                                    reference,
 *                                    names,
 *                                    interpretation);
 *         }
 * 
 *       Vector<double> mpi_owner(triangulation.n_active_cells());
 *       mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
 *       data_out.add_data_vector(mpi_owner, "owner");
 * 
 *       data_out.build_patches(mapping,
 *                              fe.degree,
 *                              DataOut<dim>::curved_inner_cells);
 * 
 *       const std::string filename =
 *         "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu";
 *       data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The EulerProblem::run() function puts all pieces together. It starts off
 * by calling the function that creates the mesh and sets up data structures,
 * and then initializing the time integrator and the two temporary vectors of
 * the low-storage integrator. We call these vectors `rk_register_1` and
 * `rk_register_2`, and use the first vector to represent the quantity
 * $\mathbf{r}_i$ and the second one for $\mathbf{k}_i$ in the formulas for
 * the Runge--Kutta scheme outlined in the introduction. Before we start the
 * time loop, we compute the time step size by the
 * `EulerOperator::compute_cell_transport_speed()` function. For reasons of
 * comparison, we compare the result obtained there with the minimal mesh
 * size and print them to screen. For velocities and speeds of sound close
 * to unity as in this tutorial program, the predicted effective mesh size
 * will be close, but they could vary if scaling were different.
 * 
 * @code
 *   template <int dim>
 *   void EulerProblem<dim>::run()
 *   {
 *     {
 *       const unsigned int n_vect_number = VectorizedArray<Number>::size();
 *       const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number;
 * 
 *       pcout << "Running with "
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
 *             << " MPI processes" << std::endl;
 *       pcout << "Vectorization over " << n_vect_number << " "
 *             << (std::is_same<Number, double>::value ? "doubles" : "floats")
 *             << " = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl;
 *     }
 * 
 *     make_grid_and_dofs();
 * 
 *     const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme);
 * 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_1;
 *     LinearAlgebra::distributed::Vector<Number> rk_register_2;
 *     rk_register_1.reinit(solution);
 *     rk_register_2.reinit(solution);
 * 
 *     euler_operator.project(ExactSolution<dim>(time), solution);
 * 
 *     double min_vertex_distance = std::numeric_limits<double>::max();
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         min_vertex_distance =
 *           std::min(min_vertex_distance, cell->minimum_vertex_distance());
 *     min_vertex_distance =
 *       Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);
 * 
 *     time_step = courant_number * integrator.n_stages() /
 *                 euler_operator.compute_cell_transport_speed(solution);
 *     pcout << "Time step size: " << time_step
 *           << ", minimal h: " << min_vertex_distance
 *           << ", initial transport scaling: "
 *           << 1. / euler_operator.compute_cell_transport_speed(solution)
 *           << std::endl
 *           << std::endl;
 * 
 *     output_results(0);
 * 
 * @endcode
 * 
 * Now we are ready to start the time loop, which we run until the time
 * has reached the desired end time. Every 5 time steps, we compute a new
 * estimate for the time step -- since the solution is nonlinear, it is
 * most effective to adapt the value during the course of the
 * simulation. In case the Courant number was chosen too aggressively, the
 * simulation will typically blow up with time step NaN, so that is easy
 * to detect here. One thing to note is that roundoff errors might
 * propagate to the leading digits due to an interaction of slightly
 * different time step selections that in turn lead to slightly different
 * solutions. To decrease this sensitivity, it is common practice to round
 * or truncate the time step size to a few digits, e.g. 3 in this case. In
 * case the current time is near the prescribed 'tick' value for output
 * (e.g. 0.02), we also write the output. After the end of the time loop,
 * we summarize the computation by printing some statistics, which is
 * mostly done by the TimerOutput::print_wall_time_statistics() function.
 * 
 * @code
 *     unsigned int timestep_number = 0;
 * 
 *     while (time < final_time - 1e-12)
 *       {
 *         ++timestep_number;
 *         if (timestep_number % 5 == 0)
 *           time_step =
 *             courant_number * integrator.n_stages() /
 *             Utilities::truncate_to_n_digits(
 *               euler_operator.compute_cell_transport_speed(solution), 3);
 * 
 *         {
 *           TimerOutput::Scope t(timer, "rk time stepping total");
 *           integrator.perform_time_step(euler_operator,
 *                                        time,
 *                                        time_step,
 *                                        solution,
 *                                        rk_register_1,
 *                                        rk_register_2);
 *         }
 * 
 *         time += time_step;
 * 
 *         if (static_cast<int>(time / output_tick) !=
 *               static_cast<int>((time - time_step) / output_tick) ||
 *             time >= final_time - 1e-12)
 *           output_results(
 *             static_cast<unsigned int>(std::round(time / output_tick)));
 *       }
 * 
 *     timer.print_wall_time_statistics(MPI_COMM_WORLD);
 *     pcout << std::endl;
 *   }
 * 
 * } // namespace Euler_DG
 * 
 * 
 * 
 * @endcode
 * 
 * The main() function is not surprising and follows what was done in all
 * previous MPI programs: As we run an MPI program, we need to call `MPI_Init()`
 * and `MPI_Finalize()`, which we do through the
 * Utilities::MPI::MPI_InitFinalize data structure. Note that we run the program
 * only with MPI, and set the thread count to 1.
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   using namespace Euler_DG;
 *   using namespace dealii;
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *   try
 *     {
 *       deallog.depth_console(0);
 * 
 *       EulerProblem<dimension> euler_problem;
 *       euler_problem.run();
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
examples/step-67/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


在一台有40个进程的机器上以默认设置运行该程序，会产生以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 147,456 ( = 4 [vars] x 1,024 [cells] x 36 [dofs/cell/var] )
Time step size: 0.00689325, minimal h: 0.3125, initial transport scaling: 0.102759


Time:       0, dt:   0.0069, error rho:   2.76e-07, rho * u:  1.259e-06, energy: 2.987e-06
Time:    1.01, dt:   0.0069, error rho:   1.37e-06, rho * u:  2.252e-06, energy: 4.153e-06
Time:    2.01, dt:   0.0069, error rho:  1.561e-06, rho * u:   2.43e-06, energy: 4.493e-06
Time:    3.01, dt:   0.0069, error rho:  1.714e-06, rho * u:  2.591e-06, energy: 4.762e-06
Time:    4.01, dt:   0.0069, error rho:  1.843e-06, rho * u:  2.625e-06, energy: 4.985e-06
Time:    5.01, dt:   0.0069, error rho:  1.496e-06, rho * u:  1.961e-06, energy: 4.142e-06
Time:       6, dt:   0.0083, error rho:  1.007e-06, rho * u:  7.119e-07, energy: 2.972e-06
Time:       7, dt:   0.0095, error rho:  9.096e-07, rho * u:  3.786e-07, energy: 2.626e-06
Time:       8, dt:   0.0096, error rho:  8.439e-07, rho * u:  3.338e-07, energy:  2.43e-06
Time:       9, dt:   0.0096, error rho:  7.822e-07, rho * u:  2.984e-07, energy: 2.248e-06
Time:      10, dt:   0.0096, error rho:  7.231e-07, rho * u:  2.666e-07, energy: 2.074e-06


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     2.249s    30 |     2.249s |     2.249s     8 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |  0.008066s    13 |   0.00952s |   0.01041s    20 |
| compute transport speed       |       258 |   0.01012s    13 |   0.05392s |   0.08574s    25 |
| output                        |        11 |    0.9597s    13 |    0.9613s |    0.9623s     6 |
| rk time stepping total        |      1283 |    0.9827s    25 |     1.015s |      1.06s    13 |
| rk_stage - integrals L_h      |      6415 |    0.8803s    26 |    0.9198s |    0.9619s    14 |
| rk_stage - inv mass + vec upd |      6415 |   0.05677s    15 |   0.06487s |   0.07597s    13 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



程序输出显示，所有的误差都很小。这是由于我们使用了一个相对较细的 $32^2$ 单元的网格，用5度的多项式来求得一个平滑的解决方案。一个有趣的模式显示在时间步长上：虽然在时间5之前是0.0069，但在后来的时间里增加到0.0096。在时间5和6.5之间，一旦在声速之上有一些运动的旋涡（因此传播速度更快）离开计算域，步长就会增加。在这之后，气流只是在同一方向上是均匀的，与之前均匀速度被漩涡覆盖的状态相比，气体的最大速度有所下降。我们的时间步长公式认识到了这种影响。

最后一块输出显示了关于程序各个部分时间的详细信息；它通过显示最快和最慢的处理器所花费的时间以及平均时间将其分解开来--这在非常大的计算中通常很有用，可以发现是否有处理器持续过热（并因此节制其时钟速度）或因其他原因持续过慢。总结显示，在1.02秒内完成了1283个时间步骤（看所有MPI进程的平均时间），而11个文件的输出又花了0.96秒。将每个时间步数和五个Runge--Kutta阶段分解开来，每次评估的计算时间为0.16毫秒。这种高性能是无矩阵评估器的典型表现，也是显式时间积分对隐式求解器非常有竞争力的原因，特别是对于大规模模拟。程序运行结束时的计算时间细分显示， $\mathcal L_h$ 中的积分评估贡献了大约0.92秒，反质量矩阵的应用贡献了0.06秒。此外，对时间步长计算的运输速度的估计又贡献了0.05秒的计算时间。

如果我们再使用三个级别的全局细化和总共940万个DoF，最终的统计数据如下（对于修改后的Lax--Friedrichs通量， $p=5$  ，和同一系统的40个核心的双插槽Intel Xeon Gold 6230）。

@code
+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     244.9s    12 |     244.9s |     244.9s    34 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |    0.4239s    12 |    0.4318s |    0.4408s     9 |
| compute transport speed       |      2053 |     3.962s    12 |     6.727s |     10.12s     7 |
| output                        |        11 |     30.35s    12 |     30.36s |     30.37s     9 |
| rk time stepping total        |     10258 |     201.7s     7 |     205.1s |     207.8s    12 |
| rk_stage - integrals L_h      |     51290 |     121.3s     6 |     126.6s |     136.3s    16 |
| rk_stage - inv mass + vec upd |     51290 |     66.19s    16 |     77.52s |     81.84s    10 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



每个时间步长，求解器现在需要0.02秒，大约是147k未知数的小问题的25倍。鉴于该问题涉及64倍的未知数，计算时间的增加并不令人惊讶。由于我们也做了8倍的时间步数，计算时间在理论上应该增加512倍。实际增加的时间是205秒/1.02秒=202。这是因为由于通信开销的原因，小问题的规模不能充分利用40个核心。如果我们研究一下每个时间步长所做操作的细节，这一点就很清楚了。带有近邻通信的微分算子 $\mathcal L_h$ 的评估时间从0.92秒到127秒，也就是说，它增加了138倍。另一方面，应用反质量矩阵和向量更新的成本，完全不需要在MPI进程之间通信，增加了1195倍。这一增长超过了理论上的512倍，因为对于较大的尺寸，操作受限于RAM内存的带宽，而对于较小的尺寸，所有的矢量都适合于CPU的缓存。数字显示，尽管使用了低存储量的Runge-Kutta积分器和合并矢量操作，但质量矩阵评估和矢量更新部分几乎消耗了Runge-Kutta阶段所花费的40%的时间。而且尽管对 $\mathcal L_h$ 算子使用了过度积分。对于更简单的微分算子和更昂贵的时间积分器，花费在质量矩阵和矢量更新部分的比例也可以达到70%。如果我们以每秒处理的DoFs和Runge--Kutta阶段计算一个吞吐量数字，我们得到@f[ \text{throughput} =
\frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} = \frac{10258 \cdot 5 \cdot
9.4\,\text{MDoFs}}{205s} = 2360\, \text{MDoFs/s} @f]这个吞吐量数字非常高，因为简单地将一个向量复制到另一个向量的运行速度只有大约10,000 MDoFs/s。

如果我们进入下一个更大的规模，有3770万个DoF，总的模拟时间是2196秒，其中1978秒用于时间步进。L_h算子的运行时间增加了9.3倍（1179秒对127秒），反质量矩阵和向量更新增加了10.3倍（797秒对77.5秒）。运行时间非最佳增长的原因可以追溯到给定硬件上的缓存效应（有40MB的二级缓存和55MB的三级缓存）。虽然不是所有的相关数据都适合940万DoF的缓存（一个向量需要75MB，我们有三个向量加上MatrixFree中的一些额外数据），但还是有能力满足一个半向量的需求。考虑到现代的缓存比天真的最近使用的策略更复杂（在这种情况下，我们几乎没有重复使用，因为数据是以类似流的方式使用的），我们可以假设，在940万DoFs的情况下，确实有相当一部分数据可以从缓存中交付。在更大的情况下，即使有最佳的缓存，也只有不到10%的数据可以放入缓存中，而且会有相关的性能损失。




<a name="Convergenceratesfortheanalyticaltestcase"></a><h3>Convergence rates for the analytical test case</h3>


对于修改后的Lax--Friedrichs通量和测量动量变量的误差，我们得到以下收敛表（密度和能量变量的速率非常相似）。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.373e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.130e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.899e-03</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">5.577e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">7.381e-03</td>
    <td>3.64</td>
    <td align="right">36,864</td>
    <td align="center">2.082e-04</td>
    <td>5.42</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.724e-03</td>
    <td>3.56</td>
    <td align="right">65,536</td>
    <td align="center">3.072e-04</td>
    <td>4.59</td>
    <td align="right">147,456</td>
    <td align="center">2.625e-06</td>
    <td>6.31</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">6.205e-04</td>
    <td>2.92</td>
    <td align="right">262,144</td>
    <td align="center">1.880e-05</td>
    <td>4.03</td>
    <td align="right">589,824</td>
    <td align="center">3.268e-08</td>
    <td>6.33</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">8.279e-05</td>
    <td>2.91</td>
    <td align="right">1,048,576</td>
    <td align="center">1.224e-06</td>
    <td>3.94</td>
    <td align="right">2,359,296</td>
    <td align="center">9.252e-10</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.105e-05</td>
    <td>2.91</td>
    <td align="right">4,194,304</td>
    <td align="center">7.871e-08</td>
    <td>3.96</td>
    <td align="right">9,437,184</td>
    <td align="center">1.369e-10</td>
    <td>2.77</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">1.615e-06</td>
    <td>2.77</td>
    <td align="right">16,777,216</td>
    <td align="center">4.961e-09</td>
    <td>3.99</td>
    <td align="right">37,748,736</td>
    <td align="center">7.091e-11</td>
    <td>0.95</td>
  </tr>
</table> 

如果我们改用Harten-Lax-van Leer通量，结果如下。   <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.339e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.037e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.849e-03</td>
    <td>3.92</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">4.204e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">9.143e-03</td>
    <td>3.31</td>
    <td align="right">36,864</td>
    <td align="center">2.501e-04</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.913e-03</td>
    <td>3.09</td>
    <td align="right">65,536</td>
    <td align="center">3.257e-04</td>
    <td>4.81</td>
    <td align="right">147,456</td>
    <td align="center">3.260e-06</td>
    <td>6.26</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">7.862e-04</td>
    <td>2.64</td>
    <td align="right">262,144</td>
    <td align="center">1.588e-05</td>
    <td>4.36</td>
    <td align="right">589,824</td>
    <td align="center">2.953e-08</td>
    <td>6.79</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">1.137e-04</td>
    <td>2.79</td>
    <td align="right">1,048,576</td>
    <td align="center">9.400e-07</td>
    <td>4.08</td>
    <td align="right">2,359,296</td>
    <td align="center">4.286e-10</td>
    <td>6.11</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.476e-05</td>
    <td>2.95</td>
    <td align="right">4,194,304</td>
    <td align="center">5.799e-08</td>
    <td>4.02</td>
    <td align="right">9,437,184</td>
    <td align="center">2.789e-11</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">2.038e-06</td>
    <td>2.86</td>
    <td align="right">16,777,216</td>
    <td align="center">3.609e-09</td>
    <td>4.01</td>
    <td align="right">37,748,736</td>
    <td align="center">5.730e-11</td>
    <td>-1.04</td>
  </tr>
</table> 

表中显示，我们对两种数值通量都得到了最佳的 $\mathcal O\left(h^{p+1}\right)$ 收敛率。对于 $p=2$ 的Lax--Friedrichs通量，误差略小，但对于 $p=3$ 的情况则相反；在任何情况下，这个测试案例的差异都相对较小。

对于 $p=5$ ，我们在最细的网格上用两种通量达到了 $10^{-11}$ 的舍入精度。还要注意的是，误差是绝对的，域长为 $10^2$ ，所以相对误差低于 $10^{-12}$ 。HLL通量对于最高度数来说要好一些，这是由于Lax--Friedrichs通量的轻微不准确造成的。Lax--Friedrichs通量对离开域的解设置了一个Dirichlet条件，这导致了一个小的人工反射，这在Lax--Friedrichs通量中被凸显出来。除此之外，我们看到数值通量的影响很小，因为元素内部的多项式部分是引起反射的主要动力。当试图用高阶DG设置来接近更具挑战性的设置时，通量的有限影响也会产生影响。以第33步的参数和网格为例，一旦高质部分接近边界，我们就会在两种通量下得到振荡（这反过来会使密度为负值，并使解决方案爆炸），这与低阶有限体积情况不同（ $p=0$ ）。因此，任何导致溶液中出现冲击的情况都需要某种形式的限制性或人工耗散。对于另一种选择，请参见step-69教程程序。




<a name="Resultsforflowinchannelaroundcylinderin2D"></a><h3>Results for flow in channel around cylinder in 2D</h3>


对于渠道中圆柱体周围的流动测试案例，我们需要将第一行代码改为

@code
  constexpr unsigned int testcase = 1;
@endcode

这个测试案例从一个马赫数为0.31的恒定速度和恒定的初始密度的背景场开始；气流必须绕过一个圆柱体形式的障碍物。由于我们对圆柱体壁施加了一个无穿透的条件，最初迎面撞上圆柱体的气流必须重新排列，这就产生了一个大的声波。下面的图片显示了二维情况下5级全局细化时0.1、0.25、0.5和1.0（左上至右下）的压力，使用了102,400个单元，多项式程度为5，所有4个求解变量的自由度为1470万。我们清楚地看到，在时间0.1的第一个快照中，不连续现象在上游方向传播缓慢，在下游方向传播较快。在时间0.25，声波已经到达顶部和底部的墙壁并反射到内部。从下壁和上壁反射波的不同距离，我们可以看到以 GridGenerator::channel_with_cylinder() 为代表的Sch&auml;fer-Turek试验案例的轻微不对称性，圆柱体上方的空间与下方相比要多一些。在后来的时间里，画面更加混乱，到处都是许多声波。

 <table align="center" class="doxtable" style="width:85%">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_010.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_025.png" alt="" width="100%">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_050.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_100.png" alt="" width="100%">
    </td>
  </tr>
</table> 

下一张图片显示了在相同分辨率下，从通道入口向出口看，时间为1.0时的压力仰角图--在这里，我们可以看到大量的反射。在该图中，可以看到两种类型的波。较大振幅的波对应于初始不连续物撞击墙壁时发生的各种反射，而与元素大小相似的小振幅波则对应于数值伪影。它们起源于方案的有限分辨率，并在不连续面通过高阶多项式的元素时出现。这种效应可以通过提高分辨率来治愈。除了这种效应之外，丰富的波浪结构是高阶DG方法的传输精度的结果。

 <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_elevated.jpg" alt="" width="40%"> 

通过2级全局细化，1,600个单元，网格及其在40个MPI进程上的划分情况如下。

 <img src="https://www.dealii.org/images/steps/developer/step-67.grid-owner.png" alt="" width="70%"> 

当我们在40个核心上运行具有4级全局细化的代码时，我们得到以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 3,686,400 ( = 4 [vars] x 25,600 [cells] x 36 [dofs/cell/var] )
Time step size: 7.39876e-05, minimal h: 0.001875, initial transport scaling: 0.00110294


Time:       0, dt:  7.4e-05, norm rho:   4.17e-16, rho * u:  1.629e-16, energy: 1.381e-15
Time:    0.05, dt:  6.3e-05, norm rho:    0.02075, rho * u:    0.03801, energy:   0.08772
Time:     0.1, dt:  5.9e-05, norm rho:    0.02211, rho * u:    0.04515, energy:   0.08953
Time:    0.15, dt:  5.7e-05, norm rho:    0.02261, rho * u:    0.04592, energy:   0.08967
Time:     0.2, dt:  5.8e-05, norm rho:    0.02058, rho * u:    0.04361, energy:   0.08222
Time:    0.25, dt:  5.9e-05, norm rho:    0.01695, rho * u:    0.04203, energy:   0.06873
Time:     0.3, dt:  5.9e-05, norm rho:    0.01653, rho * u:     0.0401, energy:   0.06604
Time:    0.35, dt:  5.7e-05, norm rho:    0.01774, rho * u:    0.04264, energy:    0.0706


...


Time:    1.95, dt:  5.8e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05185
Time:       2, dt:  5.7e-05, norm rho:    0.01432, rho * u:    0.03969, energy:   0.04889


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     273.6s    13 |     273.6s |     273.6s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.01112s    35 |    0.0672s |    0.1337s     0 |
| compute transport speed       |      6914 |     5.422s    35 |     15.96s |     29.99s     1 |
| output                        |        41 |     37.24s    35 |      37.3s |     37.37s     0 |
| rk time stepping total        |     34564 |     205.4s     1 |     219.5s |     230.1s    35 |
| rk_stage - integrals L_h      |    172820 |     153.6s     1 |     164.9s |     175.6s    27 |
| rk_stage - inv mass + vec upd |    172820 |     47.13s    13 |     53.09s |     64.05s    33 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



这里显示的各种数量的规范是对背景场（即初始条件）的偏差 $\rho'$ 、 $(\rho u)'$ 和 $E'$ 。运行时间的分布总体上与之前的测试案例相似。唯一略有不同的是，与反质量矩阵和矢量更新相比，在 $\mathcal L_h$ 中花费的时间比例较大。这是因为几何体是变形的，无矩阵框架需要从内存中加载额外的几何体阵列，这些阵列在仿生网格的情况下是被压缩的。

将全局细化的数量增加到5，输出就变成了。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 14,745,600 ( = 4 [vars] x 102,400 [cells] x 36 [dofs/cell/var] )


...


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |      2693s    32 |      2693s |      2693s    23 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.04537s    32 |     0.173s |    0.3489s     0 |
| compute transport speed       |     13858 |     40.75s    32 |     85.99s |     149.8s     0 |
| output                        |        41 |     153.8s    32 |     153.9s |     154.1s     0 |
| rk time stepping total        |     69284 |      2386s     0 |      2450s |      2496s    32 |
| rk_stage - integrals L_h      |    346420 |      1365s    32 |      1574s |      1718s    19 |
| rk_stage - inv mass + vec upd |    346420 |     722.5s    10 |     870.7s |      1125s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



对性能的影响与分析性测试案例相似--理论上，计算时间应该增加8倍，但我们实际上看到时间步骤增加了11倍（219.5秒对2450秒）。这可以追溯到缓存，小的情况下大多适合缓存。一个有趣的效果，是典型的本地通信（积分 $\mathcal L_h$ ）和全局通信（计算运输速度）混合的程序，有一些负载不平衡，可以通过查看分别遇到不同阶段的最小和最大时间的MPI等级来观察。级别0报告了 "rk时间步进总数 "部分的最快吞吐量。同时，对于 "计算传输速度 "部分，它似乎是最慢的，几乎比平均水平慢了2倍，与较快的等级相比几乎是4倍。由于后者涉及到全局通信，我们可以将这部分的缓慢归因于本地Runge--Kutta阶段在这个等级上推进得更快，需要等到其他处理器跟上。在这一点上，人们可以怀疑这种不平衡的原因。在所有的MPI进程中，单元格的数量几乎是相同的。然而，无矩阵框架在位于通道出口处的仿生和笛卡尔单元上速度更快，较低的MPI等级被分配到这些单元。另一方面，报告Runga--Kutta阶段最高运行时间的等级32拥有靠近圆柱体的弯曲单元，对于这些单元不可能有数据压缩。为了提高吞吐量，我们可以在划分 parallel::distributed::Triangulation 对象时给不同的单元类型分配不同的权重，甚至可以测量几个时间步骤的运行时间，然后尝试重新平衡。

对于1470万DoFs的测试案例，在346000个Runge--Kutta阶段中，每个Runge--Kutta阶段的吞吐量可以计算到2085 MDoFs/s，比上面报告的2360 MDoFs/s的笛卡尔网格吞吐量略慢。

最后，如果我们增加一个额外的细化，我们会记录以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 58,982,400 ( = 4 [vars] x 409,600 [cells] x 36 [dofs/cell/var] )


...


Time:    1.95, dt:  1.4e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05183
Time:       2, dt:  1.4e-05, norm rho:    0.01431, rho * u:    0.03969, energy:   0.04887


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 2.166e+04s    26 | 2.166e+04s | 2.166e+04s    24 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.1758s    30 |     0.672s |     1.376s     1 |
| compute transport speed       |     27748 |     321.3s    34 |     678.8s |      1202s     1 |
| output                        |        41 |     616.3s    32 |     616.4s |     616.4s    34 |
| rk time stepping total        |    138733 | 1.983e+04s     1 | 2.036e+04s | 2.072e+04s    34 |
| rk_stage - integrals L_h      |    693665 | 1.052e+04s    32 | 1.248e+04s | 1.387e+04s    19 |
| rk_stage - inv mass + vec upd |    693665 |      6404s    10 |      7868s | 1.018e+04s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



rk时间步数总数 "部分对应的是2010 MDoFs/s的吞吐量。执行139k时间步长的总体运行时间是20k秒（5.7小时）或每秒7个时间步长--对于有近6000万个未知数来说还不错。通过在计算中添加更多的内核，可以实现更多的吞吐量。




<a name="Resultsforflowinchannelaroundcylinderin3D"></a><h3>Results for flow in channel around cylinder in 3D</h3>


将通道测试案例切换到3D，并进行3次全局细化，输出结果是

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 221,184,000 ( = 5 [vars] x 204,800 [cells] x 216 [dofs/cell/var] )


...


Time:    1.95, dt:  0.00011, norm rho:    0.01131, rho * u:    0.03056, energy:   0.04091
Time:       2, dt:  0.00011, norm rho:     0.0119, rho * u:    0.03142, energy:   0.04425


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 1.734e+04s     4 | 1.734e+04s | 1.734e+04s    38 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.6551s    34 |     3.216s |     7.281s     0 |
| compute transport speed       |      3546 |       160s    34 |     393.2s |     776.9s     0 |
| output                        |        41 |      1350s    34 |      1353s |      1357s     0 |
| rk time stepping total        |     17723 | 1.519e+04s     0 | 1.558e+04s | 1.582e+04s    34 |
| rk_stage - integrals L_h      |     88615 | 1.005e+04s    32 | 1.126e+04s |  1.23e+04s    11 |
| rk_stage - inv mass + vec upd |     88615 |      3056s    11 |      4322s |      5759s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



物理原理与二维情况类似，由于引力的作用，在Z方向有轻微的运动。在这种情况下，每个Runge-Kutta阶段的吞吐量为

@f[
\text{throughput} = \frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} =
\frac{17723 \cdot 5 \cdot 221.2\,\text{M}}{15580s} = 1258\, \text{MDoFs/s}.


@f]



吞吐量低于二维，因为 $\mathcal L_h$ 项的计算更加昂贵。这是由于 "度+2 "点的过度积分和较大比例的面积分（更差的体积-表面比率），以及更昂贵的通量计算。如果我们只考虑反质量矩阵和矢量更新部分，我们记录到等熵涡旋的二维案例的吞吐量为4857 MDoFs/s，有3770万个未知数，而三维案例的运行速度为4535 MDoFs/s。性能是相似的，因为这两种情况实际上都受到内存带宽的限制。

如果我们进行四级全局细化，我们需要增加进程的数量以在内存中容纳所有的东西--在这种情况下，计算需要大约350GB的RAM内存。另外，通过增加额外的资源，完成35k个时间步骤所需的时间也变得更容易忍受。因此，我们使用了6个节点，每个节点有40个核心，从而使计算有240个MPI进程。

@code
Running with 240 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 1,769,472,000 ( = 5 [vars] x 1,638,400 [cells] x 216 [dofs/cell/var] )


...


Time:    1.95, dt:  5.6e-05, norm rho:    0.01129, rho * u:     0.0306, energy:   0.04086
Time:       2, dt:  5.6e-05, norm rho:    0.01189, rho * u:    0.03145, energy:   0.04417


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 5.396e+04s   151 | 5.396e+04s | 5.396e+04s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |     2.632s   178 |     7.221s |     16.56s     0 |
| compute transport speed       |      7072 |       714s   193 |      1553s |      3351s     0 |
| output                        |        41 |      8065s   176 |      8070s |      8079s     0 |
| rk time stepping total        |     35350 |  4.25e+04s     0 |  4.43e+04s | 4.515e+04s   193 |
| rk_stage - integrals L_h      |    176750 | 2.936e+04s   134 | 3.222e+04s |  3.67e+04s    99 |
| rk_stage - inv mass + vec upd |    176750 |      7004s    99 | 1.207e+04s |  1.55e+04s   132 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

这个模拟有近20亿个未知数--确实是一个相当大的计算量，而每个时间步长仍然只需要大约1.5秒。




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这里介绍的代码可以直接扩展到自适应网格，给定适当的指标来设置细化标志。在声波方程的背景下，类似求解器的大规模适应性已经由<a href="https://github.com/kronbichler/exwave">exwave
project</a>实现。然而，在目前的情况下，自适应性的好处往往只限于靠近声波起源的早期时间和效果，因为波最终会反射和衍射。这就导致了到处都是陡峭的梯度，类似于湍流，以及或多或少的全局细化网格。

我们在结果部分没有讨论的另一个话题是不同时间积分方案的比较。该程序提供了四种低存储量的Runga--Kutta积分器的变体，每一种都有轻微不同的精度和稳定性行为。在这里实现的方案中，高阶方案提供了额外的精度，但在违反CFL条件之前，每级步长的效率略低。一个有趣的扩展是将这里提出的低存储变体与标准的Runge--Kutta积分器进行比较，或者使用与质量矩阵运算分开运行的矢量运算，并比较性能。




<a name="Moreadvancednumericalfluxfunctionsandskewsymmetricformulations"></a><h4>More advanced numerical flux functions and skew-symmetric formulations</h4>


正如介绍中提到的，本程序中采用的修改的Lax--Friedrichs通量和HLL通量只是文献中关于欧拉方程的大量数值通量中的两个变种。一个例子是HLLC通量（Harten-Lax-van Leer-Contact）通量，它增加了HLL通量或Roe通量中缺少的稀疏波效应。正如介绍中提到的，数值通量对高阶DG方案的影响是有争议的（与低阶离散的情况不同）。

为了提高求解器的稳定性，一个相关的改进是也要考虑空间积分项。上面使用的相当幼稚的实现方式的一个缺点是，原始欧拉方程的能量守恒（在没有冲击的情况下）只适用于离散化误差。如果解决方案的分辨率不足，离散化误差会引起数值能量的增加，并最终导致离散化的不稳定。这是因为欧拉方程中的项的不精确数值积分，其中包含有理非线性和弯曲单元的高阶内容。摆脱这种困境的方法是所谓的倾斜对称公式，见 @cite Gassner2013 的一个简单变体。倾斜对称意味着在弱式中切换解 $\mathbf{w}$ 和检验函数 $\mathbf{v}$ 的作用，除了一些边界项外，产生原始量的精确负值。在离散设置中，挑战在于当积分只被近似计算时也要保持这种倾斜对称性（在连续情况下，倾斜对称性是部分积分的结果）。偏斜对称的数值方案平衡了保守形式的空间导数  $(\nabla \mathbf v, \mathbf{F}(\mathbf w))_{K}$  和对流形式的贡献  $(\mathbf v, \tilde{\mathbf{F}}(\mathbf w)\nabla
\mathbf{w})_{K}$  ，对于某些  $\tilde{\mathbf{F}}$  。准确的条款取决于方程和积分公式，在某些情况下可以通过特殊的倾斜对称有限差分方案来理解。

要想开始，有兴趣的读者可以看看https://github.com/kronbichler/advection_miniapp，其中用deal.II对一个简单的平流方程实现了倾斜对称的DG公式。

<a name="Equippingthecodeforsupersoniccalculations"></a><h4>Equipping the code for supersonic calculations</h4>


正如介绍中提到的，欧拉方程的解随着马赫数的增加而产生冲击，这需要额外的机制来稳定方案，例如限制器的形式。除了实际实施限制器或人工粘性方法外，主要的挑战是如何平衡计算，因为在有问题的单元中限制震荡所涉及的额外计算会使它们比没有限制的普通DG单元更昂贵。此外，更好地应对不连续情况的额外数值通量也是一种选择。

对于超音速流动来说，有一个因素也是必要的，那就是适当的边界条件。与介绍中讨论的并在程序中实现的亚音速流出边界相反，所有的特性都是超音速流出边界的外在表现，所以我们不想规定任何外部数据。

@f[
\mathbf{w}^+ = \mathbf{w}^- = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E^-\end{pmatrix} \quad
 \text{(Neumann)}.


@f]



在代码中，我们将简单地添加额外的语句

@code
            else if (supersonic_outflow_boundaries.find(boundary_id) !=
                     supersonic_outflow_boundaries.end())
              {
                w_p        = w_m;
                at_outflow = true;
              }
@endcode

在 "local_apply_boundary_face() "函数中。

<a name="ExtensiontothelinearizedEulerequations"></a><h4>Extension to the linearized Euler equations</h4>


当对欧拉解的兴趣主要在于声波的传播时，围绕一个背景状态，即一个给定的密度、速度和能量（或压力）场，将欧拉方程线性化，只计算针对这些场的变化，往往是合理的。这就是航空声学的广泛领域的设置。即使有时分辨率要求大大降低，但由于线性化引起了额外的条款，实施起来就变得有些复杂了。从代码的角度来看，在算子评估中，我们还需要为代码配备要线性化的状态。这一信息可以由分析函数（根据正交点的位置进行评估）或由类似于解决方案的矢量提供。基于该矢量，我们将创建一个额外的FEEvaluation对象，从中读取并提供正交点的场值。如果背景速度为零，密度为常数，线性化的欧拉方程进一步简化，可以等效地写成声波方程的形式。

在声音传播的背景下，一个挑战往往是边界条件的定义，因为计算域需要是有限的，而实际模拟往往跨越无限的（或至少大得多）物理域。传统的Dirichlet或Neumann边界条件会引起声波的反射，最终传播到感兴趣的区域，破坏了解决方案。因此，各种非反射边界条件或海绵层的变体，通常以<a
href="https://en.wikipedia.org/wiki/Perfectly_matched_layer">perfectly
matched layers</a>的形式出现--其中解决方案被阻尼，没有反射

--是很常见的。




<a name="ExtensiontothecompressibleNavierStokesequations"></a><h4>Extension to the compressible Navier-Stokes equations</h4>


如 @cite FehnWallKronbichler2019 所述，本教程程序中的求解器也可以通过添加粘性项扩展到可压缩的Navier-Stokes方程。为了尽量保持这里获得的性能，尽管有额外的椭圆项的成本，例如通过内部惩罚方法，我们可以像步骤59的教程程序一样，将基础从FE_DGQ切换到FE_DGQHermite。




<a name="Usingcellcentricloopsandsharedmemory"></a><h4>Using cell-centric loops and shared memory</h4>


在本教程中，我们使用了以面为中心的循环。在这里，单元和面的积分在不同的循环中处理，导致对结果向量的多次写入访问，这在现代硬件上是比较昂贵的，因为写入操作通常也会导致隐含的读操作。另一方面，以元素为中心的循环是在处理一个单元的同时直接处理其所有的2d面。虽然这种循环意味着通量必须计算两次（对于一个内部面的每一面），但结果向量只需访问一次的事实--以及由此产生的算法没有竞赛条件，因此完全适合共享内存的事实--已经带来了性能的提升。如果你对这些高级主题感兴趣，你可以看一下步骤76，在那里我们对本教程进行了修改，以便我们能够使用这些功能。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-67.cc"
*/
