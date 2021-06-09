/**
@page step_32 The step-32 tutorial program
This tutorial depends on step-31, step-55.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Usingtherightpressure"> Using the "right" pressure </a>
        <li><a href="#Thescalingofdiscretizedequations"> The scaling of discretized equations </a>
        <li><a href="#ChangestotheStokespreconditionerandsolver"> Changes to the Stokes preconditioner and solver </a>
        <li><a href="#Changestotheartificialviscositystabilization"> Changes to the artificial viscosity stabilization </a>
        <li><a href="#LocallyconservativeStokesdiscretization"> Locally conservative Stokes discretization </a>
        <li><a href="#Higherordermappingsforcurvedboundaries"> Higher order mappings for curved boundaries </a>
        <li><a href="#Parallelizationonclusters"> Parallelization on clusters </a>
        <li><a href="#Parallelizationwithinindividualnodesofacluster"> Parallelization within individual nodes of a cluster </a>
        <li><a href="#Thetestcase"> The testcase </a>
        <li><a href="#Implementationdetails"> Implementation details </a>
        <li><a href="#Outlook"> Outlook </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#PreconditioningtheStokessystem">Preconditioning the Stokes system</a>
        <li><a href="#Definitionofassemblydatastructures">Definition of assembly data structures</a>
        <li><a href="#ThecodeBoussinesqFlowProblemcodeclasstemplate">The <code>BoussinesqFlowProblem</code> class template</a>
        <li><a href="#BoussinesqFlowProblemclassimplementation">BoussinesqFlowProblem class implementation</a>
      <ul>
        <li><a href="#BoussinesqFlowProblemParameters">BoussinesqFlowProblem::Parameters</a>
        <li><a href="#BoussinesqFlowProblemBoussinesqFlowProblem">BoussinesqFlowProblem::BoussinesqFlowProblem</a>
        <li><a href="#TheBoussinesqFlowProblemhelperfunctions">The BoussinesqFlowProblem helper functions</a>
      <ul>
        <li><a href="#BoussinesqFlowProblemget_maximal_velocity">BoussinesqFlowProblem::get_maximal_velocity</a>
        <li><a href="#BoussinesqFlowProblemget_cfl_number">BoussinesqFlowProblem::get_cfl_number</a>
        <li><a href="#BoussinesqFlowProblemget_entropy_variation">BoussinesqFlowProblem::get_entropy_variation</a>
        <li><a href="#BoussinesqFlowProblemget_extrapolated_temperature_range">BoussinesqFlowProblem::get_extrapolated_temperature_range</a>
        <li><a href="#BoussinesqFlowProblemcompute_viscosity">BoussinesqFlowProblem::compute_viscosity</a>
      </ul>
        <li><a href="#TheBoussinesqFlowProblemsetupfunctions">The BoussinesqFlowProblem setup functions</a>
        <li><a href="#TheBoussinesqFlowProblemassemblyfunctions">The BoussinesqFlowProblem assembly functions</a>
      <ul>
        <li><a href="#Stokespreconditionerassembly">Stokes preconditioner assembly</a>
        <li><a href="#Stokessystemassembly">Stokes system assembly</a>
        <li><a href="#Temperaturematrixassembly">Temperature matrix assembly</a>
        <li><a href="#Temperaturerighthandsideassembly">Temperature right hand side assembly</a>
      </ul>
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
        <li><a href="#Comparisonofresultswithstep31">Comparison of results with \step-31</a>
        <li><a href="#Resultsfora2dcircularshelltestcase">Results for a 2d circular shell testcase</a>
        <li><a href="#Resultsfora3dsphericalshelltestcase">Results for a 3d spherical shell testcase</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-32/doc/intro.dox

 <br> 

<i>This program was contributed by Martin Kronbichler, Wolfgang
Bangerth, and Timo Heister.


This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology; and in a continuation by the National Science
Foundation under Award No. EAR-0949446 and The University of California
&ndash; Davis. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation, The
California Institute of Technology, or of The University of California
&ndash; Davis.


The work discussed here is also presented in the following publication:
<b>
  M. Kronbichler, T. Heister, W. Bangerth:
  <i>High Accuracy Mantle Convection Simulation through Modern Numerical
  Methods</i><b>
  M. Kronbichler, T. Heister, W. Bangerth:
  <i>High Accuracy Mantle Convection Simulation through Modern Numerical
  Methods</i>, Geophysical Journal International, 2012, 191, 12-29.
  <a href="http://dx.doi.org/10.1111/j.1365-246X.2012.05609.x">[DOI]</a><i>High Accuracy Mantle Convection Simulation through Modern Numerical
  Methods</i>, Geophysical Journal International, 2012, 191, 12-29.
  <a href="http://dx.doi.org/10.1111/j.1365-246X.2012.05609.x">[DOI]</a>
</b><a href="http://dx.doi.org/10.1111/j.1365-246X.2012.05609.x">[DOI]</a>
</b>


The continuation of development of this program has led to the much larger open
source code <i>ASPECT</i><i>ASPECT</i> (see http://aspect.geodynamics.org/) which is much
more flexible in solving many kinds of related problems.
</i>


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序所做的事情与step-31已经做的差不多：它解决了描述温度不平衡的流体运动的Boussinesq方程。因此，我们在step-31中描述的所有方程仍然成立：我们使用相同的有限元方案、相同的时间步进算法和或多或少相同的温度平流-扩散方程的稳定方法来解决相同的一般偏微分方程（只做了些许修改，以适应问题设置的更多现实性）。因此，你可能首先要了解那个程序和它的实现，然后再研究当前的程序。

step-31和当前程序的不同之处在于，在这里，我们想以%并行的方式做事，既利用集群中许多机器的可用性（基于MPI的并行化），也利用一台机器中的许多处理器核心（基于线程的并行化）。因此，本程序的主要工作是引入必要的变化，以利用这些%并行计算资源的可用性。在这方面，它建立在第40步程序的基础上，该程序首先为大部分的%并行功能介绍了必要的类，而第55步则展示了如何为一个矢量值的问题做这件事。

除了这些变化之外，我们还使用了一个略微不同的预处理程序，而且我们将不得不做出一些改变，这与我们在这里想要解决一个<i>realistic</i>问题，而不是一个模型问题有关。特别是后者，将要求我们考虑比例问题，以及所考虑的方程中所有这些参数和系数的实际含义。我们将首先讨论影响数学公式和求解器结构变化的问题，然后讨论如何将事情并行化，最后讨论我们将考虑的实际测试案例。




<a name="Usingtherightpressure"></a><h3> Using the "right" pressure </h3>


在步骤31中，我们对速度和压力场使用了以下斯托克斯模型。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&


  -\rho \; \beta \; T \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0.


@f}

第一个等式的右手边显得有点无动于衷。事情其实应该是这样的。我们需要作用在流体上的外力，我们假设这些外力只是由重力给出的。在目前的情况下，我们假设流体确实为了这个重力的目的而轻微膨胀，但还不足以让我们需要修改不可压缩性条件（第二个方程）。这意味着，为了右手边的目的，我们可以假设 $\rho=\rho(T)$  。一个可能不完全合理的假设是，我们可以假设密度作为温度的函数的变化很小，导致形式为 $\rho(T) = \rho_{\text{ref}}
[1-\beta(T-T_{\text{ref}})]$  的表达，即在参考温度下密度等于 $\rho_{\text{ref}}$ ，并且随着温度的升高（随着材料的膨胀）线性下降。然后，力平衡方程看起来正确地写成这样。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&
  \rho_{\text{ref}} [1-\beta(T-T_{\text{ref}})] \mathbf{g}.


@f}

现在注意到，引力是由重力势产生的，如 $\mathbf g=-\nabla \varphi$  ，因此我们可以将其重新写成如下。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&


  -\rho_{\text{ref}} \; \beta\; T\; \mathbf{g}


  -\rho_{\text{ref}} [1+\beta T_{\text{ref}}] \nabla\varphi.


@f}

右边的第二个项是与时间无关的，因此我们可以引入一个新的 "动态 "压力 $p_{\text{dyn}}=p+\rho_{\text{ref}}
[1+\beta T_{\text{ref}}] \varphi=p_{\text{total}}-p_{\text{static}}$ ，用它来表示斯托克斯方程。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p_{\text{dyn}} &=&


  -\rho_{\text{ref}} \; \beta \; T \; \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0.


@f}

这正是我们在第31步中使用的形式，这样做是合适的，因为流体流动的所有变化只由温度差异导致的动态压力驱动。(换句话说。任何因取标量场的梯度而导致的对右手边的贡献都对速度场没有影响）。)

另一方面，我们在这里将使用考虑总压力的斯托克斯方程的形式来代替。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&
  \rho(T)\; \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0.


@f}

这有几个好处。

- 这样我们就可以在我们的程序中绘制压力图，它实际上显示的是包括温差影响以及上覆岩石的静压力在内的总压力。由于压力没有进一步出现在任何其他方程中，因此使用一个还是另一个，更多的是口味问题，而不是正确性问题。流动场是完全相同的，但我们得到的压力现在可以与地球物理书籍中给出的数值进行比较，例如，在地幔底部的压力。

- 如果我们想让这个模型更加真实，我们就必须考虑到许多材料参数（如粘度、密度等）不仅取决于温度，而且还取决于<i>total</i>压力。

- 上面的模型假设了一个线性依赖 $\rho(T) = \rho_{\text{ref}}
  [1-\beta(T-T_{\text{ref}})]$ ，并假定 $\beta$ 很小。在实践中，情况可能并非如此。事实上，现实的模型肯定不是线性的，而且 $\beta$ 至少在部分温度范围内也可能不小，因为密度的行为不仅大大取决于热膨胀，而且取决于相变。

- 这样做的最后一个原因将在结果部分讨论，涉及到对我们在这里使用的模型的可能扩展。这与我们在这里使用的温度方程（见下文）不包括包含压力的条款这一事实有关。然而，它应该包括：岩石，像气体一样，在你压缩它的时候会升温。因此，上升的物质以绝热方式冷却，而下沉的冷物质以绝热方式升温。我们在下面进一步讨论这个问题。

 @note  然而，这个程序有一个缺点。在地球上，动压比总压要小几个数量级。如果我们使用上述方程并解决所有变量，例如，4位数的精度，那么我们可能会得到正确的速度和总压力，但如果我们通过从总压力中减去静态部分来计算动态压力，我们将完全没有精度  $p_\text{static}=\rho_{\text{ref}}
[1+\beta T_{\text{ref}}] \varphi$  。例如，如果动压比静压小六个数量级，那么我们就需要将总压解到至少七位数的精度，才能得到任何精确的结果。也就是说，在实践中，这并不是一个限制性因素。




<a name="Thescalingofdiscretizedequations"></a><h3> The scaling of discretized equations </h3>


请记住，我们要解决以下方程组。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&
  \rho(T) \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0,
  \\
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma,


@f}

用适当的边界条件和初始条件加以补充。正如第31步所讨论的，我们将通过在每个时间步长中首先求解斯托克斯问题，然后将温度方程向前移动一个时间间隔来解决这组方程。

本节所考虑的问题是斯托克斯问题：如果我们像往常一样对其进行离散化，我们会得到一个线性系统

@f{eqnarray*}
  M \; X
  =
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{c}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{c}
    F_U \\ 0
  \end{array}\right)
  =
  F


@f}

在这个程序中，我们将用FGMRES求解器来解决这些问题。这个求解器一直迭代到这些线性方程的残差低于某个公差，也就是说，直到

@f[
  \left\|
  \left(\begin{array}{c}
    F_U - A U^{(k)} - B P^{(k)}
    \\
    B^T U^{(k)}
  \end{array}\right)
  \right\|
  < \text{Tol}.


@f]

从物理单位的角度来看，这没有任何意义：这里涉及的量有物理单位，所以残差的第一部分有单位 $\frac{\text{Pa}}{\text{m}}
\text{m}^{\text{dim}}$ （通过考虑术语 $(\nabla \cdot \mathbf v, p)_{\Omega}$ 和考虑压力有单位 $\text{Pa}=\frac{\text{kg}}{\text{m}\;\text{s}^2}$ 以及积分得到的系数 $\text{m}^{\text{dim}}$ 最容易确定），而残差的第二部分有单位 $\frac{\text{m}^{\text{dim}}}{\text{s}}$  。取这个残差向量的常数将得到一个单位为  $\text{m}^{\text{dim}-1} \sqrt{\left(\text{Pa}\right)^2 +
       \left(\frac{\text{m}}{\text{s}}\right)^2}$  的量。很明显，这样做是没有意义的，而且我们不应该惊讶这样做最终会伤害到我们。

那么，为什么这在这里是个问题，而在第31步却不是呢？原因是一切都很平衡：速度是1，压力也是1，粘度是1，域的直径是 $\sqrt{2}$  。结果是，虽然不符合逻辑，但没有发生什么坏事。另一方面，正如我们将在下面解释的那样，这里的事情不会是那么简单的缩放。   $\eta$ 将在 $10^{21}$ 左右，速度在 $10^{-8}$ 的数量级，压力在 $10^8$ 左右，域的直径是 $10^7$ 。换句话说，第一个方程的数量级将是  $\eta\text{div}\varepsilon(\mathbf u) \approx 10^{21} \frac{10^{-8}}{(10^7)^2}
\approx 10^{-1}$  ，而第二个方程将是  $\text{div}{\mathbf u}\approx \frac{10^{-8}}{10^7} \approx 10^{-15}$  左右。那么，这将导致这样的结果：如果求解器想使残差变小，它几乎会完全集中在第一组方程上，因为它们大得多，而忽略描述质量守恒的发散方程。这正是发生的情况：除非我们将公差设置为极小的值，否则所得到的流场肯定不是无发散的。作为一个辅助问题，事实证明，很难找到一个始终有效的公差；在实践中，人们往往最终得到一个公差，在大多数时间步骤中需要30或40次迭代，而在其他一些时间步骤中需要10,000次。

那么，在这样的情况下，数字分析员该怎么做呢？答案是要从根本上入手，首先确保一切在数学上是一致的。在我们的例子中，这意味着如果我们想联合解决斯托克斯方程组，我们必须对它们进行缩放，使它们都有相同的物理尺寸。在我们的例子中，这意味着将第二个方程乘以具有单位 $\frac{\text{Pa}\;\text{s}}{\text{m}}$ 的东西；一种选择是乘以 $\frac{\eta}{L}$ ，其中 $L$ 是我们领域的典型长度尺度（实验表明最好选择羽流的直径&mdash；大约10公里&mdash；而不是领域的直径）。使用 $\eta$ 和 $L$ 的这些%数，这个系数约为 $10^{17}$ 。因此，我们现在对斯托克斯系统得到这个。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&
  \rho(T) \; \mathbf{g},
  \\
  \frac{\eta}{L} \nabla \cdot {\mathbf u} &=& 0.


@f}

这样做的问题是，结果不再是对称的（我们在左下方有 $\frac{\eta}{L} \nabla \cdot$ ，但在右上方没有它的转置算子）。然而，这可以通过引入一个按比例的压力 $\hat p = \frac{L}{\eta}p$ 来解决，我们得到按比例的方程式

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) +
  \nabla \left(\frac{\eta}{L} \hat p\right) &=&
  \rho(T) \; \mathbf{g},
  \\
  \frac{\eta}{L} \nabla \cdot {\mathbf u} &=& 0.


@f}

这现在是对称的。很明显，我们可以很容易地从我们作为这个程序的结果计算的比例压力 $\hat p$ 中恢复原始压力 $p$ 。

在下面的程序中，我们将引入一个与 <code>EquationData::pressure_scaling</code> 相对应的因子，我们将在系统矩阵和预处理程序的装配中使用这个因子。因为这很烦人而且容易出错，我们将在线性系统的解之后立即恢复未标定的压力，也就是说，解矢量的压力分量将立即被取消标定以检索物理压力。由于求解器使用的是我们可以通过推断以前的解来使用一个好的初始猜测，所以我们也要立即对压力进行缩放<i>before</i>求解。




<a name="ChangestotheStokespreconditionerandsolver"></a><h3> Changes to the Stokes preconditioner and solver </h3>


在这个教程程序中，我们应用了步骤31中使用的预处理程序的一个变体。该预处理程序是以块状形式对系统矩阵 $M$ 进行操作，从而使乘积矩阵

@f{eqnarray*}
  P^{-1} M
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
  \end{array}\right)
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)


@f}

其形式是基于Krylov的迭代求解器，如GMRES，可以在几次迭代中解决。然后，我们用基于矢量拉普拉斯矩阵的AMG预处理程序 $\tilde{A}$ 的作用取代了 $A$ 的精确逆，用压力空间上的质量矩阵 $M_p$ 来逼近舒尔补码 $S = B A^{-1} B^T$ ，并编写了一个<tt>InverseMatrix</tt>类，用于实现 $M_p^{-1}\approx S^{-1}$ 对矢量的作用。在InverseMatrix类中，我们使用了带有不完全Cholesky（IC）预处理的CG求解器来进行内部求解。

我们可以观察到，我们仅仅使用了预处理程序的作用来逼近速度逆 $A^{-1}$ （外部GMRES迭代处理了逆的近似特性），而我们对 $M_p^{-1}$ 使用了或多或少的<i>exact</i>逆，由完全收敛的CG解实现。这似乎是不平衡的，但这种疯狂是有系统的：几乎所有的努力都用在了左上角的区块上，我们将AMG预处理程序应用于此，而即使是压力质量矩阵的精确反转也基本上不需要花费什么。因此，如果它能帮助我们在一定程度上减少总的迭代次数，那么这种努力是值得的。

也就是说，尽管求解器对step-31工作得很好，但我们这里的问题有点复杂（细胞是变形的，压力有数量级的变化，我们要为更复杂的物理学提前做计划），所以我们要稍微改变一些东西。

- 对于更复杂的问题，事实证明，仅仅使用单一的AMG V-循环作为预处理器并不总是足够的。外围求解器在大多数时候都能在合理的迭代次数内收敛（例如，少于50次），但偶尔会出现突然需要700次左右的时间步骤。到底发生了什么，很难确定，但这个问题可以通过对左上角的块使用更精确的求解器来避免。因此，我们要使用CG迭代来反转预处理矩阵的左上块，并使用AMG作为CG求解器的预处理。

- 这样做的缺点是，当然，斯托克斯预处理程序变得更加昂贵（比我们只使用单个V型循环时大约昂贵10倍）。我们的策略是这样的：让我们只用V型循环作为预处理程序做多达30次的GMRES迭代，如果没有收敛，那么在这第一轮迭代后得到的斯托克斯解的最佳近似值，并将其作为我们使用具有相当宽松容忍度的完整内部求解器作为预处理程序的迭代的起始猜测。在我们所有的实验中，这只导致了少数额外迭代的收敛。

- 我们需要注意的一点是，当在前置条件器中使用具有宽松容忍度的CG时，那么 $y = \tilde A^{-1} r$ 就不再是 $r$ 的线性函数（当然，如果我们的求解器中具有非常严格的容忍度，或者我们只应用单一的V型循环，它就是如此）。这是一个问题，因为现在我们的预处理程序不再是一个线性算子；换句话说，每次GMRES使用它时，预处理程序看起来都不一样。标准的GMRES求解器无法处理这个问题，导致收敛缓慢甚至崩溃，但F-GMRES变体正是为了处理这种情况而设计的，我们因此使用了它。

- 另一方面，一旦我们确定使用F-GMRES，我们就可以放宽在倒置 $S$ 的预处理时使用的容忍度。在第31步中，我们对 $\tilde S$ 运行了一个预处理的CG方法，直到残差减少了7个数量级。在这里，我们可以再次宽松一些，因为我们知道外部预处理程序不会受到影响。

- 在第31步中，我们使用了一个左边的预处理程序，首先反转预处理矩阵的左上块，然后应用左下块（发散）的，再反转右下块。换句话说，预处理器的应用起到了左下块三角矩阵的作用。另一种选择是使用右预处理器，这里将是右上块三角化，即我们首先反转右下舒尔补码，应用右上（梯度）算子，然后反转椭圆的左上块。在某种程度上，选择哪一个是一个品味的问题。也就是说，在GMRES类型的求解器中，右预处理有一个明显的优势：我们决定是否应该停止迭代的残差是真正的残差，而不是预处理方程的规范。因此，将其与我们通常使用的停止标准，即右手边向量的规范进行比较要简单得多。在编写这段代码时，我们发现上面讨论的缩放问题也使我们难以确定适合于左预处理线性系统的停止准则，因此本程序使用了右预处理器。

- 在第31步中，我们对舒尔补码预处理中的压力质量矩阵和温度系统的解使用了IC（不完全Cholesky）预处理。在这里，我们原则上也可以这样做，但我们确实选择了一个更简单的预处理程序，即两个系统的雅可比预处理程序。这是因为在这里我们的目标是大规模的并行计算，IC/ILU的分解必须在每个处理器上对本地拥有的自由度逐块执行。这意味着，无论如何，预处理程序会变得更像一个雅可比预处理程序，所以我们宁愿直接从这个变体开始。请注意，我们只对有质量矩阵的CG求解器使用Jacobi预处理，无论如何它们都能提供最佳的（<i>h</i>独立的）收敛性，尽管它们通常需要两倍于IC预处理的迭代次数。

最后，让我们指出，在第31步中，我们通过逼近 $-\text{div}(-\eta\Delta)^{-1}\nabla \approx \frac 1{\eta} \mathbf{1}$ 来计算舒尔补数 $S=B A^{-1} B^T$ 。然而现在，我们已经对 $B$ 和 $B^T$ 算子进行了重新缩放。所以 $S$ 现在应该近似于 $-\frac{\eta}{L}\text{div}(-\eta\Delta)^{-1}\nabla \frac{\eta}{L} \approx
\left(\frac{\eta}{L}\right)^2 \frac 1{\eta} \mathbf{1}$  。我们用这个的右手边的离散形式作为我们对 $\tilde S$ 的近似 $S$ 。




<a name="Changestotheartificialviscositystabilization"></a><h3> Changes to the artificial viscosity stabilization </h3>


与第31步类似，我们将使用一个基于方程残差的人工黏度进行稳定。  作为与步骤-31的不同之处，我们将提供两个略有不同的稳定参数的定义。对于 $\alpha=1$ ，我们使用与步骤31相同的定义。

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \nu_1(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
  h_K
  \min\left\{
    1,
    \frac{\|R_1(T)\|_{L^\infty(K)}}{c(\mathbf{u},T)}
  \right\}


@f}

我们从方程的残差 $\|R_1(T)\|_{L^\infty(K)}$ 中计算粘度，在残差较大的区域（陡峭的梯度周围），由与网格大小 $h_K$ 成比例的扩散来限制粘度。这个定义已被证明对给定的情况， $\alpha = 1$ 在step-31中效果很好，但它通常不如 $\alpha=2$ 的扩散有效。对于这种情况，我们选择一个稍微可读的粘度定义。

@f{eqnarray*}
  \nu_2(T)|_K = \min (\nu_h^\mathrm{max}|_K,\nu_h^\mathrm{E}|_K)


@f}

其中第一项又给出了最大耗散量（类似于一阶上风方案）。

@f{eqnarray*}
  \nu^\mathrm{max}_h|_K = \beta h_K \|\mathbf {u}\|_{L^\infty(K)}


@f}

而熵粘度的定义为

@f{eqnarray*}
  \nu^\mathrm{E}_h|_K = c_R \frac{h_K^2 \|R_\mathrm{2,E}(T)\|_{L^\infty(K)}}
  {\|E(T) - \bar{E}(T)\|_{L^\infty(\Omega)} }.


@f}



这个公式在<i>J.-L. Guermond, R. Pasquetti, \&
B. Popov, 2011.  Entropy viscosity method for nonlinear conservation laws, J.
Comput. Phys., 230, 4248--4267.</i>一文中有描述。与 $\alpha = 1$ 的情况相比，残差是由温度熵计算出来的， $T_m$ 是平均温度（我们在计算中选择最高和最低温度之间的平均值），这就得到了以下公式

@f{eqnarray*}
 R_\mathrm{E}(T) = \frac{\partial E(T)}{\partial t} +
    (T-T_\mathrm{m}) \left(\mathbf{u} \cdot \nabla T -  \kappa \nabla^2 T - \gamma\right).


@f}

 $\nu^\mathrm{E}_h|_K$ 公式中的分母被计算为熵与空间平均熵的整体偏差  $\bar{E}(T) =
\int_\Omega E(T) d\mathbf{x}/\int_\Omega d\mathbf{x}$  。如同在步骤31中，我们根据前两个时间层次的温度和速度来评估人工黏度，以避免其定义中的非线性。

上述粘度的定义很简单，但取决于两个参数，即  $\beta$  和  $c_R$  。  对于目前的程序，我们想在 $\alpha =1$ 的情况下对这两个参数更系统地去解决这个问题，使用我们在步骤31的结果部分选择离散化的另外两个参数 $c_k$ 和 $\beta$ 的相同推理。特别是，请记住，我们希望使人工粘度尽可能小，同时保持必要的大。在下文中，让我们描述一下人们可能遵循的一般策略。这里显示的计算是用程序的早期版本完成的，因此你在运行程序时得到的实际数值可能不再与这里显示的数值一致；尽管如此，一般的方法仍然有效，并已被用于寻找程序中实际使用的参数值。

为了了解发生了什么，请注意，下面我们将对973和4273开尔文之间的温度施加边界条件，初始条件也选择在这个范围内；出于这些考虑，我们在没有%内部热源或散热器的情况下运行程序，因此温度应该总是在这个范围内，排除任何%内部振荡。如果最低温度下降到973开尔文以下，那么我们需要通过增加 $\beta$ 或减少 $c_R$ 来增加稳定度。

正如我们在第31步所做的那样，我们首先通过使用 "传统 "公式确定 $\beta$ 的最佳值

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
    h_K,


@f}

我们知道，只要 $\beta$ 足够大，它就是稳定的。在2d中做几百个时间步数（在比程序中显示的网格更粗的网格上，用不同的粘度影响传输速度，从而影响时间步数大小），将产生以下图表。

 <img src="https://www.dealii.org/images/steps/developer/step-32.beta.2d.png" alt=""> 

可以看出， $\beta \le 0.05$ 的数值太小，而 $\beta=0.052$ 似乎是有效的，至少在这里显示的时间范围内。顺便说一句，这里至少有两个问题是人们可能想知道的。首先，当解决方案变得不稳定时，会发生什么？看一下图形输出，我们可以看到，在这些实验所选择的不合理的粗大网格下，大约在 $t=10^{15}$ 秒的时间里，一直向冷的外部边界上升，然后向侧面扩散的热物质羽流开始相互靠近，将中间的冷物质挤出去。这就形成了一个细胞层，流体从两个相对的侧面流入，并向第三个侧面流出，显然，这种情况会在没有足够稳定的情况下产生这些不稳定性。第二：在步骤31中，我们使用了 $\beta=0.015\cdot\text{dim}$ ；为什么这在这里不起作用？这个问题的答案并不完全清楚--稳定参数肯定取决于单元格的形状等因素，在第31步中我们使用的是正方形，而在当前程序中则是梯形。不管具体原因是什么，我们至少有一个 $\beta$ 的值，即2d的0.052，对当前程序有效。在3d中也可以做类似的实验，我们发现 $\beta=0.078$ 是一个很好的选择&mdash; 整齐地引出公式 $\beta=0.026 \cdot \textrm{dim}$  。

有了这个值，我们就可以回到粘度的原始公式 $\nu$ ，并玩弄常数 $c_R$ ，使其尽可能大，以便使 $\nu$ 尽可能小。这样我们就得到了这样的画面。

 <img src="https://www.dealii.org/images/steps/developer/step-32.beta_cr.2d.png" alt=""> 

因此， $c_R=0.1$ 似乎是这里的正确值。虽然这个图形是针对指数 $\alpha=1$ 得到的，但在程序中我们用 $\alpha=2$ 代替，在这种情况下，必须重新调整参数（并观察到 $c_R$ 出现在分子中而不是分母中）。事实证明， $c_R=1$ 与 $\alpha=2$ 一起工作。




<a name="LocallyconservativeStokesdiscretization"></a><h3> Locally conservative Stokes discretization </h3>


Stokes的标准Taylor-Hood离散化，使用 $Q_{k+1}^d
\times Q_k$ 元素，是全局保守的，即 $\int_{\partial\Omega}
\mathbf n \cdot \mathbf u_h = 0$  。这很容易看出：发散方程的弱形式为  $(q_h, \textrm{div}\; \mathbf u_h)=0, \forall
q_h\in Q_h$  。因为压力空间确实包含函数  $q_h=1$  ，所以我们得到

@f{align*}
  0 = (1, \textrm{div}\; \mathbf u_h)_\Omega
  = \int_\Omega \textrm{div}\; \mathbf u_h
  = \int_{\partial\Omega} \mathbf n \cdot \mathbf u_h


@f}

由发散定理决定。这个性质很重要：如果我们想用速度场 $u_h$ 沿途输送其他量（如电流方程中的温度，但也可以是化学物质的浓度或完全是人为的示踪量），那么守恒性质保证所输送的量保持恒定。

也就是说，在有些应用中，这个<i>global</i>属性是不够的。相反，我们希望它在每个单元上都持有<i>locally</i>。这可以通过使用空间 $Q_{k+1}^d \times DGP_k$ 进行离散化来实现，我们用相同程度的完整多项式的<i>discontinuous</i>空间代替压力的张量积多项式 $k$ 空间。(注意，2d中的张量积多项式包含函数 $1, x, y, xy$ ，而完全多项式只包含函数 $1,x,y$ ) 。这个空间对斯托克斯方程来说是稳定的。

因为空间是不连续的，我们现在可以特别选择测试函数  $q_h(\mathbf x)=\chi_K(\mathbf x)$  ，即单元格  $K$  的特征函数。然后我们以类似于上面的方式得到

@f{align*}
  0
  = (q_h, \textrm{div}\; \mathbf u_h)_\Omega
  = (1, \textrm{div}\; \mathbf u_h)_K
  = \int_K \textrm{div}\; \mathbf u_h
  = \int_{\partial K} \mathbf n \cdot \mathbf u_h,


@f}

显示了单元格 $K$ 的保存属性。这显然对每个细胞都是成立的。

使用这种离散化是有充分理由的。如上所述，这个元素保证了每个单元上平流量的守恒。第二个优点是，我们用作预处理的压力质量矩阵代替了Schur补码，成为块状对角线，因此非常容易反转。然而，也有缺点。首先，现在有更多的压力变量，增加了问题的总体规模，尽管这在实践中似乎没有造成太大的影响。但更重要的是，现在每个单元上的发散是零，而以前不是，这并不能保证发散是点状的小。事实上，我们可以很容易地验证，与标准Taylor-Hood离散化相比，这个离散化的 $L_2$ 准则是<i>larger</i>。然而，两者都以相同的速度收敛到零，因为很容易看到 $\|\textrm{div}\; u_h\|=
\|\textrm{div}\; (u-u_h)\|=
\|\textrm{trace}\; \nabla (u-u_h)\|\le
\|\nabla (u-u_h)\|={\cal O}(h^{k+2})$  。因此，并不是先验的，仅仅因为我们现在有更多的自由度，误差就真的小了。

鉴于这些考虑，目前还不清楚应该选择哪种离散化方式。因此，我们把这个问题留给用户，并在输入文件中规定使用哪个参数。




<a name="Higherordermappingsforcurvedboundaries"></a><h3> Higher order mappings for curved boundaries </h3>


在程序中，我们将使用一个球壳作为域。这意味着域的内部和外部边界不再是 "直的"（我们通常指它们是可以用FlatManifold类表示的双线性表面）。相反，它们是弯曲的，如果我们已经使用高阶有限元来计算速度，那么在程序中使用一个弯曲的近似值似乎是谨慎的。因此，我们将引入一个MappingQ类型的成员变量，表示这样的映射（步骤10和步骤11首次引入这样的映射），我们将在与边界相邻的单元的所有计算中使用。由于这只影响到相对较小的一部分单元格，额外的努力并不是很大，我们将对这些单元格使用四分法映射。ls.




<a name="Parallelizationonclusters"></a><h3> Parallelization on clusters </h3>


在三维空间中运行具有显著雷利数的对流代码需要大量的计算；在整个地球模拟的情况下，需要一或几亿个未知数的数量。这显然不能用一台机器来完成（至少在2010年我们开始编写这段代码时不能）。因此，我们需要将其并行化。科学代码在计算机集群的多台机器上的并行化几乎总是使用消息传递接口（MPI）来完成。这个程序也不例外，它遵循了第17步和第18步程序的一般精神，尽管在实践中它更多地借用了第40步，在该步中我们首先介绍了当我们想<i>completely</i>分布所有计算时使用的类和策略，而第55步则展示了如何为 @ref vector_valued  "向量值问题"：包括，例如，将网格分割成若干部分，使每个处理器只存储自己的份额和一些幽灵单元，以及使用任何处理器都不可能有足够的内存在本地保存组合解向量的条目的策略。我们的目标是以合理的可扩展性在数百甚至数千台处理器上运行这段代码。

 @note  即使它有一个较大的数字，步骤40在逻辑上是在当前程序之前。第55步的情况也是如此。在你试图理解我们在这里所做的事情之前，你可能会想看看这些程序。

MPI是一个相当笨拙的编程接口。它是一套半面向对象的函数，虽然人们用它在网络上发送数据，但需要明确地描述数据类型，因为MPI函数坚持以 <code>void*</code> 对象的形式获得数据的地址，而不是通过重载或模板自动推断数据类型。我们已经在第17步和第18步中看到，如何通过将所有必要的通信放到deal.II库中，或者在这些程序中放到PETSc中，来避免几乎所有的MPI。我们将在这里做一些类似的事情：就像第40步和第55步一样，deal.II和底层的p4est库负责分配网格所需的所有通信，而我们将让Trilinos库（以及命名空间TrilinosWrappers中的包装器）处理线性代数组件的并行化问题。我们已经在step-31中使用了Trilinos，在这里也会这样做，不同的是我们将使用它的%并行能力。

Trilinos由大量的包组成，实现了基本的%并行线性代数操作（Epetra包），不同的求解器和预处理包，以及对deal.II不太重要的东西（例如。deal.II的Trilinos接口封装了Trilinos提供的许多与PDE求解器相关的东西，并提供了封装类（在命名空间TrilinosWrappers中），使Trilinos的矩阵、向量、求解器和预处理器类看起来与deal.II自己对这些功能的实现非常相同。然而，与deal.II的类相比，如果我们给它们提供必要的信息，它们可以在%并行中使用。因此，有两个Trilinos类我们必须直接处理（而不是通过包装器），这两个类都是Trilinos的Epetra基本线性代数和工具类库的一部分。   <ul>   <li>  Epetra_Comm类是MPI "通信器 "的抽象，也就是说，它描述了多少台机器和哪些机器可以相互通信。   每个分布式对象，如稀疏矩阵或矢量，我们可能想在不同的机器上存储部分，需要有一个通信器对象来知道有多少部分，在哪里可以找到它们，以及如何访问它们。

  在这个程序中，我们只真正使用了一个通信器对象--基于MPI变量 <code>MPI_COMM_WORLD</code> --它包含了<i>all</i>个一起工作的进程。在 $N$ 机器上启动一个进程，但只在其中的一个子集上存储向量，产生一个只包括这个子集的机器的通信器对象是完全合法的；不过，在这里确实没有令人信服的理由这样做。

 <li>  IndexSet类用于描述一个向量的哪些元素或一个矩阵的哪些行应该驻留在作为通信器一部分的当前机器上。要创建这样一个对象，你需要知道（i）元素或行的总数，（ii）你想在本地存储的元素的索引。我们将在下面的 <code>partitioners</code> 函数中设置这些 <code>BoussinesqFlowProblem::setup_dofs</code> ，然后把它交给我们创建的每个%parallel对象。

  与PETSc不同，Trilinos没有假设矢量的元素需要被分割成连续的小块。至少在原则上，我们可以在一个处理器上存储所有偶数索引的元素，在另一个处理器上存储所有奇数索引的元素。当然，这不是很有效率，但这是可能的。此外，这些分区的元素不一定是相互排斥的。这一点很重要，因为在对解决方案进行后处理时，我们需要访问所有本地相关的或至少是本地活跃的自由度（定义见 @ref distributed 上的模块，以及步骤40中的讨论）。那么Trilinos矢量认为哪些元素是本地拥有的，对我们来说并不重要。我们所关心的是，它在本地存储了我们所需要的那些元素。   </ul> 

还有一些与将网格分布到若干处理器上有关的概念；在尝试理解这个程序之前，你可能想看一下 @ref
distributed 模块和步骤40或步骤55。  程序的其余部分几乎完全不知道我们没有完全在本地存储所有对象的事实。有几个地方我们必须将所有单元的循环限制在本地拥有的单元上，或者我们需要区分只存储本地拥有的元素的向量和存储本地相关的所有元素的向量（见 @ref GlossLocallyRelevantDof "这个词汇表条目"），但总的来说，使程序在%parallel中运行所需的大量繁重工作都很好地隐藏在这个程序赖以建立的库中。在任何情况下，当我们在程序代码中看到这些位置时，我们会对它们进行评论。




<a name="Parallelizationwithinindividualnodesofacluster"></a><h3> Parallelization within individual nodes of a cluster </h3>


使程序并行化的第二个策略是利用这样一个事实，即今天大多数计算机都有一个以上的处理器，它们都可以访问相同的内存。换句话说，在这个模型中，我们不需要明确地说哪块数据在哪里，我们需要的所有数据都可以直接访问，我们要做的就是在可用的处理器之间分割<i>processing</i>这些数据。然后，我们将把它与上述的MPI并行化结合起来，也就是说，我们将让一台机器上的所有处理器一起工作，例如，为这台机器实际 "拥有 "的单元汇集对全局矩阵的局部贡献，而不是为那些被其他机器拥有的单元。我们将把这种策略用于本程序中经常进行的四种操作：组装斯托克斯和温度矩阵，组装形成斯托克斯预处理的矩阵，以及组装温度系统的右手边。

所有这些操作基本上都是这样的：我们需要在 <code>cell-@>subdomain_id()</code> 等于我们机器在用于所有通信的通信器对象中的索引（即 <code>MPI_COMM_WORLD</code>  ，如上所述）的所有单元中循环。我们实际要使用的测试，简明扼要地描述了我们为什么要测试这个条件，是  <code>cell-@>is_locally_owned()</code>  。在每一个这样的单元上，我们需要集合对全局矩阵或向量的局部贡献，然后我们必须将每个单元的贡献复制到全局矩阵或向量中。请注意，第一部分（循环）定义了一个必须发生的迭代器的范围。第二部分，本地贡献的组装是在这个步骤序列中花费大部分CPU时间的事情，也是一个可以在%并行中完成的典型例子：每个单元的贡献完全独立于所有其他单元的贡献。第三部分，复制到全局矩阵中，不能在%parallel中进行，因为我们正在修改一个对象，所以几个线程不能同时读取一个现有的矩阵元素，增加他们的贡献，并将总和写回内存而不产生<a
href="http://en.wikipedia.org/wiki/Race_condition">race condition</a>危险。

deal.II有一个类，正是为这个工作流程而生的。WorkStream，首先在步骤9和步骤13中讨论。它的使用在 @ref threads 模块中也有大量的记录（在 @ref MTWorkStream "WorkStream类 "一节），我们不会在这里重复那里阐述的原理和详细说明，尽管你会想通读这个模块以了解从头开始的空间和每单元数据之间的区别。我只想说，我们需要以下条件。

- 迭代器的范围是我们要处理的那些单元格。这是由FilteredIterator类提供的，它的作用就像deal.II中的其他单元格迭代器一样，只是它跳过了所有不满足特定谓词（即，一个评估为真或假的标准）的单元。在我们的例子中，该谓词是一个单元格是否为本地所有。

- 一个为上面确定的每项任务在每个单元上做工作的函数，即集合对斯托克斯矩阵和预调节器、温度矩阵和温度右侧的局部贡献的函数。这些是下面代码中的 <code>BoussinesqFlowProblem::local_assemble_stokes_system</code> 、 <code>BoussinesqFlowProblem::local_assemble_stokes_preconditioner</code> 、 <code>BoussinesqFlowProblem::local_assemble_temperature_matrix</code> 和 <code>BoussinesqFlowProblem::local_assemble_temperature_rhs</code> 函数。这四个函数都可以有几个实例同时并行运行。

- 将前一个函数的结果复制到全局对象中的函数，并按顺序运行以避免竞赛条件。这些是 <code>BoussinesqFlowProblem::copy_local_to_global_stokes_system</code> 、 <code>BoussinesqFlowProblem::copy_local_to_global_stokes_preconditioner</code> 、 <code>BoussinesqFlowProblem::copy_local_to_global_temperature_matrix</code> 、和 <code>BoussinesqFlowProblem::copy_local_to_global_temperature_rhs</code> 函数。

我们将在实际代码中再评论一些要点，但总的来说，它们的结构应该从  @ref threads  的讨论中清楚。

WorkStream的底层技术识别需要处理的 "任务"（例如，在一个单元上组装本地贡献），并将这些任务自动安排到可用的处理器上。WorkStream通过将迭代器范围分割成合适的小块，自动创建这些任务。

 @note  在每个MPI进程中使用多个线程，只有当你在集群的每个节点上运行的MPI进程少于这台机器上的处理器核心时才有意义。否则，MPI已经让你的处理器很忙了，你不会从使用线程中获得任何额外的速度。例如，如果你的集群节点有8个内核，就像在写这篇文章的时候经常有的那样，如果你的批处理调度程序在每个节点上放8个MPI进程，那么使用线程并不能使程序更快。因此，你可能想在运行之前，要么配置你的deal.II不使用线程，要么将 Utilities::MPI::MPI_InitFinalize 中的线程数设置为1（第三个参数），或者 "export DEAL_II_NUM_THREADS=1"。也就是说，在写这篇文章的时候，我们只用WorkStream类来组装（部分）线性系统，而程序的75%或更多的运行时间是在没有并行化的线性求解器中度过的&mdash;换句话说，我们最好的希望是将剩下的25%并行化。




<a name="Thetestcase"></a><h3> The testcase </h3>


这个程序的设置稍微让人想起我们当初想解决的问题（见步骤31的介绍）：地幔的对流。因此，我们选择了以下数据，所有这些数据在程序中都是以米和秒为单位（国际单位制）出现的，即使我们在这里以其他单位列出它们。然而，我们注意到，这些选择基本上仍然只是示范性的，而不是要形成对地幔对流的完全现实的描述：为此，必须实现更多、更困难的物理学，而且目前这个程序中也缺少其他几个方面。我们将在结果部分再次讨论这个问题，但现在要说明的是，在写这篇文章时，提供真实的描述是正在开发的<i>ASPECT</i>代码的一个目标。

作为提醒，让我们再次说明我们要解决的方程是这些。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) +
  \nabla \left( \frac{\eta}{L} \hat p\right) &=&
  \rho(T) \mathbf{g},
  \\
  \frac{\eta}{L} \nabla \cdot {\mathbf u} &=& 0,
  \\
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma,


@f}

用边界条件和初始条件增强。然后我们必须选择以下数量的数据。   <ul>   <li>  域是一个环形（2D）或一个球壳（3D），其内外半径与地球的半径一致：地球的总半径为6371km，地幔从大约35km的深度开始（就在由<a target="_top"
  href="http://en.wikipedia.org/wiki/Continental_crust">continental</a>和<a
  target="_top" href="http://en.wikipedia.org/wiki/Oceanic_crust">oceanic
  plates</a>组成的固体地球<a target="_top"
  href="http://en.wikipedia.org/wiki/Crust_(geology)">crust</a>之下）到2890km深度（<a target="_top" href="http://en.wikipedia.org/wiki/Outer_core">outer earth
  core</a>开始）。因此半径为 $R_0=(6371-2890)\text{km},
  R_1=(6371-35)\text{km}$  。这个领域是使用 GridGenerator::hyper_shell() 函数方便地生成的。

    <li>  在地壳和地幔的界面，温度在500到900摄氏度之间，而在其底部则是4000摄氏度左右（例如，见<a target="_top"
  href="http://en.wikipedia.org/wiki/Mantle_(geology)">this Wikipedia
  entry</a>）。因此，在开尔文中，我们选择 $T_0=(4000+273)\text{K}$  ， $T_1=(500+273)\text{K}$ 作为内外边缘的边界条件。

  除此以外，我们还必须为温度场指定一些初始条件。由于已经持续了40多亿年的对流，地球的真实温度场是相当复杂的--事实上，我们正是想通过这样的程序来探索这种温度分布的特性。因此，我们在这里并没有什么有用的东西可以提供，但是我们可以希望，如果我们从一些东西开始，让事情运行一段时间，确切的初始条件就不再那么重要了&mdash; 事实上，通过查看<a href="#Results">results section
  below</a>中显示的图片就可以看出。我们在这里使用的初始温度场是由@f{align*}
    s &= \frac{\|\mathbf x\|-R_0}{R_1-R_0}, \\
    \varphi &= \arctan \frac{y}{x}, \\
    \tau &= s + \frac 15 s(1-s) \sin(6\varphi) q(z), \\
    T(\mathbf x) &= T_0(1-\tau) + T_1\tau,
  @f}给出半径的。

  其中@f{align*}
    q(z) = \left\{
    \begin{array}{ll}
      1 & \text{in 2d} \\
      \max\{0, \cos(\pi |z/R_1|)\} & \text{in 3d}
    \end{array}
    \right. .
  @f}

  这个复杂的函数本质上是内部和外部温度之间的线性轮廓的扰动。在2D中，函数 $\tau=\tau(\mathbf x)$ 看起来是这样的（我从<a
  href="http://www.wolframalpha.com/input/?i=plot+%28sqrt%28x^2%2By^2%29%2B0.2*%28sqrt%28x^2%2By^2%29*%281-sqrt%28x^2%2By^2%29%29*sin%286*atan2%28x%2Cy%29%29%29%2C+x%3D-1+to+1%2C+y%3D-1+to+1">this
  page</a>得到的图片）。

    <img src="https://www.dealii.org/images/steps/developer/step-32.2d-initial.png" alt=""> 

  这个剖面的重点是，如果我们在 $T(\mathbf x)$ 的定义中使用 $s$ 而不是 $\tau$ ，那么它将只是一个线性内插。   $\tau$ 在内部和外部边界具有与 $s$ 相同的函数值（分别为0和1），但它根据角度和3D中的 $z$ 值将温度曲线拉长一些，产生线性内插场的角度依赖性扰动。我们将在结果部分看到，这是一个完全不实际的温度场（尽管它将会产生有趣的图像），因为温度的平衡状态将是一个几乎恒定的温度，在内部和外部边界有边界层。

    <li>  温度方程的右边包含了内部加热%的速率  $\gamma$  。地球确实通过几种机制自然升温：放射性衰变、化学分离（较重的元素沉到底部，较轻的元素升到顶部；逆流耗散的能量相当于这一分离过程中的势能损失）；随着地球内部固体核心的增长，液态金属结晶释放热量；以及流体运动时粘性摩擦产生的热量耗散。

  化学分离很难建模，因为它需要将地幔物质建模为多个相；它也是一个相对较小的效应。结晶热就更难了，因为它只限于温度和压力允许相变的区域，也就是一个不连续的过程。鉴于对这两种现象进行建模的困难，我们将忽略它们。

  另外两个很容易处理，考虑到我们对温度方程进行缩放的方式，可得出方程@f[
    \gamma(\mathbf x)
     =
     \frac{\rho q+2\eta \varepsilon(\mathbf u):\varepsilon(\mathbf u)}
     {\rho c_p},
  @f]

  其中 $q$ 是 $\frac{W}{kg}$ 中的辐射性加热，列举器中的第二项是粘性摩擦加热。   $\rho$  是密度， $c_p$  是比热。文献中提供了以下近似值。   $c_p=1250 \frac{J}{kg\; K}, q=7.4\cdot 10^{-12}\frac{W}{kg}$  .   其他参数将在本节的其他地方讨论。

  我们在这里忽略了一个内部热源，即绝热加热，这将导致一个令人惊讶的温度场。这一点将在下面的结果部分进行详细评论。

    <li> 对于速度，我们在内半径处选择 $\mathbf{v}=0$ 作为边界条件（即流体粘在地心上），在外半径处选择 $\mathbf{n}\cdot\mathbf{v}=0$ （即流体沿地壳底部切向流动）。这两种情况在物理上都不过分正确：当然，在这两个边界上，流体可以切向流动，但它们会通过与界面另一侧的介质（分别是金属核心和地壳）摩擦而产生剪切应力。这样的情况可以用切向速度的罗宾式边界条件来模拟；在这两种情况下，法向（垂直）速度将为零，尽管即使这样也不完全正确，因为大陆板块也有垂直运动（例如，见<a
  href="http://en.wikipedia.org/wiki/Postglacial_rebound">post-glacial
  rebound</a>的现象）。但是，对切向速度来说，另一侧的介质也在运动，这已经使事情变得更糟了，因此，在最简单的情况下，剪应力将与<i>velocity
  difference</i>成正比，导致边界条件的形式为@f{align*}
    \mathbf{n}\cdot [2\eta \varepsilon(\mathbf v)]
    &=
    s \mathbf{n} \times [\mathbf v - \mathbf v_0],
    \\
    \mathbf{n} \cdot \mathbf v &= 0,
  @f}

  有一个比例常数  $s$  。然而，我们没有走这条路，而是选择了零（棒）和切向流的边界条件。

  顺便提一下，我们也可以在内外边界都选择切向流动条件。然而，这有一个明显的缺点：它使速度不是唯一定义的。原因是所有对应于绕域中心旋转的固体体的速度场 $\hat{\mathbf v}$ 都满足 $\mathrm{div}\;
  \varepsilon(\hat{\mathbf v})=0, \mathrm{div} \;\hat{\mathbf v} = 0$ ，和 $\mathbf{n} \cdot \hat{\mathbf v} = 0$ 。因此，如果 $\mathbf v$ 满足方程和边界条件，那么 $\mathbf v +
  \hat{\mathbf v}$  也满足。这当然不是一个我们想避免的好情况。解决这个问题的传统方法是在边界上选一个任意的点，通过选择速度在那里的所有分量为零，将其称为你的固定点。(在三维空间中，必须选择两个点。)由于这个程序开始时并不打算太现实，我们通过简单地固定整个内部边界的速度来避免这种复杂情况。

    <li> 根据第一顺序，重力矢量总是指向下方。对于像地球这样大的物体来说，问题只是："向上 "是什么地方。天真的答案当然是 "径向向内，向地球中心"。所以在地球表面，我们有@f[
    \mathbf g
    =


    -9.81 \frac{\text{m}}{\text{s}^2} \frac{\mathbf x}{\|\mathbf x\|},
  @f]

  其中 $9.81 \frac{\text{m}}{\text{s}^2}$ 刚好是地球表面的平均重力加速度。但是在地球内部，问题变得有点复杂：例如，在地球的（轨道）中心，你有物质在各个方向上同样用力拉扯，所以 $\mathbf g=0$  。在这两者之间，净力的描述如下：让我们用<a target="_top"
  href="http://en.wikipedia.org/wiki/Potential_energy#Gravitational_potential_energy">gravity
  potential</a>来定义@f[
    \varphi(\mathbf x)
    =
    \int_{\text{earth}}


    -G \frac{\rho(\mathbf y)}{\|\mathbf x-\mathbf y\|}
    \ \text{d}y,
  @f] 。

  那么 $\mathbf g(\mathbf x) = -\nabla \varphi(\mathbf x)$  。如果我们假设密度 $\rho$ 在整个地球上是恒定的，我们可以产生一个重力矢量的分析表达式（不要试图以某种方式整合上述方程--它导致了椭圆积分；一个更简单的方法是注意到 $-\Delta\varphi(\mathbf x) = -4\pi G \rho
  \chi_{\text{earth}}(\mathbf x)$ 并利用径向对称性在所有 ${\mathbb R}^3$ 中解决这个偏微分方程）。   @f[
    \mathbf g(\mathbf x) =
    \left\{
      \begin{array}{ll}


        -\frac{4}{3}\pi G \rho \|\mathbf x\| \frac{\mathbf x}{\|\mathbf x\|}
        & \text{for} \ \|\mathbf x\|<R_1, \\


        -\frac{4}{3}\pi G \rho R^3 \frac{1}{\|\mathbf x\|^2}
        \frac{\mathbf x}{\|\mathbf x\|}
        & \text{for} \ \|\mathbf x\|\ge R_1.
      \end{array}
    \right.
  @f]

  因子 $-\frac{\mathbf x}{\|\mathbf x\|}$ 是指向径向内的单位矢量。当然，在这个问题中，我们只对与地球内部有关的分支感兴趣，即 $\|\mathbf
  x\|<R_1$ 。因此，我们将只考虑表达式@f[
    \mathbf g(\mathbf x) =


        -\frac{4}{3}\pi G \rho \|\mathbf x\| \frac{\mathbf x}{\|\mathbf x\|}
        =


        -\frac{4}{3}\pi G \rho \mathbf x
        =


        - 9.81 \frac{\mathbf x}{R_1} \frac{\text{m}}{\text{s}^2},
  @f] 。

  其中我们可以推断出最后一个表达式，因为我们知道地球在表面的重力（其中 $\|x\|=R_1$  ）。

  我们可以通过整合 $\varphi(r)$ 的微分方程，在密度分布是径向对称的情况下，即 $\rho(\mathbf
  x)=\rho(\|\mathbf x\|)=\rho(r)$ ，推导出一个更一般的表达。在这种情况下，我们将得到@f[
    \varphi(r)
    = 4\pi G \int_0^r \frac 1{s^2} \int_0^s t^2 \rho(t) \; dt \; ds.
  @f] 。




  然而，这有两个问题。(i) 地球不是均匀的，即密度 $\rho$ 取决于 $\mathbf x$ ；事实上它甚至不是一个只取决于半径 $r=\|\mathbf x\|$ 的函数。因此，在现实中，重力并不总是随着我们的深入而减少：因为地心比地幔的密度大得多，重力实际上在地心地幔边界的 $10.7
  \frac{\text{m}}{\text{s}^2}$ 左右达到峰值（见<a
  target="_top" href="http://en.wikipedia.org/wiki/Earth's_gravity">this
  article</a>）。(ii) 密度，以及由此产生的重力矢量，在时间上甚至不是恒定的：毕竟，我们要解决的问题是与时间有关的热的、密度较小的物质的上涌和冷的密度大的物质的下涌。这就导致了重力矢量随空间和时间的变化而变化，并不总是直接指向下方。

  为了不使情况变得更加复杂，我们可以使用这样的近似值：在地幔的内部边界，重力是 $10.7 \frac{\text{m}}{\text{s}^2}$ ，在外部边界，重力是 $9.81 \frac{\text{m}}{\text{s}^2}$ ，在每种情况下都是径向向内的，在两者之间，重力随着离地球中心的径向距离而线性变化。也就是说，实际上稍微现实一点，假设（就像我们下面做的那样）地幔具有恒定的密度也不是那么难。在这种情况下，上面的方程可以被整合，我们得到一个 $\|\mathbf{g}\|$ 的表达式，我们可以拟合常数以匹配地幔顶部和底部的重力，得到@f[
    \|\mathbf{g}\|
    = 1.245\cdot 10^{-6} \frac{1}{\textrm{s}^2} r + 7.714\cdot 10^{13} \frac{\textrm{m}^3}{\textrm{s}^2}\frac{1}{r^2}.
  @f]



    <li> 地幔的密度在空间上有变化，但变化幅度不大。   $\rho_{\text{ref}}=3300 \frac{\text{kg}}{\text{m}^3}$ 是参考温度 $T_{\text{ref}}=293$ 开尔文时的密度的一个相对较好的平均值。

    <li>  热膨胀系数 $\beta$ 也随深度变化（通过其对温度和压力的依赖）。在接近地表的地方，它似乎是 $\beta=45\cdot 10^{-6} \frac 1{\text{K}}$ ，而在地心地幔边界，它可能更接近 $\beta=10\cdot
  10^{-6} \frac 1{\text{K}}$ 。作为一个合理的值，让我们选择 $\beta=2\cdot 10^{-5} \frac 1{\text{K}}$ 。那么密度与温度的关系是 $\rho(T)=[1-\beta(T-T_{\text{ref}})]\rho_{\text{ref}}$  。

    <li>  我们需要指定的第二个至最后一个参数是粘度  $\eta$  。这是一个棘手的问题，因为在地幔典型的温度和压力下，岩石的流动非常缓慢，以至于在实验室里无法准确地确定粘度。那么我们如何知道地幔的粘度呢？最常用的方法是考虑在冰期和冰期之后，冰盾形成和消失的时间尺度比地幔流动的时间尺度短。因此，大陆在冰盾的附加重量下慢慢沉入地幔，而在冰盾再次消失后，它们又慢慢升起（这被称为<a target="_top"
  href="http://en.wikipedia.org/wiki/Postglacial_rebound"><i>postglacial
  rebound</i><i>postglacial
  rebound</i></a>）。通过测量这种反弹的速度，我们可以推断出流向反弹的大陆板块下腾出的区域的物质的粘度。

  使用这种技术，发现 $\eta=10^{21} \text{Pa}\;\text{s}
  = 10^{21} \frac{\text{N}\;\text{s}}{\text{m}^2}
  = 10^{21} \frac{\text{kg}}{\text{m}\;\text{s}}$ 附近的数值是最有可能的，尽管这上面的误差至少是一个数量级的。

  虽然我们将使用这个值，但我们不得不再次提醒，有许多物理原因可以假设这不是正确的值。首先，它确实应该取决于温度：较热的材料很可能比较冷的材料的粘性要小。然而，在现实中，情况甚至更为复杂。地幔中的大多数岩石随着温度和压力的变化而发生相变：根据温度和压力的不同，不同的晶体构型在热力学上比其他的更受青睐，即使地幔的化学成分是均匀的。例如，常见的地幔物质MgSiO<sub>3</sub>在整个地幔的大部分地区以其<a target="_top"
  href="http://en.wikipedia.org/wiki/Perovskite_(structure)">perovskite
  structure</a>的形式存在，但在地幔下部，同样的物质只以<a targe="_top"
  href="http://en.wikipedia.org/wiki/Postperovskite">post-perovskite</a>的形式稳定。显然，为了计算现实的粘度，我们不仅需要知道地幔的确切化学成分和所有物质的粘度，而且还必须计算所有物质在每个正交点的热力学上最稳定的配置。在编写这个程序时，这不是一个可行的建议。

    <li>  我们的最后一个材料参数是热扩散率 $\kappa$  ，其定义为 $\kappa=\frac{k}{\rho c_p}$  ，其中 $k$  是热导率， $\rho$  是密度， $c_p$  是比热。对于这一点，文献表明，它从上地幔的 $0.7$ 左右增加到下地幔的 $1.7 \frac{\text{mm}^2}{\text{s}}$ 左右，尽管确切的数值其实并不那么重要：通过对流的热传输比通过热传导的热传输要重要几个数量级。可能有兴趣知道的是，地幔中最丰富的材料--过氧化物，在超过大约120GPa的压力下似乎变得透明（例如，见J. Badro等人，《科学》305，383-386（2004年））；因此，在下地幔中，通过辐射传输的热传输可能比通过热传导更有效。

  鉴于这些考虑，让我们选择 $\kappa=1 \frac{\text{mm}^2}{\text{s}} =10^{-6} \frac{\text{m}^2}{\text{s}}$ 作为本方案的目的。   </ul> 

所有这些方程数据都在程序中定义在 <code>EquationData</code> 命名空间。当运行时，该程序产生的长期最大速度大约为每年10-40厘米（见下面的结果部分），大约是物理上正确的数量级。我们将设定结束时间为10亿年。

 @note  上述常数和材料参数的选择在很大程度上遵循了G.Schubert和D.L.Turcotte和P.Olson（剑桥，2001）的综合书籍《地球和行星的地幔对流，第一部分》。它包含了关于如何使程序更加真实的广泛讨论。




<a name="Implementationdetails"></a><h3> Implementation details </h3>


与step-31相比，这个程序有一些值得注意的区别。

-  <code>EquationData</code> 命名空间要大得多，这反映了我们现在有更多的物理学需要处理的事实。也就是说，这些额外的物理细节大部分是在这个命名空间的函数中自成一体的，并没有扩散到程序的其他部分。

- 更明显的可见性是，我们把大量的参数放入由ParameterHandler类处理的输入文件中（例如，见步骤29，关于用这个类设置运行时参数文件的方法）。当人们想避免仅仅因为想玩弄一个参数而重新编译程序时，这往往是有意义的（例如，想想确定上面讨论的稳定常数的最佳值的参数研究），特别是考虑到重新编译当前规模的程序需要花费非同小可的时间。为了仅仅概述我们从固定值移入输入文件的参数种类，这里列出了一个典型的 <code>\step-32.prm</code> 文件。   @code
# Listing of Parameters
# ---------------------
# The end time of the simulation in years.
set End time                            = 1e8


# Whether graphical output is to be generated or not. You may not want to get
# graphical output if the number of processors is large.
set Generate graphical output           = false


# The number of adaptive refinement steps performed after initial global
# refinement.
set Initial adaptive refinement         = 1


# The number of global refinement steps performed on the initial coarse mesh,
# before the problem is first solved there.
set Initial global refinement           = 1


# The number of time steps between each generation of graphical output files.
set Time steps between graphical output = 50


# The number of time steps after which the mesh is to be adapted based on
# computed error indicators.
set Time steps between mesh refinement  = 10



subsection Discretization
  # The polynomial degree to use for the velocity variables in the Stokes
  # system.
  set Stokes velocity polynomial degree       = 2


  # The polynomial degree to use for the temperature variable.
  set Temperature polynomial degree           = 2


  # Whether to use a Stokes discretization that is locally conservative at the
  # expense of a larger number of degrees of freedom, or to go with a cheaper
  # discretization that does not locally conserve mass (although it is
  # globally conservative.
  set Use locally conservative discretization = true
end



subsection Stabilization parameters
  # The exponent in the entropy viscosity stabilization.
  set alpha = 2


  # The beta factor in the artificial viscosity stabilization. An appropriate
  # value for 2d is 0.052 and 0.078 for 3d.
  set beta  = 0.078


  # The c_R factor in the entropy viscosity stabilization.
  set c_R   = 0.5
end
  @endcode



- 很明显，有很多变化是与我们想在可能非常多的机器上运行我们的程序这一事实有关的。尽管人们可能会怀疑这需要我们完全重新构建我们的代码，但事实上并非如此（尽管在deal.II中实现大部分功能的类从实现的角度来看肯定非常不同，但这并没有反映在它们的公共接口中）。相反，这些变化大多是微妙的，主类的整体结构几乎没有变化。也就是说，魔鬼在细节中：正确地进行%并行计算，没有死锁，确保正确的数据在正确的地方可用（例如，见关于全分布式向量与有鬼魂元素的向量的讨论），以及避免瓶颈是很困难的，关于这个话题的讨论将出现在本程序中的很多地方。




<a name="Outlook"></a><h3> Outlook </h3>


这是一个教程性的程序。这意味着至少它的大部分重点需要放在演示如何使用deal.II和相关的库上，而不是通过过度关注物理细节来稀释这个教学课程。尽管上面有关于物理参数选择的长篇大论，但程序中专门讨论这个问题的部分实际上是很短的，而且是自成一体的。

也就是说，第31步和目前的第32步都不是偶然出现的，而肯定是作为通向更全面的计划的路标，该计划将模拟地幔的对流。我们把这个代码称为<i>ASPECT</i>（简称<i>Advanced %Solver for Problems in Earth's
ConvecTion</i>）；它的开发是由<a href="http://www.geodynamics.org">Computational Infrastructure in
Geodynamics</a>计划资助的，得到了美国国家科学基金会的支持。关于<i>ASPECT</i>的更多信息可在其<a href="https://aspect.geodynamics.org/">homepage</a>中找到。


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
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/parameter_handler.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/solver_bicgstab.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/solver_gmres.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/block_sparsity_pattern.h>
 * #include <deal.II/lac/trilinos_parallel_block_vector.h>
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * #include <deal.II/lac/trilinos_solver.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/filtered_iterator.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_dgp.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/solution_transfer.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * #include <limits>
 * #include <locale>
 * #include <string>
 * 
 * @endcode
 * 
 * This is the only include file that is new: It introduces the
 * parallel::distributed::SolutionTransfer equivalent of the
 * dealii::SolutionTransfer class to take a solution from on mesh to the next
 * one upon mesh refinement, but in the case of parallel distributed
 * triangulations:
 * 
 * @code
 * #include <deal.II/distributed/solution_transfer.h>
 * 
 * @endcode
 * 
 * The following classes are used in parallel distributed computations and
 * have all already been introduced in step-40:
 * 
 * @code
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/distributed/tria.h>
 * #include <deal.II/distributed/grid_refinement.h>
 * 
 * 
 * @endcode
 * 
 * The next step is like in all previous tutorial programs: We put everything
 * into a namespace of its own and then import the deal.II classes and
 * functions into it:
 * 
 * @code
 * namespace Step32
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * In the following namespace, we define the various pieces of equation data
 * that describe the problem. This corresponds to the various aspects of
 * making the problem at least slightly realistic and that were exhaustively
 * discussed in the description of the testcase in the introduction.
 *   

 * 
 * We start with a few coefficients that have constant values (the comment
 * after the value indicates its physical units):
 * 
 * @code
 *   namespace EquationData
 *   {
 *     constexpr double eta                   = 1e21;    /* Pa s       */
 *     constexpr double kappa                 = 1e-6;    /* m^2 / s    */
 *     constexpr double reference_density     = 3300;    /* kg / m^3   */
 *     constexpr double reference_temperature = 293;     /* K          */
 *     constexpr double expansion_coefficient = 2e-5;    /* 1/K        */
 *     constexpr double specific_heat         = 1250;    /* J / K / kg */
 *     constexpr double radiogenic_heating    = 7.4e-12; /* W / kg     */
 * 
 * 
 *     constexpr double R0 = 6371000. - 2890000.; /* m          */
 *     constexpr double R1 = 6371000. - 35000.;   /* m          */
 * 
 *     constexpr double T0 = 4000 + 273; /* K          */
 *     constexpr double T1 = 700 + 273;  /* K          */
 * 
 * 
 * @endcode
 * 
 * The next set of definitions are for functions that encode the density
 * as a function of temperature, the gravity vector, and the initial
 * values for the temperature. Again, all of these (along with the values
 * they compute) are discussed in the introduction:
 * 
 * @code
 *     double density(const double temperature)
 *     {
 *       return (
 *         reference_density *
 *         (1 - expansion_coefficient * (temperature - reference_temperature)));
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> gravity_vector(const Point<dim> &p)
 *     {
 *       const double r = p.norm();
 *       return -(1.245e-6 * r + 7.714e13 / r / r) * p / r;
 *     }
 * 
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
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  value) const override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     double TemperatureInitialValues<dim>::value(const Point<dim> &p,
 *                                                 const unsigned int) const
 *     {
 *       const double r = p.norm();
 *       const double h = R1 - R0;
 * 
 *       const double s = (r - R0) / h;
 *       const double q =
 *         (dim == 3) ? std::max(0.0, cos(numbers::PI * abs(p(2) / R1))) : 1.0;
 *       const double phi = std::atan2(p(0), p(1));
 *       const double tau = s + 0.2 * s * (1 - s) * std::sin(6 * phi) * q;
 * 
 *       return T0 * (1.0 - tau) + T1 * tau;
 *     }
 * 
 * 
 *     template <int dim>
 *     void
 *     TemperatureInitialValues<dim>::vector_value(const Point<dim> &p,
 *                                                 Vector<double> &  values) const
 *     {
 *       for (unsigned int c = 0; c < this->n_components; ++c)
 *         values(c) = TemperatureInitialValues<dim>::value(p, c);
 *     }
 * 
 * 
 * @endcode
 * 
 * As mentioned in the introduction we need to rescale the pressure to
 * avoid the relative ill-conditioning of the momentum and mass
 * conservation equations. The scaling factor is $\frac{\eta}{L}$ where
 * $L$ was a typical length scale. By experimenting it turns out that a
 * good length scale is the diameter of plumes, which is around 10 km:
 * 
 * @code
 *     constexpr double pressure_scaling = eta / 10000;
 * 
 * @endcode
 * 
 * The final number in this namespace is a constant that denotes the
 * number of seconds per (average, tropical) year. We use this only when
 * generating screen output: internally, all computations of this program
 * happen in SI units (kilogram, meter, seconds) but writing geological
 * times in seconds yields numbers that one can't relate to reality, and
 * so we convert to years using the factor defined here:
 * 
 * @code
 *     const double year_in_seconds = 60 * 60 * 24 * 365.2425;
 * 
 *   } // namespace EquationData
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PreconditioningtheStokessystem"></a> 
 * <h3>Preconditioning the Stokes system</h3>
 * 

 * 
 * This namespace implements the preconditioner. As discussed in the
 * introduction, this preconditioner differs in a number of key portions
 * from the one used in step-31. Specifically, it is a right preconditioner,
 * implementing the matrix
 * @f{align*}
 * \left(\begin{array}{cc}A^{-1} & B^T
 * \\0 & S^{-1}
 * \end{array}\right)
 * @f}
 * where the two inverse matrix operations
 * are approximated by linear solvers or, if the right flag is given to the
 * constructor of this class, by a single AMG V-cycle for the velocity
 * block. The three code blocks of the <code>vmult</code> function implement
 * the multiplications with the three blocks of this preconditioner matrix
 * and should be self explanatory if you have read through step-31 or the
 * discussion of composing solvers in step-20.
 * 
 * @code
 *   namespace LinearSolvers
 *   {
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp>
 *     class BlockSchurPreconditioner : public Subscriptor
 *     {
 *     public:
 *       BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix &S,
 *                                const TrilinosWrappers::BlockSparseMatrix &Spre,
 *                                const PreconditionerTypeMp &Mppreconditioner,
 *                                const PreconditionerTypeA & Apreconditioner,
 *                                const bool                  do_solve_A)
 *         : stokes_matrix(&S)
 *         , stokes_preconditioner_matrix(&Spre)
 *         , mp_preconditioner(Mppreconditioner)
 *         , a_preconditioner(Apreconditioner)
 *         , do_solve_A(do_solve_A)
 *       {}
 * 
 *       void vmult(TrilinosWrappers::MPI::BlockVector &      dst,
 *                  const TrilinosWrappers::MPI::BlockVector &src) const
 *       {
 *         TrilinosWrappers::MPI::Vector utmp(src.block(0));
 * 
 *         {
 *           SolverControl solver_control(5000, 1e-6 * src.block(1).l2_norm());
 * 
 *           SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
 * 
 *           solver.solve(stokes_preconditioner_matrix->block(1, 1),
 *                        dst.block(1),
 *                        src.block(1),
 *                        mp_preconditioner);
 * 
 *           dst.block(1) *= -1.0;
 *         }
 * 
 *         {
 *           stokes_matrix->block(0, 1).vmult(utmp, dst.block(1));
 *           utmp *= -1.0;
 *           utmp.add(src.block(0));
 *         }
 * 
 *         if (do_solve_A == true)
 *           {
 *             SolverControl solver_control(5000, utmp.l2_norm() * 1e-2);
 *             TrilinosWrappers::SolverCG solver(solver_control);
 *             solver.solve(stokes_matrix->block(0, 0),
 *                          dst.block(0),
 *                          utmp,
 *                          a_preconditioner);
 *           }
 *         else
 *           a_preconditioner.vmult(dst.block(0), utmp);
 *       }
 * 
 *     private:
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix>
 *         stokes_matrix;
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix>
 *                                   stokes_preconditioner_matrix;
 *       const PreconditionerTypeMp &mp_preconditioner;
 *       const PreconditionerTypeA & a_preconditioner;
 *       const bool                  do_solve_A;
 *     };
 *   } // namespace LinearSolvers
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Definitionofassemblydatastructures"></a> 
 * <h3>Definition of assembly data structures</h3>
 *   

 * 
 * As described in the introduction, we will use the WorkStream mechanism
 * discussed in the @ref threads module to parallelize operations among the
 * processors of a single machine. The WorkStream class requires that data
 * is passed around in two kinds of data structures, one for scratch data
 * and one to pass data from the assembly function to the function that
 * copies local contributions into global objects.
 *   

 * 
 * The following namespace (and the two sub-namespaces) contains a
 * collection of data structures that serve this purpose, one pair for each
 * of the four operations discussed in the introduction that we will want to
 * parallelize. Each assembly routine gets two sets of data: a Scratch array
 * that collects all the classes and arrays that are used for the
 * calculation of the cell contribution, and a CopyData array that keeps
 * local matrices and vectors which will be written into the global
 * matrix. Whereas CopyData is a container for the final data that is
 * written into the global matrices and vector (and, thus, absolutely
 * necessary), the Scratch arrays are merely there for performance reasons
 * &mdash; it would be much more expensive to set up a FEValues object on
 * each cell, than creating it only once and updating some derivative data.
 *   

 * 
 * Step-31 had four assembly routines: One for the preconditioner matrix of
 * the Stokes system, one for the Stokes matrix and right hand side, one for
 * the temperature matrices and one for the right hand side of the
 * temperature equation. We here organize the scratch arrays and CopyData
 * objects for each of those four assembly components using a
 * <code>struct</code> environment (since we consider these as temporary
 * objects we pass around, rather than classes that implement functionality
 * of their own, though this is a more subjective point of view to
 * distinguish between <code>struct</code>s and <code>class</code>es).
 *   

 * 
 * Regarding the Scratch objects, each struct is equipped with a constructor
 * that creates an @ref FEValues object using the @ref FiniteElement,
 * Quadrature, @ref Mapping (which describes the interpolation of curved
 * boundaries), and @ref UpdateFlags instances. Moreover, we manually
 * implement a copy constructor (since the FEValues class is not copyable by
 * itself), and provide some additional vector fields that are used to hold
 * intermediate data during the computation of local contributions.
 *   

 * 
 * Let us start with the scratch arrays and, specifically, the one used for
 * assembly of the Stokes preconditioner:
 * 
 * @code
 *   namespace Assembly
 *   {
 *     namespace Scratch
 *     {
 *       template <int dim>
 *       struct StokesPreconditioner
 *       {
 *         StokesPreconditioner(const FiniteElement<dim> &stokes_fe,
 *                              const Quadrature<dim> &   stokes_quadrature,
 *                              const Mapping<dim> &      mapping,
 *                              const UpdateFlags         update_flags);
 * 
 *         StokesPreconditioner(const StokesPreconditioner &data);
 * 
 * 
 *         FEValues<dim> stokes_fe_values;
 * 
 *         std::vector<Tensor<2, dim>> grad_phi_u;
 *         std::vector<double>         phi_p;
 *       };
 * 
 *       template <int dim>
 *       StokesPreconditioner<dim>::StokesPreconditioner(
 *         const FiniteElement<dim> &stokes_fe,
 *         const Quadrature<dim> &   stokes_quadrature,
 *         const Mapping<dim> &      mapping,
 *         const UpdateFlags         update_flags)
 *         : stokes_fe_values(mapping, stokes_fe, stokes_quadrature, update_flags)
 *         , grad_phi_u(stokes_fe.n_dofs_per_cell())
 *         , phi_p(stokes_fe.n_dofs_per_cell())
 *       {}
 * 
 * 
 * 
 *       template <int dim>
 *       StokesPreconditioner<dim>::StokesPreconditioner(
 *         const StokesPreconditioner &scratch)
 *         : stokes_fe_values(scratch.stokes_fe_values.get_mapping(),
 *                            scratch.stokes_fe_values.get_fe(),
 *                            scratch.stokes_fe_values.get_quadrature(),
 *                            scratch.stokes_fe_values.get_update_flags())
 *         , grad_phi_u(scratch.grad_phi_u)
 *         , phi_p(scratch.phi_p)
 *       {}
 * 
 * 
 * 
 * @endcode
 * 
 * The next one is the scratch object used for the assembly of the full
 * Stokes system. Observe that we derive the StokesSystem scratch class
 * from the StokesPreconditioner class above. We do this because all the
 * objects that are necessary for the assembly of the preconditioner are
 * also needed for the actual matrix system and right hand side, plus
 * some extra data. This makes the program more compact. Note also that
 * the assembly of the Stokes system and the temperature right hand side
 * further down requires data from temperature and velocity,
 * respectively, so we actually need two FEValues objects for those two
 * cases.
 * 
 * @code
 *       template <int dim>
 *       struct StokesSystem : public StokesPreconditioner<dim>
 *       {
 *         StokesSystem(const FiniteElement<dim> &stokes_fe,
 *                      const Mapping<dim> &      mapping,
 *                      const Quadrature<dim> &   stokes_quadrature,
 *                      const UpdateFlags         stokes_update_flags,
 *                      const FiniteElement<dim> &temperature_fe,
 *                      const UpdateFlags         temperature_update_flags);
 * 
 *         StokesSystem(const StokesSystem<dim> &data);
 * 
 * 
 *         FEValues<dim> temperature_fe_values;
 * 
 *         std::vector<Tensor<1, dim>>          phi_u;
 *         std::vector<SymmetricTensor<2, dim>> grads_phi_u;
 *         std::vector<double>                  div_phi_u;
 * 
 *         std::vector<double> old_temperature_values;
 *       };
 * 
 * 
 *       template <int dim>
 *       StokesSystem<dim>::StokesSystem(
 *         const FiniteElement<dim> &stokes_fe,
 *         const Mapping<dim> &      mapping,
 *         const Quadrature<dim> &   stokes_quadrature,
 *         const UpdateFlags         stokes_update_flags,
 *         const FiniteElement<dim> &temperature_fe,
 *         const UpdateFlags         temperature_update_flags)
 *         : StokesPreconditioner<dim>(stokes_fe,
 *                                     stokes_quadrature,
 *                                     mapping,
 *                                     stokes_update_flags)
 *         , temperature_fe_values(mapping,
 *                                 temperature_fe,
 *                                 stokes_quadrature,
 *                                 temperature_update_flags)
 *         , phi_u(stokes_fe.n_dofs_per_cell())
 *         , grads_phi_u(stokes_fe.n_dofs_per_cell())
 *         , div_phi_u(stokes_fe.n_dofs_per_cell())
 *         , old_temperature_values(stokes_quadrature.size())
 *       {}
 * 
 * 
 *       template <int dim>
 *       StokesSystem<dim>::StokesSystem(const StokesSystem<dim> &scratch)
 *         : StokesPreconditioner<dim>(scratch)
 *         , temperature_fe_values(
 *             scratch.temperature_fe_values.get_mapping(),
 *             scratch.temperature_fe_values.get_fe(),
 *             scratch.temperature_fe_values.get_quadrature(),
 *             scratch.temperature_fe_values.get_update_flags())
 *         , phi_u(scratch.phi_u)
 *         , grads_phi_u(scratch.grads_phi_u)
 *         , div_phi_u(scratch.div_phi_u)
 *         , old_temperature_values(scratch.old_temperature_values)
 *       {}
 * 
 * 
 * @endcode
 * 
 * After defining the objects used in the assembly of the Stokes system,
 * we do the same for the assembly of the matrices necessary for the
 * temperature system. The general structure is very similar:
 * 
 * @code
 *       template <int dim>
 *       struct TemperatureMatrix
 *       {
 *         TemperatureMatrix(const FiniteElement<dim> &temperature_fe,
 *                           const Mapping<dim> &      mapping,
 *                           const Quadrature<dim> &   temperature_quadrature);
 * 
 *         TemperatureMatrix(const TemperatureMatrix &data);
 * 
 * 
 *         FEValues<dim> temperature_fe_values;
 * 
 *         std::vector<double>         phi_T;
 *         std::vector<Tensor<1, dim>> grad_phi_T;
 *       };
 * 
 * 
 *       template <int dim>
 *       TemperatureMatrix<dim>::TemperatureMatrix(
 *         const FiniteElement<dim> &temperature_fe,
 *         const Mapping<dim> &      mapping,
 *         const Quadrature<dim> &   temperature_quadrature)
 *         : temperature_fe_values(mapping,
 *                                 temperature_fe,
 *                                 temperature_quadrature,
 *                                 update_values | update_gradients |
 *                                   update_JxW_values)
 *         , phi_T(temperature_fe.n_dofs_per_cell())
 *         , grad_phi_T(temperature_fe.n_dofs_per_cell())
 *       {}
 * 
 * 
 *       template <int dim>
 *       TemperatureMatrix<dim>::TemperatureMatrix(
 *         const TemperatureMatrix &scratch)
 *         : temperature_fe_values(
 *             scratch.temperature_fe_values.get_mapping(),
 *             scratch.temperature_fe_values.get_fe(),
 *             scratch.temperature_fe_values.get_quadrature(),
 *             scratch.temperature_fe_values.get_update_flags())
 *         , phi_T(scratch.phi_T)
 *         , grad_phi_T(scratch.grad_phi_T)
 *       {}
 * 
 * 
 * @endcode
 * 
 * The final scratch object is used in the assembly of the right hand
 * side of the temperature system. This object is significantly larger
 * than the ones above because a lot more quantities enter the
 * computation of the right hand side of the temperature equation. In
 * particular, the temperature values and gradients of the previous two
 * time steps need to be evaluated at the quadrature points, as well as
 * the velocities and the strain rates (i.e. the symmetric gradients of
 * the velocity) that enter the right hand side as friction heating
 * terms. Despite the number of terms, the following should be rather
 * self explanatory:
 * 
 * @code
 *       template <int dim>
 *       struct TemperatureRHS
 *       {
 *         TemperatureRHS(const FiniteElement<dim> &temperature_fe,
 *                        const FiniteElement<dim> &stokes_fe,
 *                        const Mapping<dim> &      mapping,
 *                        const Quadrature<dim> &   quadrature);
 * 
 *         TemperatureRHS(const TemperatureRHS &data);
 * 
 * 
 *         FEValues<dim> temperature_fe_values;
 *         FEValues<dim> stokes_fe_values;
 * 
 *         std::vector<double>         phi_T;
 *         std::vector<Tensor<1, dim>> grad_phi_T;
 * 
 *         std::vector<Tensor<1, dim>> old_velocity_values;
 *         std::vector<Tensor<1, dim>> old_old_velocity_values;
 * 
 *         std::vector<SymmetricTensor<2, dim>> old_strain_rates;
 *         std::vector<SymmetricTensor<2, dim>> old_old_strain_rates;
 * 
 *         std::vector<double>         old_temperature_values;
 *         std::vector<double>         old_old_temperature_values;
 *         std::vector<Tensor<1, dim>> old_temperature_grads;
 *         std::vector<Tensor<1, dim>> old_old_temperature_grads;
 *         std::vector<double>         old_temperature_laplacians;
 *         std::vector<double>         old_old_temperature_laplacians;
 *       };
 * 
 * 
 *       template <int dim>
 *       TemperatureRHS<dim>::TemperatureRHS(
 *         const FiniteElement<dim> &temperature_fe,
 *         const FiniteElement<dim> &stokes_fe,
 *         const Mapping<dim> &      mapping,
 *         const Quadrature<dim> &   quadrature)
 *         : temperature_fe_values(mapping,
 *                                 temperature_fe,
 *                                 quadrature,
 *                                 update_values | update_gradients |
 *                                   update_hessians | update_quadrature_points |
 *                                   update_JxW_values)
 *         , stokes_fe_values(mapping,
 *                            stokes_fe,
 *                            quadrature,
 *                            update_values | update_gradients)
 *         , phi_T(temperature_fe.n_dofs_per_cell())
 *         , grad_phi_T(temperature_fe.n_dofs_per_cell())
 *         ,
 * 
 *         old_velocity_values(quadrature.size())
 *         , old_old_velocity_values(quadrature.size())
 *         , old_strain_rates(quadrature.size())
 *         , old_old_strain_rates(quadrature.size())
 *         ,
 * 
 *         old_temperature_values(quadrature.size())
 *         , old_old_temperature_values(quadrature.size())
 *         , old_temperature_grads(quadrature.size())
 *         , old_old_temperature_grads(quadrature.size())
 *         , old_temperature_laplacians(quadrature.size())
 *         , old_old_temperature_laplacians(quadrature.size())
 *       {}
 * 
 * 
 *       template <int dim>
 *       TemperatureRHS<dim>::TemperatureRHS(const TemperatureRHS &scratch)
 *         : temperature_fe_values(
 *             scratch.temperature_fe_values.get_mapping(),
 *             scratch.temperature_fe_values.get_fe(),
 *             scratch.temperature_fe_values.get_quadrature(),
 *             scratch.temperature_fe_values.get_update_flags())
 *         , stokes_fe_values(scratch.stokes_fe_values.get_mapping(),
 *                            scratch.stokes_fe_values.get_fe(),
 *                            scratch.stokes_fe_values.get_quadrature(),
 *                            scratch.stokes_fe_values.get_update_flags())
 *         , phi_T(scratch.phi_T)
 *         , grad_phi_T(scratch.grad_phi_T)
 *         ,
 * 
 *         old_velocity_values(scratch.old_velocity_values)
 *         , old_old_velocity_values(scratch.old_old_velocity_values)
 *         , old_strain_rates(scratch.old_strain_rates)
 *         , old_old_strain_rates(scratch.old_old_strain_rates)
 *         ,
 * 
 *         old_temperature_values(scratch.old_temperature_values)
 *         , old_old_temperature_values(scratch.old_old_temperature_values)
 *         , old_temperature_grads(scratch.old_temperature_grads)
 *         , old_old_temperature_grads(scratch.old_old_temperature_grads)
 *         , old_temperature_laplacians(scratch.old_temperature_laplacians)
 *         , old_old_temperature_laplacians(scratch.old_old_temperature_laplacians)
 *       {}
 *     } // namespace Scratch
 * 
 * 
 * @endcode
 * 
 * The CopyData objects are even simpler than the Scratch objects as all
 * they have to do is to store the results of local computations until
 * they can be copied into the global matrix or vector objects. These
 * structures therefore only need to provide a constructor, a copy
 * operation, and some arrays for local matrix, local vectors and the
 * relation between local and global degrees of freedom (a.k.a.
 * <code>local_dof_indices</code>). Again, we have one such structure for
 * each of the four operations we will parallelize using the WorkStream
 * class:
 * 
 * @code
 *     namespace CopyData
 *     {
 *       template <int dim>
 *       struct StokesPreconditioner
 *       {
 *         StokesPreconditioner(const FiniteElement<dim> &stokes_fe);
 *         StokesPreconditioner(const StokesPreconditioner &data);
 *         StokesPreconditioner &operator=(const StokesPreconditioner &) = default;
 * 
 *         FullMatrix<double>                   local_matrix;
 *         std::vector<types::global_dof_index> local_dof_indices;
 *       };
 * 
 *       template <int dim>
 *       StokesPreconditioner<dim>::StokesPreconditioner(
 *         const FiniteElement<dim> &stokes_fe)
 *         : local_matrix(stokes_fe.n_dofs_per_cell(), stokes_fe.n_dofs_per_cell())
 *         , local_dof_indices(stokes_fe.n_dofs_per_cell())
 *       {}
 * 
 *       template <int dim>
 *       StokesPreconditioner<dim>::StokesPreconditioner(
 *         const StokesPreconditioner &data)
 *         : local_matrix(data.local_matrix)
 *         , local_dof_indices(data.local_dof_indices)
 *       {}
 * 
 * 
 * 
 *       template <int dim>
 *       struct StokesSystem : public StokesPreconditioner<dim>
 *       {
 *         StokesSystem(const FiniteElement<dim> &stokes_fe);
 * 
 *         Vector<double> local_rhs;
 *       };
 * 
 *       template <int dim>
 *       StokesSystem<dim>::StokesSystem(const FiniteElement<dim> &stokes_fe)
 *         : StokesPreconditioner<dim>(stokes_fe)
 *         , local_rhs(stokes_fe.n_dofs_per_cell())
 *       {}
 * 
 * 
 * 
 *       template <int dim>
 *       struct TemperatureMatrix
 *       {
 *         TemperatureMatrix(const FiniteElement<dim> &temperature_fe);
 * 
 *         FullMatrix<double>                   local_mass_matrix;
 *         FullMatrix<double>                   local_stiffness_matrix;
 *         std::vector<types::global_dof_index> local_dof_indices;
 *       };
 * 
 *       template <int dim>
 *       TemperatureMatrix<dim>::TemperatureMatrix(
 *         const FiniteElement<dim> &temperature_fe)
 *         : local_mass_matrix(temperature_fe.n_dofs_per_cell(),
 *                             temperature_fe.n_dofs_per_cell())
 *         , local_stiffness_matrix(temperature_fe.n_dofs_per_cell(),
 *                                  temperature_fe.n_dofs_per_cell())
 *         , local_dof_indices(temperature_fe.n_dofs_per_cell())
 *       {}
 * 
 * 
 * 
 *       template <int dim>
 *       struct TemperatureRHS
 *       {
 *         TemperatureRHS(const FiniteElement<dim> &temperature_fe);
 * 
 *         Vector<double>                       local_rhs;
 *         std::vector<types::global_dof_index> local_dof_indices;
 *         FullMatrix<double>                   matrix_for_bc;
 *       };
 * 
 *       template <int dim>
 *       TemperatureRHS<dim>::TemperatureRHS(
 *         const FiniteElement<dim> &temperature_fe)
 *         : local_rhs(temperature_fe.n_dofs_per_cell())
 *         , local_dof_indices(temperature_fe.n_dofs_per_cell())
 *         , matrix_for_bc(temperature_fe.n_dofs_per_cell(),
 *                         temperature_fe.n_dofs_per_cell())
 *       {}
 *     } // namespace CopyData
 *   }   // namespace Assembly
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
 * This is the declaration of the main class. It is very similar to step-31
 * but there are a number differences we will comment on below.
 *   

 * 
 * The top of the class is essentially the same as in step-31, listing the
 * public methods and a set of private functions that do the heavy
 * lifting. Compared to step-31 there are only two additions to this
 * section: the function <code>get_cfl_number()</code> that computes the
 * maximum CFL number over all cells which we then compute the global time
 * step from, and the function <code>get_entropy_variation()</code> that is
 * used in the computation of the entropy stabilization. It is akin to the
 * <code>get_extrapolated_temperature_range()</code> we have used in step-31
 * for this purpose, but works on the entropy instead of the temperature
 * instead.
 * 
 * @code
 *   template <int dim>
 *   class BoussinesqFlowProblem
 *   {
 *   public:
 *     struct Parameters;
 *     BoussinesqFlowProblem(Parameters &parameters);
 *     void run();
 * 
 *   private:
 *     void   setup_dofs();
 *     void   assemble_stokes_preconditioner();
 *     void   build_stokes_preconditioner();
 *     void   assemble_stokes_system();
 *     void   assemble_temperature_matrix();
 *     void   assemble_temperature_system(const double maximal_velocity);
 *     double get_maximal_velocity() const;
 *     double get_cfl_number() const;
 *     double get_entropy_variation(const double average_temperature) const;
 *     std::pair<double, double> get_extrapolated_temperature_range() const;
 *     void                      solve();
 *     void                      output_results();
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
 *       const std::vector<SymmetricTensor<2, dim>> &old_strain_rates,
 *       const std::vector<SymmetricTensor<2, dim>> &old_old_strain_rates,
 *       const double                                global_u_infty,
 *       const double                                global_T_variation,
 *       const double                                average_temperature,
 *       const double                                global_entropy_variation,
 *       const double                                cell_diameter) const;
 * 
 *   public:
 * @endcode
 * 
 * The first significant new component is the definition of a struct for
 * the parameters according to the discussion in the introduction. This
 * structure is initialized by reading from a parameter file during
 * construction of this object.
 * 
 * @code
 *     struct Parameters
 *     {
 *       Parameters(const std::string &parameter_filename);
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 *       void        parse_parameters(ParameterHandler &prm);
 * 
 *       double end_time;
 * 
 *       unsigned int initial_global_refinement;
 *       unsigned int initial_adaptive_refinement;
 * 
 *       bool         generate_graphical_output;
 *       unsigned int graphical_output_interval;
 * 
 *       unsigned int adaptive_refinement_interval;
 * 
 *       double stabilization_alpha;
 *       double stabilization_c_R;
 *       double stabilization_beta;
 * 
 *       unsigned int stokes_velocity_degree;
 *       bool         use_locally_conservative_discretization;
 * 
 *       unsigned int temperature_degree;
 *     };
 * 
 *   private:
 *     Parameters &parameters;
 * 
 * @endcode
 * 
 * The <code>pcout</code> (for <i>%parallel <code>std::cout</code></i>)
 * object is used to simplify writing output: each MPI process can use
 * this to generate output as usual, but since each of these processes
 * will (hopefully) produce the same output it will just be replicated
 * many times over; with the ConditionalOStream class, only the output
 * generated by one MPI process will actually be printed to screen,
 * whereas the output by all the other threads will simply be forgotten.
 * 
 * @code
 *     ConditionalOStream pcout;
 * 
 * @endcode
 * 
 * The following member variables will then again be similar to those in
 * step-31 (and to other tutorial programs). As mentioned in the
 * introduction, we fully distribute computations, so we will have to use
 * the parallel::distributed::Triangulation class (see step-40) but the
 * remainder of these variables is rather standard with two exceptions:
 *     

 * 
 * - The <code>mapping</code> variable is used to denote a higher-order
 * polynomial mapping. As mentioned in the introduction, we use this
 * mapping when forming integrals through quadrature for all cells that
 * are adjacent to either the inner or outer boundaries of our domain
 * where the boundary is curved.
 *     

 * 
 * - In a bit of naming confusion, you will notice below that some of the
 * variables from namespace TrilinosWrappers are taken from namespace
 * TrilinosWrappers::MPI (such as the right hand side vectors) whereas
 * others are not (such as the various matrices). This is due to legacy
 * reasons. We will frequently have to query velocities
 * and temperatures at arbitrary quadrature points; consequently, rather
 * than importing ghost information of a vector whenever we need access
 * to degrees of freedom that are relevant locally but owned by another
 * processor, we solve linear systems in %parallel but then immediately
 * initialize a vector including ghost entries of the solution for further
 * processing. The various <code>*_solution</code> vectors are therefore
 * filled immediately after solving their respective linear system in
 * %parallel and will always contain values for all
 * @ref GlossLocallyRelevantDof "locally relevant degrees of freedom";
 * the fully distributed vectors that we obtain from the solution process
 * and that only ever contain the
 * @ref GlossLocallyOwnedDof "locally owned degrees of freedom" are
 * destroyed immediately after the solution process and after we have
 * copied the relevant values into the member variable vectors.
 * 
 * @code
 *     parallel::distributed::Triangulation<dim> triangulation;
 *     double                                    global_Omega_diameter;
 * 
 *     const MappingQ<dim> mapping;
 * 
 *     const FESystem<dim>       stokes_fe;
 *     DoFHandler<dim>           stokes_dof_handler;
 *     AffineConstraints<double> stokes_constraints;
 * 
 *     TrilinosWrappers::BlockSparseMatrix stokes_matrix;
 *     TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;
 * 
 *     TrilinosWrappers::MPI::BlockVector stokes_solution;
 *     TrilinosWrappers::MPI::BlockVector old_stokes_solution;
 *     TrilinosWrappers::MPI::BlockVector stokes_rhs;
 * 
 * 
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
 *     std::shared_ptr<TrilinosWrappers::PreconditionAMG>    Amg_preconditioner;
 *     std::shared_ptr<TrilinosWrappers::PreconditionJacobi> Mp_preconditioner;
 *     std::shared_ptr<TrilinosWrappers::PreconditionJacobi> T_preconditioner;
 * 
 *     bool rebuild_stokes_matrix;
 *     bool rebuild_stokes_preconditioner;
 *     bool rebuild_temperature_matrices;
 *     bool rebuild_temperature_preconditioner;
 * 
 * @endcode
 * 
 * The next member variable, <code>computing_timer</code> is used to
 * conveniently account for compute time spent in certain "sections" of
 * the code that are repeatedly entered. For example, we will enter (and
 * leave) sections for Stokes matrix assembly and would like to accumulate
 * the run time spent in this section over all time steps. Every so many
 * time steps as well as at the end of the program (through the destructor
 * of the TimerOutput class) we will then produce a nice summary of the
 * times spent in the different sections into which we categorize the
 * run-time of this program.
 * 
 * @code
 *     TimerOutput computing_timer;
 * 
 * @endcode
 * 
 * After these member variables we have a number of auxiliary functions
 * that have been broken out of the ones listed above. Specifically, there
 * are first three functions that we call from <code>setup_dofs</code> and
 * then the ones that do the assembling of linear systems:
 * 
 * @code
 *     void setup_stokes_matrix(
 *       const std::vector<IndexSet> &stokes_partitioning,
 *       const std::vector<IndexSet> &stokes_relevant_partitioning);
 *     void setup_stokes_preconditioner(
 *       const std::vector<IndexSet> &stokes_partitioning,
 *       const std::vector<IndexSet> &stokes_relevant_partitioning);
 *     void setup_temperature_matrices(
 *       const IndexSet &temperature_partitioning,
 *       const IndexSet &temperature_relevant_partitioning);
 * 
 * 
 * @endcode
 * 
 * Following the @ref MTWorkStream "task-based parallelization" paradigm,
 * we split all the assembly routines into two parts: a first part that
 * can do all the calculations on a certain cell without taking care of
 * other threads, and a second part (which is writing the local data into
 * the global matrices and vectors) which can be entered by only one
 * thread at a time. In order to implement that, we provide functions for
 * each of those two steps for all the four assembly routines that we use
 * in this program. The following eight functions do exactly this:
 * 
 * @code
 *     void local_assemble_stokes_preconditioner(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       Assembly::Scratch::StokesPreconditioner<dim> &        scratch,
 *       Assembly::CopyData::StokesPreconditioner<dim> &       data);
 * 
 *     void copy_local_to_global_stokes_preconditioner(
 *       const Assembly::CopyData::StokesPreconditioner<dim> &data);
 * 
 * 
 *     void local_assemble_stokes_system(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       Assembly::Scratch::StokesSystem<dim> &                scratch,
 *       Assembly::CopyData::StokesSystem<dim> &               data);
 * 
 *     void copy_local_to_global_stokes_system(
 *       const Assembly::CopyData::StokesSystem<dim> &data);
 * 
 * 
 *     void local_assemble_temperature_matrix(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       Assembly::Scratch::TemperatureMatrix<dim> &           scratch,
 *       Assembly::CopyData::TemperatureMatrix<dim> &          data);
 * 
 *     void copy_local_to_global_temperature_matrix(
 *       const Assembly::CopyData::TemperatureMatrix<dim> &data);
 * 
 * 
 * 
 *     void local_assemble_temperature_rhs(
 *       const std::pair<double, double> global_T_range,
 *       const double                    global_max_velocity,
 *       const double                    global_entropy_variation,
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       Assembly::Scratch::TemperatureRHS<dim> &              scratch,
 *       Assembly::CopyData::TemperatureRHS<dim> &             data);
 * 
 *     void copy_local_to_global_temperature_rhs(
 *       const Assembly::CopyData::TemperatureRHS<dim> &data);
 * 
 * @endcode
 * 
 * Finally, we forward declare a member class that we will define later on
 * and that will be used to compute a number of quantities from our
 * solution vectors that we'd like to put into the output files for
 * visualization.
 * 
 * @code
 *     class Postprocessor;
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
 * <a name="BoussinesqFlowProblemParameters"></a> 
 * <h4>BoussinesqFlowProblem::Parameters</h4>
 *   

 * 
 * Here comes the definition of the parameters for the Stokes problem. We
 * allow to set the end time for the simulation, the level of refinements
 * (both global and adaptive, which in the sum specify what maximum level
 * the cells are allowed to have), and the interval between refinements in
 * the time stepping.
 *   

 * 
 * Then, we let the user specify constants for the stabilization parameters
 * (as discussed in the introduction), the polynomial degree for the Stokes
 * velocity space, whether to use the locally conservative discretization
 * based on FE_DGP elements for the pressure or not (FE_Q elements for
 * pressure), and the polynomial degree for the temperature interpolation.
 *   

 * 
 * The constructor checks for a valid input file (if not, a file with
 * default parameters for the quantities is written), and eventually parses
 * the parameters.
 * 
 * @code
 *   template <int dim>
 *   BoussinesqFlowProblem<dim>::Parameters::Parameters(
 *     const std::string &parameter_filename)
 *     : end_time(1e8)
 *     , initial_global_refinement(2)
 *     , initial_adaptive_refinement(2)
 *     , adaptive_refinement_interval(10)
 *     , stabilization_alpha(2)
 *     , stabilization_c_R(0.11)
 *     , stabilization_beta(0.078)
 *     , stokes_velocity_degree(2)
 *     , use_locally_conservative_discretization(true)
 *     , temperature_degree(2)
 *   {
 *     ParameterHandler prm;
 *     BoussinesqFlowProblem<dim>::Parameters::declare_parameters(prm);
 * 
 *     std::ifstream parameter_file(parameter_filename);
 * 
 *     if (!parameter_file)
 *       {
 *         parameter_file.close();
 * 
 *         std::ofstream parameter_out(parameter_filename);
 *         prm.print_parameters(parameter_out, ParameterHandler::Text);
 * 
 *         AssertThrow(
 *           false,
 *           ExcMessage(
 *             "Input parameter file <" + parameter_filename +
 *             "> not found. Creating a template file of the same name."));
 *       }
 * 
 *     prm.parse_input(parameter_file);
 *     parse_parameters(prm);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Next we have a function that declares the parameters that we expect in
 * the input file, together with their data types, default values and a
 * description:
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::Parameters::declare_parameters(
 *     ParameterHandler &prm)
 *   {
 *     prm.declare_entry("End time",
 *                       "1e8",
 *                       Patterns::Double(0),
 *                       "The end time of the simulation in years.");
 *     prm.declare_entry("Initial global refinement",
 *                       "2",
 *                       Patterns::Integer(0),
 *                       "The number of global refinement steps performed on "
 *                       "the initial coarse mesh, before the problem is first "
 *                       "solved there.");
 *     prm.declare_entry("Initial adaptive refinement",
 *                       "2",
 *                       Patterns::Integer(0),
 *                       "The number of adaptive refinement steps performed after "
 *                       "initial global refinement.");
 *     prm.declare_entry("Time steps between mesh refinement",
 *                       "10",
 *                       Patterns::Integer(1),
 *                       "The number of time steps after which the mesh is to be "
 *                       "adapted based on computed error indicators.");
 *     prm.declare_entry("Generate graphical output",
 *                       "false",
 *                       Patterns::Bool(),
 *                       "Whether graphical output is to be generated or not. "
 *                       "You may not want to get graphical output if the number "
 *                       "of processors is large.");
 *     prm.declare_entry("Time steps between graphical output",
 *                       "50",
 *                       Patterns::Integer(1),
 *                       "The number of time steps between each generation of "
 *                       "graphical output files.");
 * 
 *     prm.enter_subsection("Stabilization parameters");
 *     {
 *       prm.declare_entry("alpha",
 *                         "2",
 *                         Patterns::Double(1, 2),
 *                         "The exponent in the entropy viscosity stabilization.");
 *       prm.declare_entry("c_R",
 *                         "0.11",
 *                         Patterns::Double(0),
 *                         "The c_R factor in the entropy viscosity "
 *                         "stabilization.");
 *       prm.declare_entry("beta",
 *                         "0.078",
 *                         Patterns::Double(0),
 *                         "The beta factor in the artificial viscosity "
 *                         "stabilization. An appropriate value for 2d is 0.052 "
 *                         "and 0.078 for 3d.");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Discretization");
 *     {
 *       prm.declare_entry(
 *         "Stokes velocity polynomial degree",
 *         "2",
 *         Patterns::Integer(1),
 *         "The polynomial degree to use for the velocity variables "
 *         "in the Stokes system.");
 *       prm.declare_entry(
 *         "Temperature polynomial degree",
 *         "2",
 *         Patterns::Integer(1),
 *         "The polynomial degree to use for the temperature variable.");
 *       prm.declare_entry(
 *         "Use locally conservative discretization",
 *         "true",
 *         Patterns::Bool(),
 *         "Whether to use a Stokes discretization that is locally "
 *         "conservative at the expense of a larger number of degrees "
 *         "of freedom, or to go with a cheaper discretization "
 *         "that does not locally conserve mass (although it is "
 *         "globally conservative.");
 *     }
 *     prm.leave_subsection();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * And then we need a function that reads the contents of the
 * ParameterHandler object we get by reading the input file and puts the
 * results into variables that store the values of the parameters we have
 * previously declared:
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::Parameters::parse_parameters(
 *     ParameterHandler &prm)
 *   {
 *     end_time                  = prm.get_double("End time");
 *     initial_global_refinement = prm.get_integer("Initial global refinement");
 *     initial_adaptive_refinement =
 *       prm.get_integer("Initial adaptive refinement");
 * 
 *     adaptive_refinement_interval =
 *       prm.get_integer("Time steps between mesh refinement");
 * 
 *     generate_graphical_output = prm.get_bool("Generate graphical output");
 *     graphical_output_interval =
 *       prm.get_integer("Time steps between graphical output");
 * 
 *     prm.enter_subsection("Stabilization parameters");
 *     {
 *       stabilization_alpha = prm.get_double("alpha");
 *       stabilization_c_R   = prm.get_double("c_R");
 *       stabilization_beta  = prm.get_double("beta");
 *     }
 *     prm.leave_subsection();
 * 
 *     prm.enter_subsection("Discretization");
 *     {
 *       stokes_velocity_degree =
 *         prm.get_integer("Stokes velocity polynomial degree");
 *       temperature_degree = prm.get_integer("Temperature polynomial degree");
 *       use_locally_conservative_discretization =
 *         prm.get_bool("Use locally conservative discretization");
 *     }
 *     prm.leave_subsection();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemBoussinesqFlowProblem"></a> 
 * <h4>BoussinesqFlowProblem::BoussinesqFlowProblem</h4>
 *   

 * 
 * The constructor of the problem is very similar to the constructor in
 * step-31. What is different is the %parallel communication: Trilinos uses
 * a message passing interface (MPI) for data distribution. When entering
 * the BoussinesqFlowProblem class, we have to decide how the parallelization
 * is to be done. We choose a rather simple strategy and let all processors
 * that are running the program work together, specified by the communicator
 * <code>MPI_COMM_WORLD</code>. Next, we create the output stream (as we
 * already did in step-18) that only generates output on the first MPI
 * process and is completely forgetful on all others. The implementation of
 * this idea is to check the process number when <code>pcout</code> gets a
 * true argument, and it uses the <code>std::cout</code> stream for
 * output. If we are one processor five, for instance, then we will give a
 * <code>false</code> argument to <code>pcout</code>, which means that the
 * output of that processor will not be printed. With the exception of the
 * mapping object (for which we use polynomials of degree 4) all but the
 * final member variable are exactly the same as in step-31.
 *   

 * 
 * This final object, the TimerOutput object, is then told to restrict
 * output to the <code>pcout</code> stream (processor 0), and then we
 * specify that we want to get a summary table at the end of the program
 * which shows us wallclock times (as opposed to CPU times). We will
 * manually also request intermediate summaries every so many time steps in
 * the <code>run()</code> function below.
 * 
 * @code
 *   template <int dim>
 *   BoussinesqFlowProblem<dim>::BoussinesqFlowProblem(Parameters &parameters_)
 *     : parameters(parameters_)
 *     , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0))
 *     ,
 * 
 *     triangulation(MPI_COMM_WORLD,
 *                   typename Triangulation<dim>::MeshSmoothing(
 *                     Triangulation<dim>::smoothing_on_refinement |
 *                     Triangulation<dim>::smoothing_on_coarsening))
 *     ,
 * 
 *     global_Omega_diameter(0.)
 *     ,
 * 
 *     mapping(4)
 *     ,
 * 
 *     stokes_fe(FE_Q<dim>(parameters.stokes_velocity_degree),
 *               dim,
 *               (parameters.use_locally_conservative_discretization ?
 *                  static_cast<const FiniteElement<dim> &>(
 *                    FE_DGP<dim>(parameters.stokes_velocity_degree - 1)) :
 *                  static_cast<const FiniteElement<dim> &>(
 *                    FE_Q<dim>(parameters.stokes_velocity_degree - 1))),
 *               1)
 *     ,
 * 
 *     stokes_dof_handler(triangulation)
 *     ,
 * 
 *     temperature_fe(parameters.temperature_degree)
 *     , temperature_dof_handler(triangulation)
 *     ,
 * 
 *     time_step(0)
 *     , old_time_step(0)
 *     , timestep_number(0)
 *     , rebuild_stokes_matrix(true)
 *     , rebuild_stokes_preconditioner(true)
 *     , rebuild_temperature_matrices(true)
 *     , rebuild_temperature_preconditioner(true)
 *     ,
 * 
 *     computing_timer(MPI_COMM_WORLD,
 *                     pcout,
 *                     TimerOutput::summary,
 *                     TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheBoussinesqFlowProblemhelperfunctions"></a> 
 * <h4>The BoussinesqFlowProblem helper functions</h4>
 * 
 * <a name="BoussinesqFlowProblemget_maximal_velocity"></a> 
 * <h5>BoussinesqFlowProblem::get_maximal_velocity</h5>
 * 

 * 
 * Except for two small details, the function to compute the global maximum
 * of the velocity is the same as in step-31. The first detail is actually
 * common to all functions that implement loops over all cells in the
 * triangulation: When operating in %parallel, each processor can only work
 * on a chunk of cells since each processor only has a certain part of the
 * entire triangulation. This chunk of cells that we want to work on is
 * identified via a so-called <code>subdomain_id</code>, as we also did in
 * step-18. All we need to change is hence to perform the cell-related
 * operations only on cells that are owned by the current process (as
 * opposed to ghost or artificial cells), i.e. for which the subdomain id
 * equals the number of the process ID. Since this is a commonly used
 * operation, there is a shortcut for this operation: we can ask whether the
 * cell is owned by the current processor using
 * <code>cell-@>is_locally_owned()</code>.
 *   

 * 
 * The second difference is the way we calculate the maximum value. Before,
 * we could simply have a <code>double</code> variable that we checked
 * against on each quadrature point for each cell. Now, we have to be a bit
 * more careful since each processor only operates on a subset of
 * cells. What we do is to first let each processor calculate the maximum
 * among its cells, and then do a global communication operation
 * <code>Utilities::MPI::max</code> that computes the maximum value among
 * all the maximum values of the individual processors. MPI provides such a
 * call, but it's even simpler to use the respective function in namespace
 * Utilities::MPI using the MPI communicator object since that will do the
 * right thing even if we work without MPI and on a single machine only. The
 * call to <code>Utilities::MPI::max</code> needs two arguments, namely the
 * local maximum (input) and the MPI communicator, which is MPI_COMM_WORLD
 * in this example.
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::get_maximal_velocity() const
 *   {
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(),
 *                                             parameters.stokes_velocity_degree);
 *     const unsigned int   n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim>               fe_values(mapping,
 *                             stokes_fe,
 *                             quadrature_formula,
 *                             update_values);
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     double max_local_velocity = 0;
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit(cell);
 *           fe_values[velocities].get_function_values(stokes_solution,
 *                                                     velocity_values);
 * 
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             max_local_velocity =
 *               std::max(max_local_velocity, velocity_values[q].norm());
 *         }
 * 
 *     return Utilities::MPI::max(max_local_velocity, MPI_COMM_WORLD);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_cfl_number"></a> 
 * <h5>BoussinesqFlowProblem::get_cfl_number</h5>
 * 

 * 
 * The next function does something similar, but we now compute the CFL
 * number, i.e., maximal velocity on a cell divided by the cell
 * diameter. This number is necessary to determine the time step size, as we
 * use a semi-explicit time stepping scheme for the temperature equation
 * (see step-31 for a discussion). We compute it in the same way as above:
 * Compute the local maximum over all locally owned cells, then exchange it
 * via MPI to find the global maximum.
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::get_cfl_number() const
 *   {
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(),
 *                                             parameters.stokes_velocity_degree);
 *     const unsigned int   n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim>               fe_values(mapping,
 *                             stokes_fe,
 *                             quadrature_formula,
 *                             update_values);
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     double max_local_cfl = 0;
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit(cell);
 *           fe_values[velocities].get_function_values(stokes_solution,
 *                                                     velocity_values);
 * 
 *           double max_local_velocity = 1e-10;
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             max_local_velocity =
 *               std::max(max_local_velocity, velocity_values[q].norm());
 *           max_local_cfl =
 *             std::max(max_local_cfl, max_local_velocity / cell->diameter());
 *         }
 * 
 *     return Utilities::MPI::max(max_local_cfl, MPI_COMM_WORLD);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_entropy_variation"></a> 
 * <h5>BoussinesqFlowProblem::get_entropy_variation</h5>
 * 

 * 
 * Next comes the computation of the global entropy variation
 * $\|E(T)-\bar{E}(T)\|_\infty$ where the entropy $E$ is defined as
 * discussed in the introduction.  This is needed for the evaluation of the
 * stabilization in the temperature equation as explained in the
 * introduction. The entropy variation is actually only needed if we use
 * $\alpha=2$ as a power in the residual computation. The infinity norm is
 * computed by the maxima over quadrature points, as usual in discrete
 * computations.
 *   

 * 
 * In order to compute this quantity, we first have to find the
 * space-average $\bar{E}(T)$ and then evaluate the maximum. However, that
 * means that we would need to perform two loops. We can avoid the overhead
 * by noting that $\|E(T)-\bar{E}(T)\|_\infty =
 * \max\big(E_{\textrm{max}}(T)-\bar{E}(T),
 * \bar{E}(T)-E_{\textrm{min}}(T)\big)$, i.e., the maximum out of the
 * deviation from the average entropy in positive and negative
 * directions. The four quantities we need for the latter formula (maximum
 * entropy, minimum entropy, average entropy, area) can all be evaluated in
 * the same loop over all cells, so we choose this simpler variant.
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::get_entropy_variation(
 *     const double average_temperature) const
 *   {
 *     if (parameters.stabilization_alpha != 2)
 *       return 1.;
 * 
 *     const QGauss<dim>  quadrature_formula(parameters.temperature_degree + 1);
 *     const unsigned int n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim>       fe_values(temperature_fe,
 *                             quadrature_formula,
 *                             update_values | update_JxW_values);
 *     std::vector<double> old_temperature_values(n_q_points);
 *     std::vector<double> old_old_temperature_values(n_q_points);
 * 
 * @endcode
 * 
 * In the two functions above we computed the maximum of numbers that were
 * all non-negative, so we knew that zero was certainly a lower bound. On
 * the other hand, here we need to find the maximum deviation from the
 * average value, i.e., we will need to know the maximal and minimal
 * values of the entropy for which we don't a priori know the sign.
 *     

 * 
 * To compute it, we can therefore start with the largest and smallest
 * possible values we can store in a double precision number: The minimum
 * is initialized with a bigger and the maximum with a smaller number than
 * any one that is going to appear. We are then guaranteed that these
 * numbers will be overwritten in the loop on the first cell or, if this
 * processor does not own any cells, in the communication step at the
 * latest. The following loop then computes the minimum and maximum local
 * entropy as well as keeps track of the area/volume of the part of the
 * domain we locally own and the integral over the entropy on it:
 * 
 * @code
 *     double min_entropy = std::numeric_limits<double>::max(),
 *            max_entropy = -std::numeric_limits<double>::max(), area = 0,
 *            entropy_integrated = 0;
 * 
 *     for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit(cell);
 *           fe_values.get_function_values(old_temperature_solution,
 *                                         old_temperature_values);
 *           fe_values.get_function_values(old_old_temperature_solution,
 *                                         old_old_temperature_values);
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             {
 *               const double T =
 *                 (old_temperature_values[q] + old_old_temperature_values[q]) / 2;
 *               const double entropy =
 *                 ((T - average_temperature) * (T - average_temperature));
 * 
 *               min_entropy = std::min(min_entropy, entropy);
 *               max_entropy = std::max(max_entropy, entropy);
 *               area += fe_values.JxW(q);
 *               entropy_integrated += fe_values.JxW(q) * entropy;
 *             }
 *         }
 * 
 * @endcode
 * 
 * Now we only need to exchange data between processors: we need to sum
 * the two integrals (<code>area</code>, <code>entropy_integrated</code>),
 * and get the extrema for maximum and minimum. We could do this through
 * four different data exchanges, but we can it with two:
 * Utilities::MPI::sum also exists in a variant that takes an array of
 * values that are all to be summed up. And we can also utilize the
 * Utilities::MPI::max function by realizing that forming the minimum over
 * the minimal entropies equals forming the negative of the maximum over
 * the negative of the minimal entropies; this maximum can then be
 * combined with forming the maximum over the maximal entropies.
 * 
 * @code
 *     const double local_sums[2]   = {entropy_integrated, area},
 *                  local_maxima[2] = {-min_entropy, max_entropy};
 *     double global_sums[2], global_maxima[2];
 * 
 *     Utilities::MPI::sum(local_sums, MPI_COMM_WORLD, global_sums);
 *     Utilities::MPI::max(local_maxima, MPI_COMM_WORLD, global_maxima);
 * 
 * @endcode
 * 
 * Having computed everything this way, we can then compute the average
 * entropy and find the $L^\infty$ norm by taking the larger of the
 * deviation of the maximum or minimum from the average:
 * 
 * @code
 *     const double average_entropy = global_sums[0] / global_sums[1];
 *     const double entropy_diff    = std::max(global_maxima[1] - average_entropy,
 *                                          average_entropy - (-global_maxima[0]));
 *     return entropy_diff;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_extrapolated_temperature_range"></a> 
 * <h5>BoussinesqFlowProblem::get_extrapolated_temperature_range</h5>
 * 

 * 
 * The next function computes the minimal and maximal value of the
 * extrapolated temperature over the entire domain. Again, this is only a
 * slightly modified version of the respective function in step-31. As in
 * the function above, we collect local minima and maxima and then compute
 * the global extrema using the same trick as above.
 *   

 * 
 * As already discussed in step-31, the function needs to distinguish
 * between the first and all following time steps because it uses a higher
 * order temperature extrapolation scheme when at least two previous time
 * steps are available.
 * 
 * @code
 *   template <int dim>
 *   std::pair<double, double>
 *   BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range() const
 *   {
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(),
 *                                             parameters.temperature_degree);
 *     const unsigned int   n_q_points = quadrature_formula.size();
 * 
 *     FEValues<dim>       fe_values(mapping,
 *                             temperature_fe,
 *                             quadrature_formula,
 *                             update_values);
 *     std::vector<double> old_temperature_values(n_q_points);
 *     std::vector<double> old_old_temperature_values(n_q_points);
 * 
 *     double min_local_temperature = std::numeric_limits<double>::max(),
 *            max_local_temperature = -std::numeric_limits<double>::max();
 * 
 *     if (timestep_number != 0)
 *       {
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *           if (cell->is_locally_owned())
 *             {
 *               fe_values.reinit(cell);
 *               fe_values.get_function_values(old_temperature_solution,
 *                                             old_temperature_values);
 *               fe_values.get_function_values(old_old_temperature_solution,
 *                                             old_old_temperature_values);
 * 
 *               for (unsigned int q = 0; q < n_q_points; ++q)
 *                 {
 *                   const double temperature =
 *                     (1. + time_step / old_time_step) *
 *                       old_temperature_values[q] -
 *                     time_step / old_time_step * old_old_temperature_values[q];
 * 
 *                   min_local_temperature =
 *                     std::min(min_local_temperature, temperature);
 *                   max_local_temperature =
 *                     std::max(max_local_temperature, temperature);
 *                 }
 *             }
 *       }
 *     else
 *       {
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators())
 *           if (cell->is_locally_owned())
 *             {
 *               fe_values.reinit(cell);
 *               fe_values.get_function_values(old_temperature_solution,
 *                                             old_temperature_values);
 * 
 *               for (unsigned int q = 0; q < n_q_points; ++q)
 *                 {
 *                   const double temperature = old_temperature_values[q];
 * 
 *                   min_local_temperature =
 *                     std::min(min_local_temperature, temperature);
 *                   max_local_temperature =
 *                     std::max(max_local_temperature, temperature);
 *                 }
 *             }
 *       }
 * 
 *     double local_extrema[2] = {-min_local_temperature, max_local_temperature};
 *     double global_extrema[2];
 *     Utilities::MPI::max(local_extrema, MPI_COMM_WORLD, global_extrema);
 * 
 *     return std::make_pair(-global_extrema[0], global_extrema[1]);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemcompute_viscosity"></a> 
 * <h5>BoussinesqFlowProblem::compute_viscosity</h5>
 * 

 * 
 * The function that calculates the viscosity is purely local and so needs
 * no communication at all. It is mostly the same as in step-31 but with an
 * updated formulation of the viscosity if $\alpha=2$ is chosen:
 * 
 * @code
 *   template <int dim>
 *   double BoussinesqFlowProblem<dim>::compute_viscosity(
 *     const std::vector<double> &                 old_temperature,
 *     const std::vector<double> &                 old_old_temperature,
 *     const std::vector<Tensor<1, dim>> &         old_temperature_grads,
 *     const std::vector<Tensor<1, dim>> &         old_old_temperature_grads,
 *     const std::vector<double> &                 old_temperature_laplacians,
 *     const std::vector<double> &                 old_old_temperature_laplacians,
 *     const std::vector<Tensor<1, dim>> &         old_velocity_values,
 *     const std::vector<Tensor<1, dim>> &         old_old_velocity_values,
 *     const std::vector<SymmetricTensor<2, dim>> &old_strain_rates,
 *     const std::vector<SymmetricTensor<2, dim>> &old_old_strain_rates,
 *     const double                                global_u_infty,
 *     const double                                global_T_variation,
 *     const double                                average_temperature,
 *     const double                                global_entropy_variation,
 *     const double                                cell_diameter) const
 *   {
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
 *         const SymmetricTensor<2, dim> strain_rate =
 *           (old_strain_rates[q] + old_old_strain_rates[q]) / 2;
 * 
 *         const double T = (old_temperature[q] + old_old_temperature[q]) / 2;
 *         const double dT_dt =
 *           (old_temperature[q] - old_old_temperature[q]) / old_time_step;
 *         const double u_grad_T =
 *           u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2;
 * 
 *         const double kappa_Delta_T =
 *           EquationData::kappa *
 *           (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) /
 *           2;
 *         const double gamma =
 *           ((EquationData::radiogenic_heating * EquationData::density(T) +
 *             2 * EquationData::eta * strain_rate * strain_rate) /
 *            (EquationData::density(T) * EquationData::specific_heat));
 * 
 *         double residual = std::abs(dT_dt + u_grad_T - kappa_Delta_T - gamma);
 *         if (parameters.stabilization_alpha == 2)
 *           residual *= std::abs(T - average_temperature);
 * 
 *         max_residual = std::max(residual, max_residual);
 *         max_velocity = std::max(std::sqrt(u * u), max_velocity);
 *       }
 * 
 *     const double max_viscosity =
 *       (parameters.stabilization_beta * max_velocity * cell_diameter);
 *     if (timestep_number == 0)
 *       return max_viscosity;
 *     else
 *       {
 *         Assert(old_time_step > 0, ExcInternalError());
 * 
 *         double entropy_viscosity;
 *         if (parameters.stabilization_alpha == 2)
 *           entropy_viscosity =
 *             (parameters.stabilization_c_R * cell_diameter * cell_diameter *
 *              max_residual / global_entropy_variation);
 *         else
 *           entropy_viscosity =
 *             (parameters.stabilization_c_R * cell_diameter *
 *              global_Omega_diameter * max_velocity * max_residual /
 *              (global_u_infty * global_T_variation));
 * 
 *         return std::min(max_viscosity, entropy_viscosity);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheBoussinesqFlowProblemsetupfunctions"></a> 
 * <h4>The BoussinesqFlowProblem setup functions</h4>
 * 

 * 
 * The following three functions set up the Stokes matrix, the matrix used
 * for the Stokes preconditioner, and the temperature matrix. The code is
 * mostly the same as in step-31, but it has been broken out into three
 * functions of their own for simplicity.
 *   

 * 
 * The main functional difference between the code here and that in step-31
 * is that the matrices we want to set up are distributed across multiple
 * processors. Since we still want to build up the sparsity pattern first
 * for efficiency reasons, we could continue to build the <i>entire</i>
 * sparsity pattern as a BlockDynamicSparsityPattern, as we did in
 * step-31. However, that would be inefficient: every processor would build
 * the same sparsity pattern, but only initialize a small part of the matrix
 * using it. It also violates the principle that every processor should only
 * work on those cells it owns (and, if necessary the layer of ghost cells
 * around it).
 *   

 * 
 * Rather, we use an object of type TrilinosWrappers::BlockSparsityPattern,
 * which is (obviously) a wrapper around a sparsity pattern object provided
 * by Trilinos. The advantage is that the Trilinos sparsity pattern class
 * can communicate across multiple processors: if this processor fills in
 * all the nonzero entries that result from the cells it owns, and every
 * other processor does so as well, then at the end after some MPI
 * communication initiated by the <code>compress()</code> call, we will have
 * the globally assembled sparsity pattern available with which the global
 * matrix can be initialized.
 *   

 * 
 * There is one important aspect when initializing Trilinos sparsity
 * patterns in parallel: In addition to specifying the locally owned rows
 * and columns of the matrices via the @p stokes_partitioning index set, we
 * also supply information about all the rows we are possibly going to write
 * into when assembling on a certain processor. The set of locally relevant
 * rows contains all such rows (possibly also a few unnecessary ones, but it
 * is difficult to find the exact row indices before actually getting
 * indices on all cells and resolving constraints). This additional
 * information allows to exactly determine the structure for the
 * off-processor data found during assembly. While Trilinos matrices are
 * able to collect this information on the fly as well (when initializing
 * them from some other reinit method), it is less efficient and leads to
 * problems when assembling matrices with multiple threads. In this program,
 * we pessimistically assume that only one processor at a time can write
 * into the matrix while assembly (whereas the computation is parallel),
 * which is fine for Trilinos matrices. In practice, one can do better by
 * hinting WorkStream at cells that do not share vertices, allowing for
 * parallelism among those cells (see the graph coloring algorithms and
 * WorkStream with colored iterators argument). However, that only works
 * when only one MPI processor is present because Trilinos' internal data
 * structures for accumulating off-processor data on the fly are not thread
 * safe. With the initialization presented here, there is no such problem
 * and one could safely introduce graph coloring for this algorithm.
 *   

 * 
 * The only other change we need to make is to tell the
 * DoFTools::make_sparsity_pattern() function that it is only supposed to
 * work on a subset of cells, namely the ones whose
 * <code>subdomain_id</code> equals the number of the current processor, and
 * to ignore all other cells.
 *   

 * 
 * This strategy is replicated across all three of the following functions.
 *   

 * 
 * Note that Trilinos matrices store the information contained in the
 * sparsity patterns, so we can safely release the <code>sp</code> variable
 * once the matrix has been given the sparsity structure.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::setup_stokes_matrix(
 *     const std::vector<IndexSet> &stokes_partitioning,
 *     const std::vector<IndexSet> &stokes_relevant_partitioning)
 *   {
 *     stokes_matrix.clear();
 * 
 *     TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning,
 *                                               stokes_partitioning,
 *                                               stokes_relevant_partitioning,
 *                                               MPI_COMM_WORLD);
 * 
 *     Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *     for (unsigned int c = 0; c < dim + 1; ++c)
 *       for (unsigned int d = 0; d < dim + 1; ++d)
 *         if (!((c == dim) && (d == dim)))
 *           coupling[c][d] = DoFTools::always;
 *         else
 *           coupling[c][d] = DoFTools::none;
 * 
 *     DoFTools::make_sparsity_pattern(stokes_dof_handler,
 *                                     coupling,
 *                                     sp,
 *                                     stokes_constraints,
 *                                     false,
 *                                     Utilities::MPI::this_mpi_process(
 *                                       MPI_COMM_WORLD));
 *     sp.compress();
 * 
 *     stokes_matrix.reinit(sp);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::setup_stokes_preconditioner(
 *     const std::vector<IndexSet> &stokes_partitioning,
 *     const std::vector<IndexSet> &stokes_relevant_partitioning)
 *   {
 *     Amg_preconditioner.reset();
 *     Mp_preconditioner.reset();
 * 
 *     stokes_preconditioner_matrix.clear();
 * 
 *     TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning,
 *                                               stokes_partitioning,
 *                                               stokes_relevant_partitioning,
 *                                               MPI_COMM_WORLD);
 * 
 *     Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
 *     for (unsigned int c = 0; c < dim + 1; ++c)
 *       for (unsigned int d = 0; d < dim + 1; ++d)
 *         if (c == d)
 *           coupling[c][d] = DoFTools::always;
 *         else
 *           coupling[c][d] = DoFTools::none;
 * 
 *     DoFTools::make_sparsity_pattern(stokes_dof_handler,
 *                                     coupling,
 *                                     sp,
 *                                     stokes_constraints,
 *                                     false,
 *                                     Utilities::MPI::this_mpi_process(
 *                                       MPI_COMM_WORLD));
 *     sp.compress();
 * 
 *     stokes_preconditioner_matrix.reinit(sp);
 *   }
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::setup_temperature_matrices(
 *     const IndexSet &temperature_partitioner,
 *     const IndexSet &temperature_relevant_partitioner)
 *   {
 *     T_preconditioner.reset();
 *     temperature_mass_matrix.clear();
 *     temperature_stiffness_matrix.clear();
 *     temperature_matrix.clear();
 * 
 *     TrilinosWrappers::SparsityPattern sp(temperature_partitioner,
 *                                          temperature_partitioner,
 *                                          temperature_relevant_partitioner,
 *                                          MPI_COMM_WORLD);
 *     DoFTools::make_sparsity_pattern(temperature_dof_handler,
 *                                     sp,
 *                                     temperature_constraints,
 *                                     false,
 *                                     Utilities::MPI::this_mpi_process(
 *                                       MPI_COMM_WORLD));
 *     sp.compress();
 * 
 *     temperature_matrix.reinit(sp);
 *     temperature_mass_matrix.reinit(sp);
 *     temperature_stiffness_matrix.reinit(sp);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The remainder of the setup function (after splitting out the three
 * functions above) mostly has to deal with the things we need to do for
 * parallelization across processors. Because setting all of this up is a
 * significant compute time expense of the program, we put everything we do
 * here into a timer group so that we can get summary information about the
 * fraction of time spent in this part of the program at its end.
 *   

 * 
 * At the top as usual we enumerate degrees of freedom and sort them by
 * component/block, followed by writing their numbers to the screen from
 * processor zero. The DoFHandler::distributed_dofs() function, when applied
 * to a parallel::distributed::Triangulation object, sorts degrees of
 * freedom in such a way that all degrees of freedom associated with
 * subdomain zero come before all those associated with subdomain one,
 * etc. For the Stokes part, this entails, however, that velocities and
 * pressures become intermixed, but this is trivially solved by sorting
 * again by blocks; it is worth noting that this latter operation leaves the
 * relative ordering of all velocities and pressures alone, i.e. within the
 * velocity block we will still have all those associated with subdomain
 * zero before all velocities associated with subdomain one, etc. This is
 * important since we store each of the blocks of this matrix distributed
 * across all processors and want this to be done in such a way that each
 * processor stores that part of the matrix that is roughly equal to the
 * degrees of freedom located on those cells that it will actually work on.
 *   

 * 
 * When printing the numbers of degrees of freedom, note that these numbers
 * are going to be large if we use many processors. Consequently, we let the
 * stream put a comma separator in between every three digits. The state of
 * the stream, using the locale, is saved from before to after this
 * operation. While slightly opaque, the code works because the default
 * locale (which we get using the constructor call
 * <code>std::locale("")</code>) implies printing numbers with a comma
 * separator for every third digit (i.e., thousands, millions, billions).
 *   

 * 
 * In this function as well as many below, we measure how much time
 * we spend here and collect that in a section called "Setup dof
 * systems" across function invocations. This is done using an
 * TimerOutput::Scope object that gets a timer going in the section
 * with above name of the `computing_timer` object upon construction
 * of the local variable; the timer is stopped again when the
 * destructor of the `timing_section` variable is called.  This, of
 * course, happens either at the end of the function, or if we leave
 * the function through a `return` statement or when an exception is
 * thrown somewhere -- in other words, whenever we leave this
 * function in any way. The use of such "scope" objects therefore
 * makes sure that we do not have to manually add code that tells
 * the timer to stop at every location where this function may be
 * left.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::setup_dofs()
 *   {
 *     TimerOutput::Scope timing_section(computing_timer, "Setup dof systems");
 * 
 *     stokes_dof_handler.distribute_dofs(stokes_fe);
 * 
 *     std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
 *     stokes_sub_blocks[dim] = 1;
 *     DoFRenumbering::component_wise(stokes_dof_handler, stokes_sub_blocks);
 * 
 *     temperature_dof_handler.distribute_dofs(temperature_fe);
 * 
 *     const std::vector<types::global_dof_index> stokes_dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(stokes_dof_handler, stokes_sub_blocks);
 * 
 *     const unsigned int n_u = stokes_dofs_per_block[0],
 *                        n_p = stokes_dofs_per_block[1],
 *                        n_T = temperature_dof_handler.n_dofs();
 * 
 *     std::locale s = pcout.get_stream().getloc();
 *     pcout.get_stream().imbue(std::locale(""));
 *     pcout << "Number of active cells: " << triangulation.n_global_active_cells()
 *           << " (on " << triangulation.n_levels() << " levels)" << std::endl
 *           << "Number of degrees of freedom: " << n_u + n_p + n_T << " (" << n_u
 *           << '+' << n_p << '+' << n_T << ')' << std::endl
 *           << std::endl;
 *     pcout.get_stream().imbue(s);
 * 
 * 
 * @endcode
 * 
 * After this, we have to set up the various partitioners (of type
 * <code>IndexSet</code>, see the introduction) that describe which parts
 * of each matrix or vector will be stored where, then call the functions
 * that actually set up the matrices, and at the end also resize the
 * various vectors we keep around in this program.
 * 
 * @code
 *     std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
 *     IndexSet              temperature_partitioning(n_T),
 *       temperature_relevant_partitioning(n_T);
 *     IndexSet stokes_relevant_set;
 *     {
 *       IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
 *       stokes_partitioning.push_back(stokes_index_set.get_view(0, n_u));
 *       stokes_partitioning.push_back(stokes_index_set.get_view(n_u, n_u + n_p));
 * 
 *       DoFTools::extract_locally_relevant_dofs(stokes_dof_handler,
 *                                               stokes_relevant_set);
 *       stokes_relevant_partitioning.push_back(
 *         stokes_relevant_set.get_view(0, n_u));
 *       stokes_relevant_partitioning.push_back(
 *         stokes_relevant_set.get_view(n_u, n_u + n_p));
 * 
 *       temperature_partitioning = temperature_dof_handler.locally_owned_dofs();
 *       DoFTools::extract_locally_relevant_dofs(
 *         temperature_dof_handler, temperature_relevant_partitioning);
 *     }
 * 
 * @endcode
 * 
 * Following this, we can compute constraints for the solution vectors,
 * including hanging node constraints and homogeneous and inhomogeneous
 * boundary values for the Stokes and temperature fields. Note that as for
 * everything else, the constraint objects can not hold <i>all</i>
 * constraints on every processor. Rather, each processor needs to store
 * only those that are actually necessary for correctness given that it
 * only assembles linear systems on cells it owns. As discussed in the
 * @ref distributed_paper "this paper", the set of constraints we need to
 * know about is exactly the set of constraints on all locally relevant
 * degrees of freedom, so this is what we use to initialize the constraint
 * objects.
 * 
 * @code
 *     {
 *       stokes_constraints.clear();
 *       stokes_constraints.reinit(stokes_relevant_set);
 * 
 *       DoFTools::make_hanging_node_constraints(stokes_dof_handler,
 *                                               stokes_constraints);
 * 
 *       FEValuesExtractors::Vector velocity_components(0);
 *       VectorTools::interpolate_boundary_values(
 *         stokes_dof_handler,
 *         0,
 *         Functions::ZeroFunction<dim>(dim + 1),
 *         stokes_constraints,
 *         stokes_fe.component_mask(velocity_components));
 * 
 *       std::set<types::boundary_id> no_normal_flux_boundaries;
 *       no_normal_flux_boundaries.insert(1);
 *       VectorTools::compute_no_normal_flux_constraints(stokes_dof_handler,
 *                                                       0,
 *                                                       no_normal_flux_boundaries,
 *                                                       stokes_constraints,
 *                                                       mapping);
 *       stokes_constraints.close();
 *     }
 *     {
 *       temperature_constraints.clear();
 *       temperature_constraints.reinit(temperature_relevant_partitioning);
 * 
 *       DoFTools::make_hanging_node_constraints(temperature_dof_handler,
 *                                               temperature_constraints);
 *       VectorTools::interpolate_boundary_values(
 *         temperature_dof_handler,
 *         0,
 *         EquationData::TemperatureInitialValues<dim>(),
 *         temperature_constraints);
 *       VectorTools::interpolate_boundary_values(
 *         temperature_dof_handler,
 *         1,
 *         EquationData::TemperatureInitialValues<dim>(),
 *         temperature_constraints);
 *       temperature_constraints.close();
 *     }
 * 
 * @endcode
 * 
 * All this done, we can then initialize the various matrix and vector
 * objects to their proper sizes. At the end, we also record that all
 * matrices and preconditioners have to be re-computed at the beginning of
 * the next time step. Note how we initialize the vectors for the Stokes
 * and temperature right hand sides: These are writable vectors (last
 * boolean argument set to @p true) that have the correct one-to-one
 * partitioning of locally owned elements but are still given the relevant
 * partitioning for means of figuring out the vector entries that are
 * going to be set right away. As for matrices, this allows for writing
 * local contributions into the vector with multiple threads (always
 * assuming that the same vector entry is not accessed by multiple threads
 * at the same time). The other vectors only allow for read access of
 * individual elements, including ghosts, but are not suitable for
 * solvers.
 * 
 * @code
 *     setup_stokes_matrix(stokes_partitioning, stokes_relevant_partitioning);
 *     setup_stokes_preconditioner(stokes_partitioning,
 *                                 stokes_relevant_partitioning);
 *     setup_temperature_matrices(temperature_partitioning,
 *                                temperature_relevant_partitioning);
 * 
 *     stokes_rhs.reinit(stokes_partitioning,
 *                       stokes_relevant_partitioning,
 *                       MPI_COMM_WORLD,
 *                       true);
 *     stokes_solution.reinit(stokes_relevant_partitioning, MPI_COMM_WORLD);
 *     old_stokes_solution.reinit(stokes_solution);
 * 
 *     temperature_rhs.reinit(temperature_partitioning,
 *                            temperature_relevant_partitioning,
 *                            MPI_COMM_WORLD,
 *                            true);
 *     temperature_solution.reinit(temperature_relevant_partitioning,
 *                                 MPI_COMM_WORLD);
 *     old_temperature_solution.reinit(temperature_solution);
 *     old_old_temperature_solution.reinit(temperature_solution);
 * 
 *     rebuild_stokes_matrix              = true;
 *     rebuild_stokes_preconditioner      = true;
 *     rebuild_temperature_matrices       = true;
 *     rebuild_temperature_preconditioner = true;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheBoussinesqFlowProblemassemblyfunctions"></a> 
 * <h4>The BoussinesqFlowProblem assembly functions</h4>
 *   

 * 
 * Following the discussion in the introduction and in the @ref threads
 * module, we split the assembly functions into different parts:
 *   

 * 
 * <ul> <li> The local calculations of matrices and right hand sides, given
 * a certain cell as input (these functions are named
 * <code>local_assemble_*</code> below). The resulting function is, in other
 * words, essentially the body of the loop over all cells in step-31. Note,
 * however, that these functions store the result from the local
 * calculations in variables of classes from the CopyData namespace.
 *   

 * 
 * <li>These objects are then given to the second step which writes the
 * local data into the global data structures (these functions are named
 * <code>copy_local_to_global_*</code> below). These functions are pretty
 * trivial.
 *   

 * 
 * <li>These two subfunctions are then used in the respective assembly
 * routine (called <code>assemble_*</code> below), where a WorkStream object
 * is set up and runs over all the cells that belong to the processor's
 * subdomain.  </ul>
 * 

 * 
 * 
 * <a name="Stokespreconditionerassembly"></a> 
 * <h5>Stokes preconditioner assembly</h5>
 *   

 * 
 * Let us start with the functions that builds the Stokes
 * preconditioner. The first two of these are pretty trivial, given the
 * discussion above. Note in particular that the main point in using the
 * scratch data object is that we want to avoid allocating any objects on
 * the free space each time we visit a new cell. As a consequence, the
 * assembly function below only has automatic local variables, and
 * everything else is accessed through the scratch data object, which is
 * allocated only once before we start the loop over all cells:
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::local_assemble_stokes_preconditioner(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     Assembly::Scratch::StokesPreconditioner<dim> &        scratch,
 *     Assembly::CopyData::StokesPreconditioner<dim> &       data)
 *   {
 *     const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell();
 *     const unsigned int n_q_points =
 *       scratch.stokes_fe_values.n_quadrature_points;
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 *     scratch.stokes_fe_values.reinit(cell);
 *     cell->get_dof_indices(data.local_dof_indices);
 * 
 *     data.local_matrix = 0;
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *           {
 *             scratch.grad_phi_u[k] =
 *               scratch.stokes_fe_values[velocities].gradient(k, q);
 *             scratch.phi_p[k] = scratch.stokes_fe_values[pressure].value(k, q);
 *           }
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             data.local_matrix(i, j) +=
 *               (EquationData::eta *
 *                  scalar_product(scratch.grad_phi_u[i], scratch.grad_phi_u[j]) +
 *                (1. / EquationData::eta) * EquationData::pressure_scaling *
 *                  EquationData::pressure_scaling *
 *                  (scratch.phi_p[i] * scratch.phi_p[j])) *
 *               scratch.stokes_fe_values.JxW(q);
 *       }
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::copy_local_to_global_stokes_preconditioner(
 *     const Assembly::CopyData::StokesPreconditioner<dim> &data)
 *   {
 *     stokes_constraints.distribute_local_to_global(data.local_matrix,
 *                                                   data.local_dof_indices,
 *                                                   stokes_preconditioner_matrix);
 *   }
 * 
 * 
 * @endcode
 * 
 * Now for the function that actually puts things together, using the
 * WorkStream functions.  WorkStream::run needs a start and end iterator to
 * enumerate the cells it is supposed to work on. Typically, one would use
 * DoFHandler::begin_active() and DoFHandler::end() for that but here we
 * actually only want the subset of cells that in fact are owned by the
 * current processor. This is where the FilteredIterator class comes into
 * play: you give it a range of cells and it provides an iterator that only
 * iterates over that subset of cells that satisfy a certain predicate (a
 * predicate is a function of one argument that either returns true or
 * false). The predicate we use here is IteratorFilters::LocallyOwnedCell,
 * i.e., it returns true exactly if the cell is owned by the current
 * processor. The resulting iterator range is then exactly what we need.
 *   

 * 
 * With this obstacle out of the way, we call the WorkStream::run
 * function with this set of cells, scratch and copy objects, and
 * with pointers to two functions: the local assembly and
 * copy-local-to-global function. These functions need to have very
 * specific signatures: three arguments in the first and one
 * argument in the latter case (see the documentation of the
 * WorkStream::run function for the meaning of these arguments).
 * Note how we use a lambda functions to
 * create a function object that satisfies this requirement. It uses
 * function arguments for the local assembly function that specify
 * cell, scratch data, and copy data, as well as function argument
 * for the copy function that expects the
 * data to be written into the global matrix (also see the discussion in
 * step-13's <code>assemble_linear_system()</code> function). On the other
 * hand, the implicit zeroth argument of member functions (namely
 * the <code>this</code> pointer of the object on which that member
 * function is to operate on) is <i>bound</i> to the
 * <code>this</code> pointer of the current function and is captured. The
 * WorkStream::run function, as a consequence, does not need to know
 * anything about the object these functions work on.
 *   

 * 
 * When the WorkStream is executed, it will create several local assembly
 * routines of the first kind for several cells and let some available
 * processors work on them. The function that needs to be synchronized,
 * i.e., the write operation into the global matrix, however, is executed by
 * only one thread at a time in the prescribed order. Of course, this only
 * holds for the parallelization on a single MPI process. Different MPI
 * processes will have their own WorkStream objects and do that work
 * completely independently (and in different memory spaces). In a
 * distributed calculation, some data will accumulate at degrees of freedom
 * that are not owned by the respective processor. It would be inefficient
 * to send data around every time we encounter such a dof. What happens
 * instead is that the Trilinos sparse matrix will keep that data and send
 * it to the owner at the end of assembly, by calling the
 * <code>compress()</code> command.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner()
 *   {
 *     stokes_preconditioner_matrix = 0;
 * 
 *     const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree + 1);
 * 
 *     using CellFilter =
 *       FilteredIterator<typename DoFHandler<2>::active_cell_iterator>;
 * 
 *     auto worker =
 *       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *              Assembly::Scratch::StokesPreconditioner<dim> &        scratch,
 *              Assembly::CopyData::StokesPreconditioner<dim> &       data) {
 *         this->local_assemble_stokes_preconditioner(cell, scratch, data);
 *       };
 * 
 *     auto copier =
 *       [this](const Assembly::CopyData::StokesPreconditioner<dim> &data) {
 *         this->copy_local_to_global_stokes_preconditioner(data);
 *       };
 * 
 *     WorkStream::run(CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                                stokes_dof_handler.begin_active()),
 *                     CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                                stokes_dof_handler.end()),
 *                     worker,
 *                     copier,
 *                     Assembly::Scratch::StokesPreconditioner<dim>(
 *                       stokes_fe,
 *                       quadrature_formula,
 *                       mapping,
 *                       update_JxW_values | update_values | update_gradients),
 *                     Assembly::CopyData::StokesPreconditioner<dim>(stokes_fe));
 * 
 *     stokes_preconditioner_matrix.compress(VectorOperation::add);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The final function in this block initiates assembly of the Stokes
 * preconditioner matrix and then in fact builds the Stokes
 * preconditioner. It is mostly the same as in the serial case. The only
 * difference to step-31 is that we use a Jacobi preconditioner for the
 * pressure mass matrix instead of IC, as discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::build_stokes_preconditioner()
 *   {
 *     if (rebuild_stokes_preconditioner == false)
 *       return;
 * 
 *     TimerOutput::Scope timer_section(computing_timer,
 *                                      "   Build Stokes preconditioner");
 *     pcout << "   Rebuilding Stokes preconditioner..." << std::flush;
 * 
 *     assemble_stokes_preconditioner();
 * 
 *     std::vector<std::vector<bool>> constant_modes;
 *     FEValuesExtractors::Vector     velocity_components(0);
 *     DoFTools::extract_constant_modes(stokes_dof_handler,
 *                                      stokes_fe.component_mask(
 *                                        velocity_components),
 *                                      constant_modes);
 * 
 *     Mp_preconditioner =
 *       std::make_shared<TrilinosWrappers::PreconditionJacobi>();
 *     Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
 * 
 *     TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
 *     Amg_data.constant_modes        = constant_modes;
 *     Amg_data.elliptic              = true;
 *     Amg_data.higher_order_elements = true;
 *     Amg_data.smoother_sweeps       = 2;
 *     Amg_data.aggregation_threshold = 0.02;
 * 
 *     Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1, 1));
 *     Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0, 0),
 *                                    Amg_data);
 * 
 *     rebuild_stokes_preconditioner = false;
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Stokessystemassembly"></a> 
 * <h5>Stokes system assembly</h5>
 * 

 * 
 * The next three functions implement the assembly of the Stokes system,
 * again split up into a part performing local calculations, one for writing
 * the local data into the global matrix and vector, and one for actually
 * running the loop over all cells with the help of the WorkStream
 * class. Note that the assembly of the Stokes matrix needs only to be done
 * in case we have changed the mesh. Otherwise, just the
 * (temperature-dependent) right hand side needs to be calculated
 * here. Since we are working with distributed matrices and vectors, we have
 * to call the respective <code>compress()</code> functions in the end of
 * the assembly in order to send non-local data to the owner process.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::local_assemble_stokes_system(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     Assembly::Scratch::StokesSystem<dim> &                scratch,
 *     Assembly::CopyData::StokesSystem<dim> &               data)
 *   {
 *     const unsigned int dofs_per_cell =
 *       scratch.stokes_fe_values.get_fe().n_dofs_per_cell();
 *     const unsigned int n_q_points =
 *       scratch.stokes_fe_values.n_quadrature_points;
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure(dim);
 * 
 *     scratch.stokes_fe_values.reinit(cell);
 * 
 *     typename DoFHandler<dim>::active_cell_iterator temperature_cell(
 *       &triangulation, cell->level(), cell->index(), &temperature_dof_handler);
 *     scratch.temperature_fe_values.reinit(temperature_cell);
 * 
 *     if (rebuild_stokes_matrix)
 *       data.local_matrix = 0;
 *     data.local_rhs = 0;
 * 
 *     scratch.temperature_fe_values.get_function_values(
 *       old_temperature_solution, scratch.old_temperature_values);
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         const double old_temperature = scratch.old_temperature_values[q];
 * 
 *         for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *           {
 *             scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value(k, q);
 *             if (rebuild_stokes_matrix)
 *               {
 *                 scratch.grads_phi_u[k] =
 *                   scratch.stokes_fe_values[velocities].symmetric_gradient(k, q);
 *                 scratch.div_phi_u[k] =
 *                   scratch.stokes_fe_values[velocities].divergence(k, q);
 *                 scratch.phi_p[k] =
 *                   scratch.stokes_fe_values[pressure].value(k, q);
 *               }
 *           }
 * 
 *         if (rebuild_stokes_matrix == true)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               data.local_matrix(i, j) +=
 *                 (EquationData::eta * 2 *
 *                    (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]) -
 *                  (EquationData::pressure_scaling * scratch.div_phi_u[i] *
 *                   scratch.phi_p[j]) -
 *                  (EquationData::pressure_scaling * scratch.phi_p[i] *
 *                   scratch.div_phi_u[j])) *
 *                 scratch.stokes_fe_values.JxW(q);
 * 
 *         const Tensor<1, dim> gravity = EquationData::gravity_vector(
 *           scratch.stokes_fe_values.quadrature_point(q));
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           data.local_rhs(i) += (EquationData::density(old_temperature) *
 *                                 gravity * scratch.phi_u[i]) *
 *                                scratch.stokes_fe_values.JxW(q);
 *       }
 * 
 *     cell->get_dof_indices(data.local_dof_indices);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::copy_local_to_global_stokes_system(
 *     const Assembly::CopyData::StokesSystem<dim> &data)
 *   {
 *     if (rebuild_stokes_matrix == true)
 *       stokes_constraints.distribute_local_to_global(data.local_matrix,
 *                                                     data.local_rhs,
 *                                                     data.local_dof_indices,
 *                                                     stokes_matrix,
 *                                                     stokes_rhs);
 *     else
 *       stokes_constraints.distribute_local_to_global(data.local_rhs,
 *                                                     data.local_dof_indices,
 *                                                     stokes_rhs);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_system()
 *   {
 *     TimerOutput::Scope timer_section(computing_timer,
 *                                      "   Assemble Stokes system");
 * 
 *     if (rebuild_stokes_matrix == true)
 *       stokes_matrix = 0;
 * 
 *     stokes_rhs = 0;
 * 
 *     const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree + 1);
 * 
 *     using CellFilter =
 *       FilteredIterator<typename DoFHandler<2>::active_cell_iterator>;
 * 
 *     WorkStream::run(
 *       CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                  stokes_dof_handler.begin_active()),
 *       CellFilter(IteratorFilters::LocallyOwnedCell(), stokes_dof_handler.end()),
 *       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *              Assembly::Scratch::StokesSystem<dim> &                scratch,
 *              Assembly::CopyData::StokesSystem<dim> &               data) {
 *         this->local_assemble_stokes_system(cell, scratch, data);
 *       },
 *       [this](const Assembly::CopyData::StokesSystem<dim> &data) {
 *         this->copy_local_to_global_stokes_system(data);
 *       },
 *       Assembly::Scratch::StokesSystem<dim>(
 *         stokes_fe,
 *         mapping,
 *         quadrature_formula,
 *         (update_values | update_quadrature_points | update_JxW_values |
 *          (rebuild_stokes_matrix == true ? update_gradients : UpdateFlags(0))),
 *         temperature_fe,
 *         update_values),
 *       Assembly::CopyData::StokesSystem<dim>(stokes_fe));
 * 
 *     if (rebuild_stokes_matrix == true)
 *       stokes_matrix.compress(VectorOperation::add);
 *     stokes_rhs.compress(VectorOperation::add);
 * 
 *     rebuild_stokes_matrix = false;
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Temperaturematrixassembly"></a> 
 * <h5>Temperature matrix assembly</h5>
 * 

 * 
 * The task to be performed by the next three functions is to calculate a
 * mass matrix and a Laplace matrix on the temperature system. These will be
 * combined in order to yield the semi-implicit time stepping matrix that
 * consists of the mass matrix plus a time step-dependent weight factor
 * times the Laplace matrix. This function is again essentially the body of
 * the loop over all cells from step-31.
 *   

 * 
 * The two following functions perform similar services as the ones above.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::local_assemble_temperature_matrix(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     Assembly::Scratch::TemperatureMatrix<dim> &           scratch,
 *     Assembly::CopyData::TemperatureMatrix<dim> &          data)
 *   {
 *     const unsigned int dofs_per_cell =
 *       scratch.temperature_fe_values.get_fe().n_dofs_per_cell();
 *     const unsigned int n_q_points =
 *       scratch.temperature_fe_values.n_quadrature_points;
 * 
 *     scratch.temperature_fe_values.reinit(cell);
 *     cell->get_dof_indices(data.local_dof_indices);
 * 
 *     data.local_mass_matrix      = 0;
 *     data.local_stiffness_matrix = 0;
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *           {
 *             scratch.grad_phi_T[k] =
 *               scratch.temperature_fe_values.shape_grad(k, q);
 *             scratch.phi_T[k] = scratch.temperature_fe_values.shape_value(k, q);
 *           }
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             {
 *               data.local_mass_matrix(i, j) +=
 *                 (scratch.phi_T[i] * scratch.phi_T[j] *
 *                  scratch.temperature_fe_values.JxW(q));
 *               data.local_stiffness_matrix(i, j) +=
 *                 (EquationData::kappa * scratch.grad_phi_T[i] *
 *                  scratch.grad_phi_T[j] * scratch.temperature_fe_values.JxW(q));
 *             }
 *       }
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::copy_local_to_global_temperature_matrix(
 *     const Assembly::CopyData::TemperatureMatrix<dim> &data)
 *   {
 *     temperature_constraints.distribute_local_to_global(data.local_mass_matrix,
 *                                                        data.local_dof_indices,
 *                                                        temperature_mass_matrix);
 *     temperature_constraints.distribute_local_to_global(
 *       data.local_stiffness_matrix,
 *       data.local_dof_indices,
 *       temperature_stiffness_matrix);
 *   }
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::assemble_temperature_matrix()
 *   {
 *     if (rebuild_temperature_matrices == false)
 *       return;
 * 
 *     TimerOutput::Scope timer_section(computing_timer,
 *                                      "   Assemble temperature matrices");
 *     temperature_mass_matrix      = 0;
 *     temperature_stiffness_matrix = 0;
 * 
 *     const QGauss<dim> quadrature_formula(parameters.temperature_degree + 2);
 * 
 *     using CellFilter =
 *       FilteredIterator<typename DoFHandler<2>::active_cell_iterator>;
 * 
 *     WorkStream::run(
 *       CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                  temperature_dof_handler.begin_active()),
 *       CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                  temperature_dof_handler.end()),
 *       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *              Assembly::Scratch::TemperatureMatrix<dim> &           scratch,
 *              Assembly::CopyData::TemperatureMatrix<dim> &          data) {
 *         this->local_assemble_temperature_matrix(cell, scratch, data);
 *       },
 *       [this](const Assembly::CopyData::TemperatureMatrix<dim> &data) {
 *         this->copy_local_to_global_temperature_matrix(data);
 *       },
 *       Assembly::Scratch::TemperatureMatrix<dim>(temperature_fe,
 *                                                 mapping,
 *                                                 quadrature_formula),
 *       Assembly::CopyData::TemperatureMatrix<dim>(temperature_fe));
 * 
 *     temperature_mass_matrix.compress(VectorOperation::add);
 *     temperature_stiffness_matrix.compress(VectorOperation::add);
 * 
 *     rebuild_temperature_matrices       = false;
 *     rebuild_temperature_preconditioner = true;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Temperaturerighthandsideassembly"></a> 
 * <h5>Temperature right hand side assembly</h5>
 * 

 * 
 * This is the last assembly function. It calculates the right hand side of
 * the temperature system, which includes the convection and the
 * stabilization terms. It includes a lot of evaluations of old solutions at
 * the quadrature points (which are necessary for calculating the artificial
 * viscosity of stabilization), but is otherwise similar to the other
 * assembly functions. Notice, once again, how we resolve the dilemma of
 * having inhomogeneous boundary conditions, by just making a right hand
 * side at this point (compare the comments for the <code>project()</code>
 * function above): We create some matrix columns with exactly the values
 * that would be entered for the temperature stiffness matrix, in case we
 * have inhomogeneously constrained dofs. That will account for the correct
 * balance of the right hand side vector with the matrix system of
 * temperature.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::local_assemble_temperature_rhs(
 *     const std::pair<double, double> global_T_range,
 *     const double                    global_max_velocity,
 *     const double                    global_entropy_variation,
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     Assembly::Scratch::TemperatureRHS<dim> &              scratch,
 *     Assembly::CopyData::TemperatureRHS<dim> &             data)
 *   {
 *     const bool use_bdf2_scheme = (timestep_number != 0);
 * 
 *     const unsigned int dofs_per_cell =
 *       scratch.temperature_fe_values.get_fe().n_dofs_per_cell();
 *     const unsigned int n_q_points =
 *       scratch.temperature_fe_values.n_quadrature_points;
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     data.local_rhs     = 0;
 *     data.matrix_for_bc = 0;
 *     cell->get_dof_indices(data.local_dof_indices);
 * 
 *     scratch.temperature_fe_values.reinit(cell);
 * 
 *     typename DoFHandler<dim>::active_cell_iterator stokes_cell(
 *       &triangulation, cell->level(), cell->index(), &stokes_dof_handler);
 *     scratch.stokes_fe_values.reinit(stokes_cell);
 * 
 *     scratch.temperature_fe_values.get_function_values(
 *       old_temperature_solution, scratch.old_temperature_values);
 *     scratch.temperature_fe_values.get_function_values(
 *       old_old_temperature_solution, scratch.old_old_temperature_values);
 * 
 *     scratch.temperature_fe_values.get_function_gradients(
 *       old_temperature_solution, scratch.old_temperature_grads);
 *     scratch.temperature_fe_values.get_function_gradients(
 *       old_old_temperature_solution, scratch.old_old_temperature_grads);
 * 
 *     scratch.temperature_fe_values.get_function_laplacians(
 *       old_temperature_solution, scratch.old_temperature_laplacians);
 *     scratch.temperature_fe_values.get_function_laplacians(
 *       old_old_temperature_solution, scratch.old_old_temperature_laplacians);
 * 
 *     scratch.stokes_fe_values[velocities].get_function_values(
 *       stokes_solution, scratch.old_velocity_values);
 *     scratch.stokes_fe_values[velocities].get_function_values(
 *       old_stokes_solution, scratch.old_old_velocity_values);
 *     scratch.stokes_fe_values[velocities].get_function_symmetric_gradients(
 *       stokes_solution, scratch.old_strain_rates);
 *     scratch.stokes_fe_values[velocities].get_function_symmetric_gradients(
 *       old_stokes_solution, scratch.old_old_strain_rates);
 * 
 *     const double nu =
 *       compute_viscosity(scratch.old_temperature_values,
 *                         scratch.old_old_temperature_values,
 *                         scratch.old_temperature_grads,
 *                         scratch.old_old_temperature_grads,
 *                         scratch.old_temperature_laplacians,
 *                         scratch.old_old_temperature_laplacians,
 *                         scratch.old_velocity_values,
 *                         scratch.old_old_velocity_values,
 *                         scratch.old_strain_rates,
 *                         scratch.old_old_strain_rates,
 *                         global_max_velocity,
 *                         global_T_range.second - global_T_range.first,
 *                         0.5 * (global_T_range.second + global_T_range.first),
 *                         global_entropy_variation,
 *                         cell->diameter());
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q)
 *       {
 *         for (unsigned int k = 0; k < dofs_per_cell; ++k)
 *           {
 *             scratch.phi_T[k] = scratch.temperature_fe_values.shape_value(k, q);
 *             scratch.grad_phi_T[k] =
 *               scratch.temperature_fe_values.shape_grad(k, q);
 *           }
 * 
 * 
 *         const double T_term_for_rhs =
 *           (use_bdf2_scheme ?
 *              (scratch.old_temperature_values[q] *
 *                 (1 + time_step / old_time_step) -
 *               scratch.old_old_temperature_values[q] * (time_step * time_step) /
 *                 (old_time_step * (time_step + old_time_step))) :
 *              scratch.old_temperature_values[q]);
 * 
 *         const double ext_T =
 *           (use_bdf2_scheme ? (scratch.old_temperature_values[q] *
 *                                 (1 + time_step / old_time_step) -
 *                               scratch.old_old_temperature_values[q] *
 *                                 time_step / old_time_step) :
 *                              scratch.old_temperature_values[q]);
 * 
 *         const Tensor<1, dim> ext_grad_T =
 *           (use_bdf2_scheme ? (scratch.old_temperature_grads[q] *
 *                                 (1 + time_step / old_time_step) -
 *                               scratch.old_old_temperature_grads[q] * time_step /
 *                                 old_time_step) :
 *                              scratch.old_temperature_grads[q]);
 * 
 *         const Tensor<1, dim> extrapolated_u =
 *           (use_bdf2_scheme ?
 *              (scratch.old_velocity_values[q] * (1 + time_step / old_time_step) -
 *               scratch.old_old_velocity_values[q] * time_step / old_time_step) :
 *              scratch.old_velocity_values[q]);
 * 
 *         const SymmetricTensor<2, dim> extrapolated_strain_rate =
 *           (use_bdf2_scheme ?
 *              (scratch.old_strain_rates[q] * (1 + time_step / old_time_step) -
 *               scratch.old_old_strain_rates[q] * time_step / old_time_step) :
 *              scratch.old_strain_rates[q]);
 * 
 *         const double gamma =
 *           ((EquationData::radiogenic_heating * EquationData::density(ext_T) +
 *             2 * EquationData::eta * extrapolated_strain_rate *
 *               extrapolated_strain_rate) /
 *            (EquationData::density(ext_T) * EquationData::specific_heat));
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           {
 *             data.local_rhs(i) +=
 *               (T_term_for_rhs * scratch.phi_T[i] -
 *                time_step * extrapolated_u * ext_grad_T * scratch.phi_T[i] -
 *                time_step * nu * ext_grad_T * scratch.grad_phi_T[i] +
 *                time_step * gamma * scratch.phi_T[i]) *
 *               scratch.temperature_fe_values.JxW(q);
 * 
 *             if (temperature_constraints.is_inhomogeneously_constrained(
 *                   data.local_dof_indices[i]))
 *               {
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                   data.matrix_for_bc(j, i) +=
 *                     (scratch.phi_T[i] * scratch.phi_T[j] *
 *                        (use_bdf2_scheme ? ((2 * time_step + old_time_step) /
 *                                            (time_step + old_time_step)) :
 *                                           1.) +
 *                      scratch.grad_phi_T[i] * scratch.grad_phi_T[j] *
 *                        EquationData::kappa * time_step) *
 *                     scratch.temperature_fe_values.JxW(q);
 *               }
 *           }
 *       }
 *   }
 * 
 * 
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::copy_local_to_global_temperature_rhs(
 *     const Assembly::CopyData::TemperatureRHS<dim> &data)
 *   {
 *     temperature_constraints.distribute_local_to_global(data.local_rhs,
 *                                                        data.local_dof_indices,
 *                                                        temperature_rhs,
 *                                                        data.matrix_for_bc);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In the function that runs the WorkStream for actually calculating the
 * right hand side, we also generate the final matrix. As mentioned above,
 * it is a sum of the mass matrix and the Laplace matrix, times some time
 * step-dependent weight. This weight is specified by the BDF-2 time
 * integration scheme, see the introduction in step-31. What is new in this
 * tutorial program (in addition to the use of MPI parallelization and the
 * WorkStream class), is that we now precompute the temperature
 * preconditioner as well. The reason is that the setup of the Jacobi
 * preconditioner takes a noticeable time compared to the solver because we
 * usually only need between 10 and 20 iterations for solving the
 * temperature system (this might sound strange, as Jacobi really only
 * consists of a diagonal, but in Trilinos it is derived from more general
 * framework for point relaxation preconditioners which is a bit
 * inefficient). Hence, it is more efficient to precompute the
 * preconditioner, even though the matrix entries may slightly change
 * because the time step might change. This is not too big a problem because
 * we remesh every few time steps (and regenerate the preconditioner then).
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
 *     if (rebuild_temperature_preconditioner == true)
 *       {
 *         T_preconditioner =
 *           std::make_shared<TrilinosWrappers::PreconditionJacobi>();
 *         T_preconditioner->initialize(temperature_matrix);
 *         rebuild_temperature_preconditioner = false;
 *       }
 * 
 * @endcode
 * 
 * The next part is computing the right hand side vectors.  To do so, we
 * first compute the average temperature $T_m$ that we use for evaluating
 * the artificial viscosity stabilization through the residual $E(T) =
 * (T-T_m)^2$. We do this by defining the midpoint between maximum and
 * minimum temperature as average temperature in the definition of the
 * entropy viscosity. An alternative would be to use the integral average,
 * but the results are not very sensitive to this choice. The rest then
 * only requires calling WorkStream::run again, binding the arguments to
 * the <code>local_assemble_temperature_rhs</code> function that are the
 * same in every call to the correct values:
 * 
 * @code
 *     temperature_rhs = 0;
 * 
 *     const QGauss<dim> quadrature_formula(parameters.temperature_degree + 2);
 *     const std::pair<double, double> global_T_range =
 *       get_extrapolated_temperature_range();
 * 
 *     const double average_temperature =
 *       0.5 * (global_T_range.first + global_T_range.second);
 *     const double global_entropy_variation =
 *       get_entropy_variation(average_temperature);
 * 
 *     using CellFilter =
 *       FilteredIterator<typename DoFHandler<2>::active_cell_iterator>;
 * 
 *     auto worker =
 *       [this, global_T_range, maximal_velocity, global_entropy_variation](
 *         const typename DoFHandler<dim>::active_cell_iterator &cell,
 *         Assembly::Scratch::TemperatureRHS<dim> &              scratch,
 *         Assembly::CopyData::TemperatureRHS<dim> &             data) {
 *         this->local_assemble_temperature_rhs(global_T_range,
 *                                              maximal_velocity,
 *                                              global_entropy_variation,
 *                                              cell,
 *                                              scratch,
 *                                              data);
 *       };
 * 
 *     auto copier = [this](const Assembly::CopyData::TemperatureRHS<dim> &data) {
 *       this->copy_local_to_global_temperature_rhs(data);
 *     };
 * 
 *     WorkStream::run(CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                                temperature_dof_handler.begin_active()),
 *                     CellFilter(IteratorFilters::LocallyOwnedCell(),
 *                                temperature_dof_handler.end()),
 *                     worker,
 *                     copier,
 *                     Assembly::Scratch::TemperatureRHS<dim>(
 *                       temperature_fe, stokes_fe, mapping, quadrature_formula),
 *                     Assembly::CopyData::TemperatureRHS<dim>(temperature_fe));
 * 
 *     temperature_rhs.compress(VectorOperation::add);
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
 * This function solves the linear systems in each time step of the
 * Boussinesq problem. First, we work on the Stokes system and then on the
 * temperature system. In essence, it does the same things as the respective
 * function in step-31. However, there are a few changes here.
 *   

 * 
 * The first change is related to the way we store our solution: we keep the
 * vectors with locally owned degrees of freedom plus ghost nodes on each
 * MPI node. When we enter a solver which is supposed to perform
 * matrix-vector products with a distributed matrix, this is not the
 * appropriate form, though. There, we will want to have the solution vector
 * to be distributed in the same way as the matrix, i.e. without any
 * ghosts. So what we do first is to generate a distributed vector called
 * <code>distributed_stokes_solution</code> and put only the locally owned
 * dofs into that, which is neatly done by the <code>operator=</code> of the
 * Trilinos vector.
 *   

 * 
 * Next, we scale the pressure solution (or rather, the initial guess) for
 * the solver so that it matches with the length scales in the matrices, as
 * discussed in the introduction. We also immediately scale the pressure
 * solution back to the correct units after the solution is completed.  We
 * also need to set the pressure values at hanging nodes to zero. This we
 * also did in step-31 in order not to disturb the Schur complement by some
 * vector entries that actually are irrelevant during the solve stage. As a
 * difference to step-31, here we do it only for the locally owned pressure
 * dofs. After solving for the Stokes solution, each processor copies the
 * distributed solution back into the solution vector that also includes
 * ghost elements.
 *   

 * 
 * The third and most obvious change is that we have two variants for the
 * Stokes solver: A fast solver that sometimes breaks down, and a robust
 * solver that is slower. This is what we already discussed in the
 * introduction. Here is how we realize it: First, we perform 30 iterations
 * with the fast solver based on the simple preconditioner based on the AMG
 * V-cycle instead of an approximate solve (this is indicated by the
 * <code>false</code> argument to the
 * <code>LinearSolvers::BlockSchurPreconditioner</code> object). If we
 * converge, everything is fine. If we do not converge, the solver control
 * object will throw an exception SolverControl::NoConvergence. Usually,
 * this would abort the program because we don't catch them in our usual
 * <code>solve()</code> functions. This is certainly not what we want to
 * happen here. Rather, we want to switch to the strong solver and continue
 * the solution process with whatever vector we got so far. Hence, we catch
 * the exception with the C++ try/catch mechanism. We then simply go through
 * the same solver sequence again in the <code>catch</code> clause, this
 * time passing the @p true flag to the preconditioner for the strong
 * solver, signaling an approximate CG solve.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::solve()
 *   {
 *     {
 *       TimerOutput::Scope timer_section(computing_timer,
 *                                        "   Solve Stokes system");
 * 
 *       pcout << "   Solving Stokes system... " << std::flush;
 * 
 *       TrilinosWrappers::MPI::BlockVector distributed_stokes_solution(
 *         stokes_rhs);
 *       distributed_stokes_solution = stokes_solution;
 * 
 *       distributed_stokes_solution.block(1) /= EquationData::pressure_scaling;
 * 
 *       const unsigned int
 *         start = (distributed_stokes_solution.block(0).size() +
 *                  distributed_stokes_solution.block(1).local_range().first),
 *         end   = (distributed_stokes_solution.block(0).size() +
 *                distributed_stokes_solution.block(1).local_range().second);
 *       for (unsigned int i = start; i < end; ++i)
 *         if (stokes_constraints.is_constrained(i))
 *           distributed_stokes_solution(i) = 0;
 * 
 * 
 *       PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;
 * 
 *       unsigned int  n_iterations     = 0;
 *       const double  solver_tolerance = 1e-8 * stokes_rhs.l2_norm();
 *       SolverControl solver_control(30, solver_tolerance);
 * 
 *       try
 *         {
 *           const LinearSolvers::BlockSchurPreconditioner<
 *             TrilinosWrappers::PreconditionAMG,
 *             TrilinosWrappers::PreconditionJacobi>
 *             preconditioner(stokes_matrix,
 *                            stokes_preconditioner_matrix,
 *                            *Mp_preconditioner,
 *                            *Amg_preconditioner,
 *                            false);
 * 
 *           SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(
 *             solver_control,
 *             mem,
 *             SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
 *               30));
 *           solver.solve(stokes_matrix,
 *                        distributed_stokes_solution,
 *                        stokes_rhs,
 *                        preconditioner);
 * 
 *           n_iterations = solver_control.last_step();
 *         }
 * 
 *       catch (SolverControl::NoConvergence &)
 *         {
 *           const LinearSolvers::BlockSchurPreconditioner<
 *             TrilinosWrappers::PreconditionAMG,
 *             TrilinosWrappers::PreconditionJacobi>
 *             preconditioner(stokes_matrix,
 *                            stokes_preconditioner_matrix,
 *                            *Mp_preconditioner,
 *                            *Amg_preconditioner,
 *                            true);
 * 
 *           SolverControl solver_control_refined(stokes_matrix.m(),
 *                                                solver_tolerance);
 *           SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(
 *             solver_control_refined,
 *             mem,
 *             SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(
 *               50));
 *           solver.solve(stokes_matrix,
 *                        distributed_stokes_solution,
 *                        stokes_rhs,
 *                        preconditioner);
 * 
 *           n_iterations =
 *             (solver_control.last_step() + solver_control_refined.last_step());
 *         }
 * 
 * 
 *       stokes_constraints.distribute(distributed_stokes_solution);
 * 
 *       distributed_stokes_solution.block(1) *= EquationData::pressure_scaling;
 * 
 *       stokes_solution = distributed_stokes_solution;
 *       pcout << n_iterations << " iterations." << std::endl;
 *     }
 * 
 * 
 * @endcode
 * 
 * Now let's turn to the temperature part: First, we compute the time step
 * size. We found that we need smaller time steps for 3D than for 2D for
 * the shell geometry. This is because the cells are more distorted in
 * that case (it is the smallest edge length that determines the CFL
 * number). Instead of computing the time step from maximum velocity and
 * minimal mesh size as in step-31, we compute local CFL numbers, i.e., on
 * each cell we compute the maximum velocity times the mesh size, and
 * compute the maximum of them. Hence, we need to choose the factor in
 * front of the time step slightly smaller.
 *     

 * 
 * After temperature right hand side assembly, we solve the linear system
 * for temperature (with fully distributed vectors without any ghosts),
 * apply constraints and copy the vector back to one with ghosts.
 *     

 * 
 * In the end, we extract the temperature range similarly to step-31 to
 * produce some output (for example in order to help us choose the
 * stabilization constants, as discussed in the introduction). The only
 * difference is that we need to exchange maxima over all processors.
 * 
 * @code
 *     {
 *       TimerOutput::Scope timer_section(computing_timer,
 *                                        "   Assemble temperature rhs");
 * 
 *       old_time_step = time_step;
 * 
 *       const double scaling = (dim == 3 ? 0.25 : 1.0);
 *       time_step            = (scaling / (2.1 * dim * std::sqrt(1. * dim)) /
 *                    (parameters.temperature_degree * get_cfl_number()));
 * 
 *       const double maximal_velocity = get_maximal_velocity();
 *       pcout << "   Maximal velocity: "
 *             << maximal_velocity * EquationData::year_in_seconds * 100
 *             << " cm/year" << std::endl;
 *       pcout << "   "
 *             << "Time step: " << time_step / EquationData::year_in_seconds
 *             << " years" << std::endl;
 * 
 *       temperature_solution = old_temperature_solution;
 *       assemble_temperature_system(maximal_velocity);
 *     }
 * 
 *     {
 *       TimerOutput::Scope timer_section(computing_timer,
 *                                        "   Solve temperature system");
 * 
 *       SolverControl solver_control(temperature_matrix.m(),
 *                                    1e-12 * temperature_rhs.l2_norm());
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);
 * 
 *       TrilinosWrappers::MPI::Vector distributed_temperature_solution(
 *         temperature_rhs);
 *       distributed_temperature_solution = temperature_solution;
 * 
 *       cg.solve(temperature_matrix,
 *                distributed_temperature_solution,
 *                temperature_rhs,
 *                *T_preconditioner);
 * 
 *       temperature_constraints.distribute(distributed_temperature_solution);
 *       temperature_solution = distributed_temperature_solution;
 * 
 *       pcout << "   " << solver_control.last_step()
 *             << " CG iterations for temperature" << std::endl;
 * 
 *       double temperature[2] = {std::numeric_limits<double>::max(),
 *                                -std::numeric_limits<double>::max()};
 *       double global_temperature[2];
 * 
 *       for (unsigned int i =
 *              distributed_temperature_solution.local_range().first;
 *            i < distributed_temperature_solution.local_range().second;
 *            ++i)
 *         {
 *           temperature[0] =
 *             std::min<double>(temperature[0],
 *                              distributed_temperature_solution(i));
 *           temperature[1] =
 *             std::max<double>(temperature[1],
 *                              distributed_temperature_solution(i));
 *         }
 * 
 *       temperature[0] *= -1.0;
 *       Utilities::MPI::max(temperature, MPI_COMM_WORLD, global_temperature);
 *       global_temperature[0] *= -1.0;
 * 
 *       pcout << "   Temperature range: " << global_temperature[0] << ' '
 *             << global_temperature[1] << std::endl;
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemoutput_results"></a> 
 * <h4>BoussinesqFlowProblem::output_results</h4>
 * 

 * 
 * Next comes the function that generates the output. The quantities to
 * output could be introduced manually like we did in step-31. An
 * alternative is to hand this task over to a class PostProcessor that
 * inherits from the class DataPostprocessor, which can be attached to
 * DataOut. This allows us to output derived quantities from the solution,
 * like the friction heating included in this example. It overloads the
 * virtual function DataPostprocessor::evaluate_vector_field(),
 * which is then internally called from DataOut::build_patches(). We have to
 * give it values of the numerical solution, its derivatives, normals to the
 * cell, the actual evaluation points and any additional quantities. This
 * follows the same procedure as discussed in step-29 and other programs.
 * 
 * @code
 *   template <int dim>
 *   class BoussinesqFlowProblem<dim>::Postprocessor
 *     : public DataPostprocessor<dim>
 *   {
 *   public:
 *     Postprocessor(const unsigned int partition, const double minimal_pressure);
 * 
 *     virtual void evaluate_vector_field(
 *       const DataPostprocessorInputs::Vector<dim> &inputs,
 *       std::vector<Vector<double>> &computed_quantities) const override;
 * 
 *     virtual std::vector<std::string> get_names() const override;
 * 
 *     virtual std::vector<
 *       DataComponentInterpretation::DataComponentInterpretation>
 *     get_data_component_interpretation() const override;
 * 
 *     virtual UpdateFlags get_needed_update_flags() const override;
 * 
 *   private:
 *     const unsigned int partition;
 *     const double       minimal_pressure;
 *   };
 * 
 * 
 *   template <int dim>
 *   BoussinesqFlowProblem<dim>::Postprocessor::Postprocessor(
 *     const unsigned int partition,
 *     const double       minimal_pressure)
 *     : partition(partition)
 *     , minimal_pressure(minimal_pressure)
 *   {}
 * 
 * 
 * @endcode
 * 
 * Here we define the names for the variables we want to output. These are
 * the actual solution values for velocity, pressure, and temperature, as
 * well as the friction heating and to each cell the number of the processor
 * that owns it. This allows us to visualize the partitioning of the domain
 * among the processors. Except for the velocity, which is vector-valued,
 * all other quantities are scalar.
 * 
 * @code
 *   template <int dim>
 *   std::vector<std::string>
 *   BoussinesqFlowProblem<dim>::Postprocessor::get_names() const
 *   {
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     solution_names.emplace_back("p");
 *     solution_names.emplace_back("T");
 *     solution_names.emplace_back("friction_heating");
 *     solution_names.emplace_back("partition");
 * 
 *     return solution_names;
 *   }
 * 
 * 
 *   template <int dim>
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *   BoussinesqFlowProblem<dim>::Postprocessor::get_data_component_interpretation()
 *     const
 *   {
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       interpretation(dim,
 *                      DataComponentInterpretation::component_is_part_of_vector);
 * 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
 * 
 *     return interpretation;
 *   }
 * 
 * 
 *   template <int dim>
 *   UpdateFlags
 *   BoussinesqFlowProblem<dim>::Postprocessor::get_needed_update_flags() const
 *   {
 *     return update_values | update_gradients | update_quadrature_points;
 *   }
 * 
 * 
 * @endcode
 * 
 * Now we implement the function that computes the derived quantities. As we
 * also did for the output, we rescale the velocity from its SI units to
 * something more readable, namely cm/year. Next, the pressure is scaled to
 * be between 0 and the maximum pressure. This makes it more easily
 * comparable -- in essence making all pressure variables positive or
 * zero. Temperature is taken as is, and the friction heating is computed as
 * $2 \eta \varepsilon(\mathbf{u}) \cdot \varepsilon(\mathbf{u})$.
 *   

 * 
 * The quantities we output here are more for illustration, rather than for
 * actual scientific value. We come back to this briefly in the results
 * section of this program and explain what one may in fact be interested in.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::Postprocessor::evaluate_vector_field(
 *     const DataPostprocessorInputs::Vector<dim> &inputs,
 *     std::vector<Vector<double>> &               computed_quantities) const
 *   {
 *     const unsigned int n_quadrature_points = inputs.solution_values.size();
 *     Assert(inputs.solution_gradients.size() == n_quadrature_points,
 *            ExcInternalError());
 *     Assert(computed_quantities.size() == n_quadrature_points,
 *            ExcInternalError());
 *     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
 * 
 *     for (unsigned int q = 0; q < n_quadrature_points; ++q)
 *       {
 *         for (unsigned int d = 0; d < dim; ++d)
 *           computed_quantities[q](d) = (inputs.solution_values[q](d) *
 *                                        EquationData::year_in_seconds * 100);
 * 
 *         const double pressure =
 *           (inputs.solution_values[q](dim) - minimal_pressure);
 *         computed_quantities[q](dim) = pressure;
 * 
 *         const double temperature        = inputs.solution_values[q](dim + 1);
 *         computed_quantities[q](dim + 1) = temperature;
 * 
 *         Tensor<2, dim> grad_u;
 *         for (unsigned int d = 0; d < dim; ++d)
 *           grad_u[d] = inputs.solution_gradients[q][d];
 *         const SymmetricTensor<2, dim> strain_rate = symmetrize(grad_u);
 *         computed_quantities[q](dim + 2) =
 *           2 * EquationData::eta * strain_rate * strain_rate;
 * 
 *         computed_quantities[q](dim + 3) = partition;
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * The <code>output_results()</code> function has a similar task to the one
 * in step-31. However, here we are going to demonstrate a different
 * technique on how to merge output from different DoFHandler objects. The
 * way we're going to achieve this recombination is to create a joint
 * DoFHandler that collects both components, the Stokes solution and the
 * temperature solution. This can be nicely done by combining the finite
 * elements from the two systems to form one FESystem, and let this
 * collective system define a new DoFHandler object. To be sure that
 * everything was done correctly, we perform a sanity check that ensures
 * that we got all the dofs from both Stokes and temperature even in the
 * combined system. We then combine the data vectors. Unfortunately, there
 * is no straight-forward relation that tells us how to sort Stokes and
 * temperature vector into the joint vector. The way we can get around this
 * trouble is to rely on the information collected in the FESystem. For each
 * dof on a cell, the joint finite element knows to which equation component
 * (velocity component, pressure, or temperature) it belongs – that's the
 * information we need! So we step through all cells (with iterators into
 * all three DoFHandlers moving in sync), and for each joint cell dof, we
 * read out that component using the FiniteElement::system_to_base_index
 * function (see there for a description of what the various parts of its
 * return value contain). We also need to keep track whether we're on a
 * Stokes dof or a temperature dof, which is contained in
 * joint_fe.system_to_base_index(i).first.first. Eventually, the dof_indices
 * data structures on either of the three systems tell us how the relation
 * between global vector and local dofs looks like on the present cell,
 * which concludes this tedious work. We make sure that each processor only
 * works on the subdomain it owns locally (and not on ghost or artificial
 * cells) when building the joint solution vector. The same will then have
 * to be done in DataOut::build_patches(), but that function does so
 * automatically.
 *   

 * 
 * What we end up with is a set of patches that we can write using the
 * functions in DataOutBase in a variety of output formats. Here, we then
 * have to pay attention that what each processor writes is really only its
 * own part of the domain, i.e. we will want to write each processor's
 * contribution into a separate file. This we do by adding an additional
 * number to the filename when we write the solution. This is not really
 * new, we did it similarly in step-40. Note that we write in the compressed
 * format @p .vtu instead of plain vtk files, which saves quite some
 * storage.
 *   

 * 
 * All the rest of the work is done in the PostProcessor class.
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::output_results()
 *   {
 *     TimerOutput::Scope timer_section(computing_timer, "Postprocessing");
 * 
 *     const FESystem<dim> joint_fe(stokes_fe, 1, temperature_fe, 1);
 * 
 *     DoFHandler<dim> joint_dof_handler(triangulation);
 *     joint_dof_handler.distribute_dofs(joint_fe);
 *     Assert(joint_dof_handler.n_dofs() ==
 *              stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(),
 *            ExcInternalError());
 * 
 *     TrilinosWrappers::MPI::Vector joint_solution;
 *     joint_solution.reinit(joint_dof_handler.locally_owned_dofs(),
 *                           MPI_COMM_WORLD);
 * 
 *     {
 *       std::vector<types::global_dof_index> local_joint_dof_indices(
 *         joint_fe.n_dofs_per_cell());
 *       std::vector<types::global_dof_index> local_stokes_dof_indices(
 *         stokes_fe.n_dofs_per_cell());
 *       std::vector<types::global_dof_index> local_temperature_dof_indices(
 *         temperature_fe.n_dofs_per_cell());
 * 
 *       typename DoFHandler<dim>::active_cell_iterator
 *         joint_cell       = joint_dof_handler.begin_active(),
 *         joint_endc       = joint_dof_handler.end(),
 *         stokes_cell      = stokes_dof_handler.begin_active(),
 *         temperature_cell = temperature_dof_handler.begin_active();
 *       for (; joint_cell != joint_endc;
 *            ++joint_cell, ++stokes_cell, ++temperature_cell)
 *         if (joint_cell->is_locally_owned())
 *           {
 *             joint_cell->get_dof_indices(local_joint_dof_indices);
 *             stokes_cell->get_dof_indices(local_stokes_dof_indices);
 *             temperature_cell->get_dof_indices(local_temperature_dof_indices);
 * 
 *             for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i)
 *               if (joint_fe.system_to_base_index(i).first.first == 0)
 *                 {
 *                   Assert(joint_fe.system_to_base_index(i).second <
 *                            local_stokes_dof_indices.size(),
 *                          ExcInternalError());
 * 
 *                   joint_solution(local_joint_dof_indices[i]) = stokes_solution(
 *                     local_stokes_dof_indices[joint_fe.system_to_base_index(i)
 *                                                .second]);
 *                 }
 *               else
 *                 {
 *                   Assert(joint_fe.system_to_base_index(i).first.first == 1,
 *                          ExcInternalError());
 *                   Assert(joint_fe.system_to_base_index(i).second <
 *                            local_temperature_dof_indices.size(),
 *                          ExcInternalError());
 *                   joint_solution(local_joint_dof_indices[i]) =
 *                     temperature_solution(
 *                       local_temperature_dof_indices
 *                         [joint_fe.system_to_base_index(i).second]);
 *                 }
 *           }
 *     }
 * 
 *     joint_solution.compress(VectorOperation::insert);
 * 
 *     IndexSet locally_relevant_joint_dofs(joint_dof_handler.n_dofs());
 *     DoFTools::extract_locally_relevant_dofs(joint_dof_handler,
 *                                             locally_relevant_joint_dofs);
 *     TrilinosWrappers::MPI::Vector locally_relevant_joint_solution;
 *     locally_relevant_joint_solution.reinit(locally_relevant_joint_dofs,
 *                                            MPI_COMM_WORLD);
 *     locally_relevant_joint_solution = joint_solution;
 * 
 *     Postprocessor postprocessor(Utilities::MPI::this_mpi_process(
 *                                   MPI_COMM_WORLD),
 *                                 stokes_solution.block(1).min());
 * 
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(joint_dof_handler);
 *     data_out.add_data_vector(locally_relevant_joint_solution, postprocessor);
 *     data_out.build_patches();
 * 
 *     static int out_index = 0;
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", out_index, MPI_COMM_WORLD, 5);
 * 
 *     out_index++;
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
 * This function isn't really new either. Since the <code>setup_dofs</code>
 * function that we call in the middle has its own timer section, we split
 * timing this function into two sections. It will also allow us to easily
 * identify which of the two is more expensive.
 *   

 * 
 * One thing of note, however, is that we only want to compute error
 * indicators on the locally owned subdomain. In order to achieve this, we
 * pass one additional argument to the KellyErrorEstimator::estimate
 * function. Note that the vector for error estimates is resized to the
 * number of active cells present on the current process, which is less than
 * the total number of active cells on all processors (but more than the
 * number of locally owned active cells); each processor only has a few
 * coarse cells around the locally owned ones, as also explained in step-40.
 *   

 * 
 * The local error estimates are then handed to a %parallel version of
 * GridRefinement (in namespace parallel::distributed::GridRefinement, see
 * also step-40) which looks at the errors and finds the cells that need
 * refinement by comparing the error values across processors. As in
 * step-31, we want to limit the maximum grid level. So in case some cells
 * have been marked that are already at the finest level, we simply clear
 * the refine flags.
 * 
 * @code
 *   template <int dim>
 *   void
 *   BoussinesqFlowProblem<dim>::refine_mesh(const unsigned int max_grid_level)
 *   {
 *     parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
 *       temperature_trans(temperature_dof_handler);
 *     parallel::distributed::SolutionTransfer<dim,
 *                                             TrilinosWrappers::MPI::BlockVector>
 *       stokes_trans(stokes_dof_handler);
 * 
 *     {
 *       TimerOutput::Scope timer_section(computing_timer,
 *                                        "Refine mesh structure, part 1");
 * 
 *       Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *       KellyErrorEstimator<dim>::estimate(
 *         temperature_dof_handler,
 *         QGauss<dim - 1>(parameters.temperature_degree + 1),
 *         std::map<types::boundary_id, const Function<dim> *>(),
 *         temperature_solution,
 *         estimated_error_per_cell,
 *         ComponentMask(),
 *         nullptr,
 *         0,
 *         triangulation.locally_owned_subdomain());
 * 
 *       parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
 *         triangulation, estimated_error_per_cell, 0.3, 0.1);
 * 
 *       if (triangulation.n_levels() > max_grid_level)
 *         for (typename Triangulation<dim>::active_cell_iterator cell =
 *                triangulation.begin_active(max_grid_level);
 *              cell != triangulation.end();
 *              ++cell)
 *           cell->clear_refine_flag();
 * 
 * @endcode
 * 
 * With all flags marked as necessary, we can then tell the
 * parallel::distributed::SolutionTransfer objects to get ready to
 * transfer data from one mesh to the next, which they will do when
 * notified by
 * Triangulation as part of the @p execute_coarsening_and_refinement() call.
 * The syntax is similar to the non-%parallel solution transfer (with the
 * exception that here a pointer to the vector entries is enough). The
 * remainder of the function further down below is then concerned with
 * setting up the data structures again after mesh refinement and
 * restoring the solution vectors on the new mesh.
 * 
 * @code
 *       std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature(2);
 *       x_temperature[0] = &temperature_solution;
 *       x_temperature[1] = &old_temperature_solution;
 *       std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes(2);
 *       x_stokes[0] = &stokes_solution;
 *       x_stokes[1] = &old_stokes_solution;
 * 
 *       triangulation.prepare_coarsening_and_refinement();
 * 
 *       temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
 *       stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);
 * 
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * 
 *     setup_dofs();
 * 
 *     {
 *       TimerOutput::Scope timer_section(computing_timer,
 *                                        "Refine mesh structure, part 2");
 * 
 *       {
 *         TrilinosWrappers::MPI::Vector distributed_temp1(temperature_rhs);
 *         TrilinosWrappers::MPI::Vector distributed_temp2(temperature_rhs);
 * 
 *         std::vector<TrilinosWrappers::MPI::Vector *> tmp(2);
 *         tmp[0] = &(distributed_temp1);
 *         tmp[1] = &(distributed_temp2);
 *         temperature_trans.interpolate(tmp);
 * 
 * @endcode
 * 
 * enforce constraints to make the interpolated solution conforming on
 * the new mesh:
 * 
 * @code
 *         temperature_constraints.distribute(distributed_temp1);
 *         temperature_constraints.distribute(distributed_temp2);
 * 
 *         temperature_solution     = distributed_temp1;
 *         old_temperature_solution = distributed_temp2;
 *       }
 * 
 *       {
 *         TrilinosWrappers::MPI::BlockVector distributed_stokes(stokes_rhs);
 *         TrilinosWrappers::MPI::BlockVector old_distributed_stokes(stokes_rhs);
 * 
 *         std::vector<TrilinosWrappers::MPI::BlockVector *> stokes_tmp(2);
 *         stokes_tmp[0] = &(distributed_stokes);
 *         stokes_tmp[1] = &(old_distributed_stokes);
 * 
 *         stokes_trans.interpolate(stokes_tmp);
 * 
 * @endcode
 * 
 * enforce constraints to make the interpolated solution conforming on
 * the new mesh:
 * 
 * @code
 *         stokes_constraints.distribute(distributed_stokes);
 *         stokes_constraints.distribute(old_distributed_stokes);
 * 
 *         stokes_solution     = distributed_stokes;
 *         old_stokes_solution = old_distributed_stokes;
 *       }
 *     }
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
 * This is the final and controlling function in this class. It, in fact,
 * runs the entire rest of the program and is, once more, very similar to
 * step-31. The only substantial difference is that we use a different mesh
 * now (a GridGenerator::hyper_shell instead of a simple cube geometry).
 * 
 * @code
 *   template <int dim>
 *   void BoussinesqFlowProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_shell(triangulation,
 *                                Point<dim>(),
 *                                EquationData::R0,
 *                                EquationData::R1,
 *                                (dim == 3) ? 96 : 12,
 *                                true);
 * 
 *     global_Omega_diameter = GridTools::diameter(triangulation);
 * 
 *     triangulation.refine_global(parameters.initial_global_refinement);
 * 
 *     setup_dofs();
 * 
 *     unsigned int pre_refinement_step = 0;
 * 
 *   start_time_iteration:
 * 
 *     {
 *       TrilinosWrappers::MPI::Vector solution(
 *         temperature_dof_handler.locally_owned_dofs());
 * @endcode
 * 
 * VectorTools::project supports parallel vector classes with most
 * standard finite elements via deal.II's own native MatrixFree framework:
 * since we use standard Lagrange elements of moderate order this function
 * works well here.
 * 
 * @code
 *       VectorTools::project(temperature_dof_handler,
 *                            temperature_constraints,
 *                            QGauss<dim>(parameters.temperature_degree + 2),
 *                            EquationData::TemperatureInitialValues<dim>(),
 *                            solution);
 * @endcode
 * 
 * Having so computed the current temperature field, let us set the member
 * variable that holds the temperature nodes. Strictly speaking, we really
 * only need to set <code>old_temperature_solution</code> since the first
 * thing we will do is to compute the Stokes solution that only requires
 * the previous time step's temperature field. That said, nothing good can
 * come from not initializing the other vectors as well (especially since
 * it's a relatively cheap operation and we only have to do it once at the
 * beginning of the program) if we ever want to extend our numerical
 * method or physical model, and so we initialize
 * <code>old_temperature_solution</code> and
 * <code>old_old_temperature_solution</code> as well. The assignment makes
 * sure that the vectors on the left hand side (which where initialized to
 * contain ghost elements as well) also get the correct ghost elements. In
 * other words, the assignment here requires communication between
 * processors:
 * 
 * @code
 *       temperature_solution         = solution;
 *       old_temperature_solution     = solution;
 *       old_old_temperature_solution = solution;
 *     }
 * 
 *     timestep_number = 0;
 *     time_step = old_time_step = 0;
 * 
 *     double time = 0;
 * 
 *     do
 *       {
 *         pcout << "Timestep " << timestep_number
 *               << ":  t=" << time / EquationData::year_in_seconds << " years"
 *               << std::endl;
 * 
 *         assemble_stokes_system();
 *         build_stokes_preconditioner();
 *         assemble_temperature_matrix();
 * 
 *         solve();
 * 
 *         pcout << std::endl;
 * 
 *         if ((timestep_number == 0) &&
 *             (pre_refinement_step < parameters.initial_adaptive_refinement))
 *           {
 *             refine_mesh(parameters.initial_global_refinement +
 *                         parameters.initial_adaptive_refinement);
 *             ++pre_refinement_step;
 *             goto start_time_iteration;
 *           }
 *         else if ((timestep_number > 0) &&
 *                  (timestep_number % parameters.adaptive_refinement_interval ==
 *                   0))
 *           refine_mesh(parameters.initial_global_refinement +
 *                       parameters.initial_adaptive_refinement);
 * 
 *         if ((parameters.generate_graphical_output == true) &&
 *             (timestep_number % parameters.graphical_output_interval == 0))
 *           output_results();
 * 
 * @endcode
 * 
 * In order to speed up linear solvers, we extrapolate the solutions
 * from the old time levels to the new one. This gives a very good
 * initial guess, cutting the number of iterations needed in solvers
 * by more than one half. We do not need to extrapolate in the last
 * iteration, so if we reached the final time, we stop here.
 *         

 * 
 * As the last thing during a time step (before actually bumping up
 * the number of the time step), we check whether the current time
 * step number is divisible by 100, and if so we let the computing
 * timer print a summary of CPU times spent so far.
 * 
 * @code
 *         if (time > parameters.end_time * EquationData::year_in_seconds)
 *           break;
 * 
 *         TrilinosWrappers::MPI::BlockVector old_old_stokes_solution;
 *         old_old_stokes_solution      = old_stokes_solution;
 *         old_stokes_solution          = stokes_solution;
 *         old_old_temperature_solution = old_temperature_solution;
 *         old_temperature_solution     = temperature_solution;
 *         if (old_time_step > 0)
 *           {
 * @endcode
 * 
 * Trilinos sadd does not like ghost vectors even as input. Copy
 * into distributed vectors for now:
 * 
 * @code
 *             {
 *               TrilinosWrappers::MPI::BlockVector distr_solution(stokes_rhs);
 *               distr_solution = stokes_solution;
 *               TrilinosWrappers::MPI::BlockVector distr_old_solution(stokes_rhs);
 *               distr_old_solution = old_old_stokes_solution;
 *               distr_solution.sadd(1. + time_step / old_time_step,
 *                                   -time_step / old_time_step,
 *                                   distr_old_solution);
 *               stokes_solution = distr_solution;
 *             }
 *             {
 *               TrilinosWrappers::MPI::Vector distr_solution(temperature_rhs);
 *               distr_solution = temperature_solution;
 *               TrilinosWrappers::MPI::Vector distr_old_solution(temperature_rhs);
 *               distr_old_solution = old_old_temperature_solution;
 *               distr_solution.sadd(1. + time_step / old_time_step,
 *                                   -time_step / old_time_step,
 *                                   distr_old_solution);
 *               temperature_solution = distr_solution;
 *             }
 *           }
 * 
 *         if ((timestep_number > 0) && (timestep_number % 100 == 0))
 *           computing_timer.print_summary();
 * 
 *         time += time_step;
 *         ++timestep_number;
 *       }
 *     while (true);
 * 
 * @endcode
 * 
 * If we are generating graphical output, do so also for the last time
 * step unless we had just done so before we left the do-while loop
 * 
 * @code
 *     if ((parameters.generate_graphical_output == true) &&
 *         !((timestep_number - 1) % parameters.graphical_output_interval == 0))
 *       output_results();
 *   }
 * } // namespace Step32
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
 * The main function is short as usual and very similar to the one in
 * step-31. Since we use a parameter file which is specified as an argument in
 * the command line, we have to read it in here and pass it on to the
 * Parameters class for parsing. If no filename is given in the command line,
 * we simply use the <code>\step-32.prm</code> file which is distributed
 * together with the program.
 * 

 * 
 * Because 3d computations are simply very slow unless you throw a lot of
 * processors at them, the program defaults to 2d. You can get the 3d version
 * by changing the constant dimension below to 3.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step32;
 *       using namespace dealii;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *         argc, argv, numbers::invalid_unsigned_int);
 * 
 *       std::string parameter_filename;
 *       if (argc >= 2)
 *         parameter_filename = argv[1];
 *       else
 *         parameter_filename = "step-32.prm";
 * 
 *       const int                              dim = 2;
 *       BoussinesqFlowProblem<dim>::Parameters parameters(parameter_filename);
 *       BoussinesqFlowProblem<dim>             flow_problem(parameters);
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
examples/step-32/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当运行时，该程序以与step-31相同的方式模拟三维对流，尽管有一个完全不同的测试案例。




<a name="Comparisonofresultswithstep31"></a><h3>Comparison of results with \step-31</h3>


然而，在我们讨论这个测试案例之前，让我们展示一下这个程序稍早的版本的一些结果，该版本正是在解决我们在第31步中使用的测试案例，只是我们现在以并行方式解决它，而且分辨率要高很多。我们展示这些结果主要是为了比较。

下面是两张图片，如果我们选择 <code>main()</code> 中的3d计算，以及设置 <code>initial_refinement=3</code> 和 <code>n_pre_refinement_steps=4</code> ，则可以看到这种更高的分辨率。在所示的时间步骤中，网格有大约72,000和236,000个单元，分别为2,680,000和8,250,000个自由度，比我们在步骤31中的可用度多了一个数量级。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-32.3d.cube.0.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-32.3d.cube.1.png" alt="">
    </td>
  </tr>
</table> 

计算是在德克萨斯A&amp;M大学Brazos集群的50个处理器的子集上完成的。




<a name="Resultsfora2dcircularshelltestcase"></a><h3>Results for a 2d circular shell testcase</h3>


接下来，我们将用目录中的参数文件运行step-32，但有一个变化：我们将最终时间增加到1e9。这里我们使用的是16个处理器。启动的命令是（注意，step-32.prm是默认的）。

<code> <pre>  $ mpirun -np 16 ./step-32
</pre>
</code>


Note that running a job on a cluster typically requires going through a job
scheduler, which we won't discuss here. The output will look roughly like
this:


<code>
<pre>
\$  mpirun -np 16 ./step-32 活动单元的数量：12,288（在6层） 自由度的数量：186,624（99,840+36,864+49,920）。

时间步数0：t=0年

   重建斯托克斯预处理程序...    解决斯托克斯系统...41次迭代。    最大速度：60.4935厘米/年 时间步长：18166.9年 温度的17次CG迭代 温度范围：973 4273.16

活动单元的数量：15,921（在7层） 自由度的数量：252,723（136,640+47,763+68,320）。

时间步数0：t=0年

   重建斯托克斯预处理程序...    解决斯托克斯系统...50次迭代。    最大速度：60.3223厘米/年 时间步长：10557.6年 温度的19次CG迭代 温度范围：973 4273.16

活动单元的数量：19,926（在8层） 自由度的数量：321,246（174,312+59,778+87,156）。

时间步数0：t=0年

   重建斯托克斯预处理程序...    解决斯托克斯系统...50次迭代。    最大速度：57.8396厘米/年 时间步长：5453.78年 温度的18次CG迭代 温度范围：973 4273.16

时间步数1：t=5453.78年

   解决斯托克斯系统...49次迭代。    最大速度：59.0231厘米/年 时间步长：5345.86年 温度的18次CG迭代 温度范围：973 4273.16

时间步数2：t=10799.6年

   解决斯托克斯系统...24次迭代。    最大速度：60.2139厘米/年 时间步长：5241.51年 温度的17次CG迭代 温度范围：973 4273.16

[...]

时间步数100：t=272151年

   解决斯托克斯系统......21次迭代。    最大速度：161.546厘米/年 时间步长：1672.96年 温度的17次CG迭代 温度范围：973 4282.57

活动单元的数量：56,085（在8层） 自由度的数量：903,408（490,102+168,255+245,051）。




+---------------------------------------------+------------+------------+ | 从开始到现在，总的壁挂时间经过了115s构建斯托克斯预调节器 | 12 | 2.09s | 1.8% | 解算斯托克斯系统 | 103 | 90.4s | 79% | 解算温度系统 | 103 | 1.53s | 1.3% | 后处理 | 3 | 0.532s | 0.完善网格结构，第一部分 | 12 | 0.93s | 0.81% | 完善网格结构，第二部分 | 12 | 0.384s | 0.33% | 设置阻尼系统 | 13 | 2.96s | 2.6% | +---------------------------------+-----------+------------+------------+

[...]

+---------------------------------------------+------------+------------+ | 从开始到现在总共经过了多少壁挂时间 | 9.14e+04s | | | | 部分 | 调用次数 | 壁挂时间 | 占总数的百分比 | +---------------------------------+-----------+------------+------------+ | 组装斯托克斯系统 | 47045 | 2.05e+03s | 2.2% | 组装温度矩阵 | 4707 | 310s | 0.34% | 组装温度rhs | 47045 | 8.7e+03s | 9.4707 | 1.48e+03s | 1.6% | 解决斯托克斯系统 | 47045 | 7.34e+04s | 80% | 解决温度系统 | 47045 | 1.46e+03s | 1.6% | 后处理 | 1883 | 222s | 0.24% | | 完善网格结构，第一部分 | 4706 | 641s | 0.7% | 完善网格结构，第二部分 | 4706 | 259s | 0.28% | 设置阻尼系统 | 4707 | 1.86e+03s | 2% | +---------------------------------+-----------+------------+------------+ </pre> </code>

当时间达到输入文件中选择的10亿年时，模拟就会终止。  你可以从中推断出不同的最终时间的模拟需要多长时间（时间步长最终确定在20,000年左右，所以计算20亿年需要100,000个时间步长，给或给20%）。  从这里可以看出，我们把大部分的计算时间花在了组装线性系统和&mdash;首先&mdash;解决斯托克斯系统。


为了演示输出，我们在这里展示了每1250个时间步骤的输出。   <table>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-000.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-050.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-100.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-150.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-200.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-250.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-300.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-350.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-400.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-450.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-500.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-550.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-time-600.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-cells.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-32-2d-partition.png" alt="">
    </td>
  </tr>
</table> 

最后两张图片显示了网格以及16个子域和16个处理器的同一计算的网格划分情况。这个模拟的全部动态只有通过看动画才能看到，例如<a
href="https://www.dealii.org/images/steps/developer/step-32-2d-temperature.webm">shown
on this site</a>。由于其艺术质量和对岩浆羽流演变的迷人描述，这个图像非常值得观看。

如果你看电影，你会看到对流模式经历了几个阶段。首先，它摆脱了不稳定的温度分层，热物质被致密的冷物质覆盖。在这个巨大的驱动力被消除后，我们有了一种稳定的情况，几个小球开始从内圈的热边界层中分离出来并上升，几个冷指也从外部边界层中掉下来。在这一阶段，解决方案仍然大部分是对称的，反映了原始网格的12倍对称性。在最后一个阶段，流体进入剧烈的混沌搅拌，其中所有的对称性都消失了。这是一个随后继续主导流动的模式。

如果我们看一下模拟中作为时间函数的最大速度，也可以确定这些不同阶段。

 <img src="https://www.dealii.org/images/steps/developer/step-32.2d.t_vs_vmax.png" alt=""> 

在这里，当温度分层不稳定时，速度（以厘米/年表示）在开始时变得非常大，达到几米/年的数量级）。然后平静下来，变成相对较小的数值，然后在混乱的搅动系统中再次回升。在那里，它保持在每年10-40厘米的范围内，完全在物理上预期的区域内。




<a name="Resultsfora3dsphericalshelltestcase"></a><h3>Results for a 3d spherical shell testcase</h3>


三维计算在计算上是非常昂贵的。此外，如上所述，有趣的行为只有在相当长的时间后才开始，需要更多的CPU时间，而不是在一个典型的集群上可用。因此，与其在这里展示一个完整的模拟，不如让我们简单地展示几张图片，我们使用这个程序的后续程序，称为<i>ASPECT</i>（简称<i>Advanced
%Solver for Problems in Earth's ConvecTion</i>），该程序正在独立于deal.II开发，已经包括了下面讨论的一些扩展。下面两张图片显示了温度的等值线和领域（连同网格）在512个处理器上的划分。

<p align="center">  <img src="https://www.dealii.org/images/steps/developer/step-32.3d-sphere.solution.png" alt=""> 

 <img src="https://www.dealii.org/images/steps/developer/step-32.3d-sphere.partition.png" alt="">  </p>  。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序有许多可以扩展的方向。正如在介绍的最后提到的，在本教程程序完成时，其中大部分正在<i>ASPECT</i>（简称<i>Advanced %Solver for Problems
in Earth's ConvecTion</i>）代码中积极开发。具体来说，下面这些肯定是人们应该解决的话题，以使程序更加有用。

 <ul>   <li>  <b>Adiabatic heating/cooling:</b> 我们在模拟中得到的温度场在一段时间后大多是恒定的，在内部和外部边界有边界层，冷和热物质的流线混合一切。然而，这并不符合我们的预期，即靠近地心的东西应该比靠近地表的东西更热。原因是我们使用的能量方程不包括一个描述绝热冷却和加热的术语：岩石，像气体一样，在你压缩它的时候会加热。因此，上升的物质以绝热方式冷却，而下沉的冷物质则以绝热方式加热。因此，正确的温度方程看起来有点像这样。   @f{eqnarray*}
    \frac{D T}{Dt}


    -
    \nabla \cdot \kappa \nabla T &=& \gamma + \tau\frac{Dp}{Dt},
  @f}

  或者，扩大平流导数  $\frac{D}{Dt} =
  \frac{\partial}{\partial t} + \mathbf u \cdot \nabla$  : @f{eqnarray*}
    \frac{\partial T}{\partial t}
    +
    {\mathbf u} \cdot \nabla T


    -
    \nabla \cdot \kappa \nabla T &=& \gamma +
    \tau\left\{\frac{\partial
    p}{\partial t} + \mathbf u \cdot \nabla p \right\}.
  @f} 。

  换句话说，随着岩石体积中压力的增加（ $\frac{Dp}{Dt}>0$ ），我们会得到一个额外的热源，反之亦然。

  压力的时间导数实施起来有点困难。如果有必要，我们可以利用导言中概述的事实进行近似，即压力可以分解为由于温差和由此产生的流动而产生的动态部分，以及仅由上层岩石的静压力产生的静态部分。由于后者要大得多，我们可以对 $p\approx p_{\text{static}}=-\rho_{\text{ref}}
  [1+\beta T_{\text{ref}}] \varphi$ 进行近似处理，从而对 $\frac{Dp}{Dt} \approx \left\{- \mathbf u \cdot \nabla \rho_{\text{ref}}
  [1+\beta T_{\text{ref}}]\varphi\right\} = \rho_{\text{ref}}
  [1+\beta T_{\text{ref}}] \mathbf u \cdot \mathbf g$ 进行处理。   换句话说，如果流体沿着重力方向（向下）运动，它将被压缩，因为在这种情况下 $\mathbf u
  \cdot \mathbf g > 0$ 我们得到一个正的热源。反之，如果流体逆着重力方向运动，它将被冷却。

 <li>  <b>Compressibility:</b> 正如在上面的温度模型中已经暗示的那样，地幔岩石不是不可压缩的。相反，鉴于地幔中的巨大压力（在地核-地幔边界，压力约为140GPa，相当于大气压力的140万倍），岩石实际上确实被压缩到它在表面压力下的密度的1.5倍左右。对这一情况进行建模会遇到很多困难。首先，质量守恒方程不再是 $\textrm{div}\;\mathbf u=0$ ，而应该是 $\textrm{div}(\rho\mathbf u)=0$ ，其中密度 $\rho$ 现在不再是空间常数，而是取决于温度和压力。一个后果是，该模型现在不再是线性的；线性化的斯托克斯方程也不再是对称的，需要我们重新考虑预处理程序，甚至可能是离散化。至于如何解决这个问题，我们在这里就不做详细介绍了。

 <li>  <b>Nonlinear material models:</b> 正如在不同地方已经暗示的那样，材料参数，如密度、粘度和各种热参数，在整个地幔中并不恒定。相反，它们非线性地依赖于压力和温度，在粘度的情况下，还依赖于应变率  $\varepsilon(\mathbf u)$  。对于复杂的模型，准确解决这些模型的唯一方法可能是在每个时间步骤中实际迭代出这种依赖关系，而不是简单地将系数冻结在从前一个（几个）时间步骤推算出来的数值上。

 <li>  <b>Checkpoint/restart:</b> 在一些处理器上以2D运行这个程序可以在一两天内解决现实的模型。然而，在3d中，计算时间非常大，以至于会遇到两个典型问题。(i) 在大多数计算集群上，排队系统将单个作业的运行时间限制在2或3天；(ii) 在数百个处理器上运行几天，由于硬件故障、错误配置或断电而丢失计算结果是一种耻辱。这两个问题都可以通过定期保存程序的状态来解决，如果有必要，在这个时候重新启动程序。这种技术通常被称为<i>checkpoint/restart</i>，它要求将程序的整个状态写到一个永久的存储位置（例如硬盘）。考虑到这个程序的数据结构的复杂性，这并不是完全微不足道的（也可能涉及到写入数千兆字节或更多的数据），但可以通过意识到可以在两个时间步骤之间保存状态，其中基本上只包括网格和解向量；在重新启动期间，然后首先以之前的方式重新列举自由度，然后重新组装矩阵。然而，考虑到这里涉及的数据结构的分布性质，保存和恢复程序的状态并不简单。一个额外的复杂性是由以下事实引入的：人们可能希望在两次运行之间改变处理器的数量，例如，因为人们可能希望在一个比用于在中间时间预计算起始温度场的网格更精细的网格上继续计算。

 <li>  <b>Predictive postprocessing:</b> 像这样的计算的重点不是简单地解决方程。相反，它通常是探索不同的物理模型，并将其与我们在地球表面可以测量到的东西进行比较，以发现哪些模型是现实的，哪些是与现实相矛盾的。为此，我们需要从我们的解决方案向量中计算出与我们可以观察到的东西有关的数量。例如，其中包括地球表面的热流，以及整个地幔的地震速度，因为这些影响到地震仪所记录的地震波。

 <li>  <b>Better refinement criteria:</b> 从上面的3D案例可以看出，3D的网格主要是沿着内部边界细化的。这是因为那里的边界层比领域中的任何其他过渡都要强，导致我们几乎只在那里细化，基本上没有沿着羽流的方向细化。我们当然需要更好的细化标准来跟踪我们真正感兴趣的部分，而不是这里使用的标准，即应用于温度的KellyErrorEstimator，能够做到。   </ul> 


还有许多其他方法来扩展当前的程序。然而，与其在这里讨论它们，不如让我们指出更大的开放源代码ASPECT（见https://aspect.geodynamics.org/），它构成了step-32的进一步发展，并且已经包括了许多这样可能的扩展。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-32.cc"
*/
