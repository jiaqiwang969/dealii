/**
@page step_18 The step-18 tutorial program
This tutorial depends on step-17.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Quasistaticelasticdeformation">Quasistatic elastic deformation</a>
      <ul>
        <li><a href="#Motivationofthemodel">Motivation of the model</a>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Updatingthestressvariable">Updating the stress variable</a>
      </ul>
        <li><a href="#Parallelgraphicaloutput">Parallel graphical output</a>
        <li><a href="#Atriangulationwithautomaticpartitioning">A triangulation with automatic partitioning</a>
        <li><a href="#Overallstructureoftheprogram">Overall structure of the program</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ThecodePointHistorycodeclass">The <code>PointHistory</code> class</a>
        <li><a href="#Thestressstraintensor">The stress-strain tensor</a>
        <li><a href="#Auxiliaryfunctions">Auxiliary functions</a>
        <li><a href="#ThecodeTopLevelcodeclass">The <code>TopLevel</code> class</a>
        <li><a href="#ThecodeBodyForcecodeclass">The <code>BodyForce</code> class</a>
        <li><a href="#ThecodeIncrementalBoundaryValuecodeclass">The <code>IncrementalBoundaryValue</code> class</a>
        <li><a href="#ImplementationofthecodeTopLevelcodeclass">Implementation of the <code>TopLevel</code> class</a>
      <ul>
        <li><a href="#Thepublicinterface">The public interface</a>
        <li><a href="#TopLevelcreate_coarse_grid">TopLevel::create_coarse_grid</a>
        <li><a href="#TopLevelsetup_system">TopLevel::setup_system</a>
        <li><a href="#TopLevelassemble_system">TopLevel::assemble_system</a>
        <li><a href="#TopLevelsolve_timestep">TopLevel::solve_timestep</a>
        <li><a href="#TopLevelsolve_linear_problem">TopLevel::solve_linear_problem</a>
        <li><a href="#TopLeveloutput_results">TopLevel::output_results</a>
        <li><a href="#TopLeveldo_initial_timestep">TopLevel::do_initial_timestep</a>
        <li><a href="#TopLeveldo_timestep">TopLevel::do_timestep</a>
        <li><a href="#TopLevelrefine_initial_grid">TopLevel::refine_initial_grid</a>
        <li><a href="#TopLevelmove_mesh">TopLevel::move_mesh</a>
        <li><a href="#TopLevelsetup_quadrature_point_history">TopLevel::setup_quadrature_point_history</a>
        <li><a href="#TopLevelupdate_quadrature_point_history">TopLevel::update_quadrature_point_history</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
      <ul>
        <li><a href="#Plasticitymodels">Plasticity models</a>
        <li><a href="#Stabilizationissues">Stabilization issues</a>
        <li><a href="#Refinementduringtimesteps">Refinement during timesteps</a>
        <li><a href="#Ensuringmeshregularity">Ensuring mesh regularity</a>
    </ul>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-18/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这个教程程序是我们已经在步骤8和步骤17中开始的弹性问题系列中的另一个。它将其扩展到两个不同的方向：首先，它用拉格朗日网格运动方法解决了大变形的准静态但与时间相关的弹性问题。其次，它又展示了一些使用PETSc的线性代数的%并行处理来解决此类问题的技术。除此之外，我们还展示了如何解决step-17的两个主要瓶颈中的一个，即我们只从一个进程中产生图形输出，而这在更多的进程和大问题上的扩展性非常差。另一个瓶颈，即每个处理器都必须持有整个网格和DoFHandler，将在第40步中解决）。最后，我们还展示了许多以前的程序中未曾展示过的各种改进和技术。

如同前面的第17步，只要你安装了PETSc，程序在单机上的运行也是一样的。关于如何告诉deal.II你的系统上安装了PETSc的信息可以在deal.II的README文件中找到，该文件可以从你安装的deal.II的<a href="../../index.html">main
documentation page</a>中链接到，或者在<a href="http://www.dealii.org/">the
deal.II webpage</a>上。




<a name="Quasistaticelasticdeformation"></a><h3>Quasistatic elastic deformation</h3>


<a name="Motivationofthemodel"></a><h4>Motivation of the model</h4>


一般来说，随时间变化的小弹性变形是由弹性波方程描述的

@f[
  \rho \frac{\partial^2 \mathbf{u}}{\partial t^2}
  + c \frac{\partial \mathbf{u}}{\partial t}


  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) = \mathbf{f}
  \qquad
  \textrm{in}\ \Omega,


@f]

其中 $\mathbf{u}=\mathbf{u} (\mathbf{x},t)$ 是身体的变形， $\rho$ 和 $c$ 是密度和衰减系数，以及 $\mathbf{f}$ 外力。此外，初始条件

@f[
  \mathbf{u}(\cdot, 0) = \mathbf{u}_0(\cdot)
  \qquad
  \textrm{on}\ \Omega,


@f]

和Dirichlet（位移）或Neumann（牵引）边界条件，需要指定一个唯一的解决方案。

@f{eqnarray*}
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega\backslash\Gamma_D.


@f}

在上述公式中， $\varepsilon(\mathbf{u})= \frac 12 (\nabla \mathbf{u} + \nabla
\mathbf{u}^T)$ 是位移的对称梯度，也称为 <em> 应变 </em>  。   $C$ 是一个等级为4的张量，称为 <em> 应力-应变张量 </em> （<a
  href="https://en.wikipedia.org/wiki/Hooke%27s_law#Hooke's_law_for_continuous_media"><em>compliance
  tensor</em></a>的逆向）。]），它包含了材料弹性强度的知识；它的对称性特性确保它将秩为2的对称张量（&ldquo;矩阵&rdquo;的维数 $d$ ，其中 $d$ 是空间维数）映射到相同秩的对称张量上。我们将在下面更多地评论应变和应力张量的作用。现在只需要说，我们将术语 $\textrm{div}\  ( C \varepsilon(\mathbf{u}))$ 解释为具有分量 $\frac \partial{\partial x_j} C_{ijkl} \varepsilon(\mathbf{u})_{kl}$ 的向量，其中对指数 $j,k,l$ 的求和是隐含的。

这个方程的准静态极限的动机如下：身体的每个小扰动，例如边界条件或强迫函数的变化，将导致身体配置的相应变化。一般来说，这将是以波的形式从扰动的位置辐射出去。由于阻尼项的存在，这些波将在例如 $\tau$ 的时间尺度上被衰减。现在，假设所有外部强制力的变化发生在比 $\tau$ 大得多的时间尺度上。在这种情况下，变化的动态性质并不重要：我们可以认为身体总是处于静态平衡状态，也就是说，我们可以假设在任何时候，身体都满足于

@f{eqnarray*}


  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega,
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N.


@f}

请注意，微分方程不再包含任何时间导数 -- 所有的时间依赖性都是通过边界条件和可能的时间变化的力函数引入的  $\mathbf{f}(\mathbf{x},t)$  。因此，配置的变化可以被认为是瞬时静止的。对此的另一种看法是， $t$ 并不是真正的时间变量，而只是一个支配问题演变的类似时间的参数。

虽然这些方程足以描述小的变形，但计算大的变形就有点复杂了，一般来说，会导致非线性方程，如步骤-44中处理的那些。在下文中，让我们考虑在模拟变形成为<i>large</i>的问题时，人们会采用的一些工具。

 @note 我们下面要考虑的模型并不是建立在任何在数学上合理的基础上：我们将考虑一个模型，在这个模型中，我们产生一个小的变形，通过这个变形使身体的物理坐标变形，然后再考虑下一个加载步骤，作为一个线性问题。这并不一致，因为线性的假设意味着变形是无限小的，所以在解决下一个线性问题之前，在我们的网格顶点周围移动一个有限的量是不一致的做法。因此，我们应该注意到，在文献中找不到下面讨论的方程，这并不奇怪。<b>The model considered here has
little to do with reality!</b>另一方面，我们所考虑的实现技术正是人们在实现<i>real</i>模型时需要使用的，我们将在步骤-44中看到。


为了回到定义我们的 "人工 "模型，让我们首先引入一个张量的应力变量 $\sigma$ ，并以应力为基础写出微分方程。

@f{eqnarray*}


  - \textrm{div}\  \sigma &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega(t),
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega(t),
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega(t)\backslash\Gamma_D.


@f}

注意这些方程是在一个随时间变化的域 $\Omega(t)$ 上提出的，边界根据边界上各点的位移 $\mathbf{u}(\mathbf{x},t)$ 而移动。为了完成这个系统，我们必须指定应力和应变之间的增量关系，如下所示。<a name="step_18.stress-strain"></a>

@f[
  \dot\sigma = C \varepsilon (\dot{\mathbf{u}}),
  \qquad
  \qquad
  \textrm{[stress-strain]}


@f]

其中点表示一个时间导数。应力 $\sigma$ 和应变 $\varepsilon(\mathbf{u})$ 都是等级2的对称张量。




<a name="Timediscretization"></a><h4>Time discretization</h4>


在数值上，该系统的求解方法如下：首先，我们使用后向欧拉方案对时间部分进行离散化。这导致了时间步长的离散平衡力  $n$  。

@f[


  -\textrm{div}\  \sigma^n = f^n,


@f]

其中

@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n),


@f]

和  $\Delta \mathbf{u}^n$  时间步长的增量位移  $n$  。此外，我们必须指定初始数据  $\mathbf{u}(\cdot,0)=\mathbf{u}_0$  。这样一来，如果我们想求解位移增量，我们必须求解以下系统。

@f{align*}


  - \textrm{div}\   C \varepsilon(\Delta\mathbf{u}^n) &= \mathbf{f} + \textrm{div}\  \sigma^{n-1}
  \qquad
  &&\textrm{in}\ \Omega(t_{n-1}),
  \\
  \Delta \mathbf{u}^n(\mathbf{x},t) &= \mathbf{d}(\mathbf{x},t_n) - \mathbf{d}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_D\subset\partial\Omega(t_{n-1}),
  \\
  \mathbf{n} \ C \varepsilon(\Delta \mathbf{u}^n(\mathbf{x},t)) &= \mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_N=\partial\Omega(t_{n-1})\backslash\Gamma_D.


@f}

这组方程的弱形式，像往常一样是有限元公式的基础，其内容如下：找到 $\Delta \mathbf{u}^n \in
\{v\in H^1(\Omega(t_{n-1}))^d: v|_{\Gamma_D}=\mathbf{d}(\cdot,t_n) - \mathbf{d}(\cdot,t_{n-1})\}$ ，使<a name="step_18.linear-system"></a>这样的方程。

@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}


  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  \\
  &\qquad
  +(\mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1}), \varphi)_{\Gamma_N}
  +(\sigma^{n-1} \mathbf{n}, \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.


@f}

利用 $\sigma^{n-1} \mathbf{n}
            = [C \varepsilon(\mathbf{u}^{n-1})] \mathbf{n}
            = \mathbf{b}(\mathbf x, t_{n-1})$ ，这些方程可以简化为

@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}


  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  +(\mathbf{b}(\mathbf{x},t_n),t_{n-1}), \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.
  \qquad
  \qquad
  \textrm{[linear-system]}


@f}



我们注意到，为了简单起见，在程序中我们总是假设没有边界力，即 $\mathbf{b} = 0$ ，并且身体的变形仅由身体力 $\mathbf{f}$ 和规定的边界位移 $\mathbf{d}$ 驱动。还值得注意的是，当通过部分积分时，我们会得到形式为 $(C \varepsilon(\Delta\mathbf{u}^n), \nabla \varphi
)_{\Omega(t_{n-1})}$ 的条款，但我们用涉及对称梯度的条款 $\varepsilon(\varphi)$ 而不是 $\nabla\varphi$ 来取代它们。由于 $C$ 的对称性，这两个项在数学上是等价的，但对称版本避免了可能出现的四舍五入错误，使得到的矩阵略显非对称性。

在时间步长 $n$ 的系统，要在旧域 $\Omega(t_{n-1})$ 上求解，其形式完全是一个静止的弹性问题，因此与我们在以前的例子程序中已经实现的类似。因此，除了说我们再次使用最低阶连续有限元之外，我们将不对空间离散化进行评论。

但也有不同之处。<ol>  <li>  我们必须在每个时间步骤之后移动（更新）网格，以便能够在新的领域上解决下一个时间步骤。

    <li> 我们需要知道 $\sigma^{n-1}$ 来计算下一个增量位移，也就是说，我们需要在时间步骤结束时计算它，以确保它可以用于下一个时间步骤。从本质上讲，应力变量是我们了解体的变形历史的窗口。   </ol>  这两个操作在程序中的 <code>move_mesh</code> 和 <code>update_quadrature_point_history</code> 函数中完成。移动网格只是一个技术问题，而更新应力则要复杂一些，将在下一节讨论。




<a name="Updatingthestressvariable"></a><h4>Updating the stress variable</h4>


如上所述，在计算时间步长 $n+1$ 时，我们需要有应力变量 $\sigma^n$ ，我们可以用<a name="step_18.stress-update"></a>来计算它。

@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update]}


@f]

尽管这个方程表面上很简单，但有两个问题我们需要讨论。第一个问题是关于我们存储 $\sigma^n$ 的方式：即使我们使用最低阶有限元计算增量更新 $\Delta\mathbf{u}^n$ ，那么其对称梯度 $\varepsilon(\Delta\mathbf{u}^n)$ 一般来说仍然是一个不容易描述的函数。特别是，它不是一个片状常数函数，在一般的网格上（单元不是平行于坐标轴的矩形）或非恒定应力-应变张量 $C$ ，它甚至不是一个双线性或三线性函数。因此，如何在计算机程序中存储 $\sigma^n$ 是先验的。

要决定这一点，我们必须看它被用在什么地方。我们需要应力的唯一地方是在术语 $(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}$ 中。在实践中，我们当然会用数值正交来代替这个项。

@f[
  (\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  =
  \sum_{K\subset {T}}
  (\sigma^{n-1},\varepsilon(\varphi))_K
  \approx
  \sum_{K\subset {T}}
  \sum_q
  w_q \ \sigma^{n-1}(\mathbf{x}_q) : \varepsilon(\varphi(\mathbf{x}_q),


@f]

其中 $w_q$ 是正交权重， $\mathbf{x}_q$ 是单元格 $K$ 上的正交点。这应该表明，我们真正需要的不是应力 $\sigma^{n-1}$ 本身，而只是所有单元上的正交点的应力值。然而，这是一个更简单的任务：我们只需要提供一个数据结构，能够为所有单元上的每个正交点（或者，由于我们是并行计算，目前的MPI进程拥有的所有单元的所有正交点&ldquo;rdquo;）容纳一个等级为2的对称张量。在每个时间步骤结束时，我们只需评估 $\varepsilon(\Delta \mathbf{u}^n(\mathbf{x}_q))$ ，将其乘以应力-应变张量 $C$ ，并使用该结果来更新正交点 $q$ 的应力 $\sigma^n(\mathbf{x}_q)$ 。

第二个复杂的问题在我们上面选择的符号中并不明显。这是由于我们在域 $\Omega(t_{n-1})$ 上计算 $\Delta u^n$ ，然后用这个位移增量来更新应力，同时移动网格节点，以达到 $\Omega(t_n)$ ，在此基础上计算下一个增量。在这种情况下，我们必须确定的是，移动网格不仅涉及到节点的移动，还涉及到应力变量的相应变化：更新的应力是一个相对于旧域中材料的坐标系而定义的变量，必须转移到新域中。其原因可以理解为：在局部，增量变形 $\Delta\mathbf{u}$ 可以分解为三个部分，线性平移（点附近的位移增量场的常数部分），扩张分量（位移场梯度中具有非零发散的那部分），以及旋转。材料的线性平移并不影响冻结在其中的应力--应力值只是沿着平移。扩张或压缩的变化产生相应的应力更新。然而，旋转分量不一定会引起非零的应力更新（想想，在2d中，例如 $\Delta\mathbf{u}=(y, -x)^T$  ，与 $\varepsilon(\Delta
\mathbf{u})=0$  的情况）。尽管如此，如果材料在某个方向上被预应力，那么这个方向将随着材料的旋转而旋转。  为此，我们必须定义一个旋转矩阵 $R(\Delta \mathbf{u}^n)$ ，描述在每一个点上由于位移增量而产生的旋转。不难看出， $R$ 对 $\Delta \mathbf{u}^n$ 的实际依赖只能是通过位移的卷曲，而不是位移本身或其全部梯度（如上所述，增量的常数分量描述平移，其发散描述扩张模式，而卷曲描述旋转模式）。由于 $R$ 的确切形式很麻烦，我们只在程序代码中说明，并注意到应力变量的正确更新公式是<a name="step_18.stress-update+rot"></a> 。

@f[
  \sigma^n
  =
  R(\Delta \mathbf{u}^n)^T
  [\sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n)]
  R(\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update+rot]}


@f]



应力更新和旋转都是在示例程序的函数 <code>update_quadrature_point_history</code> 中实现的。




<a name="Parallelgraphicaloutput"></a><h3>Parallel graphical output</h3>


在步骤17中，就运行时间而言，平行计算的主要瓶颈是只有第一个处理器产生整个领域的输出。由于生成图形输出是很昂贵的，所以当涉及到更多数量的处理器时，这并不能很好地扩展。我们将在这里解决这个问题。关于程序 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"）。

基本上，我们需要做的是让每个进程为它所拥有的那个单元子集产生图形输出，将它们写进单独的文件，并有办法同时显示某个时间步长的所有文件。这样，代码在每个时间步长的每个进程产生一个 <code>.vtu</code> 文件。两个常见的VTK文件查看器ParaView和Viscit都支持一次打开一个以上的 <code>.vtu</code> 文件。为了简化挑选正确文件的过程，并允许在时间上移动，两者都支持记录文件，以引用特定时间步长的所有文件。遗憾的是，记录文件在VisIt和Paraview之间有不同的格式，所以我们把两种格式都写出来。

代码将生成文件 <code>solution-TTTT.NNN.vtu</code> ，其中 <code>TTTT</code> 是时间步数（从1开始）， <code>NNN</code> 是进程等级（从0开始）。这些文件包含时间段和处理器的本地所有单元。文件 <code>solution-TTTT.visit</code> 是时间段的访问记录 <code>TTTT</code>, while <code>solution-TTTT.pvtu</code> 对ParaView也是如此。(较新版本的VisIt实际上也可以读取 <code>.pvtu</code> 文件，但输出两种记录文件也无妨。)最后， <code>solution.pvd</code> 文件是只有ParaView支持的特殊记录，它引用所有的时间步骤。所以在ParaView中，只需要打开solution.pvd，而在VisIt中需要选择所有的.visit文件组，才能达到同样的效果。




<a name="Atriangulationwithautomaticpartitioning"></a><h3>A triangulation with automatic partitioning</h3>


在第17步中，我们使用了一个在每个处理器上简单复制的常规三角形，以及一个相应的DoFHandler。两者都不知道它们是在%并行环境下使用的--它们只是完整地存在于每个处理器上，我们认为这最终会成为一个主要的内存瓶颈。

我们在这里不解决这个问题（我们将在第40步中解决），但使情况稍微自动化一些。在第17步中，我们创建了三角形，然后手动 "分区"，也就是说，我们给每个单元分配了 @ref GlossSubdomainId "子域ID"，以表明哪个 @ref GlossMPIProcess "MPI进程""拥有 "该单元。在这里，我们使用了一个类 parallel::shared::Triangulation ，它至少自动完成了这一部分：每当你创建或完善这样一个三角图时，它都会自动在所有参与的进程之间进行划分（它知道这些进程，因为你必须告诉它在构建三角图时连接这些进程的 @ref GlossMPICommunicator "MPI通信器"）。否则， parallel::shared::Triangulation 看起来，就所有的实际目的而言，就像一个普通的Triangulation对象。

使用这个类的便利性不仅来自于能够避免手动调用 GridTools::partition(). ，相反，DoFHandler类现在也知道你想在并行环境下使用它，并且默认情况下会自动列举自由度，使进程0拥有的所有DoF先于进程1拥有的所有DoF，等等。换句话说，你也可以避免对 DoFRenumbering::subdomain_wise(). 的调用。

还有其他好处。例如，由于三角计算知道它生活在一个%parallel universe中，它也知道它 "拥有 "某些单元（即那些子域id等于其MPI等级的单元；以前，三角计算只存储这些子域id，但没有办法使它们有意义）。因此，在汇编函数中，你可以测试一个单元是否 "本地拥有"（即由当前进程拥有，见 @ref GlossLocallyOwnedCell ），当你在所有单元上循环时，使用以下语法

@code
  if (cell->is_locally_owned())
@endcode

这种知识延伸到建立在这种三角形上的DoFHandler对象，然后它可以通过 DoFHandler::compute_n_locally_owned_dofs_per_processor() 和 DoFTools::extract_locally_relevant_dofs(). 等调用来识别哪些自由度是本地拥有的（见 @ref GlossLocallyOwnedDof ）。最后，DataOut类也知道如何处理这种三角形，并将简单地跳过在非本地拥有的单元上生成图形输出。

当然，正如在第17步的讨论中多次指出的那样，在每个进程上保持整个三角形将无法扩展：大型问题可能根本无法再适合每个进程的内存，即使我们有足够多的进程在合理的时间内解决它们。在这种情况下， parallel::shared::Triangulation 不再是一个合理的计算基础，我们将在步骤40中展示如何使用 parallel::distributed::Triangulation 类来解决这个问题，即让每个进程只存储一个<i>part</i>的三角图。




<a name="Overallstructureoftheprogram"></a><h3>Overall structure of the program</h3>


程序的整体结构可以从 <code>run()</code> 函数中推断出来，该函数首先在第一个时间步骤中调用 <code>do_initial_timestep()</code> ，然后在所有后续时间步骤中调用 <code>do_timestep()</code> 。这些函数之间的区别仅仅在于，在第一个时间步骤中，我们从一个粗略的网格开始，在其上求解，自适应地细化网格，然后在新的网格上以干净的状态重新开始。这个过程给了我们一个更好的起始网格，尽管我们当然应该在迭代过程中不断调整网格--这个程序中没有这样做，但是下面会有评论。

这两个处理时间步骤的函数的共同部分是在本网格上的以下操作序列。   <ul>   <li>   <code>assemble_system ()</code> [via <code>solve_timestep ()</code>  ] 。   这第一个函数也是最有趣的一个。它组装了对应于方程<a href="#step_18.linear-system">[linear-system]</a>离散化版本的线性系统。这导致了一个系统矩阵 $A_{ij} = \sum_K
  A^K_{ij}$ ，由每个单元 $K$ 上的局部贡献组成，其条目为@f[
    A^K_{ij} = (C \varepsilon(\varphi_j), \varepsilon(\varphi_i))_K;
  @f] 。

  在实践中， $A^K$ 是根据公式@f[
    A^K_{ij} = \sum_q w_q [\varepsilon(\varphi_i(\mathbf{x}_q)) : C :
                           \varepsilon(\varphi_j(\mathbf{x}_q))],
  @f]使用数值正交计算出来的。

  与正交点 $\mathbf{x}_q$ 和权重 $w_q$  。我们之前在步骤8和步骤17中建立了这些贡献，但在这两种情况下，我们都是通过使用等级4张量 $C$ 的组成知识，以及考虑应变张量的单个元素 $\varepsilon(\varphi_i),\varepsilon(\varphi_j)$ ，相当笨拙地完成的。这其实并不方便，特别是如果我们想考虑比各向同性的情况更复杂的弹性模型，而 $C$ 有方便的形式  $C_{ijkl}  = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik} \delta_{jl}
  + \delta_{il} \delta_{jk})$  。虽然我们在本程序中没有使用比这更复杂的形式，但我们还是希望以一种容易实现的方式来编写它。因此，很自然地要引入代表等级为2（用于应变和应力）和4（用于应力-应变张量 $C$ ）的对称张量的类。幸运的是，deal.II提供了这些： <code>SymmetricTensor<rank,dim></code> 类模板提供了等级 <code>rank</code> （需要是偶数）和维度 <code>dim</code> 的这类张量的完整实现。

  然后我们需要的是两件事：一种创建应力-应变等级4张量 $C$ 的方法，以及从形状函数 $\varphi_i$ 的梯度在给定单元上的正交点 $\mathbf{x}_q$ 创建一个等级2的对称张量（应变张量）。在这个例子程序的执行顶部，你会发现这样的函数。第一个， <code>get_stress_strain_tensor</code>  ，需要两个参数，对应于Lam&eacute; 常数 $\lambda$ 和 $\mu$ ，并返回对应于这些常数的各向同性的应力应变张量（在程序中，我们将选择对应于钢的常数）；用一个计算各向异性的张量的函数来代替这个函数是很简单的，或者考虑到晶体对称性，比如说。第二个， <code>get_strain</code> takes an object of type <code>FEValues</code> 和指数 $i$ 和 $q$ ，返回对称梯度，即应变，对应于形状函数 $\varphi_i(\mathbf{x}_q)$ ，在 <code>FEValues</code> 对象最后被重新初始化的单元上评估。

  鉴于此， <code>assemble_system</code> 的最内部循环以下列优雅的方式计算对矩阵的局部贡献（变量 <code>stress_strain_tensor</code> ，对应于张量 $C$ ，之前已经用上述第一个函数的结果初始化了）。   @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  for (unsigned int j=0; j<dofs_per_cell; ++j)
    for (unsigned int q_point=0; q_point<n_q_points;
         ++q_point)
      {
        const SymmetricTensor<2,dim>
          eps_phi_i = get_strain (fe_values, i, q_point),
          eps_phi_j = get_strain (fe_values, j, q_point);


        cell_matrix(i,j)
          += (eps_phi_i * stress_strain_tensor * eps_phi_j *
              fe_values.JxW (q_point));
      }
  @endcode

  值得注意的是这段代码的表现力，并将其与我们在以前的例子中为弹性问题所经历的复杂情况进行比较。公平地说，在写这些以前的例子时，SymmetricTensor类模板还不存在）。为了简单起见， <code>operator*</code> 在这里规定了偶数等级的对称张量之间的（双重求和）积。

  组建本地捐款@f{eqnarray*}
      f^K_i &=&
      (\mathbf{f}, \varphi_i)_K -(\sigma^{n-1},\varepsilon(\varphi_i))_K
      \\
      &\approx&
      \sum_q
      w_q \left\{
        \mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q) -
        \sigma^{n-1}_q : \varepsilon(\varphi_i(\mathbf{x}_q))
      \right\}
  @f}

  到<a href="#step_18.linear-system">[linear-system]</a>的右手边同样是直接的（注意，我们在这里不考虑任何边界牵引 $\mathbf{b}$ ）。请记住，我们只需要在单元格的正交点上存储旧的应力。在程序中，我们将提供一个变量 <code>local_quadrature_points_data</code> ，允许访问每个正交点的应力 $\sigma^{n-1}_q$ 。有了这个，右手边的代码看起来就像这样，同样相当优雅。   @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  {
    const unsigned int
      component_i = fe.system_to_component_index(i).first;


    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        const SymmetricTensor<2,dim> &old_stress
          = local_quadrature_points_data[q_point].old_stress;


        cell_rhs(i) += (body_force_values[q_point](component_i) *
                        fe_values.shape_value (i,q_point)


                        -
                        old_stress *
                        get_strain (fe_values,i,q_point)) *
                       fe_values.JxW (q_point);
      }
  }
  @endcode

  请注意，在乘法 $\mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q)$ 中，我们利用了这样一个事实：对于所选的有限元素， $\varphi_i$ 中只有一个向量分量（即 <code>component_i</code> ）是非零的，因此我们也只需要考虑 $\mathbf{f}(\mathbf{x}_q)$ 的一个分量。

  这基本上结束了我们在这个函数中提出的新材料。它后来必须处理边界条件以及悬挂节点约束，但这与我们以前在其他程序中已经要做的事情相类似。

 <li>   <code>solve_linear_problem ()</code> [via <code>solve_timestep ()</code>  ] 。   与前一个函数不同，这个函数其实并不有趣，因为它做的是以前所有教程程序中的类似函数--用CG方法求解线性系统，使用不完整的LU分解作为预处理程序（在%并行情况下，它分别使用每个处理器块的ILU）。它与第17步几乎没有变化。

 <li>   <code>update_quadrature_point_history ()</code>  [通过 <code>solve_timestep ()</code>  ] 。基于之前计算的位移场 $\Delta \mathbf{u}^n$ ，我们根据<a href="#step_18.stress-update">[stress-update]</a>和<a href="#step_18.stress-update+rot">[stress-update+rot]</a>更新所有正交点的应力值，包括坐标系的旋转。

 <li>   <code>move_mesh ()</code>  ：给定之前计算的解决方案，在这个函数中，我们通过移动每个顶点的位移矢量场来实现网格的变形。

 <li>   <code>output_results ()</code>  : 这个函数只是根据我们上面所说的输出解决方案，也就是说，每个处理器只对自己的那部分域计算输出。除了解决方案，我们还计算了每个单元上所有正交点平均的应力的规范。   </ul> 

有了这个代码的一般结构，我们只需要定义我们要解决的情况。在本程序中，我们选择模拟一个垂直圆柱体的准静态变形，其底部边界是固定的，顶部边界以规定的垂直速度被推倒。然而，顶层边界的水平速度没有被指定--我们可以把这种情况想象成一块油性良好的板从顶部推到圆柱体上，圆柱体顶层边界上的点被允许沿着板的表面水平滑动，但被板强迫向下移动。圆柱体的内部和外部边界是自由的，不受任何规定的偏转或牵引的影响。此外，重力作用于身体。

程序文本将揭示更多关于如何实现这种情况，而结果部分将显示这种模拟产生的位移模式。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * First the usual list of header files that have already been used in
 * previous example programs:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/multithread_info.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/petsc_vector.h>
 * #include <deal.II/lac/petsc_sparse_matrix.h>
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/petsc_precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparsity_tools.h>
 * #include <deal.II/distributed/shared_tria.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * And here the only three new things among the header files: an include file in
 * which symmetric tensors of rank 2 and 4 are implemented, as introduced in
 * the introduction:
 * 
 * @code
 * #include <deal.II/base/symmetric_tensor.h>
 * 
 * @endcode
 * 
 * And lastly a header that contains some functions that will help us compute
 * rotaton matrices of the local coordinate systems at specific points in the
 * domain.
 * 
 * @code
 * #include <deal.II/physics/transformations.h>
 * 
 * @endcode
 * 
 * This is then simply C++ again:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * #include <iomanip>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step18
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePointHistorycodeclass"></a> 
 * <h3>The <code>PointHistory</code> class</h3>
 * 

 * 
 * As was mentioned in the introduction, we have to store the old stress in
 * quadrature point so that we can compute the residual forces at this point
 * during the next time step. This alone would not warrant a structure with
 * only one member, but in more complicated applications, we would have to
 * store more information in quadrature points as well, such as the history
 * variables of plasticity, etc. In essence, we have to store everything
 * that affects the present state of the material here, which in plasticity
 * is determined by the deformation history variables.
 *   

 * 
 * We will not give this class any meaningful functionality beyond being
 * able to store data, i.e. there are no constructors, destructors, or other
 * member functions. In such cases of `dumb' classes, we usually opt to
 * declare them as <code>struct</code> rather than <code>class</code>, to
 * indicate that they are closer to C-style structures than C++-style
 * classes.
 * 
 * @code
 *   template <int dim>
 *   struct PointHistory
 *   {
 *     SymmetricTensor<2, dim> old_stress;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thestressstraintensor"></a> 
 * <h3>The stress-strain tensor</h3>
 * 

 * 
 * Next, we define the linear relationship between the stress and the strain
 * in elasticity. It is given by a tensor of rank 4 that is usually written
 * in the form $C_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il}
 * \delta_{jk}) + \lambda \delta_{ij} \delta_{kl}$. This tensor maps
 * symmetric tensor of rank 2 to symmetric tensors of rank 2. A function
 * implementing its creation for given values of the Lam&eacute; constants
 * $\lambda$ and $\mu$ is straightforward:
 * 
 * @code
 *   template <int dim>
 *   SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
 *                                                    const double mu)
 *   {
 *     SymmetricTensor<4, dim> tmp;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = 0; j < dim; ++j)
 *         for (unsigned int k = 0; k < dim; ++k)
 *           for (unsigned int l = 0; l < dim; ++l)
 *             tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
 *                                ((i == l) && (j == k) ? mu : 0.0) +
 *                                ((i == j) && (k == l) ? lambda : 0.0));
 *     return tmp;
 *   }
 * 
 * @endcode
 * 
 * With this function, we will define a static member variable of the main
 * class below that will be used throughout the program as the stress-strain
 * tensor. Note that in more elaborate programs, this will probably be a
 * member variable of some class instead, or a function that returns the
 * stress-strain relationship depending on other input. For example in
 * damage theory models, the Lam&eacute; constants are considered a function
 * of the prior stress/strain history of a point. Conversely, in plasticity
 * the form of the stress-strain tensor is modified if the material has
 * reached the yield stress in a certain point, and possibly also depending on
 * its prior history.
 *   

 * 
 * In the present program, however, we assume that the material is
 * completely elastic and linear, and a constant stress-strain tensor is
 * sufficient for our present purposes.
 * 

 * 
 * 

 * 
 * 

 * 
 * 
 * <a name="Auxiliaryfunctions"></a> 
 * <h3>Auxiliary functions</h3>
 * 

 * 
 * Before the rest of the program, here are a few functions that we need as
 * tools. These are small functions that are called in inner loops, so we
 * mark them as <code>inline</code>.
 *   

 * 
 * The first one computes the symmetric strain tensor for shape function
 * <code>shape_func</code> at quadrature point <code>q_point</code> by
 * forming the symmetric gradient of this shape function. We need that when
 * we want to form the matrix, for example.
 *   

 * 
 * We should note that in previous examples where we have treated
 * vector-valued problems, we have always asked the finite element object in
 * which of the vector component the shape function is actually non-zero,
 * and thereby avoided to compute any terms that we could prove were zero
 * anyway. For this, we used the <code>fe.system_to_component_index</code>
 * function that returns in which component a shape function was zero, and
 * also that the <code>fe_values.shape_value</code> and
 * <code>fe_values.shape_grad</code> functions only returned the value and
 * gradient of the single non-zero component of a shape function if this is
 * a vector-valued element.
 *   

 * 
 * This was an optimization, and if it isn't terribly time critical, we can
 * get away with a simpler technique: just ask the <code>fe_values</code>
 * for the value or gradient of a given component of a given shape function
 * at a given quadrature point. This is what the
 * <code>fe_values.shape_grad_component(shape_func,q_point,i)</code> call
 * does: return the full gradient of the <code>i</code>th component of shape
 * function <code>shape_func</code> at quadrature point
 * <code>q_point</code>. If a certain component of a certain shape function
 * is always zero, then this will simply always return zero.
 *   

 * 
 * As mentioned, using <code>fe_values.shape_grad_component</code> instead
 * of the combination of <code>fe.system_to_component_index</code> and
 * <code>fe_values.shape_grad</code> may be less efficient, but its
 * implementation is optimized for such cases and shouldn't be a big
 * slowdown. We demonstrate the technique here since it is so much simpler
 * and straightforward.
 * 
 * @code
 *   template <int dim>
 *   inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values,
 *                                             const unsigned int   shape_func,
 *                                             const unsigned int   q_point)
 *   {
 * @endcode
 * 
 * Declare a temporary that will hold the return value:
 * 
 * @code
 *     SymmetricTensor<2, dim> tmp;
 * 
 * @endcode
 * 
 * First, fill diagonal terms which are simply the derivatives in
 * direction <code>i</code> of the <code>i</code> component of the
 * vector-valued shape function:
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i)
 *       tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
 * 
 * @endcode
 * 
 * Then fill the rest of the strain tensor. Note that since the tensor is
 * symmetric, we only have to compute one half (here: the upper right
 * corner) of the off-diagonal elements, and the implementation of the
 * <code>SymmetricTensor</code> class makes sure that at least to the
 * outside the symmetric entries are also filled (in practice, the class
 * of course stores only one copy). Here, we have picked the upper right
 * half of the tensor, but the lower left one would have been just as
 * good:
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = i + 1; j < dim; ++j)
 *         tmp[i][j] =
 *           (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
 *            fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
 *           2;
 * 
 *     return tmp;
 *   }
 * 
 * 
 * @endcode
 * 
 * The second function does something very similar (and therefore is given
 * the same name): compute the symmetric strain tensor from the gradient of
 * a vector-valued field. If you already have a solution field, the
 * <code>fe_values.get_function_gradients</code> function allows you to
 * extract the gradients of each component of your solution field at a
 * quadrature point. It returns this as a vector of rank-1 tensors: one rank-1
 * tensor (gradient) per vector component of the solution. From this we have
 * to reconstruct the (symmetric) strain tensor by transforming the data
 * storage format and symmetrization. We do this in the same way as above,
 * i.e. we avoid a few computations by filling first the diagonal and then
 * only one half of the symmetric tensor (the <code>SymmetricTensor</code>
 * class makes sure that it is sufficient to write only one of the two
 * symmetric components).
 *   

 * 
 * Before we do this, though, we make sure that the input has the kind of
 * structure we expect: that is that there are <code>dim</code> vector
 * components, i.e. one displacement component for each coordinate
 * direction. We test this with the <code>Assert</code> macro that will
 * simply abort our program if the condition is not met.
 * 
 * @code
 *   template <int dim>
 *   inline SymmetricTensor<2, dim>
 *   get_strain(const std::vector<Tensor<1, dim>> &grad)
 *   {
 *     Assert(grad.size() == dim, ExcInternalError());
 * 
 *     SymmetricTensor<2, dim> strain;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       strain[i][i] = grad[i][i];
 * 
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = i + 1; j < dim; ++j)
 *         strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
 * 
 *     return strain;
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally, below we will need a function that computes the rotation matrix
 * induced by a displacement at a given point. In fact, of course, the
 * displacement at a single point only has a direction and a magnitude, it
 * is the change in direction and magnitude that induces rotations. In
 * effect, the rotation matrix can be computed from the gradients of a
 * displacement, or, more specifically, from the curl.
 *   

 * 
 * The formulas by which the rotation matrices are determined are a little
 * awkward, especially in 3d. For 2d, there is a simpler way, so we
 * implement this function twice, once for 2d and once for 3d, so that we
 * can compile and use the program in both space dimensions if so desired --
 * after all, deal.II is all about dimension independent programming and
 * reuse of algorithm thoroughly tested with cheap computations in 2d, for
 * the more expensive computations in 3d. Here is one case, where we have to
 * implement different algorithms for 2d and 3d, but then can write the rest
 * of the program in a way that is independent of the space dimension.
 *   

 * 
 * So, without further ado to the 2d implementation:
 * 
 * @code
 *   Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>> &grad_u)
 *   {
 * @endcode
 * 
 * First, compute the curl of the velocity field from the gradients. Note
 * that we are in 2d, so the rotation is a scalar:
 * 
 * @code
 *     const double curl = (grad_u[1][0] - grad_u[0][1]);
 * 
 * @endcode
 * 
 * From this, compute the angle of rotation:
 * 
 * @code
 *     const double angle = std::atan(curl);
 * 
 * @endcode
 * 
 * And from this, build the antisymmetric rotation matrix. We want this
 * rotation matrix to represent the rotation of the local coordinate system
 * with respect to the global Cartesian basis, to we construct it with a
 * negative angle. The rotation matrix therefore represents the rotation
 * required to move from the local to the global coordinate system.
 * 
 * @code
 *     return Physics::Transformations::Rotations::rotation_matrix_2d(-angle);
 *   }
 * 
 * 
 * @endcode
 * 
 * The 3d case is a little more contrived:
 * 
 * @code
 *   Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u)
 *   {
 * @endcode
 * 
 * Again first compute the curl of the velocity field. This time, it is a
 * real vector:
 * 
 * @code
 *     const Point<3> curl(grad_u[2][1] - grad_u[1][2],
 *                         grad_u[0][2] - grad_u[2][0],
 *                         grad_u[1][0] - grad_u[0][1]);
 * 
 * @endcode
 * 
 * From this vector, using its magnitude, compute the tangent of the angle
 * of rotation, and from it the actual angle of rotation with respect to
 * the Cartesian basis:
 * 
 * @code
 *     const double tan_angle = std::sqrt(curl * curl);
 *     const double angle     = std::atan(tan_angle);
 * 
 * @endcode
 * 
 * Now, here's one problem: if the angle of rotation is too small, that
 * means that there is no rotation going on (for example a translational
 * motion). In that case, the rotation matrix is the identity matrix.
 *     

 * 
 * The reason why we stress that is that in this case we have that
 * <code>tan_angle==0</code>. Further down, we need to divide by that
 * number in the computation of the axis of rotation, and we would get
 * into trouble when dividing doing so. Therefore, let's shortcut this and
 * simply return the identity matrix if the angle of rotation is really
 * small:
 * 
 * @code
 *     if (std::abs(angle) < 1e-9)
 *       {
 *         static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
 *         static const Tensor<2, 3> rot(rotation);
 *         return rot;
 *       }
 * 
 * @endcode
 * 
 * Otherwise compute the real rotation matrix. For this, again we rely on
 * a predefined function to compute the rotation matrix of the local
 * coordinate system.
 * 
 * @code
 *     const Point<3> axis = curl / tan_angle;
 *     return Physics::Transformations::Rotations::rotation_matrix_3d(axis,
 *                                                                    -angle);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeTopLevelcodeclass"></a> 
 * <h3>The <code>TopLevel</code> class</h3>
 * 

 * 
 * This is the main class of the program. Since the namespace already
 * indicates what problem we are solving, let's call it by what it does: it
 * directs the flow of the program, i.e. it is the toplevel driver.
 *   

 * 
 * The member variables of this class are essentially as before, i.e. it has
 * to have a triangulation, a DoF handler and associated objects such as
 * constraints, variables that describe the linear system, etc. There are a
 * good number of more member functions now, which we will explain below.
 *   

 * 
 * The external interface of the class, however, is unchanged: it has a
 * public constructor and destructor, and it has a <code>run</code>
 * function that initiated all the work.
 * 
 * @code
 *   template <int dim>
 *   class TopLevel
 *   {
 *   public:
 *     TopLevel();
 *     ~TopLevel();
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * The private interface is more extensive than in step-17. First, we
 * obviously need functions that create the initial mesh, set up the
 * variables that describe the linear system on the present mesh
 * (i.e. matrices and vectors), and then functions that actually assemble
 * the system, direct what has to be solved in each time step, a function
 * that solves the linear system that arises in each timestep (and returns
 * the number of iterations it took), and finally output the solution
 * vector on the correct mesh:
 * 
 * @code
 *     void create_coarse_grid();
 * 
 *     void setup_system();
 * 
 *     void assemble_system();
 * 
 *     void solve_timestep();
 * 
 *     unsigned int solve_linear_problem();
 * 
 *     void output_results() const;
 * 
 * @endcode
 * 
 * All, except for the first two, of these functions are called in each
 * timestep. Since the first time step is a little special, we have
 * separate functions that describe what has to happen in a timestep: one
 * for the first, and one for all following timesteps:
 * 
 * @code
 *     void do_initial_timestep();
 * 
 *     void do_timestep();
 * 
 * @endcode
 * 
 * Then we need a whole bunch of functions that do various things. The
 * first one refines the initial grid: we start on the coarse grid with a
 * pristine state, solve the problem, then look at it and refine the mesh
 * accordingly, and start the same process over again, again with a
 * pristine state. Thus, refining the initial mesh is somewhat simpler
 * than refining a grid between two successive time steps, since it does
 * not involve transferring data from the old to the new triangulation, in
 * particular the history data that is stored in each quadrature point.
 * 
 * @code
 *     void refine_initial_grid();
 * 
 * @endcode
 * 
 * At the end of each time step, we want to move the mesh vertices around
 * according to the incremental displacement computed in this time
 * step. This is the function in which this is done:
 * 
 * @code
 *     void move_mesh();
 * 
 * @endcode
 * 
 * Next are two functions that handle the history variables stored in each
 * quadrature point. The first one is called before the first timestep to
 * set up a pristine state for the history variables. It only works on
 * those quadrature points on cells that belong to the present processor:
 * 
 * @code
 *     void setup_quadrature_point_history();
 * 
 * @endcode
 * 
 * The second one updates the history variables at the end of each
 * timestep:
 * 
 * @code
 *     void update_quadrature_point_history();
 * 
 * @endcode
 * 
 * This is the new shared Triangulation:
 * 
 * @code
 *     parallel::shared::Triangulation<dim> triangulation;
 * 
 *     FESystem<dim> fe;
 * 
 *     DoFHandler<dim> dof_handler;
 * 
 *     AffineConstraints<double> hanging_node_constraints;
 * 
 * @endcode
 * 
 * One difference of this program is that we declare the quadrature
 * formula in the class declaration. The reason is that in all the other
 * programs, it didn't do much harm if we had used different quadrature
 * formulas when computing the matrix and the right hand side, for
 * example. However, in the present case it does: we store information in
 * the quadrature points, so we have to make sure all parts of the program
 * agree on where they are and how many there are on each cell. Thus, let
 * us first declare the quadrature formula that will be used throughout...
 * 
 * @code
 *     const QGauss<dim> quadrature_formula;
 * 
 * @endcode
 * 
 * ... and then also have a vector of history objects, one per quadrature
 * point on those cells for which we are responsible (i.e. we don't store
 * history data for quadrature points on cells that are owned by other
 * processors).
 * Note that, instead of storing and managing this data ourself, we
 * could use the CellDataStorage class like is done in step-44. However,
 * for the purpose of demonstration, in this case we manage the storage
 * manually.
 * 
 * @code
 *     std::vector<PointHistory<dim>> quadrature_point_history;
 * 
 * @endcode
 * 
 * The way this object is accessed is through a <code>user pointer</code>
 * that each cell, face, or edge holds: it is a <code>void*</code> pointer
 * that can be used by application programs to associate arbitrary data to
 * cells, faces, or edges. What the program actually does with this data
 * is within its own responsibility, the library just allocates some space
 * for these pointers, and application programs can set and read the
 * pointers for each of these objects.
 * 

 * 
 * 

 * 
 * Further: we need the objects of linear systems to be solved,
 * i.e. matrix, right hand side vector, and the solution vector. Since we
 * anticipate solving big problems, we use the same types as in step-17,
 * i.e. distributed %parallel matrices and vectors built on top of the
 * PETSc library. Conveniently, they can also be used when running on only
 * a single machine, in which case this machine happens to be the only one
 * in our %parallel universe.
 *     

 * 
 * However, as a difference to step-17, we do not store the solution
 * vector -- which here is the incremental displacements computed in each
 * time step -- in a distributed fashion. I.e., of course it must be a
 * distributed vector when computing it, but immediately after that we
 * make sure each processor has a complete copy. The reason is that we had
 * already seen in step-17 that many functions needed a complete
 * copy. While it is not hard to get it, this requires communication on
 * the network, and is thus slow. In addition, these were repeatedly the
 * same operations, which is certainly undesirable unless the gains of not
 * always having to store the entire vector outweighs it. When writing
 * this program, it turned out that we need a complete copy of the
 * solution in so many places that it did not seem worthwhile to only get
 * it when necessary. Instead, we opted to obtain the complete copy once
 * and for all, and instead get rid of the distributed copy
 * immediately. Thus, note that the declaration of
 * <code>incremental_displacement</code> does not denote a distribute
 * vector as would be indicated by the middle namespace <code>MPI</code>:
 * 
 * @code
 *     PETScWrappers::MPI::SparseMatrix system_matrix;
 * 
 *     PETScWrappers::MPI::Vector system_rhs;
 * 
 *     Vector<double> incremental_displacement;
 * 
 * @endcode
 * 
 * The next block of variables is then related to the time dependent
 * nature of the problem: they denote the length of the time interval
 * which we want to simulate, the present time and number of time step,
 * and length of present timestep:
 * 
 * @code
 *     double       present_time;
 *     double       present_timestep;
 *     double       end_time;
 *     unsigned int timestep_no;
 * 
 * @endcode
 * 
 * Then a few variables that have to do with %parallel processing: first,
 * a variable denoting the MPI communicator we use, and then two numbers
 * telling us how many participating processors there are, and where in
 * this world we are. Finally, a stream object that makes sure only one
 * processor is actually generating output to the console. This is all the
 * same as in step-17:
 * 
 * @code
 *     MPI_Comm mpi_communicator;
 * 
 *     const unsigned int n_mpi_processes;
 * 
 *     const unsigned int this_mpi_process;
 * 
 *     ConditionalOStream pcout;
 * 
 * @endcode
 * 
 * We are storing the locally owned and the locally relevant indices:
 * 
 * @code
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 * @endcode
 * 
 * Finally, we have a static variable that denotes the linear relationship
 * between the stress and strain. Since it is a constant object that does
 * not depend on any input (at least not in this program), we make it a
 * static variable and will initialize it in the same place where we
 * define the constructor of this class:
 * 
 * @code
 *     static const SymmetricTensor<4, dim> stress_strain_tensor;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBodyForcecodeclass"></a> 
 * <h3>The <code>BodyForce</code> class</h3>
 * 

 * 
 * Before we go on to the main functionality of this program, we have to
 * define what forces will act on the body whose deformation we want to
 * study. These may either be body forces or boundary forces. Body forces
 * are generally mediated by one of the four basic physical types of forces:
 * gravity, strong and weak interaction, and electromagnetism. Unless one
 * wants to consider subatomic objects (for which quasistatic deformation is
 * irrelevant and an inappropriate description anyway), only gravity and
 * electromagnetic forces need to be considered. Let us, for simplicity
 * assume that our body has a certain mass density, but is either
 * non-magnetic and not electrically conducting or that there are no
 * significant electromagnetic fields around. In that case, the body forces
 * are simply <code>rho g</code>, where <code>rho</code> is the material
 * density and <code>g</code> is a vector in negative z-direction with
 * magnitude 9.81 m/s^2.  Both the density and <code>g</code> are defined in
 * the function, and we take as the density 7700 kg/m^3, a value commonly
 * assumed for steel.
 *   

 * 
 * To be a little more general and to be able to do computations in 2d as
 * well, we realize that the body force is always a function returning a
 * <code>dim</code> dimensional vector. We assume that gravity acts along
 * the negative direction of the last, i.e. <code>dim-1</code>th
 * coordinate. The rest of the implementation of this function should be
 * mostly self-explanatory given similar definitions in previous example
 * programs. Note that the body force is independent of the location; to
 * avoid compiler warnings about unused function arguments, we therefore
 * comment out the name of the first argument of the
 * <code>vector_value</code> function:
 * 
 * @code
 *   template <int dim>
 *   class BodyForce : public Function<dim>
 *   {
 *   public:
 *     BodyForce();
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override;
 * 
 *     virtual void
 *     vector_value_list(const std::vector<Point<dim>> &points,
 *                       std::vector<Vector<double>> &  value_list) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   BodyForce<dim>::BodyForce()
 *     : Function<dim>(dim)
 *   {}
 * 
 * 
 *   template <int dim>
 *   inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/,
 *                                            Vector<double> &values) const
 *   {
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
 * 
 *     const double g   = 9.81;
 *     const double rho = 7700;
 * 
 *     values          = 0;
 *     values(dim - 1) = -rho * g;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void BodyForce<dim>::vector_value_list(
 *     const std::vector<Point<dim>> &points,
 *     std::vector<Vector<double>> &  value_list) const
 *   {
 *     const unsigned int n_points = points.size();
 * 
 *     Assert(value_list.size() == n_points,
 *            ExcDimensionMismatch(value_list.size(), n_points));
 * 
 *     for (unsigned int p = 0; p < n_points; ++p)
 *       BodyForce<dim>::vector_value(points[p], value_list[p]);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeIncrementalBoundaryValuecodeclass"></a> 
 * <h3>The <code>IncrementalBoundaryValue</code> class</h3>
 * 

 * 
 * In addition to body forces, movement can be induced by boundary forces
 * and forced boundary displacement. The latter case is equivalent to forces
 * being chosen in such a way that they induce certain displacement.
 *   

 * 
 * For quasistatic displacement, typical boundary forces would be pressure
 * on a body, or tangential friction against another body. We chose a
 * somewhat simpler case here: we prescribe a certain movement of (parts of)
 * the boundary, or at least of certain components of the displacement
 * vector. We describe this by another vector-valued function that, for a
 * given point on the boundary, returns the prescribed displacement.
 *   

 * 
 * Since we have a time-dependent problem, the displacement increment of the
 * boundary equals the displacement accumulated during the length of the
 * timestep. The class therefore has to know both the present time and the
 * length of the present time step, and can then approximate the incremental
 * displacement as the present velocity times the present timestep.
 *   

 * 
 * For the purposes of this program, we choose a simple form of boundary
 * displacement: we displace the top boundary with constant velocity
 * downwards. The rest of the boundary is either going to be fixed (and is
 * then described using an object of type
 * <code>Functions::ZeroFunction</code>) or free (Neumann-type, in which case
 * nothing special has to be done).  The implementation of the class
 * describing the constant downward motion should then be obvious using the
 * knowledge we gained through all the previous example programs:
 * 
 * @code
 *   template <int dim>
 *   class IncrementalBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     IncrementalBoundaryValues(const double present_time,
 *                               const double present_timestep);
 * 
 *     virtual void vector_value(const Point<dim> &p,
 *                               Vector<double> &  values) const override;
 * 
 *     virtual void
 *     vector_value_list(const std::vector<Point<dim>> &points,
 *                       std::vector<Vector<double>> &  value_list) const override;
 * 
 *   private:
 *     const double velocity;
 *     const double present_time;
 *     const double present_timestep;
 *   };
 * 
 * 
 *   template <int dim>
 *   IncrementalBoundaryValues<dim>::IncrementalBoundaryValues(
 *     const double present_time,
 *     const double present_timestep)
 *     : Function<dim>(dim)
 *     , velocity(.08)
 *     , present_time(present_time)
 *     , present_timestep(present_timestep)
 *   {}
 * 
 * 
 *   template <int dim>
 *   void
 *   IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/,
 *                                                Vector<double> &values) const
 *   {
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
 * 
 *     values    = 0;
 *     values(2) = -present_timestep * velocity;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void IncrementalBoundaryValues<dim>::vector_value_list(
 *     const std::vector<Point<dim>> &points,
 *     std::vector<Vector<double>> &  value_list) const
 *   {
 *     const unsigned int n_points = points.size();
 * 
 *     Assert(value_list.size() == n_points,
 *            ExcDimensionMismatch(value_list.size(), n_points));
 * 
 *     for (unsigned int p = 0; p < n_points; ++p)
 *       IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeTopLevelcodeclass"></a> 
 * <h3>Implementation of the <code>TopLevel</code> class</h3>
 * 

 * 
 * Now for the implementation of the main class. First, we initialize the
 * stress-strain tensor, which we have declared as a static const
 * variable. We chose Lam&eacute; constants that are appropriate for steel:
 * 
 * @code
 *   template <int dim>
 *   const SymmetricTensor<4, dim> TopLevel<dim>::stress_strain_tensor =
 *     get_stress_strain_tensor<dim>(/*lambda = */ 9.695e10,
 *                                   /*mu     = */ 7.617e10);
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thepublicinterface"></a> 
 * <h4>The public interface</h4>
 * 

 * 
 * The next step is the definition of constructors and destructors. There
 * are no surprises here: we choose linear and continuous finite elements
 * for each of the <code>dim</code> vector components of the solution, and a
 * Gaussian quadrature formula with 2 points in each coordinate
 * direction. The destructor should be obvious:
 * 
 * @code
 *   template <int dim>
 *   TopLevel<dim>::TopLevel()
 *     : triangulation(MPI_COMM_WORLD)
 *     , fe(FE_Q<dim>(1), dim)
 *     , dof_handler(triangulation)
 *     , quadrature_formula(fe.degree + 1)
 *     , present_time(0.0)
 *     , present_timestep(1.0)
 *     , end_time(10.0)
 *     , timestep_no(0)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
 *     , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
 *     , pcout(std::cout, this_mpi_process == 0)
 *   {}
 * 
 * 
 * 
 *   template <int dim>
 *   TopLevel<dim>::~TopLevel()
 *   {
 *     dof_handler.clear();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The last of the public functions is the one that directs all the work,
 * <code>run()</code>. It initializes the variables that describe where in
 * time we presently are, then runs the first time step, then loops over all
 * the other time steps. Note that for simplicity we use a fixed time step,
 * whereas a more sophisticated program would of course have to choose it in
 * some more reasonable way adaptively:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::run()
 *   {
 *     do_initial_timestep();
 * 
 *     while (present_time < end_time)
 *       do_timestep();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelcreate_coarse_grid"></a> 
 * <h4>TopLevel::create_coarse_grid</h4>
 * 

 * 
 * The next function in the order in which they were declared above is the
 * one that creates the coarse grid from which we start. For this example
 * program, we want to compute the deformation of a cylinder under axial
 * compression. The first step therefore is to generate a mesh for a
 * cylinder of length 3 and with inner and outer radii of 0.8 and 1,
 * respectively. Fortunately, there is a library function for such a mesh.
 *   

 * 
 * In a second step, we have to associated boundary conditions with the
 * upper and lower faces of the cylinder. We choose a boundary indicator of
 * 0 for the boundary faces that are characterized by their midpoints having
 * z-coordinates of either 0 (bottom face), an indicator of 1 for z=3 (top
 * face); finally, we use boundary indicator 2 for all faces on the inside
 * of the cylinder shell, and 3 for the outside.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::create_coarse_grid()
 *   {
 *     const double inner_radius = 0.8, outer_radius = 1;
 *     GridGenerator::cylinder_shell(triangulation, 3, inner_radius, outer_radius);
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             const Point<dim> face_center = face->center();
 * 
 *             if (face_center[2] == 0)
 *               face->set_boundary_id(0);
 *             else if (face_center[2] == 3)
 *               face->set_boundary_id(1);
 *             else if (std::sqrt(face_center[0] * face_center[0] +
 *                                face_center[1] * face_center[1]) <
 *                      (inner_radius + outer_radius) / 2)
 *               face->set_boundary_id(2);
 *             else
 *               face->set_boundary_id(3);
 *           }
 * 
 * @endcode
 * 
 * Once all this is done, we can refine the mesh once globally:
 * 
 * @code
 *     triangulation.refine_global(1);
 * 
 * @endcode
 * 
 * As the final step, we need to set up a clean state of the data that we
 * store in the quadrature points on all cells that are treated on the
 * present processor.
 * 
 * @code
 *     setup_quadrature_point_history();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_system"></a> 
 * <h4>TopLevel::setup_system</h4>
 * 

 * 
 * The next function is the one that sets up the data structures for a given
 * mesh. This is done in most the same way as in step-17: distribute the
 * degrees of freedom, then sort these degrees of freedom in such a way that
 * each processor gets a contiguous chunk of them. Note that subdivisions into
 * chunks for each processor is handled in the functions that create or
 * refine grids, unlike in the previous example program (the point where
 * this happens is mostly a matter of taste; here, we chose to do it when
 * grids are created since in the <code>do_initial_timestep</code> and
 * <code>do_timestep</code> functions we want to output the number of cells
 * on each processor at a point where we haven't called the present function
 * yet).
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 * @endcode
 * 
 * The next step is to set up constraints due to hanging nodes. This has
 * been handled many times before:
 * 
 * @code
 *     hanging_node_constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler,
 *                                             hanging_node_constraints);
 *     hanging_node_constraints.close();
 * 
 * @endcode
 * 
 * And then we have to set up the matrix. Here we deviate from step-17, in
 * which we simply used PETSc's ability to just know about the size of the
 * matrix and later allocate those nonzero elements that are being written
 * to. While this works just fine from a correctness viewpoint, it is not
 * at all efficient: if we don't give PETSc a clue as to which elements
 * are written to, it is (at least at the time of this writing) unbearably
 * slow when we set the elements in the matrix for the first time (i.e. in
 * the first timestep). Later on, when the elements have been allocated,
 * everything is much faster. In experiments we made, the first timestep
 * can be accelerated by almost two orders of magnitude if we instruct
 * PETSc which elements will be used and which are not.
 *     

 * 
 * To do so, we first generate the sparsity pattern of the matrix we are
 * going to work with, and make sure that the condensation of hanging node
 * constraints add the necessary additional entries in the sparsity
 * pattern:
 * 
 * @code
 *     DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs);
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     sparsity_pattern,
 *                                     hanging_node_constraints,
 *                                     /*keep constrained dofs*/ false);
 *     SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
 *                                                locally_owned_dofs,
 *                                                mpi_communicator,
 *                                                locally_relevant_dofs);
 * @endcode
 * 
 * Note that we have used the <code>DynamicSparsityPattern</code> class
 * here that was already introduced in step-11, rather than the
 * <code>SparsityPattern</code> class that we have used in all other
 * cases. The reason for this is that for the latter class to work we have
 * to give an initial upper bound for the number of entries in each row, a
 * task that is traditionally done by
 * <code>DoFHandler::max_couplings_between_dofs()</code>. However, this
 * function suffers from a serious problem: it has to compute an upper
 * bound to the number of nonzero entries in each row, and this is a
 * rather complicated task, in particular in 3d. In effect, while it is
 * quite accurate in 2d, it often comes up with much too large a number in
 * 3d, and in that case the <code>SparsityPattern</code> allocates much
 * too much memory at first, often several 100 MBs. This is later
 * corrected when <code>DoFTools::make_sparsity_pattern</code> is called
 * and we realize that we don't need all that much memory, but at time it
 * is already too late: for large problems, the temporary allocation of
 * too much memory can lead to out-of-memory situations.
 *     

 * 
 * In order to avoid this, we resort to the
 * <code>DynamicSparsityPattern</code> class that is slower but does
 * not require any up-front estimate on the number of nonzero entries per
 * row. It therefore only ever allocates as much memory as it needs at any
 * given time, and we can build it even for large 3d problems.
 *     

 * 
 * It is also worth noting that due to the specifics of
 * parallel::shared::Triangulation, the sparsity pattern we construct is
 * global, i.e. comprises all degrees of freedom whether they will be
 * owned by the processor we are on or another one (in case this program
 * is run in %parallel via MPI). This of course is not optimal -- it
 * limits the size of the problems we can solve, since storing the entire
 * sparsity pattern (even if only for a short time) on each processor does
 * not scale well. However, there are several more places in the program
 * in which we do this, for example we always keep the global
 * triangulation and DoF handler objects around, even if we only work on
 * part of them. At present, deal.II does not have the necessary
 * facilities to completely distribute these objects (a task that, indeed,
 * is very hard to achieve with adaptive meshes, since well-balanced
 * subdivisions of a domain tend to become unbalanced as the mesh is
 * adaptively refined).
 *     

 * 
 * With this data structure, we can then go to the PETSc sparse matrix and
 * tell it to preallocate all the entries we will later want to write to:
 * 
 * @code
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          sparsity_pattern,
 *                          mpi_communicator);
 * @endcode
 * 
 * After this point, no further explicit knowledge of the sparsity pattern
 * is required any more and we can let the <code>sparsity_pattern</code>
 * variable go out of scope without any problem.
 * 

 * 
 * The last task in this function is then only to reset the right hand
 * side vector as well as the solution vector to its correct size;
 * remember that the solution vector is a local one, unlike the right hand
 * side that is a distributed %parallel one and therefore needs to know
 * the MPI communicator over which it is supposed to transmit messages:
 * 
 * @code
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *     incremental_displacement.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelassemble_system"></a> 
 * <h4>TopLevel::assemble_system</h4>
 * 

 * 
 * Again, assembling the system matrix and right hand side follows the same
 * structure as in many example programs before. In particular, it is mostly
 * equivalent to step-17, except for the different right hand side that now
 * only has to take into account internal stresses. In addition, assembling
 * the matrix is made significantly more transparent by using the
 * <code>SymmetricTensor</code> class: note the elegance of forming the
 * scalar products of symmetric tensors of rank 2 and 4. The implementation
 * is also more general since it is independent of the fact that we may or
 * may not be using an isotropic elasticity tensor.
 *   

 * 
 * The first part of the assembly routine is as always:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::assemble_system()
 *   {
 *     system_rhs    = 0;
 *     system_matrix = 0;
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     BodyForce<dim>              body_force;
 *     std::vector<Vector<double>> body_force_values(n_q_points,
 *                                                   Vector<double>(dim));
 * 
 * @endcode
 * 
 * As in step-17, we only need to loop over all cells that belong to the
 * present processor:
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 *           fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Then loop over all indices i,j and quadrature points and assemble
 * the system matrix contributions from this cell.  Note how we
 * extract the symmetric gradients (strains) of the shape functions
 * at a given quadrature point from the <code>FEValues</code>
 * object, and the elegance with which we form the triple
 * contraction <code>eps_phi_i : C : eps_phi_j</code>; the latter
 * needs to be compared to the clumsy computations needed in
 * step-17, both in the introduction as well as in the respective
 * place in the program:
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *                 {
 *                   const SymmetricTensor<2, dim>
 *                     eps_phi_i = get_strain(fe_values, i, q_point),
 *                     eps_phi_j = get_strain(fe_values, j, q_point);
 * 
 *                   cell_matrix(i, j) += (eps_phi_i *            
 *                                         stress_strain_tensor * 
 *                                         eps_phi_j              
 *                                         ) *                    
 *                                        fe_values.JxW(q_point); 
 *                 }
 * 
 * 
 * @endcode
 * 
 * Then also assemble the local right hand side contributions. For
 * this, we need to access the prior stress value in this quadrature
 * point. To get it, we use the user pointer of this cell that
 * points into the global array to the quadrature point data
 * corresponding to the first quadrature point of the present cell,
 * and then add an offset corresponding to the index of the
 * quadrature point we presently consider:
 * 
 * @code
 *           const PointHistory<dim> *local_quadrature_points_data =
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
 * @endcode
 * 
 * In addition, we need the values of the external body forces at
 * the quadrature points on this cell:
 * 
 * @code
 *           body_force.vector_value_list(fe_values.get_quadrature_points(),
 *                                        body_force_values);
 * @endcode
 * 
 * Then we can loop over all degrees of freedom on this cell and
 * compute local contributions to the right hand side:
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               const unsigned int component_i =
 *                 fe.system_to_component_index(i).first;
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *                 {
 *                   const SymmetricTensor<2, dim> &old_stress =
 *                     local_quadrature_points_data[q_point].old_stress;
 * 
 *                   cell_rhs(i) +=
 *                     (body_force_values[q_point](component_i) *
 *                        fe_values.shape_value(i, q_point) -
 *                      old_stress * get_strain(fe_values, i, q_point)) *
 *                     fe_values.JxW(q_point);
 *                 }
 *             }
 * 
 * @endcode
 * 
 * Now that we have the local contributions to the linear system, we
 * need to transfer it into the global objects. This is done exactly
 * as in step-17:
 * 
 * @code
 *           cell->get_dof_indices(local_dof_indices);
 * 
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix,
 *                                                               cell_rhs,
 *                                                               local_dof_indices,
 *                                                               system_matrix,
 *                                                               system_rhs);
 *         }
 * 
 * @endcode
 * 
 * Now compress the vector and the system matrix:
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 * 
 * 
 * @endcode
 * 
 * The last step is to again fix up boundary values, just as we already
 * did in previous programs. A slight complication is that the
 * <code>apply_boundary_values</code> function wants to have a solution
 * vector compatible with the matrix and right hand side (i.e. here a
 * distributed %parallel vector, rather than the sequential vector we use
 * in this program) in order to preset the entries of the solution vector
 * with the correct boundary values. We provide such a compatible vector
 * in the form of a temporary vector which we then copy into the
 * sequential one.
 * 

 * 
 * We make up for this complication by showing how boundary values can be
 * used flexibly: following the way we create the triangulation, there are
 * three distinct boundary indicators used to describe the domain,
 * corresponding to the bottom and top faces, as well as the inner/outer
 * surfaces. We would like to impose boundary conditions of the following
 * type: The inner and outer cylinder surfaces are free of external
 * forces, a fact that corresponds to natural (Neumann-type) boundary
 * conditions for which we don't have to do anything. At the bottom, we
 * want no movement at all, corresponding to the cylinder being clamped or
 * cemented in at this part of the boundary. At the top, however, we want
 * a prescribed vertical downward motion compressing the cylinder; in
 * addition, we only want to restrict the vertical movement, but not the
 * horizontal ones -- one can think of this situation as a well-greased
 * plate sitting on top of the cylinder pushing it downwards: the atoms of
 * the cylinder are forced to move downward, but they are free to slide
 * horizontally along the plate.
 * 

 * 
 * The way to describe this is as follows: for boundary indicator zero
 * (bottom face) we use a dim-dimensional zero function representing no
 * motion in any coordinate direction. For the boundary with indicator 1
 * (top surface), we use the <code>IncrementalBoundaryValues</code> class,
 * but we specify an additional argument to the
 * <code>VectorTools::interpolate_boundary_values</code> function denoting
 * which vector components it should apply to; this is a vector of bools
 * for each vector component and because we only want to restrict vertical
 * motion, it has only its last component set:
 * 
 * @code
 *     FEValuesExtractors::Scalar                z_component(dim - 1);
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(dim),
 *                                              boundary_values);
 *     VectorTools::interpolate_boundary_values(
 *       dof_handler,
 *       1,
 *       IncrementalBoundaryValues<dim>(present_time, present_timestep),
 *       boundary_values,
 *       fe.component_mask(z_component));
 * 
 *     PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
 *     MatrixTools::apply_boundary_values(
 *       boundary_values, system_matrix, tmp, system_rhs, false);
 *     incremental_displacement = tmp;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_timestep"></a> 
 * <h4>TopLevel::solve_timestep</h4>
 * 

 * 
 * The next function is the one that controls what all has to happen within
 * a timestep. The order of things should be relatively self-explanatory
 * from the function names:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::solve_timestep()
 *   {
 *     pcout << "    Assembling system..." << std::flush;
 *     assemble_system();
 *     pcout << " norm of rhs is " << system_rhs.l2_norm() << std::endl;
 * 
 *     const unsigned int n_iterations = solve_linear_problem();
 * 
 *     pcout << "    Solver converged in " << n_iterations << " iterations."
 *           << std::endl;
 * 
 *     pcout << "    Updating quadrature point data..." << std::flush;
 *     update_quadrature_point_history();
 *     pcout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_linear_problem"></a> 
 * <h4>TopLevel::solve_linear_problem</h4>
 * 

 * 
 * Solving the linear system again works mostly as before. The only
 * difference is that we want to only keep a complete local copy of the
 * solution vector instead of the distributed one that we get as output from
 * PETSc's solver routines. To this end, we declare a local temporary
 * variable for the distributed vector and initialize it with the contents
 * of the local variable (remember that the
 * <code>apply_boundary_values</code> function called in
 * <code>assemble_system</code> preset the values of boundary nodes in this
 * vector), solve with it, and at the end of the function copy it again into
 * the complete local vector that we declared as a member variable. Hanging
 * node constraints are then distributed only on the local copy,
 * i.e. independently of each other on each of the processors:
 * 
 * @code
 *   template <int dim>
 *   unsigned int TopLevel<dim>::solve_linear_problem()
 *   {
 *     PETScWrappers::MPI::Vector distributed_incremental_displacement(
 *       locally_owned_dofs, mpi_communicator);
 *     distributed_incremental_displacement = incremental_displacement;
 * 
 *     SolverControl solver_control(dof_handler.n_dofs(),
 *                                  1e-16 * system_rhs.l2_norm());
 * 
 *     PETScWrappers::SolverCG cg(solver_control, mpi_communicator);
 * 
 *     PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
 * 
 *     cg.solve(system_matrix,
 *              distributed_incremental_displacement,
 *              system_rhs,
 *              preconditioner);
 * 
 *     incremental_displacement = distributed_incremental_displacement;
 * 
 *     hanging_node_constraints.distribute(incremental_displacement);
 * 
 *     return solver_control.last_step();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveloutput_results"></a> 
 * <h4>TopLevel::output_results</h4>
 * 

 * 
 * This function generates the graphical output in .vtu format as explained
 * in the introduction. Each process will only work on the cells it owns,
 * and then write the result into a file of its own. Additionally, processor
 * 0 will write the record files the reference all the .vtu files.
 *   

 * 
 * The crucial part of this function is to give the <code>DataOut</code>
 * class a way to only work on the cells that the present process owns.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 * 
 * @endcode
 * 
 * Then, just as in step-17, define the names of solution variables (which
 * here are the displacement increments) and queue the solution vector for
 * output. Note in the following switch how we make sure that if the space
 * dimension should be unhandled that we throw an exception saying that we
 * haven't implemented this case yet (another case of defensive
 * programming):
 * 
 * @code
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 1:
 *           solution_names.emplace_back("delta_x");
 *           break;
 *         case 2:
 *           solution_names.emplace_back("delta_x");
 *           solution_names.emplace_back("delta_y");
 *           break;
 *         case 3:
 *           solution_names.emplace_back("delta_x");
 *           solution_names.emplace_back("delta_y");
 *           solution_names.emplace_back("delta_z");
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 *     data_out.add_data_vector(incremental_displacement, solution_names);
 * 
 * 
 * @endcode
 * 
 * The next thing is that we wanted to output something like the average
 * norm of the stresses that we have stored in each cell. This may seem
 * complicated, since on the present processor we only store the stresses
 * in quadrature points on those cells that actually belong to the present
 * process. In other words, it seems as if we can't compute the average
 * stresses for all cells. However, remember that our class derived from
 * <code>DataOut</code> only iterates over those cells that actually do
 * belong to the present processor, i.e. we don't have to compute anything
 * for all the other cells as this information would not be touched. The
 * following little loop does this. We enclose the entire block into a
 * pair of braces to make sure that the iterator variables do not remain
 * accidentally visible beyond the end of the block in which they are
 * used:
 * 
 * @code
 *     Vector<double> norm_of_stress(triangulation.n_active_cells());
 *     {
 * @endcode
 * 
 * Loop over all the cells...
 * 
 * @code
 *       for (auto &cell : triangulation.active_cell_iterators())
 *         if (cell->is_locally_owned())
 *           {
 * @endcode
 * 
 * On these cells, add up the stresses over all quadrature
 * points...
 * 
 * @code
 *             SymmetricTensor<2, dim> accumulated_stress;
 *             for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
 *               accumulated_stress +=
 *                 reinterpret_cast<PointHistory<dim> *>(cell->user_pointer())[q]
 *                   .old_stress;
 * 
 * @endcode
 * 
 * ...then write the norm of the average to their destination:
 * 
 * @code
 *             norm_of_stress(cell->active_cell_index()) =
 *               (accumulated_stress / quadrature_formula.size()).norm();
 *           }
 * @endcode
 * 
 * And on the cells that we are not interested in, set the respective
 * value in the vector to a bogus value (norms must be positive, and a
 * large negative value should catch your eye) in order to make sure
 * that if we were somehow wrong about our assumption that these
 * elements would not appear in the output file, that we would find out
 * by looking at the graphical output:
 * 
 * @code
 *         else
 *           norm_of_stress(cell->active_cell_index()) = -1e+20;
 *     }
 * @endcode
 * 
 * Finally attach this vector as well to be treated for output:
 * 
 * @code
 *     data_out.add_data_vector(norm_of_stress, "norm_of_stress");
 * 
 * @endcode
 * 
 * As a last piece of data, let us also add the partitioning of the domain
 * into subdomains associated with the processors if this is a parallel
 * job. This works in the exact same way as in the step-17 program:
 * 
 * @code
 *     std::vector<types::subdomain_id> partition_int(
 *       triangulation.n_active_cells());
 *     GridTools::get_subdomain_association(triangulation, partition_int);
 *     const Vector<double> partitioning(partition_int.begin(),
 *                                       partition_int.end());
 *     data_out.add_data_vector(partitioning, "partitioning");
 * 
 * @endcode
 * 
 * Finally, with all this data, we can instruct deal.II to munge the
 * information and produce some intermediate data structures that contain
 * all these solution and other data vectors:
 * 
 * @code
 *     data_out.build_patches();
 * 
 * @endcode
 * 
 * Let us call a function that opens the necessary output files and writes
 * the data we have generated into them. The function automatically
 * constructs the file names from the given directory name (the first
 * argument) and file name base (second argument). It augments the resulting
 * string by pieces that result from the time step number and a "piece
 * number" that corresponds to a part of the overall domain that can consist
 * of one or more subdomains.
 *     

 * 
 * The function also writes a record files (with suffix `.pvd`) for Paraview
 * that describes how all of these output files combine into the data for
 * this single time step:
 * 
 * @code
 *     const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", timestep_no, mpi_communicator, 4);
 * 
 * @endcode
 * 
 * The record files must be written only once and not by each processor,
 * so we do this on processor 0:
 * 
 * @code
 *     if (this_mpi_process == 0)
 *       {
 * @endcode
 * 
 * Finally, we write the paraview record, that references all .pvtu
 * files and their respective time. Note that the variable
 * times_and_names is declared static, so it will retain the entries
 * from the previous timesteps.
 * 
 * @code
 *         static std::vector<std::pair<double, std::string>> times_and_names;
 *         times_and_names.push_back(
 *           std::pair<double, std::string>(present_time, pvtu_filename));
 *         std::ofstream pvd_output("solution.pvd");
 *         DataOutBase::write_pvd_record(pvd_output, times_and_names);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_initial_timestep"></a> 
 * <h4>TopLevel::do_initial_timestep</h4>
 * 

 * 
 * This and the next function handle the overall structure of the first and
 * following timesteps, respectively. The first timestep is slightly more
 * involved because we want to compute it multiple times on successively
 * refined meshes, each time starting from a clean state. At the end of
 * these computations, in which we compute the incremental displacements
 * each time, we use the last results obtained for the incremental
 * displacements to compute the resulting stress updates and move the mesh
 * accordingly. On this new mesh, we then output the solution and any
 * additional data we consider important.
 *   

 * 
 * All this is interspersed by generating output to the console to update
 * the person watching the screen on what is going on. As in step-17, the
 * use of <code>pcout</code> instead of <code>std::cout</code> makes sure
 * that only one of the parallel processes is actually writing to the
 * console, without having to explicitly code an if-statement in each place
 * where we generate output:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::do_initial_timestep()
 *   {
 *     present_time += present_timestep;
 *     ++timestep_no;
 *     pcout << "Timestep " << timestep_no << " at time " << present_time
 *           << std::endl;
 * 
 *     for (unsigned int cycle = 0; cycle < 2; ++cycle)
 *       {
 *         pcout << "  Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           create_coarse_grid();
 *         else
 *           refine_initial_grid();
 * 
 *         pcout << "    Number of active cells:       "
 *               << triangulation.n_active_cells() << " (by partition:";
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p)
 *           pcout << (p == 0 ? ' ' : '+')
 *                 << (GridTools::count_cells_with_subdomain_association(
 *                      triangulation, p));
 *         pcout << ")" << std::endl;
 * 
 *         setup_system();
 * 
 *         pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << " (by partition:";
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p)
 *           pcout << (p == 0 ? ' ' : '+')
 *                 << (DoFTools::count_dofs_with_subdomain_association(dof_handler,
 *                                                                     p));
 *         pcout << ")" << std::endl;
 * 
 *         solve_timestep();
 *       }
 * 
 *     move_mesh();
 *     output_results();
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_timestep"></a> 
 * <h4>TopLevel::do_timestep</h4>
 * 

 * 
 * Subsequent timesteps are simpler, and probably do not require any more
 * documentation given the explanations for the previous function above:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::do_timestep()
 *   {
 *     present_time += present_timestep;
 *     ++timestep_no;
 *     pcout << "Timestep " << timestep_no << " at time " << present_time
 *           << std::endl;
 *     if (present_time > end_time)
 *       {
 *         present_timestep -= (present_time - end_time);
 *         present_time = end_time;
 *       }
 * 
 * 
 *     solve_timestep();
 * 
 *     move_mesh();
 *     output_results();
 * 
 *     pcout << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelrefine_initial_grid"></a> 
 * <h4>TopLevel::refine_initial_grid</h4>
 * 

 * 
 * The following function is called when solving the first time step on
 * successively refined meshes. After each iteration, it computes a
 * refinement criterion, refines the mesh, and sets up the history variables
 * in each quadrature point again to a clean state.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::refine_initial_grid()
 *   {
 * @endcode
 * 
 * First, let each process compute error indicators for the cells it owns:
 * 
 * @code
 *     Vector<float> error_per_cell(triangulation.n_active_cells());
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       QGauss<dim - 1>(fe.degree + 1),
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       incremental_displacement,
 *       error_per_cell,
 *       ComponentMask(),
 *       nullptr,
 *       MultithreadInfo::n_threads(),
 *       this_mpi_process);
 * 
 * @endcode
 * 
 * Then set up a global vector into which we merge the local indicators
 * from each of the %parallel processes:
 * 
 * @code
 *     const unsigned int n_local_cells =
 *       triangulation.n_locally_owned_active_cells();
 * 
 *     PETScWrappers::MPI::Vector distributed_error_per_cell(
 *       mpi_communicator, triangulation.n_active_cells(), n_local_cells);
 * 
 *     for (unsigned int i = 0; i < error_per_cell.size(); ++i)
 *       if (error_per_cell(i) != 0)
 *         distributed_error_per_cell(i) = error_per_cell(i);
 *     distributed_error_per_cell.compress(VectorOperation::insert);
 * 
 * @endcode
 * 
 * Once we have that, copy it back into local copies on all processors and
 * refine the mesh accordingly:
 * 
 * @code
 *     error_per_cell = distributed_error_per_cell;
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     error_per_cell,
 *                                                     0.35,
 *                                                     0.03);
 *     triangulation.execute_coarsening_and_refinement();
 * 
 * @endcode
 * 
 * Finally, set up quadrature point data again on the new mesh, and only
 * on those cells that we have determined to be ours:
 * 
 * @code
 *     setup_quadrature_point_history();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelmove_mesh"></a> 
 * <h4>TopLevel::move_mesh</h4>
 * 

 * 
 * At the end of each time step, we move the nodes of the mesh according to
 * the incremental displacements computed in this time step. To do this, we
 * keep a vector of flags that indicate for each vertex whether we have
 * already moved it around, and then loop over all cells and move those
 * vertices of the cell that have not been moved yet. It is worth noting
 * that it does not matter from which of the cells adjacent to a vertex we
 * move this vertex: since we compute the displacement using a continuous
 * finite element, the displacement field is continuous as well and we can
 * compute the displacement of a given vertex from each of the adjacent
 * cells. We only have to make sure that we move each node exactly once,
 * which is why we keep the vector of flags.
 *   

 * 
 * There are two noteworthy things in this function. First, how we get the
 * displacement field at a given vertex using the
 * <code>cell-@>vertex_dof_index(v,d)</code> function that returns the index
 * of the <code>d</code>th degree of freedom at vertex <code>v</code> of the
 * given cell. In the present case, displacement in the k-th coordinate
 * direction corresponds to the k-th component of the finite element. Using a
 * function like this bears a certain risk, because it uses knowledge of the
 * order of elements that we have taken together for this program in the
 * <code>FESystem</code> element. If we decided to add an additional
 * variable, for example a pressure variable for stabilization, and happened
 * to insert it as the first variable of the element, then the computation
 * below will start to produce nonsensical results. In addition, this
 * computation rests on other assumptions: first, that the element we use
 * has, indeed, degrees of freedom that are associated with vertices. This
 * is indeed the case for the present Q1 element, as would be for all Qp
 * elements of polynomial order <code>p</code>. However, it would not hold
 * for discontinuous elements, or elements for mixed formulations. Secondly,
 * it also rests on the assumption that the displacement at a vertex is
 * determined solely by the value of the degree of freedom associated with
 * this vertex; in other words, all shape functions corresponding to other
 * degrees of freedom are zero at this particular vertex. Again, this is the
 * case for the present element, but is not so for all elements that are
 * presently available in deal.II. Despite its risks, we choose to use this
 * way in order to present a way to query individual degrees of freedom
 * associated with vertices.
 *   

 * 
 * In this context, it is instructive to point out what a more general way
 * would be. For general finite elements, the way to go would be to take a
 * quadrature formula with the quadrature points in the vertices of a
 * cell. The <code>QTrapezoid</code> formula for the trapezoidal rule does
 * exactly this. With this quadrature formula, we would then initialize an
 * <code>FEValues</code> object in each cell, and use the
 * <code>FEValues::get_function_values</code> function to obtain the values
 * of the solution function in the quadrature points, i.e. the vertices of
 * the cell. These are the only values that we really need, i.e. we are not
 * at all interested in the weights (or the <code>JxW</code> values)
 * associated with this particular quadrature formula, and this can be
 * specified as the last argument in the constructor to
 * <code>FEValues</code>. The only point of minor inconvenience in this
 * scheme is that we have to figure out which quadrature point corresponds
 * to the vertex we consider at present, as they may or may not be ordered
 * in the same order.
 *   

 * 
 * This inconvenience could be avoided if finite elements have support
 * points on vertices (which the one here has; for the concept of support
 * points, see @ref GlossSupport "support points"). For such a case, one
 * could construct a custom quadrature rule using
 * FiniteElement::get_unit_support_points(). The first
 * <code>cell-&gt;n_vertices()*fe.dofs_per_vertex</code>
 * quadrature points will then correspond to the vertices of the cell and
 * are ordered consistent with <code>cell-@>vertex(i)</code>, taking into
 * account that support points for vector elements will be duplicated
 * <code>fe.dofs_per_vertex</code> times.
 *   

 * 
 * Another point worth explaining about this short function is the way in
 * which the triangulation class exports information about its vertices:
 * through the <code>Triangulation::n_vertices</code> function, it
 * advertises how many vertices there are in the triangulation. Not all of
 * them are actually in use all the time -- some are left-overs from cells
 * that have been coarsened previously and remain in existence since deal.II
 * never changes the number of a vertex once it has come into existence,
 * even if vertices with lower number go away. Secondly, the location
 * returned by <code>cell-@>vertex(v)</code> is not only a read-only object
 * of type <code>Point@<dim@></code>, but in fact a reference that can be
 * written to. This allows to move around the nodes of a mesh with relative
 * ease, but it is worth pointing out that it is the responsibility of an
 * application program using this feature to make sure that the resulting
 * cells are still useful, i.e. are not distorted so much that the cell is
 * degenerated (indicated, for example, by negative Jacobians). Note that we
 * do not have any provisions in this function to actually ensure this, we
 * just have faith.
 *   

 * 
 * After this lengthy introduction, here are the full 20 or so lines of
 * code:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::move_mesh()
 *   {
 *     pcout << "    Moving mesh..." << std::endl;
 * 
 *     std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
 *     for (auto &cell : dof_handler.active_cell_iterators())
 *       for (const auto v : cell->vertex_indices())
 *         if (vertex_touched[cell->vertex_index(v)] == false)
 *           {
 *             vertex_touched[cell->vertex_index(v)] = true;
 * 
 *             Point<dim> vertex_displacement;
 *             for (unsigned int d = 0; d < dim; ++d)
 *               vertex_displacement[d] =
 *                 incremental_displacement(cell->vertex_dof_index(v, d));
 * 
 *             cell->vertex(v) += vertex_displacement;
 *           }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_quadrature_point_history"></a> 
 * <h4>TopLevel::setup_quadrature_point_history</h4>
 * 

 * 
 * At the beginning of our computations, we needed to set up initial values
 * of the history variables, such as the existing stresses in the material,
 * that we store in each quadrature point. As mentioned above, we use the
 * <code>user_pointer</code> for this that is available in each cell.
 *   

 * 
 * To put this into larger perspective, we note that if we had previously
 * available stresses in our model (which we assume do not exist for the
 * purpose of this program), then we would need to interpolate the field of
 * preexisting stresses to the quadrature points. Likewise, if we were to
 * simulate elasto-plastic materials with hardening/softening, then we would
 * have to store additional history variables like the present yield stress
 * of the accumulated plastic strains in each quadrature
 * points. Pre-existing hardening or weakening would then be implemented by
 * interpolating these variables in the present function as well.
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::setup_quadrature_point_history()
 *   {
 * @endcode
 * 
 * For good measure, we set all user pointers of all cells, whether
 * ours of not, to the null pointer. This way, if we ever access the user
 * pointer of a cell which we should not have accessed, a segmentation
 * fault will let us know that this should not have happened:
 * 

 * 
 * 
 * @code
 *     triangulation.clear_user_data();
 * 
 * @endcode
 * 
 * Next, allocate the quadrature objects that are within the responsibility
 * of this processor. This, of course, equals the number of cells that
 * belong to this processor times the number of quadrature points our
 * quadrature formula has on each cell. Since the `resize()` function does
 * not actually shrink the amount of allocated memory if the requested new
 * size is smaller than the old size, we resort to a trick to first free all
 * memory, and then reallocate it: we declare an empty vector as a temporary
 * variable and then swap the contents of the old vector and this temporary
 * variable. This makes sure that the `quadrature_point_history` is now
 * really empty, and we can let the temporary variable that now holds the
 * previous contents of the vector go out of scope and be destroyed. In the
 * next step we can then re-allocate as many elements as we need, with the
 * vector default-initializing the `PointHistory` objects, which includes
 * setting the stress variables to zero.
 * 
 * @code
 *     {
 *       std::vector<PointHistory<dim>> tmp;
 *       quadrature_point_history.swap(tmp);
 *     }
 *     quadrature_point_history.resize(
 *       triangulation.n_locally_owned_active_cells() * quadrature_formula.size());
 * 
 * @endcode
 * 
 * Finally loop over all cells again and set the user pointers from the
 * cells that belong to the present processor to point to the first
 * quadrature point objects corresponding to this cell in the vector of
 * such objects:
 * 
 * @code
 *     unsigned int history_index = 0;
 *     for (auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell->set_user_pointer(&quadrature_point_history[history_index]);
 *           history_index += quadrature_formula.size();
 *         }
 * 
 * @endcode
 * 
 * At the end, for good measure make sure that our count of elements was
 * correct and that we have both used up all objects we allocated
 * previously, and not point to any objects beyond the end of the
 * vector. Such defensive programming strategies are always good checks to
 * avoid accidental errors and to guard against future changes to this
 * function that forget to update all uses of a variable at the same
 * time. Recall that constructs using the <code>Assert</code> macro are
 * optimized away in optimized mode, so do not affect the run time of
 * optimized runs:
 * 
 * @code
 *     Assert(history_index == quadrature_point_history.size(),
 *            ExcInternalError());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelupdate_quadrature_point_history"></a> 
 * <h4>TopLevel::update_quadrature_point_history</h4>
 * 

 * 
 * At the end of each time step, we should have computed an incremental
 * displacement update so that the material in its new configuration
 * accommodates for the difference between the external body and boundary
 * forces applied during this time step minus the forces exerted through
 * preexisting internal stresses. In order to have the preexisting
 * stresses available at the next time step, we therefore have to update the
 * preexisting stresses with the stresses due to the incremental
 * displacement computed during the present time step. Ideally, the
 * resulting sum of internal stresses would exactly counter all external
 * forces. Indeed, a simple experiment can make sure that this is so: if we
 * choose boundary conditions and body forces to be time independent, then
 * the forcing terms (the sum of external forces and internal stresses)
 * should be exactly zero. If you make this experiment, you will realize
 * from the output of the norm of the right hand side in each time step that
 * this is almost the case: it is not exactly zero, since in the first time
 * step the incremental displacement and stress updates were computed
 * relative to the undeformed mesh, which was then deformed. In the second
 * time step, we again compute displacement and stress updates, but this
 * time in the deformed mesh -- there, the resulting updates are very small
 * but not quite zero. This can be iterated, and in each such iteration the
 * residual, i.e. the norm of the right hand side vector, is reduced; if one
 * makes this little experiment, one realizes that the norm of this residual
 * decays exponentially with the number of iterations, and after an initial
 * very rapid decline is reduced by roughly a factor of about 3.5 in each
 * iteration (for one testcase I looked at, other testcases, and other
 * numbers of unknowns change the factor, but not the exponential decay).
 * 

 * 
 * In a sense, this can then be considered as a quasi-timestepping scheme to
 * resolve the nonlinear problem of solving large-deformation elasticity on
 * a mesh that is moved along in a Lagrangian manner.
 *   

 * 
 * Another complication is that the existing (old) stresses are defined on
 * the old mesh, which we will move around after updating the stresses. If
 * this mesh update involves rotations of the cell, then we need to also
 * rotate the updated stress, since it was computed relative to the
 * coordinate system of the old cell.
 *   

 * 
 * Thus, what we need is the following: on each cell which the present
 * processor owns, we need to extract the old stress from the data stored
 * with each quadrature point, compute the stress update, add the two
 * together, and then rotate the result together with the incremental
 * rotation computed from the incremental displacement at the present
 * quadrature point. We will detail these steps below:
 * 
 * @code
 *   template <int dim>
 *   void TopLevel<dim>::update_quadrature_point_history()
 *   {
 * @endcode
 * 
 * First, set up an <code>FEValues</code> object by which we will evaluate
 * the incremental displacements and the gradients thereof at the
 * quadrature points, together with a vector that will hold this
 * information:
 * 
 * @code
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients);
 * 
 *     std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads(
 *       quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim));
 * 
 * @endcode
 * 
 * Then loop over all cells and do the job in the cells that belong to our
 * subdomain:
 * 
 * @code
 *     for (auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 * @endcode
 * 
 * Next, get a pointer to the quadrature point history data local to
 * the present cell, and, as a defensive measure, make sure that
 * this pointer is within the bounds of the global array:
 * 
 * @code
 *           PointHistory<dim> *local_quadrature_points_history =
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
 *           Assert(local_quadrature_points_history >=
 *                    &quadrature_point_history.front(),
 *                  ExcInternalError());
 *           Assert(local_quadrature_points_history <=
 *                    &quadrature_point_history.back(),
 *                  ExcInternalError());
 * 
 * @endcode
 * 
 * Then initialize the <code>FEValues</code> object on the present
 * cell, and extract the gradients of the displacement at the
 * quadrature points for later computation of the strains
 * 
 * @code
 *           fe_values.reinit(cell);
 *           fe_values.get_function_gradients(incremental_displacement,
 *                                            displacement_increment_grads);
 * 
 * @endcode
 * 
 * Then loop over the quadrature points of this cell:
 * 
 * @code
 *           for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
 *             {
 * @endcode
 * 
 * On each quadrature point, compute the strain increment from
 * the gradients, and multiply it by the stress-strain tensor to
 * get the stress update. Then add this update to the already
 * existing strain at this point:
 * 
 * @code
 *               const SymmetricTensor<2, dim> new_stress =
 *                 (local_quadrature_points_history[q].old_stress +
 *                  (stress_strain_tensor *
 *                   get_strain(displacement_increment_grads[q])));
 * 
 * @endcode
 * 
 * Finally, we have to rotate the result. For this, we first
 * have to compute a rotation matrix at the present quadrature
 * point from the incremental displacements. In fact, it can be
 * computed from the gradients, and we already have a function
 * for that purpose:
 * 
 * @code
 *               const Tensor<2, dim> rotation =
 *                 get_rotation_matrix(displacement_increment_grads[q]);
 * @endcode
 * 
 * Note that the result, a rotation matrix, is in general an
 * antisymmetric tensor of rank 2, so we must store it as a full
 * tensor.
 * 

 * 
 * With this rotation matrix, we can compute the rotated tensor
 * by contraction from the left and right, after we expand the
 * symmetric tensor <code>new_stress</code> into a full tensor:
 * 
 * @code
 *               const SymmetricTensor<2, dim> rotated_new_stress =
 *                 symmetrize(transpose(rotation) *
 *                            static_cast<Tensor<2, dim>>(new_stress) * rotation);
 * @endcode
 * 
 * Note that while the result of the multiplication of these
 * three matrices should be symmetric, it is not due to floating
 * point round off: we get an asymmetry on the order of 1e-16 of
 * the off-diagonal elements of the result. When assigning the
 * result to a <code>SymmetricTensor</code>, the constructor of
 * that class checks the symmetry and realizes that it isn't
 * exactly symmetric; it will then raise an exception. To avoid
 * that, we explicitly symmetrize the result to make it exactly
 * symmetric.
 * 

 * 
 * The result of all these operations is then written back into
 * the original place:
 * 
 * @code
 *               local_quadrature_points_history[q].old_stress =
 *                 rotated_new_stress;
 *             }
 *         }
 *   }
 * 
 * @endcode
 * 
 * This ends the project specific namespace <code>Step18</code>. The rest is
 * as usual and as already shown in step-17: A <code>main()</code> function
 * that initializes and terminates PETSc, calls the classes that do the
 * actual work, and makes sure that we catch all exceptions that propagate
 * up to this point:
 * 
 * @code
 * } // namespace Step18
 * 
 * 
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step18;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       TopLevel<3> elastic_problem;
 *       elastic_problem.run();
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
examples/step-18/doc/results.dox



<a name="Results"></a><h1>Results</h1>



如果使用调试模式，运行该程序需要很长时间；在我的i7台式机上需要大约11分钟。幸运的是，经过优化编译的版本要快得多；在同一台机器上用<tt>make release</tt>命令重新编译后，程序只需要大约1.5分钟，这个时间要合理得多。


如果运行，该程序会打印出以下输出，解释它在这段时间内做了什么。

@verbatim
\$ time make run
[ 66%] Built target \step-18
[100%] Run \step-18 with Release configuration
Timestep 1 at time 1
  Cycle 0:
    Number of active cells:       3712 (by partition: 3712)
    Number of degrees of freedom: 17226 (by partition: 17226)
    Assembling system... norm of rhs is 1.88062e+10
    Solver converged in 103 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       12812 (by partition: 12812)
    Number of degrees of freedom: 51738 (by partition: 51738)
    Assembling system... norm of rhs is 1.86145e+10
    Solver converged in 121 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 2 at time 2
    Assembling system... norm of rhs is 1.84169e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 3 at time 3
    Assembling system... norm of rhs is 1.82355e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 4 at time 4
    Assembling system... norm of rhs is 1.80728e+10
    Solver converged in 117 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 5 at time 5
    Assembling system... norm of rhs is 1.79318e+10
    Solver converged in 116 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 6 at time 6
    Assembling system... norm of rhs is 1.78171e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 7 at time 7
    Assembling system... norm of rhs is 1.7737e+10
    Solver converged in 112 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 8 at time 8
    Assembling system... norm of rhs is 1.77127e+10
    Solver converged in 111 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 9 at time 9
    Assembling system... norm of rhs is 1.78207e+10
    Solver converged in 113 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 10 at time 10
    Assembling system... norm of rhs is 1.83544e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...


[100%] Built target run
make run  176.82s user 0.15s system 198% cpu 1:28.94 total
@endverbatim

换句话说，它是在12,000个单元和大约52,000个未知数的情况下进行计算。不是很多，但对于一个耦合的三维问题来说，足以让计算机忙上一阵子。在一天结束的时候，这就是我们的输出。

@verbatim
\$ ls -l *vtu *visit


-rw-r--r-- 1 drwells users 1706059 Feb 13 19:36 solution-0010.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0010.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0010.visit


-rw-r--r-- 1 drwells users 1707907 Feb 13 19:36 solution-0009.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0009.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0009.visit


-rw-r--r-- 1 drwells users 1703771 Feb 13 19:35 solution-0008.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0008.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0008.visit


-rw-r--r-- 1 drwells users 1693671 Feb 13 19:35 solution-0007.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0007.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0007.visit


-rw-r--r-- 1 drwells users 1681847 Feb 13 19:35 solution-0006.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0006.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0006.visit


-rw-r--r-- 1 drwells users 1670115 Feb 13 19:35 solution-0005.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0005.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0005.visit


-rw-r--r-- 1 drwells users 1658559 Feb 13 19:35 solution-0004.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0004.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0004.visit


-rw-r--r-- 1 drwells users 1639983 Feb 13 19:35 solution-0003.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0003.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0003.visit


-rw-r--r-- 1 drwells users 1625851 Feb 13 19:35 solution-0002.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0002.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0002.visit


-rw-r--r-- 1 drwells users 1616035 Feb 13 19:34 solution-0001.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:34 solution-0001.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:34 solution-0001.visit
@endverbatim




如果我们用VisIt或Paraview将这些文件可视化，我们就能看到我们的强制压缩对圆柱体造成的灾难的全貌（图像中的颜色编码了材料中的应力规范）。


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0002.0000.png" alt="Time = 2" width="400"> </div> <div class="text" align="center"> Time = 2 </div> <div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0005.0000.png" alt="时间=5" width="400"> </div> <div class="text" align="center"> 时间=5 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0007.0000.png" alt="时间=7" width="400"> </div> <div class="text" align="center">时间=7 </div> </div> </div>


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0008.0000.png" alt="Time = 8" width="400"> </div> <div class="text" align="center"> 时间 = 8 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0009.0000.png" alt="时间=9" width="400"> </div> <div class="text" align="center"> 时间=9 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0010.0000.png" alt="Time = 10" width="400"> </div> <div class="text" align="center"> Time = 10 </div> </div> </div>


可以清楚地看到，当我们不断压缩圆柱体时，它开始在完全约束的底面附近弯曲，并在大约8个时间单位后，以方位对称的方式弯曲。


虽然这个结果对于对称几何和加载来说似乎是合理的，但计算是否完全收敛还有待确定。为了确定是否收敛，我们再次运行程序，在开始时再进行一次全局细化，并将时间步长减半。这在单机上会花费很长的时间，所以我们使用了一个合适的工作站，在16个处理器上并行运行。现在输出的开头看起来像这样。

@verbatim
Timestep 1 at time 0.5
  Cycle 0:
    Number of active cells:       29696 (by partition: 1808+1802+1894+1881+1870+1840+1884+1810+1876+1818+1870+1884+1854+1903+1816+1886)
    Number of degrees of freedom: 113100 (by partition: 6936+6930+7305+7116+7326+6869+7331+6786+7193+6829+7093+7162+6920+7280+6843+7181)
    Assembling system... norm of rhs is 1.10765e+10
    Solver converged in 209 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       102034 (by partition: 6387+6202+6421+6341+6408+6201+6428+6428+6385+6294+6506+6244+6417+6527+6299+6546)
    Number of degrees of freedom: 359337 (by partition: 23255+21308+24774+24019+22304+21415+22430+22184+22298+21796+22396+21592+22325+22553+21977+22711)
    Assembling system... norm of rhs is 1.35759e+10
    Solver converged in 268 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 2 at time 1
    Assembling system... norm of rhs is 1.34674e+10
    Solver converged in 267 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 3 at time 1.5
    Assembling system... norm of rhs is 1.33607e+10
    Solver converged in 265 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 4 at time 2
    Assembling system... norm of rhs is 1.32558e+10
    Solver converged in 263 iterations.
    Updating quadrature point data...
    Moving mesh...


[...]


Timestep 20 at time 10
    Assembling system... norm of rhs is 1.47755e+10
    Solver converged in 425 iterations.
    Updating quadrature point data...
    Moving mesh...
@endverbatim

考虑到我们是在三维空间中，这是一个相当好的未知数的数量。这个程序的输出是每个时间步骤的16个文件。

@verbatim
\$ ls -l solution-0001*


-rw-r--r-- 1 wellsd2 user 761065 Feb 13 21:09 solution-0001.000.vtu


-rw-r--r-- 1 wellsd2 user 759277 Feb 13 21:09 solution-0001.001.vtu


-rw-r--r-- 1 wellsd2 user 761217 Feb 13 21:09 solution-0001.002.vtu


-rw-r--r-- 1 wellsd2 user 761605 Feb 13 21:09 solution-0001.003.vtu


-rw-r--r-- 1 wellsd2 user 756917 Feb 13 21:09 solution-0001.004.vtu


-rw-r--r-- 1 wellsd2 user 752669 Feb 13 21:09 solution-0001.005.vtu


-rw-r--r-- 1 wellsd2 user 735217 Feb 13 21:09 solution-0001.006.vtu


-rw-r--r-- 1 wellsd2 user 750065 Feb 13 21:09 solution-0001.007.vtu


-rw-r--r-- 1 wellsd2 user 760273 Feb 13 21:09 solution-0001.008.vtu


-rw-r--r-- 1 wellsd2 user 777265 Feb 13 21:09 solution-0001.009.vtu


-rw-r--r-- 1 wellsd2 user 772469 Feb 13 21:09 solution-0001.010.vtu


-rw-r--r-- 1 wellsd2 user 760833 Feb 13 21:09 solution-0001.011.vtu


-rw-r--r-- 1 wellsd2 user 782241 Feb 13 21:09 solution-0001.012.vtu


-rw-r--r-- 1 wellsd2 user 748905 Feb 13 21:09 solution-0001.013.vtu


-rw-r--r-- 1 wellsd2 user 738413 Feb 13 21:09 solution-0001.014.vtu


-rw-r--r-- 1 wellsd2 user 762133 Feb 13 21:09 solution-0001.015.vtu


-rw-r--r-- 1 wellsd2 user   1421 Feb 13 21:09 solution-0001.pvtu


-rw-r--r-- 1 wellsd2 user    364 Feb 13 21:09 solution-0001.visit
@endverbatim




这里首先是我们计算的网格，以及16个处理器的分区。


<div class="twocolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-000mesh.png" alt="Discretization" width="400"> </div> <div class="text" align="center"> Discretization </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.p.png" alt="Parallel partitioning" width="400"> </div> <div class="text" align="center"> Parallel partitioning</div> </div> </div>


最后，这里是与我们之前展示的更小的顺序情况相同的输出。

<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.s.png" alt="Time = 2" width="400"> </div> <div class="text" align="center"> 时间 = 2 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0005.s.png" alt="时间=5" width="400"> </div> <div class="text" align="center"> 时间=5 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0007.s.png" alt="Time = 7" width="400"> </div> <div class="text" align="center"> Time = 7 </div> </div> </div>


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0008.s.png" alt="Time = 8" width="400"> </div> <div class="text" align="center"> 时间 = 8 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0009.s.png" alt="Time = 9" width="400"> </div> <div class="text" align="center"> Time = 9 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0010.s.png" alt="Time = 10" width="400"> </div> <div class="text" align="center"> Time = 10 </div> </div> </div>


和以前一样，我们观察到，在高轴向压缩时，圆柱体开始弯曲，但这一次最终是在自己身上塌陷。与我们的第一次运行相反，在模拟结束时，变形模式变得不对称（中心隆起向侧面偏转）。该模型显然没有规定这一点（我们所有的力和边界偏转都是对称的），但这种效果可能在物理上是正确的：在现实中，身体材料属性的小不均匀性会导致它向一侧弯曲以逃避强制力；在数值模拟中，小的扰动，如数值舍入或迭代求解器对线性系统的不精确求解，也会产生同样的效果。在自适应计算中，另一个典型的不对称来源是每一步只细化一定的单元，这可能导致不对称的网格，即使原来的粗网格是对称的。


如果将其与之前的运行相比较，结果在质和量上都有不同。因此，以前的计算肯定没有收敛，尽管我们不能肯定地说现在的计算有什么问题。我们需要一个更精细的计算来找出答案。然而，这一点可能是没有意义的：详细看一下最后一张图片，很明显，不仅我们选择的线性小变形模型是完全不够的，而且对于一个现实的模拟，我们还需要确保身体在变形过程中不相交（如果我们继续压缩圆柱体，我们会观察到一些自我相交）。如果没有这样的表述，我们就不能指望任何东西都有物理意义，即使它能产生漂亮的图片!




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序并没有真正解决一个在实践中有很多应用的方程：基于纯弹性规律的准静态材料变形几乎是无聊的。然而，该程序可以作为更有趣的实验的起点，而这确实是编写该程序的最初动机。这里有一些建议，说明这个程序缺少什么，以及它可以在什么方向上进行扩展。

<a name="Plasticitymodels"></a><h5>Plasticity models</h5>


最明显的扩展是使用一个更现实的材料模型来处理大规模的静态变形。这方面的自然选择是塑性，其中应力和应变之间的非线性关系取代了方程<a href="#step_18.stress-strain">[stress-strain]</a>。塑性模型的编程通常相当复杂，因为应力-应变关系通常是非平滑的。可以认为材料只能承受一个最大的应力（屈服应力），之后它就会开始&ldquo;流动&rdquo;。这方面的数学描述可以以变分不等式的形式给出，也可以将其视为弹性能量的最小化

@f[
  E(\mathbf{u}) =
  (\varepsilon(\mathbf{u}), C\varepsilon(\mathbf{u}))_{\Omega}


  - (\mathbf{f}, \mathbf{u})_{\Omega} - (\mathbf{b}, \mathbf{u})_{\Gamma_N},


@f]

受制于约束条件

@f[
  f(\sigma(\mathbf{u})) \le 0


@f]

对应力的影响。这种扩展使得在每个时间步长中要解决的问题是非线性的，所以我们需要在每个时间步长中的另一个循环。

在不进一步了解这个模型的细节的情况下，我们可以参考Simo和Hughes关于&ldquo;计算非弹性&rdquo;的优秀书籍，以全面了解解决塑性模型的计算策略。另外，在S. Commend, A. Truty, and Th. Zimmermann的文章中，对塑性的算法做了简单而简洁的描述。Zimmermann;  @cite CTZ04  。




<a name="Stabilizationissues"></a><h5>Stabilization issues</h5>


我们选择的公式，即对位移矢量的所有分量使用分片（双，三）线性元素，并将应力视为依赖于位移的变量，对于大多数材料是合适的。然而，对于不可压缩或几乎不可压缩的材料，这种所谓的基于位移的公式变得不稳定，并表现出虚假的模式。虽然流体通常不是弹性的（在大多数情况下，应力取决于速度梯度，而不是位移梯度，但也有例外，如电流变流体），但也有少数固体是几乎不可压缩的，如橡胶。另一种情况是，许多塑性模型最终让材料变得不可压缩，尽管这不在本方案的范围之内。

不可压缩性是由泊松比来表征的

@f[
  \nu = \frac{\lambda}{2(\lambda+\mu)},


@f]

其中 $\lambda,\mu$ 是材料的Lam&eacute; 常数。物理约束表明 $-1\le \nu\le \frac 12$ （该条件也来自于数学稳定性考虑）。如果 $\nu$ 接近 $\frac 12$ ，则材料变得不可压缩。在这种情况下，纯粹的基于位移的公式不再适合于解决这类问题，必须采用稳定化技术以获得稳定和准确的解决方案。上面引用的书和论文给出了如何做到这一点的指示，但在这个问题上也有大量的文献；在H.-Y. Duan和Q. Lin的论文的参考文献中可以找到一个获得该主题概述的良好开端。H.-Y. Duan and Q. Lin;  @cite DL05  。




<a name="Refinementduringtimesteps"></a><h5>Refinement during timesteps</h5>


在目前的形式下，程序只对初始网格进行若干次细化，然后就不再进行细化。对于任何一种现实的模拟，我们都希望将其扩展到每隔几步就对网格进行细化和粗化。事实上，这并不难做到，但如果你愿意的话，可以留待将来的教程程序或作为练习。

我们必须克服的主要复杂问题是，我们必须将存储在旧网格单元的正交点中的数据转移到新网格中，最好是通过某种投影方案。这方面的一般方法是这样的。

- 开始时，数据只在各个单元的正交点上可用，而不是作为一个到处定义的有限元场。

- 所以让我们找到一个<i>is</i>处处定义的有限元场，这样我们以后就可以把它插到新网格的正交点上。一般来说，要找到一个与正交点中的数值完全匹配的连续有限元场是很困难的，因为这些场的自由度数与正交点的数量不匹配，这个全局场的节点值要么是过定的，要么是欠定的。但是找到一个与正交点数值相匹配的不连续场通常不是很困难；例如，如果你有一个QGauss(2)正交公式（即2d中每个单元4个点，3d中8个点），那么就可以使用FE_DGQ(1)类型的有限元，即双/三线性函数，因为这些函数在2d中每个单元有4个自由度，在3d中有8个自由度。

- 有一些函数可以使这种从单个点到全局场的转换更简单。如果你使用QGauss(2)正交公式，下面这段伪代码应该会有所帮助。请注意，下面的投影矩阵的乘法需要一个标量分量的向量，也就是说，我们一次只能将一组标量从正交点转换成自由度，反之亦然。所以我们需要分别存储每个应力分量，这需要 <code>dim*dim</code> 个向量。我们将把这组向量存储在一个二维数组中，以便于用读出应力张量的方式来读出分量。   因此，我们将对每个单元的应力分量进行循环，并将这些值存储在全局历史域中。(前缀 <code>history_</code> 表示我们的工作与正交点中定义的历史变量有关。)   @code
    FE_DGQ<dim>     history_fe (1);
    DoFHandler<dim> history_dof_handler (triangulation);
    history_dof_handler.distribute_dofs (history_fe);


    std::vector< std::vector< Vector<double> > >
                 history_field (dim, std::vector< Vector<double> >(dim)),
                 local_history_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_history_fe_values (dim, std::vector< Vector<double> >(dim));


    for (unsigned int i=0; i<dim; i++)
      for (unsigned int j=0; j<dim; j++)
      {
        history_field[i][j].reinit(history_dof_handler.n_dofs());
        local_history_values_at_qpoints[i][j].reinit(quadrature.size());
        local_history_fe_values[i][j].reinit(history_fe.n_dofs_per_cell());
      }


    FullMatrix<double> qpoint_to_dof_matrix (history_fe.dofs_per_cell,
                                             quadrature.size());
    FETools::compute_projection_from_quadrature_points_matrix
              (history_fe,
               quadrature, quadrature,
               qpoint_to_dof_matrix);


    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();


    for (; cell!=endc; ++cell, ++dg_cell)
      {


        PointHistory<dim> *local_quadrature_points_history
          = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());


        Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
                ExcInternalError());
        Assert (local_quadrature_points_history < &quadrature_point_history.back(),
                ExcInternalError());


        for (unsigned int i=0; i<dim; i++)
          for (unsigned int j=0; j<dim; j++)
          {
            for (unsigned int q=0; q<quadrature.size(); ++q)
              local_history_values_at_qpoints[i][j](q)
                = local_quadrature_points_history[q].old_stress[i][j];


            qpoint_to_dof_matrix.vmult (local_history_fe_values[i][j],
                                        local_history_values_at_qpoints[i][j]);


            dg_cell->set_dof_values (local_history_fe_values[i][j],
                                     history_field[i][j]);
          }
      }
  @endcode



- 现在我们有了一个全局场，我们可以像往常一样使用SolutionTransfer类来细化网格并转移history_field向量。这将把所有的东西从旧的网格插值到新的网格。

- 在最后一步，我们必须将数据从现在插值的全局场返回到新网格上的正交点。下面的代码将做到这一点。   @code
    FullMatrix<double> dof_to_qpoint_matrix (quadrature.size(),
                                             history_fe.n_dofs_per_cell());
    FETools::compute_interpolation_to_quadrature_points_matrix
              (history_fe,
               quadrature,
               dof_to_qpoint_matrix);


    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();


    for (; cell != endc; ++cell, ++dg_cell)
    {
      PointHistory<dim> *local_quadrature_points_history
       = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());


      Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
              ExcInternalError());
      Assert (local_quadrature_points_history < &quadrature_point_history.back(),
              ExcInternalError());


      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
        {
          dg_cell->get_dof_values (history_field[i][j],
                                   local_history_fe_values[i][j]);


          dof_to_qpoint_matrix.vmult (local_history_values_at_qpoints[i][j],
                                      local_history_fe_values[i][j]);


          for (unsigned int q=0; q<quadrature.size(); ++q)
            local_quadrature_points_history[q].old_stress[i][j]
              = local_history_values_at_qpoints[i][j](q);
      }
  @endcode



一旦我们并行运行程序，情况就变得有点复杂了，因为那时每个进程只为它在旧网格上拥有的单元存储这些数据。也就是说，如果你在正交点转移到全局向量之后，使用 <code>history_field</code> 的并行向量就可以做到这一点。




<a name="Ensuringmeshregularity"></a><h5>Ensuring mesh regularity</h5>


目前，程序没有尝试确保一个单元在时间步数结束时移动其顶点后，仍然具有有效的几何形状（即它的雅各布行列式是正的，并且在任何地方都远离零的界限）。事实上，设置边界值和强迫项并不难，这样就可以很快得到扭曲和倒置的单元。当然，在某些大变形的情况下，这在有限网格的情况下是不可避免的，但在其他一些情况下，通过适当的网格细化和/或减少时间步长，这应该是可以避免的。这个程序没有做到这一点，但是一个更复杂的版本肯定应该采用某种启发式方法来定义哪些单元的变形量是可以接受的，哪些是不可以的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-18.cc"
*/
