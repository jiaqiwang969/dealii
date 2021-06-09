/**
@page step_61 The step-61 tutorial program
This tutorial depends on step-51.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#WeakGalerkinfiniteelementmethods"> Weak Galerkin finite element methods </a>
        <li><a href="#Theequationtosolve"> The equation to solve </a>
        <li><a href="#WeakGalerkinscheme"> Weak Galerkin scheme </a>
        <li><a href="#Representingtheweakgradient"> Representing the weak gradient </a>
        <li><a href="#Assemblingthelinearsystem"> Assembling the linear system </a>
        <li><a href="#PostprocessingandiLsub2subiiLsub2subierrors"> Post-processing and <i>L<sub>2</sub></i><i>L<sub>2</sub></i>-errors </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#TheWGDarcyEquationclasstemplate">The WGDarcyEquation class template</a>
        <li><a href="#Righthandsideboundaryvaluesandexactsolution">Right hand side, boundary values, and exact solution</a>
        <li><a href="#WGDarcyEquationclassimplementation">WGDarcyEquation class implementation</a>
      <ul>
        <li><a href="#WGDarcyEquationWGDarcyEquation">WGDarcyEquation::WGDarcyEquation</a>
        <li><a href="#WGDarcyEquationmake_grid">WGDarcyEquation::make_grid</a>
        <li><a href="#WGDarcyEquationsetup_system">WGDarcyEquation::setup_system</a>
        <li><a href="#WGDarcyEquationassemble_system">WGDarcyEquation::assemble_system</a>
        <li><a href="#WGDarcyEquationdimsolve">WGDarcyEquation<dim>::solve</a>
        <li><a href="#WGDarcyEquationdimcompute_postprocessed_velocity">WGDarcyEquation<dim>::compute_postprocessed_velocity</a>
        <li><a href="#WGDarcyEquationdimcompute_pressure_error">WGDarcyEquation<dim>::compute_pressure_error</a>
        <li><a href="#WGDarcyEquationdimcompute_velocity_error">WGDarcyEquation<dim>::compute_velocity_error</a>
        <li><a href="#WGDarcyEquationoutput_results">WGDarcyEquation::output_results</a>
        <li><a href="#WGDarcyEquationrun">WGDarcyEquation::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestresultsoniWGQsub0subQsub0subRTsub0subiiWGQsub0subQsub0subRTsub0subi">Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i><i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik0iik0i">Convergence table for <i>k=0</i><i>k=0</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub1subQsub1subRTsub1subiiWGQsub1subQsub1subRTsub1subi">Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i><i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik1iik1i">Convergence table for <i>k=1</i><i>k=1</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub2subQsub2subRTsub2subiiWGQsub2subQsub2subRTsub2subi">Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i><i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik2iik2i">Convergence table for <i>k=2</i><i>k=2</i></a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-61/doc/intro.dox

 <br> 

<i>
This program was contributed by Zhuoran Wang.
Some more information about this program, as well as more numerical
results, are presented in @cite Wang2019 .
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序介绍了泊松方程的 "弱加勒金 "有限元方法的实现。从某种意义上说，考虑这种方法的动机与步骤51中的动机相同：我们想考虑不连续的形状函数，但又需要解决这样一个事实：与通常的连续Galerkin方法相比，所产生的问题有更多的自由度（因为，例如，每个顶点携带的自由度与相邻单元一样多）。我们还必须解决这样一个事实：与连续Galerkin方法不同，<i>every</i>一个单元上的自由度与它的每个面邻单元上的所有自由度相耦合。因此，从 "传统的 "非连续Galerkin方法得到的矩阵既大又相对密集。

step-51中的混合非连续Galerkin方法（HDG）和本教程中的弱Galerkin（WG）方法都是通过引入额外的自由度来解决耦合问题的，这些自由度的形状函数只存在于单元间的一个面上（即网格的 "骨架 "上），因此它们将相邻单元上的自由度相互 "隔离"：单元自由度只与同一单元上的其他单元自由度以及面自由度耦合，而与相邻单元上的单元自由度不耦合。因此，这些细胞自由度的形状函数的耦合确实正好耦合在一个细胞和定义在其面上的自由度上。

对于一个给定的方程，例如二阶泊松方程，HDG和WG方法的区别在于如何精确地制定连接所有这些不同形状函数的问题。事实上，对于某些WG和HDG的表述，有可能表明它们是等价的）。HDG的做法是用一阶方程系统重新表述二阶问题，然后在概念上把面的自由度看作是这个一阶系统的 "通量"。相比之下，WG方法保持二阶形式，并将面的自由度视为与主解变量相同的类型，只是限制在低维的面。为了方程的目的，在定义对其应用微分算子的含义时，人们需要以某种方式将这些形状函数 "扩展 "到单元的内部。与HDG相比，该方法的优势在于它不会因为将方程重写为一阶系统而导致未知数的增加，但它也不太容易实现。然而，正如我们在下文中所看到的，这种额外的努力并不可怕。




<a name="WeakGalerkinfiniteelementmethods"></a><h3> Weak Galerkin finite element methods </h3>


弱加勒金有限元方法（WGFEMs）使用离散的弱函数来近似标量未知数，使用离散的弱梯度来近似经典梯度。该方法最初是由王俊平和叶秀在论文<a href="https://doi.org/10.1016/j.cam.2012.10.003">
<i>A weak Galerkin finite element method for second order elliptic problems</i><i>A weak Galerkin finite element method for second order elliptic problems</i>,
J. Comput. Appl. Math., 103-115, 2013</a>中提出。与连续Galerkin方法相比，弱Galerkin方法满足重要的物理特性，即局部质量守恒和体法通量连续。它的结果是一个SPD线性系统，并且通过网格细化可以获得最佳收敛率。




<a name="Theequationtosolve"></a><h3> The equation to solve </h3> 该程序使用弱加尔金有限元法求解泊松方程。


@f{align*}{
  \nabla \cdot \left( -\mathbf{K} \nabla p \right)
    &= f,
    \qquad \mathbf{x} \in \Omega, \\
  p &=  p_D,\qquad \mathbf{x} \in \Gamma^D, \\
  \mathbf{u} \cdot \mathbf{n} &= u_N,
  \qquad \mathbf{x} \in \Gamma^N,


@f}

其中 $\Omega \subset \mathbb{R}^n (n=2,3)$ 是一个有界域。在流体流经多孔介质的背景下， $p$ 是压力， $\mathbf{K}$ 是渗透性张量， $f$ 是源项， $p_D, u_N$ 代表Dirichlet和Neumann边界条件。我们可以引入一个通量， $\mathbf{u} = -\mathbf{K} \nabla p$ ，对应于达西速度（以我们在步骤20中的方式），这个变量在下面的考虑中很重要。

在这个程序中，我们将考虑一个测试案例，即在单位平方域上的确切压力为 $p = \sin \left( \pi x\right)\sin\left(\pi y \right)$ ，具有同质Dirichelet边界条件和 $\mathbf{K}$ 身份矩阵。然后我们将计算压力、速度和通量的 $L_2$ 误差。




<a name="WeakGalerkinscheme"></a><h3> Weak Galerkin scheme </h3>


上面的泊松方程有一个解 $p$ ，需要满足问题的弱表述。

@f{equation*}
\mathcal{A}\left(p,q \right) = \mathcal{F} \left(q \right),


@f}

为所有测试函数  $q$  ，其中

@f{equation*}
\mathcal{A}\left(p,q\right)
  \dealcoloneq \int_\Omega \left(\mathbf{K} \nabla p\right) \cdot \nabla q \;\mathrm{d}x,


@f}

和

@f{equation*}
\mathcal{F}\left(q\right)
  \dealcoloneq \int_\Omega f \, q \;\mathrm{d}x


  - \int_{\Gamma^N} u_N q \; \mathrm{d}x.


@f}

在这里，我们以双线性形式进行了部分积分，我们在内部评估 $p,p$ 的梯度，在域的边界评估 $q$ 的值。所有这些都是很好的定义，因为我们假设解是在 $H^1$ 中，对它来说，取梯度和评估边界值是有效的操作。

弱Galerkin方法的想法是用一个<i>discontinuous function</i> $p_h$ 来近似精确的 $p$ 解。这个函数可能只在单元格之间的界面上不连续，由于我们也想沿着界面评估这个函数，我们不仅要规定它在单元格内部应该有什么值，还要规定它在界面上的值。我们通过说 $p_h$ 实际上是一个元组， $p_h=(p^\circ,p^\partial)$ ，尽管它实际上只是一个单一的函数，它要么等于 $p^\circ(x)$ ，要么等于 $p^\partial(x)$ ，这取决于它是在位于细胞内部还是在细胞界面的某一点 $x$ 上被评估。

然后我们想把这个近似值简单地贴到上面的双线性表格中。这适用于我们必须在边界上评估测试函数 $q_h$ 的情况（我们只需取其界面部分 $q_h^\partial$ ），但我们必须小心处理梯度，因为它只在单元格内部定义。因此，泊松方程的弱Galerkin方案被定义为

@f{equation*}
\mathcal{A}_h\left(p_h,q \right) = \mathcal{F} \left(q_h \right),


@f}

对于所有离散测试函数  $q_h$  ，其中

@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,


@f}

和

@f{equation*}
\mathcal{F}\left(q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}} \int_K f \, q_h^\circ \;\mathrm{d}x


  - \sum_{\gamma \in \Gamma_h^N} \int_\gamma u_N q_h^\partial \;\mathrm{d}x,


@f}

关键的一点是，在这里，我们用<i>discrete weak gradient</i>算子 $\nabla_{w,d} p_h$ 代替了梯度 $\nabla p_h$ ，这对于我们特殊定义的近似 $p_h$ 是有意义的。

那么问题是该算子如何工作。为此，让我们首先说说我们是如何看待压力的离散近似值 $p_h$ 的。如上所述，"函数" $p_h$ 实际上由两部分组成：单元内部的值 $p_h^\circ$ 和界面上的 $p_h^\partial$ 。我们必须为这两部分定义离散的（有限维）函数空间；在这个程序中，我们将用FE_DGQ来表示 $p_h^\circ$ 作为细胞内部的空间（在每个细胞上定义，但一般沿界面是不连续的），用FE_FaceQ表示 $p_h^\partial$ 作为界面上的空间。

那么让我们只考虑一个单元（因为上面的积分都是逐个单元定义的，而且弱离散梯度是逐个单元定义的）。 $p_h$ 对 $K$ , $p_h|_K$ 的限制由一对 $(p_h^\circ|_K,p_h^\partial|_{\partial K})$ 组成。从本质上讲，我们可以认为 $\nabla_{w,d} p_h$ 是定义在 $K$ 上的某个函数，它近似于梯度；特别是，如果 $p_h|_K$ 是一个可微函数的限制（对 $K$ 的内部和边界--这将使它在内部和边界之间连续），那么 $\nabla_{w,d} p_h$  将只是精确梯度 $\nabla p_h$  。但是，由于 $p_h|_K$ 在 $K$ 的内部和边界之间不连续，我们需要一个更一般的定义；此外，我们不能处理任意函数，因此要求 $\nabla_{w,d} p_h$ 也在一个有限元空间中（由于梯度是一个矢量，必须是矢量值，而且由于弱梯度是在每个单元上单独定义的，因此在单元之间也将是不连续的）。

这样做的方法是以下列方式定义这个弱梯度算子 $\nabla_{w,d}|_K :
DGQ_k(K) \times DGQ_r(\partial K) \rightarrow RT_s(K)$ （其中 $RT_s(K)$ 是单元格 $K$ 上阶为 $s$ 的矢量值Raviart-Thomas空间）。

@f{equation*}{
  \int_K \mathbf v_h \cdot (\nabla_{w,d} p_h)
  =


  -\int_K (\nabla \cdot \mathbf v_h) p_h^\circ
  +\int_{\partial K} (\mathbf v_h \cdot \mathbf n) p_h^\partial,


@f}

为所有测试函数  $\mathbf v_h \in RT_s(K)$  。从本质上讲，这只是一个逐部积分公式的应用。换句话说，对于一个给定的 $p_h=(p^\circ_h,p^\partial_h)$ ，我们需要把 $\nabla_{w,d} p_h|_K$ 看作是度数为 $s$ 的Raviart-Thomas函数，对于这个函数，左手边和右手边在所有测试函数中是相等的。

那么，需要说明的一个关键点是以下几点。通常的梯度 $\nabla$ 是一个*本地*算子，它仅仅根据一个函数在某一点及其（无限小）邻域的值来计算导数，而弱离散梯度 $\nabla_{w,d}$ 却没有这个特性。它取决于它所应用的函数在整个单元上的值，包括单元的边界。然而，两者都是线性算子，从上面 $\nabla_{w,d}$ 的定义可以看出，这将允许我们在下面的讨论中通过矩阵来表示 $\nabla_{w,d}$ 。

 @note  值得指出的是，虽然弱的离散梯度是Raviart-Thomas空间 $RT_s(K)$ 在每个单元 $K$ 的一个元素，但它在单元之间是不连续的。另一方面，定义在整个网格上并由FE_RaviartThomas类实现的Raviart-Thomas空间 $RT_s=RT_s({\mathbb T})$ 代表在单元间界面上具有连续法线分量的函数。这意味着<i>globally</i>,  $\nabla_{w,d} p_h$ 不在 $RT_s$ 中，尽管它在 $K$ 中的每个单元上。   相反，它是在一个 "破碎的 "拉维-托马斯空间中，下面我们将用符号 $DGRT_s$ 来表示。 这里的术语 "破碎 "指的是 "把东西打碎 "的过程，而不是表达 "没有功能 "的同义词。因此，人们可能会（理所当然地）争辩说，在弱加尔金文献中使用的符号有点误导，但这往往取决于使用某种符号的背景--在目前的背景下，对Raviart-Thomas空间或元素的提及总是被理解为对 "破碎 "空间的提及。

 @note  deal.II恰好有一个实现了这个破碎的Raviart-Thomas空间。FE_DGRT类。因此，在本教程中，我们将简单地一直使用FE_DGRT类，尽管在所有那些我们必须计算单元格本地矩阵和向量的地方，它没有任何区别。




<a name="Representingtheweakgradient"></a><h3> Representing the weak gradient </h3>


由于 $p_h$ 是有限元空间的一个元素，我们可以像往常一样在一个基础上展开它，也就是说，我们可以写出

@f{equation*}{
  p_h(\mathbf x) = \sum_j P_j \varphi_j(\mathbf x).


@f}

这里，由于 $p_h$ 有两个分量（内部分量和界面分量），对于基函数 $\varphi_j(\mathbf x)$ 也必须如此，我们可以写成 $\varphi_j = (\varphi_j^\circ,\varphi_j^\partial)$  。如果你按照步骤8、步骤20和 @ref vector_valued "向量值问题文件模块 "中的描述，就不会感到奇怪，对于 $j$ 的某些值， $\varphi_j^\circ$ 将为零，而对于 $j$ 的其他值， $\varphi_j^\partial$ 将为零--也就是说，形状函数将是一种或另一种类型。然而，这在这里并不重要。重要的是，我们需要思考如何表示 $\nabla_{w,d} \varphi_j$ ，因为当我们想实现双线性形式时，这显然是问题中会出现的东西

@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  = \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,


@f}



关键的一点是，已知 $\nabla_{w,d} \varphi_j$ 是 "破碎的 "Raviart-Thomas空间 $DGRT_s$ 的一个成员。这意味着我们可以（在每个单元 $K$ 上分别表示

@f{equation*}
\nabla_{w,d} \varphi_j|_K
  = \sum_k C_{jk}^K \mathbf v_k|_K


@f}

其中，函数 $\mathbf v_k \in DGRT_s$ ，以及 $C^K$ 是一个维数的矩阵

@f{align*}{
 \text{dim}\left(DGQ_k(K) \times DGQ_r(K)\right) &\times \text{dim}\left(RT_s(K)\right)
  \\
 &=
 \left(\text{dim}(DGQ_k(K)) + \text{dim}(DGQ_r(K))\right) \times \text{dim}\left(RT_s(K)\right).


@f}

弱离散梯度可以被表示为一个矩阵，这不应该是一个惊喜：它是一个从一个有限维空间到另一个有限维空间的线性算子。如果为这两个空间都选择基数，那么<i>every linear operator</i>当然可以写成一个矩阵，将与算子的域空间的基数有关的扩展系数向量映射到与图像空间的基数有关的扩展系数向量）。)

利用这个扩展，我们可以很容易地使用上面的弱离散梯度的定义来定义矩阵要做什么。

@f{equation*}{
  \int_K \mathbf v_i \cdot \left(\sum_k C_{jk}^K \mathbf v_k\right)
  =


  -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
  +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,


@f}

对于所有的测试功能  $\mathbf v_i \in DGRT_s$  。

这显然导致了一个线性系统，其形式为

@f{equation*}{
  \sum_k M_{ik}^K C_{jk}^K
  =
  G_{ij}^K


@f}

与

@f{equation*}{
  M_{ik}^K = \int_K \mathbf v_i \cdot \mathbf v_k,
  \qquad\qquad
  G_{ij}^K = -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
             +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,


@f}

因此

@f{equation*}{
  \left(C^K\right)^T = \left(M^K\right)^{-1} G^K.


@f}

(在这最后一步中，我们假设指数 $i,j,k$ 只涉及在单元 $K$ 上活动的自由度，从而确保空间 $RT_s(K)$ 上的质量矩阵是可逆的。)等价地，利用矩阵 $M$ 的对称性，我们可以看到

@f{equation*}{
  C^K = \left(G^K\right)^{T} \left(M^K\right)^{-1}.


@f}

另外值得指出的是，矩阵 $C^K$ 和 $G^K$ 当然不是正方形而是长方形。




<a name="Assemblingthelinearsystem"></a><h3> Assembling the linear system </h3>


在解释了弱离散梯度是如何定义的之后，我们现在可以回到有关方程的线性系统应该如何组装的问题上。具体来说，利用上面显示的双线性形式 ${\cal A}_h$ 的定义，我们就需要计算局部对全局矩阵的贡献元素。

@f{equation*}{
  A^K_{ij} = \int_K \left({\mathbf K} \nabla_{w,d} \varphi_i\right) \cdot \nabla_{w,d} \varphi_j.


@f}

如上所述，我们可以用Raviart-Thomas基础在每个单元格上展开 $\nabla_{w,d} \varphi_i$ ，同样，对于 $\nabla_{w,d} \varphi_j$ 也是如此。

@f{equation*}{
  A^K_{ij} = \int_K
    \left(
      {\mathbf K}
      \sum_k C_{ik}^K \mathbf v_k|_K
    \right)
    \cdot
    \sum_l C_{jl}^K \mathbf v_l|_K.


@f}

通过重新排列和，可以得到以下表达式。

@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
     \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K.


@f}

因此，如果我们有每个单元格 $K$ 的矩阵 $C^K$ ，那么我们可以很容易地计算出单元格 $K$ 对矩阵 $A$ 的贡献 $A^K$ ，如下所示。

@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
    H^K_{kl}
    =
    \sum_k \sum_l C_{ik}^K H^K_{kl} C_{jl}^K
    =
    \left(C^K H^K (C^K)^T \right)_{ij}.


@f}

在这里。

@f{equation*}{
  H^K_{kl} =
  \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K,


@f}

这实际上只是单元 $K$ 上的质量矩阵，使用Raviart-Thomas基础并通过渗透性张量 $\mathbf K$ 加权。这里的推导表明，弱加尔金法实际上只需要我们计算每个单元 $C^K$ 和 $H^K$ 的矩阵，然后再计算 $A^K = C^K H^K (C^K)^T$ ，这很容易计算出来。下面要显示的代码正是这样做的。

在计算出单元格 $A^K$ 对全局矩阵的贡献后，我们要做的就是将这些局部贡献 "分配 "到全局矩阵中。如何做到这一点，首先显示在步骤3和步骤4中。在目前的程序中，这将通过调用 AffineConstraints::distribute_local_to_global(). 来促进。

一个线性系统当然也需要一个右手边。除了我们只需要对每个形状函数 $\varphi_i^\circ$ 使用单元格内部部分外，这里没有与计算右手边有关的困难。




<a name="PostprocessingandiLsub2subiiLsub2subierrors"></a><h3> Post-processing and <i>L<sub>2</sub></i><i>L<sub>2</sub></i>-errors </h3> 。


前面几节的讨论已经给了我们一个线性系统，我们可以求解数值压力 $p_h$  。我们可以用它来计算变量 $\mathbf u = -{\mathbf K}\nabla p$ 的近似值，如果这是我们要解决的模型，它对应于介质在多孔介质中的流动速度。这种步骤--从离散问题的解中计算一个派生量--通常被称为 "后处理"。

这里，我们不使用 $p_h$ 的精确梯度，而是使用 $p_h$ 的离散弱梯度来计算每个元素上的速度。如上所述，在每个元素上，数值压力 $\nabla p$ 的梯度可以用离散弱梯度 $ \nabla_{w,d}\phi_i$ 来近似。

@f{equation*}
\nabla_{w,d} p_h
= \nabla_{w,d} \left(\sum_{i} P_i \phi_i\right)
= \sum_{i} P_i \nabla_{w,d}\phi_i.


@f}



在单元格 $K$ 上，数值速度 $ \mathbf{u}_h = -\mathbf{K} \nabla_{w,d}p_h$ 可写为

@f{align*}{
  \mathbf{u}_h
  &= -\mathbf{K} \nabla_{w,d} p_h
   = -\mathbf{K}\sum_{i} \sum_{j} P_i C^K_{ij}\mathbf{v}_j,


@f}

其中 $C^K$ 是上面的扩展矩阵， $\mathbf{v}_j$ 是 $RT$ 空间在一个单元上的基函数。

不幸的是， $\mathbf{K} \mathbf{v}_j$ 可能不在 $RT$ 空间中（当然，除非如果 $\mathbf K$ 是常数乘以身份矩阵）。因此，为了在有限元程序中表示它，我们需要把它投射回我们可以处理的有限维空间。在这里，我们将使用 $L_2$ 投影法将其投影回（破碎的） $RT$ 空间。

我们将每个单元格  $K$  上的投影定义为  $ \mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right) =
\sum_{k} d_{jk}\mathbf{v}_k$  。对于任何  $j$  ,  $\left( \mathbf{Q}_h \left( \mathbf{Kv}_j \right),\mathbf{v}_k \right)_K =
\left( \mathbf{Kv}_j,\mathbf{v}_k \right)_K.$  所以，与其说是上面的公式，不如说是  $K$  单元上的数字速度变成了

@f{equation*}
\mathbf{u}_h = \mathbf{Q}_h \left( -\mathbf{K}\nabla_{w,d}p_h \right) =


-\sum_i \sum_j P_i B^K_{ij}\mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right),


@f}

我们有以下系统来解决系数问题  $d_{jk}$  。

@f{equation*}
 \sum_j
  \left(\mathbf{v}_i,\mathbf{v}_j\right)
   d_{jk}
   =
    \left( \mathbf{Kv}_j,\mathbf{v}_k \right).


@f}

在下面的实现中，元素为 $
   d_{jk}
$ 的矩阵被称为 <code>cell_matrix_D</code>  ，而元素为 $
      \left( \mathbf{Kv}_j,\mathbf{v}_k \right)
$ 的矩阵被称为 <code>cell_matrix_E</code> 。

那么元素速度为

@f{equation*}
\mathbf{u}_h = -\sum_{i} \sum_{j}P_ic_{ij}\sum_{k}d_{jk}\mathbf{v}_k =
\sum_{k}- \left(\sum_{j} \sum_{i} P_ic_{ij}d_{jk} \right)\mathbf{v}_k,


@f}

其中 $-\sum_{j} \sum_{i} P_ic_{ij}d_{jk}$ 在代码中被称为`细胞速度'。

利用这个通过 "后处理 "得到的速度，我们可以通过以下公式定义压力、速度和通量的 $L_2$ 误差。

@f{align*}{
\|p-p_h^\circ\|^2
  &= \sum_{K \in \mathbb{T}} \|p-p_h^\circ\|_{L_2(K)}^2, \\
 \|\mathbf{u}-\mathbf{u}_h\|^2
  &= \sum_{K \in \mathbb{T}} \|\mathbf{u}-\mathbf{u}_h\|_{L_2(K)^2}^d,\\
\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|^2
  &= \sum_{K \in \mathbb{T}} \sum_{\gamma \subset \partial K}
    \frac{|K|}{|\gamma|} \|\mathbf{u} \cdot \mathbf{n} - \mathbf{u}_h \cdot \mathbf{n}\|_{L_2(\gamma)}^2,


@f}

其中 $| K |$ 为元素的面积， $\gamma$ 为元素的面， $\mathbf{n}$ 为每个面的单位法向量。这些规范中的最后一条衡量了网格单元之间界面上速度向量的法向分量的精度。缩放因子 $|K|/|\gamma|$ 的选择是为了随着网格大小的变化，缩放出界面集合的长度（或面积）的差异。

上面的第一个错误很容易用 VectorTools::integrate_difference. 计算出来，其他的需要多做一些工作，在下面的代码中实现。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * This program is based on step-7, step-20 and step-51,
 * so most of the following header files are familiar. We
 * need the following, of which only the one that
 * imports the FE_DGRaviartThomas class (namely, `deal.II/fe/fe_dg_vector.h`)
 * is really new; the FE_DGRaviartThomas implements the "broken" Raviart-Thomas
 * space discussed in the introduction:
 * 
 * @code
 * #include <deal.II/base/quadrature.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/tensor_function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/point.h>
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_dgq.h>
 * #include <deal.II/fe/fe_raviart_thomas.h>
 * #include <deal.II/fe/fe_dg_vector.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/fe_face.h>
 * #include <deal.II/fe/component_mask.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/data_out_faces.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * Our first step, as always, is to put everything related to this tutorial
 * program into its own namespace:
 * 
 * @code
 * namespace Step61
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="TheWGDarcyEquationclasstemplate"></a> 
 * <h3>The WGDarcyEquation class template</h3>
 * 

 * 
 * This is the main class of this program. We will solve for the numerical
 * pressure in the interior and on faces using the weak Galerkin (WG) method,
 * and calculate the $L_2$ error of pressure. In the post-processing step, we
 * will also calculate $L_2$-errors of the velocity and flux.
 *   

 * 
 * The structure of the class is not fundamentally different from that of
 * previous tutorial programs, so there is little need to comment on the
 * details with one exception: The class has a member variable `fe_dgrt`
 * that corresponds to the "broken" Raviart-Thomas space mentioned in the
 * introduction. There is a matching `dof_handler_dgrt` that represents a
 * global enumeration of a finite element field created from this element, and
 * a vector `darcy_velocity` that holds nodal values for this field. We will
 * use these three variables after solving for the pressure to compute a
 * postprocessed velocity field for which we can then evaluate the error
 * and which we can output for visualization.
 * 
 * @code
 *   template <int dim>
 *   class WGDarcyEquation
 *   {
 *   public:
 *     WGDarcyEquation(const unsigned int degree);
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void compute_postprocessed_velocity();
 *     void compute_velocity_errors();
 *     void compute_pressure_error();
 *     void output_results() const;
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     FESystem<dim>   fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     FE_DGRaviartThomas<dim> fe_dgrt;
 *     DoFHandler<dim>         dof_handler_dgrt;
 *     Vector<double>          darcy_velocity;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandexactsolution"></a> 
 * <h3>Right hand side, boundary values, and exact solution</h3>
 * 

 * 
 * Next, we define the coefficient matrix $\mathbf{K}$ (here, the
 * identity matrix), Dirichlet boundary conditions, the right-hand
 * side $f = 2\pi^2 \sin(\pi x) \sin(\pi y)$, and the exact solution
 * that corresponds to these choices for $K$ and $f$, namely $p =
 * \sin(\pi x) \sin(\pi y)$.
 * 
 * @code
 *   template <int dim>
 *   class Coefficient : public TensorFunction<2, dim>
 *   {
 *   public:
 *     Coefficient()
 *       : TensorFunction<2, dim>()
 *     {}
 * 
 *     virtual void value_list(const std::vector<Point<dim>> &points,
 *                             std::vector<Tensor<2, dim>> &values) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points,
 *                                     std::vector<Tensor<2, dim>> &  values) const
 *   {
 *     Assert(points.size() == values.size(),
 *            ExcDimensionMismatch(points.size(), values.size()));
 *     for (unsigned int p = 0; p < points.size(); ++p)
 *       values[p] = unit_symmetric_tensor<dim>();
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class BoundaryValues : public Function<dim>
 *   {
 *   public:
 *     BoundaryValues()
 *       : Function<dim>(2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/,
 *                                     const unsigned int /*component*/) const
 *   {
 *     return 0;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> &p,
 *                                    const unsigned int /*component*/) const
 *   {
 *     return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) *
 *             std::sin(numbers::PI * p[1]));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The class that implements the exact pressure solution has an
 * oddity in that we implement it as a vector-valued one with two
 * components. (We say that it has two components in the constructor
 * where we call the constructor of the base Function class.) In the
 * `value()` function, we do not test for the value of the
 * `component` argument, which implies that we return the same value
 * for both components of the vector-valued function. We do this
 * because we describe the finite element in use in this program as
 * a vector-valued system that contains the interior and the
 * interface pressures, and when we compute errors, we will want to
 * use the same pressure solution to test both of these components.
 * 
 * @code
 *   template <int dim>
 *   class ExactPressure : public Function<dim>
 *   {
 *   public:
 *     ExactPressure()
 *       : Function<dim>(2)
 *     {}
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   double ExactPressure<dim>::value(const Point<dim> &p,
 *                                    const unsigned int /*component*/) const
 *   {
 *     return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class ExactVelocity : public TensorFunction<1, dim>
 *   {
 *   public:
 *     ExactVelocity()
 *       : TensorFunction<1, dim>()
 *     {}
 * 
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override;
 *   };
 * 
 * 
 * 
 *   template <int dim>
 *   Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const
 *   {
 *     Tensor<1, dim> return_value;
 *     return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) *
 *                       std::sin(numbers::PI * p[1]);
 *     return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) *
 *                       std::cos(numbers::PI * p[1]);
 *     return return_value;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationclassimplementation"></a> 
 * <h3>WGDarcyEquation class implementation</h3>
 * 

 * 
 * 
 * <a name="WGDarcyEquationWGDarcyEquation"></a> 
 * <h4>WGDarcyEquation::WGDarcyEquation</h4>
 * 

 * 
 * In this constructor, we create a finite element space for vector valued
 * functions, which will here include the ones used for the interior and
 * interface pressures, $p^\circ$ and $p^\partial$.
 * 
 * @code
 *   template <int dim>
 *   WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree)
 *     : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1)
 *     , dof_handler(triangulation)
 *     , fe_dgrt(degree)
 *     , dof_handler_dgrt(triangulation)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationmake_grid"></a> 
 * <h4>WGDarcyEquation::make_grid</h4>
 * 

 * 
 * We generate a mesh on the unit square domain and refine it.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0, 1);
 *     triangulation.refine_global(5);
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells()
 *               << std::endl
 *               << "   Total number of cells: " << triangulation.n_cells()
 *               << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationsetup_system"></a> 
 * <h4>WGDarcyEquation::setup_system</h4>
 * 

 * 
 * After we have created the mesh above, we distribute degrees of
 * freedom and resize matrices and vectors. The only piece of
 * interest in this function is how we interpolate the boundary
 * values for the pressure. Since the pressure consists of interior
 * and interface components, we need to make sure that we only
 * interpolate onto that component of the vector-valued solution
 * space that corresponds to the interface pressures (as these are
 * the only ones that are defined on the boundary of the domain). We
 * do this via a component mask object for only the interface
 * pressures.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler_dgrt.distribute_dofs(fe_dgrt);
 * 
 *     std::cout << "   Number of pressure degrees of freedom: "
 *               << dof_handler.n_dofs() << std::endl;
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 * 
 *     {
 *       constraints.clear();
 *       const FEValuesExtractors::Scalar interface_pressure(1);
 *       const ComponentMask              interface_pressure_mask =
 *         fe.component_mask(interface_pressure);
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                BoundaryValues<dim>(),
 *                                                constraints,
 *                                                interface_pressure_mask);
 *       constraints.close();
 *     }
 * 
 * 
 * @endcode
 * 
 * In the bilinear form, there is no integration term over faces
 * between two neighboring cells, so we can just use
 * <code>DoFTools::make_sparsity_pattern</code> to calculate the sparse
 * matrix.
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationassemble_system"></a> 
 * <h4>WGDarcyEquation::assemble_system</h4>
 * 

 * 
 * This function is more interesting. As detailed in the
 * introduction, the assembly of the linear system requires us to
 * evaluate the weak gradient of the shape functions, which is an
 * element in the Raviart-Thomas space. As a consequence, we need to
 * define a Raviart-Thomas finite element object, and have FEValues
 * objects that evaluate it at quadrature points. We then need to
 * compute the matrix $C^K$ on every cell $K$, for which we need the
 * matrices $M^K$ and $G^K$ mentioned in the introduction.
 *   

 * 
 * A point that may not be obvious is that in all previous tutorial
 * programs, we have always called FEValues::reinit() with a cell
 * iterator from a DoFHandler. This is so that one can call
 * functions such as FEValuesBase::get_function_values() that
 * extract the values of a finite element function (represented by a
 * vector of DoF values) on the quadrature points of a cell. For
 * this operation to work, one needs to know which vector elements
 * correspond to the degrees of freedom on a given cell -- i.e.,
 * exactly the kind of information and operation provided by the
 * DoFHandler class.
 *   

 * 
 * We could create a DoFHandler object for the "broken" Raviart-Thomas space
 * (using the FE_DGRT class), but we really don't want to here: At
 * least in the current function, we have no need for any globally defined
 * degrees of freedom associated with this broken space, but really only
 * need to reference the shape functions of such a space on the current
 * cell. As a consequence, we use the fact that one can call
 * FEValues::reinit() also with cell iterators into Triangulation
 * objects (rather than DoFHandler objects). In this case, FEValues
 * can of course only provide us with information that only
 * references information about cells, rather than degrees of freedom
 * enumerated on these cells. So we can't use
 * FEValuesBase::get_function_values(), but we can use
 * FEValues::shape_value() to obtain the values of shape functions
 * at quadrature points on the current cell. It is this kind of
 * functionality we will make use of below. The variable that will
 * give us this information about the Raviart-Thomas functions below
 * is then the `fe_values_rt` (and corresponding `fe_face_values_rt`)
 * object.
 *   

 * 
 * Given this introduction, the following declarations should be
 * pretty obvious:
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::assemble_system()
 *   {
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim>     fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values);
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     FEValues<dim>     fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell();
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size();
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();
 * 
 *     RightHandSide<dim>  right_hand_side;
 *     std::vector<double> right_hand_side_values(n_q_points);
 * 
 *     const Coefficient<dim>      coefficient;
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * 
 * @endcode
 * 
 * Next, let us declare the various cell matrices discussed in the
 * introduction:
 * 
 * @code
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell);
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt);
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 *     Vector<double>     cell_solution(dofs_per_cell);
 * 
 * @endcode
 * 
 * We need <code>FEValuesExtractors</code> to access the @p interior and
 * @p face component of the shape functions.
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure_interior(0);
 *     const FEValuesExtractors::Scalar pressure_face(1);
 * 
 * @endcode
 * 
 * This finally gets us in position to loop over all cells. On
 * each cell, we will first calculate the various cell matrices
 * used to construct the local matrix -- as they depend on the
 * cell in question, they need to be re-computed on each cell. We
 * need shape functions for the Raviart-Thomas space as well, for
 * which we need to create first an iterator to the cell of the
 * triangulation, which we can obtain by assignment from the cell
 * pointing into the DoFHandler.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 * 
 *         const typename Triangulation<dim>::active_cell_iterator cell_dgrt =
 *           cell;
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(),
 *                                    right_hand_side_values);
 *         coefficient.value_list(fe_values.get_quadrature_points(),
 *                                coefficient_values);
 * 
 * @endcode
 * 
 * The first cell matrix we will compute is the mass matrix
 * for the Raviart-Thomas space.  Hence, we need to loop over
 * all the quadrature points for the velocity FEValues object.
 * 
 * @code
 *         cell_matrix_M = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q);
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *                 {
 *                   const Tensor<1, dim> v_k =
 *                     fe_values_dgrt[velocities].value(k, q);
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q));
 *                 }
 *             }
 * @endcode
 * 
 * Next we take the inverse of this matrix by using
 * FullMatrix::gauss_jordan(). It will be used to calculate
 * the coefficient matrix $C^K$ later. It is worth recalling
 * later that `cell_matrix_M` actually contains the *inverse*
 * of $M^K$ after this call.
 * 
 * @code
 *         cell_matrix_M.gauss_jordan();
 * 
 * @endcode
 * 
 * From the introduction, we know that the right hand side
 * $G^K$ of the equation that defines $C^K$ is the difference
 * between a face integral and a cell integral. Here, we
 * approximate the negative of the contribution in the
 * interior. Each component of this matrix is the integral of
 * a product between a basis function of the polynomial space
 * and the divergence of a basis function of the
 * Raviart-Thomas space. These basis functions are defined in
 * the interior.
 * 
 * @code
 *         cell_matrix_G = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const double div_v_i =
 *                 fe_values_dgrt[velocities].divergence(i, q);
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const double phi_j_interior =
 *                     fe_values[pressure_interior].value(j, q);
 * 
 *                   cell_matrix_G(i, j) -=
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q));
 *                 }
 *             }
 * 
 * 
 * @endcode
 * 
 * Next, we approximate the integral on faces by quadrature.
 * Each component is the integral of a product between a basis function
 * of the polynomial space and the dot product of a basis function of
 * the Raviart-Thomas space and the normal vector. So we loop over all
 * the faces of the element and obtain the normal vector.
 * 
 * @code
 *         for (const auto &face : cell->face_iterators())
 *           {
 *             fe_face_values.reinit(cell, face);
 *             fe_face_values_dgrt.reinit(cell_dgrt, face);
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 const Tensor<1, dim> &normal = fe_face_values.normal_vector(q);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *                   {
 *                     const Tensor<1, dim> v_i =
 *                       fe_face_values_dgrt[velocities].value(i, q);
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                       {
 *                         const double phi_j_face =
 *                           fe_face_values[pressure_face].value(j, q);
 * 
 *                         cell_matrix_G(i, j) +=
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
 *                       }
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * @p cell_matrix_C is then the matrix product between the
 * transpose of $G^K$ and the inverse of the mass matrix
 * (where this inverse is stored in @p cell_matrix_M):
 * 
 * @code
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);
 * 
 * @endcode
 * 
 * Finally we can compute the local matrix $A^K$.  Element
 * $A^K_{ij}$ is given by $\int_{E} \sum_{k,l} C_{ik} C_{jl}
 * (\mathbf{K} \mathbf{v}_k) \cdot \mathbf{v}_l
 * \mathrm{d}x$. We have calculated the coefficients $C$ in
 * the previous step, and so obtain the following after
 * suitably re-arranging the loops:
 * 
 * @code
 *         local_matrix = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           {
 *             for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *               {
 *                 const Tensor<1, dim> v_k =
 *                   fe_values_dgrt[velocities].value(k, q);
 *                 for (unsigned int l = 0; l < dofs_per_cell_dgrt; ++l)
 *                   {
 *                     const Tensor<1, dim> v_l =
 *                       fe_values_dgrt[velocities].value(l, q);
 * 
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                       for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                         local_matrix(i, j) +=
 *                           (coefficient_values[q] * cell_matrix_C[i][k] * v_k) *
 *                           cell_matrix_C[j][l] * v_l * fe_values_dgrt.JxW(q);
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * Next, we calculate the right hand side, $\int_{K} f q \mathrm{d}x$:
 * 
 * @code
 *         cell_rhs = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               cell_rhs(i) += (fe_values[pressure_interior].value(i, q) *
 *                               right_hand_side_values[q] * fe_values.JxW(q));
 *             }
 * 
 * @endcode
 * 
 * The last step is to distribute components of the local
 * matrix into the system matrix and transfer components of
 * the cell right hand side into the system right hand side:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(
 *           local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimsolve"></a> 
 * <h4>WGDarcyEquation<dim>::solve</h4>
 * 

 * 
 * This step is rather trivial and the same as in many previous
 * tutorial programs:
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> solver(solver_control);
 *     solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_postprocessed_velocity"></a> 
 * <h4>WGDarcyEquation<dim>::compute_postprocessed_velocity</h4>
 * 

 * 
 * In this function, compute the velocity field from the pressure
 * solution previously computed. The
 * velocity is defined as $\mathbf{u}_h = \mathbf{Q}_h \left(
 * -\mathbf{K}\nabla_{w,d}p_h \right)$, which requires us to compute
 * many of the same terms as in the assembly of the system matrix.
 * There are also the matrices $E^K,D^K$ we need to assemble (see
 * the introduction) but they really just follow the same kind of
 * pattern.
 *   

 * 
 * Computing the same matrices here as we have already done in the
 * `assemble_system()` function is of course wasteful in terms of
 * CPU time. Likewise, we copy some of the code from there to this
 * function, and this is also generally a poor idea. A better
 * implementation might provide for a function that encapsulates
 * this duplicated code. One could also think of using the classic
 * trade-off between computing efficiency and memory efficiency to
 * only compute the $C^K$ matrices once per cell during the
 * assembly, storing them somewhere on the side, and re-using them
 * here. (This is what step-51 does, for example, where the
 * `assemble_system()` function takes an argument that determines
 * whether the local matrices are recomputed, and a similar approach
 * -- maybe with storing local matrices elsewhere -- could be
 * adapted for the current program.)
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_postprocessed_velocity()
 *   {
 *     darcy_velocity.reinit(dof_handler_dgrt.n_dofs());
 * 
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_quadrature_points |
 *                               update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_normal_vectors |
 *                                        update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell();
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell();
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size();
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size();
 * 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices_dgrt(
 *       dofs_per_cell_dgrt);
 * 
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell);
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_D(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 *     FullMatrix<double> cell_matrix_E(dofs_per_cell_dgrt, dofs_per_cell_dgrt);
 * 
 *     Vector<double> cell_solution(dofs_per_cell);
 *     Vector<double> cell_velocity(dofs_per_cell_dgrt);
 * 
 *     const Coefficient<dim>      coefficient;
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points_dgrt);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 *     const FEValuesExtractors::Scalar pressure_interior(0);
 *     const FEValuesExtractors::Scalar pressure_face(1);
 * 
 * @endcode
 * 
 * In the introduction, we explained how to calculate the numerical velocity
 * on the cell. We need the pressure solution values on each cell,
 * coefficients of the Gram matrix and coefficients of the $L_2$ projection.
 * We have already calculated the global solution, so we will extract the
 * cell solution from the global solution. The coefficients of the Gram
 * matrix have been calculated when we assembled the system matrix for the
 * pressures. We will do the same way here. For the coefficients of the
 * projection, we do matrix multiplication, i.e., the inverse of the Gram
 * matrix times the matrix with $(\mathbf{K} \mathbf{w}, \mathbf{w})$ as
 * components. Then, we multiply all these coefficients and call them beta.
 * The numerical velocity is the product of beta and the basis functions of
 * the Raviart-Thomas space.
 * 
 * @code
 *     typename DoFHandler<dim>::active_cell_iterator
 *       cell = dof_handler.begin_active(),
 *       endc = dof_handler.end(), cell_dgrt = dof_handler_dgrt.begin_active();
 *     for (; cell != endc; ++cell, ++cell_dgrt)
 *       {
 *         fe_values.reinit(cell);
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 *         coefficient.value_list(fe_values_dgrt.get_quadrature_points(),
 *                                coefficient_values);
 * 
 * @endcode
 * 
 * The component of this <code>cell_matrix_E</code> is the integral of
 * $(\mathbf{K} \mathbf{w}, \mathbf{w})$. <code>cell_matrix_M</code> is
 * the Gram matrix.
 * 
 * @code
 *         cell_matrix_M = 0;
 *         cell_matrix_E = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q);
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *                 {
 *                   const Tensor<1, dim> v_k =
 *                     fe_values_dgrt[velocities].value(k, q);
 * 
 *                   cell_matrix_E(i, k) +=
 *                     (coefficient_values[q] * v_i * v_k * fe_values_dgrt.JxW(q));
 * 
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q));
 *                 }
 *             }
 * 
 * @endcode
 * 
 * To compute the matrix $D$ mentioned in the introduction, we
 * then need to evaluate $D=M^{-1}E$ as explained in the
 * introduction:
 * 
 * @code
 *         cell_matrix_M.gauss_jordan();
 *         cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E);
 * 
 * @endcode
 * 
 * Then we also need, again, to compute the matrix $C$ that is
 * used to evaluate the weak discrete gradient. This is the
 * exact same code as used in the assembly of the system
 * matrix, so we just copy it from there:
 * 
 * @code
 *         cell_matrix_G = 0;
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *             {
 *               const double div_v_i =
 *                 fe_values_dgrt[velocities].divergence(i, q);
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 {
 *                   const double phi_j_interior =
 *                     fe_values[pressure_interior].value(j, q);
 * 
 *                   cell_matrix_G(i, j) -=
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q));
 *                 }
 *             }
 * 
 *         for (const auto &face : cell->face_iterators())
 *           {
 *             fe_face_values.reinit(cell, face);
 *             fe_face_values_dgrt.reinit(cell_dgrt, face);
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q)
 *               {
 *                 const Tensor<1, dim> &normal = fe_face_values.normal_vector(q);
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i)
 *                   {
 *                     const Tensor<1, dim> v_i =
 *                       fe_face_values_dgrt[velocities].value(i, q);
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                       {
 *                         const double phi_j_face =
 *                           fe_face_values[pressure_face].value(j, q);
 * 
 *                         cell_matrix_G(i, j) +=
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q));
 *                       }
 *                   }
 *               }
 *           }
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M);
 * 
 * @endcode
 * 
 * Finally, we need to extract the pressure unknowns that
 * correspond to the current cell:
 * 
 * @code
 *         cell->get_dof_values(solution, cell_solution);
 * 
 * @endcode
 * 
 * We are now in a position to compute the local velocity
 * unknowns (with respect to the Raviart-Thomas space we are
 * projecting the term $-\mathbf K \nabla_{w,d} p_h$ into):
 * 
 * @code
 *         cell_velocity = 0;
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_velocity(k) +=
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j));
 * 
 * @endcode
 * 
 * We compute Darcy velocity.
 * This is same as cell_velocity but used to graph Darcy velocity.
 * 
 * @code
 *         cell_dgrt->get_dof_indices(local_dof_indices_dgrt);
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k)
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               darcy_velocity(local_dof_indices_dgrt[k]) +=
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j));
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_pressure_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_pressure_error</h4>
 * 

 * 
 * This part is to calculate the $L_2$ error of the pressure.  We
 * define a vector that holds the norm of the error on each cell.
 * Next, we use VectorTool::integrate_difference() to compute the
 * error in the $L_2$ norm on each cell. However, we really only
 * care about the error in the interior component of the solution
 * vector (we can't even evaluate the interface pressures at the
 * quadrature points because these are all located in the interior
 * of cells) and consequently have to use a weight function that
 * ensures that the interface component of the solution variable is
 * ignored. This is done by using the ComponentSelectFunction whose
 * arguments indicate which component we want to select (component
 * zero, i.e., the interior pressures) and how many components there
 * are in total (two).
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_pressure_error()
 *   {
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     const ComponentSelectFunction<dim> select_interior_pressure(0, 2);
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       ExactPressure<dim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe.degree + 2),
 *                                       VectorTools::L2_norm,
 *                                       &select_interior_pressure);
 * 
 *     const double L2_error = difference_per_cell.l2_norm();
 *     std::cout << "L2_error_pressure " << L2_error << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_velocity_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_velocity_error</h4>
 * 

 * 
 * In this function, we evaluate $L_2$ errors for the velocity on
 * each cell, and $L_2$ errors for the flux on faces. The function
 * relies on the `compute_postprocessed_velocity()` function having
 * previous computed, which computes the velocity field based on the
 * pressure solution that has previously been computed.
 *   

 * 
 * We are going to evaluate velocities on each cell and calculate
 * the difference between numerical and exact velocities.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::compute_velocity_errors()
 *   {
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1);
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1);
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt,
 *                                  quadrature_formula,
 *                                  update_values | update_gradients |
 *                                    update_quadrature_points |
 *                                    update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt,
 *                                           face_quadrature_formula,
 *                                           update_values |
 *                                             update_normal_vectors |
 *                                             update_quadrature_points |
 *                                             update_JxW_values);
 * 
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size();
 *     const unsigned int n_face_q_points_dgrt =
 *       fe_face_values_dgrt.get_quadrature().size();
 * 
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points_dgrt);
 *     std::vector<Tensor<1, dim>> velocity_face_values(n_face_q_points_dgrt);
 * 
 *     const FEValuesExtractors::Vector velocities(0);
 * 
 *     const ExactVelocity<dim> exact_velocity;
 * 
 *     double L2_err_velocity_cell_sqr_global = 0;
 *     double L2_err_flux_sqr                 = 0;
 * 
 * @endcode
 * 
 * Having previously computed the postprocessed velocity, we here
 * only have to extract the corresponding values on each cell and
 * face and compare it to the exact values.
 * 
 * @code
 *     for (const auto &cell_dgrt : dof_handler_dgrt.active_cell_iterators())
 *       {
 *         fe_values_dgrt.reinit(cell_dgrt);
 * 
 * @endcode
 * 
 * First compute the $L_2$ error between the postprocessed velocity
 * field and the exact one:
 * 
 * @code
 *         fe_values_dgrt[velocities].get_function_values(darcy_velocity,
 *                                                        velocity_values);
 *         double L2_err_velocity_cell_sqr_local = 0;
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q)
 *           {
 *             const Tensor<1, dim> velocity = velocity_values[q];
 *             const Tensor<1, dim> true_velocity =
 *               exact_velocity.value(fe_values_dgrt.quadrature_point(q));
 * 
 *             L2_err_velocity_cell_sqr_local +=
 *               ((velocity - true_velocity) * (velocity - true_velocity) *
 *                fe_values_dgrt.JxW(q));
 *           }
 *         L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local;
 * 
 * @endcode
 * 
 * For reconstructing the flux we need the size of cells and
 * faces. Since fluxes are calculated on faces, we have the
 * loop over all four faces of each cell. To calculate the
 * face velocity, we extract values at the quadrature points from the
 * `darcy_velocity` which we have computed previously. Then, we
 * calculate the squared velocity error in normal direction. Finally, we
 * calculate the $L_2$ flux error on the cell by appropriately scaling
 * with face and cell areas and add it to the global error.
 * 
 * @code
 *         const double cell_area = cell_dgrt->measure();
 *         for (const auto &face_dgrt : cell_dgrt->face_iterators())
 *           {
 *             const double face_length = face_dgrt->measure();
 *             fe_face_values_dgrt.reinit(cell_dgrt, face_dgrt);
 *             fe_face_values_dgrt[velocities].get_function_values(
 *               darcy_velocity, velocity_face_values);
 * 
 *             double L2_err_flux_face_sqr_local = 0;
 *             for (unsigned int q = 0; q < n_face_q_points_dgrt; ++q)
 *               {
 *                 const Tensor<1, dim> velocity = velocity_face_values[q];
 *                 const Tensor<1, dim> true_velocity =
 *                   exact_velocity.value(fe_face_values_dgrt.quadrature_point(q));
 * 
 *                 const Tensor<1, dim> &normal =
 *                   fe_face_values_dgrt.normal_vector(q);
 * 
 *                 L2_err_flux_face_sqr_local +=
 *                   ((velocity * normal - true_velocity * normal) *
 *                    (velocity * normal - true_velocity * normal) *
 *                    fe_face_values_dgrt.JxW(q));
 *               }
 *             const double err_flux_each_face =
 *               L2_err_flux_face_sqr_local / face_length * cell_area;
 *             L2_err_flux_sqr += err_flux_each_face;
 *           }
 *       }
 * 
 * @endcode
 * 
 * After adding up errors over all cells and faces, we take the
 * square root and get the $L_2$ errors of velocity and
 * flux. These we output to screen.
 * 
 * @code
 *     const double L2_err_velocity_cell =
 *       std::sqrt(L2_err_velocity_cell_sqr_global);
 *     const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr);
 * 
 *     std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl
 *               << "L2_error_flux: " << L2_err_flux_face << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationoutput_results"></a> 
 * <h4>WGDarcyEquation::output_results</h4>
 * 

 * 
 * We have two sets of results to output: the interior solution and
 * the skeleton solution. We use <code>DataOut</code> to visualize
 * interior results. The graphical output for the skeleton results
 * is done by using the DataOutFaces class.
 *   

 * 
 * In both of the output files, both the interior and the face
 * variables are stored. For the interface output, the output file
 * simply contains the interpolation of the interior pressures onto
 * the faces, but because it is undefined which of the two interior
 * pressure variables you get from the two adjacent cells, it is
 * best to ignore the interior pressure in the interface output
 * file. Conversely, for the cell interior output file, it is of
 * course impossible to show any interface pressures $p^\partial$,
 * because these are only available on interfaces and not cell
 * interiors. Consequently, you will see them shown as an invalid
 * value (such as an infinity).
 *   

 * 
 * For the cell interior output, we also want to output the velocity
 * variables. This is a bit tricky since it lives on the same mesh
 * but uses a different DoFHandler object (the pressure variables live
 * on the `dof_handler` object, the Darcy velocity on the `dof_handler_dgrt`
 * object). Fortunately, there are variations of the
 * DataOut::add_data_vector() function that allow specifying which
 * DoFHandler a vector corresponds to, and consequently we can visualize
 * the data from both DoFHandler objects within the same file.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::output_results() const
 *   {
 *     {
 *       DataOut<dim> data_out;
 * 
 * @endcode
 * 
 * First attach the pressure solution to the DataOut object:
 * 
 * @code
 *       const std::vector<std::string> solution_names = {"interior_pressure",
 *                                                        "interface_pressure"};
 *       data_out.add_data_vector(dof_handler, solution, solution_names);
 * 
 * @endcode
 * 
 * Then do the same with the Darcy velocity field, and continue
 * with writing everything out into a file.
 * 
 * @code
 *       const std::vector<std::string> velocity_names(dim, "velocity");
 *       const std::vector<
 *         DataComponentInterpretation::DataComponentInterpretation>
 *         velocity_component_interpretation(
 *           dim, DataComponentInterpretation::component_is_part_of_vector);
 *       data_out.add_data_vector(dof_handler_dgrt,
 *                                darcy_velocity,
 *                                velocity_names,
 *                                velocity_component_interpretation);
 * 
 *       data_out.build_patches(fe.degree);
 *       std::ofstream output("solution_interior.vtu");
 *       data_out.write_vtu(output);
 *     }
 * 
 *     {
 *       DataOutFaces<dim> data_out_faces(false);
 *       data_out_faces.attach_dof_handler(dof_handler);
 *       data_out_faces.add_data_vector(solution, "Pressure_Face");
 *       data_out_faces.build_patches(fe.degree);
 *       std::ofstream face_output("solution_interface.vtu");
 *       data_out_faces.write_vtu(face_output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationrun"></a> 
 * <h4>WGDarcyEquation::run</h4>
 * 

 * 
 * This is the final function of the main class. It calls the other functions
 * of our class.
 * 
 * @code
 *   template <int dim>
 *   void WGDarcyEquation<dim>::run()
 *   {
 *     std::cout << "Solving problem in " << dim << " space dimensions."
 *               << std::endl;
 *     make_grid();
 *     setup_system();
 *     assemble_system();
 *     solve();
 *     compute_postprocessed_velocity();
 *     compute_pressure_error();
 *     compute_velocity_errors();
 *     output_results();
 *   }
 * 
 * } // namespace Step61
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * This is the main function. We can change the dimension here to run in 3d.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step61::WGDarcyEquation<2> wg_darcy(0);
 *       wg_darcy.run();
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
examples/step-61/doc/results.dox



<a name="Results"></a><h1>Results</h1>


我们在运行程序时，右手边会产生解  $p = \sin(\pi x) \sin(\pi y)$  ，并且在域  $\Omega = (0,1)^2$  中具有同质的迪里希特边界条件。此外，我们选择微分算子 $\mathbf{K}$ 中的系数矩阵作为身份矩阵。我们使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 、 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 和 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 元素组合测试这一设置，可以通过使用`main()`中`WGDarcyEquation`对象的适当构造参数来选择。然后我们将可视化单元内部和面上的压力值。随着网格的细化，压力、速度和流量的收敛率对于 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 应该是1，对于 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 是2，对于 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 是3。




<a name="TestresultsoniWGQsub0subQsub0subRTsub0subiiWGQsub0subQsub0subRTsub0subi"></a><h3>Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i><i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></h3> 。


下面的数字显示了使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 元素的内部压力和表面压力。网格分别细化了2倍（顶部）和4倍（底部）。(这个数字可以在`make_grid()`函数中调整)。当网格较粗时，可以看到面压 $p^\partial$ 整齐地位于两个相邻单元的内压 $p^\circ$ 的数值之间。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_2.png" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_4.png" alt=""></td>
  </tr>
</table> 

从图中我们可以看出，随着网格的细化，最大和最小的压力值正在接近我们的预期值。由于网格是一个矩形网格，每个方向的单元数是偶数，所以我们有对称的解决方案。从右边的三维图中，我们可以看到在 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 上，压力在单元的内部是一个常数，正如预期的那样。

<a name="Convergencetableforik0iik0i"></a><h4>Convergence table for <i>k=0</i><i>k=0</i></h4> 。


我们用不同的细化网格（在 "make_grid() "函数中选择）运行代码，得到压力、速度和通量（如引言中定义的）的以下收敛率。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    1.587e-01        </td><td>        5.113e-01               </td><td>   7.062e-01 </td>
  </tr>
  <tr>
   <td>   3                  </td><td>    8.000e-02        </td><td>        2.529e-01               </td><td>   3.554e-01 </td>
  </tr>
  <tr>
   <td>   4                  </td><td>    4.006e-02        </td><td>        1.260e-01               </td><td>   1.780e-01 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>    2.004e-02        </td><td>        6.297e-02               </td><td>   8.902e-02 </td>
  </tr>
  <tr>
   <th>Conv.rate             </th><th>      1.00           </th><th>          1.00                  </th><th>      1.00   </th>
  </tr>
</table> 

我们可以看到， $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 的收敛率在1左右。当然，这与我们的理论预期相符。




<a name="TestresultsoniWGQsub1subQsub1subRTsub1subiiWGQsub1subQsub1subRTsub1subi"></a><h3>Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i><i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></h3> 。


我们可以用下一个更高的多项式度数重复上面的实验。下面的数字是使用 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 实现的内部压力和表面压力。网格被细化了4次。  与之前使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 的数字相比，在每个单元上，解决方案不再是恒定的，因为我们现在使用双线性多项式来做近似。因此，在一个内部有4个压力值，在每个面上有2个压力值。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_3d_4.png" alt=""></td>
  </tr>
</table> 

与 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 组合的相应图像相比，现在的解决方案大大增加了准确性，特别是在界面上如此接近于连续，以至于我们不再能够区分相邻单元上的界面压力 $p^\partial$ 和内部压力 $p^\circ$ 。

<a name="Convergencetableforik1iik1i"></a><h4>Convergence table for <i>k=1</i><i>k=1</i></h4> 。


以下是我们使用 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 元素组合得到的压力、速度和通量的收敛率。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
    <td>  2           </td><td>           1.613e-02      </td><td>          5.093e-02     </td><td>             7.167e-02   </td>
  </tr>
  <tr>
    <td>  3           </td><td>           4.056e-03      </td><td>          1.276e-02     </td><td>             1.802e-02    </td>
  </tr>
  <tr>
    <td>  4           </td><td>           1.015e-03      </td><td>          3.191e-03     </td><td>             4.512e-03  </td>
  </tr>
  <tr>
    <td>  5           </td><td>           2.540e-04      </td><td>          7.979e-04     </td><td>             1.128e-03  </td>
  </tr>
  <tr>
    <th>Conv.rate     </th><th>              2.00        </th><th>             2.00       </th><th>                 2.00    </th>
  </tr>
</table> 

 $WG(Q_1,Q_1;RT_{[1]})$ 的收敛率在2左右，符合预期。




<a name="TestresultsoniWGQsub2subQsub2subRTsub2subiiWGQsub2subQsub2subRTsub2subi"></a><h3>Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i><i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></h3> 。


让我们再提高一个多项式等级。以下是使用 $WG(Q_2,Q_2;RT_{[2]})$ 实现的内部压力和表面压力，网格大小为 $h = 1/32$ （即5个全局网格细化步骤）。在程序中，我们在生成图形输出时使用`data_out_face.build_patches(fe.degree)`（参见 DataOut::build_patches()), 的文档，这里意味着我们将每个2d单元内部分成4个子单元，以便提供更好的二次多项式的可视化。   <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_2d_5.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_3d_5.png" alt=""></td>
  </tr>
</table> 




<a name="Convergencetableforik2iik2i"></a><h4>Convergence table for <i>k=2</i><i>k=2</i></h4> 。


和以前一样，我们可以使用 $L_2$ 组合生成压力、速度和流量的 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 误差的收敛数据。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
     <td>  2               </td><td>       1.072e-03       </td><td>         3.375e-03       </td><td>           4.762e-03   </td>
  </tr>
  <tr>
    <td>   3               </td><td>       1.347e-04       </td><td>         4.233e-04       </td><td>           5.982e-04    </td>
  </tr>
  <tr>
    <td>   4               </td><td>       1.685e-05      </td><td>          5.295e-05       </td><td>           7.487e-05  </td>
  </tr>
  <tr>
    <td>   5               </td><td>       2.107e-06      </td><td>          6.620e-06       </td><td>           9.362e-06  </td>
  </tr>
  <tr>
    <th>Conv.rate          </th><th>         3.00         </th><th>            3.00          </th><th>              3.00    </th>
  </tr>
</table> 

再一次， $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 的收敛率符合预期，其数值在3左右。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-61.cc"
*/
