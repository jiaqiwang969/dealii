/**
@page step_60 The step-60 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#NonmatchinggridconstraintsthroughdistributedLagrangemultipliers">Non-matching grid constraints through distributed Lagrange multipliers</a>
        <li><a href="#Thetestcase">The testcase</a>
        <li><a href="#References">References</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#DistributedLagrangeProblem">DistributedLagrangeProblem</a>
        <li><a href="#DistributedLagrangeProblemParameters">DistributedLagrangeProblem::Parameters</a>
        <li><a href="#Setup">Set up</a>
        <li><a href="#Assembly">Assembly</a>
        <li><a href="#Solve">Solve</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Testcase1"> Test case 1: </a>
        <li><a href="#Testcase2and3"> Test case 2 and 3: </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Runningwithspacedimequaltothree"> Running with `spacedim` equal to three</a>
        <li><a href="#Moregeneraldomains"> More general domains </a>
        <li><a href="#Preconditioner"> Preconditioner</a>
        <li><a href="#ParallelCode"> Parallel Code </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-60/doc/intro.dox

 <br> 

<i>This program was contributed by Luca Heltai and Giovanni Alzetta, SISSA, Trieste.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.1243280,https://zenodo.org/badge/DOI/10.5281/zenodo.1243280.svg} 




<a name="Introduction"></a><h1>Introduction</h1>


<a name="NonmatchinggridconstraintsthroughdistributedLagrangemultipliers"></a><h3>Non-matching grid constraints through distributed Lagrange multipliers</h3>



在本教程中，我们考虑两个域的情况， $\Omega$ 在 $R^{\text{spacedim}}$ 和 $\Gamma$ 在 $R^{\text{dim}}$ ，其中 $\Gamma$ 嵌入在 $\Omega$ （ $\Gamma \subseteq \Omega$  ）。我们想在 $\Omega$ 上解决一个偏微分方程，对问题的解决*在嵌入域* $\Gamma$ 上强制执行一些条件。

有两种有趣的情况。

- 嵌入域 $\Gamma$ 的几何维度`dim`与域 $\Omega$ 相同（`spacedim`），也就是说， $\Gamma$ 的spacedim维度不为零，或

- 嵌入域 $\Gamma$ 的内在维度`dim`小于 $\Omega$ 的维度（`spacedim`），因此其spacedim维度为零；例如，它是一条嵌入二维域的曲线，或一个嵌入三维域的曲面。

在这两种情况下，定义限制算子 $\gamma$ 为算子，给定 $\Omega$ 上的一个连续函数，返回其在 $\Gamma$ 上的（连续）限制，即：。

\f[
\gamma : C^0(\Omega) \mapsto C^0(\Gamma), \quad \text{ s.t. } \gamma u = u|_{\Gamma} \in C^0(\Gamma),
\quad \forall u \in C^0(\Omega).
\f]

众所周知，当 $\gamma$ 的内在维度与 $\Omega$ 相同时，算子 $\gamma$ 可以扩展为 $H^1(\Omega)$ 上的连续算子，将 $H^1(\Omega)$ 的函数映射为 $H^1(\Gamma)$ 的函数。

同样的道理，在一个不太规则的范围空间（即 $H^{1/2}(\Gamma)$ ）中，当 $\Gamma$ 的维度相对于 $\Omega$ 少一个，并且 $\Gamma$ 没有边界。在这第二种情况下，算子 $\gamma$ 也被称为*轨迹*算子，对于嵌入 $\Omega$ 中的Lipschitz同维度曲线和曲面 $\Gamma$ ，它有很好的定义（阅读<a
href="https://en.wikipedia.org/wiki/Trace_operator">this wikipedia article</a>了解关于轨迹算子的进一步细节）。

同维度的情况要复杂一些，一般来说，不可能构造一个连续的跟踪算子，甚至不可能从 $H^1(\Omega)$ 到 $L^2(\Gamma)$ ，当 $\Gamma$ 的维度在二维和三维中分别为零或一的时候。

在本教程中，我们对 $\gamma$ 的进一步细节不感兴趣：我们认为扩展 $\gamma$ 是理所当然的，假设嵌入域的尺寸（`dim`）总是比嵌入域的尺寸（`spacedim`）小一或相等。

我们要解决以下微分问题：给定 $g$ 上的一个足够规则的函数 $\Gamma$ ，找到 $u$ 的解。

@f{eqnarray*}{


- \Delta u + \gamma^T \lambda &=& 0  \text{ in } \Omega\\
\gamma u &=& g  \text{ in } \Gamma \\
u & = & 0 \text{ on } \partial\Omega.


@f}



这是一个约束问题，我们正在寻找一个谐波函数 $u$ ，满足 $\partial\Omega$ 上的同质边界条件，受制于使用拉格朗日乘法器的约束 $\gamma u = g$ 。

这个问题有一个物理解释：谐波函数，即满足拉普拉斯方程的函数，可以被认为是边界值被规定的膜的位移。那么，目前的情况相当于找到一个膜的形状，对于这个膜来说，不仅边界上的位移，而且 $\Gamma$ 上的位移也是规定的。例如，如果 $\Gamma$ 是二维空间中的一条封闭曲线，那么这将是一个肥皂膜的模型，它被沿 $\partial \Omega$ 的一个线环以及沿 $\Gamma$ 的第二个线环固定住。在 $\Gamma$ 是整个区域的情况下，你可以把它看成是在障碍物上伸展的膜，其中 $\Gamma$ 是接触区域。如果接触面积不知道，我们就有一个不同的问题--称为 "障碍物问题"--在步骤41中进行建模）。

作为第一个例子，我们研究 $\partial\Omega$ 上的零迪里切特边界条件。如果我们在 $\partial\Omega$ 上应用零诺伊曼边界条件或两者的混合，同样的方程也适用。

通过引入两个无限维空间 $V(\Omega)$ 和 $Q^*(\Gamma)$ ，可以得出变分公式，分别用于解 $u$ 和拉格朗日乘子 $\lambda$ 。

将第一个方程乘以 $v \in V(\Omega)$ ，第二个方程乘以 $q \in
Q(\Gamma)$ ，在可能的情况下进行部分积分，并利用 $\partial\Omega$ 的边界条件，我们得到以下变量问题。

给出 $g$ 上的一个足够规则的函数 $\Gamma$ ，求 $u$ 的解

@f{eqnarray*}{
(\nabla u, \nabla v)_{\Omega} + (\lambda, \gamma v)_{\Gamma} &=& 0 \qquad \forall v \in V(\Omega) \\
(\gamma u, q)_{\Gamma} &=& (g,q)_{\Gamma} \qquad \forall q \in Q(\Gamma),


@f}



其中 $(\cdot, \cdot)_{\Omega}$ 和 $(\cdot, \cdot)_{\Gamma}$ 分别代表 $L^2$ 中的标量积和 $\Gamma$ 中的标量积 。

对变量公式的检查告诉我们，空间 $V(\Omega)$ 可以被认为是 $H^1_0(\Omega)$ 。空间 $Q(\Gamma)$ ，在同维度为零的情况下，应取为 $H^1(\Gamma)$ ，而在同维度为一的情况下应取为 $H^{1/2}(\Gamma)$ 。

因此函数 $g$ 应该在 $H^1(\Gamma)$ （对于同维度零的情况）或者 $H^{1/2}(\Gamma)$ （对于同维度一的情况）。这使得我们在 $Q^*(\Gamma)$ 中有一个拉格朗日乘数 $\lambda$ ，它是 $H^{-1}(\Gamma)$ 或 $H^{-1/2}(\Gamma)$ 。

对于上述问题的离散化，有两种选择。可以选择匹配的离散化，即 $\Gamma$ 的三角化与 $\Omega$ 的三角化一致，也可以选择以完全独立的方式离散化这两个域。

对于我们上面提出的简单问题，第一种选择显然更有意义：对 $\Omega$ 使用一个单一的三角形就足够了，然后根据 $\Gamma$ 施加某些约束。在步骤40中研究了这种方法的一个例子，解决方案必须保持在一个障碍物之上，这是在 $\Omega$ 上施加约束实现的。

为了解决更复杂的问题，例如域 $\Gamma$ 与时间有关的问题，第二个选项可能是一个更可行的解决方案。处理不对齐的网格本身就很复杂：为了说明如何做，我们研究一个简单的问题。

我们在此描述的技术在文献中使用了许多名称之一：<b>immersed finite element method</b>、<b>fictitious boundary method</b>、<b>distributed Lagrange multiplier method</b>等。其主要原理是，两个网格的离散化和两个有限元空间的离散化保持完整。完全独立。这种技术对于模拟流体与结构的相互作用问题特别有效，其中嵌入结构的配置是问题本身的一部分，人们要解决一个（可能是非线性）弹性问题，以确定 $\Gamma$ 的（与时间有关的）配置，以及 $\Omega
\setminus \Gamma$ 的（可能是非线性）流动问题，加上流体和固体之间界面上的耦合条件。

在这个教程程序中，我们把事情弄得简单一些，我们假设嵌入式领域的配置是以两种可能的方式之一给出的。

- 作为一个变形映射 $\psi: \Gamma_0 \mapsto \Gamma \subseteq \Omega$ ，定义在 $\Gamma_0$ 的连续有限维空间上，对于任何一个点 $x \in \Gamma_0$ ，代表其在 $\Omega$ 的坐标 $\psi(x)$ 。

- 作为 $x\in \Gamma_0$ 的位移映射 $\delta \psi(x) = \psi(x)-x$ ，代表任何一点 $x$ 的位移矢量，以使 $x$ 变形为其实际配置 $\psi(x) = x +\delta\psi(x)$ 。

我们定义嵌入式参考域 $\Gamma_0$  `embedded_grid`：在这个三角形上，我们构建一个有限维空间（`embedded_configuration_dh`），通过FE_Q对象的有限元系统（`embedded_configuration_fe`）描述变形或位移。这个有限维度空间仅用于插值用户提供的函数（`embedded_configuration_function`），代表 $\psi$ （如果参数`use_displacement`被设置为 @p false) 或 $\delta\psi$ （如果参数`use_displacement`被设置为 @p true). 

拉格朗日乘数 $\lambda$ 和用户提供的函数 $g$ 是通过另一个有限维度空间`embedded_dh`和另一个有限元素`embedded_fe`定义的，使用相同的参考域。为了考虑到域的变形，MappingFEField或MappingQEulerian对象被初始化为`embedded_configuration`向量。

在嵌入空间中，一个标准的有限维空间`space_dh`被构建在嵌入网格`space_grid`上，使用有限元素`space_fe`，几乎逐字逐句地遵循步骤6中的方法。

我们用以下方法表示空间 $V$ 和 $Q$ 的离散化

\f[
V_h(\Omega) = \text{span} \{v_i\}_{i=1}^n
\f] 和

\f[
Q_h(\Gamma) = \text{span} \{q_i\}_{i=1}^m
\f]，其中 $n$ 是`空间_dh`的尺寸， $m$ 是`嵌入_dh`的尺寸。

一旦所有的有限维空间都被定义，上述问题的变异表述给我们留下了以下有限维方程组。

\f[
\begin{pmatrix}
K & C^T \\
C & 0
\end{pmatrix}
\begin{pmatrix}
u \\
\lambda
\end{pmatrix}
=
\begin{pmatrix}
0 \\
G
\end{pmatrix}
\f]

其中

@f{eqnarray*}{
K_{ij} &\dealcoloneq& (\nabla v_j, \nabla v_i)_\Omega   \qquad i,j=1,\dots,n \\
C_{\alpha j} &\dealcoloneq& (v_j, q_\alpha)_\Gamma  \qquad j=1,\dots,n, \alpha = 1,\dots, m \\\\
G_{\alpha} &\dealcoloneq& (g, q_\alpha)_\Gamma \qquad \alpha = 1,\dots, m.


@f}



虽然矩阵 $K$ 是 $\Omega$ 上泊松问题的标准刚度矩阵，而向量 $G$ 是 $g$ 上带有强制项的有限元问题的标准右手向量。矩阵 $C$ 或其转置 $C^T$ 是非标准的，因为它们是两个不匹配的网格上的信息。

特别是，在计算 $C$ 的一个条目时出现的积分，是在 $\Gamma$ 上计算的。在有限元中，我们通常将这个积分分成来自用于离散化 $\Gamma$ 的三角形的所有单元的贡献，我们将 $K$ 上的积分转换为参考元素 $\hat K$ 上的积分，其中 $F_{K}$ 是从 $\hat K$ 到 $K$ 的映射，并且使用正交公式计算 $\hat K$ 上的积分。

\f[
C_{\alpha j} \dealcoloneq (v_j, q_\alpha)_\Gamma  = \sum_{K\in \Gamma} \int_{\hat K}
\hat q_\alpha(\hat x) (v_j \circ F_{K}) (\hat x) J_K (\hat x) \mathrm{d} \hat x =
\sum_{K\in \Gamma} \sum_{i=1}^{n_q}  \big(\hat q_\alpha(\hat x_i)  (v_j \circ F_{K}) (\hat x_i) J_K (\hat x_i) w_i \big)
\f]

计算这个和是不容易的，因为我们必须评估 $(v_j \circ F_{K})
(\hat x_i)$  。一般来说，如果 $\Gamma$ 和 $\Omega$ 没有对齐，那么 $F_{K}(\hat x_i)$ 这个点相对于 $\Omega$ 来说是完全任意的，除非我们想出一个办法，在 $\Omega$ 上的任意点上插值 $V_h(\Omega)$ 的所有基函数，否则我们无法计算出矩阵 $C$ 的一个条目需要的积分。

要评估 $(v_j \circ F_{K}) (\hat x_i)$ ，需要采取以下步骤（如下图所示）。

- 对于 $\Gamma$ 中的一个给定单元 $K$ ，计算实点 $y_i \dealcoloneq F_{K} (\hat
x_i)$ ，其中 $x_i$ 是用于 $K
\subseteq \Gamma$ 上的积分的正交点之一。

- 找到 $\Omega$ 中 $y_i$ 所在的单元。我们将称这个元素为 $T$  。

- 为了评估基函数，使用映射 $G_T$ 的逆映射，将参考元素 $\hat T$ 转换为元素 $T$  :  $v_j(y_i) = \hat
v_j \circ G^{-1}_{T} (y_i)$  。

<p align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.C_interpolation.png" alt="">  </p> 

上述三个步骤可以通过依次调用来计算。

-  GridTools::find_active_cell_around_point(),  后面是

-  Mapping::transform_real_to_unit_cell().  然后我们

- 构建一个自定义的正交公式，包含参考单元格中的点，然后

- 构建一个FEValues对象，具有给定的正交公式，并以第一步中获得的单元格为初始化。

这就是deal.II函数 VectorTools::point_value() 在任意点上评估有限元场（而不仅仅是单个形状函数）时的做法；但在这种情况下，这将是低效的。

一个更好的解决方案是使用一个方便的包装器来对一个点的集合执行前三个步骤。   GridTools::compute_point_locations(). 如果人们实际上对计算完整的耦合矩阵感兴趣，那么可以调用方法 NonMatching::create_coupling_mass_matrix(), ，该方法以有效的方式执行上述步骤，重复使用所有可能的数据结构，并将昂贵的步骤聚集在一起。这就是我们在本教程后面要使用的函数。

我们通过迭代求解器来解决最终的鞍点问题，应用于Schur补数 $S$ （其构造例如在步骤20中描述），我们使用LinearOperator类构造 $S$ 。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们在这里解决的问题与步骤4相同，不同的是我们对一个嵌入域施加了一些约束  $\Gamma$  。该教程是以独立于维度的方式编写的，在结果部分我们展示了如何改变`dim`和`spacedim`。

本教程是为`dim`等于1和`spacedim`等于2而编译的。如果你想在嵌入维度`spacedim`等于3的情况下运行程序，你很可能想改变 $\Gamma$ 的参考域，例如，你从文件中读到的东西，或者你后来变形为更有趣的封闭球体。

在默认情况下， $\Gamma$ 的共维为1，本教程程序实现的是虚构边界法。事实证明，同样的技术也用在了变量沉浸式有限元方法中，上面定义的耦合算子 $C$ 在几乎所有这些非匹配方法中都是一样的。

嵌入域被假定包括在 $\Omega$ 中，我们把它当作单位平方 $[0,1]^2$ 。虚域 $\Gamma$ 的定义可以通过参数文件修改，可以给出从参考区间 $[0,1]$ 到 $\Omega$ 中的一条曲线的映射。

如果曲线是封闭的，那么结果将类似于在边界为 $\Gamma$ 的网格上运行同一问题。在非封闭的 $\Gamma$ 的情况下，程序也能愉快地运行，尽管在这些情况下，问题的数学表述更加困难，因为 $\Gamma$ 本身就有一个边界，相对于域 $\Omega$ 来说是二维的。




<a name="References"></a><h3>References</h3>


 <ul>   <li>  Glowinski, R., T.-W.Pan, T.I. Hesla, and D.D. Joseph.1999."分布式拉格朗日乘数/虚构域方法用于颗粒物流"。   International Journal of Multiphase Flow 25 (5).Pergamon: 755-94.

 <li>  Boffi, D., L. Gastaldi, L. Heltai, and C.S. Peskin。2008."关于沉浸边界法的超弹性公式"。应用力学和工程中的计算机方法197（25-28）。

 <li>  Heltai, L., and F. Costanzo.2012."浸没式有限元方法的变量实现"。应用力学和工程中的计算机方法》229-232。   </ul> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * Most of these have been introduced elsewhere, we'll comment only on the new
 * ones.
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/timer.h>
 * 
 * @endcode
 * 
 * The parameter acceptor class is the first novelty of this tutorial program:
 * in general parameter files are used to steer the execution of a program at
 * run time. While even a simple approach saves compile time, as the same
 * executable can be run with different parameter settings, it can become
 * difficult to handle hundreds of parameters simultaneously while maintaining
 * compatibility between different programs. This is where the class
 * ParameterAcceptor proves useful.
 * 

 * 
 * This class is used to define a public interface for classes that want to use
 * a single global ParameterHandler to handle parameters. The class provides a
 * static ParameterHandler member, namely ParameterAcceptor::prm, and
 * implements the "Command design pattern" (see, for example, E. Gamma, R. Helm,
 * R. Johnson, J. Vlissides, Design Patterns: Elements of Reusable
 * Object-Oriented Software, Addison-Wesley Professional, 1994.
 * https://goo.gl/FNYByc).
 * 

 * 
 * ParameterAcceptor provides a global subscription mechanism. Whenever an
 * object of a class derived from ParameterAcceptor is constructed, a pointer
 * to that object-of-derived-type is registered, together with a section entry
 * in the parameter file. Such registry is traversed upon invocation of the
 * single function ParameterAcceptor::initialize("file.prm") which in turn makes
 * sure that all classes stored in the global registry declare the parameters
 * they will be using, and after having declared them, it reads the content of
 * `file.prm` to parse the actual parameters.
 * 

 * 
 * If you call the method ParameterHandler::add_parameter for each of the
 * parameters you want to use in your code, there is nothing else you need to
 * do. If you are using an already existing class that provides the two
 * functions `declare_parameters` and `parse_parameters`, you can still use
 * ParameterAcceptor, by encapsulating the existing class into a
 * ParameterAcceptorProxy class.
 * 

 * 
 * In this example, we'll use both strategies, using ParameterAcceptorProxy for
 * deal.II classes, and deriving our own parameter classes directly from
 * ParameterAcceptor.
 * 
 * @code
 * #include <deal.II/base/parameter_acceptor.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * @endcode
 * 
 * The other new include file is the one that contains the GridTools::Cache
 * class. The structure of deal.II, as many modern numerical libraries, is
 * organized following a Directed Acyclic Graph (DAG). A DAG is a directed graph
 * with topological ordering: each node structurally represents an object, and
 * is connected to non-root nodes by one (or more) oriented edges, from the
 * parent to the child. The most significant example of this structure is the
 * Triangulation and its Triangulation::cell_iterator structure. From a
 * Triangulation (the main node), we can access each cell (children nodes of the
 * triangulation). From the cells themselves we can access over all vertices of
 * the cell. In this simple example, the DAG structure can be represented as
 * three node types (the triangulation, the cell iterator, and the vertex)
 * connected by oriented edges from the triangulation to the cell iterators, and
 * from the cell iterator to the vertices. This has several advantages, but it
 * intrinsically creates “asymmetries”, making certain operations fast and their
 * inverse very slow: finding the vertices of a cell has low computational cost,
 * and can be done by simply traversing the DAG, while finding all the cells
 * that share a vertex requires a non-trivial computation unless a new DAG data
 * structure is added that represents the inverse search.
 * 

 * 
 * Since inverse operations are usually not needed in a finite element code,
 * these are implemented in GridTools without the use of extra data structures
 * related to the Triangulation which would make them much faster. One such data
 * structure, for example, is a map from the vertices of a Triangulation to all
 * cells that share those vertices, which would reduce the computations needed
 * to answer to the previous question.
 * 

 * 
 * Some methods, for example GridTools::find_active_cell_around_point, make
 * heavy usage of these non-standard operations. If you need to call these
 * methods more than once, it becomes convenient to store those data structures
 * somewhere. GridTools::Cache does exactly this, giving you access to
 * previously computed objects, or computing them on the fly (and then storing
 * them inside the class for later use), and making sure that whenever the
 * Triangulation is updated, also the relevant data structures are recomputed.
 * 
 * @code
 * #include <deal.II/grid/grid_tools_cache.h>
 * 
 * #include <deal.II/fe/fe.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * 
 * @endcode
 * 
 * In this example, we will be using a reference domain to describe an embedded
 * Triangulation, deformed through a finite element vector field.
 * 

 * 
 * The next two include files contain the definition of two classes that can be
 * used in these cases. MappingQEulerian allows one to describe a domain through
 * a *displacement* field, based on a FESystem[FE_Q(p)^spacedim] finite element
 * space. The second is a little more generic, and allows you to use arbitrary
 * vector FiniteElement spaces, as long as they provide a *continuous*
 * description of your domain. In this case, the description is done through the
 * actual *deformation* field, rather than a *displacement* field.
 * 

 * 
 * Which one is used depends on how the user wants to specify the reference
 * domain, and/or the actual configuration. We'll provide both options, and
 * experiment a little in the results section of this tutorial program.
 * 
 * @code
 * #include <deal.II/fe/mapping_q_eulerian.h>
 * #include <deal.II/fe/mapping_fe_field.h>
 * 
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * The parsed function class is another new entry. It allows one to create a
 * Function object, starting from a string in a parameter file which is parsed
 * into an object that you can use anywhere deal.II accepts a Function (for
 * example, for interpolation, boundary conditions, etc.).
 * 
 * @code
 * #include <deal.II/base/parsed_function.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * @endcode
 * 
 * This is the last new entry for this tutorial program. The namespace
 * NonMatching contains a few methods that are useful when performing
 * computations on non-matching grids, or on curves that are not aligned with
 * the underlying mesh.
 * 

 * 
 * We'll discuss its use in detail later on in the `setup_coupling` method.
 * 
 * @code
 * #include <deal.II/non_matching/coupling.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_direct.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/linear_operator.h>
 * #include <deal.II/lac/linear_operator_tools.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * namespace Step60
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="DistributedLagrangeProblem"></a> 
 * <h3>DistributedLagrangeProblem</h3>
 *   

 * 
 * In the DistributedLagrangeProblem, we need two parameters describing the
 * dimensions of the domain $\Gamma$ (`dim`) and of the domain $\Omega$
 * (`spacedim`).
 *   

 * 
 * These will be used to initialize a Triangulation<dim,spacedim> (for
 * $\Gamma$) and a Triangulation<spacedim,spacedim> (for $\Omega$).
 *   

 * 
 * A novelty with respect to other tutorial programs is the heavy use of
 * std::unique_ptr. These behave like classical pointers, with the advantage
 * of doing automatic house-keeping: the contained object is automatically
 * destroyed as soon as the unique_ptr goes out of scope, even if it is inside
 * a container or there's an exception. Moreover it does not allow for
 * duplicate pointers, which prevents ownership problems. We do this, because
 * we want to be able to i) construct the problem, ii) read the parameters,
 * and iii) initialize all objects according to what is specified in a
 * parameter file.
 *   

 * 
 * We construct the parameters of our problem in the internal class
 * `Parameters`, derived from ParameterAcceptor. The
 * `DistributedLagrangeProblem` class takes a const reference to a
 * `Parameters` object, so that it is not possible
 * to modify the parameters from within the DistributedLagrangeProblem class
 * itself.
 *   

 * 
 * We could have initialized the parameters first, and then pass the
 * parameters to the DistributedLagrangeProblem assuming all entries are set
 * to the desired values, but this has two disadvantages:
 *   

 * 
 * - We should not make assumptions on how the user initializes a class that
 * is not under our direct control. If the user fails to initialize the
 * class, we should notice and throw an exception;
 *   

 * 
 * - Not all objects that need to read parameters from a parameter file may
 * be available when we construct the Parameters;
 * this is often the case for complex programs, with multiple physics, or
 * where we reuse existing code in some external classes. We simulate this by
 * keeping some "complex" objects, like ParsedFunction objects, inside the
 * `DistributedLagrangeProblem` instead of inside the
 * `Parameters`.
 *   

 * 
 * Here we assume that upon construction, the classes that build up our
 * problem are not usable yet. Parsing the parameter file is what ensures we
 * have all ingredients to build up our classes, and we design them so that if
 * parsing fails, or is not executed, the run is aborted.
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim = dim>
 *   class DistributedLagrangeProblem
 *   {
 *   public:
 * @endcode
 * 
 * The `Parameters` class is derived from ParameterAcceptor. This allows us
 * to use the ParameterAcceptor::add_parameter() method in its constructor.
 *     

 * 
 * The members of this function are all non-const, but the
 * `DistributedLagrangeProblem` class takes a const reference to a
 * `Parameters` object: this ensures that
 * parameters are not modified from within the `DistributedLagrangeProblem`
 * class.
 * 
 * @code
 *     class Parameters : public ParameterAcceptor
 *     {
 *     public:
 *       Parameters();
 * 
 * @endcode
 * 
 * The parameters now described can all be set externally using a
 * parameter file: if no parameter file is present when running the
 * executable, the program will create a "parameters.prm" file with the
 * default values defined here, and then abort to give the user a chance
 * to modify the parameters.prm file.
 * 

 * 
 * Initial refinement for the embedding grid, corresponding to the domain
 * $\Omega$.
 * 
 * @code
 *       unsigned int initial_refinement = 4;
 * 
 * @endcode
 * 
 * The interaction between the embedded grid $\Omega$ and the embedding
 * grid $\Gamma$ is handled through the computation of $C$, which
 * involves all cells of $\Omega$ overlapping with parts of $\Gamma$:
 * a higher refinement of such cells might improve quality of our
 * computations.
 * For this reason we define `delta_refinement`: if it is greater
 * than zero, then we mark each cell of the space grid that contains
 * a vertex of the embedded grid and its neighbors, execute the
 * refinement, and repeat this process `delta_refinement` times.
 * 
 * @code
 *       unsigned int delta_refinement = 3;
 * 
 * @endcode
 * 
 * Starting refinement of the embedded grid, corresponding to the domain
 * $\Gamma$.
 * 
 * @code
 *       unsigned int initial_embedded_refinement = 8;
 * 
 * @endcode
 * 
 * The list of boundary ids where we impose homogeneous Dirichlet boundary
 * conditions. On the remaining boundary ids (if any), we impose
 * homogeneous Neumann boundary conditions.
 * As a default problem we have zero Dirichlet boundary conditions on
 * $\partial \Omega$
 * 
 * @code
 *       std::list<types::boundary_id> homogeneous_dirichlet_ids{0, 1, 2, 3};
 * 
 * @endcode
 * 
 * FiniteElement degree of the embedding space: $V_h(\Omega)$
 * 
 * @code
 *       unsigned int embedding_space_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * FiniteElement degree of the embedded space: $Q_h(\Gamma)$
 * 
 * @code
 *       unsigned int embedded_space_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * FiniteElement degree of the space used to describe the deformation
 * of the embedded domain
 * 
 * @code
 *       unsigned int embedded_configuration_finite_element_degree = 1;
 * 
 * @endcode
 * 
 * Order of the quadrature formula used to integrate the coupling
 * 
 * @code
 *       unsigned int coupling_quadrature_order = 3;
 * 
 * @endcode
 * 
 * If set to true, then the embedded configuration function is
 * interpreted as a displacement function
 * 
 * @code
 *       bool use_displacement = false;
 * 
 * @endcode
 * 
 * Level of verbosity to use in the output
 * 
 * @code
 *       unsigned int verbosity_level = 10;
 * 
 * @endcode
 * 
 * A flag to keep track if we were initialized or not
 * 
 * @code
 *       bool initialized = false;
 *     };
 * 
 *     DistributedLagrangeProblem(const Parameters &parameters);
 * 
 * @endcode
 * 
 * Entry point for the DistributedLagrangeProblem
 * 
 * @code
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * Object containing the actual parameters
 * 
 * @code
 *     const Parameters &parameters;
 * 
 * @endcode
 * 
 * The following functions are similar to all other tutorial programs, with
 * the exception that we now need to set up things for two different
 * families of objects, namely the ones related to the *embedding* grids,
 * and the ones related to the *embedded* one.
 * 

 * 
 * 
 * @code
 *     void setup_grids_and_dofs();
 * 
 *     void setup_embedding_dofs();
 * 
 *     void setup_embedded_dofs();
 * 
 * @endcode
 * 
 * The only unconventional function we have here is the `setup_coupling()`
 * method, used to generate the sparsity patter for the coupling matrix $C$.
 * 

 * 
 * 
 * @code
 *     void setup_coupling();
 * 
 *     void assemble_system();
 * 
 *     void solve();
 * 
 *     void output_results();
 * 
 * 
 * @endcode
 * 
 * first we gather all the objects related to the embedding space geometry
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<Triangulation<spacedim>> space_grid;
 *     std::unique_ptr<GridTools::Cache<spacedim, spacedim>>
 *                                              space_grid_tools_cache;
 *     std::unique_ptr<FiniteElement<spacedim>> space_fe;
 *     std::unique_ptr<DoFHandler<spacedim>>    space_dh;
 * 
 * @endcode
 * 
 * Then the ones related to the embedded grid, with the DoFHandler
 * associated to the Lagrange multiplier `lambda`
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<Triangulation<dim, spacedim>> embedded_grid;
 *     std::unique_ptr<FiniteElement<dim, spacedim>> embedded_fe;
 *     std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_dh;
 * 
 * @endcode
 * 
 * And finally, everything that is needed to *deform* the embedded
 * triangulation
 * 
 * @code
 *     std::unique_ptr<FiniteElement<dim, spacedim>> embedded_configuration_fe;
 *     std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_configuration_dh;
 *     Vector<double>                                embedded_configuration;
 * 
 * @endcode
 * 
 * The ParameterAcceptorProxy class is a "transparent" wrapper derived
 * from both ParameterAcceptor and the type passed as its template
 * parameter. At construction, the arguments are split into two parts: the
 * first argument is an std::string, forwarded to the ParameterAcceptor
 * class, and containing the name of the section that should be used for
 * this class, while all the remaining arguments are forwarded to the
 * constructor of the templated type, in this case, to the
 * Functions::ParsedFunction constructor.
 *     

 * 
 * This class allows you to use existing classes in conjunction with the
 * ParameterAcceptor registration mechanism, provided that those classes
 * have the members `declare_parameters()` and `parse_parameters()`.
 *     

 * 
 * This is the case here, making it fairly easy to exploit the
 * Functions::ParsedFunction class: instead of requiring users to create new
 * Function objects in their code for the RHS, boundary functions, etc.,
 * (like it is done in most of the other tutorials), here we allow the user
 * to use deal.II interface to muParser (http://muparser.beltoforion.de),
 * where the specification of the function is not done at compile time, but
 * at run time, using a string that is parsed into an actual Function
 * object.
 *     

 * 
 * In this case, the `embedded_configuration_function` is a vector valued
 * Function that can be interpreted as either a *deformation* or a
 * *displacement* according to the boolean value of
 * `parameters.use_displacement`. The number of components is specified
 * later on in the construction.
 * 

 * 
 * 
 * @code
 *     ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
 *       embedded_configuration_function;
 * 
 *     std::unique_ptr<Mapping<dim, spacedim>> embedded_mapping;
 * 
 * @endcode
 * 
 * We do the same thing to specify the value of the function $g$,
 * which is what we want our solution to be in the embedded space.
 * In this case the Function is a scalar one.
 * 
 * @code
 *     ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
 *       embedded_value_function;
 * 
 * @endcode
 * 
 * Similarly to what we have done with the Functions::ParsedFunction class,
 * we repeat the same for the ReductionControl class, allowing us to
 * specify all possible stopping criteria for the Schur complement
 * iterative solver we'll use later on.
 * 
 * @code
 *     ParameterAcceptorProxy<ReductionControl> schur_solver_control;
 * 
 * @endcode
 * 
 * Next we gather all SparsityPattern, SparseMatrix, and Vector objects
 * we'll need
 * 
 * @code
 *     SparsityPattern stiffness_sparsity;
 *     SparsityPattern coupling_sparsity;
 * 
 *     SparseMatrix<double> stiffness_matrix;
 *     SparseMatrix<double> coupling_matrix;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     Vector<double> solution;
 *     Vector<double> rhs;
 * 
 *     Vector<double> lambda;
 *     Vector<double> embedded_rhs;
 *     Vector<double> embedded_value;
 * 
 * @endcode
 * 
 * The TimerOutput class is used to provide some statistics on
 * the performance of our program.
 * 
 * @code
 *     TimerOutput monitor;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="DistributedLagrangeProblemParameters"></a> 
 * <h3>DistributedLagrangeProblem::Parameters</h3>
 *   

 * 
 * At construction time, we initialize also the ParameterAcceptor class, with
 * the section name we want our problem to use when parsing the parameter
 * file.
 *   

 * 
 * Parameter files can be organized into section/subsection/etc.:
 * this has the advantage that defined objects share parameters when
 * sharing the same section/subsection/etc. ParameterAcceptor allows
 * to specify the section name using Unix conventions on paths.
 * If the section name starts with a slash ("/"), then the section is
 * interpreted as an *absolute path*, ParameterAcceptor enters a subsection
 * for each directory in the path, using the last name it encountered as
 * the landing subsection for the current class.
 *   

 * 
 * For example, if you construct your class using
 * `ParameterAcceptor("/first/second/third/My Class")`, the parameters will be
 * organized as follows:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *        ... # all the parameters
 *       end
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * Internally, the *current path* stored in ParameterAcceptor is now
 * considered to be "/first/second/third/", i.e. when you specify an
 * absolute path, ParameterAcceptor *changes* the current section to the
 * current path, i.e. to the path of the section name until the *last* "/".
 *   

 * 
 * You can now construct another class derived from ParameterAcceptor using a
 * relative path (e.g., `ParameterAcceptor("My Other Class")`) instead of the
 * absolute one (e.g. `ParameterAcceptor("/first/second/third/My Other
 * Class")`), obtaining:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *         ... # all the parameters
 *       end
 *       subsection My Other Class
 *         ... # all the parameters of MyOtherClass
 *       end
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * If the section name *ends* with a slash then subsequent classes will
 * interpret this as a full path: for example, similar to the one above, if
 * we have two classes, one initialized with
 * `ParameterAcceptor("/first/second/third/My Class/")`
 * and the other with `ParameterAcceptor("My Other Class")`, then the
 * resulting parameter file will look like:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * # Example parameter file
 * subsection first
 *   subsection second
 *     subsection third
 *       subsection My Class
 *         ... # all the parameters of MyClass
 *         ... # notice My Class subsection does not end here
 *         subsection My Other Class
 *           ... # all the parameters of MyOtherClass
 *         end # of subsection My Other Class
 *       end # of subsection My Class
 *     end
 *   end
 * end
 * @endcode
 * </div>
 *   

 * 
 * We are going to exploit this, by making our
 * `Parameters` the *parent* of all subsequently
 * constructed classes. Since most of the other classes are members of
 * `DistributedLagrangeProblem` this allows, for example, to construct two
 * `DistributedLagrangeProblem` for two different dimensions, without having
 * conflicts in the parameters for the two problems.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   DistributedLagrangeProblem<dim, spacedim>::Parameters::Parameters()
 *     : ParameterAcceptor("/Distributed Lagrange<" +
 *                         Utilities::int_to_string(dim) + "," +
 *                         Utilities::int_to_string(spacedim) + ">/")
 *   {
 * @endcode
 * 
 * The ParameterAcceptor::add_parameter() function does a few things:
 *     

 * 
 * - enters the subsection specified at construction time to
 * ParameterAcceptor
 *     

 * 
 * - calls the ParameterAcceptor::prm.add_parameter() function
 *     

 * 
 * - calls any signal you may have attached to
 * ParameterAcceptor::declare_parameters_call_back
 *     

 * 
 * - leaves the subsection
 *     

 * 
 * In turn, ParameterAcceptor::prm.add_parameter
 *     

 * 
 * - declares an entry in the parameter handler for the given variable;
 *     

 * 
 * - takes the current value of the variable
 *     

 * 
 * - transforms it to a string, used as the default value for the parameter
 * file
 *     

 * 
 * - attaches an *action* to ParameterAcceptor::prm that monitors when a
 * file is parsed, or when an entry is set, and when this happens, it
 * updates the value of the variable passed to `add_parameter()` by setting
 * it to whatever was specified in the input file (of course, after the
 * input file has been parsed and the text representation converted to the
 * type of the variable).
 * 
 * @code
 *     add_parameter("Initial embedding space refinement", initial_refinement);
 * 
 *     add_parameter("Initial embedded space refinement",
 *                   initial_embedded_refinement);
 * 
 *     add_parameter("Local refinements steps near embedded domain",
 *                   delta_refinement);
 * 
 *     add_parameter("Homogeneous Dirichlet boundary ids",
 *                   homogeneous_dirichlet_ids);
 * 
 *     add_parameter("Use displacement in embedded interface", use_displacement);
 * 
 *     add_parameter("Embedding space finite element degree",
 *                   embedding_space_finite_element_degree);
 * 
 *     add_parameter("Embedded space finite element degree",
 *                   embedded_space_finite_element_degree);
 * 
 *     add_parameter("Embedded configuration finite element degree",
 *                   embedded_configuration_finite_element_degree);
 * 
 *     add_parameter("Coupling quadrature order", coupling_quadrature_order);
 * 
 *     add_parameter("Verbosity level", verbosity_level);
 * 
 * @endcode
 * 
 * Once the parameter file has been parsed, then the parameters are good to
 * go. Set the internal variable `initialized` to true.
 * 
 * @code
 *     parse_parameters_call_back.connect([&]() -> void { initialized = true; });
 *   }
 * 
 * @endcode
 * 
 * The constructor is pretty standard, with the exception of the
 * `ParameterAcceptorProxy` objects, as explained earlier.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   DistributedLagrangeProblem<dim, spacedim>::DistributedLagrangeProblem(
 *     const Parameters &parameters)
 *     : parameters(parameters)
 *     , embedded_configuration_function("Embedded configuration", spacedim)
 *     , embedded_value_function("Embedded value")
 *     , schur_solver_control("Schur solver control")
 *     , monitor(std::cout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
 *   {
 * @endcode
 * 
 * Here is a way to set default values for a ParameterAcceptor class
 * that was constructed using ParameterAcceptorProxy.
 *     

 * 
 * In this case, we set the default deformation of the embedded grid to be a
 * circle with radius $R$ and center $(Cx, Cy)$, we set the default value
 * for the embedded_value_function to be the constant one, and specify some
 * sensible values for the SolverControl object.
 *     

 * 
 * It is fundamental for $\Gamma$ to be embedded: from the definition of
 * $C_{\alpha j}$ is clear that, if $\Gamma \not\subseteq \Omega$, certain
 * rows of the matrix $C$ will be zero. This would be a problem, as the
 * Schur complement method requires $C$ to have full column rank.
 * 
 * @code
 *     embedded_configuration_function.declare_parameters_call_back.connect(
 *       []() -> void {
 *         ParameterAcceptor::prm.set("Function constants", "R=.3, Cx=.4, Cy=.4");
 * 
 * 
 *         ParameterAcceptor::prm.set("Function expression",
 *                                    "R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy");
 *       });
 * 
 *     embedded_value_function.declare_parameters_call_back.connect(
 *       []() -> void { ParameterAcceptor::prm.set("Function expression", "1"); });
 * 
 *     schur_solver_control.declare_parameters_call_back.connect([]() -> void {
 *       ParameterAcceptor::prm.set("Max steps", "1000");
 *       ParameterAcceptor::prm.set("Reduction", "1.e-12");
 *       ParameterAcceptor::prm.set("Tolerance", "1.e-12");
 *     });
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Setup"></a> 
 * <h3>Set up</h3>
 *   

 * 
 * The function `DistributedLagrangeProblem::setup_grids_and_dofs()` is used
 * to set up the finite element spaces. Notice how `std::make_unique` is
 * used to create objects wrapped inside `std::unique_ptr` objects.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_grids_and_dofs()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Setup grids and dofs");
 * 
 * @endcode
 * 
 * Initializing $\Omega$: constructing the Triangulation and wrapping it
 * into a `std::unique_ptr` object
 * 
 * @code
 *     space_grid = std::make_unique<Triangulation<spacedim>>();
 * 
 * @endcode
 * 
 * Next, we actually create the triangulation using
 * GridGenerator::hyper_cube(). The last argument is set to true: this
 * activates colorization (i.e., assigning different boundary indicators to
 * different parts of the boundary), which we use to assign the Dirichlet
 * and Neumann conditions.
 * 
 * @code
 *     GridGenerator::hyper_cube(*space_grid, 0, 1, true);
 * 
 * @endcode
 * 
 * Once we constructed a Triangulation, we refine it globally according to
 * the specifications in the parameter file, and construct a
 * GridTools::Cache with it.
 * 
 * @code
 *     space_grid->refine_global(parameters.initial_refinement);
 *     space_grid_tools_cache =
 *       std::make_unique<GridTools::Cache<spacedim, spacedim>>(*space_grid);
 * 
 * @endcode
 * 
 * The same is done with the embedded grid. Since the embedded grid is
 * deformed, we first need to setup the deformation mapping. We do so in the
 * following few lines:
 * 
 * @code
 *     embedded_grid = std::make_unique<Triangulation<dim, spacedim>>();
 *     GridGenerator::hyper_cube(*embedded_grid);
 *     embedded_grid->refine_global(parameters.initial_embedded_refinement);
 * 
 *     embedded_configuration_fe = std::make_unique<FESystem<dim, spacedim>>(
 *       FE_Q<dim, spacedim>(
 *         parameters.embedded_configuration_finite_element_degree),
 *       spacedim);
 * 
 *     embedded_configuration_dh =
 *       std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);
 * 
 *     embedded_configuration_dh->distribute_dofs(*embedded_configuration_fe);
 *     embedded_configuration.reinit(embedded_configuration_dh->n_dofs());
 * 
 * @endcode
 * 
 * Once we have defined a finite dimensional space for the deformation, we
 * interpolate the `embedded_configuration_function` defined in the
 * parameter file:
 * 
 * @code
 *     VectorTools::interpolate(*embedded_configuration_dh,
 *                              embedded_configuration_function,
 *                              embedded_configuration);
 * 
 * @endcode
 * 
 * Now we can interpret it according to what the user has specified in the
 * parameter file: as a displacement, in which case we construct a mapping
 * that *displaces* the position of each support point of our configuration
 * finite element space by the specified amount on the corresponding
 * configuration vector, or as an absolution position.
 *     

 * 
 * In the first case, the class MappingQEulerian offers its services, while
 * in the second one, we'll use the class MappingFEField. They are in fact
 * very similar. MappingQEulerian will only work for systems of FE_Q finite
 * element spaces, where the displacement vector is stored in the first
 * `spacedim` components of the FESystem, and the degree given as a
 * parameter at construction time, must match the degree of the first
 * `spacedim` components.
 *     

 * 
 * The class MappingFEField is slightly more general, in that it allows you
 * to select arbitrary FiniteElement types when constructing your
 * approximation. Naturally some choices may (or may not) make sense,
 * according to the type of FiniteElement you choose. MappingFEField
 * implements the pure iso-parametric concept, and can be used, for example,
 * to implement iso-geometric analysis codes in deal.II, by combining it
 * with the FE_Bernstein finite element class. In this example, we'll use
 * the two interchangeably, by taking into account the fact that one
 * configuration will be a `displacement`, while the other will be an
 * absolute `deformation` field.
 * 

 * 
 * 
 * @code
 *     if (parameters.use_displacement == true)
 *       embedded_mapping =
 *         std::make_unique<MappingQEulerian<dim, Vector<double>, spacedim>>(
 *           parameters.embedded_configuration_finite_element_degree,
 *           *embedded_configuration_dh,
 *           embedded_configuration);
 *     else
 *       embedded_mapping =
 *         std::make_unique<MappingFEField<dim, spacedim, Vector<double>>>(
 *           *embedded_configuration_dh, embedded_configuration);
 * 
 *     setup_embedded_dofs();
 * 
 * @endcode
 * 
 * In this tutorial program we not only refine $\Omega$ globally,
 * but also allow a local refinement depending on the position of $\Gamma$,
 * according to the value of `parameters.delta_refinement`, that we use to
 * decide how many rounds of local refinement we should do on $\Omega$,
 * corresponding to the position of $\Gamma$.
 *     

 * 
 * With the mapping in place, it is now possible to query what is the
 * location of all support points associated with the `embedded_dh`, by
 * calling the method DoFTools::map_dofs_to_support_points.
 *     

 * 
 * This method has two variants. One that does *not* take a Mapping, and
 * one that takes a Mapping. If you use the second type, like we are doing
 * in this case, the support points are computed through the specified
 * mapping, which can manipulate them accordingly.
 *     

 * 
 * This is precisely what the `embedded_mapping` is there for.
 * 
 * @code
 *     std::vector<Point<spacedim>> support_points(embedded_dh->n_dofs());
 *     if (parameters.delta_refinement != 0)
 *       DoFTools::map_dofs_to_support_points(*embedded_mapping,
 *                                            *embedded_dh,
 *                                            support_points);
 * 
 * @endcode
 * 
 * Once we have the support points of the embedded finite element space, we
 * would like to identify what cells of the embedding space contain what
 * support point, to get a chance at refining the embedding grid where it is
 * necessary, i.e., where the embedded grid is. This can be done manually,
 * by looping over each support point, and then calling the method
 * Mapping::transform_real_to_unit_cell for each cell of the embedding
 * space, until we find one that returns points in the unit reference cell,
 * or it can be done in a more intelligent way.
 *     

 * 
 * The GridTools::find_active_cell_around_point is a possible option that
 * performs the above task in a cheaper way, by first identifying the
 * closest vertex of the embedding Triangulation to the target point, and
 * then by calling Mapping::transform_real_to_unit_cell only for those cells
 * that share the found vertex.
 *     

 * 
 * In fact, there are algorithms in the GridTools namespace that exploit a
 * GridTools::Cache object, and possibly a KDTree object to speed up these
 * operations as much as possible.
 *     

 * 
 * The simplest way to exploit the maximum speed is by calling a
 * specialized method, GridTools::compute_point_locations, that will store a
 * lot of useful information and data structures during the first point
 * search, and then reuse all of this for subsequent points.
 *     

 * 
 * GridTools::compute_point_locations returns a tuple where the first
 * element is a vector of cells containing the input points, in this
 * case support_points. For refinement, this is the only information we
 * need, and this is exactly what happens now.
 *     

 * 
 * When we need to assemble a coupling matrix, however, we'll also need the
 * reference location of each point to evaluate the basis functions of the
 * embedding space. The other elements of the tuple returned by
 * GridTools::compute_point_locations allow you to reconstruct, for each
 * point, what cell contains it, and what is the location in the reference
 * cell of the given point. Since this information is better grouped into
 * cells, then this is what the algorithm returns: a tuple, containing a
 * vector of all cells that have at least one point in them, together with a
 * list of all reference points and their corresponding index in the
 * original vector.
 *     

 * 
 * In the following loop, we will be ignoring all returned objects except
 * the first, identifying all cells contain at least one support point of
 * the embedded space. This allows for a simple adaptive refinement
 * strategy: refining these cells and their neighbors.
 *     

 * 
 * Notice that we need to do some sanity checks, in the sense that we want
 * to have an embedding grid which is well refined around the embedded grid,
 * but where two consecutive support points lie either in the same cell, or
 * in neighbor embedding cells.
 *     

 * 
 * This is only possible if we ensure that the smallest cell size of the
 * embedding grid is nonetheless bigger than the largest cell size of the
 * embedded grid. Since users can modify both levels of refinements, as well
 * as the amount of local refinement they want around the embedded grid, we
 * make sure that the resulting meshes satisfy our requirements, and if this
 * is not the case, we bail out with an exception.
 * 
 * @code
 *     for (unsigned int i = 0; i < parameters.delta_refinement; ++i)
 *       {
 *         const auto point_locations =
 *           GridTools::compute_point_locations(*space_grid_tools_cache,
 *                                              support_points);
 *         const auto &cells = std::get<0>(point_locations);
 *         for (auto &cell : cells)
 *           {
 *             cell->set_refine_flag();
 *             for (const auto face_no : cell->face_indices())
 *               if (!cell->at_boundary(face_no))
 *                 cell->neighbor(face_no)->set_refine_flag();
 *           }
 *         space_grid->execute_coarsening_and_refinement();
 *       }
 * 
 * @endcode
 * 
 * In order to construct a well posed coupling interpolation operator $C$,
 * there are some constraints on the relative dimension of the grids between
 * the embedding and the embedded domains. The coupling operator $C$ and the
 * spaces $V$ and $Q$ have to satisfy an inf-sup condition in order for the
 * problem to have a solution. It turns out that the non-matching $L^2$
 * projection satisfies such inf-sup, provided that the spaces $V$ and $Q$
 * are compatible between each other (for example, provided that they are
 * chosen to be the ones described in the introduction).
 *     

 * 
 * However, the *discrete* inf-sup condition must also hold. No
 * complications arise here, but it turns out that the discrete inf-sup
 * constant deteriorates when the non-matching grids have local diameters
 * that are too far away from each other. In particular, it turns out that
 * if you choose an embedding grid which is *finer* with respect to the
 * embedded grid, the inf-sup constant deteriorates much more than if you
 * let the embedded grid be finer.
 *     

 * 
 * In order to avoid issues, in this tutorial we will throw an exception if
 * the parameters chosen by the user are such that the maximal diameter of
 * the embedded grid is greater than the minimal diameter of the embedding
 * grid.
 *     

 * 
 * This choice guarantees that almost every cell of the embedded grid spans
 * no more than two cells of the embedding grid, with some rare exceptions,
 * that are negligible in terms of the resulting inf-sup.
 * 
 * @code
 *     const double embedded_space_maximal_diameter =
 *       GridTools::maximal_cell_diameter(*embedded_grid, *embedded_mapping);
 *     double embedding_space_minimal_diameter =
 *       GridTools::minimal_cell_diameter(*space_grid);
 * 
 *     deallog << "Embedding minimal diameter: "
 *             << embedding_space_minimal_diameter
 *             << ", embedded maximal diameter: "
 *             << embedded_space_maximal_diameter << ", ratio: "
 *             << embedded_space_maximal_diameter /
 *                  embedding_space_minimal_diameter
 *             << std::endl;
 * 
 *     AssertThrow(embedded_space_maximal_diameter <
 *                   embedding_space_minimal_diameter,
 *                 ExcMessage(
 *                   "The embedding grid is too refined (or the embedded grid "
 *                   "is too coarse). Adjust the parameters so that the minimal "
 *                   "grid size of the embedding grid is larger "
 *                   "than the maximal grid size of the embedded grid."));
 * 
 * @endcode
 * 
 * $\Omega$ has been refined and we can now set up its DoFs
 * 
 * @code
 *     setup_embedding_dofs();
 *   }
 * 
 * @endcode
 * 
 * We now set up the DoFs of $\Omega$ and $\Gamma$: since they are
 * fundamentally independent (except for the fact that $\Omega$'s mesh is more
 * refined "around"
 * $\Gamma$) the procedure is standard.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_embedding_dofs()
 *   {
 *     space_dh = std::make_unique<DoFHandler<spacedim>>(*space_grid);
 *     space_fe = std::make_unique<FE_Q<spacedim>>(
 *       parameters.embedding_space_finite_element_degree);
 *     space_dh->distribute_dofs(*space_fe);
 * 
 *     DoFTools::make_hanging_node_constraints(*space_dh, constraints);
 *     for (auto id : parameters.homogeneous_dirichlet_ids)
 *       {
 *         VectorTools::interpolate_boundary_values(
 *           *space_dh, id, Functions::ZeroFunction<spacedim>(), constraints);
 *       }
 *     constraints.close();
 * 
 * @endcode
 * 
 * By definition the stiffness matrix involves only $\Omega$'s DoFs
 * 
 * @code
 *     DynamicSparsityPattern dsp(space_dh->n_dofs(), space_dh->n_dofs());
 *     DoFTools::make_sparsity_pattern(*space_dh, dsp, constraints);
 *     stiffness_sparsity.copy_from(dsp);
 *     stiffness_matrix.reinit(stiffness_sparsity);
 *     solution.reinit(space_dh->n_dofs());
 *     rhs.reinit(space_dh->n_dofs());
 * 
 *     deallog << "Embedding dofs: " << space_dh->n_dofs() << std::endl;
 *   }
 * 
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_embedded_dofs()
 *   {
 *     embedded_dh = std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid);
 *     embedded_fe = std::make_unique<FE_Q<dim, spacedim>>(
 *       parameters.embedded_space_finite_element_degree);
 *     embedded_dh->distribute_dofs(*embedded_fe);
 * 
 * @endcode
 * 
 * By definition the rhs of the system we're solving involves only a zero
 * vector and $G$, which is computed using only $\Gamma$'s DoFs
 * 
 * @code
 *     lambda.reinit(embedded_dh->n_dofs());
 *     embedded_rhs.reinit(embedded_dh->n_dofs());
 *     embedded_value.reinit(embedded_dh->n_dofs());
 * 
 *     deallog << "Embedded dofs: " << embedded_dh->n_dofs() << std::endl;
 *   }
 * 
 * @endcode
 * 
 * Creating the coupling sparsity pattern is a complex operation,
 * but it can be easily done using the
 * NonMatching::create_coupling_sparsity_pattern, which requires the
 * two DoFHandler objects, the quadrature points for the coupling,
 * a DynamicSparsityPattern (which then needs to be copied into the
 * sparsity one, as usual), the component mask for the embedding and
 * embedded Triangulation (which we leave empty) and the mappings
 * for both the embedding and the embedded Triangulation.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::setup_coupling()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Setup coupling");
 * 
 *     QGauss<dim> quad(parameters.coupling_quadrature_order);
 * 
 *     DynamicSparsityPattern dsp(space_dh->n_dofs(), embedded_dh->n_dofs());
 * 
 *     NonMatching::create_coupling_sparsity_pattern(*space_grid_tools_cache,
 *                                                   *space_dh,
 *                                                   *embedded_dh,
 *                                                   quad,
 *                                                   dsp,
 *                                                   AffineConstraints<double>(),
 *                                                   ComponentMask(),
 *                                                   ComponentMask(),
 *                                                   *embedded_mapping);
 *     coupling_sparsity.copy_from(dsp);
 *     coupling_matrix.reinit(coupling_sparsity);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Assembly"></a> 
 * <h3>Assembly</h3>
 *   

 * 
 * The following function creates the matrices: as noted before computing the
 * stiffness matrix and the rhs is a standard procedure.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::assemble_system()
 *   {
 *     {
 *       TimerOutput::Scope timer_section(monitor, "Assemble system");
 * 
 * @endcode
 * 
 * Embedding stiffness matrix $K$, and the right hand side $G$.
 * 
 * @code
 *       MatrixTools::create_laplace_matrix(
 *         *space_dh,
 *         QGauss<spacedim>(2 * space_fe->degree + 1),
 *         stiffness_matrix,
 *         static_cast<const Function<spacedim> *>(nullptr),
 *         constraints);
 * 
 *       VectorTools::create_right_hand_side(*embedded_mapping,
 *                                           *embedded_dh,
 *                                           QGauss<dim>(2 * embedded_fe->degree +
 *                                                       1),
 *                                           embedded_value_function,
 *                                           embedded_rhs);
 *     }
 *     {
 *       TimerOutput::Scope timer_section(monitor, "Assemble coupling system");
 * 
 * @endcode
 * 
 * To compute the coupling matrix we use the
 * NonMatching::create_coupling_mass_matrix tool, which works similarly to
 * NonMatching::create_coupling_sparsity_pattern.
 * 
 * @code
 *       QGauss<dim> quad(parameters.coupling_quadrature_order);
 *       NonMatching::create_coupling_mass_matrix(*space_grid_tools_cache,
 *                                                *space_dh,
 *                                                *embedded_dh,
 *                                                quad,
 *                                                coupling_matrix,
 *                                                AffineConstraints<double>(),
 *                                                ComponentMask(),
 *                                                ComponentMask(),
 *                                                *embedded_mapping);
 * 
 *       VectorTools::interpolate(*embedded_mapping,
 *                                *embedded_dh,
 *                                embedded_value_function,
 *                                embedded_value);
 *     }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solve"></a> 
 * <h3>Solve</h3>
 *   

 * 
 * All parts have been assembled: we solve the system
 * using the Schur complement method
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::solve()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Solve system");
 * 
 * @endcode
 * 
 * Start by creating the inverse stiffness matrix
 * 
 * @code
 *     SparseDirectUMFPACK K_inv_umfpack;
 *     K_inv_umfpack.initialize(stiffness_matrix);
 * 
 * @endcode
 * 
 * Initializing the operators, as described in the introduction
 * 
 * @code
 *     auto K  = linear_operator(stiffness_matrix);
 *     auto Ct = linear_operator(coupling_matrix);
 *     auto C  = transpose_operator(Ct);
 * 
 *     auto K_inv = linear_operator(K, K_inv_umfpack);
 * 
 * @endcode
 * 
 * Using the Schur complement method
 * 
 * @code
 *     auto                     S = C * K_inv * Ct;
 *     SolverCG<Vector<double>> solver_cg(schur_solver_control);
 *     auto S_inv = inverse_operator(S, solver_cg, PreconditionIdentity());
 * 
 *     lambda = S_inv * embedded_rhs;
 * 
 *     solution = K_inv * Ct * lambda;
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * @endcode
 * 
 * The following function simply generates standard result output on two
 * separate files, one for each mesh.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::output_results()
 *   {
 *     TimerOutput::Scope timer_section(monitor, "Output results");
 * 
 *     DataOut<spacedim> embedding_out;
 * 
 *     std::ofstream embedding_out_file("embedding.vtu");
 * 
 *     embedding_out.attach_dof_handler(*space_dh);
 *     embedding_out.add_data_vector(solution, "solution");
 *     embedding_out.build_patches(
 *       parameters.embedding_space_finite_element_degree);
 *     embedding_out.write_vtu(embedding_out_file);
 * 
 * @endcode
 * 
 * The only difference between the two output routines is that in the
 * second case, we want to output the data on the current configuration, and
 * not on the reference one. This is possible by passing the actual
 * embedded_mapping to the DataOut::build_patches function. The mapping will
 * take care of outputting the result on the actual deformed configuration.
 * 

 * 
 * 
 * @code
 *     DataOut<dim, spacedim> embedded_out;
 * 
 *     std::ofstream embedded_out_file("embedded.vtu");
 * 
 *     embedded_out.attach_dof_handler(*embedded_dh);
 *     embedded_out.add_data_vector(lambda, "lambda");
 *     embedded_out.add_data_vector(embedded_value, "g");
 *     embedded_out.build_patches(*embedded_mapping,
 *                                parameters.embedded_space_finite_element_degree);
 *     embedded_out.write_vtu(embedded_out_file);
 *   }
 * 
 * @endcode
 * 
 * Similar to all other tutorial programs, the `run()` function simply calls
 * all other methods in the correct order. Nothing special to note, except
 * that we check if parsing was done before we actually attempt to run our
 * program.
 * 
 * @code
 *   template <int dim, int spacedim>
 *   void DistributedLagrangeProblem<dim, spacedim>::run()
 *   {
 *     AssertThrow(parameters.initialized, ExcNotInitialized());
 *     deallog.depth_console(parameters.verbosity_level);
 * 
 *     setup_grids_and_dofs();
 *     setup_coupling();
 *     assemble_system();
 *     solve();
 *     output_results();
 *   }
 * } // namespace Step60
 * 
 * 
 * 
 * int main(int argc, char **argv)
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       using namespace Step60;
 * 
 *       const unsigned int dim = 1, spacedim = 2;
 * 
 * @endcode
 * 
 * Differently to what happens in other tutorial programs, here we use
 * ParameterAcceptor style of initialization, i.e., all objects are first
 * constructed, and then a single call to the static method
 * ParameterAcceptor::initialize is issued to fill all parameters of the
 * classes that are derived from ParameterAcceptor.
 *       

 * 
 * We check if the user has specified a parameter file name to use when
 * the program was launched. If so, try to read that parameter file,
 * otherwise, try to read the file "parameters.prm".
 *       

 * 
 * If the parameter file that was specified (implicitly or explicitly)
 * does not exist, ParameterAcceptor::initialize will create one for you,
 * and exit the program.
 * 

 * 
 * 
 * @code
 *       DistributedLagrangeProblem<dim, spacedim>::Parameters parameters;
 *       DistributedLagrangeProblem<dim, spacedim>             problem(parameters);
 * 
 *       std::string parameter_file;
 *       if (argc > 1)
 *         parameter_file = argv[1];
 *       else
 *         parameter_file = "parameters.prm";
 * 
 *       ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");
 *       problem.run();
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
 *   return 0;
 * }
 * @endcode
examples/step-60/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序的目录中默认不包含参数文件。另一方面，这个程序想从一个叫parameters.prm的文件中读取它的参数 -- 因此，当你第一次执行它时，你会得到一个异常，即找不到这样的文件。

@code


----------------------------------------------------
Exception on processing:


--------------------------------------------------------
An error occurred in line <74> of file <../source/base/parameter_acceptor.cc> in function
    static void dealii::ParameterAcceptor::initialize(const std::string &, const std::string &, const ParameterHandler::OutputStyle, dealii::ParameterHandler &)
The violated condition was:
    false
Additional information:
    You specified <parameters.prm> as input parameter file, but it does not exist. We created it for you.


--------------------------------------------------------


Aborting!


----------------------------------------------------
@endcode



然而，正如错误信息已经指出的那样，触发该异常的代码也将生成一个参数.prm文件，该文件仅仅包含该程序所关心的所有参数的默认值。通过对参数文件的检查，我们看到以下内容。

@code
# Listing of Parameters
# ---------------------
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0, 1, 2, 3
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10



  subsection Embedded configuration
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  = R=.3, Cx=.4, Cy=.4                 # default:


    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy # default: 0


    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end


  subsection Embedded value
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  =


    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = 1     # default: 0


    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end


  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 1000   # default: 100
    set Reduction     = 1.e-12 # default: 1.e-2
    set Tolerance     = 1.e-12 # default: 1.e-10
  end


end
@endcode



如果你现在运行该程序，你将得到一个名为`used_parameters.prm`的文件，其中包含上述参数的简短版本（没有注释和文档），记录了所有用于运行你的程序的参数。

@code
# Parameter file generated with
# DEAL_II_PACKAGE_VERSION = 9.0.0
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0, 1, 2, 3
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.4, Cy=.4
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = 1
    set Variable names      = x,y,t
  end
  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 1000
    set Reduction     = 1.e-12
    set Tolerance     = 1.e-12
  end
end
@endcode



首先创建`parameters.prm`文件（第一次运行程序），然后创建`used_parameters.prm`（每隔一段时间运行程序），其理由是你可能想让大多数参数保持默认值，而只修改其中的一小部分。

例如，你可以在这个教程程序中使用以下（完全有效的）参数文件。

@code
subsection Distributed Lagrange<1,2>
  set Initial embedded space refinement            = 7
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 3
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.4, Cy=.4
    set Function expression = R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = 1
    set Variable names      = x,y,t
  end
end
@endcode



你会得到与下面测试案例1完全相同的结果。

<a name="Testcase1"></a><h3> Test case 1: </h3>


对于默认问题， $u$ 在 $\Gamma$ 上的值被设置为常数 $1$ ：这就像在 $\Gamma$ 上施加了一个常数迪里希特边界条件，被视为 $\Omega$ 在 $\Gamma$ 内的部分的边界。同样，在 $\partial
\Omega$ 上，我们有零的迪里切特边界条件。


<div class="twocolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.1_no_grid.png" alt = "" width="500"> </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.1_grid.png" alt = "" width="500"> </div> </div></div>

程序的输出将如下所示。

@code
DEAL::Embedded dofs: 129
DEAL::Embedding minimal diameter: 0.0110485, embedded maximal diameter: 0.00781250, ratio: 0.707107
DEAL::Embedding dofs: 2429
DEAL:cg::Starting value 0.166266
DEAL:cg::Convergence step 108 value 7.65958e-13



+---------------------------------------------+------------+------------+
| Total CPU time elapsed since start          |     0.586s |            |
|                                             |            |            |
| Section                         | no. calls |  CPU time  | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble coupling system        |         1 |     0.132s |        23% |
| Assemble system                 |         1 |    0.0733s |        12% |
| Output results                  |         1 |     0.087s |        15% |
| Setup coupling                  |         1 |    0.0244s |       4.2% |
| Setup grids and dofs            |         1 |    0.0907s |        15% |
| Solve system                    |         1 |     0.178s |        30% |
+---------------------------------+-----------+------------+------------+





+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.301s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble coupling system        |         1 |    0.0385s |        13% |
| Assemble system                 |         1 |    0.0131s |       4.3% |
| Output results                  |         1 |    0.0736s |        24% |
| Setup coupling                  |         1 |    0.0234s |       7.7% |
| Setup grids and dofs            |         1 |    0.0679s |        23% |
| Solve system                    |         1 |    0.0832s |        28% |
+---------------------------------+-----------+------------+------------+


@endcode



你可能会注意到，就CPU时间而言，组装耦合系统的成本是组装标准泊松系统的两倍，尽管矩阵更小。这是由于离散化的非匹配性造成的。这是否可以接受，取决于应用。

如果问题被设置在三维环境中，并且浸入式网格与时间有关，那么在每一步重新创建网格要比使用我们这里介绍的技术要昂贵得多。此外，你也许可以在一个均匀细化的正方形或立方体网格上创建一个非常快速和优化的求解器，并在你想进行计算的地方嵌入这里的技术。这就要求你只需要有一个领域的表面代表（一个更便宜和更容易制作的网格）。

为了玩一玩，我们要把虚构的领域以及我们强加给它的边界条件复杂化一点。

<a name="Testcase2and3"></a><h3> Test case 2 and 3: </h3>


如果我们使用以下参数文件。

@code
subsection Distributed Lagrange<1,2>
  set Coupling quadrature order                    = 3
  set Embedded configuration finite element degree = 1
  set Embedded space finite element degree         = 1
  set Embedding space finite element degree        = 1
  set Homogeneous Dirichlet boundary ids           = 0,1,2,3
  set Initial embedded space refinement            = 8
  set Initial embedding space refinement           = 4
  set Local refinements steps near embedded domain = 4
  set Use displacement in embedded interface       = false
  set Verbosity level                              = 10
  subsection Embedded configuration
    set Function constants  = R=.3, Cx=.5, Cy=.5, r=.1, w=12
    set Function expression = (R+r*cos(w*pi*x))*cos(2*pi*x)+Cx; (R+r*cos(w*pi*x))*sin(2*pi*x)+Cy
    set Variable names      = x,y,t
  end
  subsection Embedded value
    set Function constants  =
    set Function expression = x-.5
    set Variable names      = x,y,t
  end
  subsection Schur solver control
    set Log frequency = 1
    set Log history   = false
    set Log result    = true
    set Max steps     = 100000
    set Reduction     = 1.e-12
    set Tolerance     = 1.e-12
  end
end
@endcode



我们得到了一个看起来很 "花 "的域，在这里我们施加了一个线性边界条件  $g=x-.5$  。这个测试表明，该方法在从边界条件中恢复一个完全线性的函数方面实际上是相当准确的，即使网格没有对齐，我们也得到了一个相当好的结果。

用 $2(x-.5)^2-2(y-.5)^2$ 替换 $x-.5$ ，即修改参数文件，使我们有

@code
  ...
  subsection Embedded value
    set Function constants  =
    set Function expression = 2*(x-.5)^2-2*(y-.5)^2
    set Variable names      = x,y,t
  end
@endcode

生产右边的马鞍。

<div class="twocolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.3_grid.png" alt = "" width="500"> </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.4_grid.png" alt = "" width="500"> </div> </div></div>

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Runningwithspacedimequaltothree"></a><h4> Running with `spacedim` equal to three</h4>


虽然目前的教程程序是为`spacedim`等于2而写的，但为了使程序在不同的尺寸组合中运行，你只需要做一些小的改动。

如果你想在`spacedim`等于3，`dim`等于2的情况下运行，那么你几乎肯定要进行以下修改。

- 使用不同的参考域来嵌入网格，也许可以从一个文件中读取它。不可能用一个正方形域的单一参数化来构造一个光滑的封闭表面，因此你很可能想用一个拓扑学上等同于球体边界的参考域。

- 用位移代替变形，将 $\Gamma_0$ 映射为 $\Gamma$ 。

<a name="Moregeneraldomains"></a><h4> More general domains </h4>


我们在其他教程中看到（例如在第5步和第54步）如何从输入文件中读取网格。这个教程程序的一个很好的概括是允许用户选择从参数文件本身读取网格，而不是在教程程序本身硬编码网格类型。

<a name="Preconditioner"></a><h4> Preconditioner</h4>


目前，我们没有关于Schur补数的预处理程序。这对于二维问题来说是可以的，几百次迭代就可以把残差降低到机器的精度，但在三维问题上是行不通的。

在这里，一个好的预处理程序是什么并不明显。我们用舒尔补码解决的物理问题是将Dirichlet数据 $g$ 与Lagrange乘数 $\lambda$ 的值联系起来。   $\lambda$ 可以解释为法线梯度的*跳跃，需要强加在 $u$ 上，跨越 $\Gamma$ ，以获得迪里切特数据 $g$ 。

所以 $S$ 是某种诺伊曼到迪里切特的映射，我们希望有一个迪里切特到诺伊曼映射的良好近似。一种可能性是使用 $\Gamma$ 上的问题的边界元素近似，并构建一个与 $\Gamma$ 相关的泊松问题的超星形算子的粗略近似，这正是迪里切特到诺依曼的映射。

<a name="ParallelCode"></a><h4> Parallel Code </h4>


这里提出的简单代码可以作为更复杂问题的起点，要解决这些问题，需要在并行代码上运行，可能使用分布式网格（见步骤17、步骤40，以及 parallel::shared::Triangulation 和 parallel::distributed::Triangulation). 的文档）。

当并行使用非匹配网格时，会出现一个问题：为了计算矩阵 $C$ ，一个进程需要关于实空间同一部分的两个网格的信息，但是，当使用分布式网格时，这种信息可能无法获得，因为存储在特定处理器上的 $\Omega$ 三角形的本地所有部分可能与存储在同一处理器上的 $\Gamma$ 三角形的本地所有部分不在同一地点。

可以实施各种策略来解决这个问题。

- 分布这两个网格，以便满足这个约束条件。

- 对现实空间中不满足约束的部分使用通信。

- 对嵌入空间使用分布式三角法，对模拟配置使用共享三角法。

后一种策略显然是最容易实现的，因为本教程程序中使用的大多数函数在并行情况下也能不变地工作。当然，我们可以使用反转策略（即有一个分布式嵌入三角法和一个共享嵌入三角法）。

然而，这种策略很可能会更加昂贵，因为根据定义，嵌入网格比嵌入网格要大，而且分配两个网格中最大的网格更有意义，保持最小的网格由所有处理器共享。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-60.cc"
*/
