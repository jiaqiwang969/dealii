/**
@page step_27 The step-27 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Finiteelementcollections">Finite element collections</a>
        <li><a href="#TheDoFHandlerclassinihpiihpimodeassociatingcellswithfiniteelementsandconstraints">The DoFHandler class in <i>hp</i><i>hp</i>-mode, associating cells with finite elements, and constraints</a>
        <li><a href="#Assemblingmatricesandvectorswithhpobjects">Assembling matrices and vectors with hp-objects</a>
        <li><a href="#Asimpleindicatorforhprefinementandestimatingsmoothness">A simple indicator for hp-refinement and estimating smoothness</a>
      <ul>
        <li><a href="#Theidea">The idea</a>
        <li><a href="#Whatwehavetodo">What we have to do</a>
        <li><a href="#Compensatingforanisotropy">Compensating for anisotropy</a>
        <li><a href="#Questionsaboutcellsizes">Questions about cell sizes</a>
      </ul>
        <li><a href="#Complicationswithlinearsystemsforhpdiscretizations">Complications with linear systems for hp-discretizations</a>
      <ul>
        <li><a href="#Creatingthesparsitypattern">Creating the sparsity pattern</a>
        <li><a href="#Eliminatingconstraineddegreesoffreedom">Eliminating constrained degrees of freedom</a>
      </ul>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Themainclass">The main class</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Implementationofthemainclass">Implementation of the main class</a>
      <ul>
        <li><a href="#LaplaceProblemLaplaceProblemconstructor">LaplaceProblem::LaplaceProblem constructor</a>
        <li><a href="#LaplaceProblemLaplaceProblemdestructor">LaplaceProblem::~LaplaceProblem destructor</a>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblempostprocess">LaplaceProblem::postprocess</a>
        <li><a href="#LaplaceProblemcreate_coarse_grid">LaplaceProblem::create_coarse_grid</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Differenthpdecisionstrategies">Different hp-decision strategies</a>
        <li><a href="#Parallelhpadaptivefiniteelements">Parallel hp-adaptive finite elements</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-27/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个教程程序试图展示如何使用 $hp$ -finite element方法与deal.II。它解决的是拉普拉斯方程，因此只建立在前几个教程程序的基础上，特别是步骤4的独立维度编程和步骤6的自适应网格细化。

Babu&scaron;ka和Guo在20世纪80年代初提出了 $hp$ -有限元方法，作为(i)网格细化（即减少有限元计算中的网格参数 $h$ ）或(ii)增加用于形状函数的多项式程度 $p$ 的替代方法。它是基于这样的观察：如果解足够光滑，增加形状函数的多项式程度可以减少近似误差。另一方面，众所周知，即使对于一般良好的椭圆问题，在边界附近、角落或系数不连续的地方也不能保证较高的规则度；因此，在这些地方不能通过增加多项式度来改善近似  $p$  ，而只能通过细化网格，即减少网格大小  $h$  。这些减少误差的不同手段导致了 $hp$ 有限元的概念，即在解足够光滑的地方，近似的有限元空间被调整为具有高多项式度 $p$ ，而在解缺乏规则性的地方，网格宽度 $h$ 被减少。在关于这种方法的第一篇论文中已经意识到， $hp$ -无限元素可以是一个强有力的工具，它可以保证误差不仅以自由度数的某个负数减少，而且实际上是以指数形式减少。

为了实现这个方法，我们需要一些高于一般有限元程序所需的东西，特别是高于我们在步骤6之前的教程程序中所介绍的东西。特别是，我们将不得不讨论以下几个方面。   <ul>   <li>  我们现在不是在所有单元上使用相同的有限元，而是要有一个有限元对象的集合，并将每个单元与这个集合中的一个对象相关联。 </li> 

    <li> 然后，自由度将不得不根据与该特定单元相关的有限元来分配给每个单元。约束条件必须以与悬挂节点相同的方式生成，但我们现在还必须处理两个相邻单元分配不同有限元的情况。 </li> 

    <li>  我们将需要能够集合单元和面对全局矩阵和右手向量的贡献。 </li> 

    <li> 在求解所得到的线性系统后，我们将想分析该解决方案。特别是，我们要计算误差指标，告诉我们是否应该细化一个给定的单元和/或是否应该增加对其使用的形状函数的多项式程度。 </li>   </ul> 。

我们将在本介绍的以下几个小节中讨论所有这些方面。这些任务中的大部分已经由deal.II提供的功能很好地支持了，我们只需要提供程序应该做什么的逻辑，而不是确切地说明这一切将如何发生，这不会让人感到很惊讶。

在deal.II中， $hp$ 的功能大部分被打包到hp-namespace中。这个命名空间提供了处理 $hp$ 分解、集合矩阵和向量以及其他任务的类。我们将在下文中进一步了解它们中的许多。此外，DoFTools和VectorTools命名空间中的大多数函数除了接受非 $hp$ 的对象外，还接受 $hp$ 的对象。许多 $hp$ 的实现也在 @ref hp 文档模块和那里的链接中讨论。

也许值得在这第一部分介绍的最后给出一个稍大的观点。   $hp$ -功能已经在许多不同的有限元软件包中实现（例如，见 @ref hp_paper "hp-paper "中引用的参考文献列表）。然而，总的来说，这些软件包大多只在(i)2d情况下和/或(ii)不连续的Galerkin方法中实现了它。后者是一个重要的简化，因为根据定义，不连续的有限元不要求单元之间的面的连续性，因此，当不同的多项式程度的有限元在一个共同的面相遇时，不需要特殊处理。相比之下，deal.II实现了最普遍的情况，即它允许1d、2d和3d的连续和不连续元素，并自动处理由此产生的复杂性。特别是，它可以处理不同程度的元素在一个面或边缘相遇时的约束（类似于悬挂节点约束）的计算。许多必要的算法和数据结构技术在 @ref hp_paper "hp-paper "中描述，供那些对这些细节感兴趣的人参考。

我们希望，提供这样一个通用的实现方式将有助于进一步探索 $hp$ -方法的潜力。




<a name="Finiteelementcollections"></a><h3>Finite element collections</h3>


现在再来看看如何在deal.II中使用 $hp$ -功能的细节。我们要处理的第一个方面是，现在我们不再只有一个用于所有单元的有限元，而是有许多不同的单元可以选择使用的元素。为此，deal.II引入了<i>finite element
collection</i>的概念，在 hp::FECollection. 类中实现。实质上，这样的集合就像 <code>std::vector@<FiniteElement@></code> 类型的对象，但多了一些小功能和内存管理，更适合手头的任务。正如我们以后所看到的，我们还将使用类似的正交集合，以及&mdash; 尽管我们在这里没有使用它们&mdash; 还有映射集合的概念。所有这些类都在 @ref hpcollection 概述中有所描述。

在这个教程程序中，我们将使用阶数为2到7（在2D中）或2到5（在3D中）的连续Lagrange元素。然后可以按以下方式创建所用元素的集合。

@code
  hp::FECollection<dim> fe_collection;
  for (unsigned int degree = 2; degree <= max_degree; ++degree)
    fe_collection.push_back(FE_Q<dim>(degree));
@endcode






<a name="TheDoFHandlerclassinihpiihpimodeassociatingcellswithfiniteelementsandconstraints"></a><h3>The DoFHandler class in <i>hp</i><i>hp</i>-mode, associating cells with finite elements, and constraints</h3> 。


我们要考虑的下一个任务是如何处理我们要使用的有限元对象的列表。在之前的教程程序中，从第2步开始，我们已经看到DoFHandler类负责在网格（由三角形对象描述）和有限元之间建立联系，为网格的每个顶点、面、边和单元分配正确的自由度数。

这里的情况有点复杂，因为我们不只是有一个单一的有限元对象，而是可能想在不同的单元上使用不同的元素。因此，我们需要两样东西。(i) 一个能够处理这种情况的DoFHandler类的版本，以及(ii) 告诉DoFHandler在哪个单元上使用哪个元素的方法。

这两件事中的第一件是在DoFHandler类的<i>hp</i>模式下实现的：不是将其与一个三角形和一个有限元对象相关联，而是与一个三角形和一个有限元集合相关联。第二部分是通过在DoFHandler的所有单元上的循环来实现的，并为每个单元设置集合中的有限元的索引，该索引将被用于该单元。我们将集合中用于某个单元的有限元对象的索引称为该单元的<i>active FE
index</i>，以表明这是活跃在该单元上的有限元，而集合中的所有其他元素在该单元上是不活跃的。其大致内容是这样的。

@code
  DoFHandler<dim> dof_handler(triangulation);
  for (auto &cell: dof_handler.active_cell_iterators())
    cell->set_active_fe_index(...);
  dof_handler.distribute_dofs(fe_collection);
@endcode



调用 <code>set_active_fe_index()</code> 中的点表示我们以后必须有某种策略来决定在哪个单元格上使用哪个元素；我们以后会再来讨论这个。这里的重点是，这个代码片段的第一行和最后一行与非  $hp$  的情况几乎完全相同。

另一个复杂的情况是，这次我们不只是有来自局部网格细化的悬挂节点，我们还必须处理这样的情况：如果有两个具有不同活动有限元指数的单元在一个面上相遇（例如一个Q2和一个Q3单元），那么我们必须计算有限元场上的额外约束，以确保它是连续的。这在概念上与我们计算悬挂节点约束的方式非常相似，事实上，代码看起来也完全一样。

@code
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
@endcode

换句话说， DoFTools::make_hanging_node_constraints 不仅处理悬挂节点约束，而且同时处理 $hp$ -约束。




<a name="Assemblingmatricesandvectorswithhpobjects"></a><h3>Assembling matrices and vectors with hp-objects</h3>


在这之后，我们必须为正确大小的线性系统设置矩阵和向量，并将它们组合起来。设置它们的方法与非 $hp$ 的情况完全相同。组装则需要更多的思考。

当然，主要的想法是不变的：我们必须在所有单元中循环，集合局部贡献，然后将它们复制到全局对象中。正如在第3步中详细讨论的那样，deal.II有一个FEValues类，它将有限元描述、映射和正交公式拉到一起，帮助评估形状函数的值和梯度，以及映射到单元实际位置的每个正交点的其他信息。每当我们转到一个新的单元时，我们就会重新初始化这个FEValues对象，从而要求它重新计算从一个单元到另一个单元的那部分信息。然后，它可以被用来总结本地对双线性形式和右手边的贡献。

在 $hp$ 有限元方法的背景下，我们必须处理这样一个事实：我们不在每个单元上使用相同的有限元对象。事实上，我们甚至不应该对所有单元使用相同的正交对象，而应该对使用高阶有限元的单元使用高阶正交公式。同样地，我们可能也想在这些单元上使用高阶映射。

为了方便这些考虑，deal.II有一个类 hp::FEValues ，可以做我们在当前背景下需要的事情。不同的是，它不是一个单一的有限元、正交公式和映射，而是这些对象的集合。它的使用非常类似于常规的FEValues类，也就是说，在所有单元中循环的有趣部分看起来像这样。

@code
  hp::FEValues<dim> hp_fe_values(mapping_collection,
                                 fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      hp_fe_values.reinit(cell);


      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();


      ...  // assemble local contributions and copy them into global object
    }
@endcode



在这个教程程序中，我们将始终使用Q1映射，所以 hp::FEValues 结构的映射集合参数将被省略。在循环中，我们首先初始化当前单元的 hp::FEValues 对象。第二个、第三个和第四个参数表示我们希望在该单元上使用的正交、映射和有限元对象在各自集合中的索引。这些参数可以省略（在下面的程序中也是如此），在这种情况下， <code>cell-@>active_fe_index()</code> 被用于这个索引。之所以这样选择这些参数的顺序，是因为有时人们可能想从各自的集合中挑选不同的正交或映射对象，但几乎不可能从这个单元上使用的有限元中挑选不同的有限元，即索引与 <code>cell-@>active_fe_index()</code> 不同。因此，有限元集合的索引是最后一个默认参数，这样可以方便地省略它。

这个 <code>reinit</code> 调用的作用如下： hp::FEValues 类检查它之前是否已经为这个有限元、正交和映射对象的组合分配了一个非 $hp$ -FEValues对象。如果没有，它将分配一个。然后，它为当前单元重新初始化这个对象，之后，现在有一个FEValues对象可用于当前单元上选定的有限元、正交和映射。然后通过调用 <code>hp_fe_values.get_present_fe_values()</code> 获得对该对象的引用，并将以通常的方式用于集合本地贡献。




<a name="Asimpleindicatorforhprefinementandestimatingsmoothness"></a><h3>A simple indicator for hp-refinement and estimating smoothness</h3>


自适应有限元方法的核心之一是，我们用一个指标来检查计算出的解（后验），告诉我们哪些是误差最大的单元，然后再对它们进行细化。在其他许多教程程序中，我们使用KellyErrorEstimator类来获得一个单元上的误差大小的指示，尽管我们也在一些程序中讨论了更复杂的策略，最重要的是在步骤14。

在任何情况下，只要决定只是 "细化这个单元 "或 "不细化这个单元"，实际的细化步骤就不是特别具有挑战性。然而，在这里，我们有一个能够进行hp细化的代码，也就是说，每当我们检测到某个单元上的误差太大，我们突然有两个选择：我们可以通过把它分割成几个小单元来细化这个单元，或者我们可以增加在它上面使用的形状函数的多项式程度。我们如何知道哪个是更有希望的策略？回答这个问题是本文写作时 $hp$ -无限元研究的核心问题。

简而言之，这个问题目前在文献中似乎还没有解决。有许多或多或少复杂的方案来解决这个问题，但没有任何方案像KellyErrorEstimator那样被普遍接受为一个好的、即使不是最佳的误差指标。大多数建议采用这样的事实：只要解是局部光滑的，增加多项式的度数是有益的，而只要网格是粗糙的，就应该细化。然而，如何确定解的局部光滑度以及决定一个解何时光滑到允许增加 $p$ 的问题无疑是很大很重要的问题。

在下文中，我们提出了一个简单的解决方案的局部平滑性的估计方法。正如我们将在结果部分看到的，这个估计器有缺陷，特别是就有局部悬空节点的单元而言。因此，我们不打算把下面的想法作为问题的完整解决方案。相反，它是作为一个值得进一步研究和调查的想法来处理的。换句话说，我们不打算在关于一般问题的答案的争论中进入一个复杂的建议。然而，为了证明我们对 $hp$ -无限元素的方法，我们需要一个简单的指标，它确实产生一些有用的信息，能够驱动本教程程序将进行的简单计算。




<a name="Theidea"></a><h4>The idea</h4>


我们在这里的方法很简单：对于一个函数 $u({\bf x})$ 来说，它是在单元格 $K$ 上的索博列夫空间 $H^s(K)$ 中，它必须满足条件

@f[
   \int_K |\nabla^s u({\bf x})|^2 \; d{\bf x} < \infty.


@f]

假设单元格 $K$ 不是退化的，即从单元格到单元格 $K$ 的映射足够规则，上述条件当然等同于

@f[
   \int_{\hat K} |\nabla^s \hat u(\hat{\bf x})|^2 \; d\hat{\bf x} < \infty\,,


@f]

其中 $\hat u(\hat{\bf x})$ 是映射回单元格 $\hat K$ 的函数 $u({\bf x})$  。从这里，我们可以做以下工作：首先，让我们定义 $\hat u$ 的傅里叶级数为

@f[
   \hat u(\hat{\bf x})
   = \sum_{\bf k} \hat U_{\bf k}\,e^{-i {\bf k}\cdot \hat{\bf x}},


@f]

傅里叶向量 ${\bf k}=(k_x,k_y)$ 在2d中， ${\bf k}=(k_x,k_y,k_z)$ 在3d中，等等，以及 $k_x,k_y,k_z=0,2\pi,4\pi,\ldots$  。扩张 $\hat U_{\bf k}$ 的系数可以用 $L^2$ 得到--指数基的正交性

@f[
\int_{\hat K} e^{-i {\bf m}\cdot \hat{\bf x}} e^{i {\bf n}\cdot \hat{\bf x}} d\hat{\bf x} = \delta_{\bf m \bf n},


@f]

导致以下表达式

@f[
   \hat U_{\bf k}
   = \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat u(\hat{\bf x}) d\hat{\bf x} \,.


@f]

很明显，我们可以将 $H^s$ 的 $\hat u$ 准则写成

@f[
  \int_{\hat K} |\nabla^s \hat u(\hat{\bf x})|^2 \; d\hat{\bf x}
  =
  \int_{\hat K}
  \left|
    \sum_{\bf k} |{\bf k}|^s e^{-i{\bf k}\cdot \hat{\bf x}} \hat U_{\bf k}
  \right|^2 \; d\hat{\bf x}
  =
  \sum_{\bf k}
    |{\bf k}|^{2s}
    |\hat U_{\bf k}|^2.


@f]

换句话说，如果这个规范是有限的（即，对于 $\hat u(\hat{\bf x})$ 在 $H^s(\hat K)$ 中），我们需要

@f[
   |\hat U_{\bf k}| = {\cal O}\left(|{\bf k}|^{-\left(s+1/2+\frac{d-1}{2}+\epsilon\right)}\right).


@f]

换句话说：我们想要的规则性 $s$ 越高，傅里叶系数归零的速度就越快。如果你想知道额外的指数 $\frac{d-1}2$ 从何而来：我们想利用这样一个事实： $\sum_l a_l < \infty$ 如果序列 $a_l =
{\cal O}(l^{-1-\epsilon})$ 对于任何 $\epsilon>0$ 。问题是，我们在这里不仅有一个单一变量的求和，而且有位于 $d$ 维球内的 $2\pi$ 的所有整数倍的求和，因为我们有向量成分 $k_x, k_y,
\ldots$  。就像我们通过用整个直线上的积分代替总和来证明上面的序列 $a_l$ 收敛一样，我们可以用 $d$ -维空间上的积分来代替我们的 $d$ -维总和。现在我们必须注意到，在距离 $|{\bf k}|$ 和 $|{\bf k}|+d|{\bf k}|$ 之间，存在着多达一个常数的 $|{\bf k}|^{d-1}$ 模式，这与我们可以将体积元素 $dx\;dy$ 转化为 $2\pi r\; dr$ 的方式相同。因此，不再是 $|{\bf k}|^{2s}|\hat
U_{\bf k}|^2$ 必须衰变为 ${\cal O}(|{\bf k}|^{-1-\epsilon})$ ，而实际上是 $|{\bf k}|^{2s}|\hat U_{\bf k}|^2 |{\bf k}|^{d-1}$ 。指数的比较产生了结果。

我们可以把这个问题转过来。假设我们得到了一个未知平滑度的函数 $\hat u$ 。让我们计算它的傅里叶系数 $\hat U_{\bf k}$ ，看看它们衰减的速度。如果它们的衰减速度为

@f[
   |\hat U_{\bf k}| = {\cal O}(|{\bf k}|^{-\mu-\epsilon}),


@f]

因此，我们这里的函数是在 $H^{\mu-d/2}$ 。




<a name="Whatwehavetodo"></a><h4>What we have to do</h4>


那么，我们要做什么来估计 $u({\bf x})$ 在单元格 $K$ 上的局部光滑度呢？显然，第一步是计算我们解决方案的傅里叶系数。傅里叶级数是无限级数，我们通过只计算级数的前几项来简化我们的任务，例如， $|{\bf k}|\le 2\pi N$ 有一个截止点 $N$ 。让我们顺便说一下，我们希望选择 $N$ 足够大，这样我们至少可以捕获那些变化最大的形状函数的变化。另一方面，我们不应该把 $N$ 选得太大：显然，一个有限元函数，作为一个多项式，在任何给定的单元上都在 $C^\infty$ 中，所以系数将不得不在一个点上指数衰减；由于我们想估计这个多项式所近似的函数的平稳性，而不是多项式本身，我们需要为 $N$ 选择一个合理的截止点。无论怎样，计算这个数列并不特别困难：从定义上看

@f[
   \hat U_{\bf k}
   = \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat u(\hat{\bf x}) d\hat{\bf x}


@f]

我们看到，我们可以计算系数 $\hat U_{\bf k}$ 为

@f[
   \hat U_{\bf k}
   =
   \sum_{i=0}^{\textrm{dofs per cell}}
   \left[\int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat \varphi_i(\hat{\bf x})
   d\hat{\bf x} \right] u_i,


@f]

其中 $u_i$ 是这个单元上 $i$ 个自由度的值。换句话说，我们可以把它写成一个矩阵-向量乘积

@f[
   \hat U_{\bf k}
   = {\cal F}_{{\bf k},j} u_j,


@f]

与矩阵

@f[
   {\cal F}_{{\bf k},j}
   =
   \int_{\hat K} e^{i {\bf k}\cdot \hat{\bf x}} \hat \varphi_j(\hat{\bf x}) d\hat{\bf x}.


@f]

对于给定数量的形状函数  $\varphi_j$  和傅里叶模式  $N$  ，这个矩阵很容易计算出来。因此，寻找系数 $\hat U_{\bf k}$ 是一个相当琐碎的工作。为了进一步简化我们的生活，我们将使用 FESeries::Fourier 类，它正是这样做的。

接下来的任务是，我们必须估计这些系数随 $|{\bf k}|$ 衰减的速度。问题是，当然，我们首先只有有限的这些系数。换句话说，我们能做的最好的事情是将一个函数 $\alpha |{\bf k}|^{-\mu}$ 拟合到我们的数据点 $\hat U_{\bf k}$ ，例如通过最小二乘法程序确定 $\alpha,\mu$ 。

@f[
   \min_{\alpha,\mu}
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( |\hat U_{\bf k}| - \alpha |{\bf k}|^{-\mu}\right)^2


@f]

然而，这样做的问题是，它导致了一个非线性问题，这是我们想避免的事实。另一方面，如果我们试图将我们的系数的对数与 $\alpha |{\bf k}|^{-\mu}$ 的对数相适应，我们可以将问题转化为一个更简单的问题，就像这样。

@f[
   \min_{\alpha,\mu}
   Q(\alpha,\mu) =
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( \ln |\hat U_{\bf k}| - \ln (\alpha |{\bf k}|^{-\mu})\right)^2.


@f]

利用关于对数的一般事实，我们可以看到，这就产生了一个问题

@f[
   \min_{\beta,\mu}
   Q(\beta,\mu) =
   \frac 12 \sum_{{\bf k}, |{\bf k}|\le N}
   \left( \ln |\hat U_{\bf k}| - \beta + \mu \ln |{\bf k}|\right)^2,


@f]

其中  $\beta=\ln \alpha$  。现在这是一个问题，对于这个问题，最优性条件  $\frac{\partial Q}{\partial\beta}=0,
\frac{\partial Q}{\partial\mu}=0$  , 在  $\beta,\mu$  中是线性的。我们可以把这些条件写成如下。

@f[
   \left(\begin{array}{cc}
   \sum_{{\bf k}, |{\bf k}|\le N} 1 &
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|
   \\
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}| &
   \sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2
   \end{array}\right)
   \left(\begin{array}{c}
   \beta \\ -\mu
   \end{array}\right)
   =
   \left(\begin{array}{c}
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|
   \\
   \sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}|
   \end{array}\right)


@f]

这个线性系统很容易被倒置，从而得到

@f[
   \beta =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }


@f]

和

@f[
   \mu =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }.


@f]



这无非是线性回归拟合，为了做到这一点，我们将使用 FESeries::linear_regression(). 虽然我们对 $\beta$ 的实际值不是特别感兴趣，但上面的公式给了我们一个计算指数 $\mu$ 的平均值，然后我们可以用来确定 $\hat u(\hat{\bf x})$ 与 $s=\mu-\frac d2$ 在一起。

上面概述的这些步骤适用于许多不同的场景，这促使我们在deal.II中引入了一个通用函数 SmoothnessEstimator::Fourier::coefficient_decay() ，该函数将本节描述的所有任务结合在一个简单的函数调用中。我们将在本程序的实现中使用它。




<a name="Compensatingforanisotropy"></a><h4>Compensating for anisotropy</h4>


在上面的公式中，我们已经得出了傅里叶系数 $\hat U_{\bf
k}$  。因为 ${\bf k}$ 是一个矢量，对于相同的绝对值 $|{\bf k}|$ ，我们会得到许多傅里叶系数 $\hat U_{{\bf k}}$ ，对应于不同方向的傅里叶变换。如果我们现在考虑一个像 $|x|y^2$ 这样的函数，那么我们会发现在 $x$ 方向有很多大的傅里叶系数，因为这个方向的函数是不平滑的，但在 $y$ 方向有快速衰减的傅里叶系数，因为那里的函数是平滑的。由此产生的问题是：如果我们简单地将我们的多项式衰减 $\alpha |{\bf k}|^\mu$ 与<i>all</i>的傅里叶系数拟合，我们将把它拟合成一个光滑度<i>averaged in all spatial directions</i>。这就是我们想要的吗？还是只考虑所有 ${\bf k}$ 中幅度最大的系数 $\hat U_{{\bf k}}$ ，本质上是想确定解在那个空间方向上的平滑度，在这个方向上解显得最粗糙？

人们也许可以为这两种情况争辩。如果deal.II有能力使用各向异性的有限元，即在不同的空间方向使用不同的多项式度数的有限元，那么这个问题将更有意义，因为它们能够更好地利用方向上的可变平滑度。唉，在编写这个教程程序时，这种能力并不存在。

无论怎样，由于我们只有同位素的有限元类，我们采取的观点是，我们应该将多项式的程度调整到最低的规则性，以保持低的数值努力。因此，不使用公式

@f[
   \mu =
   \frac
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}|\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |\hat U_{{\bf k}}| \ln |{\bf k}| \right)
   }
   {
   \left(\sum_{{\bf k}, |{\bf k}|\le N} 1\right)
   \left(\sum_{{\bf k}, |{\bf k}|\le N} (\ln |{\bf k}|)^2\right)


   -
   \left(\sum_{{\bf k}, |{\bf k}|\le N} \ln |{\bf k}|\right)^2
   }.


@f]

为了计算如上所示的 $\mu$ ，我们必须对所有的和稍作修改：不是对所有的傅里叶模式求和，而是只对那些傅里叶系数是所有 $\hat U_{{\bf k}}$ 中最大的、具有相同幅度的 $|{\bf k}|$ 求和，也就是说，上面的所有和必须由以下的和来代替。

@f[
  \sum_{{\bf k}, |{\bf k}|\le N}
  \longrightarrow
  \sum_{\begin{matrix}{{\bf k}, |{\bf k}|\le N} \\ {|\hat U_{{\bf k}}| \ge |\hat U_{{\bf k}'}|
  \ \textrm{for all}\ {\bf k}'\ \textrm{with}\ |{\bf k}'|=|{\bf k}|}\end{matrix}}.


@f]

这是我们将在计划中实施的形式。




<a name="Questionsaboutcellsizes"></a><h4>Questions about cell sizes</h4>


有人可能会问，我们只在解的<i>reference cell</i>（而不是实数单元）上计算傅里叶变换，这是否是一个问题？毕竟，在变换过程中，我们将解决方案拉伸了一个系数 $\frac 1h$ ，从而使傅里叶频率移动了一个系数 $h$ 。这是一个特别值得关注的问题，因为我们可能有相邻的单元，其网格大小 $h$ 相差2倍，如果其中一个单元比另一个更精细。这个问题也是出于这样的考虑：正如我们在下面的结果部分所看到的，估计的解决方案的平滑度应该是一个或多或少的连续函数，但在网格大小跳跃的地方表现出跳跃。因此，我们似乎很自然地要问，我们是否必须对这种转换进行补偿。

简短的回答是 "不"。在上述过程中，我们试图找到系数 $\beta,\mu$ ，使条款的平方之和最小。

@f[
   \ln |\hat U_{{\bf k}}| - \beta + \mu \ln |{\bf k}|.


@f]

补偿变换意味着不试图拟合相对于傅里叶频率 ${\bf k}$ <i>on the unit
cell</i>的衰减 $|{\bf k}|^\mu$ ，而是拟合在参考单元<i>to the Fourier frequencies on the real cell $|\bf
k|h$</i>上计算的系数 $\hat U_{{\bf k}}$ ，其中 $h$ 是变换算子的规范（即类似单元直径的东西）。换句话说，我们将不得不最小化条款的平方之和

@f[
   \ln |\hat U_{{\bf k}}| - \beta + \mu \ln (|{\bf k}|h).


@f]

来代替。然而，利用对数的基本属性，这只是相当于最小化了

@f[
   \ln |\hat U_{{\bf k}}| - (\beta - \mu \ln h) + \mu \ln (|{\bf k}|).


@f]

换句话说，这个问题和原来的最小二乘法问题将产生相同的最佳拟合指数 $\mu$ ，尽管偏移量在一种情况下是 $\beta$ ，在另一种情况下是 $\beta-\mu \ln h$  。然而，由于我们对偏移量根本不感兴趣，而只对指数感兴趣，所以我们是否对傅里叶频率进行缩放以考虑网格大小的影响并不重要，在两种情况下估计的平滑度指数都是一样的。




<a name="Complicationswithlinearsystemsforhpdiscretizations"></a><h3>Complications with linear systems for hp-discretizations</h3>


<a name="Creatingthesparsitypattern"></a><h4>Creating the sparsity pattern</h4>


 $hp$ -方法的问题之一是，形状函数的高多项式程度与大量受限自由度一起导致矩阵的某些行有大量非零条目。同时，因为有些地方我们使用的是低多项式度，因此矩阵行的非零项相对较少。因此，为这些矩阵分配稀疏性模式是一个挑战：我们不能简单地从带宽的估计开始组装一个SparsityPattern，而不使用大量的额外内存。

我们为底层线性系统创建SparsityPattern的方式与我们用来执行约束的策略紧密相连。 deal.II支持以两种方式处理线性系统中的约束。<ol>  <li>  在不考虑约束条件的情况下组装矩阵，之后用 AffineConstraints::condense, 或 </li> 来应用约束条件  <li>  在我们用 AffineConstraints::distribute_local_to_global.</li>  </ol> 来组装系统时应用约束条件。然后系统矩阵使用从DynamicSparsityPattern复制过来的SparsityPattern。这个方法在步骤2中解释，并在大多数教程程序中使用。

早期的教程程序使用一阶或二阶有限元，因此去除与受限自由度相对应的稀疏模式中的条目不会对矩阵明确存储的零的总体数量产生很大影响。然而，由于多达三分之一的自由度在hp微分中可能受到约束（对于高阶元素，这些约束可以将一个自由度与多达10个或20个其他自由度相联系），值得考虑这些约束，因为所产生的矩阵将更加稀疏（因此，矩阵-向量乘积或因子化也将大大加快）。




<a name="Eliminatingconstraineddegreesoffreedom"></a><h4>Eliminating constrained degrees of freedom</h4>


 $hp$ 方法的第二个问题是，我们有如此多的受限自由度：通常有三分之一的自由度（在三维中）是受限的，因为它们要么属于有悬空节点的单元，要么位于与具有更高或更低多项式度的单元相邻的单元上。事实上，这并不比非 $hp$ 模式中受约束自由度的比例高多少，但不同的是，每个受约束的悬空节点不仅受制于相邻的两个自由度，而且还受制于更多的自由度。

事实证明，在步骤6中首先提出的在用 AffineConstraints::distribute_local_to_global 计算元素矩阵和向量时消除约束的策略，对于这种情况也是最有效的方法。另一种策略是首先建立没有约束的矩阵，然后 "浓缩 "掉有约束的自由度，这种策略要昂贵得多。事实证明，用这种低效率的算法建立稀疏模式至少需要 ${\cal O}(N \log N)$ 个未知数，而理想的有限元程序当然只有与未知数成线性的算法。对稀疏模式的创建以及矩阵的装配进行计时显示，步骤6中提出的算法（并在下面的代码中使用）确实更快。

在我们的程序中，我们也将把边界条件作为（可能是不均匀的）约束条件，并把矩阵的行和列也消除。为此我们要做的就是在设置阶段调用插值Dirichlet边界条件的函数，以便告诉AffineConstraints对象关于它们的情况，然后同时在矩阵和向量上做从局部到全局的数据转移。这正是我们在步骤6中所展示的。




<a name="Thetestcase"></a><h3>The test case</h3>


我们要用这个程序解决的测试案例是我们在第14步中已经看过的一个案例的重考：我们要解决拉普拉斯方程

@f[


   -\Delta u = f


@f]

2d，与 $f=(x+1)(y+1)$ ，以及 $u$ 的零Dirichlet边界值。我们在域 $[-1,1]^2\backslash[-\frac 12,\frac 12]^2$ 上这样做，即一个中间有一个方孔的正方形。

当然，与第14步的不同之处在于，我们使用 $hp$ -无限元素来求解。这个测试案例是有意义的，因为它在洞的四角有重入角，在这些地方的解有奇异性。因此，我们期望解在域的内部是平滑的，而在奇点附近是粗糙的。希望我们的细化和光滑度指标能够看到这种行为，并在远离奇点的地方细化网格，同时提高多项式的度数。正如我们将在结果部分看到的，情况确实如此。


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
 * The first few files have already been covered in previous examples and will
 * thus not be further commented on.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * These are the new files we need. The first and second provide the
 * FECollection and the <i>hp</i> version of the FEValues class as described in
 * the introduction of this program. The next one provides the functionality
 * for automatic $hp$-adaptation, for which we will use the estimation
 * algorithms based on decaying series expansion coefficients that are part of
 * the last two files.
 * 
 * @code
 * #include <deal.II/hp/fe_collection.h>
 * #include <deal.II/hp/fe_values.h>
 * #include <deal.II/hp/refinement.h>
 * #include <deal.II/fe/fe_series.h>
 * #include <deal.II/numerics/smoothness_estimator.h>
 * 
 * @endcode
 * 
 * The last set of include files are standard C++ headers.
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * Finally, this is as in previous programs:
 * 
 * @code
 * namespace Step27
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * The main class of this program looks very much like the one already used
 * in the first few tutorial programs, for example the one in step-6. The
 * main difference is that we have merged the refine_grid and output_results
 * functions into one since we will also want to output some of the
 * quantities used in deciding how to refine the mesh (in particular the
 * estimated smoothness of the solution).
 *   

 * 
 * As far as member variables are concerned, we use the same structure as
 * already used in step-6, but we need collections instead of
 * individual finite element, quadrature, and face quadrature objects. We
 * will fill these collections in the constructor of the class. The last
 * variable, <code>max_degree</code>, indicates the maximal polynomial
 * degree of shape functions used.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem();
 *     ~LaplaceProblem();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void create_coarse_grid();
 *     void postprocess(const unsigned int cycle);
 * 
 *     Triangulation<dim> triangulation;
 * 
 *     DoFHandler<dim>          dof_handler;
 *     hp::FECollection<dim>    fe_collection;
 *     hp::QCollection<dim>     quadrature_collection;
 *     hp::QCollection<dim - 1> face_quadrature_collection;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     const unsigned int max_degree;
 *   };
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
 * Next, let us define the right hand side function for this problem. It is
 * $x+1$ in 1d, $(x+1)(y+1)$ in 2d, and so on.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> &p,
 *                                    const unsigned int /*component*/) const
 *   {
 *     double product = 1;
 *     for (unsigned int d = 0; d < dim; ++d)
 *       product *= (p[d] + 1);
 *     return product;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Implementationofthemainclass"></a> 
 * <h3>Implementation of the main class</h3>
 * 

 * 
 * 
 * <a name="LaplaceProblemLaplaceProblemconstructor"></a> 
 * <h4>LaplaceProblem::LaplaceProblem constructor</h4>
 * 

 * 
 * The constructor of this class is fairly straightforward. It associates
 * the DoFHandler object with the triangulation, and then sets the
 * maximal polynomial degree to 7 (in 1d and 2d) or 5 (in 3d and higher). We
 * do so because using higher order polynomial degrees becomes prohibitively
 * expensive, especially in higher space dimensions.
 *   

 * 
 * Following this, we fill the collections of finite element, and cell and
 * face quadrature objects. We start with quadratic elements, and each
 * quadrature formula is chosen so that it is appropriate for the matching
 * finite element in the hp::FECollection object.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem()
 *     : dof_handler(triangulation)
 *     , max_degree(dim <= 2 ? 7 : 5)
 *   {
 *     for (unsigned int degree = 2; degree <= max_degree; ++degree)
 *       {
 *         fe_collection.push_back(FE_Q<dim>(degree));
 *         quadrature_collection.push_back(QGauss<dim>(degree + 1));
 *         face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemLaplaceProblemdestructor"></a> 
 * <h4>LaplaceProblem::~LaplaceProblem destructor</h4>
 * 

 * 
 * The destructor is unchanged from what we already did in step-6:
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::~LaplaceProblem()
 *   {
 *     dof_handler.clear();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 *   

 * 
 * This function is again a verbatim copy of what we already did in
 * step-6. Despite function calls with exactly the same names and arguments,
 * the algorithms used internally are different in some aspect since the
 * dof_handler variable here is in $hp$-mode.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe_collection);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
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
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system</h4>
 * 

 * 
 * This is the function that assembles the global matrix and right hand side
 * vector from the local contributions of each cell. Its main working is as
 * has been described in many of the tutorial programs before. The
 * significant deviations are the ones necessary for <i>hp</i> finite
 * element methods. In particular, that we need to use a collection of
 * FEValues object (implemented through the hp::FEValues class), and that we
 * have to eliminate constrained degrees of freedom already when copying
 * local contributions into global objects. Both of these are explained in
 * detail in the introduction of this program.
 *   

 * 
 * One other slight complication is the fact that because we use different
 * polynomial degrees on different cells, the matrices and vectors holding
 * local contributions do not have the same size on all cells. At the
 * beginning of the loop over all cells, we therefore each time have to
 * resize them to the correct size (given by <code>dofs_per_cell</code>).
 * Because these classes are implemented in such a way that reducing the size
 * of a matrix or vector does not release the currently allocated memory
 * (unless the new size is zero), the process of resizing at the beginning of
 * the loop will only require re-allocation of memory during the first few
 * iterations. Once we have found in a cell with the maximal finite element
 * degree, no more re-allocations will happen because all subsequent
 * <code>reinit</code> calls will only set the size to something that fits the
 * currently allocated memory. This is important since allocating memory is
 * expensive, and doing so every time we visit a new cell would take
 * significant compute time.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_system()
 *   {
 *     hp::FEValues<dim> hp_fe_values(fe_collection,
 *                                    quadrature_collection,
 *                                    update_values | update_gradients |
 *                                      update_quadrature_points |
 *                                      update_JxW_values);
 * 
 *     RightHandSide<dim> rhs_function;
 * 
 *     FullMatrix<double> cell_matrix;
 *     Vector<double>     cell_rhs;
 * 
 *     std::vector<types::global_dof_index> local_dof_indices;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
 * 
 *         cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 *         cell_matrix = 0;
 * 
 *         cell_rhs.reinit(dofs_per_cell);
 *         cell_rhs = 0;
 * 
 *         hp_fe_values.reinit(cell);
 * 
 *         const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
 * 
 *         std::vector<double> rhs_values(fe_values.n_quadrature_points);
 *         rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values);
 * 
 *         for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
 *              ++q_point)
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             {
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                 cell_matrix(i, j) +=
 *                   (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q)
 *                    fe_values.shape_grad(j, q_point) * // grad phi_j(x_q)
 *                    fe_values.JxW(q_point));           // dx
 * 
 *               cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
 *                               rhs_values[q_point] *               // f(x_q)
 *                               fe_values.JxW(q_point));            // dx
 *             }
 * 
 *         local_dof_indices.resize(dofs_per_cell);
 *         cell->get_dof_indices(local_dof_indices);
 * 
 *         constraints.distribute_local_to_global(
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * The function solving the linear system is entirely unchanged from
 * previous examples. We simply try to reduce the initial residual (which
 * equals the $l_2$ norm of the right hand side) by a certain factor:
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(system_rhs.size(),
 *                                  1e-12 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblempostprocess"></a> 
 * <h4>LaplaceProblem::postprocess</h4>
 * 

 * 
 * After solving the linear system, we will want to postprocess the
 * solution. Here, all we do is to estimate the error, estimate the local
 * smoothness of the solution as described in the introduction, then write
 * graphical output, and finally refine the mesh in both $h$ and $p$
 * according to the indicators computed before. We do all this in the same
 * function because we want the estimated error and smoothness indicators
 * not only for refinement, but also include them in the graphical output.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::postprocess(const unsigned int cycle)
 *   {
 * @endcode
 * 
 * Let us start with computing estimated error and smoothness indicators,
 * which each are one number for each active cell of our
 * triangulation. For the error indicator, we use the KellyErrorEstimator
 * class as always.
 * 
 * @code
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 *     KellyErrorEstimator<dim>::estimate(
 *       dof_handler,
 *       face_quadrature_collection,
 *       std::map<types::boundary_id, const Function<dim> *>(),
 *       solution,
 *       estimated_error_per_cell);
 * 
 * @endcode
 * 
 * Estimating the smoothness is performed with the method of decaying
 * expansion coefficients as outlined in the introduction. We will first
 * need to create an object capable of transforming the finite element
 * solution on every single cell into a sequence of Fourier series
 * coefficients. The SmoothnessEstimator namespace offers a factory function
 * for such a FESeries::Fourier object that is optimized for the process of
 * estimating smoothness. The actual determination of the decay of Fourier
 * coefficients on every individual cell then happens in the last function.
 * 
 * @code
 *     Vector<float> smoothness_indicators(triangulation.n_active_cells());
 *     FESeries::Fourier<dim> fourier =
 *       SmoothnessEstimator::Fourier::default_fe_series(fe_collection);
 *     SmoothnessEstimator::Fourier::coefficient_decay(fourier,
 *                                                     dof_handler,
 *                                                     solution,
 *                                                     smoothness_indicators);
 * 
 * @endcode
 * 
 * Next we want to generate graphical output. In addition to the two
 * estimated quantities derived above, we would also like to output the
 * polynomial degree of the finite elements used on each of the elements
 * on the mesh.
 *     

 * 
 * The way to do that requires that we loop over all cells and poll the
 * active finite element index of them using
 * <code>cell-@>active_fe_index()</code>. We then use the result of this
 * operation and query the finite element collection for the finite
 * element with that index, and finally determine the polynomial degree of
 * that element. The result we put into a vector with one element per
 * cell. The DataOut class requires this to be a vector of
 * <code>float</code> or <code>double</code>, even though our values are
 * all integers, so that is what we use:
 * 
 * @code
 *     {
 *       Vector<float> fe_degrees(triangulation.n_active_cells());
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         fe_degrees(cell->active_cell_index()) =
 *           fe_collection[cell->active_fe_index()].degree;
 * 
 * @endcode
 * 
 * With now all data vectors available -- solution, estimated errors and
 * smoothness indicators, and finite element degrees --, we create a
 * DataOut object for graphical output and attach all data:
 * 
 * @code
 *       DataOut<dim> data_out;
 * 
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "solution");
 *       data_out.add_data_vector(estimated_error_per_cell, "error");
 *       data_out.add_data_vector(smoothness_indicators, "smoothness");
 *       data_out.add_data_vector(fe_degrees, "fe_degree");
 *       data_out.build_patches();
 * 
 * @endcode
 * 
 * The final step in generating output is to determine a file name, open
 * the file, and write the data into it (here, we use VTK format):
 * 
 * @code
 *       const std::string filename =
 *         "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk";
 *       std::ofstream output(filename);
 *       data_out.write_vtk(output);
 *     }
 * 
 * @endcode
 * 
 * After this, we would like to actually refine the mesh, in both $h$ and
 * $p$. The way we are going to do this is as follows: first, we use the
 * estimated error to flag those cells for refinement that have the
 * largest error. This is what we have always done:
 * 
 * @code
 *     {
 *       GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.3,
 *                                                       0.03);
 * 
 * @endcode
 * 
 * Next we would like to figure out which of the cells that have been
 * flagged for refinement should actually have $p$ increased instead of
 * $h$ decreased. The strategy we choose here is that we look at the
 * smoothness indicators of those cells that are flagged for refinement,
 * and increase $p$ for those with a smoothness larger than a certain
 * relative threshold. In other words, for every cell for which (i) the
 * refinement flag is set, (ii) the smoothness indicator is larger than
 * the threshold, and (iii) we still have a finite element with a
 * polynomial degree higher than the current one in the finite element
 * collection, we will assign a future FE index that corresponds to a
 * polynomial with degree one higher than it currently is. The following
 * function is capable of doing exactly this. Absent any better
 * strategies, we will set the threshold via interpolation between the
 * minimal and maximal smoothness indicators on cells flagged for
 * refinement. Since the corner singularities are strongly localized, we
 * will favor $p$- over $h$-refinement quantitatively. We achieve this
 * with a low threshold by setting a small interpolation factor of 0.2. In
 * the same way, we deal with cells that are going to be coarsened and
 * decrease their polynomial degree when their smoothness indicator is
 * below the corresponding threshold determined on cells to be coarsened.
 * 
 * @code
 *       hp::Refinement::p_adaptivity_from_relative_threshold(
 *         dof_handler, smoothness_indicators, 0.2, 0.2);
 * 
 * @endcode
 * 
 * The above function only determines whether the polynomial degree will
 * change via future FE indices, but does not manipulate the
 * $h$-refinement flags. So for cells that are flagged for both refinement
 * categories, we prefer $p$- over $h$-refinement. The following function
 * call ensures that only one of $p$- or $h$-refinement is imposed, and
 * not both at once.
 * 
 * @code
 *       hp::Refinement::choose_p_over_h(dof_handler);
 * 
 * @endcode
 * 
 * For grid adaptive refinement, we ensure a 2:1 mesh balance by limiting
 * the difference of refinement levels of neighboring cells to one by
 * calling Triangulation::prepare_coarsening_and_refinement(). We would
 * like to achieve something similar for the p-levels of neighboring
 * cells: levels of future finite elements are not allowed to differ by
 * more than a specified difference. With its default parameters, a call
 * of hp::Refinement::limit_p_level_difference() ensures that their level
 * difference is limited to one. This will not necessarily decrease the
 * number of hanging nodes in the domain, but makes sure that high order
 * polynomials are not constrained to much lower polynomials on faces,
 * e.g., fifth order to second order polynomials.
 * 
 * @code
 *       triangulation.prepare_coarsening_and_refinement();
 *       hp::Refinement::limit_p_level_difference(dof_handler);
 * 
 * @endcode
 * 
 * At the end of this procedure, we then refine the mesh. During this
 * process, children of cells undergoing bisection inherit their mother
 * cell's finite element index. Further, future finite element indices
 * will turn into active ones, so that the new finite elements will be
 * assigned to cells after the next call of DoFHandler::distribute_dofs().
 * 
 * @code
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemcreate_coarse_grid"></a> 
 * <h4>LaplaceProblem::create_coarse_grid</h4>
 * 

 * 
 * The following function is used when creating the initial grid. The grid we
 * would like to create is actually similar to the one from step-14, i.e., the
 * square domain with the square hole in the middle. It can be generated by
 * exactly the same function. However, since its implementation is only a
 * specialization of the 2d case, we will present a different way of creating
 * this domain which is dimension independent.
 *   

 * 
 * We first create a hypercube triangulation with enough cells so that it
 * already holds our desired domain $[-1,1]^d$, subdivided into $4^d$ cells.
 * We then remove those cells in the center of the domain by testing the
 * coordinate values of the vertices on each cell. In the end, we refine the
 * so created grid globally as usual.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::create_coarse_grid()
 *   {
 *     Triangulation<dim> cube;
 *     GridGenerator::subdivided_hyper_cube(cube, 4, -1., 1.);
 * 
 *     std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
 *     for (const auto &cell : cube.active_cell_iterators())
 *       for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
 *         if (cell->vertex(v).square() < .1)
 *           cells_to_remove.insert(cell);
 * 
 *     GridGenerator::create_triangulation_with_removed_cells(cube,
 *                                                            cells_to_remove,
 *                                                            triangulation);
 * 
 *     triangulation.refine_global(3);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * This function implements the logic of the program, as did the respective
 * function in most of the previous programs already, see for example step-6.
 *   

 * 
 * Basically, it contains the adaptive loop: in the first iteration create a
 * coarse grid, and then set up the linear system, assemble it, solve, and
 * postprocess the solution including mesh refinement. Then start over
 * again. In the meantime, also output some information for those staring at
 * the screen trying to figure out what the program does:
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle)
 *       {
 *         std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           create_coarse_grid();
 * 
 *         setup_system();
 * 
 *         std::cout << "   Number of active cells      : "
 *                   << triangulation.n_active_cells() << std::endl
 *                   << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *                   << std::endl
 *                   << "   Number of constraints       : "
 *                   << constraints.n_constraints() << std::endl;
 * 
 *         assemble_system();
 *         solve();
 *         postprocess(cycle);
 *       }
 *   }
 * } // namespace Step27
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * The main function is again verbatim what we had before: wrap creating and
 * running an object of the main class into a <code>try</code> block and catch
 * whatever exceptions are thrown, thereby producing meaningful output if
 * anything should go wrong:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step27;
 * 
 *       LaplaceProblem<2> laplace_problem;
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
examples/step-27/doc/results.dox



<a name="Results"></a><h1>Results</h1>


在这一节中，我们讨论了运行当前教程程序所产生的一些结果。更多的结果，特别是对三维计算的扩展和确定程序的各个部分需要多少计算时间，在 @ref hp_paper "hp-paper "中给出。

当运行时，这是该程序产生的结果。

@code
> make run
[ 66%] Built target step-27
[100%] Run step-27 with Release configuration
Cycle 0:
   Number of active cells      : 768
   Number of degrees of freedom: 3264
   Number of constraints       : 384
Cycle 1:
   Number of active cells      : 807
   Number of degrees of freedom: 4764
   Number of constraints       : 756
Cycle 2:
   Number of active cells      : 927
   Number of degrees of freedom: 8226
   Number of constraints       : 1856
Cycle 3:
   Number of active cells      : 978
   Number of degrees of freedom: 12146
   Number of constraints       : 2944
Cycle 4:
   Number of active cells      : 1104
   Number of degrees of freedom: 16892
   Number of constraints       : 3998
Cycle 5:
   Number of active cells      : 1149
   Number of degrees of freedom: 22078
   Number of constraints       : 5230
@endcode



我们从中了解到的第一件事是，受限自由度的数量是总自由度的20-25%，至少在后来的网格上，当我们有相对高阶的元素时（在三维中，受限自由度的比例可以达到30%）。事实上，这与非 $hp$ 分化的数量级相同。例如，在第6步程序的最后一步，我们有18353个自由度，其中4432个是受约束的。不同的是，在后一个程序中，每个受约束的悬挂节点只对相邻的两个自由度进行约束，而在 $hp$ -案例中，受约束的节点对许多自由度进行约束。还要注意的是，目前的程序在约束列表中还包括受迪里希特边界条件约束的节点。在第0周期中，所有的约束实际上都是因为边界条件。

也许更感兴趣的是看一下图形输出。首先，这是该问题的解决方案。

<img src="https://www.dealii.org/images/steps/developer/step-27-solution.png" alt="解决方案的立面图，显示出在内部（再入）角落附近缺乏规则性。" width="200" height="200">

其次，让我们看看生成网格的顺序。

<div class="threecolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-00.svg" alt="包含无适应性细化的再入角的三角图。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-01.svg" alt="包含有一级细化的再入角的三角图。新的单元格被放置在角落附近。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-02.svg" alt="包含具有两级细化的重入角的三角图。新的单元格被放置在角落附近。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-03.svg" alt="包含三层细化的重入角的三角图。新的单元格被放置在角落附近。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-04.svg" alt="含有四级细化的重入角的三角图。新的单元格被放置在角落附近。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.mesh-05.svg" alt="包含有五级细化的重入角的三角图。新的单元被放置在角附近。" width="200" height="200"> </div> </div> </div>

我们可以清楚地看到，在角部奇点附近的网格是如何被细化的，正如我们所期望的那样。更有趣的是，我们应该好奇地看看这些网格单元的有限元多项式的度数分布，其中最浅的颜色对应于度数2，最深的对应于度数7。

<div class="threecolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-00.svg" alt="初始网格，所有单元格只包含双二次函数。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-01.svg" alt="一次细化后的局部近似度描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-02.svg" alt="两次细化后的局部近似度描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-03.svg" alt="三次细化后的局部近似度描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-04.svg" alt="四次细化后的局部近似度描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.fe_degree-05.svg" alt="五次细化后的局部近似度描述。" width="200" height="200"> </div> </div>

虽然这肯定不是一个完美的安排，但它确实有一定的意义：我们在靠近边界和角落的地方使用低阶元素，那里的规则性很低。另一方面，在以下情况下使用高阶元素：(i) 误差一度相当大，即主要在角部奇点周围的一般区域和右上角解大的地方；(ii) 解是平滑的，即远离边界的地方。

这种多项式程度的安排当然是由我们的平滑度估计器得出的。这里是对解决方案的平滑度的估计，深色表示最不平滑，浅色表示最平滑的区域。

<div class="threecolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-00.svg" alt="初始网格上每个单元的估计规则性。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-01.svg" alt="经过一次细化后每个单元的估计规则性的描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-02.svg" alt="经过两次细化后每个单元估计规则性的描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-03.svg" alt="经过三次细化后每个单元格的估计规则性的描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-04.svg" alt="经过四次细化后每个单元格的估计规则性的描述。" width="200" height="200"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-27.smoothness-05.svg" alt="经过五次细化后每个单元格的估计规则性的描述。" width="200" height="200"> </div> </div>

从这里可以得出的主要结论是，内角的规则性损失是一个高度局部的现象；它似乎只影响到与角本身相邻的单元，所以当我们细化网格时，黑色的颜色不再可见。除了角部，这一连串的图示意味着平滑度估计在某种程度上与网格细化无关，特别是当我们远离边界时。同样明显的是，平滑度估计值与解的实际尺寸无关（见上面的解的图片），这也是应该的。然而，更值得关注的一点是，人们在仔细观察后发现，我们的估计器高估了有悬空节点的单元格的解决方案的平滑度。这反过来又导致了这些区域的多项式度数更高，使有限元在单元上的分配出现偏差。

对于这种效果，我们目前还没有很好的解释。一种理论是，在有悬挂节点的单元上的数值解当然是受限制的，因此不能完全自由地探索函数空间以接近精确解。这种自由度的缺乏可能表现为在这些单元上产生具有抑制振荡的数值解，意味着更高的平滑度。估计器会捕捉到这个信号，估计的平滑度会高估实际值。然而，这个程序的作者目前还没有得到关于发生了什么的明确答案。

当然，更大的问题是如何避免这个问题。可能的做法包括不是在单个细胞上估计平滑度，而是在每个细胞周围的细胞集合体或斑块上估计平滑度。也可能为每个细胞找到简单的校正因子，这取决于它所拥有的受限自由度的数量。无论哪种情况，都有大量的机会来进一步研究寻找好的 $hp$ -精化标准。另一方面，目前方案的主要内容是在deal.II中演示使用 $hp$ -技术，这不受我们使用可能的次优细化标准的影响。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Differenthpdecisionstrategies"></a><h4>Different hp-decision strategies</h4>


本教程只演示了决定 $h$ -和 $p$ -适应的一个特殊策略。事实上，还有许多自动决定适应类型的方法，其中一些已经在deal.II中实现：  <ul>   <li>  <i>Fourier coefficient decay:</i> 这是本教程中目前实现的策略。有关该策略的更多信息，请参见 SmoothnessEstimator::Fourier 命名空间的一般文档。 </li> 

    <li>  <i>Legendre coefficient decay:</i> 这个策略与目前的策略很相似，但使用了Legendre级数扩展，而不是傅里叶级数：这个策略使用Legendre多项式，而不是正弦波作为基础函数。当然，由于我们在每个单元上使用有限维度的多项式来近似解，因此解在Legendre多项式中的扩展也是有限的，因此，当我们谈论这个扩展的 "衰减 "时，我们只能考虑这个扩展的有限多个非零系数，而不是用渐进的方式来思考。   但是，如果我们有足够多的这些系数，我们当然可以把这些系数的衰减看作是精确解的系数衰减的特征（一般来说，精确解不是多项式的，所以会有一个无限的Legendre扩展），考虑我们拥有的系数应该可以揭示一些关于精确解的特性。

  从傅里叶策略过渡到Legendre策略是非常简单的。   你只需要改变序列扩展类和相应的平滑度估计函数，使其成为适当命名空间 FESeries::Legendre 和 SmoothnessEstimator::Legendre. 的一部分，这个策略在step-75中使用。关于这个策略的理论背景，请参考 SmoothnessEstimator::Legendre 命名空间的一般文档，以及 @cite mavriplis1994hp 、 @cite eibner2007hp 和 @cite davydov2017hp  。 </li> 

    <li>  <i>Refinement history:</i> 最后一种策略与其他两种截然不同。在理论上，我们知道在改变函数空间离散化后，误差将如何收敛。用 $h$ -细化，解会以代数方式收敛，正如在步骤7中已经指出的。不过，如果解足够平滑，我们预计解将随着有限元的多项式程度的增加而指数级收敛。我们可以在下面的步骤中把对误差的正确预测与实际误差进行比较，看看我们对适应类型的选择是否合理。

  向这一策略的过渡要复杂一些。为此，我们需要一个纯 $h$ -或 $p$ -细化的初始化步骤，我们需要将预测的误差转移到适应的网格上。 hp::Refinement::predict_error() 函数的大量文档不仅描述了这种方法的理论细节，而且还介绍了如何在你的代码中实现这种策略的蓝图。欲了解更多信息，请参见  @cite melenk2001hp  。

  请注意，用这个特殊的函数，你无法预测时间依赖性问题的下一个时间步骤的误差。因此，这种策略不能不加思索地应用于这种类型的问题。另外，也可以使用下面的方法，这也适用于所有其他的策略：从每个时间步长的粗网格开始，不断细化，直到对结果满意，然后才进入下一个时间步长。 </li>   </ul> 。

试着在本教程中实施这些策略之一，观察结果的微妙变化。你会注意到，所有的策略都能够识别重心角附近的奇点，并在这些区域执行 $h$ -精简，而在体域中更倾向于 $p$ -精简。这些策略的详细比较见 @cite fehling2020  。




<a name="Parallelhpadaptivefiniteelements"></a><h4>Parallel hp-adaptive finite elements</h4>


本教程中介绍的所有功能都已经适用于顺序和并行应用。不费吹灰之力就可以改成 parallel::shared::Triangulation 或 parallel::distributed::Triangulation 类。如果你觉得急于尝试，我们建议先阅读前者的第18步和后者的第40步，以进一步了解该主题的背景信息，然后再回到本教程来尝试你新获得的技能。

我们在第75步中更进一步：在这里，我们将hp-adapative和MatrixFree方法与 parallel::distributed::Triangulation 对象结合在一起。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-27.cc"
*/
