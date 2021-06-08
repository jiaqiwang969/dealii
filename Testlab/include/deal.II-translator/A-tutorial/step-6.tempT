//include/deal.II-translator/A-tutorial/step-6_0.txt
/**
  @page step_6 The step-6 tutorial program 
* 本教程依赖于  step-5  。
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Whatadaptivelyrefinedmesheslooklike"> What adaptively refined meshes look like </a><a href="#Whatadaptivelyrefinedmesheslooklike"> What adaptively refined meshes look like </a>
        <li><a href="#Whyadapativelyrefinedmeshes"> Why adapatively refined meshes? </a><a href="#Whyadapativelyrefinedmeshes"> Why adapatively refined meshes? </a>
        <li><a href="#Howtodealwithhangingnodesintheory"> How to deal with hanging nodes in theory </a><a href="#Howtodealwithhangingnodesintheory"> How to deal with hanging nodes in theory </a>
        <li><a href="#Howtodealwithhangingnodesinpractice"> How to deal with hanging nodes in practice </a><a href="#Howtodealwithhangingnodesinpractice"> How to deal with hanging nodes in practice </a>
        <li><a href="#Howweobtainlocallyrefinedmeshes"> How we obtain locally refined meshes </a><a href="#Howweobtainlocallyrefinedmeshes"> How we obtain locally refined meshes </a>
        <li><a href="#Boundaryconditions"> Boundary conditions </a><a href="#Boundaryconditions"> Boundary conditions </a>
        <li><a href="#Otherthingsthisprogramshows"> Other things this program shows </a><a href="#Otherthingsthisprogramshows"> Other things this program shows </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeStep6codeclasstemplate">The <code>Step6</code> class template</a><a href="#ThecodeStep6codeclasstemplate">The <code>Step6</code> class template</a>
        <li><a href="#Nonconstantcoefficients">Nonconstant coefficients</a><a href="#Nonconstantcoefficients">Nonconstant coefficients</a>
        <li><a href="#ThecodeStep6codeclassimplementation">The <code>Step6</code> class implementation</a><a href="#ThecodeStep6codeclassimplementation">The <code>Step6</code> class implementation</a>
      <ul>
        <li><a href="#Step6Step6">Step6::Step6</a><a href="#Step6Step6">Step6::Step6</a>
        <li><a href="#Step6setup_system">Step6::setup_system</a><a href="#Step6setup_system">Step6::setup_system</a>
        <li><a href="#Step6assemble_system">Step6::assemble_system</a><a href="#Step6assemble_system">Step6::assemble_system</a>
        <li><a href="#Step6solve">Step6::solve</a><a href="#Step6solve">Step6::solve</a>
        <li><a href="#Step6refine_grid">Step6::refine_grid</a><a href="#Step6refine_grid">Step6::refine_grid</a>
        <li><a href="#Step6output_results">Step6::output_results</a><a href="#Step6output_results">Step6::output_results</a>
        <li><a href="#Step6run">Step6::run</a><a href="#Step6run">Step6::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Solversandpreconditioners">Solvers and preconditioners</a><a href="#Solversandpreconditioners">Solvers and preconditioners</a>
        <li><a href="#Abettermesh">A better mesh</a><a href="#Abettermesh">A better mesh</a>
        <li><a href="#Playingwiththeregularityofthesolution">Playing with the regularity of the solution</a><a href="#Playingwiththeregularityofthesolution">Playing with the regularity of the solution</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
*<a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1> 。


*  @dealiiVideoLecture{15,16,17,17.25,17.5,17.75} 
* 这个程序最终是关于deal.II的主要特征之一：使用自适应（局部）细化网格。这个程序仍然基于 step-4 和 step-5 ，正如你将看到的，它实际上不需要太多的代码来实现自适应性。事实上，虽然我们做了大量的解释，但自适应网格可以被添加到一个现有的程序中，只需十几行的额外代码。该程序显示了这些行是什么，以及自适应网格细化（AMR）的另一个重要成分：一个标准，可以用来确定是否有必要细化一个单元，因为它的误差很大，是否可以粗化该单元，因为它的误差特别小，或者我们是否应该保持该单元的原状。我们将在下文中讨论所有这些问题。
*

*<a name="Whatadaptivelyrefinedmesheslooklike"></a><h3> What adaptively refined meshes look like </h3>


* 有许多方法可以适应性地细化网格。整个算法的基本结构总是相同的，包括以下步骤的循环。
*
* 解决当前网格上的PDE。
*
* - 使用一些指示性的标准来估计每个单元的误差。
*
* - 标记那些有大误差的单元进行细化，标记那些有特别小误差的单元进行粗化，其余的不做处理。
*
* - 细化和粗化那些被标记的单元，以获得一个新的网格。
*
* - 在新的网格上重复上述步骤，直到整体误差足够小。
* 由于一些可能被历史遗忘的原因(可能是这些函数曾经用FORTRAN语言实现，这种语言并不关心某个东西是用小写字母还是大写字母拼写，程序员往往习惯性地选择大写字母)，上述循环在关于网格适应性的出版物中经常被引用为SOLVE-ESTIMATE-MARK-REFINE循环(用这种拼法)。
* 然而，除了这个结构之外，还有各种方法来实现这一目标。从根本上说，它们的区别在于如何从上一个网格生成一个网格。
* 如果使用三角形(deal.II不做)，那么有两种基本的可能性。
*
* - 最长边细化。在这个策略中，通过从最长边的中点到对面的顶点引入一条新的边，将标记为细化的三角形切割成两个。当然，来自最长边的中点必须以某种方式通过也*细化该边另一侧的单元格（如果有的话）来平衡。如果有问题的边也是相邻单元的最长边，那么我们就可以直接通过相邻单元运行一条新的边到对面的顶点；否则就需要一个稍微复杂的构造，在相邻单元的至少一条其他边上增加更多的新顶点，然后可能传播到相邻单元的邻居，直到算法终止。这很难用语言描述，而且因为deal.II不使用三角形，不值得在这里花时间。  但如果你很好奇，你可以随时在本介绍顶部显示的链接中观看视频讲座15。
*
* - 红绿色的细化。另一个选择是所谓的 "红绿细化"。  这种策略甚至更难描述（但也在视频讲座中讨论过），其优点是细化不会传播到我们想要细化的单元的近邻之外。然而，它的实施难度大大增加。
* 这些方法还有其他的变化，但重要的一点是，它们总是产生一个网格，其中两个单元的接触线是两个相邻单元的整个边缘。只要稍加努力，这种策略就可以很容易地适用于由四面体构成的三维网格。
* 这两种方法对二维的四边形和三维的六面体都不起作用，或者至少不容易。原因是，从要精化的四边形单元的四边形邻居中创建的过渡元素将是三角形，而我们并不希望如此。因此，在deal.II中选择的适应性方法是使用网格，其中相邻的单元在细化水平上可能相差一个。这就导致在单元的界面上出现了属于一方的节点，但在另一方是不平衡的。这些节点的通用术语是 "悬空节点"，这些网格在一个非常简单的情况下看起来是这样的。
*  @image html hanging_nodes.png "A simple mesh with hanging nodes" 
* 一个更复杂的二维网格会是这样的（在下面的 "结果 "部分会讨论）。
* <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg" alt="第五个自适应细化的Ladutenko网格：单元格沿着内圈聚集。" width="300" height="300">。
* 最后，一个具有这种悬挂节点的三维网格（来自 step-43 ）显示在这里。
* <img src="https://www.dealii.org/images/steps/developer/  step-43  .3d.mesh.png" alt="" width="300" height="300">。
* 第一个和第三个网格当然是基于一个正方形和一个立方体，但正如第二个网格所显示的，这并不是必要的。重要的一点是，我们可以独立于其邻居来细化一个网格（受制于一个单元只能比其邻居多细化一次的约束），但是如果我们这样做，我们就会出现这些&ldquo;悬空节点&rdquo;。
*

*<a name="Whyadapativelyrefinedmeshes"></a><h3> Why adapatively refined meshes? </h3>


* 现在你已经看到了这些自适应细化网格的样子，你应该问<i>why</i>我们为什么要这样做。毕竟，我们从理论上知道，如果我们对网格进行全局细化，误差会下降到零，因为
* @f{align*}{
  \|\nabla(u-u_h)\|_{\Omega} \le C h_\text{max}^p \| \nabla^{p+1} u \|_{\Omega},
@f}
* 其中 $C$ 是独立于 $h$ 和 $u$ 的某个常数， $p$ 是使用中的有限元的多项式程度， $h_\text{max}$ 是最大单元的直径。那么，如果<i>largest</i>单元很重要，为什么我们要在域的某些部分而不是所有的部分制作精细的网格？
* 答案是上面的公式并不是最优的。事实上，一些更多的工作表明，以下是一个更好的估计（你应该与上述估计的平方进行比较）。
* @f{align*}{
  \|\nabla(u-u_h)\|_{\Omega}^2 \le C \sum_K h_K^{2p} \| \nabla^{p+1} u \|^2_K.
@f}
* (因为 $h_K\le h_\text{max}$ ，如果你只是把网格大小从总和中拉出来，这个公式立刻就暗示了前面的公式。)这个公式所暗示的是，没有必要把<i>largest</i>单元做小，而是单元真的只需要做小<i>where $\| \nabla^{p+1} u \|_K$ is large</i>！换句话说。例如，如果我们使用一个线性元素 $p=1$ ，那么即使网格很粗，那些几乎是线性的地方（如 $\nabla^2 u$ 很小）也会得到很好的解决。只有那些二阶导数大的地方才会被大的元素解决得很差，因此这就是我们应该把网格变小的地方。
* 当然，这个<i>a priori estimate</i>在实践中不是很有用，因为我们不知道问题的精确解 $u$ ，因此，我们不能计算 $\nabla^{p+1}u$ 。但是，这也是通常采取的方法，我们可以只根据我们之前计算的离散解 $u_h$ 来计算 $\nabla^{p+1}u$ 的数字近似值。我们将在下面稍微详细地讨论这个问题。这将有助于我们确定哪些单元具有较大的 $p+1$ st导数，然后这些单元将成为细化网格的候选单元。
*

*<a name="Howtodealwithhangingnodesintheory"></a><h3> How to deal with hanging nodes in theory </h3>


* 上面提到的使用三角形网格的方法要花很大力气确保每个顶点都是所有相邻单元的顶点。
*
* 也就是说，没有悬空的节点。这就自动地确保了我们能够以这样的方式定义形状函数，即它们是全局连续的（如果我们使用迄今为止在教程程序中一直使用的常见的 $Q_p$ Lagrange有限元方法，如FE_Q类所代表的）。
* 另一方面，如果我们在有悬挂节点的网格上定义形状函数，我们最终可能得到不连续的形状函数。请看这一点，想想上面的情况，即右上角的单元没有被细化，并考虑一下使用双线性有限元。在这种情况下，与悬挂节点相关的形状函数是以明显的方式定义在与每个悬挂节点相邻的两个小单元上。但我们如何将它们扩展到相邻的大单元上呢？显然，函数对大单元格的扩展不能是双线性的，因为这样它需要沿着大单元格的每条边线性化，这意味着它在整条边上需要为零，因为它在该边上的大单元格的两个顶点需要为零。但是，当从小单元一侧看时，它在悬空节点本身不是零。
*
* - 所以它不是连续的。以下三幅图显示了沿着有关边的三个形状函数，当以通常的方式简单地根据它们相邻的单元格来定义时，这些形状函数被证明是不连续的。
*  <div class="threecolumn" style="width: 80%">
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_1.png "A discontinuous shape function adjacent to a hanging node"
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_2.png "A discontinuous shape function at a hanging node"
    </div>
  </div>
  <div class="parent">
    <div class="img" align="center">
      @image html hanging_nodes_shape_functions_3.png "A discontinuous shape function adjacent to a hanging node"
    </div>
  </div>
</div> 
*

*但我们确实希望有限元解是连续的，这样我们就有一个&ldquo;符合要求的有限元方法&rdquo;，其中离散有限元空间是我们寻求拉普拉斯方程解的 $H^1$ 函数空间的一个适当子集。诀窍是要认识到，虽然上面显示的形状函数是不连续的（因此它们的<i>arbitrary</i>线性组合也是不连续的），但形状函数加起来为 $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$ 的线性组合可以是连续的<i>if the coefficients $U_j$ satisfy certain relationships</i>。换句话说，系数 $U_j$ 不能任意选择，而必须满足某些约束条件，使函数 $u_h$ 实际上是连续的。这些约束条件在概念上比较容易理解，但在软件中的实现是复杂的，需要几千行代码。另一方面，在用户代码中，当处理悬空节点时，你只需要添加大约半打的行。
* 在下面的程序中，我们将展示如何从deal.II中获得约束条件，以及如何在解决线性方程组中使用它们。在讨论下面程序的细节之前，你可能想看一下 @ref constraints 文档模块，它解释了这些约束如何计算，以及deal.II中的哪些类对它们起作用。
*

*<a name="Howtodealwithhangingnodesinpractice"></a><h3> How to deal with hanging nodes in practice </h3> 。


* 悬空节点约束的实践比我们上述的理论要简单得多。在现实中，你只需要在 step-4 这样的程序中添加半打额外的代码，就可以使它与具有悬挂节点的自适应网格一起工作。有趣的是，这与你要解决的方程完全无关。这些约束的代数性质与方程无关，只取决于有限元的选择。因此，处理这些约束的代码完全包含在deal.II库本身中，你不需要担心细节问题。
* 你需要使其工作的步骤基本上是这样的。
*
* 你必须创建一个AffineConstraints对象，（顾名思义）它将存储有限元空间的所有约束。在目前的情况下，这些约束是由于我们希望在有悬空节点的情况下也能保持解空间的连续。(下面我们还将简要地提到，我们还将把边界值放到这个对象中，但这是一个单独的问题)。
*
* - 你必须使用函数 DoFTools::make_hanging_node_constraints() 来填充这个对象，以确保有限元空间的元素的连续性。
*
* - 当你通过使用 AffineConstraints::distribute_local_to_global(). 将矩阵和右手边的局部贡献复制到全局对象中时，你必须使用这个对象。 到目前为止，我们已经自己完成了这些，但现在有了约束，这就是神奇的地方，我们将约束应用到线性系统中。这个函数所做的是确保位于悬空节点的自由度事实上不是真正的自由。相反，通过将它们的行和列设置为零，并在对角线上放置一些东西以确保矩阵保持可反转，它们实际上被从线性系统中消除了。  对于我们在这里解决的拉普拉斯方程来说，这个过程产生的矩阵仍然是对称和正定的，所以我们可以继续使用共轭梯度法来解决。
*
* - 然后你像往常一样求解线性系统，但在这一步结束时，你需要确保位于悬挂节点上的 "自由度 "得到正确的（约束的）值，这样你随后可视化的或以其他方式评估的解决方案实际上是连续的。这可以通过在求解后立即调用 AffineConstraints::distribute() 来实现。
* 这四个步骤实际上是所有必要的
*
* -从用户的角度来看，就是这么简单。事实上，在上面提到的函数调用中，你将通过几千行并不复杂的代码来运行，这一点完全不重要。在用户代码中，实际上只有四个额外的步骤。
*

*<a name="Howweobtainlocallyrefinedmeshes"></a><h3> How we obtain locally refined meshes </h3>


* 下一个问题，既然我们知道如何<i>deal</i>处理有这些悬挂节点的网格，那么我们如何<i>obtain</i>它们。
* 一个简单的方法已经在 step-1 中展示过了：如果你<i>know</i>有必要细化网格，那么你可以手工创建一个。但是在现实中，我们并不知道这一点。我们不知道前面的PDE的解（因为，如果我们知道，我们就不必使用有限元方法），因此我们不知道哪里需要增加局部网格细化，以更好地解决解有强烈变化的地方。但是上面的讨论表明，也许我们可以用一个网格上的离散解 $u_h$ 来估计导数 $\nabla^{p+1} u$ ，然后用它来确定哪些单元太大，哪些已经足够小。我们可以使用局部网格细化法从当前的网格中生成一个新的网格。
*
* - 或者，更常见的是，直到我们耗尽了计算资源或耐心。
* 所以这正是我们要做的。局部细化的网格是用一个<i>error estimator</i>来产生的，它估计了拉普拉斯算子的数值解的能量误差。由于它是由Kelly和他的同事开发的，所以我们经常在库、文档和邮件列表中把它称为&ldquo;Kelly refinementindicator&rdquo; 。实现它的类被称为KellyErrorEstimator，在该类的文档中可以找到大量的信息，这里不需要重复。然而，摘要是，该类计算一个有 @ref GlossActive 个 "活动单元 "的条目的向量，其中每个条目包含对该单元的误差估计。我们不需要用手来做这个。一旦我们获得了误差估计矢量，GridRefinement中的函数将为我们完成这一切。
* 值得注意的是，虽然Kelly误差估计器是为拉普拉斯方程开发的，但它已被证明是一个合适的工具，可以为广泛的方程生成局部细化网格，甚至不限于椭圆问题。尽管它将为其他方程创建非最优网格，但它通常是快速生成适应解的特征的网格的好方法，如大变化区域或不连续性。
*


*<a name="Boundaryconditions"></a><h3> Boundary conditions </h3>


* 事实证明，我们可以把迪里希特边界条件看作是对自由度的另一种约束。这确实是一个特别简单的条件。如果 $j$ 是边界上的一个自由度，其位置为 $\mathbf x_j$ ，那么在 $\partial\Omega$ 上施加边界条件 $u=g$ 就会产生约束 $U_j=g({\mathbf x}_j)$ 。
*AffineConstraints类也可以处理这样的约束，这使得我们很方便地让我们用于悬挂节点约束的同一个对象也处理这些Dirichlet边界条件。这样，我们不需要在装配后应用边界条件（像我们在前面的步骤中做的那样）。
*

* <a name="Otherthingsthisprogramshows"></a> <h3> Other things this program shows </h3>。


*
* 由于用于局部细化网格的概念非常重要，我们在这个例子中没有展示很多其他材料。最重要的例外是，我们展示了如何使用双二次元，而不是之前所有例子中使用的双线性元素。事实上，使用高阶元素只需替换程序中的三行，即初始化本程序主类构造函数中的 <code>fe</code> 成员变量，以及在两处使用适当的正交公式。程序的其他部分没有变化。
* 唯一的新东西是在 <code>main</code> 函数中捕捉异常的方法，以便在程序因某种原因崩溃时输出一些信息。下面将详细讨论这个问题。
*

* <a name="CommProg"></a> <h1> The commented program</h1>。
* <a name="Includefiles"></a> <h3>Include files</h3>。
*

*
* 前面的几个文件已经在前面的例子中涉及到了，因此将不再进一步评论。
*


* @code
 #include <deal.II/base/quadrature_lib.h>
* 
 #include <deal.II/dofs/dof_handler.h>
 #include <deal.II/dofs/dof_tools.h>
* 
 #include <deal.II/fe/fe_values.h>
* 
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/grid_generator.h>
* 
 #include <deal.II/lac/dynamic_sparsity_pattern.h>
 #include <deal.II/lac/full_matrix.h>
 #include <deal.II/lac/precondition.h>
 #include <deal.II/lac/solver_cg.h>
 #include <deal.II/lac/sparse_matrix.h>
 #include <deal.II/lac/vector.h>
* 
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/numerics/vector_tools.h>
* 
 #include <fstream>
* 
 @endcode
*
* 从下面的include文件中我们将导入H1-conforming有限元形状函数的声明。这个系列的有限元被称为  <code>FE_Q</code>  ，在之前的所有例子中已经用于定义通常的双线性或三线性元素，但我们现在将把它用于双二次元元素。
*


* @code
 #include <deal.II/fe/fe_q.h>
 @endcode
*
* 我们不会像前面的例子那样从文件中读取网格，而是使用库的一个函数来生成网格。然而，我们将希望在每一步中写出局部精炼的网格（只是网格，而不是解决方案），因此我们需要以下的include文件来代替 <code>grid_in.h</code>  。
*


* @code
 #include <deal.II/grid/grid_out.h>
* 
 
 @endcode
*
* 当使用局部细化网格时，我们会得到所谓的<code>悬空节点</code>。然而，标准的有限元方法假定离散解空间是连续的，所以我们需要确保悬挂节点上的自由度符合一些约束条件，从而使全局解是连续的。我们也要在这个对象中存储边界条件。下面的文件包含一个用来处理这些约束条件的类。
*


* @code
 #include <deal.II/lac/affine_constraints.h>
* 
 @endcode
*
* 为了在本地细化我们的网格，我们需要一个来自库的函数，根据我们计算出的误差指标来决定哪些单元需要细化或粗化的标志。这个函数被定义在这里。
*


* @code
 #include <deal.II/grid/grid_refinement.h>
* 
 @endcode
*
* 最后，我们需要一个简单的方法来实际计算基于某种误差估计的细化指标。虽然在一般情况下，适应性是非常具体的问题，但以下文件中的误差指标通常会对一类广泛的问题产生相当好的适应网格。
*


* @code
 #include <deal.II/numerics/error_estimator.h>
* 
 @endcode
*
* 最后，这和以前的程序一样。
*


* @code
 using namespace dealii;
* 
 
 @endcode
*
* <a name="ThecodeStep6codeclasstemplate"></a> <h3>The <code>Step6</code> class template</h3>。


*
* 主类又是几乎没有变化的。然而，我们增加了两项内容：我们增加了 <code>refine_grid</code> 函数，它被用来自适应地细化网格（而不是之前例子中的全局细化），还有一个变量，它将保存约束。
*


* @code
 template <int dim>
 class Step6
 {
 public:
   Step6();
* 
   void run();
* 
 private:
   void setup_system();
   void assemble_system();
   void solve();
   void refine_grid();
   void output_results(const unsigned int cycle) const;
* 
   Triangulation<dim> triangulation;
* 
   FE_Q<dim>       fe;
   DoFHandler<dim> dof_handler;
* 
 
 @endcode
*
* 这是主类中的新变量。我们需要一个对象，它持有一个约束条件的列表，以保持悬挂节点和边界条件。
*


* @code
   AffineConstraints<double> constraints;
* 
   SparseMatrix<double> system_matrix;
   SparsityPattern      sparsity_pattern;
* 
   Vector<double> solution;
   Vector<double> system_rhs;
 };
* 
 
 @endcode
*
* <a name="Nonconstantcoefficients"></a> <h3>Nonconstant coefficients</h3>。


*
* 非恒定系数的实现是逐字复制自  step-5  。
*


* @code
 template <int dim>
 double coefficient(const Point<dim> &p)
 {
   if (p.square() < 0.5 0.5)
     return 20;
   else
     return 1;
 }
* 
 
* 
 @endcode
*
* <a name="ThecodeStep6codeclassimplementation"></a> <h3>The <code>Step6</code> class implementation</h3>。



* <a name="Step6Step6"></a> <h4>Step6::Step6</h4>。


*
* 这个类的构造函数与之前的基本相同，但是这一次我们要使用二次元。为此，我们只需用所需的多项式度数（这里是 <code>2</code> ）来替换构造函数的参数（在以前的例子中是 <code>1</code> ）。
*


* @code
 template <int dim>
 Step6<dim>::Step6()
   : fe(2)
   , dof_handler(triangulation)
 {}
* 
 
* 
 @endcode
*
* <a name="Step6setup_system"></a> <h4>Step6::setup_system</h4>。


*
* 下一个函数设置了所有描述线性有限元问题的变量，如DoFHandler、矩阵和向量。与我们在 step-5 中所做的不同的是，我们现在还必须照顾到悬挂节点的约束。这些约束几乎完全由库来处理，也就是说，你只需要知道它们的存在以及如何获得它们，但你不必知道它们是如何形成的，也不必知道对它们到底做了什么。
*

*
* 在函数的开头，你会发现所有与 step-5 中相同的东西：设置自由度（这次我们有二次元，但从用户代码的角度看与线性的没有区别
*
* - 或任何其他程度，对于这个问题
*
* 情况下），生成稀疏模式，并初始化解和右手边的向量。请注意，现在每行的稀疏模式将有更多的条目，因为现在每个单元有9个自由度（而不是只有4个），它们可以相互耦合。
*


* @code
 template <int dim>
 void Step6<dim>::setup_system()
 {
   dof_handler.distribute_dofs(fe);
* 
   solution.reinit(dof_handler.n_dofs());
   system_rhs.reinit(dof_handler.n_dofs());
* 
 @endcode
*
* 我们现在可以用悬挂节点的约束来填充AffineConstraints对象。由于我们将在一个循环中调用这个函数，我们首先要清除上一个系统中的当前约束集，然后计算新的约束。
*


* @code
   constraints.clear();
   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
* 
 
 @endcode
*
* 现在我们准备用指标0（整个边界）来插值边界值，并将产生的约束存储在我们的 <code>constraints</code> 对象中。请注意，我们并不像在前面的步骤中那样，在装配后应用边界条件：相反，我们将所有的约束条件放在AffineConstraints对象中的我们的函数空间。我们可以以任何顺序向AffineConstraints对象添加约束：如果两个约束发生冲突，那么约束矩阵要么中止，要么通过Assert宏抛出一个异常。
*


* @code
   VectorTools::interpolate_boundary_values(dof_handler,
                                            0,
                                            Functions::ZeroFunction<dim>(),
                                            constraints);
* 
 @endcode
*
* 在所有约束条件被添加后，需要对它们进行排序和重新排列，以便更有效地执行一些操作。这种后处理是用 <code>close()</code> 函数完成的，之后就不能再添加任何约束了。
*


* @code
   constraints.close();
* 
 @endcode
*
* 现在我们首先建立我们的压缩稀疏度模式，就像我们在前面的例子中做的那样。尽管如此，我们并没有立即将其复制到最终的稀疏度模式中。  请注意，我们调用了make_sparsity_pattern的一个变体，它把AffineConstraints对象作为第三个参数。我们通过将参数 <code>keep_constrained_dofs</code> 设置为false（换句话说，我们永远不会写入对应于受限自由度的矩阵条目），让该例程知道我们永远不会写入 <code>constraints</code> 所给的位置。如果我们要在组装后压缩约束，我们就必须通过 <code>true</code> 来代替，因为这样我们就会先写进这些位置，然后在压缩过程中再将它们设置为零。
*


* @code
   DynamicSparsityPattern dsp(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler,
                                   dsp,
                                   constraints,
                                    /*keep_constrained_dofs = */  false);
* 
 @endcode
*
* 现在，矩阵的所有非零条目都是已知的（即那些来自定期组装矩阵的条目和那些通过消除约束引入的条目）。我们可以将我们的中间对象复制到稀疏模式中。
*


* @code
   sparsity_pattern.copy_from(dsp);
* 
 @endcode
*
* 我们现在可以，最后，初始化稀疏矩阵。
*


* @code
   system_matrix.reinit(sparsity_pattern);
 }
* 
 
 @endcode
*
* <a name="Step6assemble_system"></a> <h4>Step6::assemble_system</h4>。


*
* 接下来，我们要对矩阵进行组装。然而，为了将每个单元上的本地矩阵和向量复制到全局系统中，我们不再使用手写的循环。相反，我们使用 AffineConstraints::distribute_local_to_global() ，在内部执行这个循环，同时对对应于受限自由度的行和列进行高斯消除。
*

*
* 形成本地贡献的其他代码保持不变。然而，值得注意的是，在引擎盖下，有几处地方与以前不同。首先，变量 <code>dofs_per_cell</code> 和返回值 <code>quadrature_formula.size()</code> 现在各为9，以前是4。引入这样的变量作为缩写是一个很好的策略，可以使代码与不同的元素一起工作，而不需要改变太多的代码。其次， <code>fe_values</code> 对象当然也需要做其他事情，因为现在的形状函数是二次的，而不是线性的，在每个坐标变量中。不过，这也是完全由库来处理的事情。
*


* @code
 template <int dim>
 void Step6<dim>::assemble_system()
 {
   const QGauss<dim> quadrature_formula(fe.degree + 1);
* 
   FEValues<dim> fe_values(fe,
                           quadrature_formula,
                           update_values | update_gradients |
                             update_quadrature_points | update_JxW_values);
* 
   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
* 
   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double>     cell_rhs(dofs_per_cell);
* 
   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
* 
   for (const auto &cell : dof_handler.active_cell_iterators())
     {
       cell_matrix = 0;
       cell_rhs    = 0;
* 
       fe_values.reinit(cell);
* 
       for (const unsigned int q_index : fe_values.quadrature_point_indices())
         {
           const double current_coefficient =
             coefficient(fe_values.quadrature_point(q_index));
           for (const unsigned int i : fe_values.dof_indices())
             {
               for (const unsigned int j : fe_values.dof_indices())
                 cell_matrix(i, j) +=
                   (current_coefficient              // a(x_q)
                    fe_values.shape_grad(i, q_index) // grad phi_i(x_q)
                    fe_values.shape_grad(j, q_index) // grad phi_j(x_q)
                    fe_values.JxW(q_index));           // dx
* 
               cell_rhs(i) += (1.0                               // f(x)
                               fe_values.shape_value(i, q_index) // phi_i(x_q)
                               fe_values.JxW(q_index));            // dx
             }
         }
* 
 @endcode
*
* 最后，将 @p cell_matrix 和 @p cell_rhs 中的贡献转入全局对象。
*


* @code
       cell->get_dof_indices(local_dof_indices);
       constraints.distribute_local_to_global(
         cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
     }
 @endcode

* 现在我们已经完成了线性系统的组装。约束矩阵照顾到了边界条件的应用，也消除了悬挂的节点约束。受约束的节点仍然在线性系统中（在矩阵的对角线上有一个非零条目，选择的方式是使矩阵具有良好的条件，并且这一行的所有其他条目都被设置为零），但是计算出来的值是无效的（也就是说， <code>system_rhs</code> 中的相应条目目前是没有意义的）。我们在 <code>solve</code> 函数的最后为这些节点计算出正确的值。
*


* @code
 }
* 
 
 @endcode
*
* <a name="Step6solve"></a> <h4>Step6::solve</h4>。


*
* 我们继续逐步改进。解决线性系统的函数再次使用了SSOR预处理程序，而且除了我们必须加入悬空节点约束外，再次没有任何变化。如上所述，通过对矩阵的行和列进行特殊处理，从AffineConstraints对象中删除了对应于悬挂节点约束和边界值的自由度。这样一来，这些自由度的值在求解线性系统后就有了错误的、但定义明确的值。然后我们要做的就是利用约束条件给它们分配它们应该有的值。这个过程被称为 <code>distributing</code> 约束，从无约束的节点的值中计算出约束节点的值，只需要一个额外的函数调用，你可以在这个函数的末尾找到。
*


*


* @code
 template <int dim>
 void Step6<dim>::solve()
 {
   SolverControl            solver_control(1000, 1e-12);
   SolverCG<Vector<double>> solver(solver_control);
* 
   PreconditionSSOR<SparseMatrix<double>> preconditioner;
   preconditioner.initialize(system_matrix, 1.2);
* 
   solver.solve(system_matrix, solution, system_rhs, preconditioner);
* 
   constraints.distribute(solution);
 }
* 
 
 @endcode
*
* <a name="Step6refine_grid"></a> <h4>Step6::refine_grid</h4>。


*
* 我们使用一个复杂的误差估计方案来细化网格，而不是全局细化。我们将使用KellyErrorEstimator类，该类实现了拉普拉斯方程的误差估计器；原则上它可以处理可变系数，但我们不会使用这些高级功能，而是使用其最简单的形式，因为我们对定量结果不感兴趣，只对生成局部细化网格的快速方法感兴趣。
*

*
* 尽管Kelly等人推导的误差估计器最初是为拉普拉斯方程开发的，但我们发现它也很适合为一类广泛的问题快速生成局部细化网格。这个误差估计器使用了解梯度在单元面上的跳跃（这是一个测量二阶导数的方法），并将其按单元的大小进行缩放。因此，它是对每个单元的解的局部平滑性的测量，因此可以理解，它对双曲运输问题或波浪方程也能产生合理的网格，尽管这些网格与专门为该问题定制的方法相比肯定是次优的。因此，这个误差估计器可以理解为测试自适应程序的一种快速方法。
*

*
* 估算器的工作方式是将描述自由度的 <code>DoFHandler</code> 对象和每个自由度的数值向量作为输入，并为三角形的每个活动单元计算一个指标值（即每个活动单元一个数值）。为此，它需要两个额外的信息：一个面部正交公式，即 <code>dim-1</code> 维物体上的正交公式。我们再次使用3点高斯法则，这个选择与本程序中的双二次方有限元形状函数是一致和合适的。当然，什么是合适的正交规则取决于对误差估计器评估解场的方式的了解。如上所述，梯度的跳跃在每个面上都是集成的，对于本例中使用的二次元元素来说，这将是每个面上的二次元函数。然而，事实上，它是梯度跳跃的平方，正如该类文件中所解释的那样，这是一个二次函数，对于它来说，3点高斯公式就足够了，因为它可以精确地整合5阶以下的多项式）。)
*

*
* 其次，该函数需要一个边界指标列表，用于那些我们施加了 $\partial_n u(\mathbf x) = h(\mathbf x)$ 类型的诺伊曼值的边界，以及每个此类边界的一个函数 $h(\mathbf
 x)$ 。这些信息通过一个从边界指标到描述诺伊曼边界值的函数对象的映射来表示。在本例程序中，我们不使用诺伊曼边界值，所以这个映射是空的，实际上是在函数调用期望得到相应函数参数的地方使用映射的默认构造器构造的。
*

*
* 输出是所有活动单元的值的向量。虽然非常精确地计算一个解的自由度的<b>value</b>可能是有意义的，但通常没有必要特别精确地计算一个单元上对应于解的<b>error
 indicator</b>。因此，我们通常使用一个浮点数的矢量而不是一个双数的矢量来表示误差指标。
*


* @code
 template <int dim>
 void Step6<dim>::refine_grid()
 {
   Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
* 
   KellyErrorEstimator<dim>::estimate(dof_handler,
                                      QGauss<dim
* 
- 1>(fe.degree + 1),
                                      {},
                                      solution,
                                      estimated_error_per_cell);
* 
 @endcode
*
* 上述函数为 <code>estimated_error_per_cell</code> 数组中的每个单元格返回一个错误指标值。现在的细化工作如下：细化那些误差值最高的30%的单元格，粗化那些误差值最低的3%的单元格。   
* 人们可以很容易地验证，如果第二个数字为零，这大约会导致在两个空间维度上的每一步中的细胞增加一倍，因为对于每一个30%的细胞，四个新的将被替换，而其余70%的细胞保持不动。在实践中，通常会产生一些更多的单元，因为不允许一个单元被精炼两次而相邻的单元没有被精炼；在这种情况下，相邻的单元也会被精炼。   
* 在许多应用中，被粗化的单元格数量将被设置为大于百分之三的数值。一个非零的值是有用的，特别是当由于某种原因，初始（粗）网格已经相当精细了。在这种情况下，可能有必要在某些区域进行细化，而在另一些区域进行粗化是有用的。在我们这里，初始网格是非常粗的，所以粗化只需要在一些可能发生过度细化的区域。因此，一个小的、非零的值在这里是合适的。   
* 下面的函数现在采用这些细化指标，并使用上述方法对三角形的一些单元进行细化或粗化标记。它来自一个类，该类实现了几种不同的算法，根据单元的误差指标来细化三角形。
*


* @code
   GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                   estimated_error_per_cell,
                                                   0.3,
                                                   0.03);
* 
 @endcode
*
* 在前一个函数退出后，一些单元被标记为细化，另一些单元被标记为粗化。然而，细化或粗化本身并没有被执行，因为在有些情况下，进一步修改这些标志是有用的。在这里，我们不想做任何这样的事情，所以我们可以告诉三角计算执行单元格被标记的动作。
*


* @code
   triangulation.execute_coarsening_and_refinement();
 }
* 
 
 @endcode
*
* <a name="Step6output_results"></a> <h4>Step6::output_results</h4>。


*
* 在每个网格的计算结束后，就在我们继续下一个网格细化周期之前，我们要输出这个周期的结果。
*

*
* 我们已经在 step-1 中看到了如何对网格本身实现这一点。在这里，我们改变了一些东西。   <ol>   <li>  我们使用两种不同的格式。gnuplot和VTU。 </li>   <li>  我们将周期号嵌入到输出文件名中。 </li>   <li>  对于gnuplot输出，我们设置了一个 GridOutFlags::Gnuplot 对象，以提供一些额外的可视化参数，使边缘看起来是弯曲的。这在  step-10  中有进一步的详细解释。 </li>   </ol> 。
*


* @code
 template <int dim>
 void Step6<dim>::output_results(const unsigned int cycle) const
 {
   {
     GridOut               grid_out;
     std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot");
     GridOutFlags::Gnuplot gnuplot_flags(false, 5);
     grid_out.set_flags(gnuplot_flags);
     MappingQGeneric<dim> mapping(3);
     grid_out.write_gnuplot(triangulation, output, &mapping);
   }
* 
   {
     DataOut<dim> data_out;
     data_out.attach_dof_handler(dof_handler);
     data_out.add_data_vector(solution, "solution");
     data_out.build_patches();
* 
     std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
     data_out.write_vtu(output);
   }
 }
* 
 
 @endcode
*
* <a name="Step6run"></a> <h4>Step6::run</h4>。


*
*  <code>main()</code> 之前的最后一个函数仍然是该类的主要驱动， <code>run()</code>  。它与  step-5  的函数类似，只是我们在程序中再次生成一个文件，而不是从磁盘中读取，我们自适应地而不是全局地细化网格，并且我们在本函数中输出最终网格上的解。
*

*
* 该函数主循环的第一个块是处理网格的生成。如果这是该程序的第一个循环，而不是像前面的例子那样从磁盘上的文件中读取网格，我们现在再次使用一个库函数来创建它。域还是一个圆，中心在原点，半径为1（这是函数的两个隐藏参数，有默认值）。
*

*
* 通过观察粗略的网格，你会发现它的质量比我们在前面的例子中从文件中读出的网格要差：单元格的形成不太平均。然而，使用库函数，这个程序可以在任何空间维度上工作，这在以前是不可能的。
*


* 万一我们发现这不是第一个周期，我们要细化网格。与上一个例子程序中采用的全局细化不同，我们现在使用上述的适应性程序。
*

*
* 循环的其余部分看起来和以前一样。
*


* @code
 template <int dim>
 void Step6<dim>::run()
 {
   for (unsigned int cycle = 0; cycle < 8; ++cycle)
     {
       std::cout << "Cycle " << cycle << ':' << std::endl;
* 
       if (cycle == 0)
         {
           GridGenerator::hyper_ball(triangulation);
           triangulation.refine_global(1);
         }
       else
         refine_grid();
* 
 
       std::cout << "   Number of active cells:       "
                 << triangulation.n_active_cells() << std::endl;
* 
       setup_system();
* 
       std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                 << std::endl;
* 
       assemble_system();
       solve();
       output_results(cycle);
     }
 }
* 
 
 @endcode
*
* <a name="Thecodemaincodefunction"></a> <h3>The <code>main</code> function</h3>。


*
* 主函数的功能与前面的例子相比没有改变，但我们采取了额外的谨慎措施。有时，会出现一些问题（比如写输出文件时磁盘空间不足，试图分配向量或矩阵时内存不足，或者由于某种原因我们无法从文件中读出或写入文件），在这些情况下，库会抛出异常。由于这些都是运行时的问题，而不是可以一劳永逸的编程错误，所以这种异常在优化模式下是不会被关闭的，这与我们用来测试编程错误的 <code>Assert</code> 宏相反。如果没有被捕获，这些异常会传播到 <code>main</code> 函数的调用树上，如果它们在那里也没有被捕获，程序就会被中止。在很多情况下，比如内存或磁盘空间不足，我们什么也做不了，但我们至少可以打印一些文字，试图解释程序失败的原因。下面显示了一种方法。以这种方式编写任何较大的程序当然是有用的，你可以通过或多或少地复制这个函数来做到这一点，除了 <code>try</code> 块，它实际上编码了本应用程序所特有的功能。
*


* @code
 int main()
 {
 @endcode
*
* 这个函数布局的总体思路如下：让我们试着像以前那样运行程序......
*


* @code
   try
     {
       Step6<2> laplace_problem_2d;
       laplace_problem_2d.run();
     }
 @endcode
*
* ...如果这应该是失败的，尽量收集尽可能多的信息。具体来说，如果被抛出的异常是一个从C++标准类 <code>exception</code>, then we can use the <code>what</code> 派生出来的类的对象，可以通过成员函数得到一个字符串，这个字符串描述了异常被抛出的原因。   
* deal.II的异常类都是从标准类派生出来的，特别是 <code>exc.what()</code> 函数将返回与使用 <code>Assert</code> 宏抛出的异常所产生的字符串大致相同。在前面的例子中，你已经看到了这种异常的输出，然后你知道它包含了异常发生的文件和行号，以及其他一些信息。这也是下面的语句会打印的内容。   
* 除此以外，除了用错误代码退出程序（这就是 <code>return 1;</code> 的作用），我们能做的并不多。
*


* @code
   catch (std::exception &exc)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Exception on processing: " << std::endl
                 << exc.what() << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
* 
       return 1;
     }
 @endcode
*
* 如果在某处抛出的异常不是从标准 <code>exception</code> 类派生出来的对象，那么我们根本就不能做任何事情。我们就简单地打印一个错误信息并退出。
*


* @code
   catch (...)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Unknown exception!" << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       return 1;
     }
* 
 @endcode
*
* 如果我们走到这一步，就没有任何异常传播到主函数上（可能有异常，但它们在程序或库的某个地方被捕获）。因此，程序按预期执行，我们可以无误返回。
*


* @code
   return 0;
 }
 @endcode
* <a name="Results"></a><h1>Results</h1> 。



* 程序的输出看起来如下。
* @code
Cycle 0:
   Number of active cells:       20
   Number of degrees of freedom: 89
Cycle 1:
   Number of active cells:       44
   Number of degrees of freedom: 209
Cycle 2:
   Number of active cells:       92
   Number of degrees of freedom: 449
Cycle 3:
   Number of active cells:       200
   Number of degrees of freedom: 921
Cycle 4:
   Number of active cells:       440
   Number of degrees of freedom: 2017
Cycle 5:
   Number of active cells:       956
   Number of degrees of freedom: 4425
Cycle 6:
   Number of active cells:       1916
   Number of degrees of freedom: 8993
Cycle 7:
   Number of active cells:       3860
   Number of degrees of freedom: 18353
@endcode
*
*

*
* 正如预期的那样，在每个周期中，细胞的数量大约增加了一倍。度数是单元数的四倍多；在无限网格的两个空间维度上，我们可以期待一个正好是四的系数（因为自由度之间的间距是单元宽度的一半：每个边缘有一个额外的自由度，每个单元的中间有一个），但是由于网格的有限尺寸和由悬挂节点和局部细化引入的额外自由度，它比这个系数大。
*

*
* 该程序在每个循环中输出解决方案和网格。解决方案看起来如下。
*  <img src="https://www.dealii.org/images/steps/developer/step-6.solution.9.2.png" alt=""> 
* 有趣的是，我们可以看到程序是如何得出最终的网格的。
*  <div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0.svg"
         alt="Initial grid: the five-cell circle grid with one global refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1.svg"
         alt="First grid: the five-cell circle grid with two global refinements."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2.svg"
         alt="Second grid: the five-cell circle grid with one adaptive refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3.svg"
         alt="Third grid: the five-cell circle grid with two adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4.svg"
         alt="Fourth grid: the five-cell circle grid with three adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5.svg"
         alt="Fifth grid: the five-cell circle grid with four adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_6.svg"
         alt="Sixth grid: the five-cell circle grid with five adaptive
         refinements, showing clustering around the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_7.svg"
         alt="Last grid: the five-cell circle grid with six adaptive
         refinements, showing that most cells are clustered around the inner circle."
         width="300" height="300">
  </div>
</div> 


* 可以清楚地看到，在解有扭结的区域，即距中心0.5的径向距离的圆，被精炼得最多。此外，解非常光滑且几乎平坦的中心区域几乎完全没有被细化，但这是由于我们没有考虑到那里的系数很大的事实。外面的区域被任意地细化了，因为那里的二阶导数是恒定的，因此细化主要是基于单元的大小和它们与最优方格的偏差。
*


* <a name="extensions"></a><a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3> 。


*<a name="Solversandpreconditioners"></a><h4>Solvers and preconditioners</h4> 。


*
* 如果一个人要解决相当大的问题（比我们这里的问题大得多），有一件事总是值得一试的，那就是尝试不同的求解器或预处理器。在目前的情况下，线性系统是对称的和正定的，这使得CG算法几乎成了解题的典型选择。然而，我们在 <code>solve()</code> 函数中使用的SSOR预处理器是有待商榷的。
* 在deal.II中，改变预处理程序是比较简单的。例如，通过改变现有的几行代码
* @code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
@endcode
* 变成
* @code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);
@endcode
* 我们可以尝试SSOR的不同松弛参数。通过使用
* @code
  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode
* 我们可以使用Jacobi作为预处理程序。而通过使用
* @code
  SparseILU<double> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode
* 我们可以使用一个简单的不完整的LU分解，没有任何对角线的阈值或加强（要使用这个预处理程序，你还必须把头文件 <code>deal.II/lac/sparse_ilu.h</code> 加入到文件顶部的include列表中）。
* 使用这些不同的预处理程序，我们可以比较所需的CG迭代次数（可通过 <code>solver_control.last_step()</code> 调用，见 step-4 ）以及所需的CPU时间（使用Timer类，例如在 step-28 中讨论），得到以下结果（左边：迭代次数；右边：CPU时间）。
*  <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table> 
* 我们可以看到，所有的预处理程序在这个简单问题上的表现都差不多，迭代次数以 ${\cal
O}(N^{1/2})$ 的方式增长，由于每次迭代需要大约 ${\cal
O}(N)$ 次操作，总的CPU时间以 ${\cal
O}(N^{3/2})$ 的方式增长（对于几个最小的网格，CPU时间非常小，以至于没有记录）。请注意，即使是简单的方法，Jacobi对于这个问题也是最快的。
* 当有限元不是本程序构造函数中设定的二次元，而是二次元时，情况会有一些变化。如果作出这种改变，结果如下。
*  <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table> 
* 换句话说，虽然迭代次数和CPU时间的增加和以前一样，但雅可比现在是需要最多迭代的方法；不过，由于它必须执行的操作很简单，所以它仍然是最快的方法。这并不是说Jacobi实际上是一个好的预处理程序。
*
* 对于相当大的问题，它肯定不是，其他方法会好得多。
*
* 但实际上只是因为它的实现非常简单，所以它可以补偿更多的迭代次数。
* 从这里得到的信息并不是说预设条件器的简单性总是最好的。虽然这对当前的问题可能是正确的，但一旦我们转向更复杂的问题（弹性或斯托克斯，例如 step-8 或 step-22 ），就绝对不是这样了。其次，所有这些预处理程序仍然导致迭代次数随着自由度 $N$ 的增加而增加，例如 ${\cal O}(N^\alpha)$ ；这反过来又导致总工作量的增加，因为每次迭代都需要 ${\cal O}(N)$ 工作。这种行为是不可取的：我们真的希望用 $N$ 个未知数解决线性系统，总工作量为 ${\cal O}(N)$ 个；有一类预处理程序可以实现这一点，即几何（ step-16 ,  step-37 ,  step-39 ）或代数多网格（ step-31 ,  step-40  ，以及其他几个）预处理程序。然而，它们要比上述的预处理程序复杂得多。
* 最后，要带回家的最后一个信息是，当上面显示的数据产生时（2018年），有10万个未知数的线性系统在台式机上很容易在大约一秒钟内解决，使相对简单的2D问题的解决甚至达到非常高的精度，即使在过去也不是一个大任务。当时，三维问题的情况完全不同，但即使是这样，在过去的时间里也发生了很大的变化。
*
* - 尽管在三维中解决高精确度的问题仍然是一个挑战。
*

*<a name="Abettermesh"></a><h4>A better mesh</h4>


* 如果你看一下上面的网格，你会发现尽管领域是一个圆盘，而且系数的跳动是沿着一个圆的，但是组成网格的单元并没有很好地跟踪这个几何形状。原因在 step-1 中已经提示过了，在没有其他信息的情况下，三角剖分类只看到一堆粗大的网格单元，但当然不知道它们在一起时可能代表什么样的几何。出于这个原因，我们需要告诉Triangulation当一个单元被细化时应该做什么：边缘中点和单元中点的新顶点应该位于哪里，以便子单元比父单元更好地代表所需的几何形状。
* 为了直观地了解三角测量对几何体的实际了解，仅仅输出顶点的位置和为每条边画一条直线是不够的；相反，我们必须将内部线和边界线都输出为多条线段，使它们看起来是弯曲的。我们可以通过对 <code>output_results</code> 的gnuplot部分做一个改动来实现这一点。
* @code
{
  GridOut       grid_out;
  std::ofstream output("grid-" + std::to_string(cycle) + ".gnuplot");
  GridOutFlags::Gnuplot gnuplot_flags(false, 5,  /*curved_interior_cells*/ true);
  grid_out.set_flags(gnuplot_flags);
  MappingQGeneric<dim> mapping(3);
  grid_out.write_gnuplot(triangulation, output, &mapping);
}
@endcode
*
* 在上面的代码中，我们已经对位于边界的面做了这个处理：由于我们使用了 GridGenerator::hyper_ball, ，它将SphericalManifold附加到域的边界上，所以这就自动发生了。为了使网格<i>interior</i>也能追踪到一个圆形域，我们需要更努力一些。首先，回顾一下我们的粗略网格是由一个中心的正方形单元和周围的四个单元组成的。现在首先考虑一下，如果我们将球形网格对象不仅连接到四个外部面，而且连接到周边的四个单元以及它们的所有面，会发生什么。我们可以通过添加下面的片段来实现这一点（测试一个单元格的中心是否大于单元格直径的一个小倍数，比如说十分之一，只对网格中心的正方形无效）。
* @code
GridGenerator::hyper_ball(triangulation);
* 


* 
const Point<dim> mesh_center;
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance (cell->center()) > cell->diameter()/10)
    cell->set_all_manifold_ids(0);
* 
triangulation.refine_global(1);
@endcode
*
* 经过几个全局细化步骤后，这将导致以下类型的网格。
*

*  <div class="onecolumn" style="width: 80%">
    <div>
      <img src="https://www.dealii.org/images/steps/developer/step_6_bad_grid_4.svg"
           alt="Grid where some central cells are nearly triangular."
           width="300" height="300">
    </div>
  </div>
* 
This is not a good mesh: the central cell has been refined in such a way that
the children located in the four corners of the original central cell
<i>degenerate</i>: they all tend towards triangles as mesh refinement
continues. This means that the Jacobian matrix of the transformation from
reference cell to actual cell degenerates for these cells, and because
all error estimates for finite element solutions contain the norm of the
inverse of the Jacobian matrix, you will get very large errors on these
cells and, in the limit as mesh refinement, a loss of convergence order because
the cells in these corners become worse and worse under mesh refinement.
* 
So we need something smarter. To this end, consider the following solution
originally developed by Konstantin Ladutenko. We will use the following code:
* 
@code
GridGenerator::hyper_ball(triangulation);
* 
const Point<dim> mesh_center;
const double core_radius  = 1.0/5.0,
             inner_radius = 1.0/3.0;
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) < 1e-5)
    {
      for (const auto v : cell->vertex_indices())
        cell->vertex(v)= core_radius/mesh_center.distance(cell->vertex(v));
    }
* 


* 
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) >= 1e-5)
    cell->set_refine_flag();
triangulation.execute_coarsening_and_refinement();
* 


* 
* 


* 
* 


* 
* 


* 
* 
for (const auto &cell : triangulation.active_cell_iterators())
  for (const auto v : cell->vertex_indices())
    {
      const double dist = mesh_center.distance(cell->vertex(v));
      if (dist > core_radius*1.0001 && dist < 0.9999)
        cell->vertex(v)= inner_radius/dist;
    }
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
for (const auto &cell : triangulation.active_cell_iterators())
  {
    bool is_in_inner_circle = false;
    for (const auto v : cell->vertex_indices())
      if (mesh_center.distance(cell->vertex(v)) < inner_radius)
        {
          is_in_inner_circle = true;
          break;
        }
* 
    if (is_in_inner_circle == false)
    // The Triangulation already has a SphericalManifold with
    // manifold id 0 (see the documentation of
    // GridGenerator::hyper_ball) so we just attach it to the outer
    // ring here:
      cell->set_all_manifold_ids(0);
  }
@endcode
* 
* 
This code then generates the following, much better sequence of meshes:
* 
<div class="twocolumn" style="width: 80%">
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0_ladutenko.svg"
         alt="Initial grid: the Ladutenko grid with one global refinement."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1_ladutenko.svg"
         alt="First adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2_ladutenko.svg"
         alt="Second adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3_ladutenko.svg"
         alt="Third adaptively refined Ladutenko grid."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4_ladutenko.svg"
         alt="Fourth adaptively refined Ladutenko grid. The cells are clustered
         along the inner circle."
         width="300" height="300">
  </div>
  <div>
    <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg"
         alt="Fifth adaptively refined Ladutenko grid: the cells are clustered
         along the inner circle."
         width="300" height="300">
  </div>
</div> 
* 创建好的网格，特别是使它们适合你想要的几何形状，本身就是一个复杂的话题。你可以在 step-49 、 step-53 和 step-54 中找到更多关于这个问题的内容，以及其他涉及这个问题的教程程序。   step-65 展示了另一种不太需要手工操作的方法来实现对这里的问题的网格井拟合。关于弯曲域的信息也可以在 @ref manifold 的文档模块 "流形描述 "中找到。
* 为什么选择一个跟踪内部界面的网格是有意义的？有很多原因，但最重要的原因是我们在双线型中实际集成了什么。从概念上讲，我们想把 $A_{ij}^K=\int_K
a(\mathbf x) \nabla \varphi_i(\mathbf x) \nabla \varphi_j(\mathbf x) ; dx$ 作为单元格 $K$ 对矩阵条目 $A_{ij}$ 的贡献进行积分。我们无法精确计算，只能采用正交法。我们知道，如果积分是平滑的，正交是准确的。这是因为正交法实质上是计算积分的多项式近似值，与正交点的积分重合，然后计算该多项式下的体积，作为原始积分下的体积的近似值。如果积分在一个单元上是平滑的，那么这个多项式插值是准确的，但如果积分在一个单元上是不连续的，那么它通常是相当不准确的。
* 因此，值得注意的是，在排列单元时，系数不连续的界面要与单元的界面对齐。这样，系数在每个单元上都是恒定的，这样积分就会很平滑，其多项式近似和积分的正交近似都会很准确。请注意，这样的排列方式在实际情况中很常见，所以deal.II提供了一些函数（如 @ref GlossMaterialId "material_id"）来帮助管理这样的情况。请参考 step-28 和 step-46 关于如何应用material ids的例子。
* 最后，让我们考虑一个在空间中具有平滑和非均匀分布的系数的情况。我们可以再一次重复上面关于用正交函数表示这种函数的所有讨论。因此，为了准确地模拟它，有几个现成的选择：你可以减少单元格的大小，增加正交公式中使用的多项式的阶数，选择一个更合适的正交公式，或执行这些步骤的组合。关键是，用正交多项式提供系数的空间依赖性的最佳拟合，将导致PDE的更精确的有限元解。
* 最后说明：前面的讨论表明，我们在这里有一个非常具体的方法来说明我们认为的好的网格
*
*-它应该与系数的跳动保持一致。但是我们也可以在一个更普遍的环境中提出这样的问题。给出一些具有光滑解和光滑系数的方程，我们能说一个好的网格是什么样子的吗？这个问题的答案在直觉上比数学上更容易表述。一个好的网格所包含的单元大体上都像正方形（或立方体，在三维空间）。一个不好的网格会包含一些在某些方向上非常细长的单元，或者，更广泛地说，有一些单元的边缘既短又长。有很多方法可以给每个单元分配一个数字质量指数，以衡量该单元是 "好 "还是 "坏"；其中一些方法经常被选择，因为它们便宜且容易计算，而另一些方法则基于收敛性的证明。前者的一个例子是一个单元格的最长边与最短边的比率。在理想的情况下，这个比率应该是1；不好的单元格的数值远远大于1。后者的例子是考虑从参考单元 $\hat K=[0,1]^d$ 到实际单元 $K$ 的映射的梯度（"雅各布"）；这个梯度是一个矩阵，进入误差估计的一个量是参考单元上所有点的最大比率，即这个矩阵的最大和最小的特征值。同样不难看出，如果单元格 $K$ 是 $\hat K$ 的仿生图像，这个比率是恒定的，对于正方形和立方体来说，这个比率是一个。
* 在实践中，将这种质量度量可视化可能是很有趣的。函数 GridTools::compute_aspect_ratio_of_cells() 提供了获得这种信息的途径。甚至更好的是，可视化工具，如VisIt，通常允许你在可视化软件中对各种措施进行可视化；在VisIt的情况下，只需添加一个 "伪色 "图，并选择其中一个网格质量措施，而不是解域。
*

*<a name="Playingwiththeregularityofthesolution"></a><h4>Playing with the regularity of the solution</h4>。


*从数学的角度来看，拉普拉斯方程@f[
* 


* 
* 


* 
* 
-\Delta u = f
@f]在光滑有界的凸域上的解，已知本身是光滑的。平滑性的确切程度，即解所处的函数空间，取决于域的边界到底有多平滑，以及右手边的平滑程度。特别是如果右手边满足 $f\in C^\infty(\Omega)$ ，那么 $u \in C^\infty(\Omega_i)$ ，其中 $\Omega_i$ 是 $\Omega$ 的任何紧凑子集（ $\Omega$ 是一个开放域，所以紧凑子集需要与 $\partial\Omega$ 保持一个正距离），解的一些规律性可能会在边界上丢失，但一般来说，解是比右手边多一倍的可微积分。
* 然而，我们为目前的例子选择的情况是不同的：我们看一个具有非常数系数的方程 $a(\mathbf x)$ ：@f[
* 


* 
* 


* 
* 
-\nabla \cdot (a \nabla u) = f.
@f]这里，如果 $a$ 不光滑，那么无论 $f$ 如何，解也不会光滑。特别是，我们希望在 $a$ 沿直线（或沿三维平面）不连续的地方，解会有一个结点。这很容易看出来：例如，如果 $f$ 是连续的，那么 $f=-\nabla \cdot (a \nabla u)$ 也需要是连续的。这意味着 $a \nabla u$ 必须是可连续微分的（不存在扭结）。因此，如果 $a$ 有一个不连续，那么 $\nabla u$ 必须有一个相反的不连续，这样两者正好抵消，它们的乘积得到一个没有不连续的函数。但要使 $\nabla u$ 有一个不连续， $u$ 必须有一个结点。当然，这正是当前例子中所发生的情况，在解的图片中也很容易观察到。
* 一般来说，如果系数 $a(\mathbf x)$ 沿着2D的直线或3D的平面是不连续的，那么解可能会有一个结点，但解的梯度不会到无穷大。这意味着，该解至少仍在<a href="https://en.wikipedia.org/wiki/Sobolev_space">Sobolev space</a> $W^{1,\infty}$ 中（即，大致上，在导数有界的函数空间中）。另一方面，我们知道，在最极端的情况下
*
*--即域有重合角的情况下，右手边只满足 $f\in H^{-1}$ ，或者系数 $a$ 只在 $L^\infty$ 。
*
* - 我们所能期望的是， $u\in H^1$ （即导数是可平方整数的函数的<a
href="https://en.wikipedia.org/wiki/Sobolev_space#Sobolev_spaces_with_integer_k">Sobolev
space</a>），是一个比 $W^{1,\infty}$ 大得多的空间。创造一些案例并不难，在这些案例中，解在一个空间 $H^{1+s}$ 中，我们可以让 $s$ 变得像我们想要的那样小。这种情况经常被用来测试自适应有限元方法，因为网格必须解决导致解不在 $W^{1,\infty}$ 中的奇点。
* 典型的例子是<i>Kellogg problem</i>（指 @cite Kel74 ），在常用的形式中，它的系数 $a(\mathbf x)$ 在平面的四个象限内有不同的值（或在 ${\mathbb R}^3$ 的八个象限内有不同的值）。确切的规则性程度（上述索博列夫空间索引中的 $s$ ）取决于 $a(\mathbf x)$ 的值在原点处聚集，通过选择足够大的跳跃，可以使解的规则性尽可能地接近 $H^1$ 。
* 为了实现类似的目的，我们可以用下面的方法来代替系数函数（这里只显示2d的情况）。
* @code
template <int dim>
double coefficient (const Point<dim> &p)
{
  if ((p[0] < 0) && (p[1] < 0))           // lower left quadrant
    return 1;
  else if ((p[0] >= 0) && (p[1] < 0))     // lower right quadrant
    return 10;
  else if ((p[0] < 0) && (p[1] >= 0))     // upper left quadrant
    return 100;
  else if ((p[0] >= 0) && (p[1] >= 0))    // upper right quadrant
    return 1000;
  else
    {
      Assert(false, ExcInternalError());
      return 0;
    }
}
@endcode
*（在结尾处添加 <code>Assert</code> 可以确保在我们到达那个点时，要么抛出一个异常，要么程序中止。
*
* 当然我们不应该这样做，但这是一个很好的方法来保证自己：我们都会犯错误，有时没有想到所有的情况，例如检查 <code>p[0]</code> 是否小于和大于零，而不是大于或等于零，从而忘记了一些情况，否则会导致很难发现的错误。最后的 <code>return 0;</code> 只是为了避免编译器对函数没有以 <code>return</code> 语句结束的警告。
*
* - 编译器无法看到，由于前面的 <code>Assert</code> 语句，函数实际上永远不会到达这一点）。)
* 通过玩弄这种四个或更多扇形组合在一起，并且在这些扇形上的系数有不同的值的情况，我们可以构建出解在原点有奇异点的情况。我们还可以看到在这种情况下网格是如何被细化的。
*

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-6.cc"  。
*
*/


