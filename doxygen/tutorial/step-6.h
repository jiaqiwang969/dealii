/**
@page step_6 The step-6 tutorial program
This tutorial depends on step-5.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Whatadaptivelyrefinedmesheslooklike"> What adaptively refined meshes look like </a>
        <li><a href="#Whyadapativelyrefinedmeshes"> Why adapatively refined meshes? </a>
        <li><a href="#Howtodealwithhangingnodesintheory"> How to deal with hanging nodes in theory </a>
        <li><a href="#Howtodealwithhangingnodesinpractice"> How to deal with hanging nodes in practice </a>
        <li><a href="#Howweobtainlocallyrefinedmeshes"> How we obtain locally refined meshes </a>
        <li><a href="#Boundaryconditions"> Boundary conditions </a>
        <li><a href="#Otherthingsthisprogramshows"> Other things this program shows </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeStep6codeclasstemplate">The <code>Step6</code> class template</a>
        <li><a href="#Nonconstantcoefficients">Nonconstant coefficients</a>
        <li><a href="#ThecodeStep6codeclassimplementation">The <code>Step6</code> class implementation</a>
      <ul>
        <li><a href="#Step6Step6">Step6::Step6</a>
        <li><a href="#Step6setup_system">Step6::setup_system</a>
        <li><a href="#Step6assemble_system">Step6::assemble_system</a>
        <li><a href="#Step6solve">Step6::solve</a>
        <li><a href="#Step6refine_grid">Step6::refine_grid</a>
        <li><a href="#Step6output_results">Step6::output_results</a>
        <li><a href="#Step6run">Step6::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Solversandpreconditioners">Solvers and preconditioners</a>
        <li><a href="#Abettermesh">A better mesh</a>
        <li><a href="#Playingwiththeregularityofthesolution">Playing with the regularity of the solution</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-6/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{15,16,17,17.25,17.5,17.75} 

这个程序最后是关于deal.II的主要特征之一：使用自适应（局部）细化网格。这个程序仍然是基于步骤4和步骤5的，而且，正如你将看到的，实际上不需要花太多的代码来实现自适应性。事实上，虽然我们做了大量的解释，但自适应网格可以被添加到一个现有的程序中，几乎不需要十几行额外的代码。该程序显示了这些行是什么，以及自适应网格细化（AMR）的另一个重要成分：一个标准，可以用来确定是否有必要细化一个单元，因为它上面的误差很大，是否可以粗化这个单元，因为它上面的误差特别小，或者我们是否应该让这个单元保持原样。我们将在下文中讨论所有这些问题。




<a name="Whatadaptivelyrefinedmesheslooklike"></a><h3> What adaptively refined meshes look like </h3>


有许多方法可以自适应地细化网格。整个算法的基本结构总是相同的，由以下步骤的循环组成。

- 在当前网格上求解PDE。

- 用一些能说明误差的标准来估计每个单元格的误差。

- 把那些误差大的单元格标记为细化，把那些误差特别小的单元格标记为粗化，其余的就不用管了。

- 细化和粗化如此标记的单元，得到一个新的网格。

- 在新的网格上重复上述步骤，直到整体误差足够小。

由于一些可能被历史遗忘的原因（也许是这些函数过去是用FORTRAN语言实现的，这种语言并不关心某个东西是用小写字母还是大写字母拼写的，程序员经常习惯性地选择大写字母），上述循环在关于网格适应性的出版物中经常被称为SOLVE-ESTIMATE-MARK-REFINE循环（用这种拼法）。

然而，在这个结构之外，有多种方法可以实现这一点。从根本上说，它们的区别在于究竟如何从前一个网格中生成一个网格。

如果要使用三角形（deal.II没有这样做），那么就有两种基本的可能性。

- 最长边细化。在这个策略中，通过从最长边的中点到对面的顶点引入一条新的边，将一个标记为细化的三角形切成两段。当然，来自最长边的中点必须以某种方式通过*也*完善该边另一侧的单元格（如果有的话）来平衡。如果有问题的边也是相邻单元的最长边，那么我们可以直接运行一条新的边穿过相邻单元到对面的顶点；否则就需要一个稍微复杂的结构，在相邻单元的至少一条其他边上增加更多的新顶点，然后可能传播到相邻单元的邻居，直到算法终止。这很难用语言描述，而且因为deal.II不使用三角形，不值得在这里花时间。   但如果你很好奇，你可以随时在本介绍顶部显示的链接中观看视频讲座15。

- 红-绿细化。另一个选择是所谓的 "红绿细化"。   这种策略甚至更难描述（但在视频讲座中也讨论过），其优点是细化不会传播到我们想要细化的单元的近邻之外。然而，它的实施难度要大得多。

这些方法还有其他的变化，但重要的一点是，它们总是产生一个网格，其中两个单元的接触线是两个相邻单元的整个边缘。只要稍加努力，这种策略就可以很容易地适用于由四面体构成的三维网格。

这两种方法对2D的四边形和3D的六面体都不起作用，或者至少不容易。原因是要精化的四边形单元的四边形邻居所产生的过渡元素将是三角形，而我们不希望这样。因此，在deal.II中选择的适应性方法是使用网格，其中相邻的单元在细化水平上可能相差一个。这就导致在单元的界面上出现属于一方的节点，但在另一方是不平衡的。这些节点的通用术语是&ldquo;悬挂节点&rdquo;，这些网格在非常简单的情况下看起来是这样的。

 @image html hanging_nodes.png "A simple mesh with hanging nodes" 

一个更复杂的二维网格看起来是这样的（并在下面的 "结果 "部分讨论）。

<img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg" alt="第五个自适应细化的拉杜腾科网格：单元格沿着内圈聚拢。" width="300" height="300">

最后，这里展示了一个具有这种悬挂节点的三维网格（来自步骤-43）。

<img src="https://www.dealii.org/images/steps/developer/step-43.3d.mesh.png" alt="" width="300" height="300">

第一个和第三个网格当然是基于一个正方形和一个立方体，但正如第二个网格所显示的，这不是必要的。重要的一点是，我们可以独立于其邻居来细化一个网格（受制于一个单元只能比其邻居多细化一次的约束），但如果我们这样做，最终会出现这些&ldquo;悬空节点&rdquo;。




<a name="Whyadapativelyrefinedmeshes"></a><h3> Why adapatively refined meshes? </h3>


现在你已经看到了这些自适应细化网格的样子，你应该问<i>why</i>我们为什么要这样做。毕竟，我们从理论上知道，如果我们对网格进行全局细化，误差会下降到零，因为

@f{align*}{
  \|\nabla(u-u_h)\|_{\Omega} \le C h_\text{max}^p \| \nabla^{p+1} u \|_{\Omega},


@f}

其中 $C$ 是独立于 $h$ 和 $u$ 的一些常数， $p$ 是使用中的有限元的多项式程度， $h_\text{max}$ 是最大单元的直径。那么，如果<i>largest</i>单元很重要，那么为什么我们要在域的某些部分将网格做得很细，而不是全部？

答案在于观察到上面的公式不是最佳的。事实上，一些更多的工作表明，以下是一个更好的估计（你应该与上述估计的平方进行比较）。

@f{align*}{
  \|\nabla(u-u_h)\|_{\Omega}^2 \le C \sum_K h_K^{2p} \| \nabla^{p+1} u \|^2_K.


@f}

(因为 $h_K\le h_\text{max}$ ，如果你只是把网格大小从总和中拉出来，这个公式立即暗示了前一个公式)。这个公式所暗示的是，没有必要把<i>largest</i>单元格做得很小，而单元格真正只需要做小的<i>where $\| \nabla^{p+1} u \|_K$ is large</i>!换句话说。网格实际上只需要在解有较大变化的地方做得很细，正如 $p+1$ st导数所表明的。这是有直观意义的：例如，如果我们使用一个线性元素 $p=1$ ，那么即使网格很粗，那些解几乎是线性的地方（如 $\nabla^2 u$ 所示的小地方）也会被很好地解决。只有那些二阶导数大的地方才会被大元素解决得很差，因此我们应该把网格做得很小。

当然，这个<i>a priori estimate</i>在实践中不是很有用，因为我们不知道问题的精确解 $u$ ，因此，我们不能计算 $\nabla^{p+1}u$  。但是，这也是通常采取的方法，我们可以只根据之前计算的离散解 $u_h$ 来计算 $\nabla^{p+1}u$ 的数值近似值。我们将在下面稍微详细地讨论这个问题。这将有助于我们确定哪些单元具有较大的 $p+1$ st导数，然后这些单元将成为细化网格的候选单元。




<a name="Howtodealwithhangingnodesintheory"></a><h3> How to deal with hanging nodes in theory </h3>


上面提到的使用三角形网格的方法，都是为了确保每个顶点都是所有相邻单元的顶点--也就是说，没有悬空节点。这就自动确保了我们能够以这样的方式定义形状函数，即它们是全局连续的（如果我们使用到目前为止在教程程序中一直使用的常见的 $Q_p$ 拉格朗日有限元方法，如FE_Q类所代表的）。

另一方面，如果我们在有悬挂节点的网格上定义形状函数，我们最终可能得到不连续的形状函数。要看到这一点，请想一下上面的情况，即右上角的单元没有被细化，并考虑一下使用双线性有限元的情况。在这种情况下，与悬挂节点相关的形状函数是以明显的方式定义在与每个悬挂节点相邻的两个小单元上。但我们如何将它们扩展到相邻的大单元呢？显然，函数对大单元的扩展不能是双线性的，因为那样的话，它需要沿着大单元的每条边线性化，这意味着它在整条边上需要为零，因为它需要在大单元的两个顶点上为零。但从小单元一侧看，它在悬挂节点本身并不是零--所以它不是连续的。下面三幅图显示了沿着有关边缘的三个形状函数，当以通常的方式简单地根据它们相邻的单元格来定义时，这些形状函数变成了不连续的。

<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center">  @image html hanging_nodes_shape_functions_1.png "A discontinuous shape function adjacent to a hanging node"  </div> </div> <div class="parent"> <div class="img" align="center">  @image html hanging_nodes_shape_functions_2.png "A discontinuous shape function at a hanging node"  </div> </div> <div class="parent"> <div class="img" align="center">  @image html hanging_nodes_shape_functions_3.png "A discontinuous shape function adjacent to a hanging node"  </div> </div></div>


但我们确实希望有限元解是连续的，这样我们就有了&ldquo;符合要求的有限元方法&rdquo;，其中离散有限元空间是我们寻求拉普拉斯方程解的 $H^1$ 函数空间的一个适当子集。为了保证全局解在这些节点上也是连续的，我们必须对这些节点上的解的值提出一些额外的约束。诀窍是要认识到，虽然上面显示的形状函数是不连续的（因此它们的<i>arbitrary</i>线性组合也是不连续的），但形状函数加起来为 $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf x)$ 的线性组合可以是连续的<i>if the coefficients $U_j$ satisfy certain relationships</i>。换句话说，系数 $U_j$ 不能任意选择，而必须满足某些约束条件，这样，函数 $u_h$ 实际上是连续的。这些约束条件在概念上相对容易理解，但在软件中的实现却很复杂，需要几千行的代码。另一方面，在用户代码中，在处理挂起的节点时，你只需要添加大约半打的行。

在下面的程序中，我们将展示如何从deal.II中获得这些约束，以及如何在线性方程组的求解中使用它们。在了解下面程序的细节之前，你可能想看看 @ref constraints 文件模块，它解释了这些约束如何计算以及deal.II中哪些类对它们起作用。




<a name="Howtodealwithhangingnodesinpractice"></a><h3> How to deal with hanging nodes in practice </h3>


悬挂节点约束的实践比我们上面概述的理论更简单。实际上，你只需要在step-4这样的程序中增加半打额外的代码，就可以使它在有悬挂节点的自适应网格中工作。有趣的是，这与你要解决的方程完全无关。这些约束的代数性质与方程无关，只取决于对有限元的选择。因此，处理这些约束的代码完全包含在deal.II库本身，你不需要担心细节问题。

你需要使其发挥作用的步骤基本上是这样的。

- 你必须创建一个AffineConstraints对象，（顾名思义）它将存储有限元空间的所有约束。在目前的情况下，这些约束是由于我们希望保持解空间的连续，甚至在有悬空节点的情况下。(下面我们还将简要地提到，我们还将把边界值放到这个对象中，但这是一个单独的问题)。

- 你必须使用函数 DoFTools::make_hanging_node_constraints() 来填充这个对象，以确保有限元空间的元素的连续性。

- 当你通过使用 AffineConstraints::distribute_local_to_global(). 将矩阵和右手边的局部贡献复制到全局对象时，你必须使用这个对象。 到目前为止，我们已经自己完成了这个工作，但现在有了约束，这就是神奇的地方，我们将约束应用到线性系统中。这个函数所做的是确保位于悬空节点的自由度事实上不是真正的自由。相反，通过将它们的行和列设置为零，并在对角线上放置一些东西以确保矩阵保持可反转，它们实际上被从线性系统中消除了。   对于我们在这里解决的拉普拉斯方程来说，这个过程产生的矩阵仍然是对称和正定的，所以我们可以继续使用共轭梯度法来解决。

- 然后你像往常一样求解线性系统，但在这一步结束时，你需要确保位于悬挂节点上的 "自由度 "得到正确的（约束的）值，这样你随后可视化的或以其他方式评估的解决方案实际上是连续的。这可以通过在求解后立即调用 AffineConstraints::distribute() 来实现。

这四个步骤实际上是所有必要的--从用户的角度来看就是这么简单。事实上，在上面提到的函数调用中，你将运行几千行并不复杂的代码，这一点完全不重要。在用户代码中，实际上只有四个额外的步骤。




<a name="Howweobtainlocallyrefinedmeshes"></a><h3> How we obtain locally refined meshes </h3>


下一个问题是，既然我们知道如何<i>deal</i>处理有这些悬挂节点的网格，那么我们如何<i>obtain</i>它们。

一个简单的方法已经在步骤1中展示过了：如果你<i>know</i>哪里需要细化网格，那么你可以手工创建一个。但是在现实中，我们并不知道这些。我们不知道PDE的解在前面（因为，如果我们知道，我们就不必使用有限元方法），因此，我们不知道哪里需要增加局部网格细化来更好地解决解有强烈变化的区域。但是上面的讨论表明，也许我们可以用一个网格上的离散解 $u_h$ 来估计导数 $\nabla^{p+1} u$ ，然后用这个来确定哪些单元太大，哪些已经足够小。然后，我们可以使用局部网格细化技术从当前的网格中生成一个新的网格。如果有必要，这个步骤会重复进行，直到我们对我们的数值解决方案感到满意--或者，更常见的是，直到我们耗尽了计算资源或耐心。

所以这正是我们要做的。局部细化网格是使用一个<i>error estimator</i>产生的，它可以估计拉普拉斯算子的数值解的能量误差。由于它是由Kelly和他的同事开发的，我们经常在库、文档和邮件列表中把它称为&ldquo;Kelly细化指标&rdquo;。实现它的类被称为KellyErrorEstimator，在该类的文档中可以找到大量的信息，这里不需要重复。然而，总结起来就是，该类计算出一个具有与 @ref GlossActive "活动单元 "一样多的条目的向量，其中每个条目包含对该单元的误差估计。这个估计值然后被用来细化网格的单元：那些有大误差的单元将被标记为细化，那些有特别小估计值的单元将被标记为粗化。我们不需要用手去做这些。一旦我们获得了误差估计矢量，命名空间GridRefinement中的函数将为我们完成这一切。

值得注意的是，虽然Kelly误差估计器是为拉普拉斯方程开发的，但它已被证明是为广泛的方程生成局部细化网格的合适工具，甚至不限于只针对椭圆问题。尽管它对其他方程会产生非最优网格，但它往往是快速产生网格的好方法，能很好地适应解的特征，如大变化区域或不连续性。




<a name="Boundaryconditions"></a><h3> Boundary conditions </h3>


事实证明，人们可以把迪里希特边界条件看作是对自由度的另一种约束。这的确是一个特别简单的约束。如果 $j$ 是边界上的一个自由度，其位置为 $\mathbf x_j$ ，那么在 $\partial\Omega$ 上施加边界条件 $u=g$ 就会产生约束 $U_j=g({\mathbf x}_j)$  。

AffineConstraints类也可以处理这样的约束，这使得我们可以方便地让我们用于悬挂节点约束的同一个对象也处理这些Dirichlet边界条件。这样一来，我们就不需要在装配后应用边界条件（就像我们在前面的步骤中做的那样）。所有需要的是我们调用 VectorTools::interpolate_boundary_values() 的变体，该变体在AffineConstraints对象中返回其信息，而不是我们在以前的教程程序中使用的 `std::map` 。


<a name="Otherthingsthisprogramshows"></a><h3> Other things this program shows </h3>



由于用于局部细化网格的概念非常重要，我们在这个例子中没有展示很多其他材料。最重要的例外是，我们展示了如何使用双二次元而不是之前所有例子中使用的双线性元素。事实上，使用高阶元素只需替换程序中的三行，即在本程序主类的构造函数中初始化 <code>fe</code> 成员变量，以及在两个地方使用适当的正交公式。程序的其他部分没有变化。

其他唯一的新东西是在 <code>main</code> 函数中捕捉异常的方法，以便在程序因某种原因崩溃时输出一些信息。下面将详细讨论这个问题。


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
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * #include <fstream>
 * 
 * @endcode
 * 
 * From the following include file we will import the declaration of
 * H1-conforming finite element shape functions. This family of finite
 * elements is called <code>FE_Q</code>, and was used in all examples before
 * already to define the usual bi- or tri-linear elements, but we will now use
 * it for bi-quadratic elements:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * @endcode
 * 
 * We will not read the grid from a file as in the previous example, but
 * generate it using a function of the library. However, we will want to write
 * out the locally refined grids (just the grid, not the solution) in each
 * step, so we need the following include file instead of
 * <code>grid_in.h</code>:
 * 
 * @code
 * #include <deal.II/grid/grid_out.h>
 * 
 * 
 * @endcode
 * 
 * When using locally refined grids, we will get so-called <code>hanging
 * nodes</code>. However, the standard finite element methods assumes that the
 * discrete solution spaces be continuous, so we need to make sure that the
 * degrees of freedom on hanging nodes conform to some constraints such that
 * the global solution is continuous. We are also going to store the boundary
 * conditions in this object. The following file contains a class which is
 * used to handle these constraints:
 * 
 * @code
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * @endcode
 * 
 * In order to refine our grids locally, we need a function from the library
 * that decides which cells to flag for refinement or coarsening based on the
 * error indicators we have computed. This function is defined here:
 * 
 * @code
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * @endcode
 * 
 * Finally, we need a simple way to actually compute the refinement indicators
 * based on some error estimate. While in general, adaptivity is very
 * problem-specific, the error indicator in the following file often yields
 * quite nicely adapted grids for a wide class of problems.
 * 
 * @code
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * Finally, this is as in previous programs:
 * 
 * @code
 * using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep6codeclasstemplate"></a> 
 * <h3>The <code>Step6</code> class template</h3>
 * 

 * 
 * The main class is again almost unchanged. Two additions, however, are made:
 * we have added the <code>refine_grid</code> function, which is used to
 * adaptively refine the grid (instead of the global refinement in the
 * previous examples), and a variable which will hold the constraints.
 * 
 * @code
 * template <int dim>
 * class Step6
 * {
 * public:
 *   Step6();
 * 
 *   void run();
 * 
 * private:
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void refine_grid();
 *   void output_results(const unsigned int cycle) const;
 * 
 *   Triangulation<dim> triangulation;
 * 
 *   FE_Q<dim>       fe;
 *   DoFHandler<dim> dof_handler;
 * 
 * 
 * @endcode
 * 
 * This is the new variable in the main class. We need an object which holds
 * a list of constraints to hold the hanging nodes and the boundary
 * conditions.
 * 
 * @code
 *   AffineConstraints<double> constraints;
 * 
 *   SparseMatrix<double> system_matrix;
 *   SparsityPattern      sparsity_pattern;
 * 
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Nonconstantcoefficients"></a> 
 * <h3>Nonconstant coefficients</h3>
 * 

 * 
 * The implementation of nonconstant coefficients is copied verbatim from
 * step-5:
 * 
 * @code
 * template <int dim>
 * double coefficient(const Point<dim> &p)
 * {
 *   if (p.square() < 0.5 * 0.5)
 *     return 20;
 *   else
 *     return 1;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep6codeclassimplementation"></a> 
 * <h3>The <code>Step6</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Step6Step6"></a> 
 * <h4>Step6::Step6</h4>
 * 

 * 
 * The constructor of this class is mostly the same as before, but this time
 * we want to use the quadratic element. To do so, we only have to replace the
 * constructor argument (which was <code>1</code> in all previous examples) by
 * the desired polynomial degree (here <code>2</code>):
 * 
 * @code
 * template <int dim>
 * Step6<dim>::Step6()
 *   : fe(2)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6setup_system"></a> 
 * <h4>Step6::setup_system</h4>
 * 

 * 
 * The next function sets up all the variables that describe the linear
 * finite element problem, such as the DoFHandler, matrices, and
 * vectors. The difference to what we did in step-5 is only that we now also
 * have to take care of hanging node constraints. These constraints are
 * handled almost exclusively by the library, i.e. you only need to know
 * that they exist and how to get them, but you do not have to know how they
 * are formed or what exactly is done with them.
 * 

 * 
 * At the beginning of the function, you find all the things that are the same
 * as in step-5: setting up the degrees of freedom (this time we have
 * quadratic elements, but there is no difference from a user code perspective
 * to the linear -- or any other degree, for that matter -- case), generating
 * the sparsity pattern, and initializing the solution and right hand side
 * vectors. Note that the sparsity pattern will have significantly more
 * entries per row now, since there are now 9 degrees of freedom per cell
 * (rather than only four), that can couple with each other.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 * 
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * We may now populate the AffineConstraints object with the hanging node
 * constraints. Since we will call this function in a loop we first clear
 * the current set of constraints from the last system and then compute new
 * ones:
 * 
 * @code
 *   constraints.clear();
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 * 
 * @endcode
 * 
 * Now we are ready to interpolate the boundary values with indicator 0 (the
 * whole boundary) and store the resulting constraints in our
 * <code>constraints</code> object. Note that we do not to apply the
 * boundary conditions after assembly, like we did in earlier steps: instead
 * we put all constraints on our function space in the AffineConstraints
 * object. We can add constraints to the AffineConstraints object in either
 * order: if two constraints conflict then the constraint matrix either abort
 * or throw an exception via the Assert macro.
 * 
 * @code
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<dim>(),
 *                                            constraints);
 * 
 * @endcode
 * 
 * After all constraints have been added, they need to be sorted and
 * rearranged to perform some actions more efficiently. This postprocessing
 * is done using the <code>close()</code> function, after which no further
 * constraints may be added any more:
 * 
 * @code
 *   constraints.close();
 * 
 * @endcode
 * 
 * Now we first build our compressed sparsity pattern like we did in the
 * previous examples. Nevertheless, we do not copy it to the final sparsity
 * pattern immediately.  Note that we call a variant of
 * make_sparsity_pattern that takes the AffineConstraints object as the third
 * argument. We are letting the routine know that we will never write into
 * the locations given by <code>constraints</code> by setting the argument
 * <code>keep_constrained_dofs</code> to false (in other words, that we will
 * never write into entries of the matrix that correspond to constrained
 * degrees of freedom). If we were to condense the
 * constraints after assembling, we would have to pass <code>true</code>
 * instead because then we would first write into these locations only to
 * later set them to zero again during condensation.
 * 
 * @code
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler,
 *                                   dsp,
 *                                   constraints,
 *                                   /*keep_constrained_dofs = */ false);
 * 
 * @endcode
 * 
 * Now all non-zero entries of the matrix are known (i.e. those from
 * regularly assembling the matrix and those that were introduced by
 * eliminating constraints). We may copy our intermediate object to the
 * sparsity pattern:
 * 
 * @code
 *   sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * We may now, finally, initialize the sparse matrix:
 * 
 * @code
 *   system_matrix.reinit(sparsity_pattern);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6assemble_system"></a> 
 * <h4>Step6::assemble_system</h4>
 * 

 * 
 * Next, we have to assemble the matrix. However, to copy the local matrix and
 * vector on each cell into the global system, we are no longer using a
 * hand-written loop. Instead, we use
 * AffineConstraints::distribute_local_to_global() that internally executes
 * this loop while performing Gaussian elimination on rows and columns
 * corresponding to constrained degrees on freedom.
 * 

 * 
 * The rest of the code that forms the local contributions remains
 * unchanged. It is worth noting, however, that under the hood several things
 * are different than before. First, the variable <code>dofs_per_cell</code>
 * and return value of <code>quadrature_formula.size()</code> now are 9 each,
 * where they were 4 before. Introducing such variables as abbreviations is a
 * good strategy to make code work with different elements without having to
 * change too much code. Secondly, the <code>fe_values</code> object of course
 * needs to do other things as well, since the shape functions are now
 * quadratic, rather than linear, in each coordinate variable. Again, however,
 * this is something that is completely handled by the library.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::assemble_system()
 * {
 *   const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *   FEValues<dim> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients |
 *                             update_quadrature_points | update_JxW_values);
 * 
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 *       fe_values.reinit(cell);
 * 
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         {
 *           const double current_coefficient =
 *             coefficient(fe_values.quadrature_point(q_index));
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               for (const unsigned int j : fe_values.dof_indices())
 *                 cell_matrix(i, j) +=
 *                   (current_coefficient *              // a(x_q)
 *                    fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                    fe_values.JxW(q_index));           // dx
 * 
 *               cell_rhs(i) += (1.0 *                               // f(x)
 *                               fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                               fe_values.JxW(q_index));            // dx
 *             }
 *         }
 * 
 * @endcode
 * 
 * Finally, transfer the contributions from @p cell_matrix and
 * @p cell_rhs into the global objects.
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices);
 *       constraints.distribute_local_to_global(
 *         cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *     }
 * @endcode
 * 
 * Now we are done assembling the linear system. The constraint matrix took
 * care of applying the boundary conditions and also eliminated hanging node
 * constraints. The constrained nodes are still in the linear system (there
 * is a nonzero entry, chosen in a way that the matrix is well conditioned,
 * on the diagonal of the matrix and all other entries for this line are set
 * to zero) but the computed values are invalid (i.e., the corresponding
 * entries in <code>system_rhs</code> are currently meaningless). We compute
 * the correct values for these nodes at the end of the <code>solve</code>
 * function.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6solve"></a> 
 * <h4>Step6::solve</h4>
 * 

 * 
 * We continue with gradual improvements. The function that solves the linear
 * system again uses the SSOR preconditioner, and is again unchanged except
 * that we have to incorporate hanging node constraints. As mentioned above,
 * the degrees of freedom from the AffineConstraints object corresponding to
 * hanging node constraints and boundary values have been removed from the
 * linear system by giving the rows and columns of the matrix a special
 * treatment. This way, the values for these degrees of freedom have wrong,
 * but well-defined values after solving the linear system. What we then have
 * to do is to use the constraints to assign to them the values that they
 * should have. This process, called <code>distributing</code> constraints,
 * computes the values of constrained nodes from the values of the
 * unconstrained ones, and requires only a single additional function call
 * that you find at the end of this function:
 * 

 * 
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::solve()
 * {
 *   SolverControl            solver_control(1000, 1e-12);
 *   SolverCG<Vector<double>> solver(solver_control);
 * 
 *   PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *   preconditioner.initialize(system_matrix, 1.2);
 * 
 *   solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *   constraints.distribute(solution);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6refine_grid"></a> 
 * <h4>Step6::refine_grid</h4>
 * 

 * 
 * We use a sophisticated error estimation scheme to refine the mesh instead
 * of global refinement. We will use the KellyErrorEstimator class which
 * implements an error estimator for the Laplace equation; it can in principle
 * handle variable coefficients, but we will not use these advanced features,
 * but rather use its most simple form since we are not interested in
 * quantitative results but only in a quick way to generate locally refined
 * grids.
 * 

 * 
 * Although the error estimator derived by Kelly et al. was originally
 * developed for the Laplace equation, we have found that it is also well
 * suited to quickly generate locally refined grids for a wide class of
 * problems. This error estimator uses the solution gradient's jump at
 * cell faces (which is a measure for the second derivatives) and
 * scales it by the size of the cell. It is therefore a measure for the local
 * smoothness of the solution at the place of each cell and it is thus
 * understandable that it yields reasonable grids also for hyperbolic
 * transport problems or the wave equation as well, although these grids are
 * certainly suboptimal compared to approaches specially tailored to the
 * problem. This error estimator may therefore be understood as a quick way to
 * test an adaptive program.
 * 

 * 
 * The way the estimator works is to take a <code>DoFHandler</code> object
 * describing the degrees of freedom and a vector of values for each degree of
 * freedom as input and compute a single indicator value for each active cell
 * of the triangulation (i.e. one value for each of the active cells). To do
 * so, it needs two additional pieces of information: a face quadrature formula,
 * i.e., a quadrature formula on <code>dim-1</code> dimensional objects. We use
 * a 3-point Gauss rule again, a choice that is consistent and appropriate with
 * the bi-quadratic finite element shape functions in this program.
 * (What constitutes a suitable quadrature rule here of course depends on
 * knowledge of the way the error estimator evaluates the solution field. As
 * said above, the jump of the gradient is integrated over each face, which
 * would be a quadratic function on each face for the quadratic elements in
 * use in this example. In fact, however, it is the square of the jump of the
 * gradient, as explained in the documentation of that class, and that is a
 * quartic function, for which a 3 point Gauss formula is sufficient since it
 * integrates polynomials up to order 5 exactly.)
 * 

 * 
 * Secondly, the function wants a list of boundary indicators for those
 * boundaries where we have imposed Neumann values of the kind
 * $\partial_n u(\mathbf x) = h(\mathbf x)$, along with a function $h(\mathbf
 * x)$ for each such boundary. This information is represented by a map from
 * boundary indicators to function objects describing the Neumann boundary
 * values. In the present example program, we do not use Neumann boundary
 * values, so this map is empty, and in fact constructed using the default
 * constructor of the map in the place where the function call expects the
 * respective function argument.
 * 

 * 
 * The output is a vector of values for all active cells. While it may
 * make sense to compute the <b>value</b> of a solution degree of freedom
 * very accurately, it is usually not necessary to compute the <b>error
 * indicator</b> corresponding to the solution on a cell particularly
 * accurately. We therefore typically use a vector of floats instead of a vector
 * of doubles to represent error indicators.
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::refine_grid()
 * {
 *   Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *   KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                      QGauss<dim - 1>(fe.degree + 1),
 *                                      {},
 *                                      solution,
 *                                      estimated_error_per_cell);
 * 
 * @endcode
 * 
 * The above function returned one error indicator value for each cell in
 * the <code>estimated_error_per_cell</code> array. Refinement is now done
 * as follows: refine those 30 per cent of the cells with the highest error
 * values, and coarsen the 3 per cent of cells with the lowest values.
 *   

 * 
 * One can easily verify that if the second number were zero, this would
 * approximately result in a doubling of cells in each step in two space
 * dimensions, since for each of the 30 per cent of cells, four new would be
 * replaced, while the remaining 70 per cent of cells remain untouched. In
 * practice, some more cells are usually produced since it is disallowed
 * that a cell is refined twice while the neighbor cell is not refined; in
 * that case, the neighbor cell would be refined as well.
 *   

 * 
 * In many applications, the number of cells to be coarsened would be set to
 * something larger than only three per cent. A non-zero value is useful
 * especially if for some reason the initial (coarse) grid is already rather
 * refined. In that case, it might be necessary to refine it in some
 * regions, while coarsening in some other regions is useful. In our case
 * here, the initial grid is very coarse, so coarsening is only necessary in
 * a few regions where over-refinement may have taken place. Thus a small,
 * non-zero value is appropriate here.
 *   

 * 
 * The following function now takes these refinement indicators and flags
 * some cells of the triangulation for refinement or coarsening using the
 * method described above. It is from a class that implements several
 * different algorithms to refine a triangulation based on cell-wise error
 * indicators.
 * 
 * @code
 *   GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                   estimated_error_per_cell,
 *                                                   0.3,
 *                                                   0.03);
 * 
 * @endcode
 * 
 * After the previous function has exited, some cells are flagged for
 * refinement, and some other for coarsening. The refinement or coarsening
 * itself is not performed by now, however, since there are cases where
 * further modifications of these flags is useful. Here, we don't want to do
 * any such thing, so we can tell the triangulation to perform the actions
 * for which the cells are flagged:
 * 
 * @code
 *   triangulation.execute_coarsening_and_refinement();
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6output_results"></a> 
 * <h4>Step6::output_results</h4>
 * 

 * 
 * At the end of computations on each grid, and just before we continue the
 * next cycle with mesh refinement, we want to output the results from this
 * cycle.
 * 

 * 
 * We have already seen in step-1 how this can be achieved for the
 * mesh itself. Here, we change a few things:
 * <ol>
 * <li>We use two different formats: gnuplot and VTU.</li>
 * <li>We embed the cycle number in the output file name.</li>
 * <li>For gnuplot output, we set up a GridOutFlags::Gnuplot object to
 * provide a few extra visualization arguments so that edges appear
 * curved. This is explained in further detail in step-10.</li>
 * </ol>
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::output_results(const unsigned int cycle) const
 * {
 *   {
 *     GridOut               grid_out;
 *     std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot");
 *     GridOutFlags::Gnuplot gnuplot_flags(false, 5);
 *     grid_out.set_flags(gnuplot_flags);
 *     MappingQGeneric<dim> mapping(3);
 *     grid_out.write_gnuplot(triangulation, output, &mapping);
 *   }
 * 
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches();
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
 *     data_out.write_vtu(output);
 *   }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step6run"></a> 
 * <h4>Step6::run</h4>
 * 

 * 
 * The final function before <code>main()</code> is again the main driver of
 * the class, <code>run()</code>. It is similar to the one of step-5, except
 * that we generate a file in the program again instead of reading it from
 * disk, in that we adaptively instead of globally refine the mesh, and that
 * we output the solution on the final mesh in the present function.
 * 

 * 
 * The first block in the main loop of the function deals with mesh generation.
 * If this is the first cycle of the program, instead of reading the grid from
 * a file on disk as in the previous example, we now again create it using a
 * library function. The domain is again a circle with center at the origin and
 * a radius of one (these are the two hidden arguments to the function, which
 * have default values).
 * 

 * 
 * You will notice by looking at the coarse grid that it is of inferior
 * quality than the one which we read from the file in the previous example:
 * the cells are less equally formed. However, using the library function this
 * program works in any space dimension, which was not the case before.
 * 

 * 
 * In case we find that this is not the first cycle, we want to refine the
 * grid. Unlike the global refinement employed in the last example program, we
 * now use the adaptive procedure described above.
 * 

 * 
 * The rest of the loop looks as before:
 * 
 * @code
 * template <int dim>
 * void Step6<dim>::run()
 * {
 *   for (unsigned int cycle = 0; cycle < 8; ++cycle)
 *     {
 *       std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *       if (cycle == 0)
 *         {
 *           GridGenerator::hyper_ball(triangulation);
 *           triangulation.refine_global(1);
 *         }
 *       else
 *         refine_grid();
 * 
 * 
 *       std::cout << "   Number of active cells:       "
 *                 << triangulation.n_active_cells() << std::endl;
 * 
 *       setup_system();
 * 
 *       std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *                 << std::endl;
 * 
 *       assemble_system();
 *       solve();
 *       output_results(cycle);
 *     }
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * The main function is unaltered in its functionality from the previous
 * example, but we have taken a step of additional caution. Sometimes,
 * something goes wrong (such as insufficient disk space upon writing an
 * output file, not enough memory when trying to allocate a vector or a
 * matrix, or if we can't read from or write to a file for whatever reason),
 * and in these cases the library will throw exceptions. Since these are
 * run-time problems, not programming errors that can be fixed once and for
 * all, this kind of exceptions is not switched off in optimized mode, in
 * contrast to the <code>Assert</code> macro which we have used to test
 * against programming errors. If uncaught, these exceptions propagate the
 * call tree up to the <code>main</code> function, and if they are not caught
 * there either, the program is aborted. In many cases, like if there is not
 * enough memory or disk space, we can't do anything but we can at least print
 * some text trying to explain the reason why the program failed. A way to do
 * so is shown in the following. It is certainly useful to write any larger
 * program in this way, and you can do so by more or less copying this
 * function except for the <code>try</code> block that actually encodes the
 * functionality particular to the present application.
 * 
 * @code
 * int main()
 * {
 * @endcode
 * 
 * The general idea behind the layout of this function is as follows: let's
 * try to run the program as we did before...
 * 
 * @code
 *   try
 *     {
 *       Step6<2> laplace_problem_2d;
 *       laplace_problem_2d.run();
 *     }
 * @endcode
 * 
 * ...and if this should fail, try to gather as much information as
 * possible. Specifically, if the exception that was thrown is an object of
 * a class that is derived from the C++ standard class
 * <code>exception</code>, then we can use the <code>what</code> member
 * function to get a string which describes the reason why the exception was
 * thrown.
 *   

 * 
 * The deal.II exception classes are all derived from the standard class,
 * and in particular, the <code>exc.what()</code> function will return
 * approximately the same string as would be generated if the exception was
 * thrown using the <code>Assert</code> macro. You have seen the output of
 * such an exception in the previous example, and you then know that it
 * contains the file and line number of where the exception occurred, and
 * some other information. This is also what the following statements would
 * print.
 *   

 * 
 * Apart from this, there isn't much that we can do except exiting the
 * program with an error code (this is what the <code>return 1;</code>
 * does):
 * 
 * @code
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
 * @endcode
 * 
 * If the exception that was thrown somewhere was not an object of a class
 * derived from the standard <code>exception</code> class, then we can't do
 * anything at all. We then simply print an error message and exit.
 * 
 * @code
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
 * @endcode
 * 
 * If we got to this point, there was no exception which propagated up to
 * the main function (there may have been exceptions, but they were caught
 * somewhere in the program or the library). Therefore, the program
 * performed as was expected and we can return without error.
 * 
 * @code
 *   return 0;
 * }
 * @endcode
examples/step-6/doc/results.dox



<a name="Results"></a><h1>Results</h1>



程序的输出看起来如下。

@code
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






正如预期的那样，在每个周期中，单元格的数量大约增加了一倍。度数略多于单元数的四倍；人们期望在无限网格的两个空间维度上的系数正好是四（因为自由度之间的间隔是单元宽度的一半：每个边缘有一个额外的自由度，每个单元的中间有一个），但由于网格的有限尺寸和由悬挂节点和局部细化引入的额外自由度，它大于这个系数。




程序在细化循环的每个周期都输出解决方案和网格。解决方案看起来如下。

 <img src="https://www.dealii.org/images/steps/developer/step-6.solution.9.2.png" alt=""> 

关注该程序如何得出最终的网格是很有趣的。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0.svg" alt="初始网格：具有一个全局细化的五格圆形网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1.svg" alt="第一个网格：具有两个全局细化的五格圆形网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2.svg" alt="第二个网格：有一个自适应细化的五格圆形网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3.svg" alt="第三个网格：有两个自适应细化的五格圆形网格，显示围绕内圆的聚类。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4.svg" alt="第四个网格：具有三个自适应细化的五格圆形网格，显示了围绕内圈的聚类。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5.svg" alt="第五个网格：具有四个自适应细化的五格圆形网格，显示了围绕内圈的聚类。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_6.svg" alt="第六个网格：具有五个自适应细化的五格圆形网格，显示了围绕内圈的聚类。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_7.svg" alt="最后一个网格：有六个自适应细化的五格圆形网格，显示大多数单元都聚集在内圆周围。" width="300" height="300"> </div> </div>


可以清楚地看到，在解有扭结的区域，也就是离中心0.5的径向距离的圆，被精炼得最多。此外，解非常光滑和几乎平坦的中心区域几乎完全没有被细化，但这是由于我们没有考虑到那里的系数很大的事实。外面的区域被任意细化，因为那里的二阶导数是恒定的，因此细化主要是基于单元的大小和它们与最佳方形的偏差。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Solversandpreconditioners"></a><h4>Solvers and preconditioners</h4>



如果一个人要解决相当大的问题（比我们这里的问题大得多），有一件事总是值得一试的，那就是尝试不同的求解器或预处理器。在目前的情况下，线性系统是对称的和正定的，这使得CG算法几乎成了求解的典型选择。然而，我们在 <code>solve()</code> 函数中使用的SSOR预处理器是可以争夺的。

在deal.II中，改变预处理程序是比较简单的。例如，通过改变现有的几行代码

@code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
@endcode

进入

@code
  PreconditionSSOR<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);
@endcode

我们可以尝试SSOR的不同放松参数。通过使用

@code
  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode

我们可以使用Jacobi作为预处理程序。而通过使用

@code
  SparseILU<double> preconditioner;
  preconditioner.initialize(system_matrix);
@endcode

我们可以使用一个简单的不完全LU分解，不需要任何阈值处理或加强对角线（要使用这个预处理程序，你还必须把头文件 <code>deal.II/lac/sparse_ilu.h</code> 添加到文件顶部的包含列表中）。

使用这些不同的预处理程序，我们可以比较所需的CG迭代次数（可通过 <code>solver_control.last_step()</code> 调用，见步骤4）以及所需的CPU时间（使用Timer类，例如在步骤28中讨论的），得到如下结果（左：迭代次数；右：CPU时间）。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q2.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table> 

我们可以看到，在这个简单的问题上，所有的预处理程序的表现都差不多，迭代次数的增长是 ${\cal
O}(N^{1/2})$ ，由于每次迭代需要大约 ${\cal
O}(N)$ 次操作，总的CPU时间增长是 ${\cal
O}(N^{3/2})$ （对于几个最小的网格，CPU时间小到没有记录）。请注意，尽管它是最简单的方法，但对于这个问题，雅可比是最快的。

当有限元不是本程序构造函数中设定的双二次元，而是双线性的时候，情况会有一些变化。如果做此改变，结果如下。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_iterations.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-6.q1.dofs_vs_time.png" alt="">
    </td>
  </tr>
</table> 

换句话说，虽然迭代次数和CPU时间的增加与以前一样，但Jacobi现在是需要迭代次数最多的方法；不过，由于它必须执行的操作很简单，所以它仍然是最快的方法。这并不是说Jacobi实际上是一个好的预处理方法--对于规模可观的问题来说，它绝对不是，其他方法会好得多--而实际上只是因为它的实现非常简单，可以补偿更多的迭代次数，所以它的速度很快。

从这里得到的信息并不是预处理程序的简单性总是最好的。虽然这对目前的问题可能是正确的，但一旦我们转向更复杂的问题（弹性或斯托克斯，例如第8步或第22步），就绝对不是这样了。其次，所有这些预处理程序仍然会导致迭代次数随着自由度数 $N$ 的增加而增加，例如 ${\cal O}(N^\alpha)$ ；这反过来又会导致总工作量增加为 ${\cal O}(N^{1+\alpha})$ ，因为每次迭代都需要 ${\cal O}(N)$ 的工作。这种行为是不可取的：我们真的希望用 $N$ 个未知数解决线性系统，总工作量为 ${\cal O}(N)$ 个；有一类预处理程序可以实现这一点，即几何（step-16、step-37、step-39）或代数多网格（step-31、step-40和其他几个）预处理程序。然而，它们要比上述的预处理程序复杂得多。

最后，要带回家的最后一个信息是，当上面显示的数据产生时（2018年），有10万个未知数的线性系统在台式机上很容易在大约一秒钟内解决，使相对简单的2d问题的解决甚至达到非常高的精度，甚至在过去也不是一个大任务。当时，三维问题的情况完全不同，但即使是这样，在过去的时间里也发生了很大的变化--尽管在三维中解决高精度的问题仍然是一个挑战。




<a name="Abettermesh"></a><h4>A better mesh</h4>


如果你看一下上面的网格，你会发现即使域是单位盘，系数的跳动是沿着圆的，构成网格的单元也不能很好地跟踪这个几何体。原因在步骤1中已经暗示过了，在没有其他信息的情况下，Triangulation类只看到一堆粗略的网格单元，但当然不知道它们在一起看时可能代表什么样的几何形状。出于这个原因，我们需要告诉Triangulation在一个单元被细化时应该做什么：边缘中点和单元中点的新顶点应该位于哪里，以便子单元比父单元更好地代表所需的几何图形。

为了直观地了解三角计算对几何体的实际了解，仅仅输出顶点的位置和为每条边画一条直线是不够的；相反，我们必须将内部线和边界线都输出为多段线，使它们看起来是弯曲的。我们可以通过对 <code>output_results</code> 的gnuplot部分做一个改变来做到这一点。

@code
{
  GridOut       grid_out;
  std::ofstream output("grid-" + std::to_string(cycle) + ".gnuplot");
  GridOutFlags::Gnuplot gnuplot_flags(false, 5, /*curved_interior_cells*/true);
  grid_out.set_flags(gnuplot_flags);
  MappingQGeneric<dim> mapping(3);
  grid_out.write_gnuplot(triangulation, output, &mapping);
}
@endcode



在上面的代码中，我们已经对位于边界的面做了这个处理：由于我们使用了 GridGenerator::hyper_ball, ，它将一个SphericalManifold附着在域的边界上，所以这是自动发生的。为了使网格<i>interior</i>也能追踪到一个圆形域，我们需要更努力一些。首先，回顾一下我们的粗略网格由一个中心的方形单元和周围的四个单元组成。现在首先考虑一下，如果我们不仅将SphericalManifold对象连接到四个外部面，而且还连接到周边的四个单元以及它们的所有面，会发生什么。我们可以通过添加下面的片段来实现（测试一个单元的中心是否大于单元直径的一个小倍数，比如说十分之一，远离网格中心的单元直径，只对网格中心的正方形失效）。

@code
GridGenerator::hyper_ball(triangulation);
// after GridGenerator::hyper_ball is called the Triangulation has
// a SphericalManifold with id 0. We can use it again on the interior.
const Point<dim> mesh_center;
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance (cell->center()) > cell->diameter()/10)
    cell->set_all_manifold_ids(0);


triangulation.refine_global(1);
@endcode



经过几个全局细化的步骤，这将导致以下类型的网格。


  <div class="onecolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_bad_grid_4.svg" alt="一些中央单元格几乎呈三角形的网格。" width="300" height="300"> </div> </div>

这不是一个好的网格：中心单元已经被细化，位于原始中心单元<i>degenerate</i>的四个角的子单元：随着网格细化的继续，它们都倾向于三角形。这意味着从参考单元到实际单元的转换的Jacobian矩阵对这些单元来说是退化的，由于有限元解的所有误差估计都包含Jacobian矩阵的反值，你会在这些单元上得到非常大的误差，而且随着网格细化的极限，收敛顺序的损失，因为这些角落的单元在网格细化下变得越来越差。

所以我们需要更聪明的东西。为此，考虑以下最初由Konstantin Ladutenko开发的解决方案。我们将使用以下代码。

@code
GridGenerator::hyper_ball(triangulation);


const Point<dim> mesh_center;
const double core_radius  = 1.0/5.0,
             inner_radius = 1.0/3.0;


// Step 1: Shrink the inner cell
//
// We cannot get a circle out of the inner cell because of
// the degeneration problem mentioned above. Rather, shrink
// the inner cell to a core radius of 1/5 that stays
// sufficiently far away from the place where the
// coefficient will have a discontinuity and where we want
// to have cell interfaces that actually lie on a circle.
// We do this shrinking by just scaling the location of each
// of the vertices, given that the center of the circle is
// simply the origin of the coordinate system.
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) < 1e-5)
    {
      for (const auto v : cell->vertex_indices())
        cell->vertex(v) *= core_radius/mesh_center.distance(cell->vertex(v));
    }


// Step 2: Refine all cells except the central one
for (const auto &cell : triangulation.active_cell_iterators())
  if (mesh_center.distance(cell->center()) >= 1e-5)
    cell->set_refine_flag();
triangulation.execute_coarsening_and_refinement();


// Step 3: Resize the inner children of the outer cells
//
// The previous step replaced each of the four outer cells
// by its four children, but the radial distance at which we
// have intersected is not what we want to later refinement
// steps. Consequently, move the vertices that were just
// created in radial direction to a place where we need
// them.
for (const auto &cell : triangulation.active_cell_iterators())
  for (const auto v : cell->vertex_indices())
    {
      const double dist = mesh_center.distance(cell->vertex(v));
      if (dist > core_radius*1.0001 && dist < 0.9999)
        cell->vertex(v) *= inner_radius/dist;
    }


// Step 4: Apply curved manifold description
//
// As discussed above, we can not expect to subdivide the
// inner four cells (or their faces) onto concentric rings,
// but we can do so for all other cells that are located
// outside the inner radius. To this end, we loop over all
// cells and determine whether it is in this zone. If it
// isn't, then we set the manifold description of the cell
// and all of its bounding faces to the one that describes
// the spherical manifold already introduced above and that
// will be used for all further mesh refinement.
for (const auto &cell : triangulation.active_cell_iterators())
  {
    bool is_in_inner_circle = false;
    for (const auto v : cell->vertex_indices())
      if (mesh_center.distance(cell->vertex(v)) < inner_radius)
        {
          is_in_inner_circle = true;
          break;
        }


    if (is_in_inner_circle == false)
    // The Triangulation already has a SphericalManifold with
    // manifold id 0 (see the documentation of
    // GridGenerator::hyper_ball) so we just attach it to the outer
    // ring here:
      cell->set_all_manifold_ids(0);
  }
@endcode



然后，这段代码生成了以下更好的网格序列。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_0_ladutenko.svg" alt="初始网格：带有一个全局细化的Ladutenko网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_1_ladutenko.svg" alt="第一个自适应细化的Ladutenko网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_2_ladutenko.svg" alt="第二个自适应细化的Ladutenko网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_3_ladutenko.svg" alt="第三个自适应细化Ladutenko网格。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_4_ladutenko.svg" alt="第四个自适应细化Ladutenko网格。细胞沿着内圈聚集。" width="300" height="300"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_6_grid_5_ladutenko.svg" alt="第五个自适应改进的拉杜腾科网格：细胞沿着内圈聚集。" width="300" height="300"> </div> </div>

创建好的网格，特别是使它们适合你想要的几何形状，本身就是一个复杂的话题。你可以在步骤49、步骤53和步骤54中找到更多关于这个问题的内容，以及其他涉及这个问题的教程程序。步骤65展示了另一个不那么手动的方法来实现一个很适合这里的问题的网格。关于弯曲域的信息也可以在 @ref manifold "流形描述 "的文档模块中找到。

为什么选择一个跟踪内部界面的网格是有意义的？有很多原因，但最重要的原因是我们在双线性表格中实际整合的内容。从概念上讲，我们想把 $A_{ij}^K=\int_K
a(\mathbf x) \nabla \varphi_i(\mathbf x) \nabla \varphi_j(\mathbf x) ; dx$ 作为单元格 $K$ 对矩阵条目 $A_{ij}$ 的贡献来整合。我们无法精确计算，只能求助于正交法。我们知道，如果积分是平滑的，正交是准确的。这是因为正交法实质上是计算积分的多项式近似值，与积分在正交点上重合，然后计算这个多项式下的体积，作为原始积分下体积的近似值。如果积分在一个单元上是平滑的，这个多项式插值是准确的，但如果积分在一个单元上是不连续的，它通常是相当不准确的。

因此，值得将细胞排列成这样的方式，使系数不连续的界面与细胞界面对齐。这样一来，系数在每个单元上都是恒定的，之后积分将是平滑的，其多项式近似和积分的正交近似都将是准确的。请注意，这样的排列方式在许多实际案例中很常见，因此deal.II提供了一些函数（如 @ref GlossMaterialId "material_id"）来帮助管理这种情况。请参考步骤28和步骤46，了解如何应用material ids的例子。

最后，让我们考虑一个在空间中具有平滑和非均匀分布的系数的情况。我们可以再一次重复上面关于用正交法表示这种函数的所有讨论。所以，为了准确地模拟它，有几个现成的选择：你可以减少单元格的大小，增加正交公式中使用的多项式的阶数，选择一个更合适的正交公式，或进行这些步骤的组合。关键是，用正交多项式提供系数的空间依赖性的最佳拟合将导致PDE的更精确的有限元解。

最后说明一下：前面几段的讨论表明，我们在这里有一种非常具体的方式来说明我们认为的好的网格--它应该与系数中的跳跃相一致。但人们也可以在一个更普遍的环境中提出这样的问题。给定一些具有光滑解和光滑系数的方程，我们能说一个好的网格是什么样子的吗？这个问题的答案在直觉上比数学上更容易表述。一个好的网格，其单元格大体上都像正方形（或立方体，在三维空间）。一个不好的网格会包含一些在某些方向上非常细长的单元，或者，更广泛地说，其中有一些单元的边缘既短又长。有很多方法可以给每个单元分配一个数字质量指数，以衡量该单元是 "好 "还是 "坏"；其中一些经常被选择，因为它们便宜且容易计算，而另一些则是基于收敛性证明中的内容。前者的一个例子是一个单元格的最长边与最短边的比率。在理想的情况下，这个比率是1；不好的单元格的值远远大于1。后者的例子是考虑从参考单元 $\hat K=[0,1]^d$ 到实际单元 $K$ 的映射的梯度（"Jacobian"）；这个梯度是一个矩阵，进入误差估计的一个量是参考单元上所有点的最大值，即这个矩阵的最大和最小的特征值的比率。同样不难看出，如果单元格 $K$ 是 $\hat K$ 的仿生图像，这个比率是恒定的，对于正方形和长方体来说是一个。

在实践中，将这种质量衡量标准可视化可能是很有趣的。函数 GridTools::compute_aspect_ratio_of_cells() 提供了一种获得这种信息的方法。更好的是，可视化工具，如VisIt，通常允许你在可视化软件中对各种措施进行可视化；在VisIt的情况下，只需添加一个 "伪色 "图，并选择一个网格质量措施，而不是解决方案领域。




<a name="Playingwiththeregularityofthesolution"></a><h4>Playing with the regularity of the solution</h4>


从数学的角度来看，拉普拉斯方程的解决方案

@f[


  -\Delta u = f


@f]

在光滑有界的凸域上，已知其本身是光滑的。确切的光滑程度，即解所处的函数空间，取决于域的边界到底有多光滑，以及右手边的光滑程度。边界处可能会失去解的某些规律性，但一般来说，在域的紧凑子集中，解的可微性比右手边多一倍。特别是如果右手边满足 $f\in C^\infty(\Omega)$ ，那么 $u \in C^\infty(\Omega_i)$ ，其中 $\Omega_i$ 是 $\Omega$ 的任何紧凑子集（ $\Omega$ 是一个开放域，所以紧凑子集需要与 $\partial\Omega$ 保持一个正距离）。

然而，我们为目前的例子选择的情况是不同的：我们看的是一个具有非常数系数的方程  $a(\mathbf x)$  。

@f[


  -\nabla \cdot (a \nabla u) = f.


@f]

这里，如果 $a$ 不光滑，那么无论 $f$ 如何，解也不会光滑。特别是，我们希望在 $a$ 沿直线（或沿三维平面）不连续的地方，解会有一个结点。这很容易看出来：例如，如果 $f$ 是连续的，那么 $f=-\nabla \cdot (a \nabla u)$ 也需要是连续的。这意味着 $a \nabla u$ 必须是可连续微分的（不存在扭结）。因此，如果 $a$ 有一个不连续，那么 $\nabla u$ 必须有一个相反的不连续，从而使两者完全抵消，它们的乘积得到一个没有不连续的函数。但是要使 $\nabla u$ 有一个不连续， $u$ 必须有一个结点。当然，这正是当前例子中所发生的情况，在解的图片中也很容易观察到。

一般来说，如果系数 $a(\mathbf x)$ 沿着2D的直线或3D的平面是不连续的，那么解可能有一个结点，但解的梯度不会到无限大。这意味着，解至少还在<a href="https://en.wikipedia.org/wiki/Sobolev_space">Sobolev space</a> $W^{1,\infty}$ 中（也就是说，大致上是在导数有界的函数空间中）。另一方面，我们知道，在最极端的情况下--即域有重入角，右手边只满足 $f\in H^{-1}$ ，或者系数 $a$ 只在 $L^\infty$ 中--我们所能期望的是， $u\in H^1$ （即导数是可平方整除的函数的<a
href="https://en.wikipedia.org/wiki/Sobolev_space#Sobolev_spaces_with_integer_k">Sobolev
space</a>），是比 $W^{1,\infty}$ 大很多的空间 。要创造出解在空间 $H^{1+s}$ 中的案例并不十分困难，我们可以让 $s$ 变得像我们想要的那样小。这样的情况经常被用来测试自适应有限元方法，因为网格要解决导致解不再在 $W^{1,\infty}$ 中的奇异点。

人们为此使用的典型例子叫做<i>Kellogg problem</i>（指 @cite Kel74 ），在常用的形式中，它的系数 $a(\mathbf x)$ 在平面的四个象限有不同的值（或在 ${\mathbb R}^3$ 的八个象限有不同的值）。确切的规则性程度（上述索博列夫空间索引中的 $s$ ）取决于 $a(\mathbf x)$ 的值在原点处聚集，通过选择足够大的跳跃，可以使解的规则性尽可能地接近 $H^1$  。

为了实现这样的东西，可以用以下方法来代替系数函数（这里只显示2d情况）。

@code
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

(在结尾处添加 <code>Assert</code> ，以确保在我们到达那个点时，要么抛出一个异常，要么程序中止。

--当然我们不应该这样做，但这是给自己上保险的好方法：我们都会犯错，因为有时没有想到所有的情况，例如检查 <code>p[0]</code> 是否小于和大于零，而不是大于或等于零，从而忘记了一些情况，否则会导致难以发现的错误。最后的 <code>return 0;</code> 只是为了避免编译器警告说函数没有在 <code>return</code> 语句中结束 -- 编译器无法看到由于前面的 <code>Assert</code> 语句，函数实际上永远不会到达那个点）。)

通过玩弄这种四个或更多的扇形聚集在一起，并且在这些扇形上的系数有不同的值的情况，我们可以构造出解在原点有奇异点的情况。我们还可以看到在这种情况下网格是如何被细化的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-6.cc"
*/
