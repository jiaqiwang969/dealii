/**
@page step_2 The step-2 tutorial program
This tutorial depends on step-1.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Sparsity"> Sparsity </a>
        <li><a href="#Howdegreesoffreedomareenumerated"> How degrees of freedom are enumerated </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Meshgeneration">Mesh generation</a>
        <li><a href="#CreationofaDoFHandler">Creation of a DoFHandler</a>
        <li><a href="#RenumberingofDoFs">Renumbering of DoFs</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-2/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{9} 

在前面的例子中，我们已经创建了一个网格，现在我们展示如何在这个网格上定义自由度。在这个例子中，我们将使用最低阶（ $Q_1$ ）的有限元，自由度与网格的顶点相关联。以后的例子将展示更高阶的元素，自由度不一定与顶点相关，但可以与边、面或单元相关。

术语 "自由度 "在有限元界通常用来表示两个略有不同但相关的事情。首先是我们希望将有限元解表示为形状函数的线性组合，形式为 $u_h(\mathbf x) = \sum_{j=0}^{N-1} U_j \varphi_j(\mathbf
x)$  。这里， $U_j$ 是一个膨胀系数的向量。因为我们还不知道它们的值（我们将计算它们作为线性或非线性系统的解），它们被称为 "未知数 "或 "自由度"。该术语的第二个含义可以解释如下。对有限元问题的数学描述通常是说，我们正在寻找一个满足某些方程组的有限维函数 $u_h \in V_h$ （例如， $a(u_h,\varphi_h)=(f,\varphi_h)$ 为所有测试函数 $\varphi_h\in
V_h$ ）。换句话说，我们在这里说的是，解决方案需要位于某个空间  $V_h$  中。然而，为了在计算机上实际解决这个问题，我们需要选择这个空间的一个基；这就是我们在上面用系数 $U_j$ 对 $u_h(\mathbf x)$ 进行展开的形状函数 $\varphi_j(\mathbf x)$ 的集合。当然，空间 $V_h$ 的基数有很多，但我们将特别选择由传统上在网格单元上局部定义的有限元函数描述的基数。在这种情况下描述 "自由度 "需要我们简单地 <i>enumerate</i> 空间的基函数  $V_h$  。对于 $Q_1$ 元素，这意味着简单地以某种方式列举网格的顶点，但对于高阶元素，还必须列举与网格的边、面或单元内部相关的形状函数。换句话说，自由度的枚举是完全独立于我们用于顶点的索引的。提供这种列举 $V_h$ 的基础函数的类被称为DoFHandler。

在网格上定义自由度（简称 "DoF"）是一个相当简单的任务，因为这个库为你做了所有的工作。基本上，你所要做的就是创建一个有限元对象（从deal.II已有的众多有限元类中选取，例如参见 @ref fe 文档），并通过 DoFHandler::distribute_dofs 函数将其交给DoFHandler对象（"分配DoF "是我们用来描述上文讨论的<i>enumerating</i>基函数过程的术语）。DoFHandler是一个知道哪些自由度住在哪里的类，也就是说，它可以回答 "全局有多少自由度 "和 "在这个单元上，给我住在这里的形状函数的全局索引 "这样的问题。当你决定你的系统矩阵应该有多大时，以及当把单个单元的贡献复制到全局矩阵时，你需要这种信息。

<a name="Sparsity"></a><h3> Sparsity </h3>


然后，下一步将是利用这个有限元和网格计算与特定微分方程对应的矩阵和右手。我们将为第三步程序保留这一步骤，而是谈论有限元程序的一个实际问题，即有限元矩阵总是非常稀疏的：这些矩阵中的几乎所有条目都是零。

更准确地说，如果一个矩阵中的非零项<i>per row</i>的数量与整个自由度的数量无关，我们就说该矩阵是稀疏的。例如，拉普拉斯方程的有限差分近似的简单5点模版导致了一个稀疏矩阵，因为每行的非零条目数是5，因此与矩阵的总大小无关。对于更复杂的问题--比如说步骤22的斯托克斯问题--特别是在三维中，每行的条目数可能是几百个。但重要的一点是，这个数字与问题的总体大小无关：如果你细化网格，每行未知数的最大数量保持不变。

与使用泰勒扩展和匹配系数来逼近偏微分方程的解，或使用傅里叶基相比，稀疏性是有限元方法的一个突出特点。

在实践中，正是由于矩阵的稀疏性，使我们能够解决有数百万或数十亿未知数的问题。为了理解这一点，请注意，一个有 $N$ 行的矩阵，每个非零项的数量都有固定的上限，需要 ${\cal O}(N)$ 个内存位置来存储，而矩阵-向量乘法也只需要 ${\cal O}(N)$ 次操作。因此，如果我们有一个线性求解器，只需要固定数量的矩阵向量乘法就能得出这个矩阵的线性系统的解，那么我们就会有一个能以最佳复杂度找到所有 $N$ 未知数的值的求解器，也就是说，总共只需要 ${\cal O}(N)$ 次操作。很明显，如果矩阵不是稀疏的，这是不可能的（因为那样的话，矩阵中的条目数必须是 ${\cal O}(N^s)$ 与一些 $s>1$ ，做固定数量的矩阵-向量乘积将需要 ${\cal O}(N^s)$ 次操作），但这也需要非常专业的求解器，如多网格方法，以满足求解只需要固定数量的矩阵-向量乘法的要求。我们将在本教程的剩余程序中经常研究使用什么求解器的问题。

稀疏性是由以下事实产生的：有限元形状函数是在单个单元上定义的<i>locally</i>，而不是全局的，并且双线性形式中的局部微分算子只对支持度重叠的形状函数进行耦合。一个函数的 "支持 "是指它的非零区域。对于有限元方法，形状函数的支持通常是指与它所定义的顶点、边或面相邻的单元。)换句话说，自由度 $i$ 和 $j$ 如果不是定义在同一个单元上，就不会重叠，因此，矩阵条目 $A_{ij}$ 将为零。  (在某些情况下，如非连续加尔金法，形状函数也可以通过面积分连接到相邻的单元。但是有限元方法一般不会将形状函数与定义了该函数的单元的近邻相联系）。)




<a name="Howdegreesoffreedomareenumerated"></a><h3> How degrees of freedom are enumerated </h3>


默认情况下，DoFHandler类以一种相当随机的方式枚举网格上的自由度；因此，稀疏度模式也没有为任何特定的目的进行优化。为了说明这一点，下面的代码将演示一个简单的方法来输出对应于DoFHandler的 "稀疏模式"，即一个对象代表了在网格上离散偏微分方程时可能建立的矩阵的所有潜在非零元素及其DoFHandler。这种缺乏结构的疏散模式将从我们下面展示的图片中显现出来。

对于大多数应用和算法来说，自由度的确切编号方式并不重要。例如，我们用来解决线性系统的共轭梯度方法并不关心。另一方面，有些算法确实关心：特别是一些预处理程序，如SSOR，如果它们能以特定的顺序走过自由度，就能更好地工作，如果我们能以这样的方式排序，使SSOR能以这样的顺序从零到 $N$ 迭代它们，那就太好了。其他的例子包括计算不完整的LU或Cholesky分解，或者如果我们关心矩阵的块状结构（见步骤20的例子）。因此，deal.II在命名空间DoFRenumbering中有一些算法可以以特定的方式重新列举自由度。重新编号可以被认为是选择了一个不同的、经排列的有限元空间的基础。因此，这种重新编号所产生的稀疏模式和矩阵与我们没有明确的重新编号所得到的相比，也只是行和列的排列组合。

在下面的程序中，我们将使用Cuthill和McKee的算法来完成。我们将在<a href="#Results">results section</a>中展示原始自由度列举和下面重新编号的版本的稀疏模式。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The first few includes are just like in the previous program, so do not
 * require additional comments:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * @endcode
 * 
 * However, the next file is new. We need this include file for the
 * association of degrees of freedom ("DoF"s) to vertices, lines, and cells:
 * 
 * @code
 * #include <deal.II/dofs/dof_handler.h>
 * 
 * @endcode
 * 
 * The following include contains the description of the bilinear finite
 * element, including the facts that it has one degree of freedom on each
 * vertex of the triangulation, but none on faces and none in the interior of
 * the cells.
 * 

 * 
 * (In fact, the file contains the description of Lagrange elements in
 * general, i.e. also the quadratic, cubic, etc versions, and not only for 2d
 * but also 1d and 3d.)
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * @endcode
 * 
 * In the following file, several tools for manipulating degrees of freedom
 * can be found:
 * 
 * @code
 * #include <deal.II/dofs/dof_tools.h>
 * @endcode
 * 
 * We will use a sparse matrix to visualize the pattern of nonzero entries
 * resulting from the distribution of degrees of freedom on the grid. That
 * class can be found here:
 * 
 * @code
 * #include <deal.II/lac/sparse_matrix.h>
 * @endcode
 * 
 * We will also need to use an intermediate sparsity pattern structure, which
 * is found in this file:
 * 
 * @code
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * @endcode
 * 
 * We will want to use a special algorithm to renumber degrees of freedom. It
 * is declared here:
 * 
 * @code
 * #include <deal.II/dofs/dof_renumbering.h>
 * 
 * @endcode
 * 
 * And this is again needed for C++ output:
 * 
 * @code
 * #include <fstream>
 * 
 * @endcode
 * 
 * Finally, as in step-1, we import the deal.II namespace into the global
 * scope:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Meshgeneration"></a> 
 * <h3>Mesh generation</h3>
 * 

 * 
 * This is the function that produced the circular grid in the previous step-1
 * example program with fewer refinements steps. The sole difference is that it
 * returns the grid it produces via its argument.
 * 
 * @code
 * void make_grid(Triangulation<2> &triangulation)
 * {
 *   const Point<2> center(1, 0);
 *   const double   inner_radius = 0.5, outer_radius = 1.0;
 *   GridGenerator::hyper_shell(
 *     triangulation, center, inner_radius, outer_radius, 5);
 * 
 *   for (unsigned int step = 0; step < 3; ++step)
 *     {
 *       for (auto &cell : triangulation.active_cell_iterators())
 *         for (const auto v : cell->vertex_indices())
 *           {
 *             const double distance_from_center =
 *               center.distance(cell->vertex(v));
 * 
 *             if (std::fabs(distance_from_center - inner_radius) <=
 *                 1e-6 * inner_radius)
 *               {
 *                 cell->set_refine_flag();
 *                 break;
 *               }
 *           }
 * 
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="CreationofaDoFHandler"></a> 
 * <h3>Creation of a DoFHandler</h3>
 * 

 * 
 * Up to now, we only have a grid, i.e. some geometrical (the position of the
 * vertices) and some topological information (how vertices are connected to
 * lines, and lines to cells, as well as which cells neighbor which other
 * cells). To use numerical algorithms, one needs some logic information in
 * addition to that: we would like to associate degree of freedom numbers to
 * each vertex (or line, or cell, in case we were using higher order elements)
 * to later generate matrices and vectors which describe a finite element
 * field on the triangulation.
 * 

 * 
 * This function shows how to do this. The object to consider is the
 * <code>DoFHandler</code> class template.  Before we do so, however, we first
 * need something that describes how many degrees of freedom are to be
 * associated to each of these objects. Since this is one aspect of the
 * definition of a finite element space, the finite element base class stores
 * this information. In the present context, we therefore create an object of
 * the derived class <code>FE_Q</code> that describes Lagrange elements. Its
 * constructor takes one argument that states the polynomial degree of the
 * element, which here is one (indicating a bi-linear element); this then
 * corresponds to one degree of freedom for each vertex, while there are none
 * on lines and inside the quadrilateral. A value of, say, three given to the
 * constructor would instead give us a bi-cubic element with one degree of
 * freedom per vertex, two per line, and four inside the cell. In general,
 * <code>FE_Q</code> denotes the family of continuous elements with complete
 * polynomials (i.e. tensor-product polynomials) up to the specified order.
 * 

 * 
 * We first need to create an object of this class and then pass it on to the
 * <code>DoFHandler</code> object to allocate storage for the degrees of
 * freedom (in deal.II lingo: we <i>distribute degrees of
 * freedom</i>).
 * 
 * @code
 * void distribute_dofs(DoFHandler<2> &dof_handler)
 * {
 *   const FE_Q<2> finite_element(1);
 *   dof_handler.distribute_dofs(finite_element);
 * 
 * @endcode
 * 
 * Now that we have associated a degree of freedom with a global number to
 * each vertex, we wonder how to visualize this?  There is no simple way to
 * directly visualize the DoF number associated with each vertex. However,
 * such information would hardly ever be truly important, since the
 * numbering itself is more or less arbitrary. There are more important
 * factors, of which we will demonstrate one in the following.
 *   

 * 
 * Associated with each vertex of the triangulation is a shape
 * function. Assume we want to solve something like Laplace's equation, then
 * the different matrix entries will be the integrals over the gradient of
 * each pair of such shape functions. Obviously, since the shape functions
 * are nonzero only on the cells adjacent to the vertex they are associated
 * with, matrix entries will be nonzero only if the supports of the shape
 * functions associated to that column and row %numbers intersect. This is
 * only the case for adjacent shape functions, and therefore only for
 * adjacent vertices. Now, since the vertices are numbered more or less
 * randomly by the above function (DoFHandler::distribute_dofs), the pattern
 * of nonzero entries in the matrix will be somewhat ragged, and we will
 * take a look at it now.
 *   

 * 
 * First we have to create a structure which we use to store the places of
 * nonzero elements. This can then later be used by one or more sparse
 * matrix objects that store the values of the entries in the locations
 * stored by this sparsity pattern. The class that stores the locations is
 * the SparsityPattern class. As it turns out, however, this class has some
 * drawbacks when we try to fill it right away: its data structures are set
 * up in such a way that we need to have an estimate for the maximal number
 * of entries we may wish to have in each row. In two space dimensions,
 * reasonable values for this estimate are available through the
 * DoFHandler::max_couplings_between_dofs() function, but in three
 * dimensions the function almost always severely overestimates the true
 * number, leading to a lot of wasted memory, sometimes too much for the
 * machine used, even if the unused memory can be released immediately after
 * computing the sparsity pattern. In order to avoid this, we use an
 * intermediate object of type DynamicSparsityPattern that uses a
 * different %internal data structure and that we can later copy into the
 * SparsityPattern object without much overhead. (Some more information on
 * these data structures can be found in the @ref Sparsity module.) In order
 * to initialize this intermediate data structure, we have to give it the
 * size of the matrix, which in our case will be square with as many rows
 * and columns as there are degrees of freedom on the grid:
 * 
 * @code
 *   DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
 *                                                   dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * We then fill this object with the places where nonzero elements will be
 * located given the present numbering of degrees of freedom:
 * 
 * @code
 *   DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
 * 
 * @endcode
 * 
 * Now we are ready to create the actual sparsity pattern that we could
 * later use for our matrix. It will just contain the data already assembled
 * in the DynamicSparsityPattern.
 * 
 * @code
 *   SparsityPattern sparsity_pattern;
 *   sparsity_pattern.copy_from(dynamic_sparsity_pattern);
 * 
 * @endcode
 * 
 * With this, we can now write the results to a file:
 * 
 * @code
 *   std::ofstream out("sparsity_pattern1.svg");
 *   sparsity_pattern.print_svg(out);
 * @endcode
 * 
 * The result is stored in an <code>.svg</code> file, where each nonzero entry
 * in the matrix corresponds with a red square in the image. The output will
 * be shown below.
 *   

 * 
 * If you look at it, you will note that the sparsity pattern is
 * symmetric. This should not come as a surprise, since we have not given
 * the <code>DoFTools::make_sparsity_pattern</code> any information that
 * would indicate that our bilinear form may couple shape functions in a
 * non-symmetric way. You will also note that it has several distinct
 * region, which stem from the fact that the numbering starts from the
 * coarsest cells and moves on to the finer ones; since they are all
 * distributed symmetrically around the origin, this shows up again in the
 * sparsity pattern.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="RenumberingofDoFs"></a> 
 * <h3>Renumbering of DoFs</h3>
 * 

 * 
 * In the sparsity pattern produced above, the nonzero entries extended quite
 * far off from the diagonal. For some algorithms, for example for incomplete
 * LU decompositions or Gauss-Seidel preconditioners, this is unfavorable, and
 * we will show a simple way how to improve this situation.
 * 

 * 
 * Remember that for an entry $(i,j)$ in the matrix to be nonzero, the
 * supports of the shape functions i and j needed to intersect (otherwise in
 * the integral, the integrand would be zero everywhere since either the one
 * or the other shape function is zero at some point). However, the supports
 * of shape functions intersected only if they were adjacent to each other, so
 * in order to have the nonzero entries clustered around the diagonal (where
 * $i$ equals $j$), we would like to have adjacent shape functions to be
 * numbered with indices (DoF numbers) that differ not too much.
 * 

 * 
 * This can be accomplished by a simple front marching algorithm, where one
 * starts at a given vertex and gives it the index zero. Then, its neighbors
 * are numbered successively, making their indices close to the original
 * one. Then, their neighbors, if not yet numbered, are numbered, and so on.
 * 

 * 
 * One algorithm that adds a little bit of sophistication along these lines is
 * the one by Cuthill and McKee. We will use it in the following function to
 * renumber the degrees of freedom such that the resulting sparsity pattern is
 * more localized around the diagonal. The only interesting part of the
 * function is the first call to <code>DoFRenumbering::Cuthill_McKee</code>,
 * the rest is essentially as before:
 * 
 * @code
 * void renumber_dofs(DoFHandler<2> &dof_handler)
 * {
 *   DoFRenumbering::Cuthill_McKee(dof_handler);
 * 
 *   DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
 *                                                   dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
 * 
 *   SparsityPattern sparsity_pattern;
 *   sparsity_pattern.copy_from(dynamic_sparsity_pattern);
 * 
 *   std::ofstream out("sparsity_pattern2.svg");
 *   sparsity_pattern.print_svg(out);
 * }
 * 
 * @endcode
 * 
 * Again, the output is shown below. Note that the nonzero entries are
 * clustered far better around the diagonal than before. This effect is even
 * more distinguished for larger matrices (the present one has 1260 rows and
 * columns, but large matrices often have several 100,000s).
 * 

 * 
 * It is worth noting that the <code>DoFRenumbering</code> class offers a
 * number of other algorithms as well to renumber degrees of freedom. For
 * example, it would of course be ideal if all couplings were in the lower or
 * upper triangular part of a matrix, since then solving the linear system
 * would amount to only forward or backward substitution. This is of course
 * unachievable for symmetric sparsity patterns, but in some special
 * situations involving transport equations, this is possible by enumerating
 * degrees of freedom from the inflow boundary along streamlines to the
 * outflow boundary. Not surprisingly, <code>DoFRenumbering</code> also has
 * algorithms for this.
 * 

 * 
 * 

 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * Finally, this is the main program. The only thing it does is to allocate
 * and create the triangulation, then create a <code>DoFHandler</code> object
 * and associate it to the triangulation, and finally call above two functions
 * on it:
 * 
 * @code
 * int main()
 * {
 *   Triangulation<2> triangulation;
 *   make_grid(triangulation);
 * 
 *   DoFHandler<2> dof_handler(triangulation);
 * 
 *   distribute_dofs(dof_handler);
 *   renumber_dofs(dof_handler);
 * }
 * @endcode
examples/step-2/doc/results.dox



<a name="Results"></a><h1>Results</h1>


该程序运行后，产生了两个稀疏模式。我们可以通过在网络浏览器中打开 <code>.svg</code> 文件来可视化它们。

结果是这样的（每一个点都表示一个可能为非零的条目；当然，这个条目是否真的为零取决于所考虑的方程，但矩阵中的指示位置告诉我们，在离散化局部，即微分方程时，哪些形状函数可以，哪些不可以耦合）。   <table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-1.svg" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-2.sparsity-2.svg" alt=""></td>
  </tr>
</table> 

左图中的不同区域，由线条中的扭结和左边和上面的单点表示，代表了三角法不同细化层次上的自由度。  从右图中可以看出，重新编号后，稀疏模式在矩阵的主对角线附近的聚类情况要好得多。虽然这可能不明显，但两张图片中非零项的数量当然是一样的。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


就像第1步一样，你可能想在程序中玩一下，熟悉一下deal.II。例如，在 <code>distribute_dofs</code> 函数中，我们使用线性有限元（FE_Q对象的参数 "1 "就是如此）。探索一下如果你使用高阶元素，例如立方或五元元素（使用3和5作为各自的参数），稀疏模式会有什么变化。

你也可以通过细化网格来探索稀疏性模式的变化。你会发现，不仅矩阵的大小会发生变化，其带宽（矩阵中离对角线最远的那些非零元素与对角线的距离）也会发生变化，不过带宽与大小的比例通常会缩小，也就是说，矩阵在对角线周围聚集得更多。

实验的另一个想法是尝试DoFRenumbering命名空间中除Cuthill-McKee之外的其他重新编号策略，看看它们如何影响稀疏性模式。

你也可以使用<a
href="http://www.gnuplot.info/">GNUPLOT</a>（较简单的可视化程序之一；也许不是最容易使用的，因为它是命令行驱动的，但在所有Linux和其他类似Unix的系统上也是普遍可用的）通过改变 <code>print_svg()</code> to <code>print_gnuplot()</code> in <code>distribute_dofs()</code> and <code>renumber_dofs()</code> 来使输出可视化。

@code
examples/\step-2> gnuplot


        G N U P L O T
        Version 3.7 patchlevel 3
        last modified Thu Dec 12 13:00:00 GMT 2002
        System: Linux 2.6.11.4-21.10-default


        Copyright(C) 1986 - 1993, 1998 - 2002
        Thomas Williams, Colin Kelley and many others


        Type `help` to access the on-line reference manual
        The gnuplot FAQ is available from
        http://www.gnuplot.info/gnuplot-faq.html


        Send comments and requests for help to <info-gnuplot@dartmouth.edu>
        Send bugs, suggestions and mods to <bug-gnuplot@dartmouth.edu>



Terminal type set to 'x11'
gnuplot> set style data points
gnuplot> plot "sparsity_pattern.1"
@endcode



另一个基于<a href="http://www.gnuplot.info/">GNUPLOT</a>的做法是尝试打印出带有支撑点位置和编号的网格。为此，你需要包含GridOut和MappingQ1的头文件。这方面的代码是。

@code
  std::ofstream out("gnuplot.gpl");
  out << "plot '-' using 1:2 with lines, "
      << "'-' with labels point pt 2 offset 1,1"
      << std::endl;
  GridOut().write_gnuplot (triangulation, out);
  out << "e" << std::endl;
  const int dim = 2;
  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points (MappingQ1<dim>(),
                                        dof_handler,
                                        support_points);
  DoFTools::write_gnuplot_dof_support_point_info(out,
                                                 support_points);
  out << "e" << std::endl;
@endcode

在我们运行该代码后，我们得到了一个名为gnuplot.gpl的文件。要查看这个文件，我们可以在命令行中运行以下代码。

@code
gnuplot -p gnuplot.gpl
@endcode.有了这个，你会得到一个类似于 @image html support_point_dofs1.png 的图片，这取决于你正在看的网格。更多信息，见 DoFTools::write_gnuplot_dof_support_point_info. 。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-2.cc"
*/
