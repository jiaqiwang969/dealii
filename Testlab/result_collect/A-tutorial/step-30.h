//include/deal.II-translator/A-tutorial/step-30_0.txt
/**
  @page step_30 The step-30 tutorial program 
* 本教程依赖于  step-12  。
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Overview">Overview</a><a href="#Overview">Overview</a>
        <li><a href="#Anisotropicrefinement">Anisotropic refinement</a><a href="#Anisotropicrefinement">Anisotropic refinement</a>
      <ul>
        <li><a href="#Motivation">Motivation</a><a href="#Motivation">Motivation</a>
      </ul>
        <li><a href="#Implementation">Implementation</a><a href="#Implementation">Implementation</a>
      <ul>
        <li><a href="#Meshsmoothing">Mesh smoothing</a><a href="#Meshsmoothing">Mesh smoothing</a>
      </ul>
        <li><a href="#Jumpindicator">Jump indicator</a><a href="#Jumpindicator">Jump indicator</a>
        <li><a href="#Theproblem">The problem</a><a href="#Theproblem">The problem</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a><a href="#Equationdata">Equation data</a>
        <li><a href="#ClassDGTransportEquation">Class: DGTransportEquation</a><a href="#ClassDGTransportEquation">Class: DGTransportEquation</a>
        <li><a href="#ClassDGMethod">Class: DGMethod</a><a href="#ClassDGMethod">Class: DGMethod</a>
      <ul>
        <li><a href="#Functionassemble_system">Function: assemble_system</a><a href="#Functionassemble_system">Function: assemble_system</a>
      </ul>
        <li><a href="#Solver">Solver</a><a href="#Solver">Solver</a>
        <li><a href="#Refinement">Refinement</a><a href="#Refinement">Refinement</a>
        <li><a href="#TheRest">The Rest</a><a href="#TheRest">The Rest</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
*<a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1> 。



*<a name="Overview"></a><h3>Overview</h3>


* 这个例子专门讨论 <em> 各向异性细化 </em> ，它扩展了局部细化的可能性。在大多数情况下，这是对 step-12 教程程序的修改，我们使用相同的DG方法来处理线性转场方程。这个程序将涵盖以下主题。 <ol>   <li>   <em>  各向异性的细化  </em>  ：各向异性的细化是什么意思？    <li>   <em>  实现  </em>  ：对代码进行必要的修改，以便与各向异性的细化网格一起工作。    <li>   <em>  跳跃指标  </em>  : 在DG方法的背景下，各向异性细化的简单指标。 </ol>  将不讨论离散化本身，也不讨论此处使用的非各向异性细化的实现技术。请参考 step-12 的内容。
* 请注意，在编写这个教程程序的时候，各向异性细化只在不连续的Galerkin有限元中完全实现。这一点以后可能会改变（或者已经改变）。
*


*  @note  虽然这个程序是对 step-12 的修改，但它是对 step-12 版本的改编，是在deal.II历史上早期写的，当时MeshWorker框架还没有出现。因此，它与现在的 step-12 没有什么相似之处，除了它用同样的离散化方法解决了同样的方程。
*


*<a name="Anisotropicrefinement"></a><h3>Anisotropic refinement</h3>


* 在前面的教程程序中，所有的适应过程都是基于 <em> 各向同性 </em> 细化单元，它将所有的边切成两半，并在这些分割的边上形成新的单元（当然，还包括一些额外的边、面和顶点）。在交易二中， <em> 各向异性细化 </em> 指的是只拆分部分边而不改变其他边的过程。以简单的方形单元为例，考虑一下。
* @code
* 


* 
* 
-------*
  |       |
  |       |
  |       |
* 


* 
* 
-------*
@endcode
* 经过通常的细化，它将由四个子单元组成，看起来像这样。
* @code
* 


* 
* 
---*---*
  |   |   |
* 


* 
* 
---*---*     RefinementCase<2>::cut_xy
  |   |   |
* 


* 
* 
---*---*
@endcode
* 新的各向异性细化可以采取两种形式：要么我们可以将平行于水平X轴的边缘分割开来，形成这两个子单元。
* @code
* 


* 
* 
---*---*
  |   |   |
  |   |   |     RefinementCase<2>::cut_x
  |   |   |
* 


* 
* 
---*---*
@endcode
* 或者我们可以分割沿Y轴运行的两条边，再次产生两个子单元，不过看起来是这样的。
* @code
* 


* 
* 
-------*
  |       |
* 


* 
* 
-------*     RefinementCase<2>::cut_y
  |       |
* 


* 
* 
-------*
@endcode
* 所有单元格的细化情况都由一个枚举 RefinementPossibilities::Possibilities, 来描述，上述各向异性情况被称为 @p cut_x 和 @p cut_y ，原因很明显。各向异性的细化情况在二维中被称为 @p cut_xy ，可以通过 RefinementCase<dim>::isotropic_refinement. 从RefinementCase类中请求。
* 在三维中，有第三个轴可以被分割，即Z轴，因此我们在这里有一个额外的细化案例 @p cut_z 。各向同性的细化现在将沿x轴、y轴和z轴细化一个单元，因此被称为 @p  cut_xyz。另外还有 @p cut_xy,  @p cut_xz 和 @p cut_yz 的情况，它们沿两个轴精炼单元，但不沿第三个轴精炼。给出一个x轴向右、y轴 "进入页面"、z轴在顶部的六面体单元。
* @code
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
-----------*
     /           /|
    /           / |
   /           /  |
* 


* 
* 
-----------*   |
  |           |   |
  |           |  
  |           |  /
  |           | /
  |           |/
* 


* 
* 
-----------*
@endcode
*我们有各向同性的细化情况。
* @code
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
-----*-----*
     /     /     /|
* 


* 
* 


* 
* 


* 
* 
-----*-----* |
   /     /     /|
* 


* 
* 
-----*-----* |/|
  |     |     | |
  |     |     |/|
* 


* 
* 
-----*-----* |/
  |     |     |
  |     |     |/
* 


* 
* 
-----*-----*
* 
  RefinementCase<3>::cut_xyz
@endcode
*三种各向异性的情况，只对一个轴进行细化。
* @code
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
-----*-----*
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
-----------*
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
-----------*
     /     /     /|            /           /|            /           /|
    /     /     / |
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
* 


* 
* 


* 
* 


* 
* 
-----------* |           /           / |
   /     /     /  |          /           /| |          /           /
* 


* 
* 


* 
* 
-----*-----*   |
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
* 


* 
* 
-----------* | |
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
* 


* 
* 
-----------*  /|
  |     |     |   |         |           | | |         |           | / |
  |     |     |           |           | |         |           |/ 
  |     |     |  /          |           | |/
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
* 


* 
* 


* 
* 
-----------*  /
  |     |     | /           |           |           |           | /
  |     |     |/            |           |/            |           |/
* 


* 
* 
-----*-----*
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
-----------*
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
-----------*
* 
  RefinementCase<3>::cut_x  RefinementCase<3>::cut_y  RefinementCase<3>::cut_z
@endcode
*还有三个细化三个轴中的两个轴的情况。
* @code
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
-----*-----*
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
-----*-----*
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
-----------*
     /     /     /|            /     /     /|            /           /|
* 


* 
* 


* 
* 


* 
* 
-----*-----* |           /     /     / |
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
* 


* 
* 


* 
* 


* 
* 
-----------* |
   /     /     /| |          /     /     /           /           /|
* 


* 
* 
-----*-----* | |
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
* 


* 
* 
-----*-----*  /|
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
* 


* 
* 
-----------* |/|
  |     |     | | |         |     |     | / |         |           | |
  |     |     | |         |     |     |/          |           |/|
  |     |     | |/
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
* 


* 
* 


* 
* 
-----*-----*  /
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
* 


* 
* 


* 
* 
-----------* |/
  |     |     |           |     |     | /           |           |
  |     |     |/            |     |     |/            |           |/
* 


* 
* 
-----*-----*
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
-----*-----*
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
-----------*
* 
  RefinementCase<3>::cut_xy RefinementCase<3>::cut_xz RefinementCase<3>::cut_yz
@endcode
* 对于一维问题，各向异性的细化不会产生任何影响，因为一个单元只有一个坐标方向，所以除了各向同性外，不可能以任何其他方式分割。
*<a name="Motivation"></a><h4>Motivation</h4>
* 自适应局部细化被用来获得精细的网格，这些网格能够很好地适应有效解决手头的问题。简而言之，产生较大误差的单元的尺寸被缩小，以获得对手头问题的更好的近似解。然而，很多问题都含有各向异性的特征。突出的例子是可压缩粘性流动中的冲击或边界层。一个有效的网格可以用较高长宽比的单元来逼近这些特征，而这些单元是根据上述特征来确定方向的。只使用各向同性的细化，原始网格单元的长宽比会被保留，因为它们会被单元的子代所继承。因此，从各向同性的网格开始，边界层将被细化，以捕捉壁面法线方向上流场的快速变化，从而导致在法线和切线方向上都具有非常小的边长的单元。通常情况下，切线方向的边长要大得多，因此可以使用更少的单元，而不会在近似精度上有明显的损失。各向异性的细化过程可以通过每个细化步骤将母细胞和子细胞的长宽比修改为2的系数。在多次细化的过程中，精细单元的长宽比可以被优化，节省了相当数量的单元和相应的自由度，从而节省了计算资源、内存和CPU时间。
*<a name="Implementation"></a><h3>Implementation</h3>
*

* 大多数时候，当我们进行有限元计算时，我们一次只考虑一个单元，例如计算单元对全局矩阵的贡献，或插值边界值。然而，有时我们不得不关注单元在我们的算法中是如何关联的。单元之间的关系有两种形式：邻居关系和母子关系。对于各向异性细化的情况，deal.II对始终保持的单元格关系使用了某些约定（不变量）。例如，一个细化的单元总是正好有 $2^{dim}$ 个孩子。而且（除了1d情况），两个相邻的单元格最多可以相差一个细化级别：它们同样经常被细化，或者其中一个正好再细化一次，在共同面上正好留下一个悬空的节点。几乎所有的时候，这些不变量只在库的内部实现中才会被关注。然而，在有些情况下，对它们的了解也与应用程序有关。
* 在目前的情况下，值得注意的是，网格细化的种类影响到一些最基本的假设。因此，在应用程序中发现的一些常规代码将需要修改，以利用使用各向异性细化创建的网格的特征。对于那些对deal.II如何演变感兴趣的人来说，可能有兴趣的是，这种不变性的松动需要一些不兼容的变化。例如，库中曾经有一个成员 GeometryInfo<dim>::children_per_cell ，指定单元格一旦被细化后有多少个孩子。对于各向同性的细化，这个数字等于 $2^{dim}$ ，如上所述。然而，对于各向异性的细化，这个数字并不存在，因为在二维中可以是2或4，在三维中可以是2、4或8，因此成员 GeometryInfo<dim>::children_per_cell 已被删除。它现在被 GeometryInfo<dim>::max_children_per_cell 所取代，后者规定了一个单元格可以拥有的<i>maximum</i>个子女。以前，被细化的单元有多少个孩子是作为静态信息提供的，但现在它取决于单元的实际细化状态，可以使用 TriaAccessor::n_children(), 来获取，这个调用对各向同性和各向异性的细化都同样有效。对于面和它们的子面也有非常类似的情况：相关的信息可以使用 GeometryInfo<dim>::max_children_per_face 或 <code>face->n_children()</code> 进行查询，这取决于上下文。
* 另一个重要的方面，也是本教程中最重要的方面，是在组装单元格之间的表面上的跳跃项时对相邻关系的处理。查看 step-12 中assemble_system函数的文档，我们注意到，我们需要决定一个相邻的单元是否更粗、更细或者与我们当前的单元处于同一（细化）水平。这些决定对于各向异性的细化并不适用，因为单元的 <em> 级 </em> 所提供的信息并不足以完全描述各向异性的单元；例如，一个二维单元的终端子女是在 $x$ 方向切割的二维单元，其子女在 $y$ 方向上被切割，是在第2层，还是在第1层，因为如果该单元被各向同性地细化一次，就会产生同一组最好的单元？
* 在各向异性细化之后，一个更粗的邻居不一定完全是比我们低一级的，而是几乎可以有相对于当前级别的任何级别；事实上，它甚至可以在一个更高的级别上，尽管它更粗。因此，决定必须在不同的基础上做出，而决定的意图却保持不变。
* 在下文中，我们将讨论当我们想要计算形式为@f[
  \int_{\partial K} \varphi_i(x) \varphi_j(x) \; dx
@f]或类似的矩阵（或右手边）的贡献时可能发生的情况；记住，我们使用FEFaceValues和FESubfaceValues类来整合这样的条款。我们还将展示如何编写适用于各向同性和各向异性细化的代码。
*  <ul> 
*  <li>   <em>  细化邻居  </em>  ：如果我们在一个活动单元上，想在一个面 $f\subset \partial K$ 上进行整合，第一种可能性是这个面后面的邻居更细化，即有孩子只占据了公共面的一部分。在这种情况下，所考虑的面必须是一个精致的面，这可以通过询问  <code>if (face->has_children())</code>  来确定。如果这是真的，我们需要循环所有的子面，并得到这个子面后面的邻居的子，这样我们就可以用邻居重新输入一个FEFaceValues对象，用我们的单元格和相应的子面重新输入一个FESubfaceValues对象。
* 对于各向同性的细化，这种情况相当简单，因为我们知道在deal.II中各向同性细化的自适应网格的一个不变因素是，邻居只能正好相差一个细化等级。然而，对于各向异性细化的网格来说，这并不完全正确，特别是在三维中；在那里，我们感兴趣的 $f$ 另一边的活动单元实际上可能不是我们邻居的孩子，而可能是孙子或甚至更远的后代。幸运的是，这种复杂性被隐藏在库的内部。我们所要做的就是调用 CellAccessor::neighbor_child_on_subface() 函数。尽管如此，在3D中，有两种情况需要特别考虑。    <ul>   <li>  如果邻居被各向异性地细化了一次以上，可能这里需要考虑的不是两个或四个而是三个子面。想象一下我们正在考虑的（三维）邻接单元的（二维）面的以下细化过程：首先该面沿x方向细化，后来只沿y方向细化左侧子面。
* @code
* 


* 
* 


* 
* 
-------*
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
* 
---*---*
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
* 
---*---*
   |       |        |   |   |        |   |   |
   |       |
* 


* 
* 
--->  |   |   |
* 


* 
* 
--->
* 
---*   |
   |       |        |   |   |        |   |   |
* 


* 
* 


* 
* 
-------*
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
* 
---*---*
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
* 
---*---*
@endcode
* 这里子面的数量是三个。需要注意的是，对于一个面， TriaAccessor::n_children() 和 TriaAccessor::n_active_descendants(). 之间的细微差别。第一个函数返回直系子代的数量，对于上面的例子是两个，而第二个函数返回活动子代的数量（即包括子代、孙代和进一步的子代），在上面的例子中是正确的三个。使用 <code>face->n_active_descendants()</code> 对各向同性和各向异性以及二维和三维情况都有效，所以应该始终使用它。应该注意的是，如果最右边图像左侧的两个小子面后面的任何一个单元被进一步细化，那么当前的单元（即我们正在查看这个共同面的一侧）也要被细化：之所以这样，是因为否则就会违反每条边只有一个悬挂节点的不变性。
*  <li>  可能的情况是，邻居比较粗糙，但仍有比我们当前单元格更精细的子节点。如果两个同样粗糙的单元被细化，其中一个单元在所考虑的面有两个孩子，另一个有四个孩子，这种情况就会发生。下图中的单元格只是相互分离，以显示各个细化的情况。
* @code
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
-----------*
* 


* 
* 


* 
* 


* 
* 
-----------*
     /           /|    /           /|
    ############# |   +++++++++++++ |
   #           ## |  +           ++
  ############# # | +++++++++++++ +/|
  #           # # | +           + + |
  #           # # +           +++
  #           # #/  +++++++++++++ +/
  #           # #   +           + +
  #           ##    +           ++
  #############     +++++++++++++
@endcode
*
* 这里，左边的两个单元是在 $y$ 方向对母单元进行各向异性分割的结果，而右边的四个单元是在 $y$ -和 $z$ 方向同时进行各向异性细化的结果。  标有#的左边单元有两个标有+的更细的邻居，但左边单元的实际邻居是完整的右边母单元，因为标有+的两个单元更细，它们的直接母体是一个大单元。    </ul> 
* 然而，幸运的是， CellAccessor::neighbor_child_on_subface() 可以自己处理这些情况，如果你在正确的子面数量上循环，在上面的例子中，这是两个。 FESubfaceValues<dim>::reinit 函数也会照顾到这一点，所以结果的状态总是正确的。然而，有一个小的注意事项。为了重新调用邻居的FEFaceValues对象，你需要知道指向当前单元格的面的索引。通常你假设直接得到的邻居和你一样粗或一样细，如果它有孩子的话，因此这个信息可以通过 CellAccessor::neighbor_of_neighbor(). 得到。为了方便你，有一个 CellAccessor::neighbor_face_no() 可以为你做正确的事情，并返回所需的结果。
*  <li>   <em>  邻居和我们的单元格一样细  </em>  ：在我们排除了所有存在更细的子单元的情况后，我们只需要决定，这里的邻居是否更粗。为此，有一个 CellAccessor::neighbor_is_coarser() 函数，返回一个布尔值。为了得到相同粗度的邻居的相关情况，我们将使用  <code>else if (!cell->neighbor_is_coarser(face_no))</code>  。这个块里面的代码可以不动。然而，这里有一件事要提到。如果我们想使用一个规则，哪一个单元应该在一个给定的面上组合某些条款，我们可以考虑 step-12 中提出的规则。我们知道，我们必须舍弃将我们的单元格的水平与邻居的水平进行比较的部分，而代之以上面提出的对更粗的邻居的测试。然而，我们也必须考虑到具有相同粗度的相邻单元具有相同指数（在不同水平上）的可能性。因此，我们必须包括单元格具有相同指数的情况，并给出一个额外的条件，即哪一个单元格应该集合条款，例如，我们可以选择较低层次的单元格。这个概念的细节可以在下面的实现中看到。
*  <li>   <em>  较粗的邻居  </em>  : 剩下的情况很明显：如果没有精炼的邻居，而且邻居没有当前单元格那么细，那么它一定是较粗的。因此，我们可以留下旧的条件短语，简单地使用  <code>else</code>  。 CellAccessor::neighbor_of_coarser_neighbor() 函数照顾到各向异性细化的所有复杂性，结合一般三维网格上可能出现的非标准面方向、翻转和旋转。
*  </ul> 
* <a name="Meshsmoothing"></a><h4>Mesh smoothing</h4> 。
* 当一个三角形被细化时，没有被标记为细化的单元仍然可以被细化。这是由于额外的平滑算法，这些算法是必要的或明确要求的。特别是，在每条边上最多有一个悬空节点的限制，经常迫使细化那些已经很细的、被标记为要进一步细化的单元的邻近单元。
* 然而，deal.II也实现了一些算法，以确保产生的网格比最低限度的平滑，例如确保没有孤立的细化单元被非细化单元包围，因为这些土地上的额外自由度几乎都受到悬挂节点的约束。(参见Triangulation类的文档和它的 Triangulation::MeshSmoothing 成员以获得更多关于网格平滑的信息)。
* 大多数最初为各向同性开发的平滑算法已经被调整为以非常相似的方式用于各向异性和各向同性的细化。然而，有两种算法值得一提。 <ol>   <li>   <code>MeshSmoothing::limit_level_difference_at_vertices</code>  ：在各向同性的环境中，该算法试图通过减少在共同顶点相遇的单元的细化水平差异来确保良好的近似质量。然而，对于各向异性的细化没有明确的对应概念，因此该算法不能与各向异性的细化结合使用。这个限制是由一个断言强制执行的，一旦在一个已经被各向异性细化的三角形上调用该算法，就会抛出一个错误。
*  <li>   <code>MeshSmoothing::allow_anisotropic_smoothing</code>  ：如果引入细化来限制悬空节点的数量，往往不需要额外的单元来提高近似质量。这对DG方法来说尤其如此。如果你设置了标志 <code>allow_anisotropic_smoothing</code> ，平滑算法试图通过使用各向异性的细化来尽量减少可能不需要的额外单元的数量。如果你设置了这个平滑标志，你可能会得到各向异性的细化单元，即使你从未将一个细化标志设置为各向异性的细化。请注意，如果你的代码尊重各向异性网格的可能性，你只应该使用这个标志。结合一个合适的各向异性指标，这个标志可以帮助节省额外的单元，从而节省精力。 </ol> 
*

*<a name="Jumpindicator"></a><h3>Jump indicator</h3>


* 利用各向异性细化的好处，需要一个指标来捕捉溶液的各向异性特征，并利用它们来进行细化过程。一般来说，各向异性的细化过程将包括几个步骤。 <ol>   <li>  计算一个误差指标。    <li>  使用误差指标标记单元进行细化，例如，使用固定数量或分数的单元。这些单元将被自动标记为各向同性的细化。    <li>  仅在被标记的单元上评估一个明显的各向异性指标。    <li>  使用各向异性指标为合适的单元设置一个新的各向异性细化标志，否则保持标志不变。    <li>  调用 Triangulation<dim>::execute_coarsening_and_refinement 来执行要求的细化，使用要求的各向同性和各向异性标志。 </ol>  这种方法类似于我们在 step-27 中用于hp-细化的方法，具有很大的灵活性优势。任何误差指标都可以在各向异性过程中使用，也就是说，如果你有相当多的面向后方的误差指标可用，你可以像使用一个简单的Kelly误差估计器一样容易地使用它们。细化过程的各向异性部分不受这种选择的影响。此外，在deal.II或你的应用程序中没有任何各向异性变化之前，只需省去第三步和第四步，就可以得到相同的各向异性细化结果。作为最后一个优点，只在标记为细化的单元上工作，可以更快地评估各向异性指标，如果指标涉及很多单元，在更细的网格上会变得很明显。
* 在这里，我们使用一个非常简单的方法，它只适用于DG方法。一般的想法是非常简单的。DG方法允许离散解在单元面上跳跃，而在每个单元内是平滑的。当然，在极限情况下，我们希望随着我们对网格的细化和对真实解的逼近，跳跃会趋于零。因此，在某个面上的大跳跃表明该单元应该被细化（至少是）正交于该面，而小跳跃则不会导致这一结论的产生。当然，确切的解不可能是光滑的，它也有跳跃的特征。然而，在这种情况下，一个面的大跳跃表明，这个面或多或少地与跳跃平行，并与之相邻，因此，我们再次期望与所考虑的面正交的细化能够有效。
* 建议的指标计算平均跳跃 $K_j$ ，即离散解 $u$ 在两个面 $f_i^j$ 、 $i=1,2$ 、 $j=1..d$ 上正交于坐标方向 $j$ 的绝对跳跃的平均值@f[
K_j = \frac{\sum_{i=1}^2 \int_{f_i^j}|[u]| dx}{\sum_{i=1}^2 |f_i^j|} .
@f] 。如果一个方向的平均跳动比其他方向的平均跳动大一定的系数 $\kappa$ ，即如果 $K_i > \kappa \frac 1{d-1} \sum_{j=1, j\neq i}^d K_j$ ，则单元格只沿该特定方向细化 $i$ ，否则单元格是各向同性细化的。
* 这样的标准很容易被推广到方程组中：跳跃的绝对值将被矢量值跳跃的适当规范所取代。
*


*<a name="Theproblem"></a><h3>The problem</h3>


* 我们解决 step-12 中提出的线性传输方程。域被扩展到覆盖二维的 $[-1,1]\times[0,1]$ ，其中流场 $\beta$ 在域的右半部分描述了围绕原点的逆时针四分之一圆，在域的左半部分平行于X轴。流入边界同样位于 $x=1$ 处，并沿x轴的正向部分，边界条件选择如 step-12 。
*

* <a name="CommProg"></a> <h1> The commented program</h1>。
* deal.II的包含文件在以前的例子中已经介绍过了，因此不再做进一步评论。
*


* @code
 #include <deal.II/base/function.h>
 #include <deal.II/base/quadrature_lib.h>
 #include <deal.II/base/timer.h>
 #include <deal.II/lac/precondition_block.h>
 #include <deal.II/lac/solver_richardson.h>
 #include <deal.II/lac/sparse_matrix.h>
 #include <deal.II/lac/vector.h>
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/grid_out.h>
 #include <deal.II/grid/grid_refinement.h>
 #include <deal.II/dofs/dof_handler.h>
 #include <deal.II/dofs/dof_tools.h>
 #include <deal.II/fe/fe_values.h>
 #include <deal.II/fe/mapping_q1.h>
 #include <deal.II/fe/fe_dgq.h>
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/numerics/derivative_approximation.h>
* 
 @endcode

* 而这又是C++。
*


* @code
 #include <array>
 #include <iostream>
 #include <fstream>
* 
 @endcode

* 最后一步和以前的所有程序一样。



* @code
 namespace Step30
 {
   using namespace dealii;
* 
 @endcode
*
* <a name="Equationdata"></a> <h3>Equation data</h3>.
* 描述方程数据的类和单个项的实际装配几乎完全照搬自  step-12  。我们将对差异进行评论。
*


* @code
   template <int dim>
   class RHS : public Function<dim>
   {
   public:
     virtual void value_list(const std::vector<Point<dim>> &points,
                             std::vector<double> &          values,
                             const unsigned int  /*component*/  = 0) const override
     {
       (void)points;
       Assert(values.size() == points.size(),
              ExcDimensionMismatch(values.size(), points.size()));
* 
       std::fill(values.begin(), values.end(), 0.);
     }
   };
* 
 
   template <int dim>
   class BoundaryValues : public Function<dim>
   {
   public:
     virtual void value_list(const std::vector<Point<dim>> &points,
                             std::vector<double> &          values,
                             const unsigned int  /*component*/  = 0) const override
     {
       Assert(values.size() == points.size(),
              ExcDimensionMismatch(values.size(), points.size()));
* 
       for (unsigned int i = 0; i < values.size(); ++i)
         {
           if (points[i](0) < 0.5)
             values[i] = 1.;
           else
             values[i] = 0.;
         }
     }
   };
* 
 
   template <int dim>
   class Beta
   {
   public:
 @endcode
*
*流场选择为四分之一圆，流向为逆时针，原点为域的右半部分的中点，数值为正 $x$ ，而在域的左边部分，流速只是向左走，与从右边进来的流速一致。在圆形部分，流速的大小与离原点的距离成正比。这与 step-12 不同，在该定义中，到处都是1。新定义导致 $\beta$ 沿单元的每个给定面的线性变化。另一方面， $u(x,y)$ 的解决方案与之前完全相同。
*


* @code
     void value_list(const std::vector<Point<dim>> &points,
                     std::vector<Point<dim>> &      values) const
     {
       Assert(values.size() == points.size(),
              ExcDimensionMismatch(values.size(), points.size()));
* 
       for (unsigned int i = 0; i < points.size(); ++i)
         {
           if (points[i](0) > 0)
             {
               values[i](0) =
* 
-points[i](1);
               values[i](1) = points[i](0);
             }
           else
             {
               values[i]    = Point<dim>();
               values[i](0) =
* 
-points[i](1);
             }
         }
     }
   };
* 
 
* 
 @endcode
*
* <a name="ClassDGTransportEquation"></a> <h3>Class: DGTransportEquation</h3>。
* 这个类的声明完全不受我们目前的变化影响。
*


* @code
   template <int dim>
   class DGTransportEquation
   {
   public:
     DGTransportEquation();
* 
     void assemble_cell_term(const FEValues<dim> &fe_v,
                             FullMatrix<double> & ui_vi_matrix,
                             Vector<double> &     cell_vector) const;
* 
     void assemble_boundary_term(const FEFaceValues<dim> &fe_v,
                                 FullMatrix<double> &     ui_vi_matrix,
                                 Vector<double> &         cell_vector) const;
* 
     void assemble_face_term(const FEFaceValuesBase<dim> &fe_v,
                             const FEFaceValuesBase<dim> &fe_v_neighbor,
                             FullMatrix<double> &         ui_vi_matrix,
                             FullMatrix<double> &         ue_vi_matrix,
                             FullMatrix<double> &         ui_ve_matrix,
                             FullMatrix<double> &         ue_ve_matrix) const;
* 
   private:
     const Beta<dim>           beta_function;
     const RHS<dim>            rhs_function;
     const BoundaryValues<dim> boundary_function;
   };
* 
 
* 
 @endcode
*
* 同样地，该类的构造函数以及组装对应于单元格内部和边界面的术语的函数与之前没有变化。装配单元间面术语的函数也没有改变，因为它所做的只是对两个FEFaceValuesBase类型的对象（它是FEFaceValues和FESubfaceValues的基类）进行操作。这些对象从何而来，即它们是如何被初始化的，对这个函数来说并不重要：它只是假设这两个对象所代表的面或子面上的正交点对应于物理空间中的相同点。
*


* @code
   template <int dim>
   DGTransportEquation<dim>::DGTransportEquation()
     : beta_function()
     , rhs_function()
     , boundary_function()
   {}
* 
 
* 
   template <int dim>
   void DGTransportEquation<dim>::assemble_cell_term(
     const FEValues<dim> &fe_v,
     FullMatrix<double> & ui_vi_matrix,
     Vector<double> &     cell_vector) const
   {
     const std::vector<double> &JxW = fe_v.get_JxW_values();
* 
     std::vector<Point<dim>> beta(fe_v.n_quadrature_points);
     std::vector<double>     rhs(fe_v.n_quadrature_points);
* 
     beta_function.value_list(fe_v.get_quadrature_points(), beta);
     rhs_function.value_list(fe_v.get_quadrature_points(), rhs);
* 
     for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
       for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
         {
           for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
             ui_vi_matrix(i, j)
* 
-= beta[point] fe_v.shape_grad(i, point)
                                   fe_v.shape_value(j, point) JxW[point];
* 
           cell_vector(i) +=
             rhs[point] fe_v.shape_value(i, point) JxW[point];
         }
   }
* 
 
   template <int dim>
   void DGTransportEquation<dim>::assemble_boundary_term(
     const FEFaceValues<dim> &fe_v,
     FullMatrix<double> &     ui_vi_matrix,
     Vector<double> &         cell_vector) const
   {
     const std::vector<double> &        JxW     = fe_v.get_JxW_values();
     const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();
* 
     std::vector<Point<dim>> beta(fe_v.n_quadrature_points);
     std::vector<double>     g(fe_v.n_quadrature_points);
* 
     beta_function.value_list(fe_v.get_quadrature_points(), beta);
     boundary_function.value_list(fe_v.get_quadrature_points(), g);
* 
     for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
       {
         const double beta_n = beta[point] normals[point];
         if (beta_n > 0)
           for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
             for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
               ui_vi_matrix(i, j) += beta_n fe_v.shape_value(j, point)
                                     fe_v.shape_value(i, point) JxW[point];
         else
           for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
             cell_vector(i)
* 
-=
               beta_n g[point] fe_v.shape_value(i, point) JxW[point];
       }
   }
* 
 
   template <int dim>
   void DGTransportEquation<dim>::assemble_face_term(
     const FEFaceValuesBase<dim> &fe_v,
     const FEFaceValuesBase<dim> &fe_v_neighbor,
     FullMatrix<double> &         ui_vi_matrix,
     FullMatrix<double> &         ue_vi_matrix,
     FullMatrix<double> &         ui_ve_matrix,
     FullMatrix<double> &         ue_ve_matrix) const
   {
     const std::vector<double> &        JxW     = fe_v.get_JxW_values();
     const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();
* 
     std::vector<Point<dim>> beta(fe_v.n_quadrature_points);
* 
     beta_function.value_list(fe_v.get_quadrature_points(), beta);
* 
     for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
       {
         const double beta_n = beta[point] normals[point];
         if (beta_n > 0)
           {
             for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
               for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                 ui_vi_matrix(i, j) += beta_n fe_v.shape_value(j, point)
                                       fe_v.shape_value(i, point) JxW[point];
* 
             for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
               for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                 ui_ve_matrix(k, j)
* 
-= beta_n fe_v.shape_value(j, point)
                                       fe_v_neighbor.shape_value(k, point)
                                       JxW[point];
           }
         else
           {
             for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
               for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
                 ue_vi_matrix(i, l) += beta_n
                                       fe_v_neighbor.shape_value(l, point)
                                       fe_v.shape_value(i, point) JxW[point];
* 
             for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
               for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
                 ue_ve_matrix(k, l)
* 
-=
                   beta_n fe_v_neighbor.shape_value(l, point)
                   fe_v_neighbor.shape_value(k, point) JxW[point];
           }
       }
   }
* 
 
 @endcode
*
* <a name="ClassDGMethod"></a> <h3>Class: DGMethod</h3>。
* 这个声明很像  step-12  的声明。然而，我们引入了一个新的例程（set_anisotropic_flags）并修改了另一个例程（refine_grid）。
*


* @code
   template <int dim>
   class DGMethod
   {
   public:
     DGMethod(const bool anisotropic);
* 
     void run();
* 
   private:
     void setup_system();
     void assemble_system();
     void solve(Vector<double> &solution);
     void refine_grid();
     void set_anisotropic_flags();
     void output_results(const unsigned int cycle) const;
* 
     Triangulation<dim>   triangulation;
     const MappingQ1<dim> mapping;
 @endcode
*
* 我们再次希望使用1度的DG元素（但这只在构造函数中指定）。如果你想使用不同程度的DG方法，请在构造函数中用新的程度替换1。
*


* @code
     const unsigned int degree;
     FE_DGQ<dim>        fe;
     DoFHandler<dim>    dof_handler;
* 
     SparsityPattern      sparsity_pattern;
     SparseMatrix<double> system_matrix;
 @endcode
*
* 这是新的，在介绍中解释的各向异性跳跃指标的评估中使用的阈值。它的值在构造函数中被设置为3.0，但是它可以很容易地被改变为一个大于1的不同值。
*


* @code
     const double anisotropic_threshold_ratio;
 @endcode
*
* 这是一个指示是否使用各向异性细化的bool标志。它由构造函数设置，构造函数需要一个同名的参数。
*


* @code
     const bool anisotropic;
* 
     const QGauss<dim>     quadrature;
     const QGauss<dim
* 
- 1> face_quadrature;
* 
     Vector<double> solution2;
     Vector<double> right_hand_side;
* 
     const DGTransportEquation<dim> dg;
   };
* 
 
   template <int dim>
   DGMethod<dim>::DGMethod(const bool anisotropic)
     : mapping()
     ,
 @endcode

* 在这里为不同程度的DG方法进行改变。
*


* @code
     degree(1)
     , fe(degree)
     , dof_handler(triangulation)
     , anisotropic_threshold_ratio(3.)
     , anisotropic(anisotropic)
     ,
 @endcode
*
* 由于β是一个线性函数，我们可以选择所得积分正确的正交度数。因此，我们选择使用 <code>degree+1</code> 高斯点，这使我们能够准确地积分度数为 <code>2*degree+1</code> 的多项式，足以满足我们在本程序中要进行的所有积分。
*


* @code
     quadrature(degree + 1)
     , face_quadrature(degree + 1)
     , dg()
   {}
* 
 
* 
   template <int dim>
   void DGMethod<dim>::setup_system()
   {
     dof_handler.distribute_dofs(fe);
     sparsity_pattern.reinit(dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             (GeometryInfo<dim>::faces_per_cell
                                GeometryInfo<dim>::max_children_per_face +
                              1)
                               fe.n_dofs_per_cell());
* 
     DoFTools::make_flux_sparsity_pattern(dof_handler, sparsity_pattern);
* 
     sparsity_pattern.compress();
* 
     system_matrix.reinit(sparsity_pattern);
* 
     solution2.reinit(dof_handler.n_dofs());
     right_hand_side.reinit(dof_handler.n_dofs());
   }
* 
 
 @endcode
*
* <a name="Functionassemble_system"></a> <h4>Function: assemble_system</h4>。
* 我们继续使用 <code>assemble_system</code> 函数来实现DG离散化。这个函数与 step-12 中的 <code>assemble_system</code> 函数的作用相同（但没有MeshWorker）。  一个单元的邻居关系所考虑的四种情况与各向同性的情况相同，即a)单元在边界上，b)有更细的邻居单元，c)邻居既不粗也不细，d)邻居更粗。  然而，我们决定哪种情况的方式是按照介绍中描述的方式修改的。
*


* @code
   template <int dim>
   void DGMethod<dim>::assemble_system()
   {
     const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
     std::vector<types::global_dof_index> dofs(dofs_per_cell);
     std::vector<types::global_dof_index> dofs_neighbor(dofs_per_cell);
* 
     const UpdateFlags update_flags = update_values | update_gradients |
                                      update_quadrature_points |
                                      update_JxW_values;
* 
     const UpdateFlags face_update_flags =
       update_values | update_quadrature_points | update_JxW_values |
       update_normal_vectors;
* 
     const UpdateFlags neighbor_face_update_flags = update_values;
* 
     FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags);
     FEFaceValues<dim>    fe_v_face(mapping,
                                 fe,
                                 face_quadrature,
                                 face_update_flags);
     FESubfaceValues<dim> fe_v_subface(mapping,
                                       fe,
                                       face_quadrature,
                                       face_update_flags);
     FEFaceValues<dim>    fe_v_face_neighbor(mapping,
                                          fe,
                                          face_quadrature,
                                          neighbor_face_update_flags);
* 
 
     FullMatrix<double> ui_vi_matrix(dofs_per_cell, dofs_per_cell);
     FullMatrix<double> ue_vi_matrix(dofs_per_cell, dofs_per_cell);
* 
     FullMatrix<double> ui_ve_matrix(dofs_per_cell, dofs_per_cell);
     FullMatrix<double> ue_ve_matrix(dofs_per_cell, dofs_per_cell);
* 
     Vector<double> cell_vector(dofs_per_cell);
* 
     for (const auto &cell : dof_handler.active_cell_iterators())
       {
         ui_vi_matrix = 0;
         cell_vector  = 0;
* 
         fe_v.reinit(cell);
* 
         dg.assemble_cell_term(fe_v, ui_vi_matrix, cell_vector);
* 
         cell->get_dof_indices(dofs);
* 
         for (const auto face_no : cell->face_indices())
           {
             const auto face = cell->face(face_no);
* 
 @endcode
*
* 情况（a）。该面在边界上。
*


* @code
             if (face->at_boundary())
               {
                 fe_v_face.reinit(cell, face_no);
* 
                 dg.assemble_boundary_term(fe_v_face, ui_vi_matrix, cell_vector);
               }
             else
               {
                 Assert(cell->neighbor(face_no).state() == IteratorState::valid,
                        ExcInternalError());
                 const auto neighbor = cell->neighbor(face_no);
* 
 @endcode
*
* 情况（b）。这是一个内部面，邻居是精炼的（我们可以通过询问当前单元格的面是否有孩子来测试）。在这种情况下，我们需要对 "子面 "进行整合，即当前单元格的面的子女。                 
*（有一个稍微令人困惑的角落案例。如果我们是在1d
*
*--诚然，当前的程序及其各向异性细化的演示并不特别相关
*
* - 那么单元格之间的面总是相同的：它们只是顶点。换句话说，在1d中，我们不想区别对待不同层次的单元之间的面。我们在这里检查的条件`face->has_children()`确保了这一点：在1d中，这个函数总是返回`false'，因此在1d中我们永远不会进入这个`if'分支。但我们将不得不在下面的情况（c）中回到这个角落。
*


* @code
                 if (face->has_children())
                   {
 @endcode
*
* 我们需要知道，哪个邻居的面孔指向我们的单元格的方向。使用 @p  neighbor_face_no函数，我们可以得到较粗和非较粗邻居的这些信息。
*


* @code
                     const unsigned int neighbor2 =
                       cell->neighbor_face_no(face_no);
* 
 @endcode
*
* 现在我们在所有的子面中循环，也就是当前面的子代和可能的孙代。
*


* @code
                     for (unsigned int subface_no = 0;
                          subface_no < face->n_active_descendants();
                          ++subface_no)
                       {
 @endcode
*
* 要获得当前子面后面的单元，我们可以使用 @p neighbor_child_on_subface 函数。它照顾到各向异性细化和非标准面的所有复杂情况。
*


* @code
                         const auto neighbor_child =
                           cell->neighbor_child_on_subface(face_no, subface_no);
                         Assert(!neighbor_child->has_children(),
                                ExcInternalError());
* 
 @endcode
*
* 本案例的其余部分没有变化。
*


* @code
                         ue_vi_matrix = 0;
                         ui_ve_matrix = 0;
                         ue_ve_matrix = 0;
* 
                         fe_v_subface.reinit(cell, face_no, subface_no);
                         fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
* 
                         dg.assemble_face_term(fe_v_subface,
                                               fe_v_face_neighbor,
                                               ui_vi_matrix,
                                               ue_vi_matrix,
                                               ui_ve_matrix,
                                               ue_ve_matrix);
* 
                         neighbor_child->get_dof_indices(dofs_neighbor);
* 
                         for (unsigned int i = 0; i < dofs_per_cell; ++i)
                           for (unsigned int j = 0; j < dofs_per_cell; ++j)
                             {
                               system_matrix.add(dofs[i],
                                                 dofs_neighbor[j],
                                                 ue_vi_matrix(i, j));
                               system_matrix.add(dofs_neighbor[i],
                                                 dofs[j],
                                                 ui_ve_matrix(i, j));
                               system_matrix.add(dofs_neighbor[i],
                                                 dofs_neighbor[j],
                                                 ue_ve_matrix(i, j));
                             }
                       }
                   }
                 else
                   {
 @endcode
*
* 情况（c）。如果这是一个内部面，并且邻居没有进一步细化，我们就得到这里（或者，如上所述，我们是在1d中，在这种情况下，我们对每个内部面都得到这里）。然后，我们需要决定是否要对当前面进行整合。如果邻居实际上更粗，那么我们就忽略这个面，而是在访问相邻的单元并查看当前面时进行处理（除了在1d中，如上所述，这是不可能的）。
*


* @code
                     if (dim > 1 && cell->neighbor_is_coarser(face_no))
                       continue;
* 
 @endcode
*
* 另一方面，如果邻居更加精细，那么我们已经处理了上面(b)情况下的面（1d除外）。所以对于2d和3d，我们只需要决定是要处理来自当前一侧的同一层次的单元格之间的面还是来自相邻一侧的面。  我们通过引入一个平局来做到这一点。我们只取索引较小的单元格（在当前细化级别内）。在1d中，我们取较粗的单元，或者如果它们在同一层次，则取该层次中指数较小的单元。这就导致了一个复杂的条件，希望在上面的描述中可以理解。
*


* @code
                     if (((dim > 1) && (cell->index() < neighbor->index())) ||
                         ((dim == 1) && ((cell->level() < neighbor->level()) ||
                                         ((cell->level() == neighbor->level()) &&
                                          (cell->index() < neighbor->index())))))
                       {
 @endcode
*
* 这里我们知道，邻居不是更粗的，所以我们可以使用通常的 @p neighbor_of_neighbor 函数。然而，我们也可以使用更通用的 @p neighbor_face_no 函数。
*


* @code
                         const unsigned int neighbor2 =
                           cell->neighbor_of_neighbor(face_no);
* 
                         ue_vi_matrix = 0;
                         ui_ve_matrix = 0;
                         ue_ve_matrix = 0;
* 
                         fe_v_face.reinit(cell, face_no);
                         fe_v_face_neighbor.reinit(neighbor, neighbor2);
* 
                         dg.assemble_face_term(fe_v_face,
                                               fe_v_face_neighbor,
                                               ui_vi_matrix,
                                               ue_vi_matrix,
                                               ui_ve_matrix,
                                               ue_ve_matrix);
* 
                         neighbor->get_dof_indices(dofs_neighbor);
* 
                         for (unsigned int i = 0; i < dofs_per_cell; ++i)
                           for (unsigned int j = 0; j < dofs_per_cell; ++j)
                             {
                               system_matrix.add(dofs[i],
                                                 dofs_neighbor[j],
                                                 ue_vi_matrix(i, j));
                               system_matrix.add(dofs_neighbor[i],
                                                 dofs[j],
                                                 ui_ve_matrix(i, j));
                               system_matrix.add(dofs_neighbor[i],
                                                 dofs_neighbor[j],
                                                 ue_ve_matrix(i, j));
                             }
                       }
* 
 @endcode
*
* 我们不需要考虑情况(d)，因为那些面孔是'从情况(b)内的另一边处理的。
*


* @code
                   }
               }
           }
* 
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
           for (unsigned int j = 0; j < dofs_per_cell; ++j)
             system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i, j));
* 
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
           right_hand_side(dofs[i]) += cell_vector(i);
       }
   }
* 
 
 @endcode
*
* <a name="Solver"></a> <h3>Solver</h3>。
* 对于这个简单的问题，我们再次使用简单的Richardson迭代。该求解器完全不受我们各向异性变化的影响。
*


* @code
   template <int dim>
   void DGMethod<dim>::solve(Vector<double> &solution)
   {
     SolverControl                    solver_control(1000, 1e-12, false, false);
     SolverRichardson<Vector<double>> solver(solver_control);
* 
     PreconditionBlockSSOR<SparseMatrix<double>> preconditioner;
* 
     preconditioner.initialize(system_matrix, fe.n_dofs_per_cell());
* 
     solver.solve(system_matrix, solution, right_hand_side, preconditioner);
   }
* 
 
 @endcode
*
* <a name="Refinement"></a> <h3>Refinement</h3>.
* 我们根据 step-12 中使用的相同的简单细化标准来细化网格，即对解的梯度的近似。
*


* @code
   template <int dim>
   void DGMethod<dim>::refine_grid()
   {
     Vector<float> gradient_indicator(triangulation.n_active_cells());
* 
 @endcode
*
* 我们对梯度进行近似计算。
*


* @code
     DerivativeApproximation::approximate_gradient(mapping,
                                                   dof_handler,
                                                   solution2,
                                                   gradient_indicator);
* 
 @endcode

* 并对其进行缩放，以获得一个误差指标。
*


* @code
     for (const auto &cell : triangulation.active_cell_iterators())
       gradient_indicator[cell->active_cell_index()]=
         std::pow(cell->diameter(), 1 + 1.0 dim / 2);
 @endcode
*
* 然后我们用这个指标来标记误差指标最高的30%的单元格来进行精炼。
*


* @code
     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                     gradient_indicator,
                                                     0.3,
                                                     0.1);
 @endcode
*
* 现在，细化标志被设置为那些具有大误差指标的单元。如果不做任何改变，这些单元将被各向同性地细化。如果给这个函数的 @p anisotropic 标志被设置，我们现在调用set_anisotropic_flags()函数，该函数使用跳转指标将一些细化标志重置为各向异性细化。
*


* @code
     if (anisotropic)
       set_anisotropic_flags();
 @endcode
*
* 现在执行细化，考虑各向异性以及各向同性的细化标志。
*


* @code
     triangulation.execute_coarsening_and_refinement();
   }
* 
 @endcode
*
* 一旦错误指标被评估，误差最大的单元被标记为细化，我们要再次循环这些被标记的单元，以决定它们是否需要各向同性的细化，或者各向异性的细化更合适。这就是在介绍中解释的各向异性跳跃指标。
*


* @code
   template <int dim>
   void DGMethod<dim>::set_anisotropic_flags()
   {
 @endcode
*
* 我们想在有标记的单元格的面上评估跳跃，所以我们需要一些对象来评估面上的解决方案的值。
*


* @code
     UpdateFlags face_update_flags =
       UpdateFlags(update_values | update_JxW_values);
* 
     FEFaceValues<dim>    fe_v_face(mapping,
                                 fe,
                                 face_quadrature,
                                 face_update_flags);
     FESubfaceValues<dim> fe_v_subface(mapping,
                                       fe,
                                       face_quadrature,
                                       face_update_flags);
     FEFaceValues<dim>    fe_v_face_neighbor(mapping,
                                          fe,
                                          face_quadrature,
                                          update_values);
* 
 @endcode
*
* 现在我们需要在所有的活动单元中进行循环。
*


* @code
     for (const auto &cell : dof_handler.active_cell_iterators())
 @endcode
*
* 我们只需要考虑那些被标记为细化的单元。
*


* @code
       if (cell->refine_flag_set())
         {
           Point<dim> jump;
           Point<dim> area;
* 
           for (const auto face_no : cell->face_indices())
             {
               const auto face = cell->face(face_no);
* 
               if (!face->at_boundary())
                 {
                   Assert(cell->neighbor(face_no).state() ==
                            IteratorState::valid,
                          ExcInternalError());
                   const auto neighbor = cell->neighbor(face_no);
* 
                   std::vector<double> u(fe_v_face.n_quadrature_points);
                   std::vector<double> u_neighbor(fe_v_face.n_quadrature_points);
* 
 @endcode
*
* 在汇编例程中看到的四种不同的相邻关系的情况，在这里基本以相同的方式重复。
*


* @code
                   if (face->has_children())
                     {
 @endcode
*
* 邻居被细化。  首先，我们存储信息，即邻居的哪个面指向我们当前单元的方向。这个属性会继承给子代。
*


* @code
                       unsigned int neighbor2 = cell->neighbor_face_no(face_no);
 @endcode
*
* 现在我们在所有的子界面上循环。
*


* @code
                       for (unsigned int subface_no = 0;
                            subface_no < face->n_active_descendants();
                            ++subface_no)
                         {
 @endcode
*
*得到一个迭代器，指向当前子面后面的单元格...
*


* @code
                           const auto neighbor_child =
                             cell->neighbor_child_on_subface(face_no,
                                                             subface_no);
                           Assert(!neighbor_child->has_children(),
                                  ExcInternalError());
 @endcode
*
* ... 并重新启动各自的FEFaceValues和FESSubFaceValues对象。
*


* @code
                           fe_v_subface.reinit(cell, face_no, subface_no);
                           fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
 @endcode
*
* 我们获得了函数值
*


* @code
                           fe_v_subface.get_function_values(solution2, u);
                           fe_v_face_neighbor.get_function_values(solution2,
                                                                  u_neighbor);
 @endcode
*
*以及正交权重，乘以雅各布行列式。
*


* @code
                           const std::vector<double> &JxW =
                             fe_v_subface.get_JxW_values();
 @endcode
*
* 现在我们在所有的正交点上循环。
*


* @code
                           for (unsigned int x = 0;
                                x < fe_v_subface.n_quadrature_points;
                                ++x)
                             {
 @endcode
*
*并整合解决方案的跳跃的绝对值，即分别从当前单元格和相邻单元格看到的函数值的绝对值。我们知道，前两个面与单元格上的第一个坐标方向正交，后两个面与第二个坐标方向正交，以此类推，所以我们把这些值累积成有 <code>dim</code> 分量的向量。
*


* @code
                               jump[face_no / 2] +=
                                 std::abs(u[x]
* 
- u_neighbor[x]) JxW[x];
 @endcode
*
* 我们还将按比例的权重相加，得到脸部的量度。
*


* @code
                               area[face_no / 2] += JxW[x];
                             }
                         }
                     }
                   else
                     {
                       if (!cell->neighbor_is_coarser(face_no))
                         {
 @endcode
*
* 我们当前的单元和邻近的单元沿着所考虑的面有相同的细化。除此之外，我们的做法与上述案例中的一个子单元基本相同。
*


* @code
                           unsigned int neighbor2 =
                             cell->neighbor_of_neighbor(face_no);
* 
                           fe_v_face.reinit(cell, face_no);
                           fe_v_face_neighbor.reinit(neighbor, neighbor2);
* 
                           fe_v_face.get_function_values(solution2, u);
                           fe_v_face_neighbor.get_function_values(solution2,
                                                                  u_neighbor);
* 
                           const std::vector<double> &JxW =
                             fe_v_face.get_JxW_values();
* 
                           for (unsigned int x = 0;
                                x < fe_v_face.n_quadrature_points;
                                ++x)
                             {
                               jump[face_no / 2] +=
                                 std::abs(u[x]
* 
- u_neighbor[x]) JxW[x];
                               area[face_no / 2] += JxW[x];
                             }
                         }
                       else // i.e. neighbor is coarser than cell
                         {
 @endcode
*
* 现在的邻居实际上更粗了。这种情况是新的，因为它没有在汇编程序中出现。在这里，我们必须考虑它，但这并不太复杂。我们只需使用 @p  neighbor_of_coarser_neighbor函数，它再次自行处理各向异性的细化和非标准面的方向。
*


* @code
                           std::pair<unsigned int, unsigned int>
                             neighbor_face_subface =
                               cell->neighbor_of_coarser_neighbor(face_no);
                           Assert(neighbor_face_subface.first < cell->n_faces(),
                                  ExcInternalError());
                           Assert(neighbor_face_subface.second <
                                    neighbor->face(neighbor_face_subface.first)
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
* 


* 
* 


* 
* 
->n_active_descendants(),
                                  ExcInternalError());
                           Assert(neighbor->neighbor_child_on_subface(
                                    neighbor_face_subface.first,
                                    neighbor_face_subface.second) == cell,
                                  ExcInternalError());
* 
                           fe_v_face.reinit(cell, face_no);
                           fe_v_subface.reinit(neighbor,
                                               neighbor_face_subface.first,
                                               neighbor_face_subface.second);
* 
                           fe_v_face.get_function_values(solution2, u);
                           fe_v_subface.get_function_values(solution2,
                                                            u_neighbor);
* 
                           const std::vector<double> &JxW =
                             fe_v_face.get_JxW_values();
* 
                           for (unsigned int x = 0;
                                x < fe_v_face.n_quadrature_points;
                                ++x)
                             {
                               jump[face_no / 2] +=
                                 std::abs(u[x]
* 
- u_neighbor[x]) JxW[x];
                               area[face_no / 2] += JxW[x];
                             }
                         }
                     }
                 }
             }
 @endcode
*
* 现在我们分析一下平均跳跃的大小，我们用跳跃除以各自面的度量得到。
*


* @code
           std::array<double, dim> average_jumps;
           double                  sum_of_average_jumps = 0.;
           for (unsigned int i = 0; i < dim; ++i)
             {
               average_jumps[i] = jump(i) / area(i);
               sum_of_average_jumps += average_jumps[i];
             }
* 
 @endcode
*
* 现在我们在单元格的 <code>dim</code> 坐标方向上进行循环，比较与该方向正交的面的平均跳跃和与其余方向正交的面的平均跳跃。如果前者比后者大一个给定的系数，我们只沿帽轴进行细化。否则，我们不改变细化标志，导致各向同性的细化。
*


* @code
           for (unsigned int i = 0; i < dim; ++i)
             if (average_jumps[i] > anisotropic_threshold_ratio
                                      (sum_of_average_jumps
* 
- average_jumps[i]))
               cell->set_refine_flag(RefinementCase<dim>::cut_axis(i));
         }
   }
* 
 @endcode
*
* <a name="TheRest"></a> <h3>The Rest</h3>。
* 程序的其余部分非常遵循以前教程程序的方案。我们以VTU格式输出网格（就像我们在 step-1 中所做的那样，例如），并以VTU格式输出可视化，我们几乎总是这样做。
*


* @code
   template <int dim>
   void DGMethod<dim>::output_results(const unsigned int cycle) const
   {
     std::string refine_type;
     if (anisotropic)
       refine_type = ".aniso";
     else
       refine_type = ".iso";
* 
     {
       const std::string filename =
         "grid-" + std::to_string(cycle) + refine_type + ".svg";
       std::cout << "   Writing grid to <" << filename << ">..." << std::endl;
       std::ofstream svg_output(filename);
* 
       GridOut grid_out;
       grid_out.write_svg(triangulation, svg_output);
     }
* 
     {
       const std::string filename =
         "sol-" + std::to_string(cycle) + refine_type + ".vtu";
       std::cout << "   Writing solution to <" << filename << ">..."
                 << std::endl;
       std::ofstream gnuplot_output(filename);
* 
       DataOut<dim> data_out;
       data_out.attach_dof_handler(dof_handler);
       data_out.add_data_vector(solution2, "u");
* 
       data_out.build_patches(degree);
* 
       data_out.write_vtu(gnuplot_output);
     }
   }
* 
 
* 
   template <int dim>
   void DGMethod<dim>::run()
   {
     for (unsigned int cycle = 0; cycle < 6; ++cycle)
       {
         std::cout << "Cycle " << cycle << ':' << std::endl;
* 
         if (cycle == 0)
           {
 @endcode
*
* 创建矩形域。
*


* @code
             Point<dim> p1, p2;
             p1(0) = 0;
             p1(0) =
* 
-1;
             for (unsigned int i = 0; i < dim; ++i)
               p2(i) = 1.;
 @endcode
*
* 调整不同方向的单元数，以获得原始网格的完全各向同性的单元。
*


* @code
             std::vector<unsigned int> repetitions(dim, 1);
             repetitions[0] = 2;
             GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                       repetitions,
                                                       p1,
                                                       p2);
* 
             triangulation.refine_global(5
* 
- dim);
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
         Timer assemble_timer;
         assemble_system();
         std::cout << "   Time of assemble_system: " << assemble_timer.cpu_time()
                   << std::endl;
         solve(solution2);
* 
         output_results(cycle);
* 
         std::cout << std::endl;
       }
   }
 } // namespace Step30
* 
 
* 
 int main()
 {
   try
     {
       using namespace Step30;
* 
 @endcode
*
* 如果你想以3D方式运行程序，只需将下面一行改为 <code>const unsigned int dim = 3;</code>  。
*


* @code
       const unsigned int dim = 2;
* 
       {
 @endcode
*
* 首先，我们进行一次各向同性的细化运行。
*


* @code
         std::cout << "Performing a " << dim
                   << "D run with isotropic refinement..." << std::endl
                   << "------------------------------------------------"
                   << std::endl;
         DGMethod<dim> dgmethod_iso(false);
         dgmethod_iso.run();
       }
* 
       {
 @endcode
*
* 现在我们做第二次运行，这次是用各向异性的细化方法。
*


* @code
         std::cout << std::endl
                   << "Performing a " << dim
                   << "D run with anisotropic refinement..." << std::endl
                   << "--------------------------------------------------"
                   << std::endl;
         DGMethod<dim> dgmethod_aniso(true);
         dgmethod_aniso.run();
       }
     }
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
       return 1;
     }
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
     };
* 
   return 0;
 }
 @endcode
*<a name="Results"></a><h1>Results</h1>


*
* 这个程序的输出包括控制台输出，包含网格的SVG文件，以及以VTU格式给出的解决方案。
* @code
Performing a 2D run with isotropic refinement...
* 
------------------------------------------------
Cycle 0:
   Number of active cells:       128
   Number of degrees of freedom: 512
   Time of assemble_system: 0.092049
   Writing grid to <grid-0.iso.svg>...
   Writing solution to <sol-0.iso.vtu>...
* 
Cycle 1:
   Number of active cells:       239
   Number of degrees of freedom: 956
   Time of assemble_system: 0.109519
   Writing grid to <grid-1.iso.svg>...
   Writing solution to <sol-1.iso.vtu>...
* 
Cycle 2:
   Number of active cells:       491
   Number of degrees of freedom: 1964
   Time of assemble_system: 0.08303
   Writing grid to <grid-2.iso.svg>...
   Writing solution to <sol-2.iso.vtu>...
* 
Cycle 3:
   Number of active cells:       1031
   Number of degrees of freedom: 4124
   Time of assemble_system: 0.278987
   Writing grid to <grid-3.iso.svg>...
   Writing solution to <sol-3.iso.vtu>...
* 
Cycle 4:
   Number of active cells:       2027
   Number of degrees of freedom: 8108
   Time of assemble_system: 0.305869
   Writing grid to <grid-4.iso.svg>...
   Writing solution to <sol-4.iso.vtu>...
* 
Cycle 5:
   Number of active cells:       4019
   Number of degrees of freedom: 16076
   Time of assemble_system: 0.47616
   Writing grid to <grid-5.iso.svg>...
   Writing solution to <sol-5.iso.vtu>...
* 


* 
Performing a 2D run with anisotropic refinement...
* 
--------------------------------------------------
Cycle 0:
   Number of active cells:       128
   Number of degrees of freedom: 512
   Time of assemble_system: 0.052866
   Writing grid to <grid-0.aniso.svg>...
   Writing solution to <sol-0.aniso.vtu>...
* 
Cycle 1:
   Number of active cells:       171
   Number of degrees of freedom: 684
   Time of assemble_system: 0.050917
   Writing grid to <grid-1.aniso.svg>...
   Writing solution to <sol-1.aniso.vtu>...
* 
Cycle 2:
   Number of active cells:       255
   Number of degrees of freedom: 1020
   Time of assemble_system: 0.064132
   Writing grid to <grid-2.aniso.svg>...
   Writing solution to <sol-2.aniso.vtu>...
* 
Cycle 3:
   Number of active cells:       394
   Number of degrees of freedom: 1576
   Time of assemble_system: 0.119849
   Writing grid to <grid-3.aniso.svg>...
   Writing solution to <sol-3.aniso.vtu>...
* 
Cycle 4:
   Number of active cells:       648
   Number of degrees of freedom: 2592
   Time of assemble_system: 0.218244
   Writing grid to <grid-4.aniso.svg>...
   Writing solution to <sol-4.aniso.vtu>...
* 
Cycle 5:
   Number of active cells:       1030
   Number of degrees of freedom: 4120
   Time of assemble_system: 0.128121
   Writing grid to <grid-5.aniso.svg>...
   Writing solution to <sol-5.aniso.vtu>...
@endcode
*
* 这个文本输出显示了各向异性细化的连续应用所带来的单元数量的减少。在最后一步细化之后，节省的数量已经积累到几乎是各向同性情况下所需单元和自由度的四倍。装配所需的时间也以类似的因素增加。
* 第一个有趣的部分当然是看网格的样子。左边是各向同性的细化网格，右边是各向异性的网格（颜色表示单元的细化程度）。
*  <table width="80%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-0.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-0.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-1.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-1.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-2.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-2.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-3.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-3.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-4.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-4.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-5.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.grid-5.aniso.9.2.png" alt="">
    </td>
  </tr>
</table> 
*

* 另一个有趣的事情是，当然，要看到两个网格序列上的解决方案。在这里，它们是在细化周期1和4上，清楚地显示出解决方案确实是由<i>discontinuous</i>片断多项式组成的。
*  <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.sol-1.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.sol-1.aniso.9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.sol-4.iso.9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-30.sol-4.aniso.9.2.png" alt="">
    </td>
  </tr>
</table> 
* 我们看到，各向异性细化网格上的解与各向同性细化网格上的解非常相似。因此，各向异性指标似乎可以有效地选择适当的单元进行各向异性的细化。
* 图片也解释了为什么网格被细化成这样。在整个域的左边部分，细化只沿着单元的 $y$ 轴进行。在域的右边部分，细化是以各向异性细化为主，因为解决方案的各向异性特征
*
*-从一到零的跳跃
*
* - 与平流方向转弯的网格不太一致。然而，在四分之一圆的底部和最接近观察者的部分，这个跳跃又变得越来越与网格对齐，细化算法的反应是创建各向异性的单元，也就是增加长宽比。
* 似乎各向异性特征和粗大网格的必要对齐会大大降低现实世界问题的性能。这在一般情况下是不会错的。例如，如果将各向异性细化应用于出现冲击的问题（例如，在 step-69 中求解的方程），那么在许多情况下，各向异性细化的性能会大大降低。)，那么在很多情况下，冲击并不与网格对齐，各向异性的细化对其帮助不大，除非同时引入技术将网格与冲击对齐。在这些情况下，网格已经与各向异性的特征对齐，因为它当然是与边界本身对齐的，各向异性的细化几乎总是能提高这些情况下适应网格的计算效率。
*

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-30.cc" 。
* */


