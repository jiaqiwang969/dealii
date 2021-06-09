/**
@page step_1 The step-1 tutorial program
@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Aboutthetutorial"> About the tutorial </a>
        <li><a href="#Videolecturesontutorialprograms"> Video lectures on tutorial programs </a>
        <li><a href="#Whatthisprogramdoes"> What this program does </a>
        <li><a href="#Aboutscientificcomputingingeneral"> About scientific computing in general </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Creatingthefirstmesh">Creating the first mesh</a>
        <li><a href="#Creatingthesecondmesh">Creating the second mesh</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Differentadaptiverefinementstrategies"> Different adaptive refinement strategies </a>
        <li><a href="#Differentgeometries"> Different geometries </a>
        <li><a href="#Commentsaboutprogramminganddebugging"> Comments about programming and debugging </a>
        <li><a href="#Moreaboutgraphicaloutput"> More about graphical output </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-1/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Aboutthetutorial"></a><h3> About the tutorial </h3>


由于这是第一个教程程序，让我们首先评论一下这个教程和deal.II的其他文档应该如何工作。deal.II的文档基本上有三个不同的层次。

- 该教程。这是一个程序集，展示了deal.II在实践中的应用。它通常不在单个参数的层面上讨论单个函数，而是希望给出事物如何共同工作的大画面。换句话说，它讨论的是 "概念"：什么是deal.II的构件，它们如何在有限元程序中一起使用。

- 该手册。这是deal.II中每一个类和每一个（成员）函数的文档。例如，如果你点击本页顶部的 "主页 "或 "类 "选项卡，你就可以看到。在这里你可以查到 Triangulation::create_triangulation_compatibility 的第二个参数是什么意思，这只是一个略显晦涩的例子。当你知道你要做什么，但忘记了函数到底是怎么命名的，它的参数是什么，或者它的返回值是什么时，你就需要这种级别的文档。请注意，当你读完教程并点击任何一个类或函数名称时，你也会进入手册，也就是说，当你需要对某个函数或类进行更详细的描述时，教程包含了大量进入手册的链接。另一方面，手册并不是学习deal.II的好地方，因为它只给你一个微观的观点，而没有告诉你一个函数是如何融入大局的。

- 模块。这些是一起工作或具有相关功能的类和函数组。如果你点击本页顶部的 "模块 "标签，你就会进入一个列出许多此类组的页面。每个模块都讨论了这些类的基本原理；例如， @ref Sparsity 模块讨论了与存储矩阵的稀疏模式有关的各种不同问题。这就是中级水平的文档：它们给你一个特定领域的概述。例如，当你想知道存在哪些有限元类时，你会看一下 @ref fe 模块。当然，这些模块也与手册（有时也与教程）有交叉链接；如果你点击一个类的名字，比如说三角法，如果你想了解更多关于这个类的背景，也会在类名的右上方得到一个指向这个类所属模块的链接。

让我们回到教程中来，因为你正在看的是它的第一个程序（或 "步骤"）。每个教程的程序都被细分为以下几个部分。<ol>  <li>  <b>Introduction:</b> 这是讨论程序的作用，包括数学模型，以及与以前的教程程序相比有哪些新的编程技术。     <li>  <b>The commented program:</b> 广泛地记录了源代码的清单。在这里，我们经常记录个别行或代码块，并讨论它们做什么，如何做，以及为什么。评论中经常提到介绍，也就是说，你必须先了解<i>what</i>程序想要达到的目标（介绍中讨论的目标），然后才能了解<i>how</i>它想要达到的目标。     <li>  <b>Results:</b> 程序的输出，包括注释和解释。这一部分也经常有一个小节，给出如何在不同方向上扩展程序的建议；在早期的程序中，这是为了给你提供小实验的方向，旨在使你熟悉deal.II，而在后来的程序中，更多的是关于如何使用更高级的数值技术。     <li>  <b>The plain program:</b> 剥去所有注释的源代码。如果你想看到代码的 "全貌"，这很有用，因为程序的注释版本中间有很多文字，往往很难在屏幕上一次看到单个函数的全部代码。   </ol> 

教程不仅意味着是静态文档，而且你应该玩玩它们。为此，进入 <code>example/step-1</code> 目录（或任何你感兴趣的教程的编号），然后输入

@code
  cmake .
  make
  make run
@endcode

第一条命令设置了描述本教程程序所依赖的包含文件、如何编译以及如何运行的文件。这条命令应该也能找到已安装的deal.II库，这些库是在你按照<a href="../../readme.html" target="body">README</a>文件中描述的方式编译和安装一切时产生的。如果这个命令不能找到deal.II库，那么你需要用命令提供安装的路径

@code
  cmake -DDEAL_II_DIR=/path/to/installed/deal.II .
@endcode

而不是。

上述命令中的第二条将源代码编译成可执行文件，而最后一条则是执行它（严格来说，如果可执行文件还不存在， <code>make run</code> 也会编译代码，所以如果你想的话，你可以跳过第二条命令）。这就是运行代码和产生输出所需的全部内容，在教程程序的 "结果 "部分讨论。这个顺序需要在你想玩的所有教程目录中重复。

当学习这个库时，你需要玩玩它，看看会发生什么。为此，用你喜欢的编辑器打开 <code>example/step-1/step-1.cc</code> 的源文件，并以某种方式进行修改，保存它并按上述方式运行它。在这个程序的结果部分的末尾给出了一些可能的修改建议，在那里我们还提供了一些其他有用信息的链接。




<a name="Videolecturesontutorialprograms"></a><h3> Video lectures on tutorial programs </h3>


在关于deal.II和计算科学的<a
href="http://www.math.colostate.edu/~bangerth/videos.html">Wolfgang
Bangerth video lectures</a>中也讨论和演示了这个和其他几个教程程序。特别是，你可以看到他为运行这个和其他程序所执行的步骤，你会对可以用来处理deal.II的工具有一个更好的了解。特别是，第2和第4讲概述了deal.II和任何有限元代码的构建模块。(  @dealiiVideoLectureSeeAlso{2,4}) 

如果你还不熟悉使用Linux和在命令行上运行东西，你可能会有兴趣观看讲座2.9和2.91。(  @dealiiVideoLectureSeeAlso{2.9,2.91}) These give overviews over the command 行和关于编译程序时发生的事情，分别。

请注意，deal.II正在积极开发，在开发过程中，我们偶尔会对这些视频讲座中仍然引用的函数或类进行重新命名或废弃。  例如，视频讲座5中的步骤1代码使用了一个HyperShellBoundary类，后来被SphericalManifold类取代。此外，从deal.II 9.0版本开始， GridGenerator::hyper_shell() 现在自动将SphericalManifold附加到Triangulation上。否则，讲座材料的其余部分都是相关的。

<a name="Whatthisprogramdoes"></a><h3> What this program does </h3>


让我们回到步骤1，即当前的程序。在这第一个例子中，我们实际上并没有做很多事情，而是展示了两种技术：生成三角形对象的语法是什么，以及所有单元格上简单循环的一些元素。我们创建了两个网格，一个是有规律地细化的正方形（不是很刺激，但对于一些问题来说是常见的起始网格），还有一个是更多的几何尝试：一个环形域，向内边缘细化。通过这些，你将了解到每一个有限元程序都必须有的三样东西。一个用于网格的Triangulation类型的对象；对GridGenerator函数的调用以生成网格；以及涉及迭代器的所有单元的循环（迭代器是指针的泛化，在C++标准库中经常使用；在deal.II的背景下， @ref Iterators 模块谈到了它们）。

该程序在其他方面足够小，不需要大量的介绍。

 @dealiiVideoLecture{5,6} 




<a name="Aboutscientificcomputingingeneral"></a><h3> About scientific computing in general </h3>


如果你正在阅读这个教程程序，很可能你有兴趣继续使用deal.II来完成你自己的项目。因此，你即将开始一个使用大规模科学计算库的编程练习。除非你已经是大规模编程方法的资深用户，否则这对你来说可能是一个新的领域；伴随着所有的新规则，比如你将不得不处理别人写的代码，你可能不得不考虑记录自己的代码，因为你可能在一年后不记得它到底在做什么（或者因为别人也会使用它），或者想出一些方法来测试你的程序是否在做正确的事情。这些都不是我们通常训练数学家、工程师或科学家的东西，但当你开始编写超过几百行的软件时，这些就很重要了。请记住。制作软件并不等同于仅仅编写代码。

为了使你在这一旅程中生活得更轻松，让我们指出一些资源，这些资源在你开始任何大规模的编程之前是值得浏览的。

- 在<a href="https://github.com/dealii/dealii/wiki/Frequently-Asked-Questions">
  deal.II FAQ</a>中，有大量关于deal.II特定方面问题的答案，但也有一些更普遍的问题，如 "我如何调试科学计算代码？"或 "我能否训练自己写出错误更少的代码？"。

- 你将从成为一个更好的程序员中受益。为此，一个很好的资源是Steve McConnell的[Code Complete](https://en.wikipedia.org/wiki/Code_Complete)  @cite CodeComplete  。这本书已经有几年的历史了，最后一版是在2004年出版的，但它作为良好的编程实践指南的吸引力丝毫不减，一些主要的开发者把它作为他们研究小组每一代成员的集体阅读项目。

- <a href="http://software-carpentry.org/">Software Carpentry project</a>，提供了处理软件的许多重要主题的介绍，如版本控制、make文件、测试等。它是专门为科学家和工程师编写的，而不是为计算机科学家编写的，并以简短的、实用的课程为重点。

- <a href="https://bssw.io/">Better Scientific Software
  project</a>有大量的资源（和有趣的博文），涵盖了编写科学软件的许多方面。

- <a href="https://ideas-productivity.org/">IDEAS
  project</a>也有关于软件开发的资源，特别是用于并行计算。在该网站的 "活动 "部分有录制的教程和网络研讨会，涉及许多有趣的主题。

- 一篇关于<a href="http://arxiv.org/abs/1210.0530">Best
  Practices for Scientific Computing</a>的文章，介绍了许多方法，通过这些方法，你可以确保你是一个高效的程序员，写出的程序可以正常工作。

作为一个一般性建议。如果你期望在未来花几天时间来编写软件，请帮你自己的忙，学习能够使你的生活更有效率的工具，特别是调试器和集成开发环境。(  @dealiiVideoLectureSeeAlso{7,8,8.01,25})  你会发现，通过提高工作效率，你很快就会把学习这些工具的时间拿回来几倍!上面提到的几个视频讲座展示了如何使用集成开发环境或调试器等工具。


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
 * The most fundamental class in the library is the Triangulation class, which
 * is declared here:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * @endcode
 * 
 * Here are some functions to generate standard grids:
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h>
 * @endcode
 * 
 * Output of grids in various graphics formats:
 * 
 * @code
 * #include <deal.II/grid/grid_out.h>
 * 
 * @endcode
 * 
 * This is needed for C++ output:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * @endcode
 * 
 * And this for the declarations of the `std::sqrt` and `std::fabs` functions:
 * 
 * @code
 * #include <cmath>
 * 
 * @endcode
 * 
 * The final step in importing deal.II is this: All deal.II functions and
 * classes are in a namespace <code>dealii</code>, to make sure they don't
 * clash with symbols from other libraries you may want to use in conjunction
 * with deal.II. One could use these functions and classes by prefixing every
 * use of these names by <code>dealii::</code>, but that would quickly become
 * cumbersome and annoying. Rather, we simply import the entire deal.II
 * namespace for general use:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Creatingthefirstmesh"></a> 
 * <h3>Creating the first mesh</h3>
 * 

 * 
 * In the following, first function, we simply use the unit square as domain
 * and produce a globally refined grid from it.
 * 
 * @code
 * void first_grid()
 * {
 * @endcode
 * 
 * The first thing to do is to define an object for a triangulation of a
 * two-dimensional domain:
 * 
 * @code
 *   Triangulation<2> triangulation;
 * @endcode
 * 
 * Here and in many following cases, the string "<2>" after a class name
 * indicates that this is an object that shall work in two space
 * dimensions. Likewise, there are versions of the triangulation class that
 * are working in one ("<1>") and three ("<3>") space dimensions. The way
 * this works is through some template magic that we will investigate in
 * some more detail in later example programs; there, we will also see how
 * to write programs in an essentially dimension independent way.
 * 

 * 
 * Next, we want to fill the triangulation with a single cell for a square
 * domain. The triangulation is the refined four times, to yield $4^4=256$
 * cells in total:
 * 
 * @code
 *   GridGenerator::hyper_cube(triangulation);
 *   triangulation.refine_global(4);
 * 
 * @endcode
 * 
 * Now we want to write a graphical representation of the mesh to an output
 * file. The GridOut class of deal.II can do that in a number of different
 * output formats; here, we choose scalable vector graphics (SVG) format
 * that you can visualize using the web browser of your choice:
 * 
 * @code
 *   std::ofstream out("grid-1.svg");
 *   GridOut       grid_out;
 *   grid_out.write_svg(triangulation, out);
 *   std::cout << "Grid written to grid-1.svg" << std::endl;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Creatingthesecondmesh"></a> 
 * <h3>Creating the second mesh</h3>
 * 

 * 
 * The grid in the following, second function is slightly more complicated in
 * that we use a ring domain and refine the result once globally.
 * 
 * @code
 * void second_grid()
 * {
 * @endcode
 * 
 * We start again by defining an object for a triangulation of a
 * two-dimensional domain:
 * 
 * @code
 *   Triangulation<2> triangulation;
 * 
 * @endcode
 * 
 * We then fill it with a ring domain. The center of the ring shall be the
 * point (1,0), and inner and outer radius shall be 0.5 and 1. The number of
 * circumferential cells could be adjusted automatically by this function,
 * but we choose to set it explicitly to 10 as the last argument:
 * 
 * @code
 *   const Point<2> center(1, 0);
 *   const double   inner_radius = 0.5, outer_radius = 1.0;
 *   GridGenerator::hyper_shell(
 *     triangulation, center, inner_radius, outer_radius, 10);
 * @endcode
 * 
 * By default, the triangulation assumes that all boundaries are straight
 * lines, and all cells are bi-linear quads or tri-linear hexes, and that
 * they are defined by the cells of the coarse grid (which we just
 * created). Unless we do something special, when new points need to be
 * introduced the domain is assumed to be delineated by the straight
 * lines of the coarse mesh, and new points will simply be in the middle
 * of the surrounding ones. Here, however, we know that the domain is
 * curved, and we would like to have the Triangulation place new points
 * according to the underlying geometry. Fortunately, some good soul
 * implemented an object which describes a spherical domain, of which the
 * ring is a section; it only needs the center of the ring and
 * automatically figures out how to instruct the Triangulation where to
 * place the new points. The way this works in deal.II is that you tag
 * parts of the triangulation you want to be curved with a number that is
 * usually referred to as "manifold indicator" and then tell the
 * triangulation to use a particular "manifold object" for all places
 * with this manifold indicator. How exactly this works is not important
 * at this point (you can read up on it in step-53 and @ref manifold).
 * The functions in GridGenerator handle this for us in most
 * circumstances: they attach the correct manifold to a domain so that
 * when the triangulation is refined new cells are placed in the correct
 * places. In the present case GridGenerator::hyper_shell attaches a
 * SphericalManifold to all cells: this causes cells to be refined with
 * calculations in spherical coordinates (so new cells have edges that
 * are either radial or lie along concentric circles around the origin).
 *   

 * 
 * By default (i.e., for a Triangulation created by hand or without a
 * call to a GridGenerator function like GridGenerator::hyper_shell or
 * GridGenerator::hyper_ball), all cells and faces of the Triangulation
 * have their manifold_id set to numbers::flat_manifold_id, which is
 * the default if you want a manifold that produces straight edges, but
 * you can change this number for individual cells and faces. In that
 * case, the curved manifold thus associated with number zero will not
 * apply to those parts with a non-zero manifold indicator, but other
 * manifold description objects can be associated with those non-zero
 * indicators. If no manifold description is associated with a particular
 * manifold indicator, a manifold that produces straight edges is
 * implied. (Manifold indicators are a slightly complicated topic; if
 * you're confused about what exactly is happening here, you may want to
 * look at the @ref GlossManifoldIndicator "glossary entry on this
 * topic".) Since the default chosen by GridGenerator::hyper_shell is
 * reasonable we leave things alone.
 *   

 * 
 * In order to demonstrate how to write a loop over all cells, we will
 * refine the grid in five steps towards the inner circle of the domain:
 * 
 * @code
 *   for (unsigned int step = 0; step < 5; ++step)
 *     {
 * @endcode
 * 
 * Next, we need to loop over the active cells of the triangulation. You
 * can think of a triangulation as a collection of cells. If it were an
 * array, you would just get a pointer that you increment from one
 * element to the next using the operator `++`. The cells of a
 * triangulation aren't stored as a simple array, but the concept of an
 * <i>iterator</i> generalizes how pointers work to arbitrary collections
 * of objects (see <a href=
 * "http://en.wikipedia.org/wiki/Iterator#C.2B.2B">wikipedia</a> for more
 * information). Typically, any container type in C++ will return an
 * iterator pointing to the start of the collection with a method called
 * `begin`, and an iterator point to 1 past the end of the collection with
 * a method called `end`. We can increment an iterator `it` with the
 * operator `++it`, dereference it to get the underlying data with `*it`,
 * and check to see if we're done by comparing `it != collection.end()`.
 *       

 * 
 * The second important piece is that we only need the active cells.
 * Active cells are those that are not further refined, and the only
 * ones that can be marked for further refinement. deal.II provides
 * iterator categories that allow us to iterate over <i>all</i> cells
 * (including the parent cells of active ones) or only over the active
 * cells. Because we want the latter, we need to call the method
 * Triangulation::active_cell_iterators().
 *       

 * 
 * Putting all of this together, we can loop over all the active cells of
 * a triangulation with
 * <div class=CodeFragmentInTutorialComment>
 * @code{.cpp}
 *     for (auto it = triangulation.active_cell_iterators().begin();
 *          it != triangulation.active_cell_iterators().end();
 *          ++it)
 *       {
 *         auto cell = *it;
 *         // Then a miracle occurs...
 *       }
 * @endcode
 * </div>
 * In the initializer of this loop, we've used the `auto` keyword for the
 * type of the iterator `it`. The `auto` keyword means that the type of
 * the object being declared will be inferred from the context. This
 * keyword is useful when the actual type names are long or possibly even
 * redundant. If you're unsure of what the type is and want to look up
 * what operations the result supports, you can go to the documentation
 * for the method Triangulation::active_cell_iterators(). In this case,
 * the type of `it` is `Triangulation::active_cell_iterator`.
 *       

 * 
 * While the `auto` keyword can save us from having to type out long names
 * of data types, we still have to type a lot of redundant declarations
 * about the start and end iterator and how to increment it. Instead of
 * doing that, we'll use
 * <a href="http://en.cppreference.com/w/cpp/language/range-for">range-
 * based for loops</a>, which wrap up all of the syntax shown above into a
 * much shorter form:
 * 
 * @code
 *       for (auto &cell : triangulation.active_cell_iterators())
 *         {
 * @endcode
 * 
 * @note See @ref Iterators for more information about the iterator
 * classes used in deal.II, and @ref CPP11 for more information about
 * range-based for loops and the `auto` keyword.
 *           

 * 
 * Next, we loop over all vertices of the cells. For that purpose
 * we query an iterator over the vertex indices (in 2d, this is an
 * array that contains the elements `{0,1,2,3}`, but since
 * `cell->vertex_indices()` knows the dimension the cell lives in, the
 * array so returned is correct in all dimensions and this enables
 * this code to be correct whether we run it in 2d or 3d, i.e., it
 * enables "dimension-independent programming" -- a big part of what
 * we will discuss in step-4).
 * 
 * @code
 *           for (const auto v : cell->vertex_indices())
 *             {
 * @endcode
 * 
 * If this cell is at the inner boundary, then at least one of its
 * vertices must sit on the inner ring and therefore have a radial
 * distance from the center of exactly 0.5, up to floating point
 * accuracy. So we compute this distance, and if we find a vertex
 * with this property, we flag this cell for later refinement. We
 * can then also break the loop over all vertices and move on to
 * the next cell.
 *               

 * 
 * Because the distance from the center is computed as a floating
 * point number, we have to expect that whatever we compute is
 * only accurate to within
 * [round-off](https://en.wikipedia.org/wiki/Round-off_error). As
 * a consequence, we can never expect to compare the distance
 * with the inner radius by equality: A statement such as
 * `if (distance_from_center == inner_radius)` will fail
 * unless we get exceptionally lucky. Rather, we need to do this
 * comparison with a certain tolerance, and the usual way to do
 * this is to write it as `if (std::abs(distance_from_center -
 * inner_radius) <= tolerance)`
 * where `tolerance` is some small number larger
 * than round-off. The question is how to choose it: We could just
 * pick, say, `1e-10`, but this is only appropriate if the objects
 * we compare are of size one. If we had created a mesh with cells
 * of size `1e+10`, then `1e-10` would be far lower than round-off
 * and, as before, the comparison will only succeed if we get
 * exceptionally lucky. Rather, it is almost always useful to make
 * the tolerance *relative* to a typical "scale" of the objects
 * being compared. Here, the "scale" would be the inner radius, or
 * maybe the diameter of cells. We choose the former and set the
 * tolerance equal to $10^{-6}$ times the inner radius of the
 * annulus.
 * 
 * @code
 *               const double distance_from_center =
 *                 center.distance(cell->vertex(v));
 * 
 *               if (std::fabs(distance_from_center - inner_radius) <=
 *                   1e-6 * inner_radius)
 *                 {
 *                   cell->set_refine_flag();
 *                   break;
 *                 }
 *             }
 *         }
 * 
 * @endcode
 * 
 * Now that we have marked all the cells that we want refined, we let
 * the triangulation actually do this refinement. The function that does
 * so owes its long name to the fact that one can also mark cells for
 * coarsening, and the function does coarsening and refinement all at
 * once:
 * 
 * @code
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * 
 * 
 * @endcode
 * 
 * Finally, after these five iterations of refinement, we want to again
 * write the resulting mesh to a file, again in SVG format. This works just
 * as above:
 * 
 * @code
 *   std::ofstream out("grid-2.svg");
 *   GridOut       grid_out;
 *   grid_out.write_svg(triangulation, out);
 * 
 *   std::cout << "Grid written to grid-2.svg" << std::endl;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * Finally, the main function. There isn't much to do here, only to call the
 * two subfunctions, which produce the two grids.
 * 
 * @code
 * int main()
 * {
 *   first_grid();
 *   second_grid();
 * }
 * @endcode
examples/step-1/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序会产生两个网格的图形（grid-1.svg和grid-2.svg）。你可以用大多数网络浏览器打开它们--在最简单的情况下，只要在文件系统资源管理器中打开当前目录，然后点击文件。如果你喜欢在命令行上工作，你可以用该文件调用你的网络浏览器。`firefox grid-1.svg`，`google-chrome grid-1.svg`，或者任何你的浏览器的名字。如果你这样做，这两个网格应该是这样的。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-1-r9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-2-r9.2.png" alt="">
    </td>
  </tr>
</table> 

左边那个，嗯，不是很刺激。右边的是&mdash；至少是&mdash；非传统的。这些图片对每个单元的 "细化水平 "进行了颜色编码。一个粗略的网格单元要被细分多少次才能得到给定的单元。在左图中，这是无聊的，因为网格被全局细化了若干次，也就是说，<i>every</i>单元被细化的次数相同。

(虽然第二个网状结构完全是人为捏造的，当然在应用中也不太实用，但令大家惊讶的是，它已经进入了文献：见  @cite Mu05  。显然，它至少对某些事情是有好处的）。)




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="Differentadaptiverefinementstrategies"></a><h4> Different adaptive refinement strategies </h4>


这个程序显然没有太多的功能，但特别是 <code>second_grid</code> 函数有一堆你可以玩弄它的地方。例如，你可以修改我们决定细化哪些单元格的标准。一个例子是把条件改成这样。

@code
      for (auto &cell: triangulation.active_cell_iterators())
        if (cell->center()[1] > 0)
          cell->set_refine_flag ();
@endcode

这将细化所有单元中心的 $y$ 坐标大于零的单元（我们通过解除引用 <code>cell</code> 迭代器调用的 <code>TriaAccessor::center</code> 函数返回一个Point<2>对象；下标 <code>[0]</code> 将得到 $x$ 坐标，下标 <code>[1]</code> 得到 $y$  坐标）。通过查看TriaAccessor提供的函数，你也可以使用更复杂的标准进行细化。

一般来说，你能用`cell->something()`形式的操作做什么，在文档中有点困难，因为`cell`不是一个指针，而是一个迭代器。你可以在单元格上调用的函数可以在`TriaAccessor'类的文档中找到（它的函数也可以在单元格的面或更普遍的、出现在三角形中的各种几何对象上调用），以及`CellAccessor'（它增加了一些专门针对*单元格的函数）。

对整个迭代器概念的更彻底的描述可以在 @ref Iterators 文档模块中找到。




<a name="Differentgeometries"></a><h4> Different geometries </h4>


另一种可能性是生成完全不同几何形状的网格。虽然对于复杂的几何体来说，使用网格生成器获得的网格是没有办法的，但是有大量的几何体，deal.II可以使用GridGenerator命名空间的函数来创建网格。许多这样的几何体（如本例程序中使用的几何体）包含有弯曲面的单元：换句话说，我们希望放置在边界上的新顶点位于一个圆上。deal.II通过Manifold类（以及从它继承的类）处理复杂的几何体；尤其是GridGenerator中对应于非笛卡尔网格的函数（如 GridGenerator::hyper_shell 或 GridGenerator::truncated_cone) 将一个Manifold对象附加到三角网格中应该是曲线的部分（分别为SphericalManifold和CylindricalManifold），并在应该是平面的部分使用另一个Manifold（FlatManifold）。关于这些类的设计理念和接口的描述，请参见Manifold的文档或 @ref manifold "manifold模块"。看看它们提供了什么，看看如何在这样的程序中使用它们。

我们还在第49步中讨论了其他各种创建和操作网格的方法（并描述了附加Manifolds的过程）。




<a name="Commentsaboutprogramminganddebugging"></a><h4> Comments about programming and debugging </h4>


最后，我们对用deal.II修改或编写程序做一个总体的评论。当你开始使用教程程序或你自己的应用程序时，你会发现错误会发生：你的程序会包含一些代码，这些代码要么是立即中止程序，要么是一些简单地导致错误结果的bug。无论哪种情况，你都会发现知道如何使用调试器是非常有帮助的：你可能会通过把调试输出放到你的程序中，编译它，然后运行它来应付一段时间，但最终用调试器寻找错误会更快，更方便，更可靠，因为你不必总是重新编译程序，而且你可以检查变量的值和它们的变化。

与其推迟学习如何使用调试器，直到你真的看不到任何其他方法来发现一个错误，这里是我们将在这个项目中提供的一个建议：尽快学习如何使用调试器。这将是很好的时间投资。( 从顶层的<a
href="http://www.dealii.org/">deal.II webpage</a>链接到的 @dealiiVideoLectureSeeAlso{25}) The deal.II Frequently Asked 问题（FAQ）页面也提供了大量关于调试deal.II程序的提示。




<a name="Moreaboutgraphicaloutput"></a><h4> More about graphical output </h4>


在你的论文或出版物中包含网格往往是有用的。为此，按细化级别对单元格进行颜色编码，并在每个单元格上打印单元格号，可能不是很有用。但这并不意味着一定要这样做 -- GridOut类允许为每种可能的输出格式设置标志（参见GridOutFlags命名空间中的类），以控制网格的具体绘制方式。当然，你也可以选择其他的输出文件格式，如VTK或VTU；这对三维网格特别有用，因为二维格式如SVG并不特别有用，因为它固定了三维物体的特定视角。因此，你可能想探索GridOut类中的其他选项。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-1.cc"
*/
