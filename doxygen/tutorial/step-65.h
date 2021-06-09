/**
@page step_65 The step-65 tutorial program
This tutorial depends on step-49.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Workingwithmanifolds">Working with manifolds</a>
      <ul>
        <li><a href="#Whatwewant">What we want</a>
      </ul>
        <li><a href="#TheclassTransfiniteInterpolationManifold">The class TransfiniteInterpolationManifold</a>
      <ul>
        <li><a href="#Howitworks">How it works</a>
      </ul>
        <li><a href="#Transfiniteinterpolationisexpensiveandhowtodealwithit">Transfinite interpolation is expensive and how to deal with it</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Analyticalsolutionandcoefficient">Analytical solution and coefficient</a>
        <li><a href="#ThePoissonProblemclass">The PoissonProblem class</a>
        <li><a href="#Gridcreationandinitializationofthemanifolds">Grid creation and initialization of the manifolds</a>
        <li><a href="#Setupofdatastructures">Setup of data structures</a>
        <li><a href="#Assemblyofthesystemmatrixandrighthandside">Assembly of the system matrix and right hand side</a>
        <li><a href="#Solutionofthelinearsystem">Solution of the linear system</a>
        <li><a href="#Outputofthesolutionandcomputationoferrors">Output of the solution and computation of errors</a>
        <li><a href="#ThePoissonProblemrunfunction">The PoissonProblem::run() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-65/doc/intro.dox



 <br> 

<i>
This program was contributed by Martin Kronbichler.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序介绍了一个高级流形类--TransfiniteInterpolationManifold，以及如何解决其主要缺点--相对较高的成本。

<a name="Workingwithmanifolds"></a><h3>Working with manifolds</h3>


<a name="Whatwewant"></a><h4>What we want</h4>


在许多应用中，有限元网格必须能够表示一个相对复杂的几何体。在step-1、step-49和step-53教程程序中，已经介绍了deal.II库中可用的一些生成网格的技术。给定一个基础网格，deal.II能够通过将单元格均匀地或仅在计算域的选定部分细分为子单元来创建一个更细的网格。除了GridGenerator命名空间中的基本网格划分功能外，deal.II还提供了一些接口，可以使用命名空间GridIn的功能读入由（仅四边形和六边形）网格生成器生成的网格，例如在步骤5中演示。外部生成的网格的一个基本限制是，网格中生成的单元所提供的信息只包括顶点的位置和它们的连接性，而没有最初创建这个网格的网格生成器所提供的底层几何体的背景。一旦网格在deal.II中被细化，需要放置更多的点，这就成了问题。第54步教程程序显示了如何通过使用OpenCASCADE库方面的CAD曲面来克服这个限制，第53步通过在源代码内以编程方式提供同类信息。

在deal.II中，网格细化过程中新点的放置或高阶映射的定义由流形对象控制，详见 @ref manifold  "流形模块"。举个例子，考虑以下二维环形的情况（图片取自流形模块）。如果我们从10个单元的初始网格开始，在不附加任何流形的情况下全局细化网格三次，我们将得到以下网格。

 @image html hypershell-nothing.png "" 

图片看起来是这样的，因为在默认情况下，deal.II只知道通过平均父单元的顶点位置来放置子单元的顶点。这就产生了一个多边形域，其面是原始（粗略的网格）单元的边缘。很明显，我们必须对三角形的边界面附加一个曲线描述，以便在网格细化时重现圆形，就像下图一样。

 @image html hypershell-boundary-only.png "" 

这样就好多了。如果我们继续细化网格，至少现在的内外边界已经接近真实的圆了。然而，这幅图中的网格对于环形来说仍然不是最佳的，因为从一个单元到下一个单元的<i>interior</i>线在某些顶点上有扭结，人们宁愿使用下面的网格。

 @image html hypershell-all.png "" 

在这最后一种（最佳）情况下，也是由 GridGenerator::hyper_shell(), 产生的默认情况，曲线流形描述（在这种情况下是极地流形描述）不仅适用于边界面，而且适用于整个域。每当三角化要求一个新的点，例如，当它将一个单元细化为四个子单元时，边缘或单元的中点，它将沿着极坐标系统中各自的中点来放置它们。相比之下，上面那种只有边界受制于极地流形的情况，只有沿边界的中点会沿着曲线描述放置，而内部的中点会通过直角坐标系中周围点的合适平均值来计算（更多细节见 @ref manifold  "流形模块"）。

在这一点上，人们可能会认为曲线形的体积描述是一种方式。这一般来说是不会错的，尽管有时并不那么容易描述这到底应该如何工作。这里有几个例子。

- 想象一下，上面的网格实际上是一个圆盘，而不只是一个环。   在这种情况下，极地流形会在原点退化，不会产生合理的新点。事实上，为那些应该 "看起来是圆的 "但可能在原点或接近原点的东西定义流形描述是令人惊讶的非常困难的。

- 当人们试图将球形流形附加到整个体积上时，类似的事情也会发生在三维球的原点；在这种情况下，新流形点的计算会因异常而中止。

- CAD几何体通常只描述域的边界，就像我们在上面第二张图中只把流形附在边界上一样。同样地，Step-54只使用CAD几何体来生成表面网格（也许是因为这是解决相关问题所需要的），但是如果想解决那里描述的水或船周围空气中的问题，我们就需要有一个体积网格。那么问题来了，我们应该如何准确地描述域的内部应该发生什么。

这些简单的例子清楚地表明，对于许多有趣的情况，我们必须从对全体积的分析性曲率描述的愿望中退一步。将需要有<i>some</i>种导致曲率也在内部的信息，但必须有可能做到这一点，而不需要实际写下描述那种几何的明确公式。

那么，如果我们不在内部做任何事情，只把表面描述为流形，会发生什么呢？有时，如上图所示的环，结果并不可怕。但有时却很糟糕。考虑一个环状体的情况（例如用 GridGenerator::torus()) 生成的环状体，只在表面附加一个TorusManifold对象，内部单元和面没有附加流形，细化前有六个单元在环状方向。如果对网格进行一次细化，我们将得到如下的网格，图中显示的是网格的上半部分被剪切掉了。

 @image html torus_no_inner_manifold.png "" 

这显然是次优的。事实上，如果我们从少于上面显示的环形方向的六个单元开始，映射在某些区域实际上是倒置的，因为沿着内部单元放置的新点与边界相交，因为它们不是沿着环形方向的圆形。环形的简单情况仍然可以被固定下来，因为我们知道环形方向是沿着圆柱坐标系的，所以在曲面上附加一个TorusManifold，结合CylindricalManifold，在环形方向上有适当的周期性，应用于所有的内部实体，将产生一个高质量的网格，如下所示，现在有两个顶部的单元被隐藏起来。

 @image html torus_cylindrical_inner_manifold.png "" 

这个网格是相当不错的，但显然它与对体积的良好描述有关，而我们在其他情况下缺乏这种描述。实际上，在这种情况下也有一个不完美之处，因为我们可以看到域的内部的两个相邻单元的一些不自然的扭结，这些扭结被顶部的两个边界单元所隐藏，与下面的设置相反（由 GridGenerator::torus() 应用的默认流形并使用TransfiniteInterpolationManifold）。

 @image html torus_transfinite_manifold.png "" 

<a name="TheclassTransfiniteInterpolationManifold"></a><h3>The class TransfiniteInterpolationManifold</h3>


为了找到一个更好的策略，让我们再看看二维圆盘（这也是沿环形方向旋转的环形实体的基础）。正如我们上面所学到的，我们只能将弯曲的极坐标描述应用于边界（或离原点足够远的单元的边缘），但最终必须过渡到对圆盘中心的直线描述。如果我们在单元的内部使用平坦的流形（也就是说，新的顶点是通过相邻的现有顶点的平均化而产生的），而极地流形只用于圆盘的边界，那么在进行四次全局细化之后，我们会得到以下网格。

 @image html circular_mesh_only_boundary_manifold.png "" 

这并不是一个可怕的网格。同时，如果你知道最初的粗略网格是由中间的一个正方形和周围的四个盖子组成的，那么就不难看出这个网格的每一个细化步骤都是为了得到上面的图片。

虽然deal.II的三角测量类在创建新点时试图将信息从边界传播到内部，但这种算法的影响范围是有限的。

 @image html circular_mesh_boundary_cells.png "" 

上图突出了盘面上那些触及边界的单元，在当时只看一个单元的情况下，原则上可以考虑边界信息。显然，随着网格的细化，可以考虑到一些曲率的区域变得更加有限，因此在网格中形成了看似不规则的点。当计算最左边图片中任何一个边界单元的中心时，理想的位置是外圈和中间的单元之间的中点。这正是三角剖分类中第一个细化步骤所使用的内容。然而，对于第二次细化，所有的内部边缘以及内部单元格层只能根据平面流形描述来加点。

在这一点上，我们意识到需要什么来创建一个更好的网格。对于<i>all</i>子单元中的<i>any</i>新点，即在最左边的红色阴影层内创建的单元，我们要计算相对于各自粗单元覆盖区域内的曲率的插值。这可以通过在上图最左边面板的粗略网格的高亮单元中添加TransfiniteInterpolationManifold类来实现。该类遵守一般流形接口，即在其定义域内给定任何一组点，它可以计算符合流形的加权平均数（使用的公式将在一分钟后给出）。这些加权平均数在细化网格的时候，或者在给定单元上评估符合该流形的高阶映射（如MappingQGeneric或MappingC1）的时候都会用到。在圆盘粗网格的阴影单元上使用该流形（即不仅仅是最外层的单元），在经过四个全局步骤的细化后产生以下网格。

 @image html circular_mesh_transfinite_interpolation.png "" 

这个网格的线条仍有一些扭结，但它们仅限于粗大的网格单元之间的面，而网格的其他部分则如人们所希望的那样平滑。事实上，给定一个直边的中心单元，这个表示法是最好的，因为所有的网格单元都遵循一个平滑的过渡，从内部的方形块的直边到边界的圆形。(我们可以做得更好一些，在中心的方形块中也允许一些曲率，随着中心的接近，这些曲率最终会消失)。




<a name="Howitworks"></a><h4>How it works</h4>


在一个有一条曲线和三条直线的圆盘的简单情况下，我们可以明确写下如何实现形状的混合。为此，将物理单元（如顶部的单元）映射回参考坐标系 $(\xi,\eta)\in (0,1)^2$ 是很有用的，在那里我们计算某些点之间的平均数。如果我们使用一个由四个顶点 $(x_0,y_0), (x_1,y_1), (x_2,y_2), (x_3, y_3)$ 跨越的简单双线性地图，一个点 $(\xi, \eta)\in (0,1)^2$ 的图像将是

@f{align*}{
(x,y) = (1-\xi)(1-\eta) (x_0,y_0) + \xi(1-\eta) (x_1,y_1) +
       (1-\xi)\eta  (x_2,y_2) + \xi\eta  (x_3,y_3).


@f}



对于曲面的情况，我们要修改这个公式。对于圆盘粗略网格的顶部单元，我们可以假设点 $(x_0,y_0)$ 和 $(x_1,y_1)$ 位于下端的直线上，点 $(x_2,y_2)$ 和 $(x_3,y_3)$ 沿顶部由一个四分之一圆连接。然后我们将点 $(\xi, \eta)$ 映射为

@f{align*}{
(x,y) = (1-\eta) \big[(1-\xi) (x_0,y_0) + \xi (x_1,y_1)\big] +
      \eta \mathbf{c}_3(\xi),


@f}

其中 $\mathbf{c}_3(\xi)$ 是一条曲线，用一个arclength参数 $\xi\in (0,1)$ 来描述四分之一圆的 $(x,y)$ 坐标。这代表了单元格的直线下边缘和曲线上边缘之间的线性插值，也是上图的基础。

这个公式很容易被推广到所有四条边都由曲线而不是直线描述的情况。我们称这四个函数为 $\xi$ 或 $\eta$ ，在水平和垂直方向上，分别为四边形的左、右、下和上边缘的参数。然后，内插法为

@f{align*}{
(x,y) =& (1-\xi)\mathbf{c}_0(\eta) + \xi \mathbf{c}_1(\eta)
        +(1-\eta)\mathbf{c}_2(\xi) + \eta \mathbf{c}_3(\xi)\\
       &-\big[(1-\xi)(1-\eta) (x_0,y_0) + \xi(1-\eta) (x_1,y_1) +
        (1-\xi)\eta  (x_2,y_2) + \xi\eta  (x_3,y_3)\big].


@f}



这个公式假设边界曲线与顶点 $(x_0,y_0), (x_1,y_1), (x_2,y_2), (x_3, y_3)$ ，例如 $\mathbf{c}_0(0)
= (x_0,y_0)$ 或 $\mathbf{c}_0(1) = (x_2,y_2)$ 匹配和重合。公式第二行的双线性插值的减法确保了在边界上完全遵循规定的曲线：沿着四条边中的每一条，我们需要减去在角上评估的两条相邻边的贡献，这时的贡献只是一个顶点位置。很容易检查出，如果四条曲线中的三条 $\mathbf{c}_i$ 是直的，从而与双线性插值重合，那么上面的圆的公式就再现了。

这个公式被称为转折内插，由<a
href="https://doi.org/10.1002%2Fnme.1620070405">Gordon and Hall</a>在1973年提出。尽管转义插值本质上只表示边界曲线的线性混合，但插值完全遵循每个实数 $\xi\in (0,1)$ 或 $\eta\in (0,1)$ 的边界曲线，也就是说，它插值的点数是无限的，这也是Gordon和Hall将这种插值的变体称为转义插值的最初动机。另一种解释是，无限插值从左右和上下线性插值，我们需要从中减去双线性插值，以确保在域的内部有一个单位重量。

无限插值很容易被推广到三个空间维度。在这种情况下，插值允许为一个三维单元的任何一个四边形混合6种不同的表面描述，为一个单元的线混合12种边缘描述。同样，为了确保地图的一致性，有必要减去边缘的贡献，再加上顶点的贡献，使曲线遵循规定的表面或边缘描述。在三维的情况下，也可以使用从曲线边缘到相邻面和相邻单元的转折插值。

在处理.II中的转置插值是通用的，因为它可以处理任意的曲线。它将以 $d$ -维空间的原始坐标来评估曲线，但有一个（或两个，在三维的边缘情况下）坐标固定在 $0$ 或 $1$ ，以确保任何其他流形类，包括CAD文件（如果需要），都可以应用于开箱。无限插值是网格生成器的标准成分，因此在deal.II库中集成这一功能的主要优势是在自适应细化和粗化网格时启用，以及用于创建更高程度的映射，使用流形来插入网格顶点以外的额外点。

关于无限插值的最后一句话，我们提到在没有体积流形描述的情况下，deal.II中的网格细化策略也是基于无限插值的权重，在这个意义上是最优的。不同的是，默认算法一次只看到一个单元，因此将只对那些接触曲面流形的单元应用最优算法。相反，在整个<i>patches</i>个单元上使用转义映射（源于一个较粗的单元），可以用转义插值法，将信息从边界传播到远处的单元。




<a name="Transfiniteinterpolationisexpensiveandhowtodealwithit"></a><h3>Transfinite interpolation is expensive and how to deal with it</h3>


一个具有无限流形描述的网格通常分两步建立。第一步是创建一个粗略的网格（或者从文件中读入），并在一些网格实体上附加一个曲面流形。对于上面的圆盘例子，我们将极坐标流形附加到沿外圆的面（这是由 GridGenerator::hyper_ball()). 自动完成的）在我们开始细化网格之前，我们再为网格的所有内部单元和边缘分配一个TransfiniteInterpolationManifold，当然，这需要基于我们分配给这些实体的一些流形ID（除了边界上的圆以外的一切）。我们是否也给圆盘的内部正方形分配一个TransfiniteInterpolationManifold并不重要，因为对具有直边（或3D中的平坦面）的粗单元进行Transfinite插值，只是产生具有直边（平坦面）的细分子。

之后，当网格被细化或基于此网格建立高阶映射时，单元将查询底层流形对象的新点。这个过程需要一组周围的点，例如一个二维单元的四个顶点，以及一组对这些点的权重，用于定义一个新的点。对于一个单元格的中间点，四个顶点中的每一个都将得到0.25的权重。对于无限插值流形，建立加权和的过程需要一些严肃的工作。根据结构，我们要根据周围点的参考坐标 $\xi$ 和 $\eta$ （或三维的 $\xi, \eta, \zeta$ ）来组合这些点。然而，deal.II中流形类的接口并没有得到周围点的参考坐标（因为它们不是全局存储的），而只是物理坐标。因此，无限插值流形必须做的第一步是反转映射，并在无限插值的一个粗网格单元内找到参考坐标（例如，上述磁盘网格的四个阴影粗网格单元之一）。这种反演是通过牛顿迭代（或者说，基于有限差分的牛顿方案与布罗伊登方法相结合）完成的，并根据上述公式多次查询无限期插值。这些查询中的每一次都可能依次调用一个昂贵的流形，例如球的球形描述，并且本身就很昂贵。由于deal.II的Manifold接口类只提供了一组点，所以转置插值最初甚至不知道周围的点集属于哪个粗大的网格单元，需要根据一些启发式方法在几个单元中搜索。就<a
href="https://en.wikipedia.org/wiki/Atlas_(topology)#Charts">charts</a>而言，我们可以把无限插值的实现描述为基于<a
href="https://en.wikipedia.org/wiki/Atlas_(topology)">atlas</a>的实现。三角形初始粗网格的每个单元都代表一个有自己参考空间的图表，而周围的流形提供了一种从图表空间（即参考单元）到物理空间的转换方式。粗网格单元的图表的集合是一个图集，像往常一样，在图集中查找东西时，人们做的第一件事就是找到正确的图表。

一旦找到周围点的参考坐标，参考坐标系中的新点将通过简单的加权和计算出来。最后，参考点被插入到无限插值的公式中，这就得到了所需的新点。

在许多情况下，曲面流形不仅在网格细化过程中使用，而且还用于确保计算域单元内边界的曲面表示。这对于保证复杂几何上的高阶多项式的高阶收敛性是必须的，但有时对于线性形状函数也需要一个精确的几何。这通常是通过对单元的多项式描述来实现的，如果表示曲面网格元素的多项式程度与数值解的多项式程度相同，则称为等参数概念。如果几何体的度数高于或低于解的度数，则分别称为超参数或次参数的几何体表示。在deal.II中，多项式表示的标准类是MappingQGeneric。例如，如果在三维中使用这个类的多项式度数 $4$ ，总共需要125个（即 $(4+1)^3$ ）点来进行内插。在这些点中，8个是单元的顶点，已经可以从网格中获得，但其他117个点需要由流形提供。如果使用无限插值流形，我们可以想象，通过回拉到一些尚未确定的粗略单元的参考坐标，然后在117个点中的每个点上进行后续的前推，是一个大量的工作，可能非常耗时。

更糟糕的是，许多程序的结构是这样的：对于同一个单元，映射被独立地查询了几次。它的主要用途是组装线性系统，即通过FEValues对象的`mapping`参数计算系统矩阵和右手边。然而，边界值的插值、数值误差的计算、输出的写入以及误差估计器的评估也必须涉及相同的映射，以确保对解向量的一致解释。因此，即使是解决一次的线性静止问题，也会多次评估映射的点。对于上面提到的三维案例，这意味着通过昂贵的算法多次计算每个单元的117个点。对于非线性或随时间变化的问题，情况更加紧迫，这些操作要反复进行。

由于通过无限插值进行流形描述的成本比对平流形的类似查询要高几百倍，所以只计算一次额外的点并在所有后续调用中使用这些点是有意义的。deal.II库提供的MappingQCache类正是为了这个目的。与系统矩阵所消耗的内存相比，这个缓存通常不会太大，在看这个教程程序的结果时就会明白。MappingQCache的用法很简单。一旦网格被设置好（或在细化过程中被改变），我们就调用 MappingQCache::initialize() ，将所需的三角形以及所需的映射作为参数。然后，初始化会遍历网格的所有单元，并查询给定映射的附加点。这些点会被保存为单元的标识符，以便以后当映射计算一些与单元相关的量（如参考坐标和物理坐标之间的映射的Jacobian）时，它们可以被返回。

最后，我们提到，TransfiniteInterpolationManifold也使得网格的细化变得更加昂贵。在这种情况下，MappingQCache无济于事，因为它将计算随后不能再使用的点；目前在deal.II中不存在一个更有效的机制。然而，网格细化也包含许多其他昂贵的步骤，所以与其他计算相比，它并不是一个大问题。它在每个时间步长或非线性迭代中也最多只发生一次。

<a name="Thetestcase"></a><h3>The test case</h3>


在这个教程程序中，TransfiniteInterpolationManifold与MappingQCache的结合被举例说明。这个测试案例相对简单，占用了许多典型程序中涉及的解决阶段，例如步骤6的教程程序。作为一个几何体，我们选择了TransfiniteInterpolationManifold的一个使用原型，即涉及一个球形的设置，而这个球形又被一个立方体所包围。例如，这样的设置将用于嵌入背景介质中的球形包容物，如果该包容物具有不同的材料特性，需要通过元素界面来跟踪两种材料之间的界面。这里给出了一个可视化的网格。

 <img src="https://www.dealii.org/images/steps/developer/step-65-mesh.png" alt=""> 

对于这种情况，我们要在域内的表面附加一个球形描述，并使用转折插值来平滑地切换到外部立方体的直线和球中心的立方体。

在该程序中，我们将遵循有限元程序的典型流程，从DoFHandler和稀疏模式的设置开始，组装一个线性系统来解决带有跳跃系数的泊松方程，用一个简单的迭代方法解决，用 VectorTools::integrate_difference() 计算一些数值误差以及误差估计器。我们记录了每个部分的时间，并运行了两次代码。在第一次运行中，我们将一个MappingQGeneric对象分别交给程序的每个阶段，在这些阶段中，点被反复计算。在第二次运行中，我们用MappingQCache代替。


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
 * The include files for this tutorial are essentially the same as in
 * step-6. Importantly, the TransfiniteInterpolationManifold class we
 * will be using is provided by `deal.II/grid/manifold_lib.h`.
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/manifold_lib.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q_generic.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * #include <fstream>
 * 
 * @endcode
 * 
 * The only new include file is the one for the MappingQCache class.
 * 
 * @code
 * #include <deal.II/fe/mapping_q_cache.h>
 * 
 * 
 * namespace Step65
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Analyticalsolutionandcoefficient"></a> 
 * <h3>Analytical solution and coefficient</h3>
 * 

 * 
 * In this tutorial program, we want to solve the Poisson equation
 * with a coefficient that jumps along a sphere of radius 0.5, and
 * using a constant right hand side of value $f(\mathbf{x}) = -3$. (This
 * setup is similar to step-5 and step-6, but the concrete values
 * for the coefficient and the right hand side are different.)
 * Due to the jump in the
 * coefficient, the analytical solution must have a kink where the
 * coefficient switches from one value to the other. To keep things simple,
 * we select an analytical solution that is quadratic in all components,
 * i.e., $u(x,y,z) = x^2 + y^2 + z^2$ in the ball of radius 0.5 and
 * $u(x,y,z) = 0.1(x^2 + y^2 + z^2) + 0.25-0.025$ in the outer part of the
 * domain. This analytical solution is compatible with the right hand side
 * in case the coefficient is 0.5 in the inner ball and 5 outside. It is
 * also continuous along the circle of radius 0.5.
 * 
 * @code
 *   template <int dim>
 *   class ExactSolution : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/ = 0) const override
 *     {
 *       if (p.norm_square() < 0.25)
 *         return p.norm_square();
 *       else
 *         return 0.1 * p.norm_square() + (0.25 - 0.025);
 *     }
 * 
 *     virtual Tensor<1, dim>
 *     gradient(const Point<dim> &p,
 *              const unsigned int /*component*/ = 0) const override
 *     {
 *       if (p.norm_square() < 0.25)
 *         return 2. * p;
 *       else
 *         return 0.2 * p;
 *     }
 *   };
 * 
 * 
 *   template <int dim>
 *   double coefficient(const Point<dim> &p)
 *   {
 *     if (p.norm_square() < 0.25)
 *       return 0.5;
 *     else
 *       return 5.0;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePoissonProblemclass"></a> 
 * <h3>The PoissonProblem class</h3>
 *   

 * 
 * The implementation of the Poisson problem is very similar to what
 * we used in the step-5 tutorial program. The two main differences
 * are that we pass a mapping object to the various steps in the
 * program in order to switch between two mapping representations as
 * explained in the introduction, and the `timer` object (of type
 * TimerOutput) that will be used for measuring the run times in the
 * various cases. (The concept of mapping objects was first
 * introduced in step-10 and step-11, in case you want to look up
 * the use of these classes.)
 * 
 * @code
 *   template <int dim>
 *   class PoissonProblem
 *   {
 *   public:
 *     PoissonProblem();
 *     void run();
 * 
 *   private:
 *     void create_grid();
 *     void setup_system(const Mapping<dim> &mapping);
 *     void assemble_system(const Mapping<dim> &mapping);
 *     void solve();
 *     void postprocess(const Mapping<dim> &mapping);
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 *     SparsityPattern           sparsity_pattern;
 *     SparseMatrix<double>      system_matrix;
 *     Vector<double>            solution;
 *     Vector<double>            system_rhs;
 * 
 *     TimerOutput timer;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * In the constructor, we set up the timer object to record wall times but
 * be quiet during the normal execution. We will query it for timing details
 * in the `PoissonProblem::run()` function. Furthermore, we select a
 * relatively high polynomial degree of three for the finite element in use.
 * 
 * @code
 *   template <int dim>
 *   PoissonProblem<dim>::PoissonProblem()
 *     : fe(3)
 *     , dof_handler(triangulation)
 *     , timer(std::cout, TimerOutput::never, TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Gridcreationandinitializationofthemanifolds"></a> 
 * <h3>Grid creation and initialization of the manifolds</h3>
 *   

 * 
 * The next function presents the typical usage of
 * TransfiniteInterpolationManifold. The first step is to create the desired
 * grid, which can be done by composition of two grids from
 * GridGenerator. The inner ball mesh is simple enough: We run
 * GridGenerator::hyper_cube() centered at the origin with radius 0.5 (third
 * function argument). The second mesh is more interesting and constructed
 * as follows: We want to have a mesh that is spherical in the interior but
 * flat on the outer surface. Furthermore, the mesh topology of the inner
 * ball should be compatible with the outer grid in the sense that their
 * vertices coincide so as to allow the two grid to be merged. The grid coming
 * out of GridGenerator::hyper_shell fulfills the requirements on the inner
 * side in case it is created with $2d$ coarse cells (6 coarse cells in 3D
 * which we are going to use) &ndash; this is the same number of cells as
 * there are boundary faces for the ball. For the outer surface, we use the
 * fact that the 6 faces on the surface of the shell without a manifold
 * attached would degenerate to the surface of a cube. What we are still
 * missing is the radius of the outer shell boundary. Since we desire a cube
 * of extent
 * $[-1, 1]$ and the 6-cell shell puts its 8 outer vertices at the 8
 * opposing diagonals, we must translate the points $(\pm 1, \pm 1, \pm 1)$
 * into a radius: Clearly, the radius must be $\sqrt{d}$ in $d$ dimensions,
 * i.e., $\sqrt{3}$ for the three-dimensional case we want to consider.
 *   

 * 
 * Thus, we have a plan: After creating the inner triangulation for
 * the ball and the one for the outer shell, we merge those two
 * grids but remove all manifolds that the functions in
 * GridGenerator may have set from the resulting triangulation, to
 * ensure that we have full control over manifolds. In particular,
 * we want additional points added on the boundary during refinement
 * to follow a flat manifold description. To start the process of
 * adding more appropriate manifold ids, we assign the manifold id 0
 * to all mesh entities (cells, faces, lines), which will later be
 * associated with the TransfiniteInterpolationManifold. Then, we
 * must identify the faces and lines that are along the sphere of
 * radius 0.5 and mark them with a different manifold id, so as to then
 * assign a SphericalManifold to those. We will choose the manifold
 * id of 1. Since we have thrown away all manifolds that pre-existed
 * after calling GridGenerator::hyper_ball(), we manually go through
 * the cells of the mesh and all their faces. We have found a face
 * on the sphere if all four vertices have a radius of 0.5, or, as
 * we write in the program, have $r^2-0.25 \approx 0$. Note that we call
 * `cell->face(f)->set_all_manifold_ids(1)` to set the manifold id
 * both on the faces and the surrounding lines. Furthermore, we want
 * to distinguish the cells inside the ball and outside the ball by
 * a material id for visualization, corresponding to the picture in the
 * introduction.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::create_grid()
 *   {
 *     Triangulation<dim> tria_inner;
 *     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);
 * 
 *     Triangulation<dim> tria_outer;
 *     GridGenerator::hyper_shell(
 *       tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);
 * 
 *     GridGenerator::merge_triangulations(tria_inner, tria_outer, triangulation);
 * 
 *     triangulation.reset_all_manifolds();
 *     triangulation.set_all_manifold_ids(0);
 * 
 *     for (const auto &cell : triangulation.cell_iterators())
 *       {
 *         for (const auto &face : cell->face_iterators())
 *           {
 *             bool face_at_sphere_boundary = true;
 *             for (const auto v : face->vertex_indices())
 *               {
 *                 if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12)
 *                   face_at_sphere_boundary = false;
 *               }
 *             if (face_at_sphere_boundary)
 *               face->set_all_manifold_ids(1);
 *           }
 *         if (cell->center().norm_square() < 0.25)
 *           cell->set_material_id(1);
 *         else
 *           cell->set_material_id(0);
 *       }
 * 
 * @endcode
 * 
 * With all cells, faces and lines marked appropriately, we can
 * attach the Manifold objects to those numbers. The entities with
 * manifold id 1 will get a spherical manifold, whereas the other
 * entities, which have the manifold id 0, will be assigned the
 * TransfiniteInterpolationManifold. As mentioned in the
 * introduction, we must explicitly initialize the manifold with
 * the current mesh using a call to
 * TransfiniteInterpolationManifold::initialize() in order to pick
 * up the coarse mesh cells and the manifolds attached to the
 * boundaries of those cells. We also note that the manifold
 * objects we create locally in this function are allowed to go
 * out of scope (as they do at the end of the function scope),
 * because the Triangulation object internally copies them.
 *     

 * 
 * With all manifolds attached, we will finally go about and refine the
 * mesh a few times to create a sufficiently large test case.
 * 
 * @code
 *     triangulation.set_manifold(1, SphericalManifold<dim>());
 * 
 *     TransfiniteInterpolationManifold<dim> transfinite_manifold;
 *     transfinite_manifold.initialize(triangulation);
 *     triangulation.set_manifold(0, transfinite_manifold);
 * 
 *     triangulation.refine_global(9 - 2 * dim);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Setupofdatastructures"></a> 
 * <h3>Setup of data structures</h3>
 *   

 * 
 * The following function is well-known from other tutorials in that
 * it enumerates the degrees of freedom, creates a constraint object
 * and sets up a sparse matrix for the linear system. The only thing
 * worth mentioning is the fact that the function receives a
 * reference to a mapping object that we then pass to the
 * VectorTools::interpolate_boundary_values() function to ensure
 * that our boundary values are evaluated on the high-order mesh
 * used for assembly. In the present example, it does not really
 * matter because the outer surfaces are flat, but for curved outer
 * cells this leads to more accurate approximation of the boundary
 * values.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::setup_system(const Mapping<dim> &mapping)
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     std::cout << "   Number of active cells:       "
 *               << triangulation.n_global_active_cells() << std::endl;
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 *     {
 *       TimerOutput::Scope scope(timer, "Compute constraints");
 * 
 *       constraints.clear();
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *       VectorTools::interpolate_boundary_values(
 *         mapping, dof_handler, 0, ExactSolution<dim>(), constraints);
 * 
 *       constraints.close();
 *     }
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 * 
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Assemblyofthesystemmatrixandrighthandside"></a> 
 * <h3>Assembly of the system matrix and right hand side</h3>
 *   

 * 
 * The function that assembles the linear system is also well known
 * from the previous tutorial programs. One thing to note is that we
 * set the number of quadrature points to the polynomial degree plus
 * two, not the degree plus one as in most other tutorials. This is
 * because we expect some extra accuracy as the mapping also
 * involves a degree one more than the polynomials for the solution.
 *   

 * 
 * The only somewhat unusual code in the assembly is the way we compute the
 * cell matrix. Rather than using three nested loop over the quadrature
 * point index, the row, and the column of the matrix, we first collect the
 * derivatives of the shape function, multiplied by the square root of the
 * product of the coefficient and the integration factor `JxW` in a separate
 * matrix `partial_matrix`. To compute the cell matrix, we then execute
 * `cell_matrix = partial_matrix * transpose(partial_matrix)` in the line
 * `partial_matrix.mTmult(cell_matrix, partial_matrix);`. To understand why
 * this works, we realize that the matrix-matrix multiplication performs a
 * summation over the columns of `partial_matrix`. If we denote the
 * coefficient by $a(\mathbf{x}_q)$, the entries in the temporary matrix are
 * $\sqrt{\text{det}(J) w_q a(x)} \frac{\partial \varphi_i(\boldsymbol
 * \xi_q)}{\partial x_k}$. If we take the product of the <i>i</i>th row with
 * the <i>j</i>th column of that matrix, we compute a nested sum involving
 * $\sum_q \sum_{k=1}^d \sqrt{\text{det}(J) w_q a(x)} \frac{\partial
 * \varphi_i(\boldsymbol \xi_q)}{\partial x_k} \sqrt{\text{det}(J) w_q a(x)}
 * \frac{\partial \varphi_j(\boldsymbol \xi_q)}{\partial x_k} = \sum_q
 * \sum_{k=1}^d\text{det}(J) w_q a(x)\frac{\partial \varphi_i(\boldsymbol
 * \xi_q)}{\partial x_k} \frac{\partial \varphi_j(\boldsymbol
 * \xi_q)}{\partial x_k}$, which is exactly the terms needed for the
 * bilinear form of the Laplace equation.
 *   

 * 
 * The reason for choosing this somewhat unusual scheme is due to the heavy
 * work involved in computing the cell matrix for a relatively high
 * polynomial degree in 3D. As we want to highlight the cost of the mapping
 * in this tutorial program, we better do the assembly in an optimized way
 * in order to not chase bottlenecks that have been solved by the community
 * already. Matrix-matrix multiplication is one of the best optimized
 * kernels in the HPC context, and the FullMatrix::mTmult() function will
 * call into those optimized BLAS functions. If the user has provided a good
 * BLAS library when configuring deal.II (like OpenBLAS or Intel's MKL), the
 * computation of the cell matrix will execute close to the processor's peak
 * arithmetic performance. As a side note, we mention that despite an
 * optimized matrix-matrix multiplication, the current strategy is
 * sub-optimal in terms of complexity as the work to be done is proportional
 * to $(p+1)^9$ operations for degree $p$ (this also applies to the usual
 * evaluation with FEValues). One could compute the cell matrix with
 * $\mathcal O((p+1)^7)$ operations by utilizing the tensor product
 * structure of the shape functions, as is done by the matrix-free framework
 * in deal.II. We refer to step-37 and the documentation of the
 * tensor-product-aware evaluators FEEvaluation for details on how an even
 * more efficient cell matrix computation could be realized.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::assemble_system(const Mapping<dim> &mapping)
 *   {
 *     TimerOutput::Scope scope(timer, "Assemble linear system");
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 2);
 *     FEValues<dim>     fe_values(mapping,
 *                             fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 *     FullMatrix<double> partial_matrix(dofs_per_cell, dim * n_q_points);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_rhs = 0.;
 *         fe_values.reinit(cell);
 * 
 *         for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
 *           {
 *             const double current_coefficient =
 *               coefficient(fe_values.quadrature_point(q_index));
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               {
 *                 for (unsigned int d = 0; d < dim; ++d)
 *                   partial_matrix(i, q_index * dim + d) =
 *                     std::sqrt(fe_values.JxW(q_index) * current_coefficient) *
 *                     fe_values.shape_grad(i, q_index)[d];
 *                 cell_rhs(i) +=
 *                   (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                    (-dim) *                            // f(x_q)
 *                    fe_values.JxW(q_index));            // dx
 *               }
 *           }
 * 
 *         partial_matrix.mTmult(cell_matrix, partial_matrix);
 * 
 *         cell->get_dof_indices(local_dof_indices);
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
 * <a name="Solutionofthelinearsystem"></a> 
 * <h3>Solution of the linear system</h3>
 *   

 * 
 * For solving the linear system, we pick a simple Jacobi-preconditioned
 * conjugate gradient solver, similar to the settings in the early tutorials.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::solve()
 *   {
 *     TimerOutput::Scope scope(timer, "Solve linear system");
 * 
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix);
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 *     constraints.distribute(solution);
 * 
 *     std::cout << "   Number of solver iterations:  "
 *               << solver_control.last_step() << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Outputofthesolutionandcomputationoferrors"></a> 
 * <h3>Output of the solution and computation of errors</h3>
 *   

 * 
 * In the next function we do various post-processing steps with the
 * solution, all of which involve the mapping in one way or the other.
 *   

 * 
 * The first operation we do is to write the solution as well as the
 * material ids to a VTU file. This is similar to what was done in many
 * other tutorial programs. The new ingredient presented in this tutorial
 * program is that we want to ensure that the data written to the file
 * used for visualization is actually a faithful representation of what
 * is used internally by deal.II. That is because most of the visualization
 * data formats only represent cells by their vertex coordinates, but
 * have no way of representing the curved boundaries that are used
 * in deal.II when using higher order mappings -- in other words, what
 * you see in the visualization tool is not actually what you are computing
 * on. (The same, incidentally, is true when using higher order shape
 * functions: Most visualization tools only render bilinear/trilinear
 * representations. This is discussed in detail in DataOut::build_patches().)
 *   

 * 
 * So we need to ensure that a high-order representation is written
 * to the file. We need to consider two particular topics. Firstly, we tell
 * the DataOut object via the DataOutBase::VtkFlags that we intend to
 * interpret the subdivisions of the elements as a high-order Lagrange
 * polynomial rather than a collection of bilinear patches.
 * Recent visualization programs, like ParaView version 5.5
 * or newer, can then render a high-order solution (see a <a
 * href="https://github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
 * page</a> for more details).
 * Secondly, we need to make sure that the mapping is passed to the
 * DataOut::build_patches() method. Finally, the DataOut class only prints
 * curved faces for <i>boundary</i> cells by default, so we need to ensure
 * that also inner cells are printed in a curved representation via the
 * mapping.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::postprocess(const Mapping<dim> &mapping)
 *   {
 *     {
 *       TimerOutput::Scope scope(timer, "Write output");
 * 
 *       DataOut<dim> data_out;
 * 
 *       DataOutBase::VtkFlags flags;
 *       flags.write_higher_order_cells = true;
 *       data_out.set_flags(flags);
 * 
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "solution");
 * 
 *       Vector<double> material_ids(triangulation.n_active_cells());
 *       for (const auto &cell : triangulation.active_cell_iterators())
 *         material_ids[cell->active_cell_index()] = cell->material_id();
 *       data_out.add_data_vector(material_ids, "material_ids");
 * 
 *       data_out.build_patches(mapping,
 *                              fe.degree,
 *                              DataOut<dim>::curved_inner_cells);
 * 
 *       std::ofstream file(
 *         ("solution-" +
 *          std::to_string(triangulation.n_global_levels() - 10 + 2 * dim) +
 *          ".vtu")
 *           .c_str());
 * 
 *       data_out.write_vtu(file);
 *     }
 * 
 * @endcode
 * 
 * The next operation in the postprocessing function is to compute the $L_2$
 * and $H^1$ errors against the analytical solution. As the analytical
 * solution is a quadratic polynomial, we expect a very accurate result at
 * this point. If we were solving on a simple mesh with planar faces and a
 * coefficient whose jumps are aligned with the faces between cells, then
 * we would expect the numerical result to coincide with the
 * analytical solution up to roundoff accuracy. However, since we are using
 * deformed cells following a sphere, which are only tracked by
 * polynomials of degree 4 (one more than the degree for the finite
 * elements), we will see that there is an error around $10^{-7}$. We could
 * get more accuracy by increasing the polynomial degree or refining the
 * mesh.
 * 
 * @code
 *     {
 *       TimerOutput::Scope scope(timer, "Compute error norms");
 * 
 *       Vector<double> norm_per_cell_p(triangulation.n_active_cells());
 * 
 *       VectorTools::integrate_difference(mapping,
 *                                         dof_handler,
 *                                         solution,
 *                                         ExactSolution<dim>(),
 *                                         norm_per_cell_p,
 *                                         QGauss<dim>(fe.degree + 2),
 *                                         VectorTools::L2_norm);
 *       std::cout << "   L2 error vs exact solution:   "
 *                 << norm_per_cell_p.l2_norm() << std::endl;
 * 
 *       VectorTools::integrate_difference(mapping,
 *                                         dof_handler,
 *                                         solution,
 *                                         ExactSolution<dim>(),
 *                                         norm_per_cell_p,
 *                                         QGauss<dim>(fe.degree + 2),
 *                                         VectorTools::H1_norm);
 *       std::cout << "   H1 error vs exact solution:   "
 *                 << norm_per_cell_p.l2_norm() << std::endl;
 *     }
 * 
 * @endcode
 * 
 * The final post-processing operation we do here is to compute an error
 * estimate with the KellyErrorEstimator. We use the exact same settings
 * as in the step-6 tutorial program, except for the fact that we also
 * hand in the mapping to ensure that errors are evaluated along the
 * curved element, consistent with the remainder of the program. However,
 * we do not really use the result here to drive a mesh adaptation step
 * (that would refine the mesh around the material interface along the
 * sphere), as the focus here is on the cost of this operation.
 * 
 * @code
 *     {
 *       TimerOutput::Scope scope(timer, "Compute error estimator");
 * 
 *       Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 *       KellyErrorEstimator<dim>::estimate(
 *         mapping,
 *         dof_handler,
 *         QGauss<dim - 1>(fe.degree + 1),
 *         std::map<types::boundary_id, const Function<dim> *>(),
 *         solution,
 *         estimated_error_per_cell);
 *       std::cout << "   Max cell-wise error estimate: "
 *                 << estimated_error_per_cell.linfty_norm() << std::endl;
 *     }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePoissonProblemrunfunction"></a> 
 * <h3>The PoissonProblem::run() function</h3>
 *   

 * 
 * Finally, we define the `run()` function that controls how we want to
 * execute this program (which is called by the main() function in the usual
 * way). We start by calling the `create_grid()` function that sets up our
 * geometry with the appropriate manifolds. We then run two instances of a
 * solver chain, starting from the setup of the equations, the assembly of
 * the linear system, its solution with a simple iterative solver, and the
 * postprocessing discussed above. The two instances differ in the way they
 * use the mapping. The first uses a conventional MappingQGeneric mapping
 * object which we initialize to a degree one more than we use for the
 * finite element &ndash; after all, we expect the geometry representation
 * to be the bottleneck as the analytic solution is only a quadratic
 * polynomial. (In reality, things are interlinked to quite some extent
 * because the evaluation of the polynomials in real coordinates involves
 * the mapping of a higher-degree polynomials, which represent some smooth
 * rational functions. As a consequence, higher-degree polynomials still pay
 * off, so it does not make sense to increase the degree of the mapping
 * further.) Once the first pass is completed, we let the timer print a
 * summary of the compute times of the individual stages.
 * 
 * @code
 *   template <int dim>
 *   void PoissonProblem<dim>::run()
 *   {
 *     create_grid();
 * 
 *     {
 *       std::cout << std::endl
 *                 << "====== Running with the basic MappingQGeneric class ====== "
 *                 << std::endl
 *                 << std::endl;
 * 
 *       MappingQGeneric<dim> mapping(fe.degree + 1);
 *       setup_system(mapping);
 *       assemble_system(mapping);
 *       solve();
 *       postprocess(mapping);
 * 
 *       timer.print_summary();
 *       timer.reset();
 *     }
 * 
 * @endcode
 * 
 * For the second instance, we instead set up the MappingQCache class. Its
 * use is very simple: After constructing it (with the degree, given that
 * we want it to show the correct degree functionality in other contexts),
 * we fill the cache via the MappingQCache::initialize() function. At this
 * stage, we specify which mapping we want to use (obviously, the same
 * MappingQGeneric as previously in order to repeat the same computations)
 * for the cache, and then run through the same functions again, now
 * handing in the modified mapping. In the end, we again print the
 * accumulated wall times since the reset to see how the times compare to
 * the original setting.
 * 
 * @code
 *     {
 *       std::cout
 *         << "====== Running with the optimized MappingQCache class ====== "
 *         << std::endl
 *         << std::endl;
 * 
 *       MappingQCache<dim> mapping(fe.degree + 1);
 *       {
 *         TimerOutput::Scope scope(timer, "Initialize mapping cache");
 *         mapping.initialize(MappingQGeneric<dim>(fe.degree + 1), triangulation);
 *       }
 *       std::cout << "   Memory consumption cache:     "
 *                 << 1e-6 * mapping.memory_consumption() << " MB" << std::endl;
 * 
 *       setup_system(mapping);
 *       assemble_system(mapping);
 *       solve();
 *       postprocess(mapping);
 * 
 *       timer.print_summary();
 *     }
 *   }
 * } // namespace Step65
 * 
 * 
 * 
 * int main()
 * {
 *   Step65::PoissonProblem<3> test_program;
 *   test_program.run();
 *   return 0;
 * }
 * @endcode
examples/step-65/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


如果我们用三度的多项式运行这个程序的三维版本，我们会得到以下程序输出。

@code
> make run
Scanning dependencies of target \step-65
[ 33%] Building CXX object CMakeFiles/\step-65.dir/\step-65.cc.o
[ 66%] Linking CXX executable \step-65
[ 66%] Built target \step-65
[100%] Run \step-65 with Release configuration


====== Running with the basic MappingQGeneric class ======


   Number of active cells:       6656
   Number of degrees of freedom: 181609
   Number of solver iterations:  285
   L2 error vs exact solution:   8.99339e-08
   H1 error vs exact solution:   6.45341e-06
   Max cell-wise error estimate: 0.00743406



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      49.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble linear system          |         1 |       5.8s |        12% |
| Compute constraints             |         1 |     0.109s |      0.22% |
| Compute error estimator         |         1 |      16.5s |        33% |
| Compute error norms             |         1 |      9.11s |        18% |
| Solve linear system             |         1 |      9.92s |        20% |
| Write output                    |         1 |      4.85s |       9.8% |
+---------------------------------+-----------+------------+------------+


====== Running with the optimized MappingQCache class ======


   Memory consumption cache:     22.9981 MB
   Number of active cells:       6656
   Number of degrees of freedom: 181609
   Number of solver iterations:  285
   L2 error vs exact solution:   8.99339e-08
   H1 error vs exact solution:   6.45341e-06
   Max cell-wise error estimate: 0.00743406



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      18.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble linear system          |         1 |      1.44s |       7.8% |
| Compute constraints             |         1 |   0.00336s |         0% |
| Compute error estimator         |         1 |     0.476s |       2.6% |
| Compute error norms             |         1 |     0.505s |       2.7% |
| Initialize mapping cache        |         1 |      4.96s |        27% |
| Solve linear system             |         1 |      9.95s |        54% |
| Write output                    |         1 |     0.875s |       4.8% |
+---------------------------------+-----------+------------+------------+


[100%] Built target run
@endcode



在讨论时序之前，我们看一下MappingQCache对象的内存消耗。我们的程序打印出它占用了23MB的内存。如果我们把这个数字与单个（解或右侧）向量的内存消耗联系起来，即1.5MB（即181,609个元素乘以每条8字节的双精度），或者与系统矩阵和稀疏模式所消耗的内存联系起来（274MB），我们就会发现，考虑到它的好处，这并不是一个过于沉重的数据结构。

在计时器方面，我们看到程序的整体运行时间明显改善了2.7倍。如果我们不考虑迭代求解器，它在两种情况下都是一样的（考虑到我们使用的简单预处理程序，以及稀疏矩阵-向量乘积对立方多项式的操作浪费的事实，它并不是最佳的），优势是几乎5倍。这对于一个线性静止问题来说是相当令人印象深刻的，而对于时间依赖性和非线性问题来说，成本的节省确实要突出得多，因为在这些问题上，会多次调用装配。如果我们研究一下各个组件，我们就会更清楚地了解发生了什么，以及为什么缓存会如此高效。在MappingQGeneric的情况下，基本上每一个涉及到映射的操作都需要至少5秒的时间来运行。规范的计算运行了两个 VectorTools::integrate_difference() 函数，每个都需要5秒。(约束条件的计算比较便宜，因为它只在边界的单元格中评估映射，用于边界条件的插值)。如果我们将这5秒与填充MappingQCache的时间相比较，即5.2秒（对于所有单元，而不仅仅是活动单元），很明显，在MappingQGeneric的情况下，映射支持点的计算比其他一切都重要。也许最引人注目的结果是误差估计器的时间，标记为 "计算误差估计器"，其中MappingQGeneric的实现需要17.3秒，而MappingQCache的变体不到0.5秒。前者如此昂贵的原因（比如说，比装配贵三倍）是误差估计涉及到对面的量的评估，其中网格中的每个面都要求映射的额外点，而这些点又要经过非常昂贵的TransfiniteInterpolationManifold类。由于每个单元有六个面，这种情况比在装配中发生得更频繁。同样，MappingQCache很好地消除了重复评估，将所有涉及流形的昂贵步骤汇总到一个被重复使用的初始化调用中。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-65.cc"
*/
