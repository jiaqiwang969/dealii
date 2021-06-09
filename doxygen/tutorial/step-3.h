/**
@page step_3 The step-3 tutorial program
This tutorial depends on step-2.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thebasicsetupoffiniteelementmethods">The basic set up of finite element methods</a>
        <li><a href="#Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"> Should we multiply by a test function from the left or from the right? </a>
        <li><a href="#Computingthematrixandrighthandsidevector"> Computing the matrix and right hand side vector </a>
        <li><a href="#Abouttheimplementation">About the implementation</a>
        <li><a href="#Anoteontypes"> A note on types </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Manynewincludefiles">Many new include files</a>
        <li><a href="#ThecodeStep3codeclass">The <code>Step3</code> class</a>
      <ul>
        <li><a href="#Step3Step3">Step3::Step3</a>
        <li><a href="#Step3make_grid">Step3::make_grid</a>
        <li><a href="#Step3setup_system">Step3::setup_system</a>
        <li><a href="#Step3assemble_system">Step3::assemble_system</a>
        <li><a href="#Step3solve">Step3::solve</a>
        <li><a href="#Step3output_results">Step3::output_results</a>
        <li><a href="#Step3run">Step3::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
        <li><a href="#UsingHDF5tooutputthesolutionandadditionaldata">Using HDF5 to output the solution and additional data</a>
      <ul>
        <li><a href="#Changingtheoutputtoh5"> Changing the output to .h5</a>
        <li><a href="#Addingthepointvalueandthemeanseeextensionaboveintotheh5file"> Adding the point value and the mean (see extension above) into the .h5 file</a>
      </ul>
        <li><a href="#UsingRandggplot2togenerateplots"> Using R and ggplot2 to generate plots</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-3/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{10} 

<a name="Thebasicsetupoffiniteelementmethods"></a><h3>The basic set up of finite element methods</h3>


这是第一个我们实际使用有限元来计算的例子。我们将解决一个简单的泊松方程，其边界值为零，但右手边非零。

@f{align*}


  -\Delta u &= f \qquad\qquad & \text{in}\ \Omega,
  \\
  u &= 0 \qquad\qquad & \text{on}\ \partial\Omega.


@f}

我们将在正方形 $\Omega=[-1,1]^2$ 上求解这个方程，你已经在步骤1和步骤2中学习了如何生成网格。在这个程序中，我们也将只考虑 $f(\mathbf x)=1$ 这个特殊情况，并在下一个教程程序中再来讨论如何实现更一般的情况，即步骤4。

如果你学过有限元方法的基本知识，你会记得我们需要采取的步骤，用有限维度的近似方法来近似解 $u$ 。具体来说，我们首先需要推导出上述方程的弱形式，通过将方程乘以测试函数 $\varphi$ <i>from the left</i>（我们将在下面回到从左而非从右相乘的原因）并在域 $\Omega$ 上积分得到。

@f{align*}


  -\int_\Omega \varphi \Delta u = \int_\Omega \varphi f.


@f}

这可以通过部件进行整合。

@f{align*}
  \int_\Omega \nabla\varphi \cdot \nabla u


  -
  \int_{\partial\Omega} \varphi \mathbf{n}\cdot \nabla u
   = \int_\Omega \varphi f.


@f}

测试函数 $\varphi$ 必须满足同样的边界条件（用数学术语来说：它需要来自我们寻求解决方案的集合的切线空间），因此在边界上 $\varphi=0$ ，因此我们正在寻找的弱形式为

@f{align*}
  (\nabla\varphi, \nabla u)
   = (\varphi, f),


@f}

其中我们使用了常用的符号  $(a,b)=\int_\Omega a\; b$  。然后，问题要求从适当的空间（这里是空间 $H^1$ ）中找出一个函数 $u$ ，对于该函数，这一声明对于所有测试函数 $\varphi$ 都是真的。

当然，在一般情况下，我们无法在计算机上找到这样的函数，而是寻求一个近似值 $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf
x)$  ，其中 $U_j$ 是我们需要确定的未知膨胀系数（这个问题的 "自由度"）， $\varphi_i(\mathbf x)$ 是我们将使用的有限元形状函数。为了定义这些形状函数，我们需要以下内容。

- 一个用来定义形状函数的网格。你已经看到如何在步骤1和步骤2中生成和操作描述网格的对象。

- 一个描述我们想在参考单元上使用的形状函数的有限元（在deal.II中总是单位间隔 $[0,1]$ 、单位正方形 $[0,1]^2$ 或单位立方体 $[0,1]^3$ ，取决于你在哪个空间维度工作）。在步骤2中，我们已经使用了FE_Q<2>类型的对象，它表示通常的拉格朗日元素，通过对支持点的插值来定义形状函数。最简单的是FE_Q<2>(1)，它使用1度的多项式。在2D中，这些通常被称为<i>bilinear</i>，因为它们在参考单元的两个坐标中都是线性的。(在1d中，它们是<i>linear</i>，在3d中是<i>tri-linear</i>；然而，在deal.II文档中，我们经常不做这种区分，而总是简单地称这些函数为 "线性"。)

- 一个DoFHandler对象，以有限元对象提供的参考单元描述为基础，枚举网格上的所有自由度。你也已经在步骤2中看到了如何做到这一点。

- 一个映射，告诉你如何从参考单元上的有限元类定义的形状函数中获得实数单元上的形状函数。默认情况下，除非你明确说明，否则deal.II将使用（双，三）线性映射，所以在大多数情况下，你不必担心这个步骤。

通过这些步骤，我们现在有一组函数 $\varphi_i$ ，我们可以定义离散问题的弱形式：找到一个函数 $u_h$ ，即找到上面提到的扩展系数 $U_j$ ，以便

@f{align*}
  (\nabla\varphi_i, \nabla u_h)
   = (\varphi_i, f),
   \qquad\qquad
   i=0\ldots N-1.


@f}

请注意，我们在此遵循惯例，即一切从零开始计算，这在C和C++中很常见。如果你插入表示法 $u_h(\mathbf x)=\sum_j U_j
\varphi_j(\mathbf x)$ ，这个方程可以重写为一个线性系统，然后观察到

@f{align*}{
  (\nabla\varphi_i, \nabla u_h)
  &= \left(\nabla\varphi_i, \nabla \Bigl[\sum_j U_j \varphi_j\Bigr]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \left[U_j \varphi_j\right]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \varphi_j \right) U_j.


@f}

有了这个，问题就成了。找到一个向量 $U$ ，以便

@f{align*}{
  A U = F,


@f}

其中矩阵 $A$ 和右手边 $F$ 定义为

@f{align*}
  A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j),
  \\
  F_i &= (\varphi_i, f).


@f}






<a name="Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"></a><h3> Should we multiply by a test function from the left or from the right? </h3>


在我们继续描述如何计算这些数量之前，请注意，如果我们从<i>right</i>乘以一个测试函数而不是从左边乘以原方程，那么我们将得到一个形式为的线性系统

@f{align*}
  U^T A = F^T


@f}

有一个行向量  $F^T$  。通过转置这个系统，这当然等同于解决了

@f{align*}
  A^T U = F


@f}

这里与上面的 $A=A^T$ 相同。但一般来说不是，为了避免任何形式的混淆，经验表明，只要养成从左边而不是从右边乘方程的习惯（正如数学文献中经常做的那样），就可以避免一类常见的错误，因为在比较理论和实现时，矩阵会自动正确，不需要转置。本教程的第一个例子见第9步，我们有一个非对称的双线性方程，对于这个方程，我们从右面还是从左面相乘是有区别的。




<a name="Computingthematrixandrighthandsidevector"></a><h3> Computing the matrix and right hand side vector </h3>


现在我们知道我们需要什么（即：持有矩阵和向量的对象，以及计算 $A_{ij},F_i$ 的方法），我们可以看看需要什么来实现这一点。

-  $A$ 的对象是SparseMatrix类型，而 $U$ 和 $F$ 的对象则是Vector类型。我们将在下面的程序中看到哪些类是用来解决线性系统的。

- 我们需要一种方法来形成积分。在有限元方法中，最常见的是使用正交法，也就是说，积分被每个单元上的一组点的加权和所取代。也就是说，我们首先将 $\Omega$ 的积分分成所有单元的积分，@f{align*}
    A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j)
    = \sum_{K \in {\mathbb T}} \int_K \nabla\varphi_i \cdot \nabla \varphi_j,
    \\
    F_i &= (\varphi_i, f)
    = \sum_{K \in {\mathbb T}} \int_K \varphi_i f,
  @f}

  然后用正交法对每个单元的贡献进行近似。   @f{align*}
    A^K_{ij} &=
    \int_K \nabla\varphi_i \cdot \nabla \varphi_j
    \approx
    \sum_q \nabla\varphi_i(\mathbf x^K_q) \cdot \nabla
    \varphi_j(\mathbf x^K_q) w_q^K,
    \\
    F^K_i &=
    \int_K \varphi_i f
    \approx
    \sum_q \varphi_i(\mathbf x^K_q) f(\mathbf x^K_q) w^K_q,
  @f}

  其中 $\mathbf x^K_q$ 是 $q$ 单元上的第三个正交点 $K$ ， $w^K_q$ 是 $q$ 的正交权。这样做需要有不同的部分，接下来我们将依次讨论它们。

- 首先，我们需要一种方法来描述正交点的位置  $\mathbf x_q^K$  和它们的权重  $w^K_q$  。它们通常以与形状函数相同的方式从参考单元映射出来，即隐含地使用MappingQ1类，或者，如果你明确地说，通过从Mapping派生的其他类之一。参考单元上的位置和权重由派生自正交基类的对象来描述。通常，人们选择一个正交公式（即一组点和权重），使正交正好等于矩阵中的积分；这可以实现，因为积分中的所有因子都是多项式，由高斯正交公式完成，在QGauss类中实现。

- 然后我们需要一些东西来帮助我们在 $K$ 单元上评估 $\varphi_i(\mathbf x^K_q)$ 。这就是FEValues类的作用：它需要一个有限元对象来描述参考单元上的 $\varphi$ ，一个正交对象来描述正交点和权重，以及一个映射对象（或隐含地采用MappingQ1类），并在位于 $K$ 的正交点上提供形状函数的值和导数，以及积分所需的各种其他信息。

FEValues确实是装配过程中的核心类。你可以这样看待它。FiniteElement和派生类描述了形状<i>functions</i>，即无限维度的对象：函数在每一点都有值。由于理论上的原因，我们需要这样做，因为我们想用函数的积分来进行分析。然而，对于计算机来说，这是一个非常困难的概念，因为它们一般只能处理有限的信息量，所以我们用正交点上的和来代替积分，我们通过使用定义在参考单元（正交对象）上的点映射（映射对象）到真实单元上的点来获得。实质上，我们将问题简化为我们只需要有限的信息，即形状函数值和导数、正交权重、法向量等，只需要在有限的点集合上。FEValues类就是将这三个部分结合在一起，并在一个特定的单元上提供这个有限的信息集  $K$  。当我们组装下面的线性系统时，你会看到它的作用。

值得注意的是，如果你只是在应用程序中自己创建这三个对象，并自己处理这些信息，那么所有这些也都可以实现。然而，这样做既不简单（FEValues类提供的正是你实际需要的信息），也不快：FEValues类经过高度优化，只在每个单元中计算你需要的特定信息；如果有任何东西可以从上一个单元中重复使用，那么它就会这样做，而且该类中有很多代码可以确保在有利的地方进行缓存。

这个介绍的最后一块是要提到，在得到一个线性系统后，要用迭代求解器进行求解，然后进行后处理：我们用DataOut类创建一个输出文件，然后可以用一个常见的可视化程序进行可视化。

 @note  前面对任何有限元实现的所有重要步骤的概述，在deal.II中也有对应的内容：该库可以自然地归纳为若干 "模块"，涵盖刚才概述的基本概念。你可以通过本页面顶部的标签访问这些模块。在<a href="index.html">front page of the deal.II manual</a>上也有对最基本的概念组的概述。




<a name="Abouttheimplementation"></a><h3>About the implementation</h3>


虽然这是你能用有限元方法解决的最简单的方程，但这个程序显示了大多数有限元程序的基本结构，也是几乎所有下面的程序基本上都会遵循的模板。具体来说，这个程序的主类看起来像这样。

@code
class Step3
{
  public:
    Step3 ();
    void run ();


  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;


    Triangulation<2>     triangulation;
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;


    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};
@endcode



这遵循了<a
href="http://en.wikipedia.org/wiki/Encapsulation_(object-oriented_programming)">data
encapsulation</a>的面向对象编程口号，也就是说，我们尽力将这个类的几乎所有内部细节隐藏在外部无法访问的私有成员中。

让我们从成员变量开始。这些遵循我们在上面的要点中所概述的构建模块，即我们需要一个三角形和一个DoFHandler对象，以及一个描述我们想要使用的各种形状函数的有限元对象。第二组对象与线性代数有关：系统矩阵和右手边以及解向量，还有一个描述矩阵稀疏模式的对象。这就是这个类所需要的全部内容（也是任何静止PDE的求解器所需要的基本内容），并且需要在整个程序中存活。与此相反，我们在装配时需要的FEValues对象只在整个装配过程中需要，因此我们在进行装配的函数中把它作为一个局部对象来创建，并在结束时再次销毁它。

其次，让我们来看看成员函数。这些，也已经构成了几乎所有下面的教程程序都会使用的共同结构。   <ul>   <li>   <code>make_grid()</code>  : 这就是人们所说的<i>preprocessing function</i>。顾名思义，它设置了存储三角图的对象。在以后的例子中，它还可以处理边界条件、几何形状等。     <li>   <code>setup_system()</code>  : 这是一个函数，其中设置了解决问题所需的所有其他数据结构。特别是，它将初始化DoFHandler对象并正确确定与线性代数有关的各种对象的大小。这个函数通常与上面的预处理函数分开，因为在一个与时间相关的程序中，每当网格被自适应细化时（我们将在步骤6中看到如何做），它可能至少每隔几个时间步就会被调用。另一方面，在上面的预处理函数中，设置网格本身只在程序开始时进行一次，因此，它被分离成自己的函数。     <li>   <code>assemble_system()</code>  : 这就是计算矩阵和右手边的内容的地方，在上面的介绍中已经详细讨论过。由于对这个线性系统进行处理在概念上与计算其条目有很大不同，我们将其与以下函数分开。     <li>   <code>solve()</code>  : 这就是我们计算线性系统 $U$ 的解的函数。在当前的程序中，这是一个简单的任务，因为矩阵是如此简单，但只要问题不再那么微不足道，它就会成为程序规模的重要部分（例如，一旦你对库有了更多的了解，请参阅步骤20，步骤22，或步骤31）。     <li>   <code>output_results()</code>  : 最后，当你计算出一个解决方案后，你可能想用它做一些事情。例如，你可能想以可视化的格式输出它，或者你可能想计算你感兴趣的量：例如，热交换器中的热通量、机翼的空气摩擦系数、最大桥梁载荷，或者仅仅是某一点上的数值解的值。因此，这个函数是对你的解进行后处理的地方。   </ul> 所有这些都是由单一的公共函数（除构造函数外），即 <code>run()</code> 函数来支撑的。它是在创建这种类型的对象的地方被调用的，它是按正确顺序调用所有其他函数的函数。把这个操作封装到 <code>run()</code> 函数中，而不是从 <code>main()</code> 中调用所有其他函数，确保你可以改变这个类中的关注点分离的实现方式。例如，如果其中一个函数变得太大了，你可以把它分成两个，而你唯一需要关注的地方就是这个类中的变化，而不是其他地方。

如上所述，你会看到这种一般的结构&mdash；有时在函数名称的拼写上会有一些变化，但基本上是按照这种功能分离的顺序&mdash；在下面的许多教程程序中也是如此。




<a name="Anoteontypes"></a><h3> A note on types </h3>


deal.II通过命名空间 dealii::types. 中的别名定义了一些积分%类型（在前一句中，"积分 "一词被用作与名词 "整数 "相对应的<i>adjective</i>。它不应该与表示曲线或曲面下的面积或体积的<i>noun</i>"积分 "混淆起来。形容词 "积分 "在C++世界中被广泛使用，如 "积分类型"、"积分常数 "等。）特别是，在这个程序中，你会在几个地方看到 types::global_dof_index ：一个整数类型，用来表示自由度的<i>global</i>索引，即在定义在三角形之上的DoFHandler对象中特定自由度的索引（而不是特定单元中的特定自由度的索引）。对于当前的程序（以及几乎所有的教程程序），你将有几千个到几百万个全局未知数（而且，对于 $Q_1$ 元素，你将有4个<i>locally on each cell</i>的2D和8个3D）。因此，允许为全局DoF指数存储足够大的数字的数据类型是 <code>unsigned int</code> ，因为它允许存储0到略高于40亿的数字（在大多数系统中，整数是32位的）。事实上，这就是 types::global_dof_index 的作用。

那么，为什么不马上使用 <code>unsigned int</code> 呢？deal.II在7.3版本之前一直是这样做的。然而，deal.II支持非常大的计算（通过步骤40中讨论的框架），当分布在几千个处理器上时，可能有超过40亿个未知数。因此，有些情况下 <code>unsigned int</code> 不够大，我们需要一个64位的无符号积分类型。为了实现这一点，我们引入了 types::global_dof_index ，它默认被定义为<code>unsigned int</code>，而如果有必要，可以通过在配置过程中传递一个特定的标志，将其定义为<code>unsigned long long int</code>（见ReadMe文件）。

这涵盖了技术方面。但是还有一个文档的目的：在图书馆和建立在它之上的代码中，如果你看到一个地方使用数据类型 types::global_dof_index, ，你就会立即知道被引用的数量实际上是一个全局dof指数。如果我们只是使用 <code>unsigned int</code> （它也可能是一个局部索引，一个边界指示器，一个材料ID，等等），就不会有这样的意义了。立即知道一个变量指的是什么也有助于避免错误：如果你看到一个 types::global_dof_index 类型的对象被分配给 types::subdomain_id, 类型的变量，这很明显一定有一个错误，尽管它们都是用无符号整数表示，因此编译器不会抱怨。

在更实际的情况下，这种类型的存在意味着在装配过程中，我们创建一个 $4\times 4$ 矩阵（在2d中，使用 $Q_1$ 元素）来表示我们当前所在单元的贡献，然后我们需要将这个矩阵的元素添加到全局（系统）矩阵的相应元素中。为此，我们需要获得当前单元的局部自由度的全局指数，为此我们将始终使用下面这段代码。

@code
  cell->get_dof_indices (local_dof_indices);
@endcode

其中 <code>local_dof_indices</code> 被声明为

@code
  std::vector<types::global_dof_index> local_dof_indices (fe.n_dofs_per_cell());
@endcode

这个变量的名字可能有点名不副实--它代表 "在当前单元上局部定义的那些自由度的全局指数"--但持有这种信息的变量在整个库中普遍是这样命名的。

 @note   types::global_dof_index  并不是这个命名空间中定义的唯一类型。相反，有一整个系列，包括 types::subdomain_id,  types::boundary_id, 和 types::material_id. 所有这些都是整数数据类型的别名，但正如上面所解释的，它们被用于整个库，以便（i）变量的意图变得更容易辨别，以及（ii）如果有必要，可以将实际类型改为一个更大的类型，而不必翻阅整个库，找出 <code>unsigned int</code> 的特定使用是否对应于，例如，一个材料指标。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Manynewincludefiles"></a> 
 * <h3>Many new include files</h3>
 * 

 * 
 * These include files are already known to you. They declare the classes
 * which handle triangulations and enumeration of degrees of freedom:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/dofs/dof_handler.h>
 * @endcode
 * 
 * And this is the file in which the functions are declared that create grids:
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h>
 * 
 * @endcode
 * 
 * This file contains the description of the Lagrange interpolation finite
 * element:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * 
 * @endcode
 * 
 * And this file is needed for the creation of sparsity patterns of sparse
 * matrices, as shown in previous examples:
 * 
 * @code
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * The next two files are needed for assembling the matrix using quadrature on
 * each cell. The classes declared in them will be explained below:
 * 
 * @code
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * @endcode
 * 
 * The following three include files we need for the treatment of boundary
 * values:
 * 
 * @code
 * #include <deal.II/base/function.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * @endcode
 * 
 * We're now almost to the end. The second to last group of include files is
 * for the linear algebra which we employ to solve the system of equations
 * arising from the finite element discretization of the Laplace equation. We
 * will use vectors and full matrices for assembling the system of equations
 * locally on each cell, and transfer the results into a sparse matrix. We
 * will then use a Conjugate Gradient solver to solve the problem, for which
 * we need a preconditioner (in this program, we use the identity
 * preconditioner which does nothing, but we need to include the file anyway):
 * 
 * @code
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * @endcode
 * 
 * Finally, this is for output to a file and to the console:
 * 
 * @code
 * #include <deal.II/numerics/data_out.h>
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * ...and this is to import the deal.II namespace into the global scope:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep3codeclass"></a> 
 * <h3>The <code>Step3</code> class</h3>
 * 

 * 
 * Instead of the procedural programming of previous examples, we encapsulate
 * everything into a class for this program. The class consists of functions
 * which each perform certain aspects of a finite element program, a `main`
 * function which controls what is done first and what is done next, and a
 * list of member variables.
 * 

 * 
 * The public part of the class is rather short: it has a constructor and a
 * function `run` that is called from the outside and acts as something like
 * the `main` function: it coordinates which operations of this class shall be
 * run in which order. Everything else in the class, i.e. all the functions
 * that actually do anything, are in the private section of the class:
 * 
 * @code
 * class Step3
 * {
 * public:
 *   Step3();
 * 
 *   void run();
 * 
 * @endcode
 * 
 * Then there are the member functions that mostly do what their names
 * suggest and whose have been discussed in the introduction already. Since
 * they do not need to be called from outside, they are made private to this
 * class.
 * 

 * 
 * 
 * @code
 * private:
 *   void make_grid();
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void output_results() const;
 * 
 * @endcode
 * 
 * And finally we have some member variables. There are variables describing
 * the triangulation and the global numbering of the degrees of freedom (we
 * will specify the exact polynomial degree of the finite element in the
 * constructor of this class)...
 * 
 * @code
 *   Triangulation<2> triangulation;
 *   FE_Q<2>          fe;
 *   DoFHandler<2>    dof_handler;
 * 
 * @endcode
 * 
 * ...variables for the sparsity pattern and values of the system matrix
 * resulting from the discretization of the Laplace equation...
 * 
 * @code
 *   SparsityPattern      sparsity_pattern;
 *   SparseMatrix<double> system_matrix;
 * 
 * @endcode
 * 
 * ...and variables which will hold the right hand side and solution
 * vectors.
 * 
 * @code
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * @endcode
 * 
 * 
 * <a name="Step3Step3"></a> 
 * <h4>Step3::Step3</h4>
 * 

 * 
 * Here comes the constructor. It does not much more than first to specify
 * that we want bi-linear elements (denoted by the parameter to the finite
 * element object, which indicates the polynomial degree), and to associate
 * the dof_handler variable to the triangulation we use. (Note that the
 * triangulation isn't set up with a mesh at all at the present time, but the
 * DoFHandler doesn't care: it only wants to know which triangulation it will
 * be associated with, and it only starts to care about an actual mesh once
 * you try to distribute degree of freedom on the mesh using the
 * distribute_dofs() function.) All the other member variables of the Step3
 * class have a default constructor which does all we want.
 * 
 * @code
 * Step3::Step3()
 *   : fe(1)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3make_grid"></a> 
 * <h4>Step3::make_grid</h4>
 * 

 * 
 * Now, the first thing we've got to do is to generate the triangulation on
 * which we would like to do our computation and number each vertex with a
 * degree of freedom. We have seen these two steps in step-1 and step-2
 * before, respectively.
 * 

 * 
 * This function does the first part, creating the mesh.  We create the grid
 * and refine all cells five times. Since the initial grid (which is the
 * square $[-1,1] \times [-1,1]$) consists of only one cell, the final grid
 * has 32 times 32 cells, for a total of 1024.
 * 

 * 
 * Unsure that 1024 is the correct number? We can check that by outputting the
 * number of cells using the <code>n_active_cells()</code> function on the
 * triangulation.
 * 
 * @code
 * void Step3::make_grid()
 * {
 *   GridGenerator::hyper_cube(triangulation, -1, 1);
 *   triangulation.refine_global(5);
 * 
 *   std::cout << "Number of active cells: " << triangulation.n_active_cells()
 *             << std::endl;
 * }
 * 
 * @endcode
 * 
 * @note We call the Triangulation::n_active_cells() function, rather than
 * Triangulation::n_cells(). Here, <i>active</i> means the cells that aren't
 * refined any further. We stress the adjective "active" since there are more
 * cells, namely the parent cells of the finest cells, their parents, etc, up
 * to the one cell which made up the initial grid. Of course, on the next
 * coarser level, the number of cells is one quarter that of the cells on the
 * finest level, i.e. 256, then 64, 16, 4, and 1. If you called
 * <code>triangulation.n_cells()</code> instead in the code above, you would
 * consequently get a value of 1365 instead. On the other hand, the number of
 * cells (as opposed to the number of active cells) is not typically of much
 * interest, so there is no good reason to print it.
 * 

 * 
 * 

 * 
 * 
 * <a name="Step3setup_system"></a> 
 * <h4>Step3::setup_system</h4>
 * 

 * 
 * Next we enumerate all the degrees of freedom and set up matrix and vector
 * objects to hold the system data. Enumerating is done by using
 * DoFHandler::distribute_dofs(), as we have seen in the step-2 example. Since
 * we use the FE_Q class and have set the polynomial degree to 1 in the
 * constructor, i.e. bilinear elements, this associates one degree of freedom
 * with each vertex. While we're at generating output, let us also take a look
 * at how many degrees of freedom are generated:
 * 
 * @code
 * void Step3::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 *   std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 * @endcode
 * 
 * There should be one DoF for each vertex. Since we have a 32 times 32
 * grid, the number of DoFs should be 33 times 33, or 1089.
 * 

 * 
 * As we have seen in the previous example, we set up a sparsity pattern by
 * first creating a temporary structure, tagging those entries that might be
 * nonzero, and then copying the data over to the SparsityPattern object
 * that can then be used by the system matrix.
 * 
 * @code
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *   sparsity_pattern.copy_from(dsp);
 * 
 * @endcode
 * 
 * Note that the SparsityPattern object does not hold the values of the
 * matrix, it only stores the places where entries are. The entries
 * themselves are stored in objects of type SparseMatrix, of which our
 * variable system_matrix is one.
 *   

 * 
 * The distinction between sparsity pattern and matrix was made to allow
 * several matrices to use the same sparsity pattern. This may not seem
 * relevant here, but when you consider the size which matrices can have,
 * and that it may take some time to build the sparsity pattern, this
 * becomes important in large-scale problems if you have to store several
 * matrices in your program.
 * 
 * @code
 *   system_matrix.reinit(sparsity_pattern);
 * 
 * @endcode
 * 
 * The last thing to do in this function is to set the sizes of the right
 * hand side vector and the solution vector to the right values:
 * 
 * @code
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Step3assemble_system"></a> 
 * <h4>Step3::assemble_system</h4>
 * 

 * 
 * 

 * 
 * The next step is to compute the entries of the matrix and right hand side
 * that form the linear system from which we compute the solution. This is the
 * central function of each finite element program and we have discussed the
 * primary steps in the introduction already.
 * 

 * 
 * The general approach to assemble matrices and vectors is to loop over all
 * cells, and on each cell compute the contribution of that cell to the global
 * matrix and right hand side by quadrature. The point to realize now is that
 * we need the values of the shape functions at the locations of quadrature
 * points on the real cell. However, both the finite element shape functions
 * as well as the quadrature points are only defined on the reference
 * cell. They are therefore of little help to us, and we will in fact hardly
 * ever query information about finite element shape functions or quadrature
 * points from these objects directly.
 * 

 * 
 * Rather, what is required is a way to map this data from the reference cell
 * to the real cell. Classes that can do that are derived from the Mapping
 * class, though one again often does not have to deal with them directly:
 * many functions in the library can take a mapping object as argument, but
 * when it is omitted they simply resort to the standard bilinear Q1
 * mapping. We will go this route, and not bother with it for the moment (we
 * come back to this in step-10, step-11, and step-12).
 * 

 * 
 * So what we now have is a collection of three classes to deal with: finite
 * element, quadrature, and mapping objects. That's too much, so there is one
 * type of class that orchestrates information exchange between these three:
 * the FEValues class. If given one instance of each three of these objects
 * (or two, and an implicit linear mapping), it will be able to provide you
 * with information about values and gradients of shape functions at
 * quadrature points on a real cell.
 * 

 * 
 * Using all this, we will assemble the linear system for this problem in the
 * following function:
 * 
 * @code
 * void Step3::assemble_system()
 * {
 * @endcode
 * 
 * Ok, let's start: we need a quadrature formula for the evaluation of the
 * integrals on each cell. Let's take a Gauss formula with two quadrature
 * points in each direction, i.e. a total of four points since we are in
 * 2D. This quadrature formula integrates polynomials of degrees up to three
 * exactly (in 1D). It is easy to check that this is sufficient for the
 * present problem:
 * 
 * @code
 *   QGauss<2> quadrature_formula(fe.degree + 1);
 * @endcode
 * 
 * And we initialize the object which we have briefly talked about above. It
 * needs to be told which finite element we want to use, and the quadrature
 * points and their weights (jointly described by a Quadrature object). As
 * mentioned, we use the implied Q1 mapping, rather than specifying one
 * ourselves explicitly. Finally, we have to tell it what we want it to
 * compute on each cell: we need the values of the shape functions at the
 * quadrature points (for the right hand side $(\varphi_i,f)$), their
 * gradients (for the matrix entries $(\nabla \varphi_i, \nabla
 * \varphi_j)$), and also the weights of the quadrature points and the
 * determinants of the Jacobian transformations from the reference cell to
 * the real cells.
 *   

 * 
 * This list of what kind of information we actually need is given as a
 * collection of flags as the third argument to the constructor of
 * FEValues. Since these values have to be recomputed, or updated, every
 * time we go to a new cell, all of these flags start with the prefix
 * <code>update_</code> and then indicate what it actually is that we want
 * updated. The flag to give if we want the values of the shape functions
 * computed is #update_values; for the gradients it is
 * #update_gradients. The determinants of the Jacobians and the quadrature
 * weights are always used together, so only the products (Jacobians times
 * weights, or short <code>JxW</code>) are computed; since we need them, we
 * have to list #update_JxW_values as well:
 * 
 * @code
 *   FEValues<2> fe_values(fe,
 *                         quadrature_formula,
 *                         update_values | update_gradients | update_JxW_values);
 * @endcode
 * 
 * The advantage of this approach is that we can specify what kind of
 * information we actually need on each cell. It is easily understandable
 * that this approach can significantly speed up finite element computations,
 * compared to approaches where everything, including second derivatives,
 * normal vectors to cells, etc are computed on each cell, regardless of
 * whether they are needed or not.
 *   

 * 
 * @note The syntax <code>update_values | update_gradients |
 * update_JxW_values</code> is not immediately obvious to anyone not
 * used to programming bit operations in C for years already. First,
 * <code>operator|</code> is the <i>bitwise or operator</i>, i.e.,
 * it takes two integer arguments that are interpreted as bit
 * patterns and returns an integer in which every bit is set for
 * which the corresponding bit is set in at least one of the two
 * arguments. For example, consider the operation
 * <code>9|10</code>. In binary, <code>9=0b1001</code> (where the
 * prefix <code>0b</code> indicates that the number is to be
 * interpreted as a binary number) and <code>10=0b1010</code>. Going
 * through each bit and seeing whether it is set in one of the
 * argument, we arrive at <code>0b1001|0b1010=0b1011</code> or, in
 * decimal notation, <code>9|10=11</code>. The second piece of
 * information you need to know is that the various
 * <code>update_*</code> flags are all integers that have <i>exactly
 * one bit set</i>. For example, assume that
 * <code>update_values=0b00001=1</code>,
 * <code>update_gradients=0b00010=2</code>,
 * <code>update_JxW_values=0b10000=16</code>. Then
 * <code>update_values | update_gradients | update_JxW_values =
 * 0b10011 = 19</code>. In other words, we obtain a number that
 * <i>encodes a binary mask representing all of the operations you
 * want to happen</i>, where each operation corresponds to exactly
 * one bit in the integer that, if equal to one, means that a
 * particular piece should be updated on each cell and, if it is
 * zero, means that we need not compute it. In other words, even
 * though <code>operator|</code> is the <i>bitwise OR operation</i>,
 * what it really represents is <i>I want this AND that AND the
 * other</i>. Such binary masks are quite common in C programming,
 * but maybe not so in higher level languages like C++, but serve
 * the current purpose quite well.
 * 

 * 
 * For use further down below, we define a shortcut for a value that will
 * be used very frequently. Namely, an abbreviation for the number of degrees
 * of freedom on each cell (since we are in 2D and degrees of freedom are
 * associated with vertices only, this number is four, but we rather want to
 * write the definition of this variable in a way that does not preclude us
 * from later choosing a different finite element that has a different
 * number of degrees of freedom per cell, or work in a different space
 * dimension).
 *   

 * 
 * In general, it is a good idea to use a symbolic name instead of
 * hard-coding these numbers even if you know them, since for example,
 * you may want to change the finite element at some time. Changing the
 * element would have to be done in a different function and it is easy
 * to forget to make a corresponding change in another part of the program.
 * It is better to not rely on your own calculations, but instead ask
 * the right object for the information: Here, we ask the finite element
 * to tell us about the number of degrees of freedom per cell and we
 * will get the correct number regardless of the space dimension or
 * polynomial degree we may have chosen elsewhere in the program.
 *   

 * 
 * The shortcut here, defined primarily to discuss the basic concept
 * and not because it saves a lot of typing, will then make the following
 * loops a bit more readable. You will see such shortcuts in many places in
 * larger programs, and `dofs_per_cell` is one that is more or less the
 * conventional name for this kind of object.
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 * @endcode
 * 
 * Now, we said that we wanted to assemble the global matrix and vector
 * cell-by-cell. We could write the results directly into the global matrix,
 * but this is not very efficient since access to the elements of a sparse
 * matrix is slow. Rather, we first compute the contribution of each cell in
 * a small matrix with the degrees of freedom on the present cell, and only
 * transfer them to the global matrix when the computations are finished for
 * this cell. We do the same for the right hand side vector. So let's first
 * allocate these objects (these being local objects, all degrees of freedom
 * are coupling with all others, and we should use a full matrix object
 * rather than a sparse one for the local operations; everything will be
 * transferred to a global sparse matrix later on):
 * 
 * @code
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 * @endcode
 * 
 * When assembling the contributions of each cell, we do this with the local
 * numbering of the degrees of freedom (i.e. the number running from zero
 * through dofs_per_cell-1). However, when we transfer the result into the
 * global matrix, we have to know the global numbers of the degrees of
 * freedom. When we query them, we need a scratch (temporary) array for
 * these numbers (see the discussion at the end of the introduction for
 * the type, types::global_dof_index, used here):
 * 
 * @code
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Now for the loop over all cells. We have seen before how this works for a
 * triangulation. A DoFHandler has cell iterators that are exactly analogous
 * to those of a Triangulation, but with extra information about the degrees
 * of freedom for the finite element you're using. Looping over the active
 * cells of a degree-of-freedom handler works the same as for a triangulation.
 *   

 * 
 * Note that we declare the type of the cell as `const auto &` instead of
 * `auto` this time around. In step 1, we were modifying the cells of the
 * triangulation by flagging them with refinement indicators. Here we're only
 * examining the cells without modifying them, so it's good practice to
 * declare `cell` as `const` in order to enforce this invariant.
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 * @endcode
 * 
 * We are now sitting on one cell, and we would like the values and
 * gradients of the shape functions be computed, as well as the
 * determinants of the Jacobian matrices of the mapping between
 * reference cell and true cell, at the quadrature points. Since all
 * these values depend on the geometry of the cell, we have to have the
 * FEValues object re-compute them on each cell:
 * 
 * @code
 *       fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Next, reset the local cell's contributions to global matrix and
 * global right hand side to zero, before we fill them:
 * 
 * @code
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 * @endcode
 * 
 * Now it is time to start integration over the cell, which we
 * do by looping over all quadrature points, which we will
 * number by q_index.
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         {
 * @endcode
 * 
 * First assemble the matrix: For the Laplace problem, the
 * matrix on each cell is the integral over the gradients of
 * shape function i and j. Since we do not integrate, but
 * rather use quadrature, this is the sum over all
 * quadrature points of the integrands times the determinant
 * of the Jacobian matrix at the quadrature point times the
 * weight of this quadrature point. You can get the gradient
 * of shape function $i$ at quadrature point with number q_index by
 * using <code>fe_values.shape_grad(i,q_index)</code>; this
 * gradient is a 2-dimensional vector (in fact it is of type
 * Tensor@<1,dim@>, with here dim=2) and the product of two
 * such vectors is the scalar product, i.e. the product of
 * the two shape_grad function calls is the dot
 * product. This is in turn multiplied by the Jacobian
 * determinant and the quadrature point weight (that one
 * gets together by the call to FEValues::JxW() ). Finally,
 * this is repeated for all shape functions $i$ and $j$:
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices())
 *             for (const unsigned int j : fe_values.dof_indices())
 *               cell_matrix(i, j) +=
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                  fe_values.JxW(q_index));           // dx
 * 
 * @endcode
 * 
 * We then do the same thing for the right hand side. Here,
 * the integral is over the shape function i times the right
 * hand side function, which we choose to be the function
 * with constant value one (more interesting examples will
 * be considered in the following programs).
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices())
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                             1. *                                // f(x_q)
 *                             fe_values.JxW(q_index));            // dx
 *         }
 * @endcode
 * 
 * Now that we have the contribution of this cell, we have to transfer
 * it to the global matrix and right hand side. To this end, we first
 * have to find out which global numbers the degrees of freedom on this
 * cell have. Let's simply ask the cell for that information:
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices);
 * 
 * @endcode
 * 
 * Then again loop over all shape functions i and j and transfer the
 * local elements to the global matrix. The global numbers can be
 * obtained using local_dof_indices[i]:
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices())
 *         for (const unsigned int j : fe_values.dof_indices())
 *           system_matrix.add(local_dof_indices[i],
 *                             local_dof_indices[j],
 *                             cell_matrix(i, j));
 * 
 * @endcode
 * 
 * And again, we do the same thing for the right hand side vector.
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices())
 *         system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *     }
 * 
 * 
 * @endcode
 * 
 * Now almost everything is set up for the solution of the discrete
 * system. However, we have not yet taken care of boundary values (in fact,
 * Laplace's equation without Dirichlet boundary values is not even uniquely
 * solvable, since you can add an arbitrary constant to the discrete
 * solution). We therefore have to do something about the situation.
 *   

 * 
 * For this, we first obtain a list of the degrees of freedom on the
 * boundary and the value the shape function shall have there. For
 * simplicity, we only interpolate the boundary value function, rather than
 * projecting it onto the boundary. There is a function in the library which
 * does exactly this: VectorTools::interpolate_boundary_values(). Its
 * parameters are (omitting parameters for which default values exist and
 * that we don't care about): the DoFHandler object to get the global
 * numbers of the degrees of freedom on the boundary; the component of the
 * boundary where the boundary values shall be interpolated; the boundary
 * value function itself; and the output object.
 *   

 * 
 * The component of the boundary is meant as follows: in many cases, you may
 * want to impose certain boundary values only on parts of the boundary. For
 * example, you may have inflow and outflow boundaries in fluid dynamics, or
 * clamped and free parts of bodies in deformation computations of
 * bodies. Then you will want to denote these different parts of the
 * boundary by indicators, and tell the interpolate_boundary_values
 * function to only compute the boundary values on a certain part of the
 * boundary (e.g. the clamped part, or the inflow boundary). By default,
 * all boundaries have a 0 boundary indicator, unless otherwise specified. If
 * sections of the boundary have different boundary conditions, you have to
 * number those parts with different boundary indicators. The function call
 * below will then only determine boundary values for those parts of the
 * boundary for which the boundary indicator is in fact the zero specified as
 * the second argument.
 *   

 * 
 * The function describing the boundary values is an object of type Function
 * or of a derived class. One of the derived classes is
 * Functions::ZeroFunction, which describes (not unexpectedly) a function
 * which is zero everywhere. We create such an object in-place and pass it to
 * the VectorTools::interpolate_boundary_values() function.
 *   

 * 
 * Finally, the output object is a list of pairs of global degree of freedom
 * numbers (i.e. the number of the degrees of freedom on the boundary) and
 * their boundary values (which are zero here for all entries). This mapping
 * of DoF numbers to boundary values is done by the <code>std::map</code>
 * class.
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values;
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<2>(),
 *                                            boundary_values);
 * @endcode
 * 
 * Now that we got the list of boundary DoFs and their respective boundary
 * values, let's use them to modify the system of equations
 * accordingly. This is done by the following function call:
 * 
 * @code
 *   MatrixTools::apply_boundary_values(boundary_values,
 *                                      system_matrix,
 *                                      solution,
 *                                      system_rhs);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3solve"></a> 
 * <h4>Step3::solve</h4>
 * 

 * 
 * The following function simply solves the discretized equation. As the
 * system is quite a large one for direct solvers such as Gauss elimination or
 * LU decomposition, we use a Conjugate Gradient algorithm. You should
 * remember that the number of variables here (only 1089) is a very small
 * number for finite element computations, where 100.000 is a more usual
 * number.  For this number of variables, direct methods are no longer usable
 * and you are forced to use methods like CG.
 * 
 * @code
 * void Step3::solve()
 * {
 * @endcode
 * 
 * First, we need to have an object that knows how to tell the CG algorithm
 * when to stop. This is done by using a SolverControl object, and as
 * stopping criterion we say: stop after a maximum of 1000 iterations (which
 * is far more than is needed for 1089 variables; see the results section to
 * find out how many were really used), and stop if the norm of the residual
 * is below $10^{-12}$. In practice, the latter criterion will be the one
 * which stops the iteration:
 * 
 * @code
 *   SolverControl solver_control(1000, 1e-12);
 * @endcode
 * 
 * Then we need the solver itself. The template parameter to the SolverCG
 * class is the type of the vectors, and leaving the empty angle brackets
 * would indicate that we are taking the default argument (which is
 * <code>Vector@<double@></code>). However, we explicitly mention the template
 * argument:
 * 
 * @code
 *   SolverCG<Vector<double>> solver(solver_control);
 * 
 * @endcode
 * 
 * Now solve the system of equations. The CG solver takes a preconditioner
 * as its fourth argument. We don't feel ready to delve into this yet, so we
 * tell it to use the identity operation as preconditioner:
 * 
 * @code
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 * @endcode
 * 
 * Now that the solver has done its job, the solution variable contains the
 * nodal values of the solution function.
 * 
 * @code
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3output_results"></a> 
 * <h4>Step3::output_results</h4>
 * 

 * 
 * The last part of a typical finite element program is to output the results
 * and maybe do some postprocessing (for example compute the maximal stress
 * values at the boundary, or the average flux across the outflow, etc). We
 * have no such postprocessing here, but we would like to write the solution
 * to a file.
 * 
 * @code
 * void Step3::output_results() const
 * {
 * @endcode
 * 
 * To write the output to a file, we need an object which knows about output
 * formats and the like. This is the DataOut class, and we need an object of
 * that type:
 * 
 * @code
 *   DataOut<2> data_out;
 * @endcode
 * 
 * Now we have to tell it where to take the values from which it shall
 * write. We tell it which DoFHandler object to use, and the solution vector
 * (and the name by which the solution variable shall appear in the output
 * file). If we had more than one vector which we would like to look at in
 * the output (for example right hand sides, errors per cell, etc) we would
 * add them as well:
 * 
 * @code
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 * @endcode
 * 
 * After the DataOut object knows which data it is to work on, we have to
 * tell it to process them into something the back ends can handle. The
 * reason is that we have separated the frontend (which knows about how to
 * treat DoFHandler objects and data vectors) from the back end (which knows
 * many different output formats) and use an intermediate data format to
 * transfer data from the front- to the backend. The data is transformed
 * into this intermediate format by the following function:
 * 
 * @code
 *   data_out.build_patches();
 * 
 * @endcode
 * 
 * Now we have everything in place for the actual output. Just open a file
 * and write the data into it, using VTK format (there are many other
 * functions in the DataOut class we are using here that can write the
 * data in postscript, AVS, GMV, Gnuplot, or some other file
 * formats):
 * 
 * @code
 *   std::ofstream output("solution.vtk");
 *   data_out.write_vtk(output);
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step3run"></a> 
 * <h4>Step3::run</h4>
 * 

 * 
 * Finally, the last function of this class is the main function which calls
 * all the other functions of the <code>Step3</code> class. The order in which
 * this is done resembles the order in which most finite element programs
 * work. Since the names are mostly self-explanatory, there is not much to
 * comment about:
 * 
 * @code
 * void Step3::run()
 * {
 *   make_grid();
 *   setup_system();
 *   assemble_system();
 *   solve();
 *   output_results();
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
 * This is the main function of the program. Since the concept of a
 * main function is mostly a remnant from the pre-object oriented era
 * before C++ programming, it often does not do much more than
 * creating an object of the top-level class and calling its principle
 * function.
 * 

 * 
 * Finally, the first line of the function is used to enable output of
 * some diagnostics that deal.II can generate.  The @p deallog
 * variable (which stands for deal-log, not de-allog) represents a
 * stream to which some parts of the library write output. For
 * example, iterative solvers will generate diagnostics (starting
 * residual, number of solver steps, final residual) as can be seen
 * when running this tutorial program.
 * 

 * 
 * The output of @p deallog can be written to the console, to a file,
 * or both. Both are disabled by default since over the years we have
 * learned that a program should only generate output when a user
 * explicitly asks for it. But this can be changed, and to explain how
 * this can be done, we need to explain how @p deallog works: When
 * individual parts of the library want to log output, they open a
 * "context" or "section" into which this output will be placed. At
 * the end of the part that wants to write output, one exits this
 * section again. Since a function may call another one from within
 * the scope where this output section is open, output may in fact be
 * nested hierarchically into these sections. The LogStream class of
 * which @p deallog is a variable calls each of these sections a
 * "prefix" because all output is printed with this prefix at the left
 * end of the line, with prefixes separated by colons. There is always
 * a default prefix called "DEAL" (a hint at deal.II's history as the
 * successor of a previous library called "DEAL" and from which the
 * LogStream class is one of the few pieces of code that were taken
 * into deal.II).
 * 

 * 
 * By default, @p logstream only outputs lines with zero prefixes --
 * i.e., all output is disabled because the default "DEAL" prefix is
 * always there. But one can set a different maximal number of
 * prefixes for lines that should be output to something larger, and
 * indeed here we set it to two by calling
 * LogStream::depth_console(). This means that for all screen output,
 * a context that has pushed one additional prefix beyond the default
 * "DEAL" is allowed to print its output to the screen ("console"),
 * whereas all further nested sections that would have three or more
 * prefixes active would write to @p deallog, but @p deallog does not
 * forward this output to the screen. Thus, running this example (or
 * looking at the "Results" section), you will see the solver
 * statistics prefixed with "DEAL:CG", which is two prefixes. This is
 * sufficient for the context of the current program, but you will see
 * examples later on (e.g., in step-22) where solvers are nested more
 * deeply and where you may get useful information by setting the
 * depth even higher.
 * 
 * @code
 * int main()
 * {
 *   deallog.depth_console(2);
 * 
 *   Step3 laplace_problem;
 *   laplace_problem.run();
 * 
 *   return 0;
 * }
 * @endcode
examples/step-3/doc/results.dox



<a name="Results"></a><h1>Results</h1>


程序的输出看起来如下。

@code
Number of active cells: 1024
Number of degrees of freedom: 1089
DEAL:cg::Starting value 0.121094
DEAL:cg::Convergence step 48 value 5.33692e-13
@endcode



前两行是我们写给  <code>cout</code>  的内容。最后两行是CG求解器在没有我们的干预下生成的。前两行说明了迭代开始时的残差，而最后一行告诉我们求解器需要47次迭代才能使残差的规范值达到5.3e-13，即低于我们在 "solve "函数中设置的阈值1e-12。我们将在下一个程序中展示如何抑制这种输出，这种输出有时对调试很有用，但往往会使屏幕显示变得混乱。

除了上面显示的输出，该程序还生成了文件 <code>solution.vtk</code> ，该文件为VTK格式，被当今许多可视化程序广泛使用--包括两个重量级的<a href="https://www.llnl.gov/visit">VisIt</a>和<a href="https://www.paraview.org">Paraview</a>，是当今最常使用的程序。

使用VisIt，生成一张像这样的解决方案的图片并不是很困难。   <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-3.solution-3.png" alt="Visualization of the solution of step-3">
    </td>
  </tr>
</table>  它同时显示了解和网格，根据每一点的解的值提升到 $x$  -  $y$ 平面之上。当然，这里的解并不特别令人兴奋，但这是拉普拉斯方程所代表的内容和我们为这个程序选择的右手边 $f(\mathbf x)=1$ 的结果。拉普拉斯方程描述了（在许多其他用途中）受外部（也是垂直）力作用的膜的垂直变形。在目前的例子中，膜的边界被夹在一个没有垂直变化的方形框架上；因此，一个恒定的力密度将直观地导致膜简单地向上隆起--就像上图所示。

VisIt和Paraview都允许玩各种可视化的解决方案。几个视频讲座展示了如何使用这些程序。   @dealiiVideoLectureSeeAlso{11,32} 




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


如果你想用这个程序玩一玩，这里有几个建议。   </p> 

 <ul>   <li>  改变几何图形和网格。在程序中，我们通过使用 <code>GridGenerator::hyper_cube</code> 函数生成了一个方形域和网格。然而， <code>GridGenerator</code> 也有大量的其他函数。试试L形域，环形域，或其他你在那里找到的域。     </li> 

    <li>  改变边界条件。代码使用 Functions::ZeroFunction 函数来生成零边界条件。然而，你可能想用 <code>ConstantFunction&lt;2&gt;(1)</code> 而不是 <code>ZeroFunction&lt;2&gt;()</code> 尝试非零常数边界值，以获得单位Dirichlet边界值。在函数命名空间的文档中描述了更多的奇异函数，你可以选择一个来描述你的特定边界值。     </li> 

    <li>  修改边界条件的类型。目前，发生的情况是，我们在周围使用迪里希特边界值，因为默认情况是所有边界部分的边界指标为零，然后我们告诉 VectorTools::interpolate_boundary_values() 函数，在所有指标为零的边界部分上将边界值插值为零。    <p>  如果我们给边界的部分分配不同的指标，我们可以改变这种行为。例如，在调用 GridGenerator::hyper_cube(): @code
  triangulation.begin_active()->face(0)->set_boundary_id(1);
  @endcode后立即尝试这样做。



  这样做的目的是，首先要求三角剖分返回一个迭代器，指向第一个活动单元。当然，由于这是一个正方形的三角测量的粗略网格，此刻三角测量只有一个单元，而且它是活动的。接下来，我们要求单元格返回它的第一个面的迭代器，然后我们要求面将该面的边界指标重置为1。接下来的事情就是这样。当网格被细化时，子单元的面会继承其父母的边界指示器，也就是说，即使在最细的网格上，广场一侧的面的边界指示器为1。稍后，当我们要插值边界条件时， VectorTools::interpolate_boundary_values() 调用将只为那些边界指标为零的面产生边界值，而对那些具有不同边界指标的面则不予理会。这样做的目的是对前者施加Dirichlet边界条件，而对后者施加同质的Neumann条件（即解的法向导数为零，除非在变分等式的右侧添加额外的条款来处理潜在的非零Neumann条件）。如果你运行该程序，你会看到这一点。

  另一种改变边界指标的方法是根据面中心的笛卡尔坐标来标注边界。   例如，我们可以通过检查单元格中心的y坐标是否在-1和1的公差（这里是1e-12）范围内，将沿上下边界的所有单元格标记为边界指示器1。在调用 GridGenerator::hyper_cube(), 后，像以前一样立即尝试这样做。   @code
  for (auto &face : triangulation.active_face_iterators())
    if (std::fabs(face->center()(1) - (-1.0)) < 1e-12 ||
        std::fabs(face->center()(1) - (1.0)) < 1e-12)
      face->set_boundary_id(1);
  @endcode

  虽然这段代码比以前长了一些，但它对复杂的几何形状很有用，因为它不需要脸部标签的知识。

    <li> 最后一点的一个小变化是像上面那样设置不同的边界值，但随后为边界指标一使用不同的边界值函数。在实践中，你要做的是为边界指标一增加对 <code>interpolate_boundary_values</code> 的第二次调用。   @code
  VectorTools::interpolate_boundary_values(dof_handler,
					   1,
					   ConstantFunction<2>(1.),
					   boundary_values);
  @endcode

  如果你在这个函数的第一个调用之后立即进行这个调用，那么它将把边界指标为1的面的边界值内插到单位值，并将这些内插值与之前计算的边界指标为0的值合并。

    <li>  观察收敛情况。我们将只讨论第7步中规范的计算误差，但很容易检查计算在这里已经收敛了。例如，我们可以在一个点上评估解的值，并比较不同%的全局细化的值（全局细化的步骤数在上面的 <code>LaplaceProblem::make_grid</code> 中设定）。为了评估某个点的解决方案，例如在 $(\frac 13, \frac 13)$ ，我们可以在 <code>LaplaceProblem::output_results</code> 函数中加入以下代码。   @code
    std::cout << "Solution at (1/3,1/3): "
              << VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3))
              << std::endl;
  @endcode

  对于1到9个全局细化步骤，我们就会得到以下的点值序列。     <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$u_h(\frac 13,\frac13)$</th> </tr>
    <tr> <td>1</td> <td>0.166667</td> </tr>
    <tr> <td>2</td> <td>0.227381</td> </tr>
    <tr> <td>3</td> <td>0.237375</td> </tr>
    <tr> <td>4</td> <td>0.240435</td> </tr>
    <tr> <td>5</td> <td>0.241140</td> </tr>
    <tr> <td>6</td> <td>0.241324</td> </tr>
    <tr> <td>7</td> <td>0.241369</td> </tr>
    <tr> <td>8</td> <td>0.241380</td> </tr>
    <tr> <td>9</td> <td>0.241383</td> </tr>
  </table>  通过注意到每两个连续值之间的差异减少了大约4倍，我们可以猜测 "正确 "的值可能是 $u(\frac 13, \frac 13)\approx 0.241384$  。事实上，如果我们假设这是正确的值，我们可以证明上面的序列确实显示了 ${\cal
  O}(h^2)$ 的收敛&mdash；理论上，收敛顺序应该是 ${\cal O}(h^2 |\log h|)$ ，但是领域和网格的对称性可能导致了观察到的更好的收敛顺序。

  这方面的一个小变种是用二次元重复测试。你需要做的就是在构造函数中把有限元的多项式程度设置为2  <code>LaplaceProblem::LaplaceProblem</code>  。

    <li>  平均值的收敛。一个不同的方法是计算解的平均数，以了解解是否真的收敛了（收敛到什么程度&mdash；我们无法判断它是否真的是正确的值！）。为此，在 <code>LaplaceProblem::output_results</code> 中添加以下代码：@code
    std::cout << "Mean value: "
              << VectorTools::compute_mean_value (dof_handler,
						  QGauss<2>(fe.degree + 1),
						  solution,
						  0)
              << std::endl;
  @endcode

  该函数的文档解释了第二和第四个参数的含义，而第一和第三个参数应该是很明显的。再次做同样的研究，我们改变了全局细化步骤的数量，我们得到以下结果。     <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$\int_\Omega u_h(x)\; dx$</th> </tr>
    <tr> <td>0</td> <td>0.09375000</td> </tr>
    <tr> <td>1</td> <td>0.12790179</td> </tr>
    <tr> <td>2</td> <td>0.13733440</td> </tr>
    <tr> <td>3</td> <td>0.13976069</td> </tr>
    <tr> <td>4</td> <td>0.14037251</td> </tr>
    <tr> <td>5</td> <td>0.14052586</td> </tr>
    <tr> <td>6</td> <td>0.14056422</td> </tr>
    <tr> <td>7</td> <td>0.14057382</td> </tr>
    <tr> <td>8</td> <td>0.14057622</td> </tr>
  </table>  同样，两个相邻值之间的差异下降了约四倍，表明收敛为  ${\cal O}(h^2)$  。   </ul> 




<a name="UsingHDF5tooutputthesolutionandadditionaldata"></a><h3>Using %HDF5 to output the solution and additional data</h3>


%HDF5是一种常用的格式，可以被许多脚本语言（如R或Python）读取。让deal.II产生一些%HDF5文件并不困难，然后可以在外部脚本中使用，对该程序产生的一些数据进行后处理。这里有一些关于可能的想法。




<a name="Changingtheoutputtoh5"></a><h4> Changing the output to .h5</h4>


为了充分利用自动化，我们首先需要为全局细化步骤的数量引入一个私有变量 <code>unsigned int n_refinement_steps </code> ，它将被用于输出文件名。在 <code>make_grid()</code> we then replace <code>triangulation.refine_global(5);</code> 中用

@code
n_refinement_steps = 5;
triangulation.refine_global(n_refinement_steps);
@endcode

deal.II库有两个不同的%HDF5绑定，一个在HDF5命名空间（用于对接通用数据文件），另一个在DataOut（专门用于为解决方案的可视化写文件）。尽管HDF5 deal.II绑定支持串行和MPI，但%HDF5 DataOut绑定只支持并行输出。由于这个原因，我们需要初始化一个只有一个处理器的MPI通信器。这可以通过添加以下代码来实现。

@code
int main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  ...
}
@endcode

接下来我们改变 `Step3::output_results()` 的输出例程，如DataOutBase命名空间文档中所述。

@code
const std::string filename_h5 = "solution_" + std::to_string(n_refinement_steps) + ".h5";
DataOutBase::DataOutFilterFlags flags(true, true);
DataOutBase::DataOutFilter data_filter(flags);
data_out.write_filtered_data(data_filter);
data_out.write_hdf5_parallel(data_filter, filename_h5, MPI_COMM_WORLD);
@endcode

然后，产生的文件可以被可视化，就像教程的原始版本产生的VTK文件一样；但是，由于%HDF5是一种更通用的文件格式，它也可以很容易地用脚本语言处理，用于其他目的。




<a name="Addingthepointvalueandthemeanseeextensionaboveintotheh5file"></a><h4> Adding the point value and the mean (see extension above) into the .h5 file</h4>


在输出解决方案后，可以再次打开该文件以包括更多的数据集。  这使得我们可以将实验的所有必要信息保存在一个结果文件中，然后可以由一些后处理脚本来读取和处理。关于可能的输出选项，请看 HDF5::Group::write_dataset() 的进一步信息）。

为了实现这一点，我们首先将必要的头文件纳入我们的文件。

@code
#include <deal.II/base/hdf5.h>
@endcode

在我们的输出例程的末尾添加以下几行，将关于某一点的解的值，以及解的平均值的信息添加到我们的%HDF5文件中。

@code
HDF5::File data_file(filename_h5, HDF5::File::FileAccessMode::open, MPI_COMM_WORLD);
Vector<double> point_value(1);
point_value[0] = VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3));
data_file.write_dataset("point_value", point_value);
Vector<double> mean_value(1);
mean_value[0] = VectorTools::compute_mean_value(dof_handler,
                                                QGauss<2>(fe.degree + 1),
                                                solution, 0);
data_file.write_dataset("mean_value",mean_value);
@endcode






<a name="UsingRandggplot2togenerateplots"></a><h3> Using R and ggplot2 to generate plots</h3>


上述放入%HDF5文件的数据，然后可以从脚本语言中使用，进行进一步的后处理。在下文中，让我们展示一下，特别是如何用<a href="https://en.wikipedia.org/wiki/R_(programming_language)">R
programming language</a>这个在统计数据分析中广泛使用的语言来完成。(例如，类似的事情也可以在Python中完成。)如果你不熟悉R和ggplot2，你可以看看R的数据木工课程<a href="https://datacarpentry.org/R-ecology-lesson/index.html">here</a>。此外，由于大多数搜索引擎对 "R+主题 "这种形式的搜索很吃力，我们建议使用专门的服务<a
href="http://rseek.org">RSeek </a>来代替。

R和其他语言最突出的区别是，赋值运算符（`a = 5`）通常被写成`a <- 5`。由于后者被认为是标准的，我们将在我们的例子中也使用它。要在R语言中打开`.h5`文件，你必须安装<a href="https://bioconductor.org/packages/release/bioc/html/rhdf5.html">rhdf5</a>包，它是Bioconductor软件包的一部分。

首先，我们将包括所有必要的包，并看看我们文件中的数据是如何结构化的。

@code{.r}
library(rhdf5)     # library for handling HDF5 files
library(ggplot2)   # main plotting library
library(grDevices) # needed for output to PDF
library(viridis)   # contains good colormaps for sequential data


refinement <- 5
h5f <- H5Fopen(paste("solution_",refinement,".h5",sep=""))
print(h5f)
@endcode

这给出了以下输出

@code{.unparsed}
HDF5 FILE
   name /
filename


    name       otype  dclass     dim
0 cells       H5I_DATASET INTEGER  x 1024
1 mean_value  H5I_DATASET FLOAT   1
2 nodes       H5I_DATASET FLOAT    x 1089
3 point_value H5I_DATASET FLOAT   1
4 solution    H5I_DATASET FLOAT    x 1089
@endcode

数据集可以通过  <code>h5f\$name</code>  访问。函数  <code>dim(h5f\$cells)</code>  给我们提供了用于存储我们单元格的矩阵的尺寸。我们可以看到以下三个矩阵，以及我们添加的两个额外数据点。   <ul>   <li>   <code>cells</code>  ：一个4x1024的矩阵，存储每个单元的（C++）顶点指数  <li>   <code>nodes</code>  ：一个2x1089的矩阵，存储我们单元顶点的位置值（x，y）  <li>   <code>solution</code>  : 一个1x1089的矩阵，存储我们的解决方案在每个顶点的值  </ul>  现在我们可以使用这些数据来生成各种图表。用ggplot2作图通常分为两步。首先，数据需要被处理并添加到一个  <code>data.frame</code>  。之后，构建一个 <code>ggplot</code> 对象，并通过向其添加绘图元素来进行操作。

 <code>nodes</code> and <code>cells</code> 包含我们绘制网格所需的所有信息。下面的代码将所有的数据打包成一个数据框架，用于绘制我们的网格。

@code{.r}
# Counting in R starts at 1 instead of 0, so we need to increment all
# vertex indices by one:
cell_ids <- h5f$cells+1


# Store the x and y positions of each vertex in one big vector in a
# cell by cell fashion (every 4 entries belong to one cell):
cells_x <- h5f$nodes[1,][cell_ids]
cells_y <- h5f$nodes[2,][cell_ids]


# Construct a vector that stores the matching cell by cell grouping
# (1,1,1,1,2,2,2,2,...):
groups <- rep(1:ncol(cell_ids),each=4)


# Finally put everything into one dataframe:
meshdata <- data.frame(x = cells_x, y = cells_y, id = groups)
@endcode



有了完成的数据框架，我们就有了绘制网格所需的一切。

@code{.r}
pdf (paste("grid_",refinement,".pdf",sep=""),width = 5,height = 5) # Open new PDF file
plt <- ggplot(meshdata,aes(x=x,y=y,group=id))                      # Construction of our plot
                                                                   # object, at first only data


plt <- plt + geom_polygon(fill="white",colour="black")             # Actual plotting of the grid as polygons
plt <- plt + ggtitle(paste("grid at refinement level #",refinement))


print(plt)                                                         # Show the current state of the plot/add it to the pdf
dev.off()                                                          # Close PDF file
@endcode



这个文件的内容看起来如下（不是很令人兴奋，但你会明白的）。   <table width="60%" align="center">
  <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.grid_5.png" alt="Grid after 5 refinement steps of step-3">
   </td>
  </tr>
</table> 

我们还可以将解决方案本身可视化，这看起来会更有趣。为了给我们的解决方案做一个二维伪色图，我们将使用  <code>geom_raster</code>  。这个函数需要一个结构化的网格，即在x和y方向上是均匀的。幸运的是，我们在这一点上的数据是以正确的方式结构化的。下面的代码将我们的曲面的伪彩色表示法绘制成一个新的PDF。

@code{.r}
pdf (paste("pseudocolor_",refinement,".pdf",sep=""),width = 5,height = 4.2) # Open new PDF file
colordata <- data.frame(x = h5f$nodes[1,],y = h5f$nodes[2,] , solution = h5f$solution[1,])
plt <- ggplot(colordata,aes(x=x,y=y,fill=solution))
plt <- plt + geom_raster(interpolate=TRUE)
plt <- plt + scale_fill_viridis()
plt <- plt + ggtitle(paste("solution at refinement level #",refinement))


print(plt)
dev.off()
H5Fclose(h5f) # Close the HDF5 file
@endcode

现在的情况是这样的。   <table width="60%" align="center">
 <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.pseudocolor_5.png" alt="Solution after 5 refinement steps of step-3">
   </td>
 </tr>
</table> 

为了绘制收敛曲线，我们需要从1开始用不同的 <code>n_refinement_steps</code> 值多次重新运行C++代码。由于每个文件只包含一个数据点，我们需要对它们进行循环，并将结果串联成一个矢量。

@code{.r}
n_ref <- 8   # Maximum refinement level for which results are existing


# First we initiate all vectors with the results of the first level
h5f   <- H5Fopen("solution_1.h5")
dofs  <- dim(h5f$solution)[2]
mean  <- h5f$mean_value
point <- h5f$point_value
H5Fclose(h5f)


for (reflevel in 2:n_ref)
{
   h5f   <- H5Fopen(paste("solution_",reflevel,".h5",sep=""))
   dofs  <- c(dofs,dim(h5f\$solution)[2])
   mean  <- c(mean,h5f\$mean_value)
   point <- c(point,h5f\$point_value)
   H5Fclose(h5f)
}
@endcode

由于我们对数值本身不感兴趣，而是对与 "精确 "解决方案相比的误差感兴趣，我们将假设我们的最高细化水平是该解决方案，并从数据中省略它。

@code{.r}
# Calculate the error w.r.t. our maximum refinement step
mean_error  <- abs(mean[1:n_ref-1]-mean[n_ref])
point_error <- abs(point[1:n_ref-1]-point[n_ref])


# Remove the highest value from our DoF data
dofs     <- dofs[1:n_ref-1]
convdata <- data.frame(dofs = dofs, mean_value= mean_error, point_value = point_error)
@endcode

现在我们有所有的数据可以用来生成我们的图。在对数尺度上绘制误差往往是有用的，这在下面的代码中可以实现。

@code
pdf (paste("convergence.pdf",sep=""),width = 5,height = 4.2)
plt <- ggplot(convdata,mapping=aes(x = dofs, y = mean_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "mean value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)


plt <- ggplot(convdata,mapping=aes(x = dofs, y = point_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "point value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)


dev.off()
@endcode

这就产生了下面的图，显示了均值和所选点的解值的误差如何很好地收敛到零。   <table style="width:50%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_mean.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_point.png" alt=""></td>
  </tr>
</table> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-3.cc"
*/
