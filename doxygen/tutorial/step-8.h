/**
@page step_8 The step-8 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeElasticProblemcodeclasstemplate">The <code>ElasticProblem</code> class template</a>
        <li><a href="#Righthandsidevalues">Right hand side values</a>
        <li><a href="#ThecodeElasticProblemcodeclassimplementation">The <code>ElasticProblem</code> class implementation</a>
      <ul>
        <li><a href="#ElasticProblemElasticProblemconstructor">ElasticProblem::ElasticProblem constructor</a>
        <li><a href="#ElasticProblemsetup_system">ElasticProblem::setup_system</a>
        <li><a href="#ElasticProblemassemble_system">ElasticProblem::assemble_system</a>
        <li><a href="#ElasticProblemsolve">ElasticProblem::solve</a>
        <li><a href="#ElasticProblemrefine_grid">ElasticProblem::refine_grid</a>
        <li><a href="#ElasticProblemoutput_results">ElasticProblem::output_results</a>
        <li><a href="#ElasticProblemrun">ElasticProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-8/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



在现实生活中，大多数偏微分方程实际上是方程组。相应地，解通常是矢量值的。deal.II库支持这样的问题（见 @ref vector_valued 模块中的大量文档），我们将表明这大多是相当简单的。唯一比较复杂的问题是在组装矩阵和右手边，但这些也很容易理解。

 @dealiiVideoLecture{19} 

在这个教程程序中，我们将想解决<a href="https://en.wikipedia.org/wiki/Linear_elasticity">elastic equations</a>。它们是对拉普拉斯方程的扩展，有一个矢量值的解，描述了受力的刚体在每个空间方向的位移。当然，力也是矢量值的，意味着在每一个点上它都有一个方向和一个绝对值。

人们可以用多种方式来写弹性方程。以最明显的方式显示与拉普拉斯方程的对称性的是将其写成

@f[


  -
  \text{div}\,
  ({\mathbf C} \nabla \mathbf{u})
  =
  \mathbf f,


@f]

其中 $\mathbf u$ 是每一点的矢量值位移， $\mathbf f$ 是力， ${\mathbf C}$ 是一个等级4的张量（即它有四个指数），编码应力-应变关系--本质上，它代表胡克斯定律中的<a href="https://en.wikipedia.org/wiki/Hooke%27s_law">"spring constant"</a>，将位移与力联系起来。  在许多情况下，如果我们想要模拟的物体的变形是由不同的材料组成的，那么 ${\mathbf C}$ 将取决于 $\mathbf x$ 。

虽然上述方程的形式是正确的，但这并不是它们通常的推导方式。事实上，位移的梯度 $\nabla\mathbf u$ （一个矩阵）没有物理意义，而其对称版本。

@f[
\varepsilon(\mathbf u)_{kl} =\frac{1}{2}(\partial_k u_l + \partial_l u_k),


@f]

做，通常被称为 "应变"。(在这里和下文中， $\partial_k=\frac{\partial}{\partial x_k}$  。我们还将使用<a href="https://en.wikipedia.org/wiki/Einstein_notation">Einstein summation
convention</a>，即只要同一指数在方程式中出现两次，就意味着对该指数进行求和；但是，我们将不区分上下指数）。)有了这个应变的定义，弹性方程就读作

@f[


  -
  \text{div}\,
  ({\mathbf C} \varepsilon(\mathbf u))
  =
  \mathbf f,


@f]

你可以把它看作是拉普拉斯方程对矢量值问题的更自然的概括。(首先显示的形式等同于这种形式，因为张量 ${\mathbf C}$ 具有某些对称性，即 $C_{ijkl}=C_{ijlk}$  ，因此 ${\mathbf C} \varepsilon(\mathbf u)_{kl}
= {\mathbf C} \nabla\mathbf u$  。)

当然，我们也可以把这些方程写成组件形式。

@f[


  -
  \partial_j (c_{ijkl} \varepsilon_{kl})
  =
  f_i,
  \qquad
  i=1\ldots d.


@f]



在许多情况下，我们知道所考虑的材料是各向同性的，在这种情况下，通过引入两个系数 $\lambda$ 和 $\mu$ ，系数张量减少为

@f[
  c_{ijkl}
  =
  \lambda \delta_{ij} \delta_{kl} +
  \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).


@f]



然后，弹性方程可以用更简单的形式重写。

@f[


   -
   \nabla \lambda (\nabla\cdot {\mathbf u})


   -
   (\nabla \cdot \mu \nabla) {\mathbf u}


   -
   \nabla\cdot \mu (\nabla {\mathbf u})^T
   =
   {\mathbf f},


@f]

而各自的双线性形式则是

@f[
  a({\mathbf u}, {\mathbf v}) =
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega,


@f]

或将第一项写成成分之和。

@f[
  a({\mathbf u}, {\mathbf v}) =
  \sum_{k,l}
  \left(
    \lambda \partial_l u_l, \partial_k v_k
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega.


@f]



 @note 按照写法，如果位移很小，我们可以假设<a
href="http://en.wikipedia.org/wiki/Hookes_law">Hooke's law</a>是有效的，上面的方程一般被认为是对三维物体位移的正确描述。在这种情况下，上面的指数 $i,j,k,l$ 都是在集合 $\{1,2,3\}$ 上运行的（或者，在C++源中，在 $\{0,1,2\}$ 上运行）。然而，按照目前的情况，程序是在2d中运行的，虽然上面的方程在这种情况下也有数学意义，但它们只能描述一个真正的二维实体。特别是，它们不是对 $x-y$ 方向上无限大的体的横截面的适当描述；这与其他许多二维方程相反，这些方程可以通过假设体在 $z$ -方向上具有无限大的范围和解函数不依赖于 $z$ 坐标来获得。另一方面，也有二维弹性模型的方程；例如，见维基百科上的<a
href="http://en.wikipedia.org/wiki/Infinitesimal_strain_theory#Special_cases">plane
strain</a>、<a
href="http://en.wikipedia.org/wiki/Antiplane_shear">antiplane shear</a>和<a
href="http://en.wikipedia.org/wiki/Plane_stress#Plane_stress">plan stress</a>文章。

但让我们回到最初的问题上。我们如何为这样一个方程组装矩阵？在 @ref vector_valued 模块的文档中给出了一个很长的答案，其中有许多不同的选择。从历史上看，下面所示的解决方案是该库早期唯一可用的解决方案。事实证明，它也是最快的。另一方面，如果百分之几的计算时间并不重要，还有比下面讨论的更简单、更直观的方法来组装线性系统，但这些方法直到本教程首次编写后的几年才可用；如果你对它们感兴趣，可以看看  @ref vector_valued  模块。

让我们回到如何组装线性系统的问题上来。首先我们需要一些关于形状函数在矢量值有限元情况下如何工作的知识。基本上，这归结为以下几点：让 $n$ 为我们建立矢量元素的标量有限元素的形状函数的数量（例如，我们将对矢量值有限元素的每个分量使用双线性函数，所以标量有限元素是我们在以前的例子中已经使用过的 <code>FE_Q(1)</code> 元素，以及两个空间维度的 $n=4$ ）。此外，让 $N$ 为矢量元素的形状函数数量；在两个空间维度中，我们需要为矢量的每个分量提供 $n$ 个形状函数，因此 $N=2n$  。那么，矢量元素的 $i$ 个形状函数的形式为

@f[
  \Phi_i({\mathbf x}) = \varphi_{\text{base}(i)}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)},


@f]

其中 $e_l$ 是第 $l$ 个单位向量， $\text{comp}(i)$ 是告诉我们 $\Phi_i$ 的哪个分量是不为零的函数（对于每个向量形状函数，只有一个分量是不为零的，其他都是零）。   $\varphi_{\text{base}(i)}(x)$ 描述了形状函数的空间依赖性，它被认为是标量元素的第 $\text{base}(i)$ 个形状函数。当然，虽然 $i$ 的范围是 $0,\ldots,N-1$ ，但函数 $\text{comp}(i)$ 和 $\text{base}(i)$ 的范围分别为 $0,1$ （在二维）和 $0,\ldots,n-1$ 。

例如（尽管这种形状函数的顺序不被保证，你也不应该依赖它），下面的布局可以被库使用。

@f{eqnarray*}
  \Phi_0({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_0({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_1({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_0({\mathbf x})
  \end{array}\right),
  \\
  \Phi_2({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_1({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_3({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_1({\mathbf x})
  \end{array}\right),
  \ldots


@f}

在这里

@f[
  \text{comp}(0)=0, \quad  \text{comp}(1)=1, \quad  \text{comp}(2)=0, \quad  \text{comp}(3)=1, \quad  \ldots


@f]



@f[
  \text{base}(0)=0, \quad  \text{base}(1)=0, \quad  \text{base}(2)=1, \quad  \text{base}(3)=1, \quad  \ldots


@f]



除了非常罕见的情况，你不需要知道标量元素的哪个形状函数 $\varphi_{\text{base}(i)}$ 属于矢量元素的一个形状函数 $\Phi_i$ 。因此，让我们定义

@f[
  \phi_i = \varphi_{\text{base}(i)}


@f]

据此，我们可以将矢量形状函数写为

@f[
  \Phi_i({\mathbf x}) = \phi_{i}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)}.


@f]

现在你可以安全地忘记函数 $\text{base}(i)$ 了，至少在这个例子程序的其余部分。

现在使用这个矢量形状函数，我们可以将离散的有限元解写为

@f[
  {\mathbf u}_h({\mathbf x}) =
  \sum_i \Phi_i({\mathbf x})\ U_i


@f]

具有标量系数  $U_i$  。如果我们定义一个模拟函数 ${\mathbf v}_h$ 作为测试函数，我们可以将离散问题写成如下。找出系数 $U_i$ ，使得

@f[
  a({\mathbf u}_h, {\mathbf v}_h) = ({\mathbf f}, {\mathbf v}_h)
  \qquad
  \forall {\mathbf v}_h.


@f]



如果我们把双线性形式的定义和 ${\mathbf u}_h$ 和 ${\mathbf v}_h$ 的表示插入这个公式。

@f{eqnarray*}
  \sum_{i,j}
    U_i V_j
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_\Omega
  \right\}
\\
=
  \sum_j V_j
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_\Omega.


@f}

我们注意到，在这里和下文中，指数 $k,l$ 在空间方向上运行，即 $0\le k,l < d$  ，而指数 $i,j$ 在自由度上运行。

因此，单元 $K$ 上的局部刚度矩阵有以下条目。

@f[
  A^K_{ij}
  =
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_K
  \right\},


@f]

其中 $i,j$ 现在是局部自由度，因此 $0\le i,j < N$  。在这些公式中，我们总是取矢量形状函数 $\Phi_i$ 的一些分量，当然，这些分量是如下给出的（见其定义）。

@f[
  (\Phi_i)_l = \phi_i \delta_{l,\text{comp}(i)},


@f]

与克朗克符号  $\delta_{nm}$  。由于这一点，我们可以删除一些对  $k$  和  $l$  的和。

@f{eqnarray*}
  A^K_{ij}
  &=&
  \sum_{k,l}
  \Bigl\{
  \left(
    \lambda \partial_l \phi_i\ \delta_{l,\text{comp}(i)},
            \partial_k \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
\\
  &\qquad\qquad& +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_l \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
  +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_k \phi_j\ \delta_{l,\text{comp}(j)}
  \right)_K
  \Bigr\}
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \sum_l
  \left(
    \mu \partial_l \phi_i,
        \partial_l \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \left(
    \mu \nabla \phi_i,
        \nabla \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K.


@f}



同样地，单元格 $K$ 对右侧向量的贡献是

@f{eqnarray*}
  f^K_j
  &=&
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_K
\\
  &=&
  \sum_l
  \left(
    f_l,
    \phi_j \delta_{l,\text{comp}(j)}
  \right)_K
\\
  &=&
  \left(
    f_{\text{comp}(j)},
    \phi_j
  \right)_K.


@f}



这就是我们要实现局部刚度矩阵和右手边向量的形式。

作为最后的说明：在第17步的例子程序中，我们将重新审视这里提出的弹性问题，并将展示如何在一个计算机集群上以%并行的方式解决这个问题。因此，所产生的程序将能够以更高的精度解决这个问题，而且如果需要的话，效率更高。此外，在第20步， @ref step_21 "第21步"，以及其他一些后来的教程程序中，我们将重新审视一些矢量值问题，并展示一些技术，这些技术可能使实际通过上面显示的所有东西更简单，与 FiniteElement::system_to_component_index 等。


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
 * As usual, the first few include files are already known, so we will not
 * comment on them further.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/tensor.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * In this example, we need vector-valued finite elements. The support for
 * these can be found in the following include file:
 * 
 * @code
 * #include <deal.II/fe/fe_system.h>
 * @endcode
 * 
 * We will compose the vector-valued finite elements from regular Q1 elements
 * which can be found here, as usual:
 * 
 * @code
 * #include <deal.II/fe/fe_q.h>
 * 
 * @endcode
 * 
 * This again is C++:
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The last step is as in previous programs. In particular, just like in
 * step-7, we pack everything that's specific to this program into a namespace
 * of its own.
 * 
 * @code
 * namespace Step8
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclasstemplate"></a> 
 * <h3>The <code>ElasticProblem</code> class template</h3>
 * 

 * 
 * The main class is, except for its name, almost unchanged with respect to
 * the step-6 example.
 *   

 * 
 * The only change is the use of a different class for the <code>fe</code>
 * variable: Instead of a concrete finite element class such as FE_Q, we now
 * use a more generic one, FESystem. In fact, FESystem is not really a
 * finite element itself in that it does not implement shape functions of
 * its own. Rather, it is a class that can be used to stack several other
 * elements together to form one vector-valued finite element. In our case,
 * we will compose the vector-valued element of <code>FE_Q(1)</code>
 * objects, as shown below in the constructor of this class.
 * 
 * @code
 *   template <int dim>
 *   class ElasticProblem
 *   {
 *   public:
 *     ElasticProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system();
 *     void solve();
 *     void refine_grid();
 *     void output_results(const unsigned int cycle) const;
 * 
 *     Triangulation<dim> triangulation;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     FESystem<dim> fe;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsidevalues"></a> 
 * <h3>Right hand side values</h3>
 * 

 * 
 * Before going over to the implementation of the main class, we declare and
 * define the function which describes the right hand side. This time, the
 * right hand side is vector-valued, as is the solution, so we will describe
 * the changes required for this in some more detail.
 *   

 * 
 * To prevent cases where the return vector has not previously been set to
 * the right size we test for this case and otherwise throw an exception at
 * the beginning of the function. Note that enforcing that output arguments
 * already have the correct size is a convention in deal.II, and enforced
 * almost everywhere. The reason is that we would otherwise have to check at
 * the beginning of the function and possibly change the size of the output
 * vector. This is expensive, and would almost always be unnecessary (the
 * first call to the function would set the vector to the right size, and
 * subsequent calls would only have to do redundant checks). In addition,
 * checking and possibly resizing the vector is an operation that can not be
 * removed if we can't rely on the assumption that the vector already has
 * the correct size; this is in contract to the Assert call that is
 * completely removed if the program is compiled in optimized mode.
 *   

 * 
 * Likewise, if by some accident someone tried to compile and run the
 * program in only one space dimension (in which the elastic equations do
 * not make much sense since they reduce to the ordinary Laplace equation),
 * we terminate the program in the second assertion. The program will work
 * just fine in 3d, however.
 * 
 * @code
 *   template <int dim>
 *   void right_hand_side(const std::vector<Point<dim>> &points,
 *                        std::vector<Tensor<1, dim>> &  values)
 *   {
 *     Assert(values.size() == points.size(),
 *            ExcDimensionMismatch(values.size(), points.size()));
 *     Assert(dim >= 2, ExcNotImplemented());
 * 
 * @endcode
 * 
 * The rest of the function implements computing force values. We will use
 * a constant (unit) force in x-direction located in two little circles
 * (or spheres, in 3d) around points (0.5,0) and (-0.5,0), and y-force in
 * an area around the origin; in 3d, the z-component of these centers is
 * zero as well.
 *     

 * 
 * For this, let us first define two objects that denote the centers of
 * these areas. Note that upon construction of the Point objects, all
 * components are set to zero.
 * 
 * @code
 *     Point<dim> point_1, point_2;
 *     point_1(0) = 0.5;
 *     point_2(0) = -0.5;
 * 
 *     for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
 *       {
 * @endcode
 * 
 * If <code>points[point_n]</code> is in a circle (sphere) of radius
 * 0.2 around one of these points, then set the force in x-direction
 * to one, otherwise to zero:
 * 
 * @code
 *         if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
 *             ((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
 *           values[point_n][0] = 1.0;
 *         else
 *           values[point_n][0] = 0.0;
 * 
 * @endcode
 * 
 * Likewise, if <code>points[point_n]</code> is in the vicinity of the
 * origin, then set the y-force to one, otherwise to zero:
 * 
 * @code
 *         if (points[point_n].norm_square() < 0.2 * 0.2)
 *           values[point_n][1] = 1.0;
 *         else
 *           values[point_n][1] = 0.0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclassimplementation"></a> 
 * <h3>The <code>ElasticProblem</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="ElasticProblemElasticProblemconstructor"></a> 
 * <h4>ElasticProblem::ElasticProblem constructor</h4>
 * 

 * 
 * Following is the constructor of the main class. As said before, we would
 * like to construct a vector-valued finite element that is composed of
 * several scalar finite elements (i.e., we want to build the vector-valued
 * element so that each of its vector components consists of the shape
 * functions of a scalar element). Of course, the number of scalar finite
 * elements we would like to stack together equals the number of components
 * the solution function has, which is <code>dim</code> since we consider
 * displacement in each space direction. The FESystem class can handle this:
 * we pass it the finite element of which we would like to compose the
 * system of, and how often it shall be repeated:
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   ElasticProblem<dim>::ElasticProblem()
 *     : dof_handler(triangulation)
 *     , fe(FE_Q<dim>(1), dim)
 *   {}
 * @endcode
 * 
 * In fact, the FESystem class has several more constructors which can
 * perform more complex operations than just stacking together several
 * scalar finite elements of the same type into one; we will get to know
 * these possibilities in later examples.
 * 

 * 
 * 

 * 
 * 
 * <a name="ElasticProblemsetup_system"></a> 
 * <h4>ElasticProblem::setup_system</h4>
 * 

 * 
 * Setting up the system of equations is identical to the function used in
 * the step-6 example. The DoFHandler class and all other classes used here
 * are fully aware that the finite element we want to use is vector-valued,
 * and take care of the vector-valuedness of the finite element
 * themselves. (In fact, they do not, but this does not need to bother you:
 * since they only need to know how many degrees of freedom there are per
 * vertex, line and cell, and they do not ask what they represent,
 * i.e. whether the finite element under consideration is vector-valued or
 * whether it is, for example, a scalar Hermite element with several degrees
 * of freedom on each vertex).
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              0,
 *                                              Functions::ZeroFunction<dim>(dim),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ false);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemassemble_system"></a> 
 * <h4>ElasticProblem::assemble_system</h4>
 * 

 * 
 * The big changes in this program are in the creation of matrix and right
 * hand side, since they are problem-dependent. We will go through that
 * process step-by-step, since it is a bit more complicated than in previous
 * examples.
 *   

 * 
 * The first parts of this function are the same as before, however: setting
 * up a suitable quadrature formula, initializing an FEValues object for the
 * (vector-valued) finite element we use as well as the quadrature object,
 * and declaring a number of auxiliary arrays. In addition, we declare the
 * ever same two abbreviations: <code>n_q_points</code> and
 * <code>dofs_per_cell</code>. The number of degrees of freedom per cell we
 * now obviously ask from the composed finite element rather than from the
 * underlying scalar Q1 element. Here, it is <code>dim</code> times the
 * number of degrees of freedom per cell of the Q1 element, though this is
 * not explicit knowledge we need to care about:
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::assemble_system()
 *   {
 *     QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * As was shown in previous examples as well, we need a place where to
 * store the values of the coefficients at all the quadrature points on a
 * cell. In the present situation, we have two coefficients, lambda and
 * mu.
 * 
 * @code
 *     std::vector<double> lambda_values(n_q_points);
 *     std::vector<double> mu_values(n_q_points);
 * 
 * @endcode
 * 
 * Well, we could as well have omitted the above two arrays since we will
 * use constant coefficients for both lambda and mu, which can be declared
 * like this. They both represent functions always returning the constant
 * value 1.0. Although we could omit the respective factors in the
 * assemblage of the matrix, we use them here for purpose of
 * demonstration.
 * 
 * @code
 *     Functions::ConstantFunction<dim> lambda(1.), mu(1.);
 * 
 * @endcode
 * 
 * Like the two constant functions above, we will call the function
 * right_hand_side just once per cell to make things simpler.
 * 
 * @code
 *     std::vector<Tensor<1, dim>> rhs_values(n_q_points);
 * 
 * @endcode
 * 
 * Now we can begin with the loop over all cells:
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * Next we get the values of the coefficients at the quadrature
 * points. Likewise for the right hand side:
 * 
 * @code
 *         lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
 *         mu.value_list(fe_values.get_quadrature_points(), mu_values);
 *         right_hand_side(fe_values.get_quadrature_points(), rhs_values);
 * 
 * @endcode
 * 
 * Then assemble the entries of the local stiffness matrix and right
 * hand side vector. This follows almost one-to-one the pattern
 * described in the introduction of this example.  One of the few
 * comments in place is that we can compute the number
 * <code>comp(i)</code>, i.e. the index of the only nonzero vector
 * component of shape function <code>i</code> using the
 * <code>fe.system_to_component_index(i).first</code> function call
 * below.
 *         

 * 
 * (By accessing the <code>first</code> variable of the return value
 * of the <code>system_to_component_index</code> function, you might
 * already have guessed that there is more in it. In fact, the
 * function returns a <code>std::pair@<unsigned int, unsigned
 * int@></code>, of which the first element is <code>comp(i)</code>
 * and the second is the value <code>base(i)</code> also noted in the
 * introduction, i.e.  the index of this shape function within all the
 * shape functions that are nonzero in this component,
 * i.e. <code>base(i)</code> in the diction of the introduction. This
 * is not a number that we are usually interested in, however.)
 *         

 * 
 * With this knowledge, we can assemble the local matrix
 * contributions:
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             const unsigned int component_i =
 *               fe.system_to_component_index(i).first;
 * 
 *             for (const unsigned int j : fe_values.dof_indices())
 *               {
 *                 const unsigned int component_j =
 *                   fe.system_to_component_index(j).first;
 * 
 *                 for (const unsigned int q_point :
 *                      fe_values.quadrature_point_indices())
 *                   {
 *                     cell_matrix(i, j) +=
 * @endcode
 * 
 * The first term is $\lambda \partial_i u_i, \partial_j
 * v_j) + (\mu \partial_i u_j, \partial_j v_i)$. Note
 * that <code>shape_grad(i,q_point)</code> returns the
 * gradient of the only nonzero component of the i-th
 * shape function at quadrature point q_point. The
 * component <code>comp(i)</code> of the gradient, which
 * is the derivative of this only nonzero vector
 * component of the i-th shape function with respect to
 * the comp(i)th coordinate is accessed by the appended
 * brackets.
 * 
 * @code
 *                       (                                                  
 *                         (fe_values.shape_grad(i, q_point)[component_i] * 
 *                          fe_values.shape_grad(j, q_point)[component_j] * 
 *                          lambda_values[q_point])                         
 *                         +                                                
 *                         (fe_values.shape_grad(i, q_point)[component_j] * 
 *                          fe_values.shape_grad(j, q_point)[component_i] * 
 *                          mu_values[q_point])                             
 *                         +                                                
 * @endcode
 * 
 * The second term is $(\mu \nabla u_i, \nabla
 * v_j)$. We need not access a specific component of
 * the gradient, since we only have to compute the
 * scalar product of the two gradients, of which an
 * overloaded version of <tt>operator*</tt> takes
 * care, as in previous examples.
 *                         

 * 
 * Note that by using the <tt>?:</tt> operator, we only
 * do this if <tt>component_i</tt> equals
 * <tt>component_j</tt>, otherwise a zero is added
 * (which will be optimized away by the compiler).
 * 
 * @code
 *                         ((component_i == component_j) ?        
 *                            (fe_values.shape_grad(i, q_point) * 
 *                             fe_values.shape_grad(j, q_point) * 
 *                             mu_values[q_point]) :              
 *                            0)                                  
 *                         ) *                                    
 *                       fe_values.JxW(q_point);                  
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * Assembling the right hand side is also just as discussed in the
 * introduction:
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             const unsigned int component_i =
 *               fe.system_to_component_index(i).first;
 * 
 *             for (const unsigned int q_point :
 *                  fe_values.quadrature_point_indices())
 *               cell_rhs(i) += fe_values.shape_value(i, q_point) *
 *                              rhs_values[q_point][component_i] *
 *                              fe_values.JxW(q_point);
 *           }
 * 
 * @endcode
 * 
 * The transfer from local degrees of freedom into the global matrix
 * and right hand side vector does not depend on the equation under
 * consideration, and is thus the same as in all previous
 * examples.
 * 
 * @code
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
 * <a name="ElasticProblemsolve"></a> 
 * <h4>ElasticProblem::solve</h4>
 * 

 * 
 * The solver does not care about where the system of equations comes, as
 * long as it stays positive definite and symmetric (which are the
 * requirements for the use of the CG solver), which the system indeed
 * is. Therefore, we need not change anything.
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
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
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrefine_grid"></a> 
 * <h4>ElasticProblem::refine_grid</h4>
 * 

 * 
 * The function that does the refinement of the grid is the same as in the
 * step-6 example. The quadrature formula is adapted to the linear elements
 * again. Note that the error estimator by default adds up the estimated
 * obtained from all components of the finite element solution, i.e., it
 * uses the displacement in all directions with the same weight. If we would
 * like the grid to be adapted to the x-displacement only, we could pass the
 * function an additional parameter which tells it to do so and do not
 * consider the displacements in all other directions for the error
 * indicators. However, for the current problem, it seems appropriate to
 * consider all displacement components with equal weight.
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::refine_grid()
 *   {
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 * 
 *     KellyErrorEstimator<dim>::estimate(dof_handler,
 *                                        QGauss<dim - 1>(fe.degree + 1),
 *                                        {},
 *                                        solution,
 *                                        estimated_error_per_cell);
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation,
 *                                                     estimated_error_per_cell,
 *                                                     0.3,
 *                                                     0.03);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemoutput_results"></a> 
 * <h4>ElasticProblem::output_results</h4>
 * 

 * 
 * The output happens mostly as has been shown in previous examples
 * already. The only difference is that the solution function is vector
 * valued. The DataOut class takes care of this automatically, but we have
 * to give each component of the solution vector a different name.
 *   

 * 
 * To do this, the DataOut::add_vector() function wants a vector of
 * strings. Since the number of components is the same as the number
 * of dimensions we are working in, we use the <code>switch</code>
 * statement below.
 *   

 * 
 * We note that some graphics programs have restriction on what
 * characters are allowed in the names of variables. deal.II therefore
 * supports only the minimal subset of these characters that is supported
 * by all programs. Basically, these are letters, numbers, underscores,
 * and some other characters, but in particular no whitespace and
 * minus/hyphen. The library will throw an exception otherwise, at least
 * if in debug mode.
 *   

 * 
 * After listing the 1d, 2d, and 3d case, it is good style to let the
 * program die if we run upon a case which we did not consider. Remember
 * that the Assert macro generates an exception if the condition in the
 * first parameter is not satisfied. Of course, the condition
 * <code>false</code> can never be satisfied, so the program will always
 * abort whenever it gets to the default statement:
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     DataOut<dim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 * 
 *     std::vector<std::string> solution_names;
 *     switch (dim)
 *       {
 *         case 1:
 *           solution_names.emplace_back("displacement");
 *           break;
 *         case 2:
 *           solution_names.emplace_back("x_displacement");
 *           solution_names.emplace_back("y_displacement");
 *           break;
 *         case 3:
 *           solution_names.emplace_back("x_displacement");
 *           solution_names.emplace_back("y_displacement");
 *           solution_names.emplace_back("z_displacement");
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 * 
 * @endcode
 * 
 * After setting up the names for the different components of the
 * solution vector, we can add the solution vector to the list of
 * data vectors scheduled for output. Note that the following
 * function takes a vector of strings as second argument, whereas
 * the one which we have used in all previous examples accepted a
 * string there. (In fact, the function we had used before would
 * convert the single string into a vector with only one element
 * and forwards that to the other function.)
 * 
 * @code
 *     data_out.add_data_vector(solution, solution_names);
 *     data_out.build_patches();
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtk");
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrun"></a> 
 * <h4>ElasticProblem::run</h4>
 * 

 * 
 * The <code>run</code> function does the same things as in step-6, for
 * example. This time, we use the square [-1,1]^d as domain, and we refine
 * it globally four times before starting the first iteration.
 *   

 * 
 * The reason for refining is a bit accidental: we use the QGauss
 * quadrature formula with two points in each direction for integration of the
 * right hand side; that means that there are four quadrature points on each
 * cell (in 2D). If we only refine the initial grid once globally, then there
 * will be only four quadrature points in each direction on the
 * domain. However, the right hand side function was chosen to be rather
 * localized and in that case, by pure chance, it happens that all quadrature
 * points lie at points where the right hand side function is zero (in
 * mathematical terms, the quadrature points happen to be at points outside
 * the <i>support</i> of the right hand side function). The right hand side
 * vector computed with quadrature will then contain only zeroes (even though
 * it would of course be nonzero if we had computed the right hand side vector
 * exactly using the integral) and the solution of the system of
 * equations is the zero vector, i.e., a finite element function that is zero
 * everywhere. In a sense, we
 * should not be surprised that this is happening since we have chosen
 * an initial grid that is totally unsuitable for the problem at hand.
 *   

 * 
 * The unfortunate thing is that if the discrete solution is constant, then
 * the error indicators computed by the KellyErrorEstimator class are zero
 * for each cell as well, and the call to
 * Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
 * for refinement (why should it if the indicated error is zero for each
 * cell?). The grid in the next iteration will therefore consist of four
 * cells only as well, and the same problem occurs again.
 *   

 * 
 * The conclusion needs to be: while of course we will not choose the
 * initial grid to be well-suited for the accurate solution of the problem,
 * we must at least choose it such that it has the chance to capture the
 * important features of the solution. In this case, it needs to be able to
 * see the right hand side. Thus, we refine globally four times. (Any larger
 * number of global refinement steps would of course also work.)
 * 
 * @code
 *   template <int dim>
 *   void ElasticProblem<dim>::run()
 *   {
 *     for (unsigned int cycle = 0; cycle < 8; ++cycle)
 *       {
 *         std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, -1, 1);
 *             triangulation.refine_global(4);
 *           }
 *         else
 *           refine_grid();
 * 
 *         std::cout << "   Number of active cells:       "
 *                   << triangulation.n_active_cells() << std::endl;
 * 
 *         setup_system();
 * 
 *         std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *                   << std::endl;
 * 
 *         assemble_system();
 *         solve();
 *         output_results(cycle);
 *       }
 *   }
 * } // namespace Step8
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * After closing the <code>Step8</code> namespace in the last line above, the
 * following is the main function of the program and is again exactly like in
 * step-6 (apart from the changed class names, of course).
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step8::ElasticProblem<2> elastic_problem_2d;
 *       elastic_problem_2d.run();
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
examples/step-8/doc/results.dox



<a name="Results"></a><h1>Results</h1>



关于这个程序的结果，除了它们看起来很好之外，没有什么可说的。所有图片都是用VisIt从程序写入磁盘的输出文件中制作的。前两张图片显示了 $x$ -和 $y$ -位移的标量分量。

 <table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.x.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.y.png" alt="">
</td>
</tr>
</table> 


你可以清楚地看到 $x$ 周围的位移 $x=0.5$ 和 $x=-0.5$ 的来源，以及 $y$ 在原点的位移。

人们经常想做的是将位移显示为一个矢量场，也就是说，每一个点的矢量都说明了位移的方向和大小。不幸的是，这就有点麻烦了。为了理解为什么会这样，请记住，我们刚刚将我们的有限元定义为两个分量的集合（在 <code>dim=2</code> 维度）。我们没有说过这不仅仅是一个压力和一个浓度（两个标量），而是说这两个分量实际上是一个矢量值量的一部分，即位移。如果没有这方面的知识，DataOut类就会假定我们打印的所有单个变量都是独立的标量，然后VisIt和Paraview就会忠实地假定这确实是这样的。换句话说，一旦我们把数据写成标量，这些程序中就没有任何东西可以让我们把这两个标量字段粘贴到一起作为一个矢量字段。我们必须从根本上解决这个问题，即在  <code>ElasticProblem::output_results()</code>  。我们不会在这里这样做，而是让读者参考step-22程序，在那里我们展示了如何在一个更普遍的情况下这样做。话虽如此，我们还是忍不住要生成数据，以显示如果按照步骤22中讨论的方式实施，这将是什么样子。矢量场看起来是这样的（VisIt和Paraview随机选择几百个顶点来绘制矢量；从每个顶点绘制矢量会使图片无法阅读）。

 <img src="https://www.dealii.org/images/steps/developer/step-8.vectors.png" alt=""> 


我们注意到，由于 $x$ -和 $y$ -力相对于这些轴是对称的，人们可能直观地期望解是关于 $x$ -和 $y$ -轴的对称。然而，作为矢量的力是不对称的，因此解决方案也不对称。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-8.cc"
*/
