/**
@page step_38 The step-38 tutorial program
This tutorial depends on step-34.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Testcase">Testcase</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeLaplaceBeltramiProblemcodeclasstemplate">The <code>LaplaceBeltramiProblem</code> class template</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeLaplaceBeltramiProblemcodeclass">Implementation of the <code>LaplaceBeltramiProblem</code> class</a>
      <ul>
        <li><a href="#LaplaceBeltramiProblemmake_grid_and_dofs">LaplaceBeltramiProblem::make_grid_and_dofs</a>
        <li><a href="#LaplaceBeltramiProblemassemble_system">LaplaceBeltramiProblem::assemble_system</a>
        <li><a href="#LaplaceBeltramiProblemsolve">LaplaceBeltramiProblem::solve</a>
        <li><a href="#LaplaceBeltramiProblemoutput_result">LaplaceBeltramiProblem::output_result</a>
        <li><a href="#LaplaceBeltramiProblemcompute_error">LaplaceBeltramiProblem::compute_error</a>
        <li><a href="#LaplaceBeltramiProblemrun">LaplaceBeltramiProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-38/doc/intro.dox

 <br> 

<i>This program was contributed by Andrea Bonito and M. Sebastian Pauletti,
with editing and writing by Wolfgang Bangerth.
<br>
This material is based upon work supported by the National Science
Foundation under Grant No. DMS-0914977. Any opinions, findings and conclusions
or recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the National Science Foundation
(NSF).
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


在这个例子中，我们展示了如何解决由四边形组成的一维曲面 $\Gamma \subset \mathbb R^3$ 上的偏微分方程（PDE），即在三维的曲面或二维的直线上。我们重点讨论以下的椭圆二阶PDE

@f{align*}


-\Delta_\Gamma u &= f \qquad \text{on } \qquad \Gamma,\\
u  &= g \qquad \text{on} \qquad \partial \Gamma,


@f}

它概括了我们以前在几个早期教程程序中解决的拉普拉斯方程。我们的实现是基于step-4的。step-34也可以解决低维曲面上的问题；但是，在那里我们只考虑不涉及解变量导数的积分方程，而在这里我们实际上要研究只在一个（可能是弯曲的）曲面上定义的函数的导数是什么意思。

为了定义上述算子，我们首先要介绍一些符号。让 $\mathbf x_S:\hat S \rightarrow S$ 是一个由参考元素 $\hat S \subset \mathbb R^2$ 构成的曲面 $S$ 的参数化，即每个点 $\hat{\mathbf x}\in\hat S$ 诱导出一个点 ${\mathbf
  x}_S(\hat{\mathbf x}) \in S$  。那么让

@f[
G_S\dealcoloneq (D \mathbf{x}_S)^T \ D \mathbf{x}_S


@f]

表示相应的第一基本形式，其中 $D
\mathbf{x}_S=\left(\frac{\partial x_{S,i}(\hat{\mathbf x})}{\partial \hat x_j}\right)_{ij}$ 是映射的导数（雅各布）。在下文中， $S$ 将是整个表面 $\Gamma$ ，或者对有限元方法更方便的是任何面 $S \in
{\mathbb T}$ ，其中 ${\mathbb T}$ 是由四边形构成的 $\Gamma$ 的分区（三角化）。我们现在可以定义一个函数 $v : S \rightarrow \mathbb
R$ 的切向梯度，即

@f[
(\nabla_S v)\circ \mathbf x_S \dealcoloneq  D \mathbf x_S \ G_S^{-1} \ \nabla (v \circ \mathbf x_S).


@f]

表面拉普拉斯(也叫拉普拉斯-贝特拉米算子)的定义是  $\Delta_S \dealcoloneq \nabla_S \cdot \nabla_S$  。请注意，在光滑表面上计算表面梯度的另一种方法  $\Gamma$  是

@f[
\nabla_S v = \nabla \tilde v - \mathbf n (\mathbf n \cdot \nabla \tilde v),


@f]

其中 $\tilde v$ 是 $v$ 在 $\Gamma$ 的管状邻域的 "平滑 "扩展， $\mathbf n$ 是 $\Gamma$ 的法线。由于 $\Delta_S = \nabla_S \cdot \nabla_S$ ，我们推导出

@f[
\Delta_S v = \Delta \tilde v - \mathbf n^T \ D^2 \tilde v \ \mathbf n - (\mathbf n \cdot \nabla \tilde v) (\nabla \cdot \mathbf n - \mathbf n^T \ D \mathbf n \ \mathbf n ).


@f]

值得一提的是，上述表达式中出现的术语 $\nabla \cdot \mathbf n - \mathbf n \ D \mathbf n \ \mathbf n$ 是曲面的总曲率（主曲率之和）。

像往常一样，我们只对弱解感兴趣，为此我们可以使用 $C^0$ 有限元（而不是像强解那样要求 $C^1$ 的连续性）。因此，我们求助于弱的表述

@f[
\int_\Gamma \nabla_\Gamma u \cdot
\nabla_\Gamma v = \int_\Gamma f \ v  \qquad \forall v \in H^1_0(\Gamma)


@f]

并利用分区 ${\mathbb T}$ 的优势，进一步编写

@f[
\sum_{K\in  {\mathbb T}}\int_K \nabla_{K} u \cdot \nabla_{K} v = \sum_{K\in
  {\mathbb T}} \int_K f \ v  \qquad \forall v \in H^1_0(\Gamma).


@f]

此外，上述表达式中的每个积分都是在参考元素 $\hat K \dealcoloneq [0,1]^2$ 中计算的，因此

@f{align*}
\int_{K} \nabla_{K} u \cdot \nabla_{K} v
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} (D \mathbf
  x_K)^T D \mathbf x_K G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}
\\
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}


@f}

和

@f[
\int_{K} f \ v = \int_{\hat K} (f \circ \mathbf x_K) (v \circ \mathbf
x_K)  \sqrt{\det
    (G_K)}.


@f]

最后，我们使用由点 $\{p_l\}_{l=1}^N\subset
\hat K$ 和权重 $\{w_l\}_{l=1}^N \subset \mathbb R^+_*$ 定义的正交公式来评估上述积分，得到

@f[\int_{K} \nabla_{K} u \cdot \nabla_{K} v \approx \sum_{l=1}^N
 (\nabla (u \circ \mathbf x_K)(p_l))^T G^{-1}(p_l)  \nabla (v \circ \mathbf x_K)
(p_l) \sqrt{\det (G(p_l))} \ w_l


@f]

和

@f[
\int_{K} f \ v \approx \sum_{l=1}^N (f \circ \mathbf x_K)(p_l) \ (v \circ \mathbf x_K)(p_l) \sqrt{\det (G(p_l))} \ w_l.


@f]




幸运的是，deal.II已经有了所有的工具来计算上述表达式。事实上，它们与我们求解通常的拉普拉斯的方法几乎没有区别，只需要在FEValues类的构造函数中提供表面坐标映射。这个曲面描述给定，在二维曲面的情况下，两个例程 FEValues::shape_grad 和 FEValues::JxW 会返回

@f{align*}
\text{FEValues::shape\_grad}(i,l)&=D \mathbf x_K(p_l) G^{-1}(p_l)D(\varphi_i \circ \mathbf x_K)
  (p_l)
\\
\text{FEValues::JxW}(l) &=  \sqrt{\det (G(p_l))} \ w_l.


@f}

这正好提供了我们的计算所需的术语。

在更广泛的意义上，表面有限元逼近的细节可以在[Dziuk, in Partial differential equations and calculus of variations 1357, Lecture Notes in Math., 1988], [Demlow, SIAM J. Numer. Anal. 47(2), 2009] 和 [Bonito, Nochetto, and Pauletti, SIAM J. Numer. Anal. 48(5), 2010] 中找到。




<a name="Testcase"></a><h3>Testcase</h3>


一般来说，当你想在数值上测试一个算法的准确性和/或收敛性，你需要提供一个精确的解决方案。通常的技巧是选择一个我们希望成为解决方案的函数，然后对其应用微分算子，为右侧定义一个强制项。这就是我们在这个例子中所做的。在当前情况下，域的形式显然也是至关重要的。

我们为二维问题制作一个测试案例，为三维问题制作另一个测试案例。

 <ul>   <li>  在2d中，让我们选择一个半圆作为域。在这个域上，我们选择函数 $u(\mathbf x)=-2x_1x_2$ 作为解决方案。为了计算右手边，我们必须计算解函数的表面拉普拉斯。有（至少）两种方法可以做到这一点。第一种是使用 $u(\mathbf x)$ 的自然延伸（仍然用 $u$ 表示）在 $\mathbb R^d$ 上投影掉上面描述的法向导数，即计算@f[


    -\Delta_\Gamma u =  \Delta u - \mathbf n^T \ D^2 u \ \mathbf n - (\mathbf n \cdot \nabla u)\ \kappa,
  @f] 。

  其中  $\kappa$  是  $\Gamma$  的总曲率。   由于我们在单位圆上， $\mathbf n=\mathbf x$ 和 $\kappa = 1$ 所以@f[


    -\Delta_\Gamma u = -8 x_1x_2.
  @f]



  一个更简单的方法，至少对于目前二维空间的曲线的情况，是注意到我们可以用变换 $t \in
  [0,\pi]$ 将区间 $\Omega$ 映射到域 $\mathbf x(t)= \left(\begin{array}{c} \cos t \\ \sin t \end{array}\right)$ 。   在位置  $\mathbf x=\mathbf x(t)$  上，解的值是  $u(\mathbf x(t)) = -2\cos t \sin t$  。   考虑到转换是保长的，即长度为 $dt$ 的线段被映射到完全相同长度的曲线上，那么切向拉普拉斯就满足@f{align*}
    \Delta_\Gamma u
    &= \frac{d^2}{dt^2}(-2\cos t \sin t)
    = -2 \frac{d}{dt}(-\sin^2 t + \cos^2 t)
    = -2 (-2 \sin t \cos t - 2 \cos t \sin t)
    \\
    &= 8 \sin t \cos t
    \\
    &= 8 x_1x_2,
  @f} 。

  这当然和我们上面的结果是一样的。   </li>   <li>  在三维中，域又是单位球表面的一半，即半球或圆顶。我们选择 $u(\mathbf x)=-2\sin(\pi x_1)\cos(\pi x_2)e^z$ 作为解决方案。我们可以用上面的方法计算方程的右边， $f=-\Delta_\Gamma u$ ，（用 $\kappa = 2$ ），得到一个笨拙而冗长的表达。你可以在源代码中找到完整的表达式。   </li>   </ul> 。

在程序中，我们还将计算出解的 $H^1$ 半规范误差。由于解函数及其数值近似只在流形上定义，这个误差函数的明显定义是 $| e |_{H^1(\Gamma)}
  = | \nabla_\Gamma e |_{L_2(\Gamma)}
  = \left( \int_\Gamma | \nabla_\Gamma (u-u_h) |^2 \right)^{1/2}$  。这就要求我们为函数 VectorTools::integrate_difference （在步骤7中首次引入）提供<i>tangential</i>梯度 $\nabla_\Gamma u$ ，我们将通过在下面的程序中实现函数 <code>Solution::gradient</code> 来实现。




<a name="Implementation"></a><h3>Implementation</h3>


如果你已经读完了第4步，并且理解了上面关于解和右边如何对应的讨论，你也会立即熟悉这个程序。事实上，只有两件事是有意义的。

- 我们生成三角计算域的网格的方式。

- 我们使用映射对象的方式来描述，我们解决偏微分方程的领域不是平面的，实际上是弯曲的。

在第10步和第11步中已经介绍了映射对象，正如那里所解释的，只要你对边界的样子有一个有效的描述，你通常不需要知道它们是如何工作的。从本质上讲，我们将简单地声明一个适当的MappingQ类型的对象，它将自动从三角图中获得边界描述。然后，该映射对象将被传递给适当的函数，我们将得到库中预定义的半圆或半球形的边界描述。

该程序的其余部分紧跟步骤4，至于计算误差，则是步骤7。这个程序的某些方面，特别是在Triangulation、DoFHandler和类似的类上使用两个模板参数，已经在步骤34中作了详细描述；你不妨也读一读这个教程程序。


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
 * If you've read through step-4 and step-7, you will recognize that we have
 * used all of the following include files there already. Consequently, we
 * will not explain their meaning here again.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/solver_control.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * namespace Step38
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceBeltramiProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceBeltramiProblem</code> class template</h3>
 * 

 * 
 * This class is almost exactly similar to the <code>LaplaceProblem</code>
 * class in step-4.
 * 

 * 
 * The essential differences are these:
 *   

 * 
 * - The template parameter now denotes the dimensionality of the embedding
 * space, which is no longer the same as the dimensionality of the domain
 * and the triangulation on which we compute. We indicate this by calling
 * the parameter @p spacedim, and introducing a constant @p dim equal to
 * the dimensionality of the domain -- here equal to
 * <code>spacedim-1</code>.
 * - All member variables that have geometric aspects now need to know about
 * both their own dimensionality as well as that of the embedding
 * space. Consequently, we need to specify both of their template
 * parameters one for the dimension of the mesh @p dim, and the other for
 * the dimension of the embedding space, @p spacedim. This is exactly what
 * we did in step-34, take a look there for a deeper explanation.
 * - We need an object that describes which kind of mapping to use from the
 * reference cell to the cells that the triangulation is composed of. The
 * classes derived from the Mapping base class do exactly this. Throughout
 * most of deal.II, if you don't do anything at all, the library assumes
 * that you want an object of kind MappingQ1 that uses a (bi-, tri-)linear
 * mapping. In many cases, this is quite sufficient, which is why the use
 * of these objects is mostly optional: for example, if you have a
 * polygonal two-dimensional domain in two-dimensional space, a bilinear
 * mapping of the reference cell to the cells of the triangulation yields
 * an exact representation of the domain. If you have a curved domain, one
 * may want to use a higher order mapping for those cells that lie at the
 * boundary of the domain -- this is what we did in step-11, for
 * example. However, here we have a curved domain, not just a curved
 * boundary, and while we can approximate it with bilinearly mapped cells,
 * it is really only prudent to use a higher order mapping for all
 * cells. Consequently, this class has a member variable of type MappingQ;
 * we will choose the polynomial degree of the mapping equal to the
 * polynomial degree of the finite element used in the computations to
 * ensure optimal approximation, though this iso-parametricity is not
 * required.
 * 
 * @code
 *   template <int spacedim>
 *   class LaplaceBeltramiProblem
 *   {
 *   public:
 *     LaplaceBeltramiProblem(const unsigned degree = 2);
 *     void run();
 * 
 *   private:
 *     static constexpr unsigned int dim = spacedim - 1;
 * 
 *     void make_grid_and_dofs();
 *     void assemble_system();
 *     void solve();
 *     void output_results() const;
 *     void compute_error() const;
 * 
 * 
 *     Triangulation<dim, spacedim> triangulation;
 *     FE_Q<dim, spacedim>          fe;
 *     DoFHandler<dim, spacedim>    dof_handler;
 *     MappingQ<dim, spacedim>      mapping;
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
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * Next, let us define the classes that describe the exact solution and the
 * right hand sides of the problem. This is in analogy to step-4 and step-7
 * where we also defined such objects. Given the discussion in the
 * introduction, the actual formulas should be self-explanatory. A point of
 * interest may be how we define the value and gradient functions for the 2d
 * and 3d cases separately, using explicit specializations of the general
 * template. An alternative to doing it this way might have been to define
 * the general template and have a <code>switch</code> statement (or a
 * sequence of <code>if</code>s) for each possible value of the spatial
 * dimension.
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual Tensor<1, dim>
 *     gradient(const Point<dim> & p,
 *              const unsigned int component = 0) const override;
 *   };
 * 
 * 
 *   template <>
 *   double Solution<2>::value(const Point<2> &p, const unsigned int) const
 *   {
 *     return (-2. * p(0) * p(1));
 *   }
 * 
 * 
 *   template <>
 *   Tensor<1, 2> Solution<2>::gradient(const Point<2> &p,
 *                                      const unsigned int) const
 *   {
 *     Tensor<1, 2> return_value;
 *     return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0));
 *     return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1));
 * 
 *     return return_value;
 *   }
 * 
 * 
 *   template <>
 *   double Solution<3>::value(const Point<3> &p, const unsigned int) const
 *   {
 *     return (std::sin(numbers::PI * p(0)) * std::cos(numbers::PI * p(1)) *
 *             exp(p(2)));
 *   }
 * 
 * 
 *   template <>
 *   Tensor<1, 3> Solution<3>::gradient(const Point<3> &p,
 *                                      const unsigned int) const
 *   {
 *     using numbers::PI;
 * 
 *     Tensor<1, 3> return_value;
 * 
 *     return_value[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     return_value[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     return_value[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     return return_value;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 *   template <>
 *   double RightHandSide<2>::value(const Point<2> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     return (-8. * p(0) * p(1));
 *   }
 * 
 * 
 *   template <>
 *   double RightHandSide<3>::value(const Point<3> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     using numbers::PI;
 * 
 *     Tensor<2, 3> hessian;
 * 
 *     hessian[0][0] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[1][1] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[2][2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     hessian[0][1] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     hessian[1][0] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 * 
 *     hessian[0][2] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     hessian[2][0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     hessian[1][2] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     hessian[2][1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 * 
 *     Tensor<1, 3> gradient;
 *     gradient[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 *     gradient[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 *     gradient[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 * 
 *     Point<3> normal = p;
 *     normal /= p.norm();
 * 
 *     return (-trace(hessian) + 2 * (gradient * normal) +
 *             (hessian * normal) * normal);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeLaplaceBeltramiProblemcodeclass"></a> 
 * <h3>Implementation of the <code>LaplaceBeltramiProblem</code> class</h3>
 * 

 * 
 * The rest of the program is actually quite unspectacular if you know
 * step-4. Our first step is to define the constructor, setting the
 * polynomial degree of the finite element and mapping, and associating the
 * DoF handler to the triangulation:
 * 
 * @code
 *   template <int spacedim>
 *   LaplaceBeltramiProblem<spacedim>::LaplaceBeltramiProblem(
 *     const unsigned degree)
 *     : fe(degree)
 *     , dof_handler(triangulation)
 *     , mapping(degree)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemmake_grid_and_dofs"></a> 
 * <h4>LaplaceBeltramiProblem::make_grid_and_dofs</h4>
 * 

 * 
 * The next step is to create the mesh, distribute degrees of freedom, and
 * set up the various variables that describe the linear system. All of
 * these steps are standard with the exception of how to create a mesh that
 * describes a surface. We could generate a mesh for the domain we are
 * interested in, generate a triangulation using a mesh generator, and read
 * it in using the GridIn class. Or, as we do here, we generate the mesh
 * using the facilities in the GridGenerator namespace.
 *   

 * 
 * In particular, what we're going to do is this (enclosed between the set
 * of braces below): we generate a <code>spacedim</code> dimensional mesh
 * for the half disk (in 2d) or half ball (in 3d), using the
 * GridGenerator::half_hyper_ball function. This function sets the boundary
 * indicators of all faces on the outside of the boundary to zero for the
 * ones located on the perimeter of the disk/ball, and one on the straight
 * part that splits the full disk/ball into two halves. The next step is the
 * main point: The GridGenerator::extract_boundary_mesh function creates a
 * mesh that consists of those cells that are the faces of the previous mesh,
 * i.e. it describes the <i>surface</i> cells of the original (volume)
 * mesh. However, we do not want all faces: only those on the perimeter of
 * the disk or ball which carry boundary indicator zero; we can select these
 * cells using a set of boundary indicators that we pass to
 * GridGenerator::extract_boundary_mesh.
 *   

 * 
 * There is one point that needs to be mentioned. In order to refine a
 * surface mesh appropriately if the manifold is curved (similarly to
 * refining the faces of cells that are adjacent to a curved boundary), the
 * triangulation has to have an object attached to it that describes where
 * new vertices should be located. If you don't attach such a boundary
 * object, they will be located halfway between existing vertices; this is
 * appropriate if you have a domain with straight boundaries (e.g. a
 * polygon) but not when, as here, the manifold has curvature. So for things
 * to work properly, we need to attach a manifold object to our (surface)
 * triangulation, in much the same way as we've already done in 1d for the
 * boundary. We create such an object and attach it to the triangulation.
 *   

 * 
 * The final step in creating the mesh is to refine it a number of
 * times. The rest of the function is the same as in previous tutorial
 * programs.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs()
 *   {
 *     {
 *       Triangulation<spacedim> volume_mesh;
 *       GridGenerator::half_hyper_ball(volume_mesh);
 * 
 *       std::set<types::boundary_id> boundary_ids;
 *       boundary_ids.insert(0);
 * 
 *       GridGenerator::extract_boundary_mesh(volume_mesh,
 *                                            triangulation,
 *                                            boundary_ids);
 *     }
 *     triangulation.set_all_manifold_ids(0);
 *     triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());
 * 
 *     triangulation.refine_global(4);
 * 
 *     std::cout << "Surface mesh has " << triangulation.n_active_cells()
 *               << " cells." << std::endl;
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     std::cout << "Surface mesh has " << dof_handler.n_dofs()
 *               << " degrees of freedom." << std::endl;
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *     sparsity_pattern.copy_from(dsp);
 * 
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
 * <a name="LaplaceBeltramiProblemassemble_system"></a> 
 * <h4>LaplaceBeltramiProblem::assemble_system</h4>
 * 

 * 
 * The following is the central function of this program, assembling the
 * matrix that corresponds to the surface Laplacian (Laplace-Beltrami
 * operator). Maybe surprisingly, it actually looks exactly the same as for
 * the regular Laplace operator discussed in, for example, step-4. The key
 * is that the FEValues::shape_grad() function does the magic: It returns
 * the surface gradient $\nabla_K \phi_i(x_q)$ of the $i$th shape function
 * at the $q$th quadrature point. The rest then does not need any changes
 * either:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const QGauss<dim>       quadrature_formula(2 * fe.degree);
 *     FEValues<dim, spacedim> fe_values(mapping,
 *                                       fe,
 *                                       quadrature_formula,
 *                                       update_values | update_gradients |
 *                                         update_quadrature_points |
 *                                         update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<double>                  rhs_values(n_q_points);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     RightHandSide<spacedim> rhs;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 *         rhs.value_list(fe_values.get_quadrature_points(), rhs_values);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *               cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
 *                                    fe_values.shape_grad(j, q_point) *
 *                                    fe_values.JxW(q_point);
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             cell_rhs(i) += fe_values.shape_value(i, q_point) *
 *                            rhs_values[q_point] * fe_values.JxW(q_point);
 * 
 *         cell->get_dof_indices(local_dof_indices);
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *               system_matrix.add(local_dof_indices[i],
 *                                 local_dof_indices[j],
 *                                 cell_matrix(i, j));
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *           }
 *       }
 * 
 *     std::map<types::global_dof_index, double> boundary_values;
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 0, Solution<spacedim>(), boundary_values);
 * 
 *     MatrixTools::apply_boundary_values(
 *       boundary_values, system_matrix, solution, system_rhs, false);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemsolve"></a> 
 * <h4>LaplaceBeltramiProblem::solve</h4>
 * 

 * 
 * The next function is the one that solves the linear system. Here, too, no
 * changes are necessary:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::solve()
 *   {
 *     SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm());
 *     SolverCG<Vector<double>> cg(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemoutput_result"></a> 
 * <h4>LaplaceBeltramiProblem::output_result</h4>
 * 

 * 
 * This is the function that generates graphical output from the
 * solution. Most of it is boilerplate code, but there are two points worth
 * pointing out:
 *   

 * 
 * - The DataOut::add_data_vector() function can take two kinds of vectors:
 * Either vectors that have one value per degree of freedom defined by the
 * DoFHandler object previously attached via DataOut::attach_dof_handler();
 * and vectors that have one value for each cell of the triangulation, for
 * example to output estimated errors for each cell. Typically, the
 * DataOut class knows to tell these two kinds of vectors apart: there are
 * almost always more degrees of freedom than cells, so we can
 * differentiate by the two kinds looking at the length of a vector. We
 * could do the same here, but only because we got lucky: we use a half
 * sphere. If we had used the whole sphere as domain and $Q_1$ elements,
 * we would have the same number of cells as vertices and consequently the
 * two kinds of vectors would have the same number of elements. To avoid
 * the resulting confusion, we have to tell the DataOut::add_data_vector()
 * function which kind of vector we have: DoF data. This is what the third
 * argument to the function does.
 * - The DataOut::build_patches() function can generate output that subdivides
 * each cell so that visualization programs can resolve curved manifolds
 * or higher polynomial degree shape functions better. We here subdivide
 * each element in each coordinate direction as many times as the
 * polynomial degree of the finite element in use.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::output_results() const
 *   {
 *     DataOut<dim, spacedim> data_out;
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution,
 *                              "solution",
 *                              DataOut<dim, spacedim>::type_dof_data);
 *     data_out.build_patches(mapping, mapping.get_degree());
 * 
 *     const std::string filename =
 *       "solution-" + std::to_string(spacedim) + "d.vtk";
 *     std::ofstream output(filename);
 *     data_out.write_vtk(output);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemcompute_error"></a> 
 * <h4>LaplaceBeltramiProblem::compute_error</h4>
 * 

 * 
 * This is the last piece of functionality: we want to compute the error in
 * the numerical solution. It is a verbatim copy of the code previously
 * shown and discussed in step-7. As mentioned in the introduction, the
 * <code>Solution</code> class provides the (tangential) gradient of the
 * solution. To avoid evaluating the error only a superconvergence points,
 * we choose a quadrature rule of sufficiently high order.
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::compute_error() const
 *   {
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Solution<spacedim>(),
 *                                       difference_per_cell,
 *                                       QGauss<dim>(2 * fe.degree + 1),
 *                                       VectorTools::H1_norm);
 * 
 *     double h1_error = VectorTools::compute_global_error(triangulation,
 *                                                         difference_per_cell,
 *                                                         VectorTools::H1_norm);
 *     std::cout << "H1 error = " << h1_error << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemrun"></a> 
 * <h4>LaplaceBeltramiProblem::run</h4>
 * 

 * 
 * The last function provides the top-level logic. Its contents are
 * self-explanatory:
 * 
 * @code
 *   template <int spacedim>
 *   void LaplaceBeltramiProblem<spacedim>::run()
 *   {
 *     make_grid_and_dofs();
 *     assemble_system();
 *     solve();
 *     output_results();
 *     compute_error();
 *   }
 * } // namespace Step38
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * The remainder of the program is taken up by the <code>main()</code>
 * function. It follows exactly the general layout first introduced in step-6
 * and used in all following tutorial programs:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step38;
 * 
 *       LaplaceBeltramiProblem<3> laplace_beltrami;
 *       laplace_beltrami.run();
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
 * 
 *   return 0;
 * }
 * @endcode
examples/step-38/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你运行该程序时，应在屏幕上打印出以下输出。

@verbatim
Surface mesh has 1280 cells.
Surface mesh has 5185 degrees of freedom.
H1 error = 0.0217136
@endverbatim




通过在 <code>LaplaceBeltrami::make_grid_and_dofs</code> 函数中玩弄全局细化的数量，可以增加或减少网格细化。例如，多做一次细化，只运行三维曲面问题，得到的输出结果如下。

@verbatim
Surface mesh has 5120 cells.
Surface mesh has 20609 degrees of freedom.
H1 error = 0.00543481
@endverbatim



这就是我们所期望的：将网格尺寸缩小2倍，误差下降4倍（记住我们使用的是双二次元）。从一到五次细化的全部误差序列看起来是这样的，整齐地遵循理论上预测的模式。

@verbatim
0.339438
0.0864385
0.0217136
0.00543481
0.00135913
@endverbatim



最后，该程序产生图形输出，我们可以将其可视化。下面是一个结果图。

 <img src="https://www.dealii.org/images/steps/developer/step-38.solution-3d.png" alt=""> 

该程序也适用于2D中的1D曲线，而不仅仅是3D中的2D曲面。你可以通过改变 <code>main()</code> 中的模板参数来测试这一点，像这样。

@code
      LaplaceBeltramiProblem<2> laplace_beltrami;
@endcode

域是一条2D的曲线，我们可以通过使用第三维（和颜色）来表示函数 $u(x)$ 的值来可视化解决方案。这样看起来就像这样（白色的曲线是域，彩色的曲线是被挤压到第三维的解决方案，清楚地显示了当曲线从域的一个象限移动到相邻的象限时符号的变化）。

 <img src="https://www.dealii.org/images/steps/developer/step-38.solution-2d.png" alt=""> 


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


只有当表面比仅仅是一个半球体更有趣时，在表面上的计算才会变得有趣。为了达到这个目的，deal.II可以通过通常的GridIn类读取描述曲面的网格。或者，万一你有一个分析性的描述，一个简单的网格有时可以被拉伸和弯曲成我们感兴趣的形状。

让我们考虑一个相对简单的例子：我们把之前用过的半球体，在Z方向上拉伸10倍，然后把X和Y坐标拼一下。让我们先展示一下计算域和解决方案，然后再讨论下面的实现细节。

 <img src="https://www.dealii.org/images/steps/developer/step-38.warp-1.png" alt=""> 

 <img src="https://www.dealii.org/images/steps/developer/step-38.warp-2.png" alt=""> 

产生这种网格的方法是使用 GridTools::transform() 函数。它需要一个方法来转换每个单独的网格点到不同的位置。让我们在这里使用下面这个相当简单的函数（记住：在一个方向上拉伸，在另外两个方向上拼凑）。

@code
template <int spacedim>
Point<spacedim> warp(const Point<spacedim> &p)
{
  Point<spacedim> q = p;
  q[spacedim-1] *= 10;


  if (spacedim >= 2)
    q[0] += 2*std::sin(q[spacedim-1]);
  if (spacedim >= 3)
    q[1] += 2*std::cos(q[spacedim-1]);


  return q;
}
@endcode



如果我们遵循 <code>LaplaceBeltrami::make_grid_and_dofs</code> 函数，我们会像以前一样提取半球形表面网格，将其扭曲成我们想要的形状，并根据需要经常进行细化。但这并不像我们所希望的那样简单：细化需要我们有一个适当的流形对象附加到三角形上，描述细化时网格的新顶点应该位于何处。我相信可以通过简单地撤销上面的变换（重新得到球面），找到球面上新的点的位置，然后重新扭曲结果，以一种不太复杂的方式描述这个流形。但我是个懒人，既然这样做并不是真正的重点，我们还是让我们的生活变得简单一点：我们将提取半球体，根据需要对其进行细化，摆脱描述流形的对象，因为我们现在不再需要它，然后最后对网格进行扭曲。使用上面的函数，这将看起来如下。

@code
template <int spacedim>
void LaplaceBeltrami<spacedim>::make_grid_and_dofs()
{
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::half_hyper_ball(volume_mesh);


    volume_mesh.refine_global(4);


    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);


    GridGenerator::extract_boundary_mesh(volume_mesh, triangulation,
                                         boundary_ids);
    GridTools::transform(&warp<spacedim>, triangulation);       /* ** */
    std::ofstream x("x"), y("y");
    GridOut().write_gnuplot(volume_mesh, x);
    GridOut().write_gnuplot(triangulation, y);
  }


  std::cout << "Surface mesh has " << triangulation.n_active_cells()
            << " cells."
            << std::endl;
  ...
}
@endcode



请注意，唯一必要的补充是标有星号的那一行。不过值得指出的是：由于我们将流形描述从表面网格中分离出来，所以当我们在程序的其余部分使用映射对象时，它不再有曲线边界描述可言。相反，它将不得不使用隐含的FlatManifold类，该类用于域的所有未明确指定不同流形对象的部分。因此，无论我们使用MappingQ(2)、MappingQ(15)还是MappingQ1，我们的网格的每个单元都将使用双线性近似进行映射。

撇开所有这些缺点不谈，得到的图片还是很好看的。与步骤38中的内容唯一不同的是，我们把右手边改为 $f(\mathbf x)=\sin x_3$ ，把边界值（通过 <code>Solution</code> 类）改为 $u(\mathbf x)|_{\partial\Omega}=\cos x_3$  。当然，我们现在已经不知道确切的解决方案，所以在 <code>LaplaceBeltrami::run</code> 末尾的误差计算将得到一个毫无意义的数字。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-38.cc"
*/
