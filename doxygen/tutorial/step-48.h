/**
@page step_48 The step-48 tutorial program
This tutorial depends on step-25, step-37.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemstatementanddiscretization"> Problem statement and discretization </a>
        <li><a href="#Implementationofconstraints">Implementation of constraints</a>
        <li><a href="#Parallelization"> Parallelization </a>
        <li><a href="#Thetestcase"> The test case </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#SineGordonOperation">SineGordonOperation</a>
      <ul>
        <li><a href="#SineGordonOperationSineGordonOperation">SineGordonOperation::SineGordonOperation</a>
        <li><a href="#SineGordonOperationlocal_apply">SineGordonOperation::local_apply</a>
        <li><a href="#SineGordonOperationapply">SineGordonOperation::apply</a>
      </ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#SineGordonProblemclass">SineGordonProblem class</a>
      <ul>
        <li><a href="#SineGordonProblemSineGordonProblem">SineGordonProblem::SineGordonProblem</a>
        <li><a href="#SineGordonProblemmake_grid_and_dofs">SineGordonProblem::make_grid_and_dofs</a>
        <li><a href="#SineGordonProblemoutput_results">SineGordonProblem::output_results</a>
        <li><a href="#SineGordonProblemrun">SineGordonProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Comparisonwithasparsematrix">Comparison with a sparse matrix</a>
        <li><a href="#Parallelrunin2Dand3D">Parallel run in 2D and 3D</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-48/doc/intro.dox



<i>
This program was contributed by Katharina Kormann and Martin
Kronbichler.


The algorithm for the matrix-vector product is based on the article <a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic
interface for parallel cell-based finite element operator
application</a><a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic
interface for parallel cell-based finite element operator
application</a> by Martin Kronbichler and Katharina Kormann, Computers
and Fluids 63:135&ndash;147, 2012, and the paper &quot;Parallel finite element operator
application: Graph partitioning and coloring&quot; by Katharina
Kormann and Martin Kronbichler in: Proceedings of the 7th IEEE
International Conference on e-Science, 2011.  </i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序演示了如何使用基于单元的有限元算子与MatrixFree类的实现，这是在step-37中首次介绍的，用于解决非线性偏微分方程。此外，我们再看一下无矩阵框架内对约束条件的处理。最后，我们将使用显式时间步进方法来解决问题，并介绍高斯-洛巴托有限元，在这种情况下非常方便，因为它们的质量矩阵可以准确地被对角线矩阵所接近，因此是可逆的。这一特性的两个成分是：首先，根据Gauss-Lobatto正交规则的点分布，对Lagrange多项式的结点进行分布。其次，正交是用同样的Gauss-Lobatto正交规则完成的。在这个公式中，只要 $\int_K \varphi_i \varphi_j
dx\approx \sum_q \varphi_i \varphi_j \mathrm{det}(J) \big |_{x_q}$ ，积分 $i\neq j$ 就会变成零，因为在定义拉格朗日多项式的点中，正好有一个函数 $\varphi_j$ 是一，其他都是零。此外，拉格朗日多项式的节点的Gauss-Lobatto分布将节点向元素边界聚集。这就为高阶离散化方法提供了条件良好的多项式基础。事实上，具有等距节点的FE_Q元素的条件数随着度数的增加而呈指数级增长，这破坏了约5级以上的任何好处。由于这个原因，高斯-洛巴托点是FE_Q元素的默认分布（但在度数为1和2时，这些点相当于等距点）。

<a name="Problemstatementanddiscretization"></a><h3> Problem statement and discretization </h3>


作为一个例子，我们选择解决正弦-戈登孤子方程

\f{eqnarray*}
u_{tt} &=& \Delta u -\sin(u) \quad\mbox{for}\quad (x,t) \in
\Omega \times (t_0,t_f],\\
{\mathbf n} \cdot \nabla u &=& 0
\quad\mbox{for}\quad (x,t) \in \partial\Omega \times (t_0,t_f],\\
u(x,t_0) &=& u_0(x).
\f}

在步骤25中已经介绍过。作为一种简单的显式时间积分方法，我们选择使用方程的二阶表述的跃迁蛙方案。通过这个时间步长，该方案以弱的形式读取

\f{eqnarray*}
(v,u^{n+1}) = (v,2 u^n-u^{n-1} -
(\Delta t)^2 \sin(u^n)) - (\nabla v, (\Delta t)^2 \nabla u^n),
\f}其中<i> v</i>表示一个测试函数，索引<i>n</i>代表时间步数。

对于空间离散化，我们选择FE_Q元素，其基函数定义为插值高斯-洛巴托正交规则的支持点。此外，当我们计算基函数的积分以形成质量矩阵和上述方程右边的算子时，我们使用高斯-洛巴托正交规则，其支持点与有限元的节点点相同，以评估积分。由于有限元是拉格朗日的，这将产生方程左侧的对角线质量矩阵，使每个时间步长的线性系统的解变得微不足道。

使用这个正交规则，对于<i>p</i>th阶有限元，我们使用<i>(2p-1)</i>th阶精确公式来评估积分。由于在计算质量矩阵时，两个<i>p</i>阶基函数的乘积在每个方向上给出了一个具有多项式程度<i>2p</i>的函数，所以积分的计算并不精确。  然而，在具有仿生元素形状的网格上，整体收敛特性不受正交误差的干扰，L2误差与<i>h<sup>p+1</sup></i>成正比。但是请注意，当积分不再是多项式时，一些三维设置的L2误差<i>O(h<sup>p</sup>)</i>甚至<i>O(h<sup>p-1</sup>)</i>的次优收敛率的阶次减少已被报道<a href="https://dx.doi.org/10.1002/num.20353">in
literature</a>在变形（非affine）元素形状的波方程上。

除了在使用显式时间步进时我们可以避免用这种类型的元素来解决线性系统外，它们还具有另外两个优点。当我们使用和-因子化方法来评估有限元算子时（参见步骤37），我们必须在正交点评估函数。在Gauss-Lobatto元素的情况下，正交点和有限元的节点点重合，这种操作是微不足道的，因为正交点的函数值是由其一维系数给出的。这样一来，与一般的高斯正交相比，有限元算子评估的算术工作减少了大约两倍。

总结一下讨论，通过使用正确的有限元和正交规则组合，我们最终得到一个方案，我们只需要计算对应于上述公式的右手边向量，然后在每个时间步骤中乘以对角线质量矩阵的逆。当然，在实践中，我们提取对角线元素，只在程序开始时反转一次。

<a name="Implementationofconstraints"></a><h3>Implementation of constraints</h3>


在 <code>deal.II</code> 中处理约束的通常方法是使用AffineConstraints类，该类建立了一个稀疏矩阵，存储关于哪些自由度（DoF）被约束以及如何被约束的信息。这种格式使用了不必要的大量内存，因为没有那么多不同类型的约束：例如，在每个单元上使用线性有限元时，悬挂节点的情况下，大多数约束具有 $x_k = \frac 12 x_i + \frac 12 x_j$ 的形式，其中系数 $\frac 12$ 总是相同，只有 $i,j,k$ 不同。虽然存储这些多余的信息在一般情况下不是问题，因为在矩阵和右手边的装配过程中只需要一次，但在无矩阵的方法中，它成为一个瓶颈，因为在那里，每次我们应用算子时都要访问这些信息，而算子评估的其余部分是如此之快。因此，MatrixFree使用一个我们称为 <code>constraint_pool</code> 的变量来收集不同约束的权重，而不是AffineConstraints对象。然后，只需要存储网格中每个约束的标识符而不是所有的权重。此外，约束不是在前后处理步骤中应用的，而是在我们评估有限元算子时应用的。因此，约束信息被嵌入到变量 <code>indices_local_to_global</code> 中，用于从全局矢量中提取单元信息。如果一个DoF被约束， <code>indices_local_to_global</code> 变量包含它被约束的DoF的全局索引。然后，我们手头还有另一个变量 <code>constraint_indicator</code> ，对于每个单元，持有被约束的DoF的局部指数以及约束类型的标识符。幸运的是，你不会在示例程序中看到这些数据结构，因为类 <code>FEEvaluation</code> 会在没有用户互动的情况下处理这些约束。

在存在悬空节点的情况下，通过Gauss-Lobatto正交/节点点程序在元素层面获得的对角线质量矩阵并不能直接转化为对角线全局质量矩阵，因为遵循行和列的约束也会增加非对角线的条目。正如在<a href="https://dx.doi.org/10.4208/cicp.101214.021015a">Kormann
(2016)</a>中所解释的，在一个矢量上插值约束，保持质量矩阵的对角线形状，与方程一致，直到与正交误差相同大小的误差。在下面的程序中，我们将简单地组装质量矩阵的对角线，就像它是一个矢量一样，以实现这种近似。




<a name="Parallelization"></a><h3> Parallelization </h3>


MatrixFree类可以在三个层次上进行并行化。分布式节点集群上的MPI并行化，由线程积木库安排的线程并行化，以及最后通过SIMD数据类型在两个（或更多）单元的批次上工作的矢量化（有时称为跨元素或外部矢量化）。正如我们在第37步中已经讨论过的，通过使用特定于你的系统的指令集，你将得到最好的性能，例如，使用cmake变量<tt>-DCMAKE_CXX_FLAGS="-march=native"</tt>。MPI并行化已经在步骤37中被利用了。这里，我们额外考虑用TBB进行线程并行化。这相当简单，因为我们需要做的就是告诉MatrixFree对象的初始化，我们想通过变量 MatrixFree::AdditionalData::thread_parallel_scheme. 来使用线程并行方案。 在设置过程中，建立了一个类似于 @ref workstream_paper 中描述的依赖图，这允许安排 @p local_apply 函数在单元块上的工作，而没有几个线程访问同一个向量索引。相对于WorkStream循环，还应用了一些额外的巧妙技巧来避免<a
href="https://dx.doi.org/10.1109/eScience.2011.53">Kormann and Kronbichler
(2011)</a>中描述的全局同步。

请注意，这个程序是为分布式三角测量 (parallel::distributed::Triangulation), 而设计的，它要求deal.II配置<a href="http://www.p4est.org/">p4est</a>，如<a href="../../readme.html">deal.II ReadMe</a>文件中所述。然而，也支持非分布式三角法，在这种情况下，计算将以串行方式运行。

<a name="Thetestcase"></a><h3> The test case </h3>


在我们的例子中，我们选择初始值为\f{eqnarray*} u(x,t) =
\prod_{i=1}^{d} -4 \arctan \left(
\frac{m}{\sqrt{1-m^2}}\frac{\sin\left(\sqrt{1-m^2} t +c_2\right)}{\cosh(mx_i+c_1)}\right)
\f}，并在时间区间[-10,10]内解决方程。常数被选择为 $c_1=c_1=0$ 和<i> m=0.5</i>。如步骤25所述，在一维中<i>u</i>作为<i>t</i>的函数是正弦-戈登方程的精确解。然而，对于更高的维度，情况并非如此。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The necessary files from the deal.II library.
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/distributed/tria.h>
 * 
 * @endcode
 * 
 * This includes the data structures for the efficient implementation of
 * matrix-free methods.
 * 
 * @code
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * #include <iomanip>
 * 
 * 
 * namespace Step48
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * We start by defining two global variables to collect all parameters
 * subject to changes at one place: One for the dimension and one for the
 * finite element degree. The dimension is used in the main function as a
 * template argument for the actual classes (like in all other deal.II
 * programs), whereas the degree of the finite element is more crucial, as
 * it is passed as a template argument to the implementation of the
 * Sine-Gordon operator. Therefore, it needs to be a compile-time constant.
 * 
 * @code
 *   const unsigned int dimension = 2;
 *   const unsigned int fe_degree = 4;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperation"></a> 
 * <h3>SineGordonOperation</h3>
 * 

 * 
 * The <code>SineGordonOperation</code> class implements the cell-based
 * operation that is needed in each time step. This nonlinear operation can
 * be implemented straight-forwardly based on the <code>MatrixFree</code>
 * class, in the same way as a linear operation would be treated by this
 * implementation of the finite element operator application. We apply two
 * template arguments to the class, one for the dimension and one for the
 * degree of the finite element. This is a difference to other functions in
 * deal.II where only the dimension is a template argument. This is
 * necessary to provide the inner loops in @p FEEvaluation with information
 * about loop lengths etc., which is essential for efficiency. On the other
 * hand, it makes it more challenging to implement the degree as a run-time
 * parameter.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   class SineGordonOperation
 *   {
 *   public:
 *     SineGordonOperation(const MatrixFree<dim, double> &data_in,
 *                         const double                   time_step);
 * 
 *     void apply(LinearAlgebra::distributed::Vector<double> &dst,
 *                const std::vector<LinearAlgebra::distributed::Vector<double> *>
 *                  &src) const;
 * 
 *   private:
 *     const MatrixFree<dim, double> &            data;
 *     const VectorizedArray<double>              delta_t_sqr;
 *     LinearAlgebra::distributed::Vector<double> inv_mass_matrix;
 * 
 *     void local_apply(
 *       const MatrixFree<dim, double> &                                  data,
 *       LinearAlgebra::distributed::Vector<double> &                     dst,
 *       const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
 *       const std::pair<unsigned int, unsigned int> &cell_range) const;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationSineGordonOperation"></a> 
 * <h4>SineGordonOperation::SineGordonOperation</h4>
 * 

 * 
 * This is the constructor of the SineGordonOperation class. It receives a
 * reference to the MatrixFree holding the problem information and the time
 * step size as input parameters. The initialization routine sets up the
 * mass matrix. Since we use Gauss-Lobatto elements, the mass matrix is a
 * diagonal matrix and can be stored as a vector. The computation of the
 * mass matrix diagonal is simple to achieve with the data structures
 * provided by FEEvaluation: Just loop over all cell batches, i.e.,
 * collections of cells due to SIMD vectorization, and integrate over the
 * function that is constant one on all quadrature points by using the
 * <code>integrate</code> function with @p true argument at the slot for
 * values. Finally, we invert the diagonal entries to have the inverse mass
 * matrix directly available in each time step.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   SineGordonOperation<dim, fe_degree>::SineGordonOperation(
 *     const MatrixFree<dim, double> &data_in,
 *     const double                   time_step)
 *     : data(data_in)
 *     , delta_t_sqr(make_vectorized_array(time_step * time_step))
 *   {
 *     data.initialize_dof_vector(inv_mass_matrix);
 * 
 *     FEEvaluation<dim, fe_degree> fe_eval(data);
 *     const unsigned int           n_q_points = fe_eval.n_q_points;
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
 *       {
 *         fe_eval.reinit(cell);
 *         for (unsigned int q = 0; q < n_q_points; ++q)
 *           fe_eval.submit_value(make_vectorized_array(1.), q);
 *         fe_eval.integrate(EvaluationFlags::values);
 *         fe_eval.distribute_local_to_global(inv_mass_matrix);
 *       }
 * 
 *     inv_mass_matrix.compress(VectorOperation::add);
 *     for (unsigned int k = 0; k < inv_mass_matrix.locally_owned_size(); ++k)
 *       if (inv_mass_matrix.local_element(k) > 1e-15)
 *         inv_mass_matrix.local_element(k) =
 *           1. / inv_mass_matrix.local_element(k);
 *       else
 *         inv_mass_matrix.local_element(k) = 1;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationlocal_apply"></a> 
 * <h4>SineGordonOperation::local_apply</h4>
 * 

 * 
 * This operator implements the core operation of the program, the
 * integration over a range of cells for the nonlinear operator of the
 * Sine-Gordon problem. The implementation is based on the FEEvaluation
 * class as in step-37. Due to the special structure in Gauss-Lobatto
 * elements, certain operations become simpler, in particular the evaluation
 * of shape function values on quadrature points which is simply the
 * injection of the values of cell degrees of freedom. The MatrixFree class
 * detects possible structure of the finite element at quadrature points
 * when initializing, which is then automatically used by FEEvaluation for
 * selecting the most appropriate numerical kernel.
 * 

 * 
 * The nonlinear function that we have to evaluate for the time stepping
 * routine includes the value of the function at the present time @p current
 * as well as the value at the previous time step @p old. Both values are
 * passed to the operator in the collection of source vectors @p src, which
 * is simply a <tt>std::vector</tt> of pointers to the actual solution
 * vectors. This construct of collecting several source vectors into one is
 * necessary as the cell loop in @p MatrixFree takes exactly one source and
 * one destination vector, even if we happen to use many vectors like the
 * two in this case. Note that the cell loop accepts any valid class for
 * input and output, which does not only include vectors but general data
 * types.  However, only in case it encounters a
 * LinearAlgebra::distributed::Vector<Number> or a <tt>std::vector</tt>
 * collecting these vectors, it calls functions that exchange ghost data due
 * to MPI at the beginning and the end of the loop. In the loop over the
 * cells, we first have to read in the values in the vectors related to the
 * local values.  Then, we evaluate the value and the gradient of the
 * current solution vector and the values of the old vector at the
 * quadrature points. Next, we combine the terms in the scheme in the loop
 * over the quadrature points. Finally, we integrate the result against the
 * test function and accumulate the result to the global solution vector @p
 * dst.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void SineGordonOperation<dim, fe_degree>::local_apply(
 *     const MatrixFree<dim> &                                          data,
 *     LinearAlgebra::distributed::Vector<double> &                     dst,
 *     const std::vector<LinearAlgebra::distributed::Vector<double> *> &src,
 *     const std::pair<unsigned int, unsigned int> &cell_range) const
 *   {
 *     AssertDimension(src.size(), 2);
 *     FEEvaluation<dim, fe_degree> current(data), old(data);
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         current.reinit(cell);
 *         old.reinit(cell);
 * 
 *         current.read_dof_values(*src[0]);
 *         old.read_dof_values(*src[1]);
 * 
 *         current.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
 *         old.evaluate(EvaluationFlags::values);
 * 
 *         for (unsigned int q = 0; q < current.n_q_points; ++q)
 *           {
 *             const VectorizedArray<double> current_value = current.get_value(q);
 *             const VectorizedArray<double> old_value     = old.get_value(q);
 * 
 *             current.submit_value(2. * current_value - old_value -
 *                                    delta_t_sqr * std::sin(current_value),
 *                                  q);
 *             current.submit_gradient(-delta_t_sqr * current.get_gradient(q), q);
 *           }
 * 
 *         current.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
 *         current.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonOperationapply"></a> 
 * <h4>SineGordonOperation::apply</h4>
 * 

 * 
 * This function performs the time stepping routine based on the cell-local
 * strategy. Note that we need to set the destination vector to zero before
 * we add the integral contributions of the current time step (via the
 * FEEvaluation::distribute_local_to_global() call). In this tutorial, we
 * let the cell-loop do the zero operation via the fifth `true` argument
 * passed to MatrixFree::cell_loop. The loop can schedule the zero operation
 * closer to the operations on vector entries for supported vector entries,
 * thereby possibly increasing data locality (the vector entries that first
 * get zeroed are later re-used in the `distribute_local_to_global()`
 * call). The structure of the cell loop is implemented in the cell finite
 * element operator class. On each cell it applies the routine defined as
 * the <code>local_apply()</code> method of the class
 * <code>SineGordonOperation</code>, i.e., <code>this</code>. One could also
 * provide a function with the same signature that is not part of a
 * class. Finally, the result of the integration is multiplied by the
 * inverse mass matrix.
 * 
 * @code
 *   template <int dim, int fe_degree>
 *   void SineGordonOperation<dim, fe_degree>::apply(
 *     LinearAlgebra::distributed::Vector<double> &                     dst,
 *     const std::vector<LinearAlgebra::distributed::Vector<double> *> &src) const
 *   {
 *     data.cell_loop(
 *       &SineGordonOperation<dim, fe_degree>::local_apply, this, dst, src, true);
 *     dst.scale(inv_mass_matrix);
 *   }
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
 * We define a time-dependent function that is used as initial
 * value. Different solutions can be obtained by varying the starting
 * time. This function, taken from step-25, would represent an analytic
 * solution in 1D for all times, but is merely used for setting some
 * starting solution of interest here. More elaborate choices that could
 * test the convergence of this program are given in step-25.
 * 
 * @code
 *   template <int dim>
 *   class InitialCondition : public Function<dim>
 *   {
 *   public:
 *     InitialCondition(const unsigned int n_components = 1,
 *                      const double       time         = 0.)
 *       : Function<dim>(n_components, time)
 *     {}
 *     virtual double value(const Point<dim> &p,
 *                          const unsigned int /*component*/) const override
 *     {
 *       double t = this->get_time();
 * 
 *       const double m  = 0.5;
 *       const double c1 = 0.;
 *       const double c2 = 0.;
 *       const double factor =
 *         (m / std::sqrt(1. - m * m) * std::sin(std::sqrt(1. - m * m) * t + c2));
 *       double result = 1.;
 *       for (unsigned int d = 0; d < dim; ++d)
 *         result *= -4. * std::atan(factor / std::cosh(m * p[d] + c1));
 *       return result;
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemclass"></a> 
 * <h3>SineGordonProblem class</h3>
 * 

 * 
 * This is the main class that builds on the class in step-25.  However, we
 * replaced the SparseMatrix<double> class by the MatrixFree class to store
 * the geometry data. Also, we use a distributed triangulation in this
 * example.
 * 
 * @code
 *   template <int dim>
 *   class SineGordonProblem
 *   {
 *   public:
 *     SineGordonProblem();
 *     void run();
 * 
 *   private:
 *     ConditionalOStream pcout;
 * 
 *     void make_grid_and_dofs();
 *     void output_results(const unsigned int timestep_number);
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     MappingQ1<dim> mapping;
 * 
 *     AffineConstraints<double> constraints;
 *     IndexSet                  locally_relevant_dofs;
 * 
 *     MatrixFree<dim, double> matrix_free_data;
 * 
 *     LinearAlgebra::distributed::Vector<double> solution, old_solution,
 *       old_old_solution;
 * 
 *     const unsigned int n_global_refinements;
 *     double             time, time_step;
 *     const double       final_time;
 *     const double       cfl_number;
 *     const unsigned int output_timestep_skip;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemSineGordonProblem"></a> 
 * <h4>SineGordonProblem::SineGordonProblem</h4>
 * 

 * 
 * This is the constructor of the SineGordonProblem class. The time interval
 * and time step size are defined here. Moreover, we use the degree of the
 * finite element that we defined at the top of the program to initialize a
 * FE_Q finite element based on Gauss-Lobatto support points. These points
 * are convenient because in conjunction with a QGaussLobatto quadrature
 * rule of the same order they give a diagonal mass matrix without
 * compromising accuracy too much (note that the integration is inexact,
 * though), see also the discussion in the introduction. Note that FE_Q
 * selects the Gauss-Lobatto nodal points by default due to their improved
 * conditioning versus equidistant points. To make things more explicit, we
 * state the selection of the nodal points nonetheless.
 * 
 * @code
 *   template <int dim>
 *   SineGordonProblem<dim>::SineGordonProblem()
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *     ,
 * #ifdef DEAL_II_WITH_P4EST
 *     triangulation(MPI_COMM_WORLD)
 *     ,
 * #endif
 *     fe(QGaussLobatto<1>(fe_degree + 1))
 *     , dof_handler(triangulation)
 *     , n_global_refinements(10 - 2 * dim)
 *     , time(-10)
 *     , time_step(10.)
 *     , final_time(10.)
 *     , cfl_number(.1 / fe_degree)
 *     , output_timestep_skip(200)
 *   {}
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemmake_grid_and_dofs"></a> 
 * <h4>SineGordonProblem::make_grid_and_dofs</h4>
 * 

 * 
 * As in step-25 this functions sets up a cube grid in <code>dim</code>
 * dimensions of extent $[-15,15]$. We refine the mesh more in the center of
 * the domain since the solution is concentrated there. We first refine all
 * cells whose center is within a radius of 11, and then refine once more
 * for a radius 6.  This simple ad hoc refinement could be done better by
 * adapting the mesh to the solution using error estimators during the time
 * stepping as done in other example programs, and using
 * parallel::distributed::SolutionTransfer to transfer the solution to the
 * new mesh.
 * 
 * @code
 *   template <int dim>
 *   void SineGordonProblem<dim>::make_grid_and_dofs()
 *   {
 *     GridGenerator::hyper_cube(triangulation, -15, 15);
 *     triangulation.refine_global(n_global_refinements);
 *     {
 *       typename Triangulation<dim>::active_cell_iterator
 *         cell     = triangulation.begin_active(),
 *         end_cell = triangulation.end();
 *       for (; cell != end_cell; ++cell)
 *         if (cell->is_locally_owned())
 *           if (cell->center().norm() < 11)
 *             cell->set_refine_flag();
 *       triangulation.execute_coarsening_and_refinement();
 * 
 *       cell     = triangulation.begin_active();
 *       end_cell = triangulation.end();
 *       for (; cell != end_cell; ++cell)
 *         if (cell->is_locally_owned())
 *           if (cell->center().norm() < 6)
 *             cell->set_refine_flag();
 *       triangulation.execute_coarsening_and_refinement();
 *     }
 * 
 *     pcout << "   Number of global active cells: "
 * #ifdef DEAL_II_WITH_P4EST
 *           << triangulation.n_global_active_cells()
 * #else
 *           << triangulation.n_active_cells()
 * #endif
 *           << std::endl;
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl;
 * 
 * 
 * @endcode
 * 
 * We generate hanging node constraints for ensuring continuity of the
 * solution. As in step-40, we need to equip the constraint matrix with
 * the IndexSet of locally relevant degrees of freedom to avoid it to
 * consume too much memory for big problems. Next, the <code> MatrixFree
 * </code> object for the problem is set up. Note that we specify a
 * particular scheme for shared-memory parallelization (hence one would
 * use multithreading for intra-node parallelism and not MPI; we here
 * choose the standard option &mdash; if we wanted to disable shared
 * memory parallelization even in case where there is more than one TBB
 * thread available in the program, we would choose
 * MatrixFree::AdditionalData::TasksParallelScheme::none). Also note that,
 * instead of using the default QGauss quadrature argument, we supply a
 * QGaussLobatto quadrature formula to enable the desired
 * behavior. Finally, three solution vectors are initialized. MatrixFree
 * expects a particular layout of ghost indices (as it handles index
 * access in MPI-local numbers that need to match between the vector and
 * MatrixFree), so we just ask it to initialize the vectors to be sure the
 * ghost exchange is properly handled.
 * 
 * @code
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     constraints.close();
 * 
 *     typename MatrixFree<dim>::AdditionalData additional_data;
 *     additional_data.tasks_parallel_scheme =
 *       MatrixFree<dim>::AdditionalData::TasksParallelScheme::partition_partition;
 * 
 *     matrix_free_data.reinit(mapping,
 *                             dof_handler,
 *                             constraints,
 *                             QGaussLobatto<1>(fe_degree + 1),
 *                             additional_data);
 * 
 *     matrix_free_data.initialize_dof_vector(solution);
 *     old_solution.reinit(solution);
 *     old_old_solution.reinit(solution);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemoutput_results"></a> 
 * <h4>SineGordonProblem::output_results</h4>
 * 

 * 
 * This function prints the norm of the solution and writes the solution
 * vector to a file. The norm is standard (except for the fact that we need
 * to accumulate the norms over all processors for the parallel grid which
 * we do via the VectorTools::compute_global_error() function), and the
 * second is similar to what we did in step-40 or step-37. Note that we can
 * use the same vector for output as the one used during computations: The
 * vectors in the matrix-free framework always provide full information on
 * all locally owned cells (this is what is needed in the local evaluations,
 * too), including ghost vector entries on these cells. This is the only
 * data that is needed in the VectorTools::integrate_difference() function
 * as well as in DataOut. The only action to take at this point is to make
 * sure that the vector updates its ghost values before we read from
 * them, and to reset ghost values once done. This is a feature present only
 * in the LinearAlgebra::distributed::Vector class. Distributed vectors with
 * PETSc and Trilinos, on the other hand, need to be copied to special
 * vectors including ghost values (see the relevant section in step-40). If
 * we also wanted to access all degrees of freedom on ghost cells (e.g. when
 * computing error estimators that use the jump of solution over cell
 * boundaries), we would need more information and create a vector
 * initialized with locally relevant dofs just as in step-40. Observe also
 * that we need to distribute constraints for output - they are not filled
 * during computations (rather, they are interpolated on the fly in the
 * matrix-free method FEEvaluation::read_dof_values()).
 * 
 * @code
 *   template <int dim>
 *   void
 *   SineGordonProblem<dim>::output_results(const unsigned int timestep_number)
 *   {
 *     constraints.distribute(solution);
 * 
 *     Vector<float> norm_per_cell(triangulation.n_active_cells());
 *     solution.update_ghost_values();
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Functions::ZeroFunction<dim>(),
 *                                       norm_per_cell,
 *                                       QGauss<dim>(fe_degree + 1),
 *                                       VectorTools::L2_norm);
 *     const double solution_norm =
 *       VectorTools::compute_global_error(triangulation,
 *                                         norm_per_cell,
 *                                         VectorTools::L2_norm);
 * 
 *     pcout << "   Time:" << std::setw(8) << std::setprecision(3) << time
 *           << ", solution norm: " << std::setprecision(5) << std::setw(7)
 *           << solution_norm << std::endl;
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches(mapping);
 * 
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", timestep_number, MPI_COMM_WORLD, 3);
 * 
 *     solution.zero_out_ghost_values();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemrun"></a> 
 * <h4>SineGordonProblem::run</h4>
 * 

 * 
 * This function is called by the main function and steps into the
 * subroutines of the class.
 *   

 * 
 * After printing some information about the parallel setup, the first
 * action is to set up the grid and the cell operator. Then, the time step
 * is computed from the CFL number given in the constructor and the finest
 * mesh size. The finest mesh size is computed as the diameter of the last
 * cell in the triangulation, which is the last cell on the finest level of
 * the mesh. This is only possible for meshes where all elements on a level
 * have the same size, otherwise, one needs to loop over all cells. Note
 * that we need to query all the processors for their finest cell since
 * not all processors might hold a region where the mesh is at the finest
 * level. Then, we readjust the time step a little to hit the final time
 * exactly.
 * 
 * @code
 *   template <int dim>
 *   void SineGordonProblem<dim>::run()
 *   {
 *     {
 *       pcout << "Number of MPI ranks:            "
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
 *       pcout << "Number of threads on each rank: "
 *             << MultithreadInfo::n_threads() << std::endl;
 *       const unsigned int n_vect_doubles = VectorizedArray<double>::size();
 *       const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;
 *       pcout << "Vectorization over " << n_vect_doubles
 *             << " doubles = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl
 *             << std::endl;
 *     }
 *     make_grid_and_dofs();
 * 
 *     const double local_min_cell_diameter =
 *       triangulation.last()->diameter() / std::sqrt(dim);
 *     const double global_min_cell_diameter =
 *       -Utilities::MPI::max(-local_min_cell_diameter, MPI_COMM_WORLD);
 *     time_step = cfl_number * global_min_cell_diameter;
 *     time_step = (final_time - time) / (int((final_time - time) / time_step));
 *     pcout << "   Time step size: " << time_step
 *           << ", finest cell: " << global_min_cell_diameter << std::endl
 *           << std::endl;
 * 
 * @endcode
 * 
 * Next the initial value is set. Since we have a two-step time stepping
 * method, we also need a value of the solution at time-time_step. For
 * accurate results, one would need to compute this from the time
 * derivative of the solution at initial time, but here we ignore this
 * difficulty and just set it to the initial value function at that
 * artificial time.
 * 

 * 
 * We then go on by writing the initial state to file and collecting
 * the two starting solutions in a <tt>std::vector</tt> of pointers that
 * get later consumed by the SineGordonOperation::apply() function. Next,
 * an instance of the <code> SineGordonOperation class </code> based on
 * the finite element degree specified at the top of this file is set up.
 * 
 * @code
 *     VectorTools::interpolate(mapping,
 *                              dof_handler,
 *                              InitialCondition<dim>(1, time),
 *                              solution);
 *     VectorTools::interpolate(mapping,
 *                              dof_handler,
 *                              InitialCondition<dim>(1, time - time_step),
 *                              old_solution);
 *     output_results(0);
 * 
 *     std::vector<LinearAlgebra::distributed::Vector<double> *>
 *       previous_solutions({&old_solution, &old_old_solution});
 * 
 *     SineGordonOperation<dim, fe_degree> sine_gordon_op(matrix_free_data,
 *                                                        time_step);
 * 
 * @endcode
 * 
 * Now loop over the time steps. In each iteration, we shift the solution
 * vectors by one and call the `apply` function of the
 * `SineGordonOperator` class. Then, we write the solution to a file. We
 * clock the wall times for the computational time needed as wall as the
 * time needed to create the output and report the numbers when the time
 * stepping is finished.
 *     

 * 
 * Note how this shift is implemented: We simply call the swap method on
 * the two vectors which swaps only some pointers without the need to copy
 * data around, a relatively expensive operation within an explicit time
 * stepping method. Let us see what happens in more detail: First, we
 * exchange <code>old_solution</code> with <code>old_old_solution</code>,
 * which means that <code>old_old_solution</code> gets
 * <code>old_solution</code>, which is what we expect. Similarly,
 * <code>old_solution</code> gets the content from <code>solution</code>
 * in the next step. After this, <code>solution</code> holds
 * <code>old_old_solution</code>, but that will be overwritten during this
 * step.
 * 
 * @code
 *     unsigned int timestep_number = 1;
 * 
 *     Timer  timer;
 *     double wtime       = 0;
 *     double output_time = 0;
 *     for (time += time_step; time <= final_time;
 *          time += time_step, ++timestep_number)
 *       {
 *         timer.restart();
 *         old_old_solution.swap(old_solution);
 *         old_solution.swap(solution);
 *         sine_gordon_op.apply(solution, previous_solutions);
 *         wtime += timer.wall_time();
 * 
 *         timer.restart();
 *         if (timestep_number % output_timestep_skip == 0)
 *           output_results(timestep_number / output_timestep_skip);
 * 
 *         output_time += timer.wall_time();
 *       }
 *     timer.restart();
 *     output_results(timestep_number / output_timestep_skip + 1);
 *     output_time += timer.wall_time();
 * 
 *     pcout << std::endl
 *           << "   Performed " << timestep_number << " time steps." << std::endl;
 * 
 *     pcout << "   Average wallclock time per time step: "
 *           << wtime / timestep_number << "s" << std::endl;
 * 
 *     pcout << "   Spent " << output_time << "s on output and " << wtime
 *           << "s on computations." << std::endl;
 *   }
 * } // namespace Step48
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * As in step-40, we initialize MPI at the start of the program. Since we will
 * in general mix MPI parallelization with threads, we also set the third
 * argument in MPI_InitFinalize that controls the number of threads to an
 * invalid number, which means that the TBB library chooses the number of
 * threads automatically, typically to the number of available cores in the
 * system. As an alternative, you can also set this number manually if you
 * want to set a specific number of threads (e.g. when MPI-only is required).
 * 
 * @code
 * int main(int argc, char **argv)
 * {
 *   using namespace Step48;
 *   using namespace dealii;
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *     argc, argv, numbers::invalid_unsigned_int);
 * 
 *   try
 *     {
 *       SineGordonProblem<dimension> sg_problem;
 *       sg_problem.run();
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
examples/step-48/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Comparisonwithasparsematrix"></a><h3>Comparison with a sparse matrix</h3>


为了证明使用MatrixFree类而不是标准的 <code>deal.II</code> 汇编例程来评估旧时间步长的信息的好处，我们研究了代码在非自适应网格上的一个简单串行运行。由于很多时间花在评估正弦函数上，我们不仅显示了完整的正弦-戈登方程的数字，还显示了波浪方程（正弦-戈登方程中跳过的正弦项）的数字。我们同时使用二阶和四阶元素。结果总结在下表中。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3">wave equation</th>
    <th colspan="2">sine-Gordon</th>
  </tr>
  <tr>
    <th>&nbsp;</th>
    <th>MF</th>
    <th>SpMV</th>
    <th>dealii</th>
    <th>MF</th>
    <th>dealii</th>
  </tr>
  <tr>
    <td>2D, $\mathcal{Q}_2$</td>
    <td align="right"> 0.0106</td>
    <td align="right"> 0.00971</td>
    <td align="right"> 0.109</td>
    <td align="right"> 0.0243</td>
    <td align="right"> 0.124</td>
  </tr>
  <tr>
    <td>2D, $\mathcal{Q}_4$</td>
    <td align="right"> 0.0328</td>
    <td align="right"> 0.0706</td>
    <td align="right"> 0.528</td>
    <td align="right"> 0.0714</td>
    <td align="right"> 0.502</td>
   </tr>
   <tr>
    <td>3D, $\mathcal{Q}_2$</td>
    <td align="right"> 0.0151</td>
    <td align="right"> 0.0320</td>
    <td align="right"> 0.331</td>
    <td align="right"> 0.0376</td>
    <td align="right"> 0.364</td>
   </tr>
   <tr>
    <td>3D, $\mathcal{Q}_4$</td>
    <td align="right"> 0.0918</td>
    <td align="right"> 0.844</td>
    <td align="right"> 6.83</td>
    <td align="right"> 0.194</td>
    <td align="right"> 6.95</td>
   </tr>
</table> 

很明显，无矩阵代码远远超过了deal.II中的标准汇编程序。在三维和四阶元素中，一个运算符的评估速度也几乎是稀疏矩阵-向量乘积的十倍。

<a name="Parallelrunin2Dand3D"></a><h3>Parallel run in 2D and 3D</h3>


我们从一个具有12个核心/24个线程的工作站（一个英特尔至强E5-2687W v4 CPU运行在3.2 GHz，启用了超线程）上获得的程序输出开始，以发布模式运行程序。

@code
\$ make run
Number of MPI ranks:            1
Number of threads on each rank: 24
Vectorization over 4 doubles = 256 bits (AVX)


   Number of global active cells: 15412
   Number of degrees of freedom: 249065
   Time step size: 0.00292997, finest cell: 0.117188


   Time:     -10, solution norm:  9.5599
   Time:   -9.41, solution norm:  17.678
   Time:   -8.83, solution norm:  23.504
   Time:   -8.24, solution norm:    27.5
   Time:   -7.66, solution norm:  29.513
   Time:   -7.07, solution norm:  29.364
   Time:   -6.48, solution norm:   27.23
   Time:    -5.9, solution norm:  23.527
   Time:   -5.31, solution norm:  18.439
   Time:   -4.73, solution norm:  11.935
   Time:   -4.14, solution norm:  5.5284
   Time:   -3.55, solution norm:  8.0354
   Time:   -2.97, solution norm:  14.707
   Time:   -2.38, solution norm:      20
   Time:    -1.8, solution norm:  22.834
   Time:   -1.21, solution norm:  22.771
   Time:  -0.624, solution norm:  20.488
   Time: -0.0381, solution norm:  16.697
   Time:   0.548, solution norm:  11.221
   Time:    1.13, solution norm:  5.3912
   Time:    1.72, solution norm:  8.4528
   Time:    2.31, solution norm:  14.335
   Time:    2.89, solution norm:  18.555
   Time:    3.48, solution norm:  20.894
   Time:    4.06, solution norm:  21.305
   Time:    4.65, solution norm:  19.903
   Time:    5.24, solution norm:  16.864
   Time:    5.82, solution norm:  12.223
   Time:    6.41, solution norm:   6.758
   Time:    6.99, solution norm:  7.2423
   Time:    7.58, solution norm:  12.888
   Time:    8.17, solution norm:  17.273
   Time:    8.75, solution norm:  19.654
   Time:    9.34, solution norm:  19.838
   Time:    9.92, solution norm:  17.964
   Time:      10, solution norm:  17.595


   Performed 6826 time steps.
   Average wallclock time per time step: 0.0013453s
   Spent 14.976s on output and 9.1831s on computations.
@endcode



在3D中，各自的输出看起来像

@code
\$ make run
Number of MPI ranks:            1
Number of threads on each rank: 24
Vectorization over 4 doubles = 256 bits (AVX)


   Number of global active cells: 17592
   Number of degrees of freedom: 1193881
   Time step size: 0.0117233, finest cell: 0.46875


   Time:     -10, solution norm:  29.558
   Time:   -7.66, solution norm:  129.13
   Time:   -5.31, solution norm:  67.753
   Time:   -2.97, solution norm:  79.245
   Time:  -0.621, solution norm:  123.52
   Time:    1.72, solution norm:  43.525
   Time:    4.07, solution norm:  93.285
   Time:    6.41, solution norm:  97.722
   Time:    8.76, solution norm:  36.734
   Time:      10, solution norm:  94.115


   Performed 1706 time steps.
   Average wallclock time per time step: 0.0084542s
   Spent 16.766s on output and 14.423s on computations.
@endcode



一个自由度超过一百万的时间步长需要0.008秒（注意，在求解线性系统时，我们需要许多处理器来达到这样的数字）。

如果我们用一个纯粹的MPI并行化取代线程并行化，时间就会变成。

@code
\$ mpirun -n 24 ./step-48
Number of MPI ranks:            24
Number of threads on each rank: 1
Vectorization over 4 doubles = 256 bits (AVX)
...
   Performed 1706 time steps.
   Average wallclock time per time step: 0.0051747s
   Spent 2.0535s on output and 8.828s on computations.
@endcode



我们观察到输出的急剧加速（这是有道理的，因为输出的大部分代码没有通过线程并行化，而对于MPI则是如此），但低于我们从并行性中期望的12的理论系数。更有趣的是，当从纯线程变量切换到纯MPI变量时，计算也变得更快。这是MatrixFree框架的一个一般观察结果（截至2019年更新此数据）。主要原因是，为实现并行执行而做出的关于冲突单元批处理工作的决定过于悲观：虽然它们确保在不同的线程上不会同时进行相邻单元的工作，但这种保守的设置意味着在相邻单元被触及时，相邻单元的数据也会从缓存中被驱逐。此外，对于给定的具有17592个单元的网格，目前的方案无法为所有24个线程提供一个恒定的负载。

目前的程序还允许将MPI并行化与线程并行化混合起来。在有多个节点的集群上运行程序时，这是最有利的，使用MPI进行节点间并行化，使用线程进行节点内并行化。在上面使用的工作站上，我们可以在超线程区域运行线程（即为12个MPI行列中的每一个使用2个线程）。将MPI与线程混合的一个重要设置是确保将任务适当地分到CPU上。在许多集群上，放置是通过`mpirun/mpiexec`环境自动进行的，或者可以有手动设置。在这里，我们简单地报告了程序的普通版本的运行时间（注意到当适当的分档完成后，事情可以向仅有MPI的程序的时间改进）。

@code
\$ mpirun -n 12 ./step-48
Number of MPI ranks:            12
Number of threads on each rank: 2
Vectorization over 4 doubles = 256 bits (AVX)
...
   Performed 1706 time steps.
   Average wallclock time per time step: 0.0056651s
   Spent 2.5175s on output and 9.6646s on computations.
@endcode






<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序中有几处可以改进，使其更加有效（除了步骤25中讨论的改进边界条件和物理东西）。

 <ul>   <li>  <b>Faster evaluation of sine terms:</b> 从上面平波方程和正弦-戈登方程的比较中可以明显看出，正弦项的评估在有限元算子应用的总时间中占主导地位。这有几个原因。首先，VectorizedArray场的deal.II正弦计算没有被矢量化（与算子应用的其他部分相反）。这可以通过将正弦计算交给一个具有矢量化正弦计算的库来解决，比如英特尔的数学内核库（MKL）。通过使用MKL中的函数 <code>vdSin</code> ，该程序在二维中使用了一半的计算时间，在三维中使用了40%的时间。另一方面，正弦计算在结构上要比其他本地操作中的加法和乘法等简单算术操作复杂得多。

    <li>  <b>Higher order time stepping:</b> 虽然该实现允许空间部分的任意顺序（通过调整有限元的程度），但时间步进方案是一个标准的二阶跃迁方案。由于波的传播问题的解通常是非常平滑的，所以误差很可能被时间步进部分所支配。当然，这可以通过使用较小的时间步长（在固定的空间分辨率下）来解决，但如果使用高阶时间步长也会更有效率。虽然对于一阶系统来说，这样做是很简单的（使用一些高阶的Runge&ndash;Kutta方案，可能结合像<a
  href="http://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method">Dormand&ndash;Prince
  method</a>那样的自适应时间步长选择），但对于二阶公式来说，这更具挑战性。至少在有限差分社区，人们通常使用PDE来寻找改善时间误差的空间修正项。

 </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-48.cc"
*/
