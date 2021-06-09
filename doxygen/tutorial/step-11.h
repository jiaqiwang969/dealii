/**
@page step_11 The step-11 tutorial program
This tutorial depends on step-10.

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
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-11/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


我们要考虑的问题是只带诺伊曼边界条件的拉普拉斯问题的解决方案。

@f{eqnarray*}


  -\Delta u &=& f \qquad \mathrm{in}\ \Omega,
  \\
  \partial_n u &=& g \qquad \mathrm{on}\ \partial\Omega.


@f}

众所周知，如果这个问题要有一个解决方案，那么力需要满足兼容性条件

@f[
  \int_\Omega f\; dx + \int_{\partial\Omega} g\; ds = 0.


@f]

我们将考虑这样的特殊情况： $\Omega$ 是围绕原点的半径为1的圆，而 $f=-2$  ， $g=1$  。这种选择满足了兼容性条件。

兼容性条件允许上述方程的解，但它仍然保留了一个模糊性：因为只有解的导数出现在方程中，解只确定到一个常数。出于这个原因，我们必须为数字解提出另一个条件，以固定这个常数。

对于这一点，有多种可能性。<ol>  <li>  将离散化的一个节点固定为零或任何其他固定值。   这相当于一个附加条件  $u_h(x_0)=0$  。虽然这是常见的做法，但不一定是个好主意，因为我们知道拉普拉斯方程的解只在  $H^1$  中，这不允许定义点值，因为它不是连续函数的子集。因此，即使固定一个节点对离散函数来说是允许的，但对连续函数来说是不允许的，在数值解的这一点上，人们常常可以看到由此产生的错误尖峰。

 <li>  将域上的均值固定为零或任何其他值。这在连续水平上是允许的，因为 $H^1(\Omega)\subset L^1(\Omega)$ 由Sobolev不等式决定，因此在离散水平上也是允许的，因为我们那里只考虑 $H^1$ 的子集。

 <li>  将域的边界上的均值固定为零或任何其他值。这在连续水平上也是允许的，因为 $H^{1/2}(\partial\Omega)\subset L^1(\partial\Omega)$  ，同样由Sobolev的不等式。   </ol>  我们将选择最后一种可能性，因为我们想用它来演示另一种技术。

虽然这描述了要解决的问题，但我们仍然要弄清楚如何实现它。基本上，除了额外的均值约束，我们已经多次解决了这个问题，使用的是迪里希特边界值，我们只需要放弃对迪里希特边界节点的处理。高阶映射的使用也是相当琐碎的，我们会在使用它的各个地方进行解释；在几乎所有可以想象的情况下，你只会把描述映射的对象视为一个黑盒子，你不需要担心，因为它们的唯一用途似乎是被传递到库的深处，在那里函数知道如何处理它们（即在 <code>FEValues</code> 类及其后代）。

这个程序中的棘手之处在于对均值约束的使用。幸运的是，库中有一个知道如何处理这种约束的类，我们已经经常使用它了，没有提到它的通用性。请注意，如果我们假设边界节点沿边界的间隔是相等的，那么均值约束

@f[
  \int_{\partial \Omega} u(x) \; ds = 0


@f]

可写为

@f[
  \sum_{i\in\partial\Omega_h} u_i = 0,


@f]

其中总和应贯穿位于计算域边界上的所有自由度指数。让我们用 $i_0$ 表示边界上数字最小的指数（或任何其他方便选择的指数），那么这个约束也可以用以下方式表示

@f[
  u_{i_0} = \sum_{i\in\partial\Omega_h\backslash i_0} -u_i.


@f]

幸运的是，这正是AffineConstraints类所设计的约束形式。请注意，我们在之前的几个例子中使用了这个类来表示悬空节点的约束，它也有这种形式：在这里，中间的顶点应具有相邻顶点的平均值。一般来说，AffineConstraints类被设计用来处理以下形式的仿生约束

@f[
  CU = b


@f]

其中 $C$ 表示一个矩阵， $b$ 表示一个向量， $U$ 是节点值的向量。在这种情况下，由于 $C$ 代表一个同质约束， $b$ 是零向量。

在这个例子中，沿边界的平均值允许这样的表示， $C$ 是一个只有一行的矩阵（即只有一个约束条件）。在实现中，我们将创建一个AffineConstraints对象，添加一个参考第一个边界节点 $i_0$ 的约束（即给矩阵添加另一行），并插入所有其他节点贡献的权重，在这个例子中刚好是 $-1$  。

稍后，我们将使用这个对象来消除线性方程组中的第一个边界节点，将其还原为一个没有常数偏移值的解。实施过程中的一个问题是，这个节点的明确消除会导致矩阵中出现一些额外的元素，我们事先不知道这些元素的位置，也不知道矩阵的每一行中会有多少额外的条目。我们将展示我们如何使用一个中间对象来解决这个问题。

但现在开始实施解决这个问题的方案......


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * As usual, the program starts with a rather long list of include files which
 * you are probably already used to by now:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/table_handler.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * Just this one is new: it declares a class
 * DynamicSparsityPattern, which we will use and explain
 * further down below.
 * 
 * @code
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * 
 * @endcode
 * 
 * We will make use of the std::find algorithm of the C++ standard library, so
 * we have to include the following file for its declaration:
 * 
 * @code
 * #include <algorithm>
 * #include <iostream>
 * #include <iomanip>
 * #include <cmath>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step11
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Then we declare a class which represents the solution of a Laplace
 * problem. As this example program is based on step-5, the class looks
 * rather the same, with the sole structural difference that the functions
 * <code>assemble_system</code> now calls <code>solve</code> itself, and is
 * thus called <code>assemble_and_solve</code>, and that the output function
 * was dropped since the solution function is so boring that it is not worth
 * being viewed.
 *   

 * 
 * The only other noteworthy change is that the constructor takes a value
 * representing the polynomial degree of the mapping to be used later on,
 * and that it has another member variable representing exactly this
 * mapping. In general, this variable will occur in real applications at the
 * same places where the finite element is declared or used.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem(const unsigned int mapping_degree);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_and_solve();
 *     void solve();
 *     void write_high_order_mesh(const unsigned cycle);
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 *     MappingQ<dim>      mapping;
 * 
 *     SparsityPattern           sparsity_pattern;
 *     SparseMatrix<double>      system_matrix;
 *     AffineConstraints<double> mean_value_constraints;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     TableHandler output_table;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * Construct such an object, by initializing the variables. Here, we use
 * linear finite elements (the argument to the <code>fe</code> variable
 * denotes the polynomial degree), and mappings of given order. Print to
 * screen what we are about to do.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem(const unsigned int mapping_degree)
 *     : fe(1)
 *     , dof_handler(triangulation)
 *     , mapping(mapping_degree)
 *   {
 *     std::cout << "Using mapping with degree " << mapping_degree << ":"
 *               << std::endl
 *               << "============================" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The first task is to set up the variables for this problem. This includes
 * generating a valid <code>DoFHandler</code> object, as well as the
 * sparsity patterns for the matrix, and the object representing the
 * constraints that the mean value of the degrees of freedom on the boundary
 * be zero.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 * @endcode
 * 
 * The first task is trivial: generate an enumeration of the degrees of
 * freedom, and initialize solution and right hand side vector to their
 * correct sizes:
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 * @endcode
 * 
 * The next task is to construct the object representing the constraint that
 * the mean value of the degrees of freedom on the boundary shall be
 * zero. For this, we first want a list of those nodes that are actually
 * at the boundary. The <code>DoFTools</code> namespace has a function
 * that returns an IndexSet object that contains the indices of all those
 * degrees of freedom that are at the boundary.
 *     

 * 
 * Once we have this index set, we wanted to know which is the first
 * index corresponding to a degree of freedom on the boundary. We need
 * this because we wanted to constrain one of the nodes on the boundary by
 * the values of all other DoFs on the boundary. To get the index of this
 * "first" degree of freedom is easy enough using the IndexSet class:
 * 
 * @code
 *     const IndexSet boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler);
 * 
 *     const types::global_dof_index first_boundary_dof =
 *       boundary_dofs.nth_index_in_set(0);
 * 
 * @endcode
 * 
 * Then generate a constraints object with just this one constraint. First
 * clear all previous content (which might reside there from the previous
 * computation on a once coarser grid), then add this one line
 * constraining the <code>first_boundary_dof</code> to the sum of other
 * boundary DoFs each with weight -1. Finally, close the constraints
 * object, i.e. do some internal bookkeeping on it for faster processing
 * of what is to come later:
 * 
 * @code
 *     mean_value_constraints.clear();
 *     mean_value_constraints.add_line(first_boundary_dof);
 *     for (types::global_dof_index i : boundary_dofs)
 *       if (i != first_boundary_dof)
 *         mean_value_constraints.add_entry(first_boundary_dof, i, -1);
 *     mean_value_constraints.close();
 * 
 * @endcode
 * 
 * Next task is to generate a sparsity pattern. This is indeed a tricky
 * task here. Usually, we just call
 * <code>DoFTools::make_sparsity_pattern</code> and condense the result
 * using the hanging node constraints. We have no hanging node constraints
 * here (since we only refine globally in this example), but we have this
 * global constraint on the boundary. This poses one severe problem in
 * this context: the <code>SparsityPattern</code> class wants us to state
 * beforehand the maximal number of entries per row, either for all rows
 * or for each row separately. There are functions in the library which
 * can tell you this number in case you just have hanging node constraints
 * (namely DoFHandler::max_couplings_between_dofs), but how is
 * this for the present case? The difficulty arises because the
 * elimination of the constrained degree of freedom requires a number of
 * additional entries in the matrix at places that are not so simple to
 * determine. We would therefore have a problem had we to give a maximal
 * number of entries per row here.
 *     

 * 
 * Since this can be so difficult that no reasonable answer can be given
 * that allows allocation of only a reasonable amount of memory, there is
 * a class DynamicSparsityPattern, that can help us out
 * here. It does not require that we know in advance how many entries rows
 * could have, but allows just about any length. It is thus significantly
 * more flexible in case you do not have good estimates of row lengths,
 * however at the price that building up such a pattern is also
 * significantly more expensive than building up a pattern for which you
 * had information in advance. Nevertheless, as we have no other choice
 * here, we'll just build such an object by initializing it with the
 * dimensions of the matrix and calling another function
 * <code>DoFTools::make_sparsity_pattern</code> to get the sparsity
 * pattern due to the differential operator, then condense it with the
 * constraints object which adds those positions in the sparsity pattern
 * that are required for the elimination of the constraint.
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *     mean_value_constraints.condense(dsp);
 * 
 * @endcode
 * 
 * Finally, once we have the full pattern, we can initialize an object of
 * type <code>SparsityPattern</code> from it and in turn initialize the
 * matrix with it. Note that this is actually necessary, since the
 * DynamicSparsityPattern is so inefficient compared to
 * the <code>SparsityPattern</code> class due to the more flexible data
 * structures it has to use, that we can impossibly base the sparse matrix
 * class on it, but rather need an object of type
 * <code>SparsityPattern</code>, which we generate by copying from the
 * intermediate object.
 *     

 * 
 * As a further sidenote, you will notice that we do not explicitly have
 * to <code>compress</code> the sparsity pattern here. This, of course, is
 * due to the fact that the <code>copy_from</code> function generates a
 * compressed object right from the start, to which you cannot add new
 * entries anymore. The <code>compress</code> call is therefore implicit
 * in the <code>copy_from</code> call.
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp);
 *     system_matrix.reinit(sparsity_pattern);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The next function then assembles the linear system of equations, solves
 * it, and evaluates the solution. This then makes three actions, and we
 * will put them into eight true statements (excluding declaration of
 * variables, and handling of temporary vectors). Thus, this function is
 * something for the very lazy. Nevertheless, the functions called are
 * rather powerful, and through them this function uses a good deal of the
 * whole library. But let's look at each of the steps.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_and_solve()
 *   {
 * @endcode
 * 
 * First, we have to assemble the matrix and the right hand side. In all
 * previous examples, we have investigated various ways how to do this
 * manually. However, since the Laplace matrix and simple right hand sides
 * appear so frequently in applications, the library provides functions
 * for actually doing this for you, i.e. they perform the loop over all
 * cells, setting up the local matrices and vectors, and putting them
 * together for the end result.
 *     

 * 
 * The following are the two most commonly used ones: creation of the
 * Laplace matrix and creation of a right hand side vector from body or
 * boundary forces. They take the mapping object, the
 * <code>DoFHandler</code> object representing the degrees of freedom and
 * the finite element in use, a quadrature formula to be used, and the
 * output object. The function that creates a right hand side vector also
 * has to take a function object describing the (continuous) right hand
 * side function.
 *     

 * 
 * Let us look at the way the matrix and body forces are integrated:
 * 
 * @code
 *     const unsigned int gauss_degree =
 *       std::max(static_cast<unsigned int>(
 *                  std::ceil(1. * (mapping.get_degree() + 1) / 2)),
 *                2U);
 *     MatrixTools::create_laplace_matrix(mapping,
 *                                        dof_handler,
 *                                        QGauss<dim>(gauss_degree),
 *                                        system_matrix);
 *     VectorTools::create_right_hand_side(mapping,
 *                                         dof_handler,
 *                                         QGauss<dim>(gauss_degree),
 *                                         Functions::ConstantFunction<dim>(-2),
 *                                         system_rhs);
 * @endcode
 * 
 * That's quite simple, right?
 *     

 * 
 * Two remarks are in order, though: First, these functions are used in a
 * lot of contexts. Maybe you want to create a Laplace or mass matrix for
 * a vector values finite element; or you want to use the default Q1
 * mapping; or you want to assembled the matrix with a coefficient in the
 * Laplace operator. For this reason, there are quite a large number of
 * variants of these functions in the <code>MatrixCreator</code> and
 * <code>MatrixTools</code> namespaces. Whenever you need a slightly
 * different version of these functions than the ones called above, it is
 * certainly worthwhile to take a look at the documentation and to check
 * whether something fits your needs.
 *     

 * 
 * The second remark concerns the quadrature formula we use: we want to
 * integrate over bilinear shape functions, so we know that we have to use
 * at least an order two Gauss quadrature formula. On the other hand, we
 * want the quadrature rule to have at least the order of the boundary
 * approximation. Since the order of Gauss rule with $r$ points is $2r -
 * 1$, and the order of the boundary approximation using polynomials of
 * degree $p$ is $p+1$, we know that $2r \geq p$. Since r has to be an
 * integer and (as mentioned above) has to be at least $2$, this makes up
 * for the formula above computing <code>gauss_degree</code>.
 *     

 * 
 * Since the generation of the body force contributions to the right hand
 * side vector was so simple, we do that all over again for the boundary
 * forces as well: allocate a vector of the right size and call the right
 * function. The boundary function has constant values, so we can generate
 * an object from the library on the fly, and we use the same quadrature
 * formula as above, but this time of lower dimension since we integrate
 * over faces now instead of cells:
 * 
 * @code
 *     Vector<double> tmp(system_rhs.size());
 *     VectorTools::create_boundary_right_hand_side(
 *       mapping,
 *       dof_handler,
 *       QGauss<dim - 1>(gauss_degree),
 *       Functions::ConstantFunction<dim>(1),
 *       tmp);
 * @endcode
 * 
 * Then add the contributions from the boundary to those from the interior
 * of the domain:
 * 
 * @code
 *     system_rhs += tmp;
 * @endcode
 * 
 * For assembling the right hand side, we had to use two different vector
 * objects, and later add them together. The reason we had to do so is
 * that the <code>VectorTools::create_right_hand_side</code> and
 * <code>VectorTools::create_boundary_right_hand_side</code> functions
 * first clear the output vector, rather than adding up their results to
 * previous contents. This can reasonably be called a design flaw in the
 * library made in its infancy, but unfortunately things are as they are
 * for some time now and it is difficult to change such things that
 * silently break existing code, so we have to live with that.
 * 

 * 
 * Now, the linear system is set up, so we can eliminate the one degree of
 * freedom which we constrained to the other DoFs on the boundary for the
 * mean value constraint from matrix and right hand side vector, and solve
 * the system. After that, distribute the constraints again, which in this
 * case means setting the constrained degree of freedom to its proper
 * value
 * 
 * @code
 *     mean_value_constraints.condense(system_matrix);
 *     mean_value_constraints.condense(system_rhs);
 * 
 *     solve();
 *     mean_value_constraints.distribute(solution);
 * 
 * @endcode
 * 
 * Finally, evaluate what we got as solution. As stated in the
 * introduction, we are interested in the H1 semi-norm of the
 * solution. Here, as well, we have a function in the library that does
 * this, although in a slightly non-obvious way: the
 * <code>VectorTools::integrate_difference</code> function integrates the
 * norm of the difference between a finite element function and a
 * continuous function. If we therefore want the norm of a finite element
 * field, we just put the continuous function to zero. Note that this
 * function, just as so many other ones in the library as well, has at
 * least two versions, one which takes a mapping as argument (which we
 * make us of here), and the one which we have used in previous examples
 * which implicitly uses <code>MappingQ1</code>.  Also note that we take a
 * quadrature formula of one degree higher, in order to avoid
 * superconvergence effects where the solution happens to be especially
 * close to the exact solution at certain points (we don't know whether
 * this might be the case here, but there are cases known of this, and we
 * just want to make sure):
 * 
 * @code
 *     Vector<float> norm_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(mapping,
 *                                       dof_handler,
 *                                       solution,
 *                                       Functions::ZeroFunction<dim>(),
 *                                       norm_per_cell,
 *                                       QGauss<dim>(gauss_degree + 1),
 *                                       VectorTools::H1_seminorm);
 * @endcode
 * 
 * Then, the function just called returns its results as a vector of
 * values each of which denotes the norm on one cell. To get the global
 * norm, we do the following:
 * 
 * @code
 *     const double norm =
 *       VectorTools::compute_global_error(triangulation,
 *                                         norm_per_cell,
 *                                         VectorTools::H1_seminorm);
 * 
 * @endcode
 * 
 * Last task -- generate output:
 * 
 * @code
 *     output_table.add_value("cells", triangulation.n_active_cells());
 *     output_table.add_value("|u|_1", norm);
 *     output_table.add_value("error",
 *                            std::fabs(norm - std::sqrt(3.14159265358 / 2)));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The following function solving the linear system of equations is copied
 * from step-5 and is explained there in some detail:
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
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
 * Next, we write the solution as well as the
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
 * page</a> for more details). Secondly, we need to make sure that the mapping
 * is passed to the DataOut::build_patches() method. Finally, the DataOut
 * class only prints curved faces for <i>boundary</i> cells by default, so we
 * need to ensure that also inner cells are printed in a curved representation
 * via the mapping.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::write_high_order_mesh(const unsigned cycle)
 *   {
 *     DataOut<dim> data_out;
 * 
 *     DataOutBase::VtkFlags flags;
 *     flags.write_higher_order_cells = true;
 *     data_out.set_flags(flags);
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 * 
 *     data_out.build_patches(mapping,
 *                            mapping.get_degree(),
 *                            DataOut<dim>::curved_inner_cells);
 * 
 *     std::ofstream file("solution-c=" + std::to_string(cycle) +
 *                        ".p=" + std::to_string(mapping.get_degree()) + ".vtu");
 * 
 *     data_out.write_vtu(file);
 *   }
 * 
 * 
 * @endcode
 * 
 * Finally the main function controlling the different steps to be
 * performed. Its content is rather straightforward, generating a
 * triangulation of a circle, associating a boundary to it, and then doing
 * several cycles on subsequently finer grids. Note again that we have put
 * mesh refinement into the loop header; this may be something for a test
 * program, but for real applications you should consider that this implies
 * that the mesh is refined after the loop is executed the last time since
 * the increment clause (the last part of the three-parted loop header) is
 * executed before the comparison part (the second one), which may be rather
 * costly if the mesh is already quite refined. In that case, you should
 * arrange code such that the mesh is not further refined after the last
 * loop run (or you should do it at the beginning of each run except for the
 * first one).
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     GridGenerator::hyper_ball(triangulation);
 * 
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle)
 *       {
 *         setup_system();
 *         assemble_and_solve();
 *         write_high_order_mesh(cycle);
 * 
 *         triangulation.refine_global();
 *       }
 * 
 * @endcode
 * 
 * After all the data is generated, write a table of results to the
 * screen:
 * 
 * @code
 *     output_table.set_precision("|u|_1", 6);
 *     output_table.set_precision("error", 6);
 *     output_table.write_text(std::cout);
 *     std::cout << std::endl;
 *   }
 * } // namespace Step11
 * 
 * 
 * 
 * @endcode
 * 
 * Finally the main function. It's structure is the same as that used in
 * several of the previous examples, so probably needs no more explanation.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       std::cout.precision(5);
 * 
 * @endcode
 * 
 * This is the main loop, doing the computations with mappings of linear
 * through cubic mappings. Note that since we need the object of type
 * <code>LaplaceProblem@<2@></code> only once, we do not even name it,
 * but create an unnamed such object and call the <code>run</code>
 * function of it, subsequent to which it is immediately destroyed
 * again.
 * 
 * @code
 *       for (unsigned int mapping_degree = 1; mapping_degree <= 3;
 *            ++mapping_degree)
 *         Step11::LaplaceProblem<2>(mapping_degree).run();
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
 *     };
 * 
 *   return 0;
 * }
 * @endcode
examples/step-11/doc/results.dox



<a name="Results"></a><h1>Results</h1>


这是该程序的输出结果。

@code
Using mapping with degree 1:
============================
cells  |u|_1    error
    5 0.680402 0.572912
   20 1.088141 0.165173
   80 1.210142 0.043172
  320 1.242375 0.010939
 1280 1.250569 0.002745
 5120 1.252627 0.000687


Using mapping with degree 2:
============================
cells  |u|_1    error
    5 1.177062 0.076252
   20 1.228978 0.024336
   80 1.245175 0.008139
  320 1.250881 0.002433
 1280 1.252646 0.000668
 5120 1.253139 0.000175


Using mapping with degree 3:
============================
cells  |u|_1    error
    5 1.193493 0.059821
   20 1.229825 0.023489
   80 1.245221 0.008094
  320 1.250884 0.002430
 1280 1.252646 0.000668
 5120 1.253139 0.000175
@endcode

正如我们所期望的，每个不同的映射的收敛顺序显然是与网格大小成二次方的。不过 <em> 有趣的是，双线性映射（即1度）的误差比高阶映射的误差大三倍以上；因此在这种情况下，使用高阶映射显然是有利的，不是因为它提高了收敛顺序，而只是为了减少收敛顺序前的常数。另一方面，除了在非常粗的网格上，使用立方体映射只能进一步改善结果，幅度微乎其微。

我们还可以通过使用例如ParaView来可视化底层网格。下面的图片显示了不同映射度的初始网格。

 <img src="https://www.dealii.org/images/steps/developer/step-11.cycle_0.png" alt=""> 

显然，当我们从线性映射到二次映射时，这种影响是最明显的。这也反映在上表中给出的误差值中。从二次方到三次方的效果没有那么明显，但由于对圆形边界的描述更加准确，所以还是很明显的。

接下来，让我们看看三次全局细化后的网格

 <img src="https://www.dealii.org/images/steps/developer/step-11.cycle_3.png" alt=""> 

在这里，差异就不那么明显了，特别是对于高阶映射。事实上，在这个细化水平上，表格中报告的误差值在二度和三度的映射之间基本上是相同的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-11.cc"
*/
