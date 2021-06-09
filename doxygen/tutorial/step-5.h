/**
@page step_5 The step-5 tutorial program
This tutorial depends on step-4.

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
        <li><a href="#ThecodeStep5codeclasstemplate">The <code>Step5</code> class template</a>
        <li><a href="#Workingwithnonconstantcoefficients">Working with nonconstant coefficients</a>
        <li><a href="#ThecodeStep5codeclassimplementation">The <code>Step5</code> class implementation</a>
      <ul>
        <li><a href="#Step5Step5">Step5::Step5</a>
        <li><a href="#Step5setup_system">Step5::setup_system</a>
        <li><a href="#Step5assemble_system">Step5::assemble_system</a>
        <li><a href="#Step5solve">Step5::solve</a>
        <li><a href="#Step5output_resultsandsettingoutputflags">Step5::output_results and setting output flags</a>
        <li><a href="#Step5run">Step5::run</a>
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
examples/step-5/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{14} 

这个例子并没有显示出革命性的新东西，但它显示了许多比以前的例子更小的改进，也显示了许多通常可以在有限元程序中发现的小东西。其中包括。   <ul>   <li>  在连续细化的网格上进行计算。至少在数学科学中，在一个层次的网格上计算解是很常见的，以便对解的精度有一个感觉；如果你在一个网格上只有一个解，你通常无法猜测解的精度。此外，deal.II被设计用来支持自适应算法，其中在连续细化的网格上的迭代求解是算法的核心。虽然在这个例子中没有使用自适应网格，但这里为它们奠定了基础。     <li>  在实际应用中，领域经常被自动网格生成器细分为三角形。为了使用它们，从文件中读取粗大的网格是很重要的。在这个例子中，我们将读取一个UCD（非结构化单元数据）格式的粗网格。当这个程序在2000年左右首次编写时，UCD格式是AVS Explorer所使用的--这个程序在当时被合理地广泛使用，但现在已经不再重要了。尽管如此，该文件格式仍然存在，并且仍然被一些程序所理解）。     <li>  有限元程序通常会使用大量的计算时间，所以有时需要进行一些优化。我们将展示其中的一些。     <li>  另一方面，有限元程序往往是相当复杂的，所以调试是一个重要方面。我们通过使用断言来支持安全编程，断言在调试模式下检查参数和%内部状态的有效性，但在优化模式下被删除。(  @dealiiVideoLectureSeeAlso{18})   <li>  关于数学方面，我们展示了如何支持椭圆算子中的可变系数，以及如何对线性方程组使用预处理迭代求解器。   </ul> 

这里要解决的方程式如下。

@f{align*}


  -\nabla \cdot a(\mathbf x) \nabla u(\mathbf x) &= 1 \qquad\qquad & \text{in}\ \Omega,
  \\
  u &= 0 \qquad\qquad & \text{on}\ \partial\Omega.


@f}

如果 $a(\mathbf x)$ 是一个恒定的系数，这就只是泊松方程了。然而，如果它确实是空间可变的，它就是一个更复杂的方程（通常被称为 "扩展泊松方程"）。根据变量 $u$ 所指的内容，它可以模拟各种情况，具有广泛的适用性。

- 如果 $u$ 是电动势，那么 $-a\nabla u$ 是介质中的电流，系数 $a$ 是介质在任何特定点的电导率。在这种情况下，方程的右侧将是电源密度，通常为零或由局部的、类似德尔塔的函数组成）。

- 如果 $u$ 是薄膜的垂直挠度，那么 $a$ 将是对局部刚度的测量。这就是让我们解释下面结果部分所显示的图像的解释。

由于拉普拉斯/泊松方程出现在如此多的场合中，因此除了上面列出的两种解释外，还有许多其他解释。

当组装这个方程的线性系统时，我们需要弱的形式，这里的内容如下。

@f{align*}
  (a \nabla \varphi, \nabla u) &= (\varphi, 1) \qquad \qquad \forall \varphi.


@f}

 <code>assemble_system</code> 函数中的实现紧随其后。


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
 * Again, the first few include files are already known, so we won't comment
 * on them:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * 
 * @endcode
 * 
 * This one is new. We want to read a triangulation from disk, and the class
 * which does this is declared in the following file:
 * 
 * @code
 * #include <deal.II/grid/grid_in.h>
 * 
 * @endcode
 * 
 * We will use a circular domain, and the object describing the boundary of it
 * comes from this file:
 * 
 * @code
 * #include <deal.II/grid/manifold_lib.h>
 * 
 * @endcode
 * 
 * This is C++ ...
 * 
 * @code
 * #include <fstream>
 * #include <iostream>
 * 
 * 
 * @endcode
 * 
 * Finally, this has been discussed in previous tutorial programs before:
 * 
 * @code
 * using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep5codeclasstemplate"></a> 
 * <h3>The <code>Step5</code> class template</h3>
 * 

 * 
 * The main class is mostly as in the previous example. The most visible
 * change is that the function <code>make_grid</code> has been
 * removed, since creating the grid is now done in the <code>run</code>
 * function and the rest of its functionality is now in
 * <code>setup_system</code>. Apart from this, everything is as before.
 * 
 * @code
 * template <int dim>
 * class Step5
 * {
 * public:
 *   Step5();
 *   void run();
 * 
 * private:
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void output_results(const unsigned int cycle) const;
 * 
 *   Triangulation<dim> triangulation;
 *   FE_Q<dim>          fe;
 *   DoFHandler<dim>    dof_handler;
 * 
 *   SparsityPattern      sparsity_pattern;
 *   SparseMatrix<double> system_matrix;
 * 
 *   Vector<double> solution;
 *   Vector<double> system_rhs;
 * };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Workingwithnonconstantcoefficients"></a> 
 * <h3>Working with nonconstant coefficients</h3>
 * 

 * 
 * In step-4, we showed how to use non-constant boundary values and right hand
 * side.  In this example, we want to use a variable coefficient in the
 * elliptic operator instead. Since we have a function which just depends on
 * the point in space we can do things a bit more simply and use a plain
 * function instead of inheriting from Function.
 * 

 * 
 * This is the implementation of the coefficient function for a single
 * point. We let it return 20 if the distance to the origin is less than 0.5,
 * and 1 otherwise.
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
 * @endcode
 * 
 * 
 * <a name="ThecodeStep5codeclassimplementation"></a> 
 * <h3>The <code>Step5</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Step5Step5"></a> 
 * <h4>Step5::Step5</h4>
 * 

 * 
 * This function is as before.
 * 
 * @code
 * template <int dim>
 * Step5<dim>::Step5()
 *   : fe(1)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5setup_system"></a> 
 * <h4>Step5::setup_system</h4>
 * 

 * 
 * This is the function <code>make_grid</code> from the previous
 * example, minus the generation of the grid. Everything else is unchanged:
 * 
 * @code
 * template <int dim>
 * void Step5<dim>::setup_system()
 * {
 *   dof_handler.distribute_dofs(fe);
 * 
 *   std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 * 
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp);
 *   sparsity_pattern.copy_from(dsp);
 * 
 *   system_matrix.reinit(sparsity_pattern);
 * 
 *   solution.reinit(dof_handler.n_dofs());
 *   system_rhs.reinit(dof_handler.n_dofs());
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5assemble_system"></a> 
 * <h4>Step5::assemble_system</h4>
 * 

 * 
 * As in the previous examples, this function is not changed much with regard
 * to its functionality, but there are still some optimizations which we will
 * show. For this, it is important to note that if efficient solvers are used
 * (such as the preconditioned CG method), assembling the matrix and right hand
 * side can take a comparable time, and you should think about using one or
 * two optimizations at some places.
 * 

 * 
 * The first parts of the function are completely unchanged from before:
 * 
 * @code
 * template <int dim>
 * void Step5<dim>::assemble_system()
 * {
 *   QGauss<dim> quadrature_formula(fe.degree + 1);
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
 * @endcode
 * 
 * Next is the typical loop over all cells to compute local contributions
 * and then to transfer them into the global matrix and vector. The only
 * change in this part, compared to step-4, is that we will use the
 * <code>coefficient()</code> function defined above to compute the
 * coefficient value at each quadrature point.
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       cell_matrix = 0.;
 *       cell_rhs    = 0.;
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
 *               cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                               1.0 *                               // f(x_q)
 *                               fe_values.JxW(q_index));            // dx
 *             }
 *         }
 * 
 * 
 *       cell->get_dof_indices(local_dof_indices);
 *       for (const unsigned int i : fe_values.dof_indices())
 *         {
 *           for (const unsigned int j : fe_values.dof_indices())
 *             system_matrix.add(local_dof_indices[i],
 *                               local_dof_indices[j],
 *                               cell_matrix(i, j));
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i);
 *         }
 *     }
 * 
 * @endcode
 * 
 * With the matrix so built, we use zero boundary values again:
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values;
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            Functions::ZeroFunction<dim>(),
 *                                            boundary_values);
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
 * <a name="Step5solve"></a> 
 * <h4>Step5::solve</h4>
 * 

 * 
 * The solution process again looks mostly like in the previous
 * examples. However, we will now use a preconditioned conjugate gradient
 * algorithm. It is not very difficult to make this change. In fact, the only
 * thing we have to alter is that we need an object which will act as a
 * preconditioner. We will use SSOR (symmetric successive overrelaxation),
 * with a relaxation factor of 1.2. For this purpose, the
 * <code>SparseMatrix</code> class has a function which does one SSOR step,
 * and we need to package the address of this function together with the
 * matrix on which it should act (which is the matrix to be inverted) and the
 * relaxation factor into one object. The <code>PreconditionSSOR</code> class
 * does this for us. (<code>PreconditionSSOR</code> class takes a template
 * argument denoting the matrix type it is supposed to work on. The default
 * value is <code>SparseMatrix@<double@></code>, which is exactly what we need
 * here, so we simply stick with the default and do not specify anything in
 * the angle brackets.)
 * 

 * 
 * Note that for the present case, SSOR doesn't really perform much better
 * than most other preconditioners (though better than no preconditioning at
 * all). A brief comparison of different preconditioners is presented in the
 * Results section of the next tutorial program, step-6.
 * 

 * 
 * With this, the rest of the function is trivial: instead of the
 * <code>PreconditionIdentity</code> object we have created before, we now use
 * the preconditioner we have declared, and the CG solver will do the rest for
 * us:
 * 
 * @code
 * template <int dim>
 * void Step5<dim>::solve()
 * {
 *   SolverControl            solver_control(1000, 1e-12);
 *   SolverCG<Vector<double>> solver(solver_control);
 * 
 *   PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *   preconditioner.initialize(system_matrix, 1.2);
 * 
 *   solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *   std::cout << "   " << solver_control.last_step()
 *             << " CG iterations needed to obtain convergence." << std::endl;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5output_resultsandsettingoutputflags"></a> 
 * <h4>Step5::output_results and setting output flags</h4>
 * 

 * 
 * Writing output to a file is mostly the same as for the previous tutorial.
 * The only difference is that we now need to construct a different filename
 * for each refinement cycle.
 * 

 * 
 * The function writes the output in VTU format, a variation of the VTK format
 * that requires less disk space because it compresses the data. Of course,
 * there are many other formats supported by the DataOut class if you
 * desire to use a program for visualization that doesn't understand
 * VTK or VTU.
 * 
 * @code
 * template <int dim>
 * void Step5<dim>::output_results(const unsigned int cycle) const
 * {
 *   DataOut<dim> data_out;
 * 
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 * 
 *   data_out.build_patches();
 * 
 *   std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
 *   data_out.write_vtu(output);
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5run"></a> 
 * <h4>Step5::run</h4>
 * 

 * 
 * The second to last thing in this program is the definition of the
 * <code>run()</code> function. In contrast to the previous programs, we will
 * compute on a sequence of meshes that after each iteration is globally
 * refined. The function therefore consists of a loop over 6 cycles. In each
 * cycle, we first print the cycle number, and then have to decide what to do
 * with the mesh. If this is not the first cycle, we simply refine the
 * existing mesh once globally. Before running through these cycles, however,
 * we have to generate a mesh:
 * 

 * 
 * In previous examples, we have already used some of the functions from the
 * <code>GridGenerator</code> class. Here we would like to read a grid from a
 * file where the cells are stored and which may originate from someone else,
 * or may be the product of a mesh generator tool.
 * 

 * 
 * In order to read a grid from a file, we generate an object of data type
 * GridIn and associate the triangulation to it (i.e. we tell it to fill our
 * triangulation object when we ask it to read the file). Then we open the
 * respective file and initialize the triangulation with the data in the file:
 * 
 * @code
 * template <int dim>
 * void Step5<dim>::run()
 * {
 *   GridIn<dim> grid_in;
 *   grid_in.attach_triangulation(triangulation);
 *   std::ifstream input_file("circle-grid.inp");
 * @endcode
 * 
 * We would now like to read the file. However, the input file is only for a
 * two-dimensional triangulation, while this function is a template for
 * arbitrary dimension. Since this is only a demonstration program, we will
 * not use different input files for the different dimensions, but rather
 * quickly kill the whole program if we are not in 2D. Of course, since the
 * main function below assumes that we are working in two dimensions we
 * could skip this check, in this version of the program, without any ill
 * effects.
 *   

 * 
 * It turns out that more than 90 per cent of programming errors are invalid
 * function parameters such as invalid array sizes, etc, so we use
 * assertions heavily throughout deal.II to catch such mistakes. For this,
 * the <code>Assert</code> macro is a good choice, since it makes sure that
 * the condition which is given as first argument is valid, and if not
 * throws an exception (its second argument) which will usually terminate
 * the program giving information where the error occurred and what the
 * reason was. (A longer discussion of what exactly the @p Assert macro
 * does can be found in the @ref Exceptions "exception documentation module".)
 * This generally reduces the time to find programming errors
 * dramatically and we have found assertions an invaluable means to program
 * fast.
 *   

 * 
 * On the other hand, all these checks (there are over 10,000 of them in the
 * library at present) should not slow down the program too much if you want
 * to do large computations. To this end, the <code>Assert</code> macro is
 * only used in debug mode and expands to nothing if in optimized
 * mode. Therefore, while you test your program on small problems and debug
 * it, the assertions will tell you where the problems are. Once your
 * program is stable, you can switch off debugging and the program will run
 * your real computations without the assertions and at maximum speed. More
 * precisely: turning off all the checks in the library (which prevent you
 * from calling functions with wrong arguments, walking off of arrays, etc.)
 * by compiling your program in optimized mode usually makes things run
 * about four times faster. Even though optimized programs are more
 * performant, we still recommend developing in debug mode since it allows
 * the library to find lots of common programming errors automatically. For
 * those who want to try: The way to switch from debug mode to optimized
 * mode is to recompile your program with the command <code>make
 * release</code>. The output of the <code>make</code> program should now
 * indicate to you that the program is now compiled in optimized mode, and
 * it will later also be linked to libraries that have been compiled for
 * optimized mode. In order to switch back to debug mode, simply recompile
 * with the command <code>make debug</code>.
 * 
 * @code
 *   Assert(dim == 2, ExcInternalError());
 * @endcode
 * 
 * ExcInternalError is a globally defined exception, which may be thrown
 * whenever something is terribly wrong. Usually, one would like to use more
 * specific exceptions, and particular in this case one would of course try
 * to do something else if <code>dim</code> is not equal to two, e.g. create
 * a grid using library functions. Aborting a program is usually not a good
 * idea and assertions should really only be used for exceptional cases
 * which should not occur, but might due to stupidity of the programmer,
 * user, or someone else. The situation above is not a very clever use of
 * Assert, but again: this is a tutorial and it might be worth to show what
 * not to do, after all.
 * 

 * 
 * So if we got past the assertion, we know that dim==2, and we can now
 * actually read the grid. It is in UCD (unstructured cell data) format
 * (though the convention is to use the suffix <code>inp</code> for UCD
 * files):
 * 
 * @code
 *   grid_in.read_ucd(input_file);
 * @endcode
 * 
 * If you like to use another input format, you have to use one of the other
 * <code>grid_in.read_xxx</code> function. (See the documentation of the
 * <code>GridIn</code> class to find out what input formats are presently
 * supported.)
 * 

 * 
 * The grid in the file describes a circle. Therefore we have to use a
 * manifold object which tells the triangulation where to put new points on
 * the boundary when the grid is refined. Unlike step-1, since GridIn does
 * not know that the domain has a circular boundary (unlike
 * GridGenerator::hyper_shell) we have to explicitly attach a manifold to
 * the boundary after creating the triangulation to get the correct result
 * when we refine the mesh.
 * 
 * @code
 *   const SphericalManifold<dim> boundary;
 *   triangulation.set_all_manifold_ids_on_boundary(0);
 *   triangulation.set_manifold(0, boundary);
 * 
 *   for (unsigned int cycle = 0; cycle < 6; ++cycle)
 *     {
 *       std::cout << "Cycle " << cycle << ':' << std::endl;
 * 
 *       if (cycle != 0)
 *         triangulation.refine_global(1);
 * 
 * @endcode
 * 
 * Now that we have a mesh for sure, we write some output and do all the
 * things that we have already seen in the previous examples.
 * 
 * @code
 *       std::cout << "   Number of active cells: "  
 *                 << triangulation.n_active_cells() 
 *                 << std::endl                      
 *                 << "   Total number of cells: "   
 *                 << triangulation.n_cells()        
 *                 << std::endl;
 * 
 *       setup_system();
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
 * The main function looks mostly like the one in the previous example, so we
 * won't comment on it further:
 * 
 * @code
 * int main()
 * {
 *   Step5<2> laplace_problem_2d;
 *   laplace_problem_2d.run();
 *   return 0;
 * }
 * @endcode
examples/step-5/doc/results.dox



<a name="Results"></a><h1>Results</h1>



下面是控制台的输出。

@code
Cycle 0:
   Number of active cells: 20
   Total number of cells: 20
   Number of degrees of freedom: 25
   13 CG iterations needed to obtain convergence.
Cycle 1:
   Number of active cells: 80
   Total number of cells: 100
   Number of degrees of freedom: 89
   18 CG iterations needed to obtain convergence.
Cycle 2:
   Number of active cells: 320
   Total number of cells: 420
   Number of degrees of freedom: 337
   29 CG iterations needed to obtain convergence.
Cycle 3:
   Number of active cells: 1280
   Total number of cells: 1700
   Number of degrees of freedom: 1313
   52 CG iterations needed to obtain convergence.
Cycle 4:
   Number of active cells: 5120
   Total number of cells: 6820
   Number of degrees of freedom: 5185
   95 CG iterations needed to obtain convergence.
Cycle 5:
   Number of active cells: 20480
   Total number of cells: 27300
   Number of degrees of freedom: 20609
   182 CG iterations needed to obtain convergence.
@endcode






在每个周期中，单元格的数量翻了两番，CG迭代的数量大约翻了一番。另外，在每个周期中，程序写出一个VTU格式的输出图形文件。它们被描述为以下内容。

 <table width="100%">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-0-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-1-r9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-2-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-3-r9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-4-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-5-r9.2.png" alt="">
    </td>
  </tr>
</table> 




由于系数的可变性（那里的曲率减少的程度与系数增加的程度相同），溶液的顶部区域被压扁了。溶液的梯度沿着界面是不连续的，尽管这在上面的图片中不是很明显。我们将在下一个例子中更详细地研究这个问题。

图片还显示，这个程序计算出来的解在非常粗的网格上其实是相当错误的（它的大小是错误的）。这是因为没有任何数值方法能保证粗大网格上的解是特别准确的--但我们知道解<i>converges</i>是精确的解，事实上你可以看到从一个网格到下一个网格的解似乎在最后不再有太大的变化。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-5.cc"
*/
