/**
@page step_10 The step-10 tutorial program
This tutorial depends on step-7.

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
examples/step-10/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这是一个相当短的例子，只显示了使用高阶映射的一些方面。我们所说的 <em> 映射 </em> 是指单元格（即单位线、正方形或立方体）与现实空间中的单元格之间的转换。在前面所有的例子中，我们都隐含地使用了线性或d-线性映射；你根本不会注意到这一点，因为如果你不做任何特别的事情，这就是发生的情况。然而，如果你的域有弯曲的边界，在有些情况下，边界的片状线性逼近（即通过直线段）是不够的，你希望你的计算域也是一个使用弯曲边界的真实域的逼近。如果边界近似使用分片的二次抛物线来近似真实边界，那么我们说这是二次或 $Q_2$ 的近似。如果我们使用成片的立体多项式的图形，那么这是一个 $Q_3$ 的近似，以此类推。




对于某些微分方程，已知如果精确域的边界是弯曲的，那么边界的片状线性近似，即 $Q_1$ 映射，是不够的。例如，使用 $C^1$ 元素的偏谐方程，或具有弯曲反射边界的域上的气体动力学的欧拉方程。在这些情况下，有必要使用高阶映射来计算积分。如果我们不使用这样的高阶映射，那么边界的逼近顺序将主导整个数值方案的收敛顺序，而不考虑域的内部离散化的收敛顺序。




我们没有用这些更复杂的例子来证明高阶映射的使用，而是只做了一个简单的计算：用两种不同的方法计算 $\pi=3.141592653589793238462643\ldots$ 的值。




第一种方法使用单位半径的圆的三角形近似，并在其上积分一个单位幅度的常数函数（ $f = 1$ ）。当然，如果领域是精确的单位圆，那么面积将是 $\pi$ ，但由于我们只使用片状多项式段的近似，我们积分的面积值并不完全是 $\pi$  。然而，众所周知，当我们细化三角形时， $Q_p$ 映射以 $h^{p+1}$ 的阶数逼近边界，其中 $h$ 是网格大小。我们将检查圆的计算面积值，以及它们在不同映射的网格细化下向 $\pi$ 收敛的情况。我们还将发现一种收敛行为，这种行为一开始是令人惊讶的，但有一个很好的解释。




第二种方法与此类似，但这次不使用三角形单位圆的面积，而是使用其周长。   $\pi$ 然后用周长的一半来近似，因为我们选择半径等于1。




 @note  本教程实质上展示了如何为积分选择一个特定的映射，方法是将一个特定的几何体附加到三角形上（例如在步骤1中已经完成），然后将一个映射参数传递给FEValues类，该类用于deal.II中的所有积分。我们选择的几何体是一个圆，deal.II已经有一个可以使用的类（SphericalManifold）。如果你想定义你自己的几何体，例如，因为它很复杂，不能用deal.II中已有的类来描述，你会想通过步骤53来阅读。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * The first of the following include files are probably well-known by now and
 * need no further explanation.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/manifold_lib.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * @endcode
 * 
 * This include file is new. Even if we are not solving a PDE in this tutorial,
 * we want to use a dummy finite element with zero degrees of freedoms provided
 * by the FE_Nothing class.
 * 
 * @code
 * #include <deal.II/fe/fe_nothing.h>
 * 
 * @endcode
 * 
 * The following header file is also new: in it, we declare the MappingQ class
 * which we will use for polynomial mappings of arbitrary order:
 * 
 * @code
 * #include <deal.II/fe/mapping_q.h>
 * 
 * @endcode
 * 
 * And this again is C++:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * #include <cmath>
 * 
 * @endcode
 * 
 * The last step is as in previous programs:
 * 
 * @code
 * namespace Step10
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * Now, as we want to compute the value of $\pi$, we have to compare to
 * something. These are the first few digits of $\pi$, which we define
 * beforehand for later use. Since we would like to compute the difference
 * between two numbers which are quite accurate, with the accuracy of the
 * computed approximation to $\pi$ being in the range of the number of
 * digits which a double variable can hold, we rather declare the reference
 * value as a <code>long double</code> and give it a number of extra digits:
 * 
 * @code
 *   const long double pi = 3.141592653589793238462643L;
 * 
 * 
 * 
 * @endcode
 * 
 * Then, the first task will be to generate some output. Since this program
 * is so small, we do not employ object oriented techniques in it and do not
 * declare classes (although, of course, we use the object oriented features
 * of the library). Rather, we just pack the functionality into separate
 * functions. We make these functions templates on the number of space
 * dimensions to conform to usual practice when using deal.II, although we
 * will only use them for two space dimensions and throw an exception when
 * attempted to use for any other spatial dimension.
 *   

 * 
 * The first of these functions just generates a triangulation of a circle
 * (hyperball) and outputs the $Q_p$ mapping of its cells for different values
 * of <code>p</code>. Then, we refine the grid once and do so again.
 * 
 * @code
 *   template <int dim>
 *   void gnuplot_output()
 *   {
 *     std::cout << "Output of grids into gnuplot files:" << std::endl
 *               << "===================================" << std::endl;
 * 
 * @endcode
 * 
 * So first generate a coarse triangulation of the circle and associate a
 * suitable boundary description to it. By default,
 * GridGenerator::hyper_ball attaches a SphericalManifold to the boundary
 * (and uses FlatManifold for the interior) so we simply call that
 * function and move on:
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     GridGenerator::hyper_ball(triangulation);
 * 
 * @endcode
 * 
 * Then alternate between generating output on the current mesh
 * for $Q_1$, $Q_2$, and $Q_3$ mappings, and (at the end of the
 * loop body) refining the mesh once globally.
 * 
 * @code
 *     for (unsigned int refinement = 0; refinement < 2; ++refinement)
 *       {
 *         std::cout << "Refinement level: " << refinement << std::endl;
 * 
 *         std::string filename_base = "ball_" + std::to_string(refinement);
 * 
 *         for (unsigned int degree = 1; degree < 4; ++degree)
 *           {
 *             std::cout << "Degree = " << degree << std::endl;
 * 
 * @endcode
 * 
 * For this, first set up an object describing the mapping. This
 * is done using the MappingQ class, which takes as
 * argument to the constructor the polynomial degree which it
 * shall use.
 * 
 * @code
 *             const MappingQ<dim> mapping(degree);
 * @endcode
 * 
 * As a side note, for a piecewise linear mapping, you
 * could give a value of <code>1</code> to the constructor
 * of MappingQ, but there is also a class MappingQ1 that
 * achieves the same effect. Historically, it did a lot of
 * things in a simpler way than MappingQ but is today just
 * a wrapper around the latter. It is, however, still the
 * class that is used implicitly in many places of the
 * library if you do not specify another mapping
 * explicitly.
 * 

 * 
 * 

 * 
 * In order to actually write out the present grid with this
 * mapping, we set up an object which we will use for output. We
 * will generate Gnuplot output, which consists of a set of lines
 * describing the mapped triangulation. By default, only one line
 * is drawn for each face of the triangulation, but since we want
 * to explicitly see the effect of the mapping, we want to have
 * the faces in more detail. This can be done by passing the
 * output object a structure which contains some flags. In the
 * present case, since Gnuplot can only draw straight lines, we
 * output a number of additional points on the faces so that each
 * face is drawn by 30 small lines instead of only one. This is
 * sufficient to give us the impression of seeing a curved line,
 * rather than a set of straight lines.
 * 
 * @code
 *             GridOut               grid_out;
 *             GridOutFlags::Gnuplot gnuplot_flags(false, 60);
 *             grid_out.set_flags(gnuplot_flags);
 * 
 * @endcode
 * 
 * Finally, generate a filename and a file for output:
 * 
 * @code
 *             std::string filename =
 *               filename_base + "_mapping_q_" + std::to_string(degree) + ".dat";
 *             std::ofstream gnuplot_file(filename);
 * 
 * @endcode
 * 
 * Then write out the triangulation to this file. The last
 * argument of the function is a pointer to a mapping object. This
 * argument has a default value, and if no value is given a simple
 * MappingQ1 object is taken, which we briefly
 * described above. This would then result in a piecewise linear
 * approximation of the true boundary in the output.
 * 
 * @code
 *             grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping);
 *           }
 *         std::cout << std::endl;
 * 
 * @endcode
 * 
 * At the end of the loop, refine the mesh globally.
 * 
 * @code
 *         triangulation.refine_global();
 *       }
 *   }
 * 
 * @endcode
 * 
 * Now we proceed with the main part of the code, the approximation of
 * $\pi$. The area of a circle is of course given by $\pi r^2$, so having a
 * circle of radius 1, the area represents just the number that is searched
 * for. The numerical computation of the area is performed by integrating
 * the constant function of value 1 over the whole computational domain,
 * i.e. by computing the areas $\int_K 1 dx=\int_{\hat K} 1
 * \ \textrm{det}\ J(\hat x) d\hat x \approx \sum_i \textrm{det}
 * \ J(\hat x_i)w(\hat x_i)$,
 * where the sum extends over all quadrature points on all active cells in
 * the triangulation, with $w(x_i)$ being the weight of quadrature point
 * $x_i$. The integrals on each cell are approximated by numerical
 * quadrature, hence the only additional ingredient we need is to set up a
 * FEValues object that provides the corresponding `JxW` values of each
 * cell. (Note that `JxW` is meant to abbreviate <i>Jacobian determinant
 * times weight</i>; since in numerical quadrature the two factors always
 * occur at the same places, we only offer the combined quantity, rather
 * than two separate ones.) We note that here we won't use the FEValues
 * object in its original purpose, i.e. for the computation of values of
 * basis functions of a specific finite element at certain quadrature
 * points. Rather, we use it only to gain the `JxW` at the quadrature
 * points, irrespective of the (dummy) finite element we will give to the
 * constructor of the FEValues object. The actual finite element given to
 * the FEValues object is not used at all, so we could give any.
 * 
 * @code
 *   template <int dim>
 *   void compute_pi_by_area()
 *   {
 *     std::cout << "Computation of Pi by the area:" << std::endl
 *               << "==============================" << std::endl;
 * 
 * @endcode
 * 
 * For the numerical quadrature on all cells we employ a quadrature rule
 * of sufficiently high degree. We choose QGauss that is of order 8 (4
 * points), to be sure that the errors due to numerical quadrature are of
 * higher order than the order (maximal 6) that will occur due to the
 * order of the approximation of the boundary, i.e. the order of the
 * mappings employed. Note that the integrand, the Jacobian determinant,
 * is not a polynomial function (rather, it is a rational one), so we do
 * not use Gauss quadrature in order to get the exact value of the
 * integral as done often in finite element computations, but could as
 * well have used any quadrature formula of like order instead.
 * 
 * @code
 *     const QGauss<dim> quadrature(4);
 * 
 * @endcode
 * 
 * Now start by looping over polynomial mapping degrees=1..4:
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree)
 *       {
 *         std::cout << "Degree = " << degree << std::endl;
 * 
 * @endcode
 * 
 * First generate the triangulation, the boundary and the mapping
 * object as already seen.
 * 
 * @code
 *         Triangulation<dim> triangulation;
 *         GridGenerator::hyper_ball(triangulation);
 * 
 *         const MappingQ<dim> mapping(degree);
 * 
 * @endcode
 * 
 * We now create a finite element. Unlike the rest of the example
 * programs, we do not actually need to do any computations with shape
 * functions; we only need the `JxW` values from an FEValues
 * object. Hence we use the special finite element class FE_Nothing
 * which has exactly zero degrees of freedom per cell (as the name
 * implies, the local basis on each cell is the empty set). A more
 * typical usage of FE_Nothing is shown in step-46.
 * 
 * @code
 *         const FE_Nothing<dim> fe;
 * 
 * @endcode
 * 
 * Likewise, we need to create a DoFHandler object. We do not actually
 * use it, but it will provide us with `active_cell_iterators` that
 * are needed to reinitialize the FEValues object on each cell of the
 * triangulation.
 * 
 * @code
 *         DoFHandler<dim> dof_handler(triangulation);
 * 
 * @endcode
 * 
 * Now we set up the FEValues object, giving the Mapping, the dummy
 * finite element and the quadrature object to the constructor,
 * together with the update flags asking for the `JxW` values at the
 * quadrature points only. This tells the FEValues object that it
 * needs not compute other quantities upon calling the
 * <code>reinit</code> function, thus saving computation time.
 *         

 * 
 * The most important difference in the construction of the FEValues
 * object compared to previous example programs is that we pass a
 * mapping object as first argument, which is to be used in the
 * computation of the mapping from unit to real cell. In previous
 * examples, this argument was omitted, resulting in the implicit use
 * of an object of type MappingQ1.
 * 
 * @code
 *         FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values);
 * 
 * @endcode
 * 
 * We employ an object of the ConvergenceTable class to store all
 * important data like the approximated values for $\pi$ and the error
 * with respect to the true value of $\pi$. We will also use functions
 * provided by the ConvergenceTable class to compute convergence rates
 * of the approximations to $\pi$.
 * 
 * @code
 *         ConvergenceTable table;
 * 
 * @endcode
 * 
 * Now we loop over several refinement steps of the triangulation.
 * 
 * @code
 *         for (unsigned int refinement = 0; refinement < 6;
 *              ++refinement, triangulation.refine_global(1))
 *           {
 * @endcode
 * 
 * In this loop we first add the number of active cells of the
 * current triangulation to the table. This function automatically
 * creates a table column with superscription `cells`, in case
 * this column was not created before.
 * 
 * @code
 *             table.add_value("cells", triangulation.n_active_cells());
 * 
 * @endcode
 * 
 * Then we distribute the degrees of freedom for the dummy finite
 * element. Strictly speaking we do not need this function call in
 * our special case but we call it to make the DoFHandler happy --
 * otherwise it would throw an assertion in the FEValues::reinit
 * function below.
 * 
 * @code
 *             dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * We define the variable area as `long double` like we did for
 * the `pi` variable before.
 * 
 * @code
 *             long double area = 0;
 * 
 * @endcode
 * 
 * Now we loop over all cells, reinitialize the FEValues object
 * for each cell, and add up all the `JxW` values for this cell to
 * `area`...
 * 
 * @code
 *             for (const auto &cell : dof_handler.active_cell_iterators())
 *               {
 *                 fe_values.reinit(cell);
 *                 for (unsigned int i = 0; i < fe_values.n_quadrature_points; ++i)
 *                   area += static_cast<long double>(fe_values.JxW(i));
 *               }
 * 
 * @endcode
 * 
 * ...and store the resulting area values and the errors in the
 * table. We need a static cast to double as there is no
 * add_value(string, long double) function implemented. Note that
 * this also concerns the second call as the <code>fabs</code>
 * function in the <code>std</code> namespace is overloaded on its
 * argument types, so there exists a version taking and returning
 * a <code>long double</code>, in contrast to the global namespace
 * where only one such function is declared (which takes and
 * returns a double).
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(area));
 *             table.add_value("error", static_cast<double>(std::fabs(area - pi)));
 *           }
 * 
 * @endcode
 * 
 * We want to compute the convergence rates of the `error`
 * column. Therefore we need to omit the other columns from the
 * convergence rate evaluation before calling
 * `evaluate_all_convergence_rates`
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells");
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi");
 *         table.evaluate_all_convergence_rates(
 *           ConvergenceTable::reduction_rate_log2);
 * 
 * @endcode
 * 
 * Finally we set the precision and scientific mode for output of some
 * of the quantities...
 * 
 * @code
 *         table.set_precision("eval.pi", 16);
 *         table.set_scientific("error", true);
 * 
 * @endcode
 * 
 * ...and write the whole table to std::cout.
 * 
 * @code
 *         table.write_text(std::cout);
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * The following, second function also computes an approximation of $\pi$
 * but this time via the perimeter $2\pi r$ of the domain instead of the
 * area. This function is only a variation of the previous function. So we
 * will mainly give documentation for the differences.
 * 
 * @code
 *   template <int dim>
 *   void compute_pi_by_perimeter()
 *   {
 *     std::cout << "Computation of Pi by the perimeter:" << std::endl
 *               << "===================================" << std::endl;
 * 
 * @endcode
 * 
 * We take the same order of quadrature but this time a `dim-1`
 * dimensional quadrature as we will integrate over (boundary) lines
 * rather than over cells.
 * 
 * @code
 *     const QGauss<dim - 1> quadrature(4);
 * 
 * @endcode
 * 
 * We loop over all degrees, create the triangulation, the boundary, the
 * mapping, the dummy finite element and the DoFHandler object as seen
 * before.
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree)
 *       {
 *         std::cout << "Degree = " << degree << std::endl;
 *         Triangulation<dim> triangulation;
 *         GridGenerator::hyper_ball(triangulation);
 * 
 *         const MappingQ<dim>   mapping(degree);
 *         const FE_Nothing<dim> fe;
 * 
 *         DoFHandler<dim> dof_handler(triangulation);
 * 
 * @endcode
 * 
 * Then we create a FEFaceValues object instead of a FEValues object
 * as in the previous function. Again, we pass a mapping as first
 * argument.
 * 
 * @code
 *         FEFaceValues<dim> fe_face_values(mapping,
 *                                          fe,
 *                                          quadrature,
 *                                          update_JxW_values);
 *         ConvergenceTable  table;
 * 
 *         for (unsigned int refinement = 0; refinement < 6;
 *              ++refinement, triangulation.refine_global(1))
 *           {
 *             table.add_value("cells", triangulation.n_active_cells());
 * 
 *             dof_handler.distribute_dofs(fe);
 * 
 * @endcode
 * 
 * Now we run over all cells and over all faces of each cell. Only
 * the contributions of the `JxW` values on boundary faces are
 * added to the long double variable `perimeter`.
 * 
 * @code
 *             long double perimeter = 0;
 *             for (const auto &cell : dof_handler.active_cell_iterators())
 *               for (const auto &face : cell->face_iterators())
 *                 if (face->at_boundary())
 *                   {
 * @endcode
 * 
 * We reinit the FEFaceValues object with the cell
 * iterator and the number of the face.
 * 
 * @code
 *                     fe_face_values.reinit(cell, face);
 *                     for (unsigned int i = 0;
 *                          i < fe_face_values.n_quadrature_points;
 *                          ++i)
 *                       perimeter +=
 *                         static_cast<long double>(fe_face_values.JxW(i));
 *                   }
 * @endcode
 * 
 * Then store the evaluated values in the table...
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(perimeter / 2.0L));
 *             table.add_value(
 *               "error", static_cast<double>(std::fabs(perimeter / 2.0L - pi)));
 *           }
 * 
 * @endcode
 * 
 * ...and end this function as we did in the previous one:
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells");
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi");
 *         table.evaluate_all_convergence_rates(
 *           ConvergenceTable::reduction_rate_log2);
 * 
 *         table.set_precision("eval.pi", 16);
 *         table.set_scientific("error", true);
 * 
 *         table.write_text(std::cout);
 * 
 *         std::cout << std::endl;
 *       }
 *   }
 * } // namespace Step10
 * 
 * 
 * @endcode
 * 
 * The following main function just calls the above functions in the order of
 * their appearance. Apart from this, it looks just like the main functions of
 * previous tutorial programs.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       std::cout.precision(16);
 * 
 *       const unsigned int dim = 2;
 * 
 *       Step10::gnuplot_output<dim>();
 * 
 *       Step10::compute_pi_by_area<dim>();
 *       Step10::compute_pi_by_perimeter<dim>();
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
examples/step-10/doc/results.dox



<a name="Results"></a><h1>Results</h1>



该程序执行两个任务，第一个是生成映射域的可视化，第二个是通过所述的两种方法计算π。让我们先看一下生成的图形。它们是以Gnuplot格式生成的，可以用以下命令查看

@code
set style data lines
set size ratio -1
unset key
unset xtics
unset ytics
plot [-1:1][-1:1] "ball_0_mapping_q_1.dat" lw 4 lt rgb "black"
@endcode

或使用其他文件名之一。第二行确保生成的输出的长宽比实际上是1:1，也就是说，一个圆在你的屏幕上被画成一个圆，而不是一个椭圆。第三行关闭了图形中的键，因为那只会打印信息（文件名），而这些信息现在并不那么重要。同样地，第四行和第五行关闭了刻度线。然后生成具有特定线宽（"lw"，这里设置为4）和线型（"lt"，这里选择线应使用RGB颜色 "黑色"）的图。

下表显示了 $Q_1$ 、 $Q_2$ 和 $Q_3$ 映射的三角计算域，为原始粗网格（左）和一次均匀细化网格（右）。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q1.svg" alt="磁盘的五单元离散化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q1.svg" alt="磁盘的20单元离散化（即。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q2.svg" alt="五格离散化的圆盘，边缘为二次方。边界与实际的圆几乎没有区别。" width="400" height="400" > </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q2.svg" alt="带有二次方边缘的20个单元离散化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q3.svg" alt="带有三次方边缘的圆的五单元离散化。边界与实际的圆几乎没有区别。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q3.svg" alt="具有立方体边缘的20单元离散化。" width="400" height="400"> </div> </div>

这些图片显示了高阶映射的明显优势：它们在相当粗的网格上也能相当好地接近真实边界。为了进一步证明这一点，这里是使用 $Q_2$ 和 $Q_3$ 映射的粗网格的右上角四分之一圈的一部分，其中红色虚线标志着实际的圆。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q2.svg" alt="二次离散化的特写。二次插值和实际圆之间的距离很小。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q3.svg" alt="立方离散化的特写。立体插值和实际圆之间的距离非常小。" width="400" height="400"> </div> </div> </div>

很明显，二次映射很好地接近了边界，而对于三次映射来说，对于粗略的网格来说，近似域和真实域之间的差异已经很难看到了。你还可以看到，映射只在三角形的外部边界有一些变化。在内部，所有的线条仍然是由线性函数表示的，这就导致了只在边界的单元上进行额外的计算。因此，高阶映射通常不会比低阶映射明显地慢，因为额外的计算只在所有单元格的一个小子集上执行。




该程序的第二个目的是计算π的值，以达到良好的准确性。这是程序的这一部分的输出。

@code
Output of grids into gnuplot files:
===================================
Refinement level: 0
Degree = 1
Degree = 2
Degree = 3


Refinement level: 1
Degree = 1
Degree = 2
Degree = 3


Computation of Pi by the area:
==============================
Degree = 1
cells      eval.pi            error
    5 1.9999999999999993 1.1416e+00    -
   20 2.8284271247461890 3.1317e-01 1.87
   80 3.0614674589207174 8.0125e-02 1.97
  320 3.1214451522580511 2.0148e-02 1.99
 1280 3.1365484905459380 5.0442e-03 2.00
 5120 3.1403311569547516 1.2615e-03 2.00


Degree = 2
cells      eval.pi            error
    5 3.1045694996615860 3.7023e-02    -
   20 3.1391475703122267 2.4451e-03 3.92
   80 3.1414377167038290 1.5494e-04 3.98
  320 3.1415829366419006 9.7169e-06 4.00
 1280 3.1415920457576898 6.0783e-07 4.00
 5120 3.1415926155921117 3.7998e-08 4.00


Degree = 3
cells      eval.pi            error
    5 3.1410033851499288 5.8927e-04    -
   20 3.1415830393583839 9.6142e-06 5.94
   80 3.1415925017363797 1.5185e-07 5.98
  320 3.1415926512106696 2.3791e-09 6.00
 1280 3.1415926535525927 3.7200e-11 6.00
 5120 3.1415926535892100 5.8302e-13 6.00


Degree = 4
cells      eval.pi            error
    5 3.1415871927401131 5.4608e-06    -
   20 3.1415926314742428 2.2116e-08 7.95
   80 3.1415926535026202 8.7173e-11 7.99
  320 3.1415926535894498 3.4350e-13 7.99
 1280 3.1415926535897896 3.4671e-15 6.63
 5120 3.1415926535897909 2.4009e-15 0.53


Computation of Pi by the perimeter:
===================================
Degree = 1
cells      eval.pi            error
    5 2.8284271247461898 3.1317e-01    -
   20 3.0614674589207178 8.0125e-02 1.97
   80 3.1214451522580520 2.0148e-02 1.99
  320 3.1365484905459389 5.0442e-03 2.00
 1280 3.1403311569547525 1.2615e-03 2.00
 5120 3.1412772509327724 3.1540e-04 2.00


Degree = 2
cells      eval.pi            error
    5 3.1248930668550594 1.6700e-02    -
   20 3.1404050605605449 1.1876e-03 3.81
   80 3.1415157631807009 7.6890e-05 3.95
  320 3.1415878042798613 4.8493e-06 3.99
 1280 3.1415923498174534 3.0377e-07 4.00
 5120 3.1415926345931995 1.8997e-08 4.00


Degree = 3
cells      eval.pi            error
    5 3.1414940401456048 9.8613e-05    -
   20 3.1415913432549156 1.3103e-06 6.23
   80 3.1415926341726910 1.9417e-08 6.08
  320 3.1415926532906897 2.9910e-10 6.02
 1280 3.1415926535851355 4.6578e-12 6.00
 5120 3.1415926535897190 7.4216e-14 5.97


Degree = 4
cells      eval.pi            error
    5 3.1415921029432572 5.5065e-07     -
   20 3.1415926513737595 2.2160e-09  7.96
   80 3.1415926535810712 8.7222e-12  7.99
  320 3.1415926535897576 3.5525e-14  7.94
 1280 3.1415926535897936 4.6729e-16  6.25
 5120 3.1415926535897918 1.4929e-15 -1.68
@endcode






从输出的一个直接观察结果是，在所有情况下，数值都迅速收敛到 $\pi=3.141592653589793238462643$ 的真实值。请注意，对于 $Q_4$ 的映射，我们已经进入了四舍五入误差的制度，收敛率趋于平缓，这已经是相当多的了。然而，也请注意，对于 $Q_1$ 映射，即使在最细的网格上，精度也明显比 $Q_3$ 映射的粗网格上要差得多!




输出的最后一列显示了收敛顺序，以网格宽度的幂为单位  $h$  。在介绍中，我们曾说过  $Q_p$  映射的收敛顺序应该是  $h^{p+1}$  。然而，在所示的例子中，顺序是 $h^{2p}$  !这个起初令人惊讶的事实可以用 $Q_p$ 映射的特性来解释。在<i>p</i>阶时，它使用的支持点是基于<i>p</i>+1点的Gauss-Lobatto正交规则，以这样的方式选择支持点，使正交规则在2<i>p</i>阶时收敛。尽管这些点在这里只用于插值<i>p</i>阶多项式，但我们在数值评估积分时得到了超收敛效应，导致观察到的高阶收敛。这种效应在下面的出版物中也有详细讨论。A. Bonito, A. Demlow, and J. Owen:"拉普拉斯-贝特拉米算子的特征值和特征函数的有限元近似的先验误差估计"，提交，2018年）。)


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-10.cc"
*/
