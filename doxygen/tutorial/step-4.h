/**
@page step_4 The step-4 tutorial program
This tutorial depends on step-3.

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
        <li><a href="#ThecodeStep4codeclasstemplate">The <code>Step4</code> class template</a>
        <li><a href="#Righthandsideandboundaryvalues">Right hand side and boundary values</a>
        <li><a href="#ImplementationofthecodeStep4codeclass">Implementation of the <code>Step4</code> class</a>
      <ul>
        <li><a href="#Step4Step4">Step4::Step4</a>
        <li><a href="#Step4make_grid">Step4::make_grid</a>
        <li><a href="#Step4setup_system">Step4::setup_system</a>
        <li><a href="#Step4assemble_system">Step4::assemble_system</a>
        <li><a href="#Step4solve">Step4::solve</a>
        <li><a href="#Step4output_results">Step4::output_results</a>
        <li><a href="#Step4run">Step4::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-4/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{12,13} 

deal.II有一个独特的功能，我们称之为 "无维度编程"。你可能已经注意到，在前面的例子中，许多类的后缀都是角括号中的数字。这是为了表明，例如，二维和三维空间的三角形是不同的，但是相关的数据%类型。我们完全可以把它们称为 <code>Triangulation2d</code> and <code>Triangulation3d</code> 而不是 <code>Triangulation@<2@></code> 和 <code>Triangulation@<3@></code> 来命名这两个类，但这有一个重要的缺点：假设你有一个功能完全相同的函数，但在2D或3D三角形上，取决于我们目前想在哪个维度上解方程（如果你不相信一个函数在所有维度上都做同样的事情是常见的情况，看看下面的代码就知道了，2D和3D之间几乎没有区别！）。我们将不得不把同一个函数写两次，一次在 <code>Triangulation2d</code> 上工作，一次在 <code>Triangulation3d</code> 上工作。这在编程中是一个不必要的障碍，并且导致了保持两个函数同步的麻烦（最好是），或者在两个版本不同步时难以发现错误（最坏的情况是；这可能是更常见的情况）。





这种障碍可以通过使用C++语言提供的一些模板魔法来规避：模板化的类和函数并不是真正的类或函数，而只是取决于一个尚未定义的数据类型参数或在定义时也未知的数值的一种模式。然而，如果你向它提供了所需的信息，编译器可以从这些模板中建立适当的类或函数。当然，模板的部分内容可以依赖于模板参数，它们将在编译时被解析为特定的模板参数。例如，考虑下面这段代码。

@code
  template <int dim>
  void make_grid (Triangulation<dim> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode






在编译器看到这个函数的时候，它对 <code>dim</code> 的实际值并不了解。编译器唯一拥有的是一个模板，即蓝图，如果给定 <code>make_grid</code> 的特定值有一个未知的值，编译器暂时没有可以生成的代码。




然而，如果以后下来，编译器会遇到一些代码，例如，看起来像这样。

@code
  Triangulation<2> triangulation;
  make_grid (triangulation);
@endcode

那么编译器将推断出请求将函数 <code>make_grid</code> 替换为 <code>dim==2</code> ，并将上述模板编译为一个到处都用2替换了dim的函数，也就是说，它将编译该函数，就好像它被定义为

@code
  void make_grid (Triangulation<2> &triangulation)
  {
    GridGenerator::hyper_cube (triangulation, -1, 1);
  };
@endcode






然而，值得注意的是，函数 <code>GridGenerator::hyper_cube</code> 也取决于维度，所以在这种情况下，编译器将调用函数 <code>GridGenerator::hyper_cube@<2@></code> ，而如果dim是3，它将调用 <code>GridGenerator::hyper_cube@<3@></code> ，这可能是（实际上是）一个完全无关的函数。




对成员变量也可以这样做。考虑一下下面的函数，它可能反过来调用上面的函数。

@code
  template <int dim>
  void make_grid_and_dofs (Triangulation<dim> &triangulation)
  {
    make_grid (triangulation);


    DoFHandler<dim> dof_handler(triangulation);
    ...
  };
@endcode

这个函数有一个类型为  <code>DoFHandler@<dim@></code>  的成员变量。同样，编译器在知道哪个维度之前不能编译这个函数。如果你像上面那样为一个特定的维度调用这个函数，编译器将使用模板，用调用的维度替换所有出现的dim，并编译它。如果你为不同的维度多次调用该函数，它将多次编译，每次都调用正确的 <code>make_grid</code> 函数，并为成员变量保留适当的内存量；注意， <code>DoFHandler</code> 的大小可能，事实上也确实取决于空间维度。




deal.II库是围绕这个独立于维度的编程概念建立的，因此允许你以一种不需要区分空间维度的方式来编程。应该注意的是，只有在极少数的地方才有必要使用 <code>if</code>s or <code>switch</code> es来实际比较尺寸。然而，由于编译器必须为每个维度单独编译每个函数，即使在那里，它在编译时也知道 <code>dim</code> 的值，因此将能够优化掉 <code>if</code> 语句和未使用的分支。




在这个例子程序中，我们将展示如何独立编程维度（事实上，这比你必须照顾到维度还要简单），我们将把上一个例子的拉普拉斯问题扩展到一个同时在两个和三个空间维度运行的程序。其他的扩展是使用非恒定的右手边函数和非零边界值。




 @note  在使用模板时，C++强加了各种语法限制，有时让人有点难以理解为什么到底要这样写。一个典型的例子是，在很多地方都需要使用关键字 <code>typename</code> 。如果你已经不完全熟悉，那么在<a
href="http://www.dealii.org/">deal.II homepage</a>中链接的deal.II常见问题解答（FAQ）中解释了其中的几个困难。

<！--我们需要一个空行来正确结束上述块。


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
 * The first few (many?) include files have already been used in the previous
 * example, so we will not explain their meaning here again.
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * This is new, however: in the previous example we got some unwanted output
 * from the linear solvers. If we want to suppress it, we have to include this
 * file and add a single line somewhere to the program (see the main()
 * function below for that):
 * 
 * @code
 * #include <deal.II/base/logstream.h>
 * 
 * @endcode
 * 
 * The final step, as in previous programs, is to import all the deal.II class
 * and function names into the global namespace:
 * 
 * @code
 * using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep4codeclasstemplate"></a> 
 * <h3>The <code>Step4</code> class template</h3>
 * 

 * 
 * This is again the same <code>Step4</code> class as in the previous
 * example. The only difference is that we have now declared it as a class
 * with a template parameter, and the template parameter is of course the
 * spatial dimension in which we would like to solve the Laplace equation. Of
 * course, several of the member variables depend on this dimension as well,
 * in particular the Triangulation class, which has to represent
 * quadrilaterals or hexahedra, respectively. Apart from this, everything is
 * as before.
 * 
 * @code
 * template <int dim>
 * class Step4
 * {
 * public:
 *   Step4();
 *   void run();
 * 
 * private:
 *   void make_grid();
 *   void setup_system();
 *   void assemble_system();
 *   void solve();
 *   void output_results() const;
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
 * <a name="Righthandsideandboundaryvalues"></a> 
 * <h3>Right hand side and boundary values</h3>
 * 

 * 
 * In the following, we declare two more classes denoting the right hand side
 * and the non-homogeneous Dirichlet boundary values. Both are functions of a
 * dim-dimensional space variable, so we declare them as templates as well.
 * 

 * 
 * Each of these classes is derived from a common, abstract base class
 * Function, which declares the common interface which all functions have to
 * follow. In particular, concrete classes have to overload the
 * <code>value</code> function, which takes a point in dim-dimensional space
 * as parameters and returns the value at that point as a
 * <code>double</code> variable.
 * 

 * 
 * The <code>value</code> function takes a second argument, which we have here
 * named <code>component</code>: This is only meant for vector-valued
 * functions, where you may want to access a certain component of the vector
 * at the point <code>p</code>. However, our functions are scalar, so we need
 * not worry about this parameter and we will not use it in the implementation
 * of the functions. Inside the library's header files, the Function base
 * class's declaration of the <code>value</code> function has a default value
 * of zero for the component, so we will access the <code>value</code>
 * function of the right hand side with only one parameter, namely the point
 * where we want to evaluate the function. A value for the component can then
 * simply be omitted for scalar functions.
 * 

 * 
 * Function objects are used in lots of places in the library (for example, in
 * step-3 we used a Functions::ZeroFunction instance as an argument to
 * VectorTools::interpolate_boundary_values) and this is the first tutorial
 * where we define a new class that inherits from Function. Since we only ever
 * call Function::value(), we could get away with just a plain function (and
 * this is what is done in step-5), but since this is a tutorial we inherit from
 * Function for the sake of example.
 * 
 * @code
 * template <int dim>
 * class RightHandSide : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 * };
 * 
 * 
 * 
 * template <int dim>
 * class BoundaryValues : public Function<dim>
 * {
 * public:
 *   virtual double value(const Point<dim> & p,
 *                        const unsigned int component = 0) const override;
 * };
 * 
 * @endcode
 * 
 * If you are not familiar with what the keywords `virtual` and `override` in
 * the function declarations above mean, you will probably want to take a look
 * at your favorite C++ book or an online tutorial such as
 * http://www.cplusplus.com/doc/tutorial/polymorphism/ . In essence, what is
 * happening here is that Function<dim> is an "abstract" base class that
 * declares a certain "interface" -- a set of functions one can call on
 * objects of this kind. But it does not actually *implement* these functions:
 * it just says "this is how Function objects look like", but what kind of
 * function it actually is, is left to derived classes that implement
 * the `value()` function.
 * 

 * 
 * Deriving one class from another is often called an "is-a" relationship
 * function. Here, the `RightHandSide` class "is a" Function class
 * because it implements the interface described by the Function base class.
 * (The actual implementation of the `value()` function is in the code block
 * below.) The `virtual` keyword then means "Yes, the
 * function here is one that can be overridden by derived classes",
 * and the `override` keyword means "Yes, this is in fact a function we know
 * has been declared as part of the base class". The `override` keyword is not
 * strictly necessary, but is an insurance against typos: If we get the name
 * of the function or the type of one argument wrong, the compiler will warn
 * us by stating "You say that this function overrides one in a base class,
 * but I don't actually know any such function with this name and these
 * arguments."
 * 

 * 
 * But back to the concrete case here:
 * For this tutorial, we choose as right hand side the function
 * $4(x^4+y^4)$ in 2D, or $4(x^4+y^4+z^4)$ in 3D. We could write this
 * distinction using an if-statement on the space dimension, but here is a
 * simple way that also allows us to use the same function in 1D (or in 4D, if
 * you should desire to do so), by using a short loop.  Fortunately, the
 * compiler knows the size of the loop at compile time (remember that at the
 * time when you define the template, the compiler doesn't know the value of
 * <code>dim</code>, but when it later encounters a statement or declaration
 * <code>RightHandSide@<2@></code>, it will take the template, replace all
 * occurrences of dim by 2 and compile the resulting function).  In other
 * words, at the time of compiling this function, the number of times the body
 * will be executed is known, and the compiler can minimize the overhead
 * needed for the loop; the result will be as fast as if we had used the
 * formulas above right away.
 * 

 * 
 * The last thing to note is that a <code>Point@<dim@></code> denotes a point
 * in dim-dimensional space, and its individual components (i.e. $x$, $y$,
 * ... coordinates) can be accessed using the () operator (in fact, the []
 * operator will work just as well) with indices starting at zero as usual in
 * C and C++.
 * 
 * @code
 * template <int dim>
 * double RightHandSide<dim>::value(const Point<dim> &p,
 *                                  const unsigned int /*component*/) const
 * {
 *   double return_value = 0.0;
 *   for (unsigned int i = 0; i < dim; ++i)
 *     return_value += 4.0 * std::pow(p(i), 4.0);
 * 
 *   return return_value;
 * }
 * 
 * 
 * @endcode
 * 
 * As boundary values, we choose $x^2+y^2$ in 2D, and $x^2+y^2+z^2$ in 3D. This
 * happens to be equal to the square of the vector from the origin to the
 * point at which we would like to evaluate the function, irrespective of the
 * dimension. So that is what we return:
 * 
 * @code
 * template <int dim>
 * double BoundaryValues<dim>::value(const Point<dim> &p,
 *                                   const unsigned int /*component*/) const
 * {
 *   return p.square();
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeStep4codeclass"></a> 
 * <h3>Implementation of the <code>Step4</code> class</h3>
 * 

 * 
 * Next for the implementation of the class template that makes use of the
 * functions above. As before, we will write everything as templates that have
 * a formal parameter <code>dim</code> that we assume unknown at the time we
 * define the template functions. Only later, the compiler will find a
 * declaration of <code>Step4@<2@></code> (in the <code>main</code> function,
 * actually) and compile the entire class with <code>dim</code> replaced by 2,
 * a process referred to as `instantiation of a template'. When doing so, it
 * will also replace instances of <code>RightHandSide@<dim@></code> by
 * <code>RightHandSide@<2@></code> and instantiate the latter class from the
 * class template.
 * 

 * 
 * In fact, the compiler will also find a declaration <code>Step4@<3@></code>
 * in <code>main()</code>. This will cause it to again go back to the general
 * <code>Step4@<dim@></code> template, replace all occurrences of
 * <code>dim</code>, this time by 3, and compile the class a second time. Note
 * that the two instantiations <code>Step4@<2@></code> and
 * <code>Step4@<3@></code> are completely independent classes; their only
 * common feature is that they are both instantiated from the same general
 * template, but they are not convertible into each other, for example, and
 * share no code (both instantiations are compiled completely independently).
 * 

 * 
 * 

 * 
 * 
 * <a name="Step4Step4"></a> 
 * <h4>Step4::Step4</h4>
 * 

 * 
 * After this introduction, here is the constructor of the <code>Step4</code>
 * class. It specifies the desired polynomial degree of the finite elements
 * and associates the DoFHandler to the triangulation just as in the previous
 * example program, step-3:
 * 
 * @code
 * template <int dim>
 * Step4<dim>::Step4()
 *   : fe(1)
 *   , dof_handler(triangulation)
 * {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4make_grid"></a> 
 * <h4>Step4::make_grid</h4>
 * 

 * 
 * Grid creation is something inherently dimension dependent. However, as long
 * as the domains are sufficiently similar in 2D or 3D, the library can
 * abstract for you. In our case, we would like to again solve on the square
 * $[-1,1]\times [-1,1]$ in 2D, or on the cube $[-1,1] \times [-1,1] \times
 * [-1,1]$ in 3D; both can be termed GridGenerator::hyper_cube(), so we may
 * use the same function in whatever dimension we are. Of course, the
 * functions that create a hypercube in two and three dimensions are very much
 * different, but that is something you need not care about. Let the library
 * handle the difficult things.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::make_grid()
 * {
 *   GridGenerator::hyper_cube(triangulation, -1, 1);
 *   triangulation.refine_global(4);
 * 
 *   std::cout << "   Number of active cells: " << triangulation.n_active_cells()
 *             << std::endl
 *             << "   Total number of cells: " << triangulation.n_cells()
 *             << std::endl;
 * }
 * 
 * @endcode
 * 
 * 
 * <a name="Step4setup_system"></a> 
 * <h4>Step4::setup_system</h4>
 * 

 * 
 * This function looks exactly like in the previous example, although it
 * performs actions that in their details are quite different if
 * <code>dim</code> happens to be 3. The only significant difference from a
 * user's perspective is the number of cells resulting, which is much higher
 * in three than in two space dimensions!
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::setup_system()
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
 * @endcode
 * 
 * 
 * <a name="Step4assemble_system"></a> 
 * <h4>Step4::assemble_system</h4>
 * 

 * 
 * Unlike in the previous example, we would now like to use a non-constant
 * right hand side function and non-zero boundary values. Both are tasks that
 * are readily achieved with only a few new lines of code in the assemblage of
 * the matrix and right hand side.
 * 

 * 
 * More interesting, though, is the way we assemble matrix and right hand side
 * vector dimension independently: there is simply no difference to the
 * two-dimensional case. Since the important objects used in this function
 * (quadrature formula, FEValues) depend on the dimension by way of a template
 * parameter as well, they can take care of setting up properly everything for
 * the dimension for which this function is compiled. By declaring all classes
 * which might depend on the dimension using a template parameter, the library
 * can make nearly all work for you and you don't have to care about most
 * things.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::assemble_system()
 * {
 *   QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 * @endcode
 * 
 * We wanted to have a non-constant right hand side, so we use an object of
 * the class declared above to generate the necessary data. Since this right
 * hand side object is only used locally in the present function, we declare
 * it here as a local variable:
 * 
 * @code
 *   RightHandSide<dim> right_hand_side;
 * 
 * @endcode
 * 
 * Compared to the previous example, in order to evaluate the non-constant
 * right hand side function we now also need the quadrature points on the
 * cell we are presently on (previously, we only required values and
 * gradients of the shape function from the FEValues object, as well as the
 * quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
 * do for us by also giving it the #update_quadrature_points flag:
 * 
 * @code
 *   FEValues<dim> fe_values(fe,
 *                           quadrature_formula,
 *                           update_values | update_gradients |
 *                             update_quadrature_points | update_JxW_values);
 * 
 * @endcode
 * 
 * We then again define the same abbreviation as in the previous program.
 * The value of this variable of course depends on the dimension which we
 * are presently using, but the FiniteElement class does all the necessary
 * work for you and you don't have to care about the dimension dependent
 * parts:
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *   Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Next, we again have to loop over all cells and assemble local
 * contributions.  Note, that a cell is a quadrilateral in two space
 * dimensions, but a hexahedron in 3D. In fact, the
 * <code>active_cell_iterator</code> data type is something different,
 * depending on the dimension we are in, but to the outside world they look
 * alike and you will probably never see a difference. In any case, the real
 * type is hidden by using `auto`:
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       fe_values.reinit(cell);
 *       cell_matrix = 0;
 *       cell_rhs    = 0;
 * 
 * @endcode
 * 
 * Now we have to assemble the local matrix and right hand side. This is
 * done exactly like in the previous example, but now we revert the
 * order of the loops (which we can safely do since they are independent
 * of each other) and merge the loops for the local matrix and the local
 * vector as far as possible to make things a bit faster.
 *       

 * 
 * Assembling the right hand side presents the only significant
 * difference to how we did things in step-3: Instead of using a
 * constant right hand side with value 1, we use the object representing
 * the right hand side and evaluate it at the quadrature points:
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             for (const unsigned int j : fe_values.dof_indices())
 *               cell_matrix(i, j) +=
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                  fe_values.JxW(q_index));           // dx
 * 
 *             const auto &x_q = fe_values.quadrature_point(q_index);
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                             right_hand_side.value(x_q) *        // f(x_q)
 *                             fe_values.JxW(q_index));            // dx
 *           }
 * @endcode
 * 
 * As a final remark to these loops: when we assemble the local
 * contributions into <code>cell_matrix(i,j)</code>, we have to multiply
 * the gradients of shape functions $i$ and $j$ at point number
 * q_index and
 * multiply it with the scalar weights JxW. This is what actually
 * happens: <code>fe_values.shape_grad(i,q_index)</code> returns a
 * <code>dim</code> dimensional vector, represented by a
 * <code>Tensor@<1,dim@></code> object, and the operator* that
 * multiplies it with the result of
 * <code>fe_values.shape_grad(j,q_index)</code> makes sure that the
 * <code>dim</code> components of the two vectors are properly
 * contracted, and the result is a scalar floating point number that
 * then is multiplied with the weights. Internally, this operator* makes
 * sure that this happens correctly for all <code>dim</code> components
 * of the vectors, whether <code>dim</code> be 2, 3, or any other space
 * dimension; from a user's perspective, this is not something worth
 * bothering with, however, making things a lot simpler if one wants to
 * write code dimension independently.
 * 

 * 
 * With the local systems assembled, the transfer into the global matrix
 * and right hand side is done exactly as before, but here we have again
 * merged some loops for efficiency:
 * 
 * @code
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
 * As the final step in this function, we wanted to have non-homogeneous
 * boundary values in this example, unlike the one before. This is a simple
 * task, we only have to replace the Functions::ZeroFunction used there by an
 * object of the class which describes the boundary values we would like to
 * use (i.e. the <code>BoundaryValues</code> class declared above):
 *   

 * 
 * The function VectorTools::interpolate_boundary_values() will only work
 * on faces that have been marked with boundary indicator 0 (because that's
 * what we say the function should work on with the second argument below).
 * If there are faces with boundary id other than 0, then the function
 * interpolate_boundary_values will do nothing on these faces. For
 * the Laplace equation doing nothing is equivalent to assuming that
 * on those parts of the boundary a zero Neumann boundary condition holds.
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values;
 *   VectorTools::interpolate_boundary_values(dof_handler,
 *                                            0,
 *                                            BoundaryValues<dim>(),
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
 * <a name="Step4solve"></a> 
 * <h4>Step4::solve</h4>
 * 

 * 
 * Solving the linear system of equations is something that looks almost
 * identical in most programs. In particular, it is dimension independent, so
 * this function is copied verbatim from the previous example.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::solve()
 * {
 *   SolverControl            solver_control(1000, 1e-12);
 *   SolverCG<Vector<double>> solver(solver_control);
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
 * 
 * @endcode
 * 
 * We have made one addition, though: since we suppress output from the
 * linear solvers, we have to print the number of iterations by hand.
 * 
 * @code
 *   std::cout << "   " << solver_control.last_step()
 *             << " CG iterations needed to obtain convergence." << std::endl;
 * }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4output_results"></a> 
 * <h4>Step4::output_results</h4>
 * 

 * 
 * This function also does what the respective one did in step-3. No changes
 * here for dimension independence either.
 * 

 * 
 * Since the program will run both 2d and 3d versions of the Laplace solver,
 * we use the dimension in the filename to generate distinct filenames for
 * each run (in a better program, one would check whether <code>dim</code> can
 * have other values than 2 or 3, but we neglect this here for the sake of
 * brevity).
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::output_results() const
 * {
 *   DataOut<dim> data_out;
 * 
 *   data_out.attach_dof_handler(dof_handler);
 *   data_out.add_data_vector(solution, "solution");
 * 
 *   data_out.build_patches();
 * 
 *   std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk");
 *   data_out.write_vtk(output);
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Step4run"></a> 
 * <h4>Step4::run</h4>
 * 

 * 
 * This is the function which has the top-level control over everything. Apart
 * from one line of additional output, it is the same as for the previous
 * example.
 * 
 * @code
 * template <int dim>
 * void Step4<dim>::run()
 * {
 *   std::cout << "Solving problem in " << dim << " space dimensions."
 *             << std::endl;
 * 
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
 * And this is the main function. It also looks mostly like in step-3, but if
 * you look at the code below, note how we first create a variable of type
 * <code>Step4@<2@></code> (forcing the compiler to compile the class template
 * with <code>dim</code> replaced by <code>2</code>) and run a 2d simulation,
 * and then we do the whole thing over in 3d.
 * 

 * 
 * In practice, this is probably not what you would do very frequently (you
 * probably either want to solve a 2d problem, or one in 3d, but not both at
 * the same time). However, it demonstrates the mechanism by which we can
 * simply change which dimension we want in a single place, and thereby force
 * the compiler to recompile the dimension independent class templates for the
 * dimension we request. The emphasis here lies on the fact that we only need
 * to change a single place. This makes it rather trivial to debug the program
 * in 2d where computations are fast, and then switch a single place to a 3 to
 * run the much more computing intensive program in 3d for `real'
 * computations.
 * 

 * 
 * Each of the two blocks is enclosed in braces to make sure that the
 * <code>laplace_problem_2d</code> variable goes out of scope (and releases
 * the memory it holds) before we move on to allocate memory for the 3d
 * case. Without the additional braces, the <code>laplace_problem_2d</code>
 * variable would only be destroyed at the end of the function, i.e. after
 * running the 3d problem, and would needlessly hog memory while the 3d run
 * could actually use it.
 * 
 * @code
 * int main()
 * {
 *   {
 *     Step4<2> laplace_problem_2d;
 *     laplace_problem_2d.run();
 *   }
 * 
 *   {
 *     Step4<3> laplace_problem_3d;
 *     laplace_problem_3d.run();
 *   }
 * 
 *   return 0;
 * }
 * @endcode
examples/step-4/doc/results.dox



<a name="Results"></a><h1>Results</h1>



程序的输出看起来如下（迭代次数可能会有一到两次的变化，这取决于你的计算机，因为这通常取决于浮点运算的舍入精度，而这在不同的处理器之间是不同的）。

@code
Solving problem in 2 space dimensions.
   Number of active cells: 256
   Total number of cells: 341
   Number of degrees of freedom: 289
   26 CG iterations needed to obtain convergence.
Solving problem in 3 space dimensions.
   Number of active cells: 4096
   Total number of cells: 4681
   Number of degrees of freedom: 4913
   30 CG iterations needed to obtain convergence.
@endcode

很明显，在三个空间维度中，单元格的数量，因此也是自由度的数量要高得多。这里看不到的是，除了矩阵中更多的行和列之外，在三个空间维度中，矩阵的每一行也有明显更多的条目。这就导致了解方程组时需要付出更多的数值努力，当你实际运行程序时，你可以从两个求解步骤的运行时间中感受到这一点。




该程序产生两个文件。   <code>solution-2d.vtk</code> 和 <code>solution-3d.vtk</code> ，可以用VisIt或Paraview程序查看（如果你没有这些程序，你可以很容易地在程序中改变输出格式，使你更容易查看）。解决方案的可视化是一门艺术，但它也可以很有趣，所以你应该玩一玩你最喜欢的可视化工具，熟悉它的功能。下面是我想出的2D解决方案。

<p align="center">  <img src="https://www.dealii.org/images/steps/developer/step-4.solution-2d.png" alt="">   </p>  。

(  @dealiiVideoLectureSeeAlso{11,32})  图片显示了所考虑的问题的解决方案，是一个三维图。可以看出，该解在域的内部几乎是平的，而在边界附近有较高的曲率。当然，这是因为对于拉普拉斯方程来说，解的曲率等于右手边，而右手边被选为四次多项式，在内部几乎为零，只有在接近域的边界时才急剧上升；右手边函数的最大值在域的角落，在那里解的移动也最迅速。很高兴看到解沿着域的边界遵循理想的二次边界值。将计算出的解与分析出的解进行验证也是很有用的。关于这一技术的解释，请参见步骤7。

另一方面，尽管图片中没有明确显示网格线，但你可以看到它们在解决方案中的小疙瘩。这清楚地表明，解决方案还没有被计算到非常高的精度，为了得到更好的解决方案，我们可能必须在更细的网格上进行计算。

在三个空间维度上，可视化就比较困难了。左图显示了解决方案和它在域的表面上计算出来的网格。这很好，但它的缺点是完全掩盖了内部的情况。右图是通过显示解的恒定值的表面（如左上角的图例所示），试图将内部的情况也可视化。如果我们把各个表面弄得稍微透明一些，这样就有可能透过它们看到后面的东西，那么等值面图片看起来就最好了。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.solution-3d.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-4.contours-3d.png" alt="">
    </td>
  </tr>
</table> 

 @note  关于可视化的最后一句话：可视化的想法是给人以洞察力，这与显示信息是不同的。特别是，在一张图片上很容易显示过多的信息，但在显示更多的信息的同时，也使人们更难收集到洞察力。举个例子，我用来生成这些图片的程序，VisIt，默认情况下在每个轴上都有刻度线，在 $x$ 轴上贴上一个大胖标签 "X轴"，其他轴也是如此，在左上方显示提取数据的文件名，在右下方显示用户的名字以及时间和日期。这些在这里都不重要：轴同样容易辨认，因为左下方的三脚架仍然可见，而且我们从程序中知道域是 $[-1,1]^3$ ，所以不需要刻度线。因此，我关掉了图片中所有不相干的东西：可视化的艺术在于把图片缩减到那些对看清自己想看的东西很重要的部分，而不是其他。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>



基本上，玩这个程序的可能性与前一个程序相同，只是它们现在也适用于3D情况。为了获得灵感，请阅读<a href="step_3.html#extensions"
target="body">possible extensions in the documentation of step 3</a>。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-4.cc"
*/
