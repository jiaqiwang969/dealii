/**
@page step_13 The step-13 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Backgroundandpurpose">Background and purpose</a>
        <li><a href="#Whattheprogramdoes">What the program does</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Evaluationofthesolution">Evaluation of the solution</a>
      <ul>
        <li><a href="#Pointevaluation">%Point evaluation</a>
        <li><a href="#Generatingoutput">Generating output</a>
        <li><a href="#Otherevaluations">Other evaluations</a>
      </ul>
        <li><a href="#TheLaplacesolverclasses">The Laplace solver classes</a>
      <ul>
        <li><a href="#Anabstractbaseclass">An abstract base class</a>
        <li><a href="#Ageneralsolverclass">A general solver class</a>
        <li><a href="#Aprimalsolver">A primal solver</a>
        <li><a href="#Globalrefinement">Global refinement</a>
        <li><a href="#LocalrefinementbytheKellyerrorindicator">Local refinement by the Kelly error indicator</a>
      </ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Thedriverroutines">The driver routines</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-13/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Backgroundandpurpose"></a><h3>Background and purpose</h3>



在这个例子程序中，我们将不太关注描述如何使用deal.II及其设施的新方法，而是介绍编写模块化和可扩展有限元程序的方法。其主要原因是现代研究软件的规模和复杂性：实现现代误差估计概念和自适应求解方法的应用程序往往变得相当大。例如，当这个程序是在2002年写的时候，deal.II的主要作者的三个最大的应用程序，在写这个例子程序的时候。<ol>  <li>  一个用非连续加尔金有限元法求解守恒双曲方程的程序。33,775行代码；  <li>  一个参数估计程序。28,980行代码；  <li>  一个波浪方程求解器：21,020行代码。   </ol> 

这个库本身--不包括示例程序和测试套件--在2002年春天有略多于15万行的代码。当然，现在已经大了好几倍了）。)这些应用程序的规模是一个人，甚至是一个有经验的程序员，所能处理的边缘。




上面的数字相当清楚地说明了一件事：没有被分解成较小的、大部分独立的片段的单体程序是没有办法生存的，因为即使是作者也会很快失去对程序不同部分之间各种依赖关系的概述。只有数据封装，例如使用面向对象的编程方法，以及通过定义小而固定的接口来实现模块化，才能帮助构建数据流和相互依赖。如果一个以上的人在开发一个程序，这也是一个绝对的先决条件，因为否则的话，混乱会很快出现，因为如果一个开发人员需要知道另一个人是否改变了不同模块的内部结构，如果它们没有被干净地分开的话。




在前面的例子中，你已经看到了库本身是如何被分成几个复合体的，每个复合体都建立在底层复合体之上，但相对独立于其他复合体。<ol>  <li>  三角形类复合物，以及相关的迭代器类； <li>  有限元类； <li>  DoFHandler类复合物，以及相关的迭代器，建立在三角形和有限元类之上； <li>  实现单元和实数单元间映射的类； <li>  FEValues类复合物，建立在有限元和映射之上。   </ol>  除了这些，还有大量的小类，当然还有以下的 "工具 "模块。<ol>  <li>  以各种图形格式输出；  <li>  线性代数类。   </ol>  这些复数也可以在deal.II手册网站的首页上以流程图的形式找到。




这个程序的目标是给出一个例子，说明一个相对简单的有限元程序的结构，使我们最终得到一组尽可能相互独立的模块。这使得我们可以在一端改变程序，而不必担心另一端的程序会被破坏，只要我们不触及两端交流的接口。当然，C++中的接口是抽象基类的声明。




在这里，我们将（再次）实现一个拉普拉斯求解器，尽管与之前的例子程序相比有一些不同。<ol>  <li>  实现数值求解方程过程的类不再负责驱动 "求解-估计误差-再求解 "的过程，而是将其委托给外部函数。这首先允许在更大的范围内将其作为一个构件，在这里拉普拉斯方程的解可能只是其中的一部分（例如，在非线性问题中，拉普拉斯方程可能要在每个非线性步骤中解决）。它还允许围绕该类建立一个框架，允许使用其他方程的求解器（但具有相同的外部接口），以备对不同类型的偏微分方程评估一些技术。   <li>  它将评估计算出的解的过程分割成一组单独的类。原因是，人们通常对偏微分方程的解本身不感兴趣，而是对它的某些方面感兴趣。例如，人们可能希望在弹性计算中计算某一边界的牵引力，或者在某一位置的接收器上计算地震波的信号。有时，人们可能对这些方面中的几个方面感兴趣。由于解的评估是通常不影响解的过程，我们把它拆成一个单独的模块，以便独立于解算器类的开发来开发这种评估过滤器。   <li>  将实现网格细化的类与计算解的类分开。   <li>  将我们要介绍的测试案例的描述与程序的其他部分分开。   <li>  使用WorkStream设施对线性系统的装配进行并行化。这是在 @ref threads "多处理器访问共享内存的并行计算 "文档模块中可以找到的广泛描述。该实现基本上遵循步骤9中已经描述过的内容。   </ol> 




该程序所做的事情并不新鲜。事实上，这更像是以前程序的混合体，从以前的例子中拆解了各种部分和功能。读者应该关注的是它们在这个程序中的安排方式，即程序中使用的软件设计技术，以达到实现所需数学方法的目的。然而，我们必须强调，软件设计在某种程度上也是一个主观的问题：不同的人有不同的编程背景，对编程的 "正确 "风格有不同的看法；因此，这个程序只表达了作者认为有用的做法，如果你对所选择的方式感到不舒服，不一定要采用这种风格来编写成功的数值软件。不过，它应该作为一个案例研究，用启发读者的想法来达到理想的目的。




一旦你完成了这个程序，你会注意到它的结构已经有些复杂了。然而，它只有大约850行代码，没有注释。在真正的应用程序中，当然会有注释和类文件，这将使其达到1200行。然而，与上面列出的应用程序相比，这仍然是很小的，因为它们的规模是它们的20到25倍。对于这么大的程序，从一开始就进行适当的设计是不可缺少的。否则，一旦它变得过于庞大而无法管理，就必须在其生命中的某一时刻重新设计它。




尽管如此，上面列出的三个程序都经历了重大的修改，甚至是重写。例如，波浪程序，在它还明显较小的时候，曾经被完全撕成了碎片，只是为了以更模块化的形式再次组装。那时，已经不可能在不影响代码的旧部分的情况下增加功能了（代码的主要问题是数据流：在时间依赖的应用中，主要的问题是什么时候把数据存储到磁盘，什么时候再重新加载；如果这不是以一种有组织的方式进行的，那么你最终会发现数据释放得太早，加载得太晚，或者根本没有释放）。尽管本例程序因此吸取了几年的经验，但它的设计肯定不是没有缺陷的，特别是可能不适合目标不同的应用。它应该作为一个灵感，让你以模块化的方式编写自己的应用程序，以避免过于紧密耦合的代码的陷阱。




<a name="Whattheprogramdoes"></a><h3>What the program does</h3>



程序实际做什么甚至不是这个程序的重点，程序的结构更重要。然而，用几句话来描述就是：求解给定右手边的拉普拉斯方程，使其解为函数  $u(x,t)=\exp(x+\sin(10y+5x^2))$  。计算的目标是得到解在点 $x_0=(0.5,0.5)$ 处的值，并比较我们在两种细化标准下解决这个值的准确性，即全局细化和通过Kelly等人的误差指标细化，我们已经在以前的例子中使用过。




像往常一样，这些结果将在本文件的相应部分进行讨论。在这样做的过程中，我们将发现一个关于两个细化标准的相对性能的略微令人恼火的观察。在以后的例子程序中，在这个程序的基础上，我们将设计一个不同的方法，希望它能比这里讨论的技术表现更好。




现在，所有的理论和传闻背景都说了这么多。了解一个项目的最好方法是看它，所以它就在这里。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * As in all programs, we start with a list of include files from the library,
 * and as usual they are in the standard order which is <code>base</code> --
 * <code>lac</code> -- <code>grid</code> -- <code>dofs</code> --
 * <code>fe</code> -- <code>numerics</code> (as each of these categories
 * roughly builds upon previous ones), then C++ standard headers:
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/table_handler.h>
 * #include <deal.II/base/thread_management.h>
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * @endcode
 * 
 * Now for the C++ standard headers:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * #include <list>
 * 
 * @endcode
 * 
 * The last step is as in all previous programs:
 * 
 * @code
 * namespace Step13
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Evaluationofthesolution"></a> 
 * <h3>Evaluation of the solution</h3>
 * 

 * 
 * As for the program itself, we first define classes that evaluate the
 * solutions of a Laplace equation. In fact, they can evaluate every kind of
 * solution, as long as it is described by a <code>DoFHandler</code> object,
 * and a solution vector. We define them here first, even before the classes
 * that actually generate the solution to be evaluated, since we need to
 * declare an abstract base class that the solver classes can refer to.
 *   

 * 
 * From an abstract point of view, we declare a pure base class that
 * provides an evaluation operator() which will do the evaluation of the
 * solution (whatever derived classes might consider an
 * <code>evaluation</code>). Since this is the only real function of this
 * base class (except for some bookkeeping machinery), one usually terms
 * such a class that only has an <code>operator()</code> a
 * <code>functor</code> in C++ terminology, since it is used just like a
 * function object.
 *   

 * 
 * Objects of this functor type will then later be passed to the solver
 * object, which applies it to the solution just computed. The evaluation
 * objects may then extract any quantity they like from the solution. The
 * advantage of putting these evaluation functions into a separate hierarchy
 * of classes is that by design they cannot use the internals of the solver
 * object and are therefore independent of changes to the way the solver
 * works. Furthermore, it is trivial to write another evaluation class
 * without modifying the solver class, which speeds up programming (not
 * being able to use internals of another class also means that you do not
 * have to worry about them -- programming evaluators is usually a rather
 * quickly done task), as well as compilation (if solver and evaluation
 * classes are put into different files: the solver only needs to see the
 * declaration of the abstract base class, and therefore does not need to be
 * recompiled upon addition of a new evaluation class, or modification of an
 * old one).  On a related note, you can reuse the evaluation classes for
 * other projects, solving different equations.
 *   

 * 
 * In order to improve separation of code into different modules, we put the
 * evaluation classes into a namespace of their own. This makes it easier to
 * actually solve different equations in the same program, by assembling it
 * from existing building blocks. The reason for this is that classes for
 * similar purposes tend to have the same name, although they were developed
 * in different contexts. In order to be able to use them together in one
 * program, it is necessary that they are placed in different
 * namespaces. This we do here:
 * 
 * @code
 *   namespace Evaluation
 *   {
 * @endcode
 * 
 * Now for the abstract base class of evaluation classes: its main purpose
 * is to declare a pure virtual function <code>operator()</code> taking a
 * <code>DoFHandler</code> object, and the solution vector. In order to be
 * able to use pointers to this base class only, it also has to declare a
 * virtual destructor, which however does nothing. Besides this, it only
 * provides for a little bit of bookkeeping: since we usually want to
 * evaluate solutions on subsequent refinement levels, we store the number
 * of the present refinement cycle, and provide a function to change this
 * number.
 * 
 * @code
 *     template <int dim>
 *     class EvaluationBase
 *     {
 *     public:
 *       virtual ~EvaluationBase() = default;
 * 
 *       void set_refinement_cycle(const unsigned int refinement_cycle);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const = 0;
 * 
 *     protected:
 *       unsigned int refinement_cycle;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step)
 *     {
 *       refinement_cycle = step;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Pointevaluation"></a> 
 * <h4>%Point evaluation</h4>
 * 

 * 
 * The next thing is to implement actual evaluation classes. As noted in
 * the introduction, we'd like to extract a point value from the solution,
 * so the first class does this in its <code>operator()</code>. The actual
 * point is given to this class through the constructor, as well as a
 * table object into which it will put its findings.
 *     

 * 
 * Finding out the value of a finite element field at an arbitrary point
 * is rather difficult, if we cannot rely on knowing the actual finite
 * element used, since then we cannot, for example, interpolate between
 * nodes. For simplicity, we therefore assume here that the point at which
 * we want to evaluate the field is actually a node. If, in the process of
 * evaluating the solution, we find that we did not encounter this point
 * upon looping over all vertices, we then have to throw an exception in
 * order to signal to the calling functions that something has gone wrong,
 * rather than silently ignore this error.
 *     

 * 
 * In the step-9 example program, we have already seen how such an
 * exception class can be declared, using the <code>DeclExceptionN</code>
 * macros. We use this mechanism here again.
 *     

 * 
 * From this, the actual declaration of this class should be evident. Note
 * that of course even if we do not list a destructor explicitly, an
 * implicit destructor is generated from the compiler, and it is virtual
 * just as the one of the base class.
 * 
 * @code
 *     template <int dim>
 *     class PointValueEvaluation : public EvaluationBase<dim>
 *     {
 *     public:
 *       PointValueEvaluation(const Point<dim> &evaluation_point,
 *                            TableHandler &    results_table);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const override;
 * 
 *       DeclException1(
 *         ExcEvaluationPointNotFound,
 *         Point<dim>,
 *         << "The evaluation point " << arg1
 *         << " was not found among the vertices of the present grid.");
 * 
 *     private:
 *       const Point<dim> evaluation_point;
 *       TableHandler &   results_table;
 *     };
 * 
 * 
 * @endcode
 * 
 * As for the definition, the constructor is trivial, just taking data and
 * storing it in object-local ones:
 * 
 * @code
 *     template <int dim>
 *     PointValueEvaluation<dim>::PointValueEvaluation(
 *       const Point<dim> &evaluation_point,
 *       TableHandler &    results_table)
 *       : evaluation_point(evaluation_point)
 *       , results_table(results_table)
 *     {}
 * 
 * 
 * 
 * @endcode
 * 
 * Now for the function that is mainly of interest in this class, the
 * computation of the point value:
 * 
 * @code
 *     template <int dim>
 *     void PointValueEvaluation<dim>::
 *          operator()(const DoFHandler<dim> &dof_handler,
 *                const Vector<double> & solution) const
 *     {
 * @endcode
 * 
 * First allocate a variable that will hold the point value. Initialize
 * it with a value that is clearly bogus, so that if we fail to set it
 * to a reasonable value, we will note at once. This may not be
 * necessary in a function as small as this one, since we can easily see
 * all possible paths of execution here, but it proved to be helpful for
 * more complex cases, and so we employ this strategy here as well.
 * 
 * @code
 *       double point_value = 1e20;
 * 
 * @endcode
 * 
 * Then loop over all cells and all their vertices, and check whether a
 * vertex matches the evaluation point. If this is the case, then
 * extract the point value, set a flag that we have found the point of
 * interest, and exit the loop.
 * 
 * @code
 *       bool evaluation_point_found = false;
 *       for (const auto &cell : dof_handler.active_cell_iterators())
 *         if (!evaluation_point_found)
 *           for (const auto vertex : cell->vertex_indices())
 *             if (cell->vertex(vertex) == evaluation_point)
 *               {
 * @endcode
 * 
 * In order to extract the point value from the global solution
 * vector, pick that component that belongs to the vertex of
 * interest, and, in case the solution is vector-valued, take
 * the first component of it:
 * 
 * @code
 *                 point_value = solution(cell->vertex_dof_index(vertex, 0));
 * @endcode
 * 
 * Note that by this we have made an assumption that is not
 * valid always and should be documented in the class
 * declaration if this were code for a real application rather
 * than a tutorial program: we assume that the finite element
 * used for the solution we try to evaluate actually has degrees
 * of freedom associated with vertices. This, for example, does
 * not hold for discontinuous elements, were the support points
 * for the shape functions happen to be located at the vertices,
 * but are not associated with the vertices but rather with the
 * cell interior, since association with vertices would imply
 * continuity there. It would also not hold for edge oriented
 * elements, and the like.
 *                 

 * 
 * Ideally, we would check this at the beginning of the
 * function, for example by a statement like <code>Assert
 * (dof_handler.get_fe().dofs_per_vertex @> 0,
 * ExcNotImplemented())</code>, which should make it quite clear
 * what is going wrong when the exception is triggered. In this
 * case, we omit it (which is indeed bad style), but knowing
 * that that does not hurt here, since the statement
 * <code>cell-@>vertex_dof_index(vertex,0)</code> would fail if
 * we asked it to give us the DoF index of a vertex if there
 * were none.
 *                 

 * 
 * We stress again that this restriction on the allowed finite
 * elements should be stated in the class documentation.
 * 

 * 
 * Since we found the right point, we now set the respective
 * flag and exit the innermost loop. The outer loop will also be
 * terminated due to the set flag.
 * 
 * @code
 *                 evaluation_point_found = true;
 *                 break;
 *               };
 * 
 * @endcode
 * 
 * Finally, we'd like to make sure that we have indeed found the
 * evaluation point, since if that were not so we could not give a
 * reasonable value of the solution there and the rest of the
 * computations were useless anyway. So make sure through the
 * <code>AssertThrow</code> macro already used in the step-9 program
 * that we have indeed found this point. If this is not so, the macro
 * throws an exception of the type that is given to it as second
 * argument, but compared to a straightforward <code>throw</code>
 * statement, it fills the exception object with a set of additional
 * information, for example the source file and line number where the
 * exception was generated, and the condition that failed. If you have a
 * <code>catch</code> clause in your main function (as this program
 * has), you will catch all exceptions that are not caught somewhere in
 * between and thus already handled, and this additional information
 * will help you find out what happened and where it went wrong.
 * 
 * @code
 *       AssertThrow(evaluation_point_found,
 *                   ExcEvaluationPointNotFound(evaluation_point));
 * @endcode
 * 
 * Note that we have used the <code>Assert</code> macro in other example
 * programs as well. It differed from the <code>AssertThrow</code> macro
 * used here in that it simply aborts the program, rather than throwing
 * an exception, and that it did so only in debug mode. It was the right
 * macro to use to check about the size of vectors passed as arguments
 * to functions, and the like.
 *       

 * 
 * However, here the situation is different: whether we find the
 * evaluation point or not may change from refinement to refinement (for
 * example, if the four cells around point are coarsened away, then the
 * point may vanish after refinement and coarsening). This is something
 * that cannot be predicted from a few number of runs of the program in
 * debug mode, but should be checked always, also in production
 * runs. Thus the use of the <code>AssertThrow</code> macro here.
 * 

 * 
 * Now, if we are sure that we have found the evaluation point, we can
 * add the results into the table of results:
 * 
 * @code
 *       results_table.add_value("DoFs", dof_handler.n_dofs());
 *       results_table.add_value("u(x_0)", point_value);
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Generatingoutput"></a> 
 * <h4>Generating output</h4>
 * 

 * 
 * A different, maybe slightly odd kind of <code>evaluation</code> of a
 * solution is to output it to a file in a graphical format. Since in the
 * evaluation functions we are given a <code>DoFHandler</code> object and
 * the solution vector, we have all we need to do this, so we can do it in
 * an evaluation class. The reason for actually doing so instead of
 * putting it into the class that computed the solution is that this way
 * we have more flexibility: if we choose to only output certain aspects
 * of it, or not output it at all. In any case, we do not need to modify
 * the solver class, we just have to modify one of the modules out of
 * which we build this program. This form of encapsulation, as above,
 * helps us to keep each part of the program rather simple as the
 * interfaces are kept simple, and no access to hidden data is possible.
 *     

 * 
 * Since this class which generates the output is derived from the common
 * <code>EvaluationBase</code> base class, its main interface is the
 * <code>operator()</code> function. Furthermore, it has a constructor
 * taking a string that will be used as the base part of the file name to
 * which output will be sent (we will augment it by a number indicating
 * the number of the refinement cycle -- the base class has this
 * information at hand --, and a suffix), and the constructor also takes a
 * value that indicates which format is requested, i.e. for which graphics
 * program we shall generate output (from this we will then also generate
 * the suffix of the filename to which we write).
 *     

 * 
 * Regarding the output format, the DataOutBase namespace
 * provides an enumeration field
 * DataOutBase::OutputFormat which lists names for all supported output
 * formats. At the time of writing of this program, the supported graphics
 * formats are represented by the enum values <code>ucd</code>,
 * <code>gnuplot</code>, <code>povray</code>, <code>eps</code>,
 * <code>gmv</code>, <code>tecplot</code>, <code>tecplot_binary</code>,
 * <code>dx</code>, <code>vtk</code>, etc, but this list will certainly
 * grow over time. Now, within various functions of that base class, you
 * can use values of this type to get information about these graphics
 * formats (for example the default suffix used for files of each format),
 * and you can call a generic <code>write</code> function, which then
 * branches to the <code>write_gnuplot</code>, <code>write_ucd</code>, etc
 * functions which we have used in previous examples already, based on the
 * value of a second argument given to it denoting the required output
 * format. This mechanism makes it simple to write an extensible program
 * that can decide which output format to use at runtime, and it also
 * makes it rather simple to write the program in a way such that it takes
 * advantage of newly implemented output formats, without the need to
 * change the application program.
 *     

 * 
 * Of these two fields, the base name and the output format descriptor,
 * the constructor takes values and stores them for later use by the
 * actual evaluation function.
 * 
 * @code
 *     template <int dim>
 *     class SolutionOutput : public EvaluationBase<dim>
 *     {
 *     public:
 *       SolutionOutput(const std::string &             output_name_base,
 *                      const DataOutBase::OutputFormat output_format);
 * 
 *       virtual void operator()(const DoFHandler<dim> &dof_handler,
 *                               const Vector<double> & solution) const override;
 * 
 *     private:
 *       const std::string               output_name_base;
 *       const DataOutBase::OutputFormat output_format;
 *     };
 * 
 * 
 *     template <int dim>
 *     SolutionOutput<dim>::SolutionOutput(
 *       const std::string &             output_name_base,
 *       const DataOutBase::OutputFormat output_format)
 *       : output_name_base(output_name_base)
 *       , output_format(output_format)
 *     {}
 * 
 * 
 * @endcode
 * 
 * Following the description above, the function generating the actual
 * output is now relatively straightforward. The only particularly
 * interesting feature over previous example programs is the use of the
 * DataOutBase::default_suffix function, returning the usual
 * suffix for files of a given format (e.g. ".eps" for encapsulated
 * postscript files, ".gnuplot" for Gnuplot files), and of the generic
 * DataOut::write() function with a second argument, which internally
 * branches to the actual output functions for the different graphics
 * formats, based on the value of the format descriptor passed as second
 * argument.
 *     

 * 
 * Also note that we have to prefix <code>this-@></code> to access a
 * member variable of the template dependent base class. The reason here,
 * and further down in the program is the same as the one described in the
 * step-7 example program (look for <code>two-stage name lookup</code>
 * there).
 * 
 * @code
 *     template <int dim>
 *     void SolutionOutput<dim>::operator()(const DoFHandler<dim> &dof_handler,
 *                                          const Vector<double> & solution) const
 *     {
 *       DataOut<dim> data_out;
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "solution");
 *       data_out.build_patches();
 * 
 *       std::ofstream out(output_name_base + "-" +
 *                         std::to_string(this->refinement_cycle) +
 *                         data_out.default_suffix(output_format));
 * 
 *       data_out.write(out, output_format);
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Otherevaluations"></a> 
 * <h4>Other evaluations</h4>
 * 

 * 
 * In practical applications, one would add here a list of other possible
 * evaluation classes, representing quantities that one may be interested
 * in. For this example, that much shall be sufficient, so we close the
 * namespace.
 * 
 * @code
 *   } // namespace Evaluation
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheLaplacesolverclasses"></a> 
 * <h3>The Laplace solver classes</h3>
 * 

 * 
 * After defining what we want to know of the solution, we should now care
 * how to get at it. We will pack everything we need into a namespace of its
 * own, for much the same reasons as for the evaluations above.
 *   

 * 
 * Since we have discussed Laplace solvers already in considerable detail in
 * previous examples, there is not much new stuff following. Rather, we have
 * to a great extent cannibalized previous examples and put them, in
 * slightly different form, into this example program. We will therefore
 * mostly be concerned with discussing the differences to previous examples.
 *   

 * 
 * Basically, as already said in the introduction, the lack of new stuff in
 * this example is deliberate, as it is more to demonstrate software design
 * practices, rather than mathematics. The emphasis in explanations below
 * will therefore be more on the actual implementation.
 * 
 * @code
 *   namespace LaplaceSolver
 *   {
 * @endcode
 * 
 * 
 * <a name="Anabstractbaseclass"></a> 
 * <h4>An abstract base class</h4>
 * 

 * 
 * In defining a Laplace solver, we start out by declaring an abstract
 * base class, that has no functionality itself except for taking and
 * storing a pointer to the triangulation to be used later.
 *     

 * 
 * This base class is very general, and could as well be used for any
 * other stationary problem. It provides declarations of functions that
 * shall, in derived classes, solve a problem, postprocess the solution
 * with a list of evaluation objects, and refine the grid,
 * respectively. None of these functions actually does something itself in
 * the base class.
 *     

 * 
 * Due to the lack of actual functionality, the programming style of
 * declaring very abstract base classes is similar to the style used in
 * Smalltalk or Java programs, where all classes are derived from entirely
 * abstract classes <code>Object</code>, even number representations. The
 * author admits that he does not particularly like the use of such a
 * style in C++, as it puts style over reason. Furthermore, it promotes
 * the use of virtual functions for everything (for example, in Java, all
 * functions are virtual per se), which, however, has proven to be rather
 * inefficient in many applications where functions are often only
 * accessing data, not doing computations, and therefore quickly return;
 * the overhead of virtual functions can then be significant. The opinion
 * of the author is to have abstract base classes wherever at least some
 * part of the code of actual implementations can be shared and thus
 * separated into the base class.
 *     

 * 
 * Besides all these theoretical questions, we here have a good reason,
 * which will become clearer to the reader below. Basically, we want to be
 * able to have a family of different Laplace solvers that differ so much
 * that no larger common subset of functionality could be found. We
 * therefore just declare such an abstract base class, taking a pointer to
 * a triangulation in the constructor and storing it henceforth. Since
 * this triangulation will be used throughout all computations, we have to
 * make sure that the triangulation is valid until it is last used. We
 * do this by keeping a <code>SmartPointer</code> to this triangulation,
 * as explained in step-7.
 *     

 * 
 * Note that while the pointer itself is declared constant
 * (i.e. throughout the lifetime of this object, the pointer points to the
 * same object), it is not declared as a pointer to a constant
 * triangulation. In fact, by this we allow that derived classes refine or
 * coarsen the triangulation within the <code>refine_grid</code> function.
 *     

 * 
 * Finally, we have a function <code>n_dofs</code> is only a tool for the
 * driver functions to decide whether we want to go on with mesh
 * refinement or not. It returns the number of degrees of freedom the
 * present simulation has.
 * 
 * @code
 *     template <int dim>
 *     class Base
 *     {
 *     public:
 *       Base(Triangulation<dim> &coarse_grid);
 *       virtual ~Base() = default;
 * 
 *       virtual void solve_problem() = 0;
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
 *       virtual void         refine_grid()                            = 0;
 *       virtual unsigned int n_dofs() const                           = 0;
 * 
 *     protected:
 *       const SmartPointer<Triangulation<dim>> triangulation;
 *     };
 * 
 * 
 * @endcode
 * 
 * The implementation of the only two non-abstract functions is then
 * rather boring:
 * 
 * @code
 *     template <int dim>
 *     Base<dim>::Base(Triangulation<dim> &coarse_grid)
 *       : triangulation(&coarse_grid)
 *     {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Ageneralsolverclass"></a> 
 * <h4>A general solver class</h4>
 * 

 * 
 * Following now the main class that implements assembling the matrix of
 * the linear system, solving it, and calling the postprocessor objects on
 * the solution. It implements the <code>solve_problem</code> and
 * <code>postprocess</code> functions declared in the base class. It does
 * not, however, implement the <code>refine_grid</code> method, as mesh
 * refinement will be implemented in a number of derived classes.
 *     

 * 
 * It also declares a new abstract virtual function,
 * <code>assemble_rhs</code>, that needs to be overloaded in
 * subclasses. The reason is that we will implement two different classes
 * that will implement different methods to assemble the right hand side
 * vector. This function might also be interesting in cases where the
 * right hand side depends not simply on a continuous function, but on
 * something else as well, for example the solution of another discretized
 * problem, etc. The latter happens frequently in non-linear problems.
 *     

 * 
 * As we mentioned previously, the actual content of this class is not
 * new, but a mixture of various techniques already used in previous
 * examples. We will therefore not discuss them in detail, but refer the
 * reader to these programs.
 *     

 * 
 * Basically, in a few words, the constructor of this class takes pointers
 * to a triangulation, a finite element, and a function object
 * representing the boundary values. These are either passed down to the
 * base class's constructor, or are stored and used to generate a
 * <code>DoFHandler</code> object later. Since finite elements and
 * quadrature formula should match, it is also passed a quadrature object.
 *     

 * 
 * The <code>solve_problem</code> sets up the data structures for the
 * actual solution, calls the functions to assemble the linear system, and
 * solves it.
 *     

 * 
 * The <code>postprocess</code> function finally takes an evaluation
 * object and applies it to the computed solution.
 *     

 * 
 * The <code>n_dofs</code> function finally implements the pure virtual
 * function of the base class.
 * 
 * @code
 *     template <int dim>
 *     class Solver : public virtual Base<dim>
 *     {
 *     public:
 *       Solver(Triangulation<dim> &      triangulation,
 *              const FiniteElement<dim> &fe,
 *              const Quadrature<dim> &   quadrature,
 *              const Function<dim> &     boundary_values);
 *       virtual ~Solver() override;
 * 
 *       virtual void solve_problem() override;
 * 
 *       virtual void postprocess(
 *         const Evaluation::EvaluationBase<dim> &postprocessor) const override;
 * 
 *       virtual unsigned int n_dofs() const override;
 * 
 * @endcode
 * 
 * In the protected section of this class, we first have a number of
 * member variables, of which the use should be clear from the previous
 * examples:
 * 
 * @code
 *     protected:
 *       const SmartPointer<const FiniteElement<dim>> fe;
 *       const SmartPointer<const Quadrature<dim>>    quadrature;
 *       DoFHandler<dim>                              dof_handler;
 *       Vector<double>                               solution;
 *       const SmartPointer<const Function<dim>>      boundary_values;
 * 
 * @endcode
 * 
 * Then we declare an abstract function that will be used to assemble
 * the right hand side. As explained above, there are various cases for
 * which this action differs strongly in what is necessary, so we defer
 * this to derived classes:
 * 
 * @code
 *       virtual void assemble_rhs(Vector<double> &rhs) const = 0;
 * 
 * @endcode
 * 
 * Next, in the private section, we have a small class which represents
 * an entire linear system, i.e. a matrix, a right hand side, and a
 * solution vector, as well as the constraints that are applied to it,
 * such as those due to hanging nodes. Its constructor initializes the
 * various subobjects, and there is a function that implements a
 * conjugate gradient method as solver.
 * 
 * @code
 *     private:
 *       struct LinearSystem
 *       {
 *         LinearSystem(const DoFHandler<dim> &dof_handler);
 * 
 *         void solve(Vector<double> &solution) const;
 * 
 *         AffineConstraints<double> hanging_node_constraints;
 *         SparsityPattern           sparsity_pattern;
 *         SparseMatrix<double>      matrix;
 *         Vector<double>            rhs;
 *       };
 * 
 * 
 * @endcode
 * 
 * Finally, there is a set of functions which will be used to
 * assemble the actual system matrix. The main function of this
 * group, <code>assemble_linear_system()</code> computes the
 * matrix in parallel on multicore systems, using the following
 * two helper functions. The mechanism for doing so is the same
 * as in the step-9 example program and follows the WorkStream
 * concept outlined in @ref threads . The main function also
 * calls the virtual function assembling the right hand side.
 * 
 * @code
 *       struct AssemblyScratchData
 *       {
 *         AssemblyScratchData(const FiniteElement<dim> &fe,
 *                             const Quadrature<dim> &   quadrature);
 *         AssemblyScratchData(const AssemblyScratchData &scratch_data);
 * 
 *         FEValues<dim> fe_values;
 *       };
 * 
 *       struct AssemblyCopyData
 *       {
 *         FullMatrix<double>                   cell_matrix;
 *         std::vector<types::global_dof_index> local_dof_indices;
 *       };
 * 
 *       void assemble_linear_system(LinearSystem &linear_system);
 * 
 *       void local_assemble_matrix(
 *         const typename DoFHandler<dim>::active_cell_iterator &cell,
 *         AssemblyScratchData &                                 scratch_data,
 *         AssemblyCopyData &                                    copy_data) const;
 * 
 *       void copy_local_to_global(const AssemblyCopyData &copy_data,
 *                                 LinearSystem &          linear_system) const;
 *     };
 * 
 * 
 * 
 * @endcode
 * 
 * Now here comes the constructor of the class. It does not do much except
 * store pointers to the objects given, and generate
 * <code>DoFHandler</code> object initialized with the given pointer to a
 * triangulation. This causes the DoF handler to store that pointer, but
 * does not already generate a finite element numbering (we only ask for
 * that in the <code>solve_problem</code> function).
 * 
 * @code
 *     template <int dim>
 *     Solver<dim>::Solver(Triangulation<dim> &      triangulation,
 *                         const FiniteElement<dim> &fe,
 *                         const Quadrature<dim> &   quadrature,
 *                         const Function<dim> &     boundary_values)
 *       : Base<dim>(triangulation)
 *       , fe(&fe)
 *       , quadrature(&quadrature)
 *       , dof_handler(triangulation)
 *       , boundary_values(&boundary_values)
 *     {}
 * 
 * 
 * @endcode
 * 
 * The destructor is simple, it only clears the information stored in the
 * DoF handler object to release the memory.
 * 
 * @code
 *     template <int dim>
 *     Solver<dim>::~Solver()
 *     {
 *       dof_handler.clear();
 *     }
 * 
 * 
 * @endcode
 * 
 * The next function is the one which delegates the main work in solving
 * the problem: it sets up the DoF handler object with the finite element
 * given to the constructor of this object, the creates an object that
 * denotes the linear system (i.e. the matrix, the right hand side vector,
 * and the solution vector), calls the function to assemble it, and
 * finally solves it:
 * 
 * @code
 *     template <int dim>
 *     void Solver<dim>::solve_problem()
 *     {
 *       dof_handler.distribute_dofs(*fe);
 *       solution.reinit(dof_handler.n_dofs());
 * 
 *       LinearSystem linear_system(dof_handler);
 *       assemble_linear_system(linear_system);
 *       linear_system.solve(solution);
 *     }
 * 
 * 
 * @endcode
 * 
 * As stated above, the <code>postprocess</code> function takes an
 * evaluation object, and applies it to the computed solution. This
 * function may be called multiply, once for each evaluation of the
 * solution which the user required.
 * 
 * @code
 *     template <int dim>
 *     void Solver<dim>::postprocess(
 *       const Evaluation::EvaluationBase<dim> &postprocessor) const
 *     {
 *       postprocessor(dof_handler, solution);
 *     }
 * 
 * 
 * @endcode
 * 
 * The <code>n_dofs</code> function should be self-explanatory:
 * 
 * @code
 *     template <int dim>
 *     unsigned int Solver<dim>::n_dofs() const
 *     {
 *       return dof_handler.n_dofs();
 *     }
 * 
 * 
 * @endcode
 * 
 * The following function assembles matrix and right hand side of
 * the linear system to be solved in each step. We will do things
 * in parallel at a couple of levels. First, note that we need to
 * assemble both the matrix and the right hand side. These are
 * independent operations, and we should do this in parallel. To
 * this end, we use the concept of "tasks" that is discussed in
 * the @ref threads documentation module. In essence, what we want
 * to say "here is something that needs to be worked on, go do it
 * whenever a CPU core is available", then do something else, and
 * when we need the result of the first operation wait for its
 * completion. At the second level, we want to assemble the matrix
 * using the exact same strategy we have already used in step-9,
 * namely the WorkStream concept.
 *     

 * 
 * While we could consider either assembling the right hand side
 * or assembling the matrix as the thing to do in the background
 * while doing the other, we will opt for the former approach
 * simply because the call to <code>Solver::assemble_rhs</code> is
 * so much simpler to write than the call to WorkStream::run with
 * its many arguments. In any case, the code then looks like this
 * to assemble the entire linear system:
 * 
 * @code
 *     template <int dim>
 *     void Solver<dim>::assemble_linear_system(LinearSystem &linear_system)
 *     {
 *       Threads::Task<void> rhs_task =
 *         Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs);
 * 
 *       auto worker =
 *         [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                AssemblyScratchData &scratch_data,
 *                AssemblyCopyData &   copy_data) {
 *           this->local_assemble_matrix(cell, scratch_data, copy_data);
 *         };
 * 
 *       auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) {
 *         this->copy_local_to_global(copy_data, linear_system);
 *       };
 * 
 *       WorkStream::run(dof_handler.begin_active(),
 *                       dof_handler.end(),
 *                       worker,
 *                       copier,
 *                       AssemblyScratchData(*fe, *quadrature),
 *                       AssemblyCopyData());
 *       linear_system.hanging_node_constraints.condense(linear_system.matrix);
 * 
 * @endcode
 * 
 * The syntax above requires
 * some explanation. There are multiple version of
 * WorkStream::run that expect different arguments. In step-9,
 * we used one version that took a pair of iterators, a pair of
 * pointers to member functions with very specific argument
 * lists, a pointer or reference to the object on which these
 * member functions have to work, and a scratch and copy data
 * object. This is a bit restrictive since the member functions
 * called this way have to have an argument list that exactly
 * matches what WorkStream::run expects: the local assembly
 * function needs to take an iterator, a scratch object and a
 * copy object; and the copy-local-to-global function needs to
 * take exactly a copy object. But, what if we want something
 * that's slightly more general? For example, in the current
 * program, the copy-local-to-global function needs to know
 * which linear system object to write the local contributions
 * into, i.e., it also has to take a <code>LinearSystem</code>
 * argument. That won't work with the approach using member
 * function pointers.
 *       

 * 
 * Fortunately, C++ offers a way out. These are called function
 * objects. In essence, what WorkStream::run wants to do is not
 * call a member function. It wants to call some function that
 * takes an iterator, a scratch object and a copy object in the
 * first case, and a copy object in the second case. Whether
 * these are member functions, global functions, or something
 * else, is really not of much concern to
 * WorkStream. Consequently, there is a second version of the
 * function that just takes function objects -- objects that
 * have an <code>operator()</code> and that consequently can be
 * called like functions, whatever they really represent. The
 * typical way to generate such function objects is using a
 * <a href="http://en.wikipedia.org/wiki/Anonymous_function">lambda
 * function</a> that wraps the function call including the individual
 * arguments with fixed values. All the arguments that are part of the
 * outer function signature are specified as regular function arguments in
 * the lambda function. The fixed values are passed into the lambda
 * function using the capture list (`[...]`). It is possible to use a
 * capture default or to list all the variables that are to be bound to
 * the lambda explicitly. For the sake of clarity we decided to omit
 * the capture default here, but that capture list could equally well be
 * `[&]`, meaning that all used variables are copied into the lambda by
 * reference.
 *       

 * 
 * At this point, we have assembled the matrix and condensed
 * it. The right hand side may or may not have been completely
 * assembled, but we would like to condense the right hand side
 * vector next. We can only do this if the assembly of this
 * vector has finished, so we have to wait for the task to
 * finish; in computer science, waiting for a task is typically
 * called "joining" the task, explaining the name of the
 * function we call below.
 *       

 * 
 * Since that task may or may not have finished, and since we
 * may have to wait for it to finish, we may as well try to pack
 * other things that need to be done anyway into this
 * gap. Consequently, we first interpolate boundary values
 * before we wait for the right hand side. Of course, another
 * possibility would have been to also interpolate the boundary
 * values on a separate task since doing so is independent of
 * the other things we have done in this function so far. Feel
 * free to find the correct syntax to also create a task for
 * this interpolation and start it at the top of this function,
 * along with the assembly of the right hand side. (You will
 * find that this is slightly more complicated since there are
 * multiple versions of
 * VectorTools::interpolate_boundary_values(), and so simply
 * taking the address
 * <code>&VectorTools::interpolate_boundary_values</code>
 * produces a set of overloaded functions that can't be passed
 * to Threads::new_task() right away -- you have to select which
 * element of this overload set you want by casting the address
 * expression to a function pointer type that is specific to the
 * version of the function that you want to call on the task.)
 * 
 * @code
 *       std::map<types::global_dof_index, double> boundary_value_map;
 *       VectorTools::interpolate_boundary_values(dof_handler,
 *                                                0,
 *                                                *boundary_values,
 *                                                boundary_value_map);
 * 
 *       rhs_task.join();
 *       linear_system.hanging_node_constraints.condense(linear_system.rhs);
 * 
 * @endcode
 * 
 * Now that we have the complete linear system, we can also
 * treat boundary values, which need to be eliminated from both
 * the matrix and the right hand side:
 * 
 * @code
 *       MatrixTools::apply_boundary_values(boundary_value_map,
 *                                          linear_system.matrix,
 *                                          solution,
 *                                          linear_system.rhs);
 *     }
 * 
 * 
 * @endcode
 * 
 * The second half of this set of functions deals with the local
 * assembly on each cell and copying local contributions into the
 * global matrix object. This works in exactly the same way as
 * described in step-9:
 * 
 * @code
 *     template <int dim>
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData(
 *       const FiniteElement<dim> &fe,
 *       const Quadrature<dim> &   quadrature)
 *       : fe_values(fe, quadrature, update_gradients | update_JxW_values)
 *     {}
 * 
 * 
 *     template <int dim>
 *     Solver<dim>::AssemblyScratchData::AssemblyScratchData(
 *       const AssemblyScratchData &scratch_data)
 *       : fe_values(scratch_data.fe_values.get_fe(),
 *                   scratch_data.fe_values.get_quadrature(),
 *                   update_gradients | update_JxW_values)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::local_assemble_matrix(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       AssemblyScratchData &                                 scratch_data,
 *       AssemblyCopyData &                                    copy_data) const
 *     {
 *       const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
 *       const unsigned int n_q_points    = quadrature->size();
 * 
 *       copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
 * 
 *       copy_data.local_dof_indices.resize(dofs_per_cell);
 * 
 *       scratch_data.fe_values.reinit(cell);
 * 
 *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *             copy_data.cell_matrix(i, j) +=
 *               (scratch_data.fe_values.shape_grad(i, q_point) *
 *                scratch_data.fe_values.shape_grad(j, q_point) *
 *                scratch_data.fe_values.JxW(q_point));
 * 
 *       cell->get_dof_indices(copy_data.local_dof_indices);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data,
 *                                            LinearSystem &linear_system) const
 *     {
 *       for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i)
 *         for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j)
 *           linear_system.matrix.add(copy_data.local_dof_indices[i],
 *                                    copy_data.local_dof_indices[j],
 *                                    copy_data.cell_matrix(i, j));
 *     }
 * 
 * 
 * @endcode
 * 
 * Now for the functions that implement actions in the linear system
 * class. First, the constructor initializes all data elements to their
 * correct sizes, and sets up a number of additional data structures, such
 * as constraints due to hanging nodes. Since setting up the hanging nodes
 * and finding out about the nonzero elements of the matrix is
 * independent, we do that in parallel (if the library was configured to
 * use concurrency, at least; otherwise, the actions are performed
 * sequentially). Note that we start only one thread, and do the second
 * action in the main thread. Since only one task is generated, we don't
 * use the <code>Threads::TaskGroup</code> class here, but rather use
 * the one created task object directly to wait for this particular
 * task's exit.
 *     

 * 
 * Note that taking up the address of the
 * <code>DoFTools::make_hanging_node_constraints</code> function is a
 * little tricky, since there are actually three of them, one for each
 * supported space dimension. Taking addresses of overloaded functions is
 * somewhat complicated in C++, since the address-of operator
 * <code>&</code> in that case returns more like a set of values (the
 * addresses of all functions with that name), and selecting the right one
 * is then the next step. If the context dictates which one to take (for
 * example by assigning to a function pointer of known type), then the
 * compiler can do that by itself, but if this set of pointers shall be
 * given as the argument to a function that takes a template, the compiler
 * could choose all without having a preference for one. We therefore have
 * to make it clear to the compiler which one we would like to have; for
 * this, we could use a cast, but for more clarity, we assign it to a
 * temporary <code>mhnc_p</code> (short for <code>pointer to
 * make_hanging_node_constraints</code>) with the right type, and using
 * this pointer instead.
 * 
 * @code
 *     template <int dim>
 *     Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler)
 *     {
 *       hanging_node_constraints.clear();
 * 
 *       void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) =
 *         &DoFTools::make_hanging_node_constraints;
 * 
 * @endcode
 * 
 * Start a side task then continue on the main thread
 * 
 * @code
 *       Threads::Task<void> side_task =
 *         Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints);
 * 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp);
 * 
 * 
 * 
 * @endcode
 * 
 * Wait for the side task to be done before going further
 * 
 * @code
 *       side_task.join();
 * 
 *       hanging_node_constraints.close();
 *       hanging_node_constraints.condense(dsp);
 *       sparsity_pattern.copy_from(dsp);
 * 
 * 
 * @endcode
 * 
 * Finally initialize the matrix and right hand side vector
 * 
 * @code
 *       matrix.reinit(sparsity_pattern);
 *       rhs.reinit(dof_handler.n_dofs());
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * The second function of this class simply solves the linear system by a
 * preconditioned conjugate gradient method. This has been extensively
 * discussed before, so we don't dwell into it any more.
 * 
 * @code
 *     template <int dim>
 *     void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const
 *     {
 *       SolverControl            solver_control(1000, 1e-12);
 *       SolverCG<Vector<double>> cg(solver_control);
 * 
 *       PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *       preconditioner.initialize(matrix, 1.2);
 * 
 *       cg.solve(matrix, solution, rhs, preconditioner);
 * 
 *       hanging_node_constraints.distribute(solution);
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Aprimalsolver"></a> 
 * <h4>A primal solver</h4>
 * 

 * 
 * In the previous section, a base class for Laplace solvers was
 * implemented, that lacked the functionality to assemble the right hand
 * side vector, however, for reasons that were explained there. Now we
 * implement a corresponding class that can do this for the case that the
 * right hand side of a problem is given as a function object.
 *     

 * 
 * The actions of the class are rather what you have seen already in
 * previous examples already, so a brief explanation should suffice: the
 * constructor takes the same data as does that of the underlying class
 * (to which it passes all information) except for one function object
 * that denotes the right hand side of the problem. A pointer to this
 * object is stored (again as a <code>SmartPointer</code>, in order to
 * make sure that the function object is not deleted as long as it is
 * still used by this class).
 *     

 * 
 * The only functional part of this class is the <code>assemble_rhs</code>
 * method that does what its name suggests.
 * 
 * @code
 *     template <int dim>
 *     class PrimalSolver : public Solver<dim>
 *     {
 *     public:
 *       PrimalSolver(Triangulation<dim> &      triangulation,
 *                    const FiniteElement<dim> &fe,
 *                    const Quadrature<dim> &   quadrature,
 *                    const Function<dim> &     rhs_function,
 *                    const Function<dim> &     boundary_values);
 * 
 *     protected:
 *       const SmartPointer<const Function<dim>> rhs_function;
 *       virtual void assemble_rhs(Vector<double> &rhs) const override;
 *     };
 * 
 * 
 * @endcode
 * 
 * The constructor of this class basically does what it is announced to do
 * above...
 * 
 * @code
 *     template <int dim>
 *     PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &      triangulation,
 *                                     const FiniteElement<dim> &fe,
 *                                     const Quadrature<dim> &   quadrature,
 *                                     const Function<dim> &     rhs_function,
 *                                     const Function<dim> &     boundary_values)
 *       : Base<dim>(triangulation)
 *       , Solver<dim>(triangulation, fe, quadrature, boundary_values)
 *       , rhs_function(&rhs_function)
 *     {}
 * 
 * 
 * 
 * @endcode
 * 
 * ... as does the <code>assemble_rhs</code> function. Since this is
 * explained in several of the previous example programs, we leave it at
 * that.
 * 
 * @code
 *     template <int dim>
 *     void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const
 *     {
 *       FEValues<dim> fe_values(*this->fe,
 *                               *this->quadrature,
 *                               update_values | update_quadrature_points |
 *                                 update_JxW_values);
 * 
 *       const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell();
 *       const unsigned int n_q_points    = this->quadrature->size();
 * 
 *       Vector<double>                       cell_rhs(dofs_per_cell);
 *       std::vector<double>                  rhs_values(n_q_points);
 *       std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *       for (const auto &cell : this->dof_handler.active_cell_iterators())
 *         {
 *           cell_rhs = 0;
 *           fe_values.reinit(cell);
 *           rhs_function->value_list(fe_values.get_quadrature_points(),
 *                                    rhs_values);
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *               cell_rhs(i) += fe_values.shape_value(i, q_point) * 
 *                              rhs_values[q_point] *               
 *                              fe_values.JxW(q_point);
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             rhs(local_dof_indices[i]) += cell_rhs(i);
 *         };
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Globalrefinement"></a> 
 * <h4>Global refinement</h4>
 * 

 * 
 * By now, all functions of the abstract base class except for the
 * <code>refine_grid</code> function have been implemented. We will now
 * have two classes that implement this function for the
 * <code>PrimalSolver</code> class, one doing global refinement, one a
 * form of local refinement.
 *     

 * 
 * The first, doing global refinement, is rather simple: its main function
 * just calls <code>triangulation-@>refine_global (1);</code>, which does
 * all the work.
 *     

 * 
 * Note that since the <code>Base</code> base class of the
 * <code>Solver</code> class is virtual, we have to declare a constructor
 * that initializes the immediate base class as well as the abstract
 * virtual one.
 *     

 * 
 * Apart from this technical complication, the class is probably simple
 * enough to be left without further comments.
 * 
 * @code
 *     template <int dim>
 *     class RefinementGlobal : public PrimalSolver<dim>
 *     {
 *     public:
 *       RefinementGlobal(Triangulation<dim> &      coarse_grid,
 *                        const FiniteElement<dim> &fe,
 *                        const Quadrature<dim> &   quadrature,
 *                        const Function<dim> &     rhs_function,
 *                        const Function<dim> &     boundary_values);
 * 
 *       virtual void refine_grid() override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     RefinementGlobal<dim>::RefinementGlobal(
 *       Triangulation<dim> &      coarse_grid,
 *       const FiniteElement<dim> &fe,
 *       const Quadrature<dim> &   quadrature,
 *       const Function<dim> &     rhs_function,
 *       const Function<dim> &     boundary_values)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           fe,
 *                           quadrature,
 *                           rhs_function,
 *                           boundary_values)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void RefinementGlobal<dim>::refine_grid()
 *     {
 *       this->triangulation->refine_global(1);
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LocalrefinementbytheKellyerrorindicator"></a> 
 * <h4>Local refinement by the Kelly error indicator</h4>
 * 

 * 
 * The second class implementing refinement strategies uses the Kelly
 * refinement indicator used in various example programs before. Since this
 * indicator is already implemented in a class of its own inside the
 * deal.II library, there is not much t do here except cal the function
 * computing the indicator, then using it to select a number of cells for
 * refinement and coarsening, and refinement the mesh accordingly.
 *     

 * 
 * Again, this should now be sufficiently standard to allow the omission
 * of further comments.
 * 
 * @code
 *     template <int dim>
 *     class RefinementKelly : public PrimalSolver<dim>
 *     {
 *     public:
 *       RefinementKelly(Triangulation<dim> &      coarse_grid,
 *                       const FiniteElement<dim> &fe,
 *                       const Quadrature<dim> &   quadrature,
 *                       const Function<dim> &     rhs_function,
 *                       const Function<dim> &     boundary_values);
 * 
 *       virtual void refine_grid() override;
 *     };
 * 
 * 
 * 
 *     template <int dim>
 *     RefinementKelly<dim>::RefinementKelly(Triangulation<dim> &      coarse_grid,
 *                                           const FiniteElement<dim> &fe,
 *                                           const Quadrature<dim> &   quadrature,
 *                                           const Function<dim> &rhs_function,
 *                                           const Function<dim> &boundary_values)
 *       : Base<dim>(coarse_grid)
 *       , PrimalSolver<dim>(coarse_grid,
 *                           fe,
 *                           quadrature,
 *                           rhs_function,
 *                           boundary_values)
 *     {}
 * 
 * 
 * 
 *     template <int dim>
 *     void RefinementKelly<dim>::refine_grid()
 *     {
 *       Vector<float> estimated_error_per_cell(
 *         this->triangulation->n_active_cells());
 *       KellyErrorEstimator<dim>::estimate(
 *         this->dof_handler,
 *         QGauss<dim - 1>(this->fe->degree + 1),
 *         std::map<types::boundary_id, const Function<dim> *>(),
 *         this->solution,
 *         estimated_error_per_cell);
 *       GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation,
 *                                                       estimated_error_per_cell,
 *                                                       0.3,
 *                                                       0.03);
 *       this->triangulation->execute_coarsening_and_refinement();
 *     }
 * 
 *   } // namespace LaplaceSolver
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
 * As this is one more academic example, we'd like to compare exact and
 * computed solution against each other. For this, we need to declare
 * function classes representing the exact solution (for comparison and for
 * the Dirichlet boundary values), as well as a class that denotes the right
 * hand side of the equation (this is simply the Laplace operator applied to
 * the exact solution we'd like to recover).
 *   

 * 
 * For this example, let us choose as exact solution the function
 * $u(x,y)=exp(x+sin(10y+5x^2))$. In more than two dimensions, simply repeat
 * the sine-factor with <code>y</code> replaced by <code>z</code> and so
 * on. Given this, the following two classes are probably straightforward
 * from the previous examples.
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double Solution<dim>::value(const Point<dim> & p,
 *                               const unsigned int component) const
 *   {
 *     (void)component;
 *     AssertIndexRange(component, 1);
 *     double q = p(0);
 *     for (unsigned int i = 1; i < dim; ++i)
 *       q += std::sin(10 * p(i) + 5 * p(0) * p(0));
 *     const double exponential = std::exp(q);
 *     return exponential;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 *     (void)component;
 *     AssertIndexRange(component, 1);
 *     double q = p(0);
 *     for (unsigned int i = 1; i < dim; ++i)
 *       q += std::sin(10 * p(i) + 5 * p(0) * p(0));
 *     const double u  = std::exp(q);
 *     double       t1 = 1, t2 = 0, t3 = 0;
 *     for (unsigned int i = 1; i < dim; ++i)
 *       {
 *         t1 += std::cos(10 * p(i) + 5 * p(0) * p(0)) * 10 * p(0);
 *         t2 += 10 * std::cos(10 * p(i) + 5 * p(0) * p(0)) -
 *               100 * std::sin(10 * p(i) + 5 * p(0) * p(0)) * p(0) * p(0);
 *         t3 += 100 * std::cos(10 * p(i) + 5 * p(0) * p(0)) *
 *                 std::cos(10 * p(i) + 5 * p(0) * p(0)) -
 *               100 * std::sin(10 * p(i) + 5 * p(0) * p(0));
 *       };
 *     t1 = t1 * t1;
 * 
 *     return -u * (t1 + t2 + t3);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Thedriverroutines"></a> 
 * <h3>The driver routines</h3>
 * 

 * 
 * What is now missing are only the functions that actually select the
 * various options, and run the simulation on successively finer grids to
 * monitor the progress as the mesh is refined.
 *   

 * 
 * This we do in the following function: it takes a solver object, and a
 * list of postprocessing (evaluation) objects, and runs them with
 * intermittent mesh refinement:
 * 
 * @code
 *   template <int dim>
 *   void run_simulation(
 *     LaplaceSolver::Base<dim> &                          solver,
 *     const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list)
 *   {
 * @endcode
 * 
 * We will give an indicator of the step we are presently computing, in
 * order to keep the user informed that something is still happening, and
 * that the program is not in an endless loop. This is the head of this
 * status line:
 * 
 * @code
 *     std::cout << "Refinement cycle: ";
 * 
 * @endcode
 * 
 * Then start a loop which only terminates once the number of degrees of
 * freedom is larger than 20,000 (you may of course change this limit, if
 * you need more -- or less -- accuracy from your program).
 * 
 * @code
 *     for (unsigned int step = 0; true; ++step)
 *       {
 * @endcode
 * 
 * Then give the <code>alive</code> indication for this
 * iteration. Note that the <code>std::flush</code> is needed to have
 * the text actually appear on the screen, rather than only in some
 * buffer that is only flushed the next time we issue an end-line.
 * 
 * @code
 *         std::cout << step << " " << std::flush;
 * 
 * @endcode
 * 
 * Now solve the problem on the present grid, and run the evaluators
 * on it. The long type name of iterators into the list is a little
 * annoying, but could be shortened by an alias, if so desired.
 * 
 * @code
 *         solver.solve_problem();
 * 
 *         for (const auto &postprocessor : postprocessor_list)
 *           {
 *             postprocessor->set_refinement_cycle(step);
 *             solver.postprocess(*postprocessor);
 *           };
 * 
 * 
 * @endcode
 * 
 * Now check whether more iterations are required, or whether the loop
 * shall be ended:
 * 
 * @code
 *         if (solver.n_dofs() < 20000)
 *           solver.refine_grid();
 *         else
 *           break;
 *       };
 * 
 * @endcode
 * 
 * Finally end the line in which we displayed status reports:
 * 
 * @code
 *     std::cout << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The final function is one which takes the name of a solver (presently
 * "kelly" and "global" are allowed), creates a solver object out of it
 * using a coarse grid (in this case the ubiquitous unit square) and a
 * finite element object (here the likewise ubiquitous bilinear one), and
 * uses that solver to ask for the solution of the problem on a sequence of
 * successively refined grids.
 *   

 * 
 * The function also sets up two of evaluation functions, one evaluating the
 * solution at the point (0.5,0.5), the other writing out the solution to a
 * file.
 * 
 * @code
 *   template <int dim>
 *   void solve_problem(const std::string &solver_name)
 *   {
 * @endcode
 * 
 * First minor task: tell the user what is going to happen. Thus write a
 * header line, and a line with all '-' characters of the same length as
 * the first one right below.
 * 
 * @code
 *     const std::string header =
 *       "Running tests with \"" + solver_name + "\" refinement criterion:";
 *     std::cout << header << std::endl
 *               << std::string(header.size(), '-') << std::endl;
 * 
 * @endcode
 * 
 * Then set up triangulation, finite element, etc.
 * 
 * @code
 *     Triangulation<dim> triangulation;
 *     GridGenerator::hyper_cube(triangulation, -1, 1);
 *     triangulation.refine_global(2);
 *     const FE_Q<dim>    fe(1);
 *     const QGauss<dim>  quadrature(4);
 *     RightHandSide<dim> rhs_function;
 *     Solution<dim>      boundary_values;
 * 
 * @endcode
 * 
 * Create a solver object of the kind indicated by the argument to this
 * function. If the name is not recognized, throw an exception!
 * The respective solver object is stored in a `std::unique_ptr` to avoid
 * having to delete the pointer after use.
 * 
 * @code
 *     std::unique_ptr<LaplaceSolver::Base<dim>> solver;
 *     if (solver_name == "global")
 *       solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>(
 *         triangulation, fe, quadrature, rhs_function, boundary_values);
 *     else if (solver_name == "kelly")
 *       solver = std::make_unique<LaplaceSolver::RefinementKelly<dim>>(
 *         triangulation, fe, quadrature, rhs_function, boundary_values);
 *     else
 *       AssertThrow(false, ExcNotImplemented());
 * 
 * @endcode
 * 
 * Next create a table object in which the values of the numerical
 * solution at the point (0.5,0.5) will be stored, and create a respective
 * evaluation object:
 * 
 * @code
 *     TableHandler                          results_table;
 *     Evaluation::PointValueEvaluation<dim> postprocessor1(Point<dim>(0.5, 0.5),
 *                                                          results_table);
 * 
 * @endcode
 * 
 * Also generate an evaluator which writes out the solution:
 * 
 * @code
 *     Evaluation::SolutionOutput<dim> postprocessor2(std::string("solution-") +
 *                                                      solver_name,
 *                                                    DataOutBase::gnuplot);
 * 
 * @endcode
 * 
 * Take these two evaluation objects and put them in a list...
 * 
 * @code
 *     std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list;
 *     postprocessor_list.push_back(&postprocessor1);
 *     postprocessor_list.push_back(&postprocessor2);
 * 
 * @endcode
 * 
 * ... which we can then pass on to the function that actually runs the
 * simulation on successively refined grids:
 * 
 * @code
 *     run_simulation(*solver, postprocessor_list);
 * 
 * @endcode
 * 
 * When this all is done, write out the results of the point evaluations:
 * 
 * @code
 *     results_table.write_text(std::cout);
 * 
 * @endcode
 * 
 * And one blank line after all results:
 * 
 * @code
 *     std::cout << std::endl;
 *   }
 * } // namespace Step13
 * 
 * 
 * 
 * @endcode
 * 
 * There is not much to say about the main function. It follows the same
 * pattern as in all previous examples, with attempts to catch thrown
 * exceptions, and displaying as much information as possible if we should get
 * some. The rest is self-explanatory.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step13::solve_problem<2>("global");
 *       Step13::solve_problem<2>("kelly");
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
examples/step-13/doc/results.dox



<a name="Results"></a><h1>Results</h1>





这个程序的结果并不那么有趣--毕竟它的目的不是为了演示一些新的数学思想，也不是为了演示如何用deal.II编程，而是为了利用我们在前面的例子中所开发的材料，形成一些演示以模块化和可扩展的方式建立现代有限元软件的方法。




尽管如此，我们当然会展示程序的结果。最感兴趣的是点值计算，为此我们实现了相应的评估类。该程序的结果（即输出）看起来如下。

@code
Running tests with "global" refinement criterion:


-------------------------------------------------
Refinement cycle: 0 1 2 3 4 5 6
DoFs  u(x_0)
   25 1.2868
   81 1.6945
  289 1.4658
 1089 1.5679
 4225 1.5882
16641 1.5932
66049 1.5945


Running tests with "kelly" refinement criterion:


------------------------------------------------
Refinement cycle: 0 1 2 3 4 5 6 7 8 9 10 11
DoFs  u(x_0)
   25 1.2868
   47 0.8775
   89 1.5365
  165 1.2974
  316 1.6442
  589 1.5221
 1093 1.5724
 2042 1.5627
 3766 1.5916
 7124 1.5876
13111 1.5942
24838 1.5932
@endcode




这里令人惊讶的是，精确的数值是1.59491554...，而且计算该解显然出奇的复杂，甚至只达到百分之一的精度，尽管该解是平滑的（事实上是无限常可微）。这种平滑性显示在程序生成的图形输出中，这里是粗网格和凯利细化指标的前9个细化步骤。


 <table width="80%" align="center">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-0.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-1.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-3.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-4.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-5.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-6.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-7.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-8.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-13.solution-kelly-9.png" alt="">
    </td>
  </tr>
</table> 


当我们已经在观看图片时，这是第八个网格，从顶部看。


 <img src="https://www.dealii.org/images/steps/developer/step-13.grid-kelly-8.png" alt=""> 


然而，我们还没有完成对点值计算的评估。事实上，将两个细化标准的误差 $e=|u(x_0)-u_h(x_0)|$ 绘制成图，可以得到下面的图片。


 <img src="https://www.dealii.org/images/steps/developer/step-13.error.png" alt=""> 





这幅图 <em> 令人不安的是，自适应网格细化不仅没有像人们通常期望的那样比全局细化好，甚至明显更差，因为它的收敛是不规则的，在使用后续网格的值时，阻止了所有的外推技术!另一方面，全局细化提供了一个完美的 $1/N$ 或 $h^{-2}$ 收敛历史，并提供了各种机会，甚至可以通过外推法来改善点值。因此，在这个例子中，全局网格细化必须被认为是优越的!这更令人惊讶，因为评估点不是在左边的某个地方，那里的网格是粗糙的，而是在右边，自适应细化也应该细化评估点周围的网格。




因此，我们以一个问题来结束对这个例子程序的讨论。

<p align="center"> <strong>  <em>  如果适应性不比全局细化好，那么它有什么问题？ </em>  </strong>





 <em>  在这个例子的结尾处进行练习。 </em>  有一个简单的原因导致适应性网格解决方案的不良和不规则行为。通过观察每个步骤中评估点周围的网格，可以很简单地找到这个原因--这个数据在程序的输出文件中。因此，一个练习是修改网格细化程序，使问题（一旦你注意到它）得以避免。第二个练习是检查结果是否比全局细化要好，如果是的话，是否能达到更好的收敛顺序（就自由度数而言），或者只达到一个更好的常数。




(  <em>  对于没有耐心的人来说，非常简短的回答。 </em> 在误差较大的步骤中，网格在评估点上是不规则的，即一些相邻的单元有悬空的节点；这破坏了一些超级近似的效果，而全局细化的网格可以从中受益。答案2：这个快速黑客

@code
  bool refinement_indicated = false;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto v : cell->vertex_indices())
	  if (cell->vertex(v) == Point<dim>(.5,.5))
	    {
	      cell->clear_coarsen_flag();
	      refinement_indicated |= cell->refine_flag_set();
	    }
  if (refinement_indicated)
    for (const auto &cell : triangulation.active_cell_iterators())
      for (const auto v : cell->vertex_indices())
	    if (cell->vertex(v) == Point<dim>(.5,.5))
	      cell->set_refine_flag ();
@endcode

在执行细化之前，在Kelly细化类的细化函数中，将改善结果（练习：代码是怎么做的？不过，行为仍然是不规则的，所以不可能有关于收敛顺序的结果）。)


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-13.cc"
*/
