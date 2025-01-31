/**
  @page step_29 The step-29 tutorial program  
* This tutorial depends on   step-4  .
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemsetting">Problem setting</a><a href="#Problemsetting">Problem setting</a>
        <li><a href="#Thetestcase">The test case</a><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeDirichletBoundaryValuescodeclass">The <code>DirichletBoundaryValues</code> class</a><a href="#ThecodeDirichletBoundaryValuescodeclass">The <code>DirichletBoundaryValues</code> class</a>
        <li><a href="#ThecodeParameterReadercodeclass">The <code>ParameterReader</code> class</a><a href="#ThecodeParameterReadercodeclass">The <code>ParameterReader</code> class</a>
      <ul>
        <li><a href="#codeParameterReaderdeclare_parameterscode"><code>ParameterReader::declare_parameters</code></a><a href="#codeParameterReaderdeclare_parameterscode"><code>ParameterReader::declare_parameters</code></a>
        <li><a href="#codeParameterReaderread_parameterscode"><code>ParameterReader::read_parameters</code></a><a href="#codeParameterReaderread_parameterscode"><code>ParameterReader::read_parameters</code></a>
      </ul>
        <li><a href="#ThecodeComputeIntensitycodeclass">The <code>ComputeIntensity</code> class</a><a href="#ThecodeComputeIntensitycodeclass">The <code>ComputeIntensity</code> class</a>
        <li><a href="#ThecodeUltrasoundProblemcodeclass">The <code>UltrasoundProblem</code> class</a><a href="#ThecodeUltrasoundProblemcodeclass">The <code>UltrasoundProblem</code> class</a>
      <ul>
        <li><a href="#codeUltrasoundProblemmake_gridcode"><code>UltrasoundProblem::make_grid</code></a><a href="#codeUltrasoundProblemmake_gridcode"><code>UltrasoundProblem::make_grid</code></a>
        <li><a href="#codeUltrasoundProblemsetup_systemcode"><code>UltrasoundProblem::setup_system</code></a><a href="#codeUltrasoundProblemsetup_systemcode"><code>UltrasoundProblem::setup_system</code></a>
        <li><a href="#codeUltrasoundProblemassemble_systemcode"><code>UltrasoundProblem::assemble_system</code></a><a href="#codeUltrasoundProblemassemble_systemcode"><code>UltrasoundProblem::assemble_system</code></a>
        <li><a href="#codeUltrasoundProblemsolvecode"><code>UltrasoundProblem::solve</code></a><a href="#codeUltrasoundProblemsolvecode"><code>UltrasoundProblem::solve</code></a>
        <li><a href="#codeUltrasoundProblemoutput_resultscode"><code>UltrasoundProblem::output_results</code></a><a href="#codeUltrasoundProblemoutput_resultscode"><code>UltrasoundProblem::output_results</code></a>
        <li><a href="#codeUltrasoundProblemruncode"><code>UltrasoundProblem::run</code></a><a href="#codeUltrasoundProblemruncode"><code>UltrasoundProblem::run</code></a>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
*   <br>  
* <i>
This program was contributed by Moritz Allmaras at Texas A&amp;M
University. Some of the work on this tutorial program has been funded
by NSF under grant DMS-0604778.
</i>
* <b>Note:</b> In order to run this program, deal.II must be configured to usethe UMFPACK sparse direct solver. Refer to the <a
href="../../readme.html#umfpack">ReadMe</a> for instructions how to do this.
* 

* <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1>
* 

* 
* A question that comes up frequently is how to solve problems involving complexvalued functions with deal.II. For many problems, instead of working withcomplex valued finite elements directly, it is often more convenient to split complex valuedfunctions into their real and imaginary parts and use separate scalar finiteelement fields for discretizing each one of them. Basically this amounts toviewing a single complex valued equation as a system of two real valuedequations. This short example demonstrates how this can be implemented indeal.II by using an   <code>FE_system</code>   object to stack two finite elementfields representing real and imaginary parts. (The opposite approach,keeping everything complex-valued, is demonstrated in a differenttutorial program: see   step-58   for this.)When split into real and imaginary parts, the equations covered herefall into the class of vector-valued problems. A toplevel overview ofthis topic can be found in the   @ref vector_valued   module.
* In addition to this discussion, we also discuss the ParameterHandlerclass, which provides a convenient way for reading parameters from aconfiguration file at runtime without the need to recompile theprogram code.
* 

* <a name="Problemsetting"></a><h3>Problem setting</h3>
* 

* The original purpose of this program is to simulate the focusing propertiesof an ultrasound wave generated by a transducer lens with variablegeometry. Recent applications in medical imaging use ultrasound waves not onlyfor imaging purposes, but also to excite certain local effects in amaterial, like changes in optical properties, that can then be measured byother imaging techniques. A vital ingredient for these methods is the abilityto focus the intensity of the ultrasound wave in a particular part of thematerial, ideally in a point, to be able to examine the properties of thematerial at that particular location.
* To derive a model for this problem, we think of ultrasound as a pressure wavegoverned by the wave equation:@f[
	\frac{\partial^2 U}{\partial t^2}	-	c^2 \Delta U = 0
@f]where   $c$   is the wave speed (that for simplicity we assume to be constant),   $U
= U(x,t),\;x \in \Omega,\;t\in\mathrm{R}$  . The boundary  $\Gamma=\partial\Omega$   is divided into two parts   $\Gamma_1$   and  $\Gamma_2=\Gamma\setminus\Gamma_1$  , with   $\Gamma_1$   representing thetransducer lens and   $\Gamma_2$   an absorbing boundary (that is, we want tochoose boundary conditions on   $\Gamma_2$   in such a way that they imitate alarger domain). On   $\Gamma_1$  , the transducer generates a wave of constantfrequency   ${\omega}>0$   and constant amplitude (that we chose to be 1 here):@f[
U(x,t) = \cos{\omega t}, \qquad x\in \Gamma_1
@f]
* If there are no other (interior or boundary) sources, and since the onlysource has frequency   $\omega$  , then the solution admits a separation ofvariables of the form   $U(x,t) = \textrm{Re}\left(u(x)\,e^{i\omega
t})\right)$  . The complex-valued function   $u(x)$   describes the spatialdependency of amplitude and phase (relative to the source) of the waves offrequency   ${\omega}$  , with the amplitude being the quantity that we areinterested in. By plugging this form of the solution into the wave equation,we see that for   $u$   we have@f{eqnarray*}
* 
-\omega^2 u(x)
* 
- c^2\Delta u(x) &=& 0, \qquad x\in\Omega,\\
u(x) &=& 1,  \qquad x\in\Gamma_1.
@f}
* 
* For finding suitable conditions on   $\Gamma_2$   that model an absorbingboundary, consider a wave of the form   $V(x,t)=e^{i(k\cdot x
* 
-\omega t)}$   withfrequency   ${\omega}$   traveling in direction   $k\in {\mathrm{R}^2}$  . In orderfor   $V$   to solve the wave equation,   $|k|={\frac{\omega}{c}}$   musthold. Suppose that this wave hits the boundary in   $x_0\in\Gamma_2$   at a rightangle, i.e.   $n=\frac{k}{|k|}$   with   $n$   denoting the outer unit normal of  $\Omega$   in   $x_0$  . Then at   $x_0$  , this wave satisfies the equation@f[
c (n\cdot\nabla V) + \frac{\partial V}{\partial t} = (i\, c\, |k|
* 
- i\, \omega) V = 0.
@f]Hence, by enforcing the boundary condition@f[
c (n\cdot\nabla U) + \frac{\partial U}{\partial t} = 0, \qquad x\in\Gamma_2,
@f]waves that hit the boundary   $\Gamma_2$   at a right angle will be perfectlyabsorbed. On the other hand, those parts of the wave field that do not hit aboundary at a right angle do not satisfy this condition and enforcing it as aboundary condition will yield partial reflections, i.e. only parts of the wavewill pass through the boundary as if it wasn't here whereas the remainingfraction of the wave will be reflected back into the domain.
* If we are willing to accept this as a sufficient approximation to an absorbingboundary we finally arrive at the following problem for   $u$  :@f{eqnarray*}
* 
-\omega^2 u
* 
- c^2\Delta u &=& 0, \qquad x\in\Omega,\\
c (n\cdot\nabla u) + i\,\omega\,u &=&0, \qquad x\in\Gamma_2,\\
u &=& 1,  \qquad x\in\Gamma_1.
@f}
* This is a Helmholtz equation (similar to the one in   step-7  , but this time with''the bad sign'') with Dirichlet data on   $\Gamma_1$   and mixed boundaryconditions on   $\Gamma_2$  . Because of the condition on   $\Gamma_2$  , we cannot justtreat the equations for real and imaginary parts of   $u$   separately. What we cando however is to view the PDE for   $u$   as a system of two PDEs for the real andimaginary parts of   $u$  , with the boundary condition on   $\Gamma_2$   representingthe coupling terms between the two components of the system. This works alongthe following lines: Let   $v=\textrm{Re}\;u,\; w=\textrm{Im}\;u$  , then in termsof   $v$   and   $w$   we have the following system:@f{eqnarray*}
  \left.\begin{array}{ccc}
* 


* 
* 


* 
* 


* 
* 


* 
* 
-\omega^2 v
* 
- c^2\Delta v &=& 0 \quad\\
* 


* 
* 


* 
* 


* 
* 


* 
* 
-\omega^2 w
* 
- c^2\Delta w &=& 0 \quad
  \end{array}\right\} &\;& x\in\Omega,
	\\
  \left.\begin{array}{ccc}
    c (n\cdot\nabla v)
* 
- \omega\,w &=& 0 \quad\\
    c (n\cdot\nabla w) + \omega\,v &=& 0 \quad
  \end{array}\right\} &\;& x\in\Gamma_2,
	\\
	\left.\begin{array}{ccc}
    v &=& 1 \quad\\
    w &=& 0 \quad
  \end{array}\right\} &\;& x\in\Gamma_1.
@f}
* 
* For test functions   $\phi,\psi$   with   $\phi|_{\Gamma_1}=\psi|_{\Gamma_1}=0$  , afterthe usual multiplication, integration over   $\Omega$   and applying integration byparts, we get the weak formulation@f{eqnarray*}
* 
-\omega^2 \langle \phi, v \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi, \nabla v \rangle_{\mathrm{L}^2(\Omega)}
* 
- c \omega \langle \phi, w \rangle_{\mathrm{L}^2(\Gamma_2)} &=& 0, \\
* 
-\omega^2 \langle \psi, w \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi, \nabla w \rangle_{\mathrm{L}^2(\Omega)}
+ c \omega \langle \psi, v \rangle_{\mathrm{L}^2(\Gamma_2)} &=& 0.
@f}
* 
* We choose finite element spaces   $V_h$   and   $W_h$   with bases   $\{\phi_j\}_{j=1}^n,
\{\psi_j\}_{j=1}^n$   and look for approximate solutions@f[
v_h = \sum_{j=1}^n \alpha_j \phi_j, \;\; w_h = \sum_{j=1}^n \beta_j \psi_j.
@f]Plugging into the variational form yields the equation system@f[
\renewcommand{\arraystretch}{2.0}
\left.\begin{array}{ccc}
\sum_{j=1}^n
\left(
* 
-\omega^2 \langle \phi_i, \phi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi_i, \nabla \phi_j \rangle_{\mathrm{L}^2(\Omega)}
\right)
\alpha_j
* 
- \left(
c\omega \langle \phi_i,\psi_j\rangle_{\mathrm{L}^2(\Gamma_2)}\right)\beta_j
&=& 0 \\
\sum_{j=1}^n
\left(
* 
-\omega^2 \langle \psi_i, \psi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi_i, \nabla \psi_j \rangle_{\mathrm{L}^2(\Omega)}
\right)\beta_j
+ \left(
c\omega \langle
\psi_i,\phi_j\rangle_{\mathrm{L}^2(\Gamma_2)}
\right)\alpha_j
&=& 0
\end{array}\right\}\;\;\forall\; i =1,\ldots,n.
@f]In matrix notation:@f[
\renewcommand{\arraystretch}{2.0}
\left(
\begin{array}{cc}
* 
-\omega^2 \langle \phi_i, \phi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi_i, \nabla \phi_j \rangle_{\mathrm{L}^2(\Omega)}
&
* 
-c\omega \langle \phi_i,\psi_j\rangle_{\mathrm{L}^2(\Gamma_2)} \\
c\omega \langle \psi_i,\phi_j\rangle_{\mathrm{L}^2(\Gamma_2)}
&
* 
-\omega^2 \langle \psi_{i}, \psi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi_{i}, \nabla \psi_j  \rangle_{\mathrm{L}^2(\Omega)}
\end{array}
\right)
\left(
\begin{array}{c}
\alpha \\ \beta
\end{array}
\right)
=
\left(
\begin{array}{c}
0 \\ 0
\end{array}
\right)
@f](One should not be fooled by the right hand side being zero here, that isbecause we haven't included the Dirichlet boundary data yet.)Because of the alternating sign in the off-diagonal blocks, we can alreadysee that this system is non-symmetric, in fact it is even indefinite.Of course, there is no necessity to choose the spaces   $V_h$   and   $W_h$   to bethe same. However, we expect real and imaginary part of the solution tohave similar properties and will therefore indeed take   $V_h=W_h$   in theimplementation, and also use the same basis functions   $\phi_i = \psi_i$   forboth spaces. The reason for the notation using different symbols is just thatit allows us to distinguish between shape functions for   $v$   and   $w$  , as thisdistinction plays an important role in the implementation.
* 

* <a name="Thetestcase"></a><h3>The test case</h3>
* 

* For the computations, we will consider wave propagation in the unit square,with ultrasound generated by a transducer lens that is shaped like a segmentof the circle with center at   $(0.5, d)$   and aradius slightly greater than   $d$  ; this shape should lead to a focusing of the soundwave at the center of the circle. Varying   $d$   changes the "focus" of the lensand affects the spatial distribution of the intensity of   $u$  , where our mainconcern is how well   $|u|=\sqrt{v^2+w^2}$   is focused.
* In the program below, we will implement the complex-valued Helmholtz equationsusing the formulation with split real and imaginary parts. We will alsodiscuss how to generate a domain that looks like a square with a slight bulgesimulating the transducer (in the  <code>UltrasoundProblem<dim>::make_grid()</code>   function), and how togenerate graphical output that not only contains the solution components   $v$   and  $w$  , but also the magnitude   $\sqrt{v^2+w^2}$   directly in the output file (in  <code>UltrasoundProblem<dim>::output_results()</code>  ). Finally, we use theParameterHandler class to easily read parameters like the focal distance   $d$  ,wave speed   $c$  , frequency   $\omega$  , and a number of other parameters from aninput file at run-time, rather than fixing those parameters in the source codewhere we would have to re-compile every time we want to change parameters.
* 

*  <a name="CommProg"></a> <h1> The commented program</h1>
*   <a name="Includefiles"></a>  <h3>Include files</h3>
* 

* 
*  The following header files have all been discussed before:
* 

* 
*  

* 
* @code
 #include <deal.II/base/quadrature_lib.h>
 #include <deal.II/base/function.h>
 #include <deal.II/base/logstream.h>
* 
 #include <deal.II/lac/vector.h>
 #include <deal.II/lac/full_matrix.h>
 #include <deal.II/lac/sparse_matrix.h>
 #include <deal.II/lac/dynamic_sparsity_pattern.h>
* 
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/manifold_lib.h>
* 
 #include <deal.II/dofs/dof_handler.h>
 #include <deal.II/dofs/dof_tools.h>
* 
 #include <deal.II/fe/fe_q.h>
 #include <deal.II/fe/fe_values.h>
* 
 #include <deal.II/numerics/matrix_tools.h>
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/numerics/vector_tools.h>
* 
 #include <iostream>
 #include <fstream>
* 
 @endcode
* 
*  This header file contains the necessary declarations for the ParameterHandler class that we will use to read our parameters from a configuration file:
* 

* 
* @code
 #include <deal.II/base/parameter_handler.h>
* 
 @endcode
* 
*  For solving the linear system, we'll use the sparse LU-decomposition provided by UMFPACK (see the SparseDirectUMFPACK class), for which the following header file is needed.  Note that in order to compile this tutorial program, the deal.II-library needs to be built with UMFPACK support, which is enabled by default:
* 

* 
* @code
 #include <deal.II/lac/sparse_direct.h>
* 
 @endcode
* 
*  The FESystem class allows us to stack several FE-objects to one compound, vector-valued finite element field. The necessary declarations for this class are provided in this header file:
* 

* 
* @code
 #include <deal.II/fe/fe_system.h>
* 
 @endcode
* 
*  Finally, include the header file that declares the Timer class that we will use to determine how much time each of the operations of our program takes:
* 

* 
* @code
 #include <deal.II/base/timer.h>
* 
 @endcode
* 
*  As the last step at the beginning of this program, we put everything that is in this program into its namespace and, within it, make everything that is in the deal.II namespace globally available, without the need to prefix everything with   <code>dealii</code><code>::</code>  :
* 

* 
* @code
 namespace Step29
 {
   using namespace dealii;
* 
 
 @endcode
* 
*   <a name="ThecodeDirichletBoundaryValuescodeclass"></a>  <h3>The <code>DirichletBoundaryValues</code> class</h3>
* 

* 
*  First we define a class for the function representing the Dirichlet boundary values. This has been done many times before and therefore does not need much explanation.   
*   Since there are two values   $v$   and   $w$   that need to be prescribed at the boundary, we have to tell the base class that this is a vector-valued function with two components, and the   <code>vector_value</code>   function and its cousin   <code>vector_value_list</code>   must return vectors with two entries. In our case the function is very simple, it just returns 1 for the real part   $v$   and 0 for the imaginary part   $w$   regardless of the point where it is evaluated.
* 

* 
* @code
   template <int dim>
   class DirichletBoundaryValues : public Function<dim>
   {
   public:
     DirichletBoundaryValues()
       : Function<dim>(2)
     {}
* 
     virtual void vector_value(const Point<dim> &  /*p*/ ,
                               Vector<double> &values) const override
     {
       Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2));
* 
       values(0) = 1;
       values(1) = 0;
     }
* 
     virtual void
     vector_value_list(const std::vector<Point<dim>> &points,
                       std::vector<Vector<double>> &  value_list) const override
     {
       Assert(value_list.size() == points.size(),
              ExcDimensionMismatch(value_list.size(), points.size()));
* 
       for (unsigned int p = 0; p < points.size(); ++p)
         DirichletBoundaryValues<dim>::vector_value(points[p], value_list[p]);
     }
   };
* 
 @endcode
* 
*   <a name="ThecodeParameterReadercodeclass"></a>  <h3>The <code>ParameterReader</code> class</h3>
* 

* 
*  The next class is responsible for preparing the ParameterHandler object and reading parameters from an input file.  It includes a function   <code>declare_parameters</code>   that declares all the necessary parameters and a   <code>read_parameters</code>   function that is called from outside to initiate the parameter reading process.
* 

* 
* @code
   class ParameterReader : public Subscriptor
   {
   public:
     ParameterReader(ParameterHandler &);
     void read_parameters(const std::string &);
* 
   private:
     void              declare_parameters();
     ParameterHandler &prm;
   };
* 
 @endcode
* 
*  The constructor stores a reference to the ParameterHandler object that is passed to it:
* 

* 
* @code
   ParameterReader::ParameterReader(ParameterHandler &paramhandler)
     : prm(paramhandler)
   {}
* 
 @endcode
* 
*   <a name="codeParameterReaderdeclare_parameterscode"></a>  <h4><code>ParameterReader::declare_parameters</code></h4>
* 

* 
*  The   <code>declare_parameters</code>   function declares all the parameters that our ParameterHandler object will be able to read from input files, along with their types, range conditions and the subsections they appear in. We will wrap all the entries that go into a section in a pair of braces to force the editor to indent them by one level, making it simpler to read which entries together form a section:
* 

* 
* @code
   void ParameterReader::declare_parameters()
   {
 @endcode
* 
*  Parameters for mesh and geometry include the number of global refinement steps that are applied to the initial coarse mesh and the focal distance   $d$   of the transducer lens. For the number of refinement steps, we allow integer values in the range   $[0,\infty)$  , where the omitted second argument to the   Patterns::Integer   object denotes the half-open interval.  For the focal distance any number greater than zero is accepted:
* 

* 
* @code
     prm.enter_subsection("Mesh & geometry parameters");
     {
       prm.declare_entry("Number of refinements",
                         "6",
                         Patterns::Integer(0),
                         "Number of global mesh refinement steps "
                         "applied to initial coarse grid");
* 
       prm.declare_entry("Focal distance",
                         "0.3",
                         Patterns::Double(0),
                         "Distance of the focal point of the lens "
                         "to the x-axis");
     }
     prm.leave_subsection();
* 
 @endcode
* 
*  The next subsection is devoted to the physical parameters appearing in the equation, which are the frequency   $\omega$   and wave speed   $c$  . Again, both need to lie in the half-open interval   $[0,\infty)$   represented by calling the   Patterns::Double   class with only the left end-point as argument:
* 

* 
* @code
     prm.enter_subsection("Physical constants");
     {
       prm.declare_entry("c", "1.5e5", Patterns::Double(0), "Wave speed");
* 
       prm.declare_entry("omega", "5.0e7", Patterns::Double(0), "Frequency");
     }
     prm.leave_subsection();
* 
 
 @endcode
* 
*  Last but not least we would like to be able to change some properties of the output, like filename and format, through entries in the configuration file, which is the purpose of the last subsection:
* 

* 
* @code
     prm.enter_subsection("Output parameters");
     {
       prm.declare_entry("Output filename",
                         "solution",
                         Patterns::Anything(),
                         "Name of the output file (without extension)");
* 
 @endcode
* 
*  Since different output formats may require different parameters for generating output (like for example, postscript output needs viewpoint angles, line widths, colors etc), it would be cumbersome if we had to declare all these parameters by hand for every possible output format supported in the library. Instead, each output format has a   <code>FormatFlags::declare_parameters</code>   function, which declares all the parameters specific to that format in an own subsection. The following call of   DataOutInterface<1>::declare_parameters   executes   <code>declare_parameters</code>   for all available output formats, so that for each format an own subsection will be created with parameters declared for that particular output format. (The actual value of the template parameter in the call,   <code>@<1@></code>   above, does not matter here: the function does the same work independent of the dimension, but happens to be in a template-parameter-dependent class.)  To find out what parameters there are for which output format, you can either consult the documentation of the DataOutBase class, or simply run this program without a parameter file present. It will then create a file with all declared parameters set to their default values, which can conveniently serve as a starting point for setting the parameters to the values you desire.
* 

* 
* @code
       DataOutInterface<1>::declare_parameters(prm);
     }
     prm.leave_subsection();
   }
* 
 @endcode
* 
*   <a name="codeParameterReaderread_parameterscode"></a>  <h4><code>ParameterReader::read_parameters</code></h4>
* 

* 
*  This is the main function in the ParameterReader class.  It gets called from outside, first declares all the parameters, and then reads them from the input file whose filename is provided by the caller. After the call to this function is complete, the   <code>prm</code>   object can be used to retrieve the values of the parameters read in from the file:
* 

* 
* @code
   void ParameterReader::read_parameters(const std::string &parameter_file)
   {
     declare_parameters();
* 
     prm.parse_input(parameter_file);
   }
* 
 
* 
 @endcode
* 
*   <a name="ThecodeComputeIntensitycodeclass"></a>  <h3>The <code>ComputeIntensity</code> class</h3>
* 

* 
*  As mentioned in the introduction, the quantity that we are really after is the spatial distribution of the intensity of the ultrasound wave, which corresponds to   $|u|=\sqrt{v^2+w^2}$  . Now we could just be content with having   $v$   and   $w$   in our output, and use a suitable visualization or postprocessing tool to derive   $|u|$   from the solution we computed. However, there is also a way to output data derived from the solution in deal.II, and we are going to make use of this mechanism here.
* 

* 
*  So far we have always used the   DataOut::add_data_vector   function to add vectors containing output data to a DataOut object.  There is a special version of this function that in addition to the data vector has an additional argument of type DataPostprocessor. What happens when this function is used for output is that at each point where output data is to be generated, the   DataPostprocessor::evaluate_scalar_field()   or   DataPostprocessor::evaluate_vector_field()   function of the specified DataPostprocessor object is invoked to compute the output quantities from the values, the gradients and the second derivatives of the finite element function represented by the data vector (in the case of face related data, normal vectors are available as well). Hence, this allows us to output any quantity that can locally be derived from the values of the solution and its derivatives.  Of course, the ultrasound intensity   $|u|$   is such a quantity and its computation doesn't even involve any derivatives of   $v$   or   $w$  .
* 

* 
*  In practice, the DataPostprocessor class only provides an interface to this functionality, and we need to derive our own class from it in order to implement the functions specified by the interface. In the most general case one has to implement several member functions but if the output quantity is a single scalar then some of this boilerplate code can be handled by a more specialized class, DataPostprocessorScalar and we can derive from that one instead. This is what the   <code>ComputeIntensity</code>   class does:
* 

* 
* @code
   template <int dim>
   class ComputeIntensity : public DataPostprocessorScalar<dim>
   {
   public:
     ComputeIntensity();
* 
     virtual void evaluate_vector_field(
       const DataPostprocessorInputs::Vector<dim> &inputs,
       std::vector<Vector<double>> &computed_quantities) const override;
   };
* 
 @endcode
* 
*  In the constructor, we need to call the constructor of the base class with two arguments. The first denotes the name by which the single scalar quantity computed by this class should be represented in output files. In our case, the postprocessor has   $|u|$   as output, so we use "Intensity".   
*   The second argument is a set of flags that indicate which data is needed by the postprocessor in order to compute the output quantities.  This can be any subset of update_values, update_gradients and update_hessians (and, in the case of face data, also update_normal_vectors), which are documented in UpdateFlags.  Of course, computation of the derivatives requires additional resources, so only the flags for data that are really needed should be given here, just as we do when we use FEValues objects. In our case, only the function values of   $v$   and   $w$   are needed to compute   $|u|$  , so we're good with the update_values flag.
* 

* 
* @code
   template <int dim>
   ComputeIntensity<dim>::ComputeIntensity()
     : DataPostprocessorScalar<dim>("Intensity", update_values)
   {}
* 
 
 @endcode
* 
*  The actual postprocessing happens in the following function. Its input is an object that stores values of the function (which is here vector-valued) representing the data vector given to   DataOut::add_data_vector,   evaluated at all evaluation points where we generate output, and some tensor objects representing derivatives (that we don't use here since   $|u|$   is computed from just   $v$   and   $w$  ). The derived quantities are returned in the   <code>computed_quantities</code>   vector. Remember that this function may only use data for which the respective update flag is specified by   <code>get_needed_update_flags</code>  . For example, we may not use the derivatives here, since our implementation of   <code>get_needed_update_flags</code>   requests that only function values are provided.
* 

* 
* @code
   template <int dim>
   void ComputeIntensity<dim>::evaluate_vector_field(
     const DataPostprocessorInputs::Vector<dim> &inputs,
     std::vector<Vector<double>> &               computed_quantities) const
   {
     Assert(computed_quantities.size() == inputs.solution_values.size(),
            ExcDimensionMismatch(computed_quantities.size(),
                                 inputs.solution_values.size()));
* 
 @endcode
* 
*  The computation itself is straightforward: We iterate over each entry in the output vector and compute   $|u|$   from the corresponding values of   $v$   and   $w$  . We do this by creating a complex number   $u$   and then calling   `std::abs()`   on the result. (One may be tempted to call   `std::norm()`,   but in a historical quirk, the C++ committee decided that   `std::norm()`   should return the <i>square</i> of the absolute value
* 
*  -  thereby not satisfying the properties mathematicians require of something called a "norm".)
* 

* 
* @code
     for (unsigned int i = 0; i < computed_quantities.size(); i++)
       {
         Assert(computed_quantities[i].size() == 1,
                ExcDimensionMismatch(computed_quantities[i].size(), 1));
         Assert(inputs.solution_values[i].size() == 2,
                ExcDimensionMismatch(inputs.solution_values[i].size(), 2));
* 
         const std::complex<double> u(inputs.solution_values[i](0),
                                      inputs.solution_values[i](1));
* 
         computed_quantities[i](0) = std::abs(u);
       }
   }
* 
 
 @endcode
* 
*   <a name="ThecodeUltrasoundProblemcodeclass"></a>  <h3>The <code>UltrasoundProblem</code> class</h3>
* 

* 
*  Finally here is the main class of this program.  It's member functions are very similar to the previous examples, in particular   step-4  , and the list of member variables does not contain any major surprises either. The ParameterHandler object that is passed to the constructor is stored as a reference to allow easy access to the parameters from all functions of the class.  Since we are working with vector valued finite elements, the FE object we are using is of type FESystem.
* 

* 
* @code
   template <int dim>
   class UltrasoundProblem
   {
   public:
     UltrasoundProblem(ParameterHandler &);
     void run();
* 
   private:
     void make_grid();
     void setup_system();
     void assemble_system();
     void solve();
     void output_results() const;
* 
     ParameterHandler &prm;
* 
     Triangulation<dim> triangulation;
     DoFHandler<dim>    dof_handler;
     FESystem<dim>      fe;
* 
     SparsityPattern      sparsity_pattern;
     SparseMatrix<double> system_matrix;
     Vector<double>       solution, system_rhs;
   };
* 
 
* 
 @endcode
* 
*  The constructor takes the ParameterHandler object and stores it in a reference. It also initializes the DoF-Handler and the finite element system, which consists of two copies of the scalar Q1 field, one for   $v$   and one for   $w$  :
* 

* 
* @code
   template <int dim>
   UltrasoundProblem<dim>::UltrasoundProblem(ParameterHandler &param)
     : prm(param)
     , dof_handler(triangulation)
     , fe(FE_Q<dim>(1), 2)
   {}
* 
 @endcode
* 
*   <a name="codeUltrasoundProblemmake_gridcode"></a>  <h4><code>UltrasoundProblem::make_grid</code></h4>
* 

* 
*  Here we setup the grid for our domain.  As mentioned in the exposition, the geometry is just a unit square (in 2d) with the part of the boundary that represents the transducer lens replaced by a sector of a circle.
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::make_grid()
   {
 @endcode
* 
*  First we generate some logging output and start a timer so we can compute execution time when this function is done:
* 

* 
* @code
     std::cout << "Generating grid... ";
     Timer timer;
* 
 @endcode
* 
*  Then we query the values for the focal distance of the transducer lens and the number of mesh refinement steps from our ParameterHandler object:
* 

* 
* @code
     prm.enter_subsection("Mesh & geometry parameters");
* 
     const double       focal_distance = prm.get_double("Focal distance");
     const unsigned int n_refinements = prm.get_integer("Number of refinements");
* 
     prm.leave_subsection();
* 
 @endcode
* 
*  Next, two points are defined for position and focal point of the transducer lens, which is the center of the circle whose segment will form the transducer part of the boundary. Notice that this is the only point in the program where things are slightly different in 2D and 3D. Even though this tutorial only deals with the 2D case, the necessary additions to make this program functional in 3D are so minimal that we opt for including them:
* 

* 
* @code
     const Point<dim> transducer =
       (dim == 2) ? Point<dim>(0.5, 0.0) : Point<dim>(0.5, 0.5, 0.0);
     const Point<dim> focal_point = (dim == 2) ?
                                      Point<dim>(0.5, focal_distance) :
                                      Point<dim>(0.5, 0.5, focal_distance);
* 
 
 @endcode
* 
*  As initial coarse grid we take a simple unit square with 5 subdivisions in each direction. The number of subdivisions is chosen so that the line segment   $[0.4,0.6]$   that we want to designate as the transducer boundary is spanned by a single face. Then we step through all cells to find the faces where the transducer is to be located, which in fact is just the single edge from 0.4 to 0.6 on the x-axis. This is where we want the refinements to be made according to a circle shaped boundary, so we mark this edge with a different manifold indicator. Since we will Dirichlet boundary conditions on the transducer, we also change its boundary indicator.
* 

* 
* @code
     GridGenerator::subdivided_hyper_cube(triangulation, 5, 0, 1);
* 
     for (auto &cell : triangulation.cell_iterators())
       for (const auto &face : cell->face_iterators())
         if (face->at_boundary() &&
             ((face->center()
* 
- transducer).norm_square() < 0.01))
           {
             face->set_boundary_id(1);
             face->set_manifold_id(1);
           }
 @endcode
* 
*  For the circle part of the transducer lens, a SphericalManifold object is used (which, of course, in 2D just represents a circle), with center computed as above.
* 

* 
* @code
     triangulation.set_manifold(1, SphericalManifold<dim>(focal_point));
* 
 @endcode
* 
*  Now global refinement is executed. Cells near the transducer location will be automatically refined according to the circle shaped boundary of the transducer lens:
* 

* 
* @code
     triangulation.refine_global(n_refinements);
* 
 @endcode
* 
*  Lastly, we generate some more logging output. We stop the timer and query the number of CPU seconds elapsed since the beginning of the function:
* 

* 
* @code
     timer.stop();
     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
* 
     std::cout << "  Number of active cells:  " << triangulation.n_active_cells()
               << std::endl;
   }
* 
 
 @endcode
* 
*   <a name="codeUltrasoundProblemsetup_systemcode"></a>  <h4><code>UltrasoundProblem::setup_system</code></h4>   
*   Initialization of the system matrix, sparsity patterns and vectors are the same as in previous examples and therefore do not need further comment. As in the previous function, we also output the run time of what we do here:
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::setup_system()
   {
     std::cout << "Setting up system... ";
     Timer timer;
* 
     dof_handler.distribute_dofs(fe);
* 
     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
     DoFTools::make_sparsity_pattern(dof_handler, dsp);
     sparsity_pattern.copy_from(dsp);
* 
     system_matrix.reinit(sparsity_pattern);
     system_rhs.reinit(dof_handler.n_dofs());
     solution.reinit(dof_handler.n_dofs());
* 
     timer.stop();
     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
* 
     std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs()
               << std::endl;
   }
* 
 
 @endcode
* 
*   <a name="codeUltrasoundProblemassemble_systemcode"></a>  <h4><code>UltrasoundProblem::assemble_system</code></h4>
* 

* 
*  As before, this function takes care of assembling the system matrix and right hand side vector:
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::assemble_system()
   {
     std::cout << "Assembling system matrix... ";
     Timer timer;
* 
 @endcode
* 
*  First we query wavespeed and frequency from the ParameterHandler object and store them in local variables, as they will be used frequently throughout this function.
* 

* 
*  

* 
* @code
     prm.enter_subsection("Physical constants");
* 
     const double omega = prm.get_double("omega"), c = prm.get_double("c");
* 
     prm.leave_subsection();
* 
 @endcode
* 
*  As usual, for computing integrals ordinary Gauss quadrature rule is used. Since our bilinear form involves boundary integrals on   $\Gamma_2$  , we also need a quadrature rule for surface integration on the faces, which are   $dim-1$   dimensional:
* 

* 
* @code
     QGauss<dim>     quadrature_formula(fe.degree + 1);
     QGauss<dim
* 
- 1> face_quadrature_formula(fe.degree + 1);
* 
     const unsigned int n_q_points      = quadrature_formula.size(),
                        n_face_q_points = face_quadrature_formula.size(),
                        dofs_per_cell   = fe.n_dofs_per_cell();
* 
 @endcode
* 
*  The FEValues objects will evaluate the shape functions for us.  For the part of the bilinear form that involves integration on   $\Omega$  , we'll need the values and gradients of the shape functions, and of course the quadrature weights.  For the terms involving the boundary integrals, only shape function values and the quadrature weights are necessary.
* 

* 
* @code
     FEValues<dim> fe_values(fe,
                             quadrature_formula,
                             update_values | update_gradients |
                               update_JxW_values);
* 
     FEFaceValues<dim> fe_face_values(fe,
                                      face_quadrature_formula,
                                      update_values | update_JxW_values);
* 
 @endcode
* 
*  As usual, the system matrix is assembled cell by cell, and we need a matrix for storing the local cell contributions as well as an index vector to transfer the cell contributions to the appropriate location in the global system matrix after.
* 

* 
* @code
     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
* 
     for (const auto &cell : dof_handler.active_cell_iterators())
       {
 @endcode
* 
*  On each cell, we first need to reset the local contribution matrix and request the FEValues object to compute the shape functions for the current cell:
* 

* 
* @code
         cell_matrix = 0;
         fe_values.reinit(cell);
* 
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
           {
             for (unsigned int j = 0; j < dofs_per_cell; ++j)
               {
 @endcode
* 
*  At this point, it is important to keep in mind that we are dealing with a finite element system with two components. Due to the way we constructed this FESystem, namely as the Cartesian product of two scalar finite element fields, each shape function has only a single nonzero component (they are, in deal.II lingo,   @ref   GlossPrimitive "primitive").  Hence, each shape function can be viewed as one of the   $\phi$  's or   $\psi$  's from the introduction, and similarly the corresponding degrees of freedom can be attributed to either   $\alpha$   or   $\beta$  . As we iterate through all the degrees of freedom on the current cell however, they do not come in any particular order, and so we cannot decide right away whether the DoFs with index   $i$   and   $j$   belong to the real or imaginary part of our solution.  On the other hand, if you look at the form of the system matrix in the introduction, this distinction is crucial since it will determine to which block in the system matrix the contribution of the current pair of DoFs will go and hence which quantity we need to compute from the given two shape functions.  Fortunately, the FESystem object can provide us with this information, namely it has a function   FESystem::system_to_component_index,   that for each local DoF index returns a pair of integers of which the first indicates to which component of the system the DoF belongs. The second integer of the pair indicates which index the DoF has in the scalar base finite element field, but this information is not relevant here. If you want to know more about this function and the underlying scheme behind primitive vector valued elements, take a look at   step-8   or the   @ref vector_valued   module, where these topics are explained in depth.
* 

* 
* @code
                 if (fe.system_to_component_index(i).first ==
                     fe.system_to_component_index(j).first)
                   {
 @endcode
* 
*  If both DoFs   $i$   and   $j$   belong to same component, i.e. their shape functions are both   $\phi$  's or both   $\psi$  's, the contribution will end up in one of the diagonal blocks in our system matrix, and since the corresponding entries are computed by the same formula, we do not bother if they actually are   $\phi$   or   $\psi$   shape functions. We can simply compute the entry by iterating over all quadrature points and adding up their contributions, where values and gradients of the shape functions are supplied by our FEValues object.
* 

* 
*  

* 
* @code
                     for (unsigned int q_point = 0; q_point < n_q_points;
                          ++q_point)
                       cell_matrix(i, j) +=
                         (((fe_values.shape_value(i, q_point)
                            fe_values.shape_value(j, q_point))
                             (-omega omega) +
                           (fe_values.shape_grad(i, q_point)
                            fe_values.shape_grad(j, q_point))
                             c c)
                          fe_values.JxW(q_point));
* 
 @endcode
* 
*  You might think that we would have to specify which component of the shape function we'd like to evaluate when requesting shape function values or gradients from the FEValues object. However, as the shape functions are primitive, they have only one nonzero component, and the FEValues class is smart enough to figure out that we are definitely interested in this one nonzero component.
* 

* 
* @code
                   }
               }
           }
* 
 
 @endcode
* 
*  We also have to add contributions due to boundary terms. To this end, we loop over all faces of the current cell and see if first it is at the boundary, and second has the correct boundary indicator associated with   $\Gamma_2$  , the part of the boundary where we have absorbing boundary conditions:
* 

* 
* @code
         for (const auto face_no : cell->face_indices())
           if (cell->face(face_no)->at_boundary() &&
               (cell->face(face_no)->boundary_id() == 0))
             {
 @endcode
* 
*  These faces will certainly contribute to the off-diagonal blocks of the system matrix, so we ask the FEFaceValues object to provide us with the shape function values on this face:
* 

* 
* @code
               fe_face_values.reinit(cell, face_no);
* 
 
 @endcode
* 
*  Next, we loop through all DoFs of the current cell to find pairs that belong to different components and both have support on the current face_no:
* 

* 
* @code
               for (unsigned int i = 0; i < dofs_per_cell; ++i)
                 for (unsigned int j = 0; j < dofs_per_cell; ++j)
                   if ((fe.system_to_component_index(i).first !=
                        fe.system_to_component_index(j).first) &&
                       fe.has_support_on_face(i, face_no) &&
                       fe.has_support_on_face(j, face_no))
 @endcode
* 
*  The check whether shape functions have support on a face is not strictly necessary: if we don't check for it we would simply add up terms to the local cell matrix that happen to be zero because at least one of the shape functions happens to be zero. However, we can save that work by adding the checks above.
* 

* 
*  In either case, these DoFs will contribute to the boundary integrals in the off-diagonal blocks of the system matrix. To compute the integral, we loop over all the quadrature points on the face and sum up the contribution weighted with the quadrature weights that the face quadrature rule provides.  In contrast to the entries on the diagonal blocks, here it does matter which one of the shape functions is a   $\psi$   and which one is a   $\phi$  , since that will determine the sign of the entry.  We account for this by a simple conditional statement that determines the correct sign. Since we already checked that DoF   $i$   and   $j$   belong to different components, it suffices here to test for one of them to which component it belongs.
* 

* 
* @code
                     for (unsigned int q_point = 0; q_point < n_face_q_points;
                          ++q_point)
                       cell_matrix(i, j) +=
                         ((fe.system_to_component_index(i).first == 0) ?
* 
-1 :
                                                                         1)
                         fe_face_values.shape_value(i, q_point)
                         fe_face_values.shape_value(j, q_point) c omega
                         fe_face_values.JxW(q_point);
             }
* 
 @endcode
* 
*  Now we are done with this cell and have to transfer its contributions from the local to the global system matrix. To this end, we first get a list of the global indices of the this cells DoFs...
* 

* 
* @code
         cell->get_dof_indices(local_dof_indices);
* 
 
 @endcode
* 
*  ...and then add the entries to the system matrix one by one:
* 

* 
* @code
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
           for (unsigned int j = 0; j < dofs_per_cell; ++j)
             system_matrix.add(local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i, j));
       }
* 
 
 @endcode
* 
*  The only thing left are the Dirichlet boundary values on   $\Gamma_1$  , which is characterized by the boundary indicator 1. The Dirichlet values are provided by the   <code>DirichletBoundaryValues</code>   class we defined above:
* 

* 
* @code
     std::map<types::global_dof_index, double> boundary_values;
     VectorTools::interpolate_boundary_values(dof_handler,
                                              1,
                                              DirichletBoundaryValues<dim>(),
                                              boundary_values);
* 
     MatrixTools::apply_boundary_values(boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);
* 
     timer.stop();
     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
   }
* 
 
* 
 @endcode
* 
*   <a name="codeUltrasoundProblemsolvecode"></a>  <h4><code>UltrasoundProblem::solve</code></h4>
* 

* 
*  As already mentioned in the introduction, the system matrix is neither symmetric nor definite, and so it is not quite obvious how to come up with an iterative solver and a preconditioner that do a good job on this matrix.  We chose instead to go a different way and solve the linear system with the sparse LU decomposition provided by UMFPACK. This is often a good first choice for 2D problems and works reasonably well even for a large number of DoFs.  The deal.II interface to UMFPACK is given by the SparseDirectUMFPACK class, which is very easy to use and allows us to solve our linear system with just 3 lines of code.
* 

* 
*  Note again that for compiling this example program, you need to have the deal.II library built with UMFPACK support.
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::solve()
   {
     std::cout << "Solving linear system... ";
     Timer timer;
* 
 @endcode
* 
*  The code to solve the linear system is short: First, we allocate an object of the right type. The following   <code>initialize</code>   call provides the matrix that we would like to invert to the SparseDirectUMFPACK object, and at the same time kicks off the LU-decomposition. Hence, this is also the point where most of the computational work in this program happens.
* 

* 
* @code
     SparseDirectUMFPACK A_direct;
     A_direct.initialize(system_matrix);
* 
 @endcode
* 
*  After the decomposition, we can use   <code>A_direct</code>   like a matrix representing the inverse of our system matrix, so to compute the solution we just have to multiply with the right hand side vector:
* 

* 
* @code
     A_direct.vmult(solution, system_rhs);
* 
     timer.stop();
     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
   }
* 
 
* 
 @endcode
* 
*   <a name="codeUltrasoundProblemoutput_resultscode"></a>  <h4><code>UltrasoundProblem::output_results</code></h4>
* 

* 
*  Here we output our solution   $v$   and   $w$   as well as the derived quantity   $|u|$   in the format specified in the parameter file. Most of the work for deriving   $|u|$   from   $v$   and   $w$   was already done in the implementation of the   <code>ComputeIntensity</code>   class, so that the output routine is rather straightforward and very similar to what is done in the previous tutorials.
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::output_results() const
   {
     std::cout << "Generating output... ";
     Timer timer;
* 
 @endcode
* 
*  Define objects of our   <code>ComputeIntensity</code>   class and a DataOut object:
* 

* 
* @code
     ComputeIntensity<dim> intensities;
     DataOut<dim>          data_out;
* 
     data_out.attach_dof_handler(dof_handler);
* 
 @endcode
* 
*  Next we query the output-related parameters from the ParameterHandler. The   DataOut::parse_parameters   call acts as a counterpart to the   DataOutInterface<1>::declare_parameters   call in   <code>ParameterReader::declare_parameters</code>  . It collects all the output format related parameters from the ParameterHandler and sets the corresponding properties of the DataOut object accordingly.
* 

* 
* @code
     prm.enter_subsection("Output parameters");
* 
     const std::string output_filename = prm.get("Output filename");
     data_out.parse_parameters(prm);
* 
     prm.leave_subsection();
* 
 @endcode
* 
*  Now we put together the filename from the base name provided by the ParameterHandler and the suffix which is provided by the DataOut class (the default suffix is set to the right type that matches the one set in the .prm file through parse_parameters()):
* 

* 
* @code
     const std::string filename = output_filename + data_out.default_suffix();
* 
     std::ofstream output(filename);
* 
 @endcode
* 
*  The solution vectors   $v$   and   $w$   are added to the DataOut object in the usual way:
* 

* 
* @code
     std::vector<std::string> solution_names;
     solution_names.emplace_back("Re_u");
     solution_names.emplace_back("Im_u");
* 
     data_out.add_data_vector(solution, solution_names);
* 
 @endcode
* 
*  For the intensity, we just call   <code>add_data_vector</code>   again, but this with our   <code>ComputeIntensity</code>   object as the second argument, which effectively adds   $|u|$   to the output data:
* 

* 
* @code
     data_out.add_data_vector(solution, intensities);
* 
 @endcode
* 
*  The last steps are as before. Note that the actual output format is now determined by what is stated in the input file, i.e. one can change the output format without having to re-compile this program:
* 

* 
* @code
     data_out.build_patches();
     data_out.write(output);
* 
     timer.stop();
     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl;
   }
* 
 
* 
 @endcode
* 
*   <a name="codeUltrasoundProblemruncode"></a>  <h4><code>UltrasoundProblem::run</code></h4>
* 

* 
*  Here we simply execute our functions one after the other:
* 

* 
* @code
   template <int dim>
   void UltrasoundProblem<dim>::run()
   {
     make_grid();
     setup_system();
     assemble_system();
     solve();
     output_results();
   }
 } // namespace Step29
* 
 
 @endcode
* 
*   <a name="Thecodemaincodefunction"></a>  <h4>The <code>main</code> function</h4>
* 

* 
*  Finally the   <code>main</code>   function of the program. It has the same structure as in almost all of the other tutorial programs. The only exception is that we define ParameterHandler and   <code>ParameterReader</code>   objects, and let the latter read in the parameter values from a textfile called   <code>step-29.prm</code>  . The values so read are then handed over to an instance of the UltrasoundProblem class:
* 

* 
* @code
 int main()
 {
   try
     {
       using namespace dealii;
       using namespace Step29;
* 
       ParameterHandler prm;
       ParameterReader  param(prm);
       param.read_parameters("step-29.prm");
* 
       UltrasoundProblem<2> ultrasound_problem(prm);
       ultrasound_problem.run();
     }
   catch (std::exception &exc)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Exception on processing: " << std::endl
                 << exc.what() << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       return 1;
     }
   catch (...)
     {
       std::cerr << std::endl
                 << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       std::cerr << "Unknown exception!" << std::endl
                 << "Aborting!" << std::endl
                 << "----------------------------------------------------"
                 << std::endl;
       return 1;
     }
   return 0;
 }
 @endcode
* <a name="Results"></a><a name="Results"></a><h1>Results</h1>
* 

* The current program reads its run-time parameters from an input filecalled   <code>\step-29.prm</code>   that looks like this:
* @code
subsection Mesh & geometry parameters
  # Distance of the focal point of the lens to the x-axis
  set Focal distance        = 0.3
* 
  # Number of global mesh refinement steps applied to initial coarse grid
  set Number of refinements = 5
end
* 


* 
subsection Physical constants
  # Wave speed
  set c     = 1.5e5
* 
  # Frequency
  set omega = 3.0e7
end
* 


* 
subsection Output parameters
  # Name of the output file (without extension)
  set Output file   = solution
* 
  # A name for the output format to be used
  set Output format = vtu
end
@endcode
* 
* As can be seen, we set  $d=0.3$  , which amounts to a focus of the transducer lensat   $x=0.5$  ,   $y=0.3$  . The coarse mesh is refined 5 times,resulting in 160x160 cells, and the output is written in vtuformat. The parameter reader understands many more parameterspertaining in particular to the generation of output, but weneed none of these parameters here and therefore stick withtheir default values.
* Here's the console output of the program in debug mode:
* @code
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Debug configuration
Generating grid... done (0.820449s)
  Number of active cells:  25600
Setting up system... done (1.18392s)
  Number of degrees of freedom: 51842
Assembling system matrix... done (2.33291s)
Solving linear system... done (1.34837s)
Generating output... done (2.05782s)
[100%] Built target run
@endcode
* 
* (Of course, execution times will differ if you run the programlocally.) The fact that most of the time is spent on assemblingthe system matrix and generating output is due to the many assertionsthat need to be checked in debug mode. In release mode these partsof the program run much faster whereas solving the linear system ishardly sped up at all:
* @code
> make run
[ 66%] Built target step-29
Scanning dependencies of target run
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0144960s)
DEAL::  Number of active cells:  25600
DEAL::Setting up system... done (0.0356880s)
DEAL::  Number of degrees of freedom: 51842
DEAL::Assembling system matrix... done (0.0436570s)
DEAL::Solving linear system... done (1.54733s)
DEAL::Generating output... done (0.720528s)
[100%] Built target run
@endcode
* 
* The graphical output of the program looks as follows:
* 

*   <table align="center" class="doxtable">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.v.png" alt="v = Re(u)">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.w.png" alt="w = Im(u)">
    </td>
  </tr>
  <tr>
    <td colspan="2">
      <img src="https://www.dealii.org/images/steps/developer/step-29.intensity.png" alt="|u|">
    </td>
  </tr>
</table>  
* The first two pictures show the real and imaginary parts of  $u$  , whereas the last shows the intensity   $|u|$  . One can clearlysee that the intensity is focused around the focal point of thelens (0.5, 0.3), and that the focusis rather sharp in   $x$  -direction but more blurred in   $y$  -direction, which is aconsequence of the geometry of the focusing lens, its finite aperture,and the wave nature of the problem.
* Because colorful graphics are always fun, and to stress the focusingeffects some more, here is another set of images highlighting how wellthe intensity is actually focused in   $x$  -direction:
*   <table align="center" class="doxtable">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.surface.png" alt="|u|">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.contours.png" alt="|u|">
    </td>
  </tr>
</table>  
* 

* As a final note, the structure of the program makes it easy todetermine which parts of the program scale nicely as the mesh isrefined and which parts don't. Here are the run times for 5, 6, and 7global refinements:
* @code
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0135260s)
DEAL::  Number of active cells:  25600
DEAL::Setting up system... done (0.0213910s)
DEAL::  Number of degrees of freedom: 51842
DEAL::Assembling system matrix... done (0.0414300s)
DEAL::Solving linear system... done (1.56621s)
DEAL::Generating output... done (0.729605s)
[100%] Built target run
* 
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0668490s)
DEAL::  Number of active cells:  102400
DEAL::Setting up system... done (0.109694s)
DEAL::  Number of degrees of freedom: 206082
DEAL::Assembling system matrix... done (0.160784s)
DEAL::Solving linear system... done (7.86577s)
DEAL::Generating output... done (2.89320s)
[100%] Built target run
* 
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.293154s)
DEAL::  Number of active cells:  409600
DEAL::Setting up system... done (0.491301s)
DEAL::  Number of degrees of freedom: 821762
DEAL::Assembling system matrix... done (0.605386s)
DEAL::Solving linear system... done (45.1989s)
DEAL::Generating output... done (11.2292s)
[100%] Built target run
@endcode
* 
* Each time we refine the mesh once, so the number of cells and degreesof freedom roughly quadruples from each step to the next. As can be seen,generating the grid, setting up degrees of freedom, assembling thelinear system, and generating output scale pretty closely to linear,whereas solving the linear system is an operation that requires 8times more time each time the number of degrees of freedom isincreased by a factor of 4, i.e. it is   ${\cal O}(N^{3/2})$  . This canbe explained by the fact that (using optimal ordering) thebandwidth of a finite element matrix is   $B={\cal O}(N^{(dim-1)/dim})$  ,and the effort to solve a banded linear system using LU decompositionis   ${\cal O}(BN)$  . This also explains why the program does run in 3das well (after changing the dimension on the  <code>UltrasoundProblem</code>   object), but scales very badly andtakes extraordinary patience before it finishes solving the linearsystem on a mesh with appreciable resolution, even though all theother parts of the program scale very nicely.
* 

* 
* <a name="extensions"></a><a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>
* 

* An obvious possible extension for this program is to run it in 3d&mdash; after all, the world around us is three-dimensional, andultrasound beams propagate in three-dimensional media. You can trythis by simply changing the template parameter of the principal classin   <code>main()</code>   and running it. This won't get you very far,though: certainly not if you do 5 global refinement steps as set inthe parameter file. You'll simply run out of memory as both the mesh(with its   $(2^5)^3 \cdot 5^3=2^{15}\cdot 125 \approx 4\cdot 10^6$   cells)and in particular the sparse direct solver take too much memory. Youcan solve with 3 global refinement steps, however, if you have a bitof time: in early 2011, the direct solve takes about half anhour. What you'll notice, however, is that the solution is completelywrong: the mesh size is simply not small enough to resolve thesolution's waves accurately, and you can see this in plots of thesolution. Consequently, this is one of the cases where adaptivity isindispensable if you don't just want to throw a bigger (presumably%parallel) machine at the problem.
* 

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-29.cc"  
*
*/
