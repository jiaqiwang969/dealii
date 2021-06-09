/**
@page step_68 The step-68 tutorial program
This tutorial depends on step-19.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Simulationofthemotionofmasslesstracerparticlesinavorticalflow">Simulation of the motion of massless tracer particles in a vortical flow</a>
        <li><a href="#ParticlesindealII">Particles in deal.II</a>
        <li><a href="#Challengesrelatedtodistributedparticlesimulations">Challenges related to distributed particle simulations</a>
      <ul>
        <li><a href="#Parallelparticlegeneration">Parallel particle generation</a>
        <li><a href="#Particleexchange">Particle exchange</a>
        <li><a href="#Balancingmeshandparticleload">Balancing mesh and particle load</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Runtimeparameterhandling">Run-time parameter handling</a>
        <li><a href="#Velocityprofile">Velocity profile</a>
        <li><a href="#ThecodeParticleTrackingcodeclassdeclaration">The <code>ParticleTracking</code> class declaration</a>
        <li><a href="#ThecodePatricleTrackingcodeclassimplementation">The <code>PatricleTracking</code> class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#Cellweight">Cell weight</a>
        <li><a href="#Particlesgeneration">Particles generation</a>
        <li><a href="#BackgroundDOFsandinterpolation">Background DOFs and interpolation</a>
        <li><a href="#Timeintegrationofthetrajectories">Time integration of the trajectories</a>
        <li><a href="#Dataoutput">Data output</a>
        <li><a href="#Runningthesimulation">Running the simulation</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Motionoftheparticles"> Motion of the particles </a>
        <li><a href="#Dynamicloadbalancing"> Dynamic load balancing </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-68/doc/intro.dox

 <br> 

<i>
This program was contributed by
Bruno Blais (Polytechnique Montréal),
Toni El Geitani Nehme (Polytechnique Montréal),
Rene Gassmöller (University of California Davis),
and Peter Munch (Technical University of Munich and Helmholtz-Zentrum Geesthacht).
Bruno Blais was supported by NSERC Discovery grant
RGPIN-2020-04510, by Compute Canada and Calcul Québec.
</i>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Simulationofthemotionofmasslesstracerparticlesinavorticalflow"></a><h3>Simulation of the motion of massless tracer particles in a vortical flow</h3>


粒子在大量应用的数值模型中发挥着重要作用。粒子通常被用作无质量追踪器，以显示瞬时流动的动态。它们也可以作为更复杂的有限元模型的一部分发挥固有的作用，如颗粒在细胞中（PIC）方法 @cite GLHPW2018 ，或者它们甚至可以用来模拟颗粒物质的运动，如离散元素法（DEM） @cite Blais2019  。在DEM的情况下，所产生的模型不再与有限元方法有关，而只是导致了一个描述颗粒运动和它们碰撞动态的常微分方程系统。所有这些模型都可以用deal.II的粒子处理能力来建立。

在本步骤中，我们使用粒子作为无质量的追踪器来说明涡流的动态。由于粒子是无质量追踪器，每个粒子 $i$ 的位置由以下常微分方程（ODE）描述。

@f[
\frac{d \textbf{x}_i}{dt} =\textbf{u}(\textbf{x}_i)


@f]



其中 $\textbf{x}_i$ 是粒子 $i$ 的位置， $\textbf{u}(\textbf{x}_i)$ 是其位置上的流速。在本步骤中，该ODE使用显式欧拉方法进行求解。由此产生的方案是。

@f[
\textbf{x}_{i}^{n+1} = \textbf{x}_{i}^{n} + \Delta t \; \textbf{u}(\textbf{x}_{i}^{n})


@f]



其中 $\textbf{x}_{i}^{n+1}$ 和 $\textbf{x}_{i}^{n}$ 分别是粒子 $i$ 在时间 $t+\Delta t$ 和 $t$ 的位置，其中 $\Delta t$ 是时间步骤。在本步骤中，粒子位置的速度以两种不同的方式获得。

- 通过在粒子的位置评估速度函数。

- 通过在背景三角图上评估速度函数，并使用有限元支持，在粒子的位置上进行插值。

第一种方法是不实际的，因为速度曲线一般是不知道的分析。第二种方法，基于在粒子位置的内插解，完全模仿了在现实的计算流体动力学模拟中的做法，这也是我们在步骤19中对粒子位置的有限元解进行评估的方式。在这一步中，我们说明了这两种策略。

我们注意到，通过使用四阶Runge-Kutta方法或其他适当的方案对粒子的运动进行时间积分，可以获得更大的精度。  实施一个更先进的时间积分方案将是这一步骤的直接延伸。

<a name="ParticlesindealII"></a><h3>Particles in deal.II</h3>


在deal.II中， Particles::Particle 是非常简单和灵活的实体，可以用来建立PIC、DEM或任何类型的基于粒子的模型。粒子在现实空间中有一个位置，在它们所在的元素的参考空间中有一个位置，还有一个唯一的ID。在大多数情况下，包含粒子的模拟需要大量的粒子。因此，通过一个聚集所有粒子的实体来处理所有的粒子变得很有趣。在deal.II中，这是通过使用 Particles::ParticleHandler 类来实现的。

默认情况下，粒子没有直径、质量或任何其他我们通常期望的物理粒子的物理属性。然而，通过ParticleHandler，粒子可以访问一个 Particles::PropertyPool. PropertyPool是一个数组，可以用来存储与粒子相关的任意数量的属性。因此，用户可以建立自己的粒子解算器，并将所需的属性归属于粒子（例如，质量、电荷、直径、温度等）。在本教程中，这被用来存储流体速度的值和粒子所属的过程ID。

<a name="Challengesrelatedtodistributedparticlesimulations"></a><h3>Challenges related to distributed particle simulations</h3>


尽管本步骤不是计算密集型的，但包括许多粒子的模拟可能对计算要求很高，需要并行化。本步骤展示了deal.II对粒子的分布式并行能力。一般来说，在包括粒子的平行分布式模拟中，有三个主要挑战。

- 在分布式三角图上生成粒子。

- 在处理器之间交换离开本地域的粒子。

- 对模拟进行负载平衡，使每个处理器都有类似的计算负载。这些挑战及其在交易.II中的解决方案已经在 @cite GLHPW2018 中进行了更详细的讨论，但我们将在下面进行总结。

当然也有关于简单设置使用粒子的代码的问题。这些问题在第19步中已经基本解决了。一些更高级的技术也将在第70步中讨论。

<a name="Parallelparticlegeneration"></a><h4>Parallel particle generation</h4>


以可扩展的方式生成分布式粒子并不简单，因为在找到它们所在的单元之前，必须首先确定它们所属的处理器。 deal.II通过 Particles::Generator 命名空间提供了许多生成粒子的能力。  其中一些粒子生成器只在本地拥有的子域上创建粒子。例如， Particles::Generators::regular_reference_locations() 在本地子域的每个单元内的相同参考位置创建粒子， Particles::Generators::probabilistic_locations() 使用全局定义的概率密度函数来确定本地生成粒子的数量和位置。

在其他情况下，如本步骤，粒子必须在单元格上的特定位置生成，而这些单元格可能只由处理器的一个子集拥有。在大多数这些情况下，粒子的插入是在非常有限的时间步长内完成的，因此，不构成计算成本的很大一部分。对于这些情况，deal.II提供了方便的 Particles::Generators ，可以在全局范围内插入粒子，即使粒子不在启动创建粒子的调用的并行进程所拥有的单元中。生成器首先定位粒子位于哪个子域上，确定它们位于哪个单元中，并在处理器之间交换必要的信息，以确保生成的粒子具有正确的属性。因此，这种类型的粒子生成可能是通信密集型的。 Particles::Generators::dof_support_points 和 Particles::Generators::quadrature_points 分别使用三角法和相关DoFHandler或正交的点来生成粒子。用于生成粒子的三角形可以是用于背景网格的同一三角形，在这种情况下，这些函数与上一段中描述的 Particles::Generators::regular_reference_locations() 函数非常相似。然而，用于生成粒子的三角法也可以与背景网格的三角法不同（不匹配），这对于生成特定形状的粒子（如本例），或者在两个不同的计算网格之间传输信息（如步骤-70）是很有用的。  此外， Particles::ParticleHandler 类提供了 Particles::ParticleHandler::insert_global_particles() 函数，可以从任意点的矢量和边界框的全局矢量中插入粒子。在本步骤中，我们在非匹配三角形上使用 Particles::Generators::quadrature_points() 函数来插入位于圆盘形状位置的粒子。

<a name="Particleexchange"></a><h4>Particle exchange</h4>


当粒子在平行分布式计算中移动时，它们可能会离开本地拥有的子域，需要转移到它们的新主人进程中。这种情况可能以两种非常不同的方式出现。首先，如果先前拥有的进程知道丢失的粒子的新主人（例如，因为粒子从一个处理器的本地拥有的单元移动到分布式三角形的相邻的幽灵单元），那么转移可以作为每个进程和新主人之间的点对点通信有效处理。每当粒子被分类到新单元时，这种转移就会自动发生。其次，以前的所有者可能不知道粒子被转移到哪个进程。在这种情况下，粒子被默认丢弃，因为全局搜索所有者的成本很高。步骤19显示了这样一个被丢弃的粒子仍然可以被收集、解释，并可能被用户重新插入。在本例中，我们通过在时间步长上施加一个CFL准则来防止第二种情况，以确保粒子最多会移动到本地进程的幽灵层，因此可以自动发送到邻近的进程。

<a name="Balancingmeshandparticleload"></a><h4>Balancing mesh and particle load</h4>


在使用粒子的并行分布式计算中出现的最后一个挑战是平衡计算负荷，即在网格上完成的工作，例如解决有限元问题，和在粒子上完成的工作，例如移动粒子或计算粒子之间或粒子与网格之间的力。默认情况下，例如在步骤40中，deal.II尽可能均匀地将背景网格分配给可用的进程，也就是说，它平衡每个进程上的单元数量。然而，如果一些单元拥有比其他单元多得多的粒子，或者一个单元的粒子比其他单元的粒子计算成本高得多，那么这个问题就不再有效地扩展了（关于我们认为的 "可扩展 "程序的讨论，见 @ref GlossParallelScaling "这个词汇条"）。因此，我们必须应用一种 "负载平衡 "的形式，这意味着我们估计与每个单元及其粒子相关的计算负载。然后，重新划分网格就会考虑到这个综合的计算负荷，而不是单元数的简化假设  @cite GLHPW2018  。

在本节中，我们只讨论了分布式计算中针对粒子的挑战。粒子与有限元解决方案共享的并行挑战（并行输出，网格细化过程中的数据传输）可以用其他例子中已经讨论过的有限元问题的解决方案来解决。

<a name="Thetestcase"></a><h3>The testcase</h3>


在本步骤中，我们使用粒子作为无质量的追踪器来说明一个特殊的涡流的动力学：Rayleigh--Kothe涡流。这种流动模式通常被用作界面跟踪方法（如流体体积法和水平集法）的复杂测试案例，因为它导致了流体的强烈旋转和伸长  @cite Blais2013  。

这个Rayleigh-Kothe涡流的流函数 $\Psi$ 被定义为。

@f[
\Psi = \frac{1}{\pi} \sin^2 (\pi x) \sin^2 (\pi y) \cos \left( \pi \frac{t}{T} \right)


@f]

其中 $T$ 为流动的一半周期。二维的速度曲线 (  $\textbf{u}=[u,v]^T$  ) 是 :

@f{eqnarray*}
   u &=&  - \frac{\partial\Psi}{\partial y} = -2 \sin^2 (\pi x) \sin (\pi y) \cos (\pi y)  \cos \left( \pi \frac{t}{T} \right)\\
   v &=&  \frac{\partial\Psi}{\partial x} = 2 \cos(\pi x) \sin(\pi x) \sin^2 (\pi y) \cos \left( \pi \frac{t}{T} \right)


@f}



速度曲线在下面的动画中得到说明。

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/m6hQm7etji8"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly



可以看出，由于项 $\cos \left( \pi \frac{t}{T} \right)$ 的存在，这个速度会周期性地逆转，材料在每一个长度为 $t=2T$ 的周期后都会结束在其起始位置。我们将正好运行这个教程程序一个周期，并将最终的粒子位置与初始位置进行比较，以说明这一流动特性。这个例子使用测试案例产生了两个模型，它们对粒子的处理略有不同。第一个模型将精确的分析速度解作为每个粒子的速度。因此，在这个模型中，分配给粒子的速度没有错误，而在给定时间内，粒子位置与分析位置的任何偏差都是由于使用时间步进法不准确地求解粒子运动方程的错误造成的。在第二个模型中，分析速度场首先被内插到一个有限元矢量空间（以模拟速度是通过求解有限元问题得到的情况，与步骤19中每个粒子的ODE取决于有限元解的方式相同）。然后在粒子的位置上评估这个有限元 "解决方案"，以解决其运动方程。这两种情况之间的差异允许评估所选择的有限元空间是否足够精确，以所选择的粒子平流方案的最佳收敛率来平流粒子，这个问题在实践中对确定组合算法的精度很重要（例如，见 @cite Gassmoller2019  ）。


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
 * 
 * @code
 * #include <deal.II/base/bounding_box.h>
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/base/mpi.h>
 * #include <deal.II/base/parameter_acceptor.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/distributed/cell_weights.h>
 * #include <deal.II/distributed/solution_transfer.h>
 * #include <deal.II/distributed/tria.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/mapping_q.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * 
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * From the following include file we import the ParticleHandler class
 * that allows you to manage
 * a collection of particles (objects of type Particles::Particle), representing
 * a collection of points with some attached properties (e.g., an id) floating
 * on a parallel::distributed::Triangulation. The methods and classes in the
 * namespace Particles allows one to easily implement Particle-In-Cell methods
 * and particle tracing on distributed triangulations:
 * 
 * @code
 * #include <deal.II/particles/particle_handler.h>
 * 
 * @endcode
 * 
 * We import the particles generator
 * which allow us to insert the particles. In the present step, the particle
 * are globally inserted using a non-matching hyper-shell triangulation:
 * 
 * @code
 * #include <deal.II/particles/generators.h>
 * 
 * @endcode
 * 
 * Since the particles do not form a triangulation, they have their
 * own specific DataOut class which will enable us to write them
 * to commonly used parallel vtu format (or any number of other file formats):
 * 
 * @code
 * #include <deal.II/particles/data_out.h>
 * 
 * #include <cmath>
 * #include <iostream>
 * 
 * 
 * 
 * namespace Step68
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameterhandling"></a> 
 * <h3>Run-time parameter handling</h3>
 * 

 * 
 * Similarly to what is done in step-60, we set up a class that holds
 * all the parameters of our problem and derive it from the ParameterAcceptor
 * class to simplify the management and creation of parameter files.
 *   

 * 
 * The ParameterAcceptor paradigm requires all parameters to be writable by
 * the ParameterAcceptor methods. In order to avoid bugs that would be very
 * difficult to track down (such as writing things like `if (time = 0)`
 * instead of `if(time == 0)`), we declare all the parameters in an external
 * class, which is initialized before the actual `ParticleTracking` class, and
 * pass it to the main class as a `const` reference.
 *   

 * 
 * The constructor of the class is responsible for the connection between the
 * members of this class and the corresponding entries in the
 * ParameterHandler. Thanks to the use of the
 * ParameterHandler::add_parameter() method, this connection is trivial, but
 * requires all members of this class to be writable.
 * 
 * @code
 *   class ParticleTrackingParameters : public ParameterAcceptor
 *   {
 *   public:
 *     ParticleTrackingParameters();
 * 
 * @endcode
 * 
 * This class consists largely of member variables that
 * describe the details of the particle tracking simulation and its
 * discretization. The following parameters are about where output should
 * written to, the spatial discretization of the velocity (the default is
 * $Q_1$), the time step and the output frequency (how many time steps
 * should elapse before we generate graphical output again):
 * 
 * @code
 *     std::string output_directory = "./";
 * 
 *     unsigned int velocity_degree       = 1;
 *     double       time_step             = 0.002;
 *     double       final_time            = 4.0;
 *     unsigned int output_frequency      = 10;
 *     unsigned int repartition_frequency = 5;
 * 
 * @endcode
 * 
 * We allow every grid to be refined independently. In this tutorial, no
 * physics is resolved on the fluid grid, and its velocity is calculated
 * analytically.
 * 
 * @code
 *     unsigned int fluid_refinement              = 4;
 *     unsigned int particle_insertion_refinement = 3;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * There remains the task of declaring what run-time parameters we can accept
 * in input files. Since we have a very limited number of parameters, all
 * parameters are declared in the same section.
 * 
 * @code
 *   ParticleTrackingParameters::ParticleTrackingParameters()
 *     : ParameterAcceptor("Particle Tracking Problem/")
 *   {
 *     add_parameter(
 *       "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1));
 * 
 *     add_parameter("Output frequency",
 *                   output_frequency,
 *                   "Iteration frequency at which output results are written",
 *                   prm,
 *                   Patterns::Integer(1));
 * 
 *     add_parameter("Repartition frequency",
 *                   repartition_frequency,
 *                   "Iteration frequency at which the mesh is load balanced",
 *                   prm,
 *                   Patterns::Integer(1));
 * 
 *     add_parameter("Output directory", output_directory);
 * 
 *     add_parameter("Time step", time_step, "", prm, Patterns::Double());
 * 
 *     add_parameter("Final time",
 *                   final_time,
 *                   "End time of the simulation",
 *                   prm,
 *                   Patterns::Double());
 * 
 *     add_parameter("Fluid refinement",
 *                   fluid_refinement,
 *                   "Refinement level of the fluid domain",
 *                   prm,
 *                   Patterns::Integer(0));
 * 
 *     add_parameter(
 *       "Particle insertion refinement",
 *       particle_insertion_refinement,
 *       "Refinement of the volumetric mesh used to insert the particles",
 *       prm,
 *       Patterns::Integer(0));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Velocityprofile"></a> 
 * <h3>Velocity profile</h3>
 * 

 * 
 * The velocity profile is provided as a Function object.
 * This function is hard-coded within
 * the example.
 * 
 * @code
 *   template <int dim>
 *   class Vortex : public Function<dim>
 *   {
 *   public:
 *     Vortex()
 *       : Function<dim>(dim)
 *     {}
 * 
 * 
 *     virtual void vector_value(const Point<dim> &point,
 *                               Vector<double> &  values) const override;
 *   };
 * 
 * 
 * @endcode
 * 
 * The velocity profile for the Rayleigh-Kothe vertex is time-dependent.
 * Consequently, the current time in the
 * simulation (t) must be gathered from the Function object.
 * 
 * @code
 *   template <int dim>
 *   void Vortex<dim>::vector_value(const Point<dim> &point,
 *                                  Vector<double> &  values) const
 *   {
 *     const double T = 4;
 *     const double t = this->get_time();
 * 
 *     const double px = numbers::PI * point(0);
 *     const double py = numbers::PI * point(1);
 *     const double pt = numbers::PI / T * t;
 * 
 *     values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py);
 *     values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px);
 *     if (dim == 3)
 *       {
 *         values[2] = 0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeParticleTrackingcodeclassdeclaration"></a> 
 * <h3>The <code>ParticleTracking</code> class declaration</h3>
 * 

 * 
 * We are now ready to introduce the main class of our tutorial program.
 * 
 * @code
 *   template <int dim>
 *   class ParticleTracking
 *   {
 *   public:
 *     ParticleTracking(const ParticleTrackingParameters &par,
 *                      const bool                        interpolated_velocity);
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * This function is responsible for the initial
 * generation of the particles on top of the background grid.
 * 
 * @code
 *     void generate_particles();
 * 
 * @endcode
 * 
 * When the velocity profile is interpolated to the position of the
 * particles, it must first be stored using degrees of freedom.
 * Consequently, as is the case for other parallel case (e.g. step-40) we
 * initialize the degrees of freedom on the background grid.
 * 
 * @code
 *     void setup_background_dofs();
 * 
 * @endcode
 * 
 * In one of the test cases, the function is mapped to the background grid
 * and a finite element interpolation is used to calculate the velocity
 * at the particle location. This function calculates the value of the
 * function at the support point of the triangulation.
 * 
 * @code
 *     void interpolate_function_to_field();
 * 
 * @endcode
 * 
 * The next two functions are responsible for carrying out step of explicit
 * Euler time integration for the cases where the velocity field is
 * interpolated at the positions of the particles or calculated
 * analytically, respectively.
 * 
 * @code
 *     void euler_step_interpolated(const double dt);
 *     void euler_step_analytical(const double dt);
 * 
 * @endcode
 * 
 * The `cell_weight()` function indicates to the triangulation how much
 * computational work is expected to happen on this cell, and consequently
 * how the domain needs to be partitioned so that every MPI rank receives a
 * roughly equal amount of work (potentially not an equal number of cells).
 * While the function is called from the outside, it is connected to the
 * corresponding signal from inside this class, therefore it can be
 * `private`.
 * 
 * @code
 *     unsigned int cell_weight(
 *       const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *         &cell,
 *       const typename parallel::distributed::Triangulation<dim>::CellStatus
 *         status) const;
 * 
 * @endcode
 * 
 * The following two functions are responsible for outputting the simulation
 * results for the particles and for the velocity profile on the background
 * mesh, respectively.
 * 
 * @code
 *     void output_particles(const unsigned int it);
 *     void output_background(const unsigned int it);
 * 
 * @endcode
 * 
 * The private members of this class are similar to other parallel deal.II
 * examples. The parameters are stored as a `const` member. It is important
 * to note that we keep the `Vortex` class as a member since its time
 * must be modified as the simulation proceeds.
 * 

 * 
 * 
 * @code
 *     const ParticleTrackingParameters &par;
 * 
 *     MPI_Comm                                  mpi_communicator;
 *     parallel::distributed::Triangulation<dim> background_triangulation;
 *     Particles::ParticleHandler<dim>           particle_handler;
 * 
 *     DoFHandler<dim>                            fluid_dh;
 *     FESystem<dim>                              fluid_fe;
 *     MappingQ1<dim>                             mapping;
 *     LinearAlgebra::distributed::Vector<double> velocity_field;
 * 
 *     Vortex<dim> velocity;
 * 
 *     ConditionalOStream pcout;
 * 
 *     bool interpolated_velocity;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePatricleTrackingcodeclassimplementation"></a> 
 * <h3>The <code>PatricleTracking</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * The constructors and destructors are rather trivial. They are very similar
 * to what is done in step-40. We set the processors we want to work on
 * to all machines available (`MPI_COMM_WORLD`) and
 * initialize the <code>pcout</code> variable to only allow processor zero
 * to output anything to the standard output.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
 *                                           const bool interpolated_velocity)
 *     : par(par)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , background_triangulation(mpi_communicator)
 *     , fluid_dh(background_triangulation)
 *     , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *     , interpolated_velocity(interpolated_velocity)
 * 
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Cellweight"></a> 
 * <h4>Cell weight</h4>
 * 

 * 
 * This function is the key component that allow us to dynamically balance the
 * computational load for this example. The function attributes a weight to
 * every cell that represents the computational work on this cell. Here the
 * majority of work is expected to happen on the particles, therefore the
 * return value of this function (representing "work for this cell") is
 * calculated based on the number of particles in the current cell.
 * The function is
 * connected to the cell_weight() signal inside the triangulation, and will be
 * called once per cell, whenever the triangulation repartitions the domain
 * between ranks (the connection is created inside the
 * generate_particles() function of this class).
 * 
 * @code
 *   template <int dim>
 *   unsigned int ParticleTracking<dim>::cell_weight(
 *     const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *       &                                                                  cell,
 *     const typename parallel::distributed::Triangulation<dim>::CellStatus status)
 *     const
 *   {
 * @endcode
 * 
 * We do not assign any weight to cells we do not own (i.e., artificial
 * or ghost cells)
 * 
 * @code
 *     if (!cell->is_locally_owned())
 *       return 0;
 * 
 * @endcode
 * 
 * This determines how important particle work is compared to cell
 * work (by default every cell has a weight of 1000).
 * We set the weight per particle much higher to indicate that
 * the particle load is the only one that is important to distribute the
 * cells in this example. The optimal value of this number depends on the
 * application and can range from 0 (cheap particle operations,
 * expensive cell operations) to much larger than 1000 (expensive
 * particle operations, cheap cell operations, like presumed in this
 * example).
 * 
 * @code
 *     const unsigned int particle_weight = 10000;
 * 
 * @endcode
 * 
 * This example does not use adaptive refinement, therefore every cell
 * should have the status `CELL_PERSIST`. However this function can also
 * be used to distribute load during refinement, therefore we consider
 * refined or coarsened cells as well.
 * 
 * @code
 *     if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
 *         status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
 *       {
 *         const unsigned int n_particles_in_cell =
 *           particle_handler.n_particles_in_cell(cell);
 *         return n_particles_in_cell * particle_weight;
 *       }
 *     else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
 *       {
 *         unsigned int n_particles_in_cell = 0;
 * 
 *         for (unsigned int child_index = 0; child_index < cell->n_children();
 *              ++child_index)
 *           n_particles_in_cell +=
 *             particle_handler.n_particles_in_cell(cell->child(child_index));
 * 
 *         return n_particles_in_cell * particle_weight;
 *       }
 * 
 *     Assert(false, ExcInternalError());
 *     return 0;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Particlesgeneration"></a> 
 * <h4>Particles generation</h4>
 * 

 * 
 * This function generates the tracer particles and the background
 * triangulation on which these particles evolve.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::generate_particles()
 *   {
 * @endcode
 * 
 * We create a hyper cube triangulation which we globally refine. This
 * triangulation covers the full trajectory of the particles.
 * 
 * @code
 *     GridGenerator::hyper_cube(background_triangulation, 0, 1);
 *     background_triangulation.refine_global(par.fluid_refinement);
 * 
 * @endcode
 * 
 * In order to consider the particles when repartitioning the triangulation
 * the algorithm needs to know three things:
 *     

 * 
 * 1. How much weight to assign to each cell (how many particles are in
 * there);
 * 2. How to pack the particles before shipping data around;
 * 3. How to unpack the particles after repartitioning.
 *     

 * 
 * We attach the correct functions to the signals inside
 * parallel::distributed::Triangulation. These signal will be called every
 * time the repartition() function is called. These connections only need to
 * be created once, so we might as well have set them up in the constructor
 * of this class, but for the purpose of this example we want to group the
 * particle related instructions.
 * 
 * @code
 *     background_triangulation.signals.cell_weight.connect(
 *       [&](
 *         const typename parallel::distributed::Triangulation<dim>::cell_iterator
 *           &cell,
 *         const typename parallel::distributed::Triangulation<dim>::CellStatus
 *           status) -> unsigned int { return this->cell_weight(cell, status); });
 * 
 *     background_triangulation.signals.pre_distributed_repartition.connect(
 *       [this]() { this->particle_handler.register_store_callback_function(); });
 * 
 *     background_triangulation.signals.post_distributed_repartition.connect(
 *       [&]() { this->particle_handler.register_load_callback_function(false); });
 * 
 * @endcode
 * 
 * This initializes the background triangulation where the particles are
 * living and the number of properties of the particles.
 * 
 * @code
 *     particle_handler.initialize(background_triangulation, mapping, 1 + dim);
 * 
 * @endcode
 * 
 * We create a particle triangulation which is solely used to generate
 * the points which will be used to insert the particles. This
 * triangulation is a hyper shell which is offset from the
 * center of the simulation domain. This will be used to generate a
 * disk filled with particles which will allow an easy monitoring
 * of the motion due to the vortex.
 * 
 * @code
 *     Point<dim> center;
 *     center[0] = 0.5;
 *     center[1] = 0.75;
 *     if (dim == 3)
 *       center[2] = 0.5;
 * 
 *     const double outer_radius = 0.15;
 *     const double inner_radius = 0.01;
 * 
 *     parallel::distributed::Triangulation<dim> particle_triangulation(
 *       MPI_COMM_WORLD);
 * 
 *     GridGenerator::hyper_shell(
 *       particle_triangulation, center, inner_radius, outer_radius, 6);
 *     particle_triangulation.refine_global(par.particle_insertion_refinement);
 * 
 * @endcode
 * 
 * We generate the necessary bounding boxes for the particles generator.
 * These bounding boxes are required to quickly identify in which
 * process's subdomain the inserted particle lies, and which cell owns it.
 * 
 * @code
 *     const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
 *       background_triangulation, IteratorFilters::LocallyOwnedCell());
 *     const auto global_bounding_boxes =
 *       Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);
 * 
 * @endcode
 * 
 * We generate an empty vector of properties. We will attribute the
 * properties to the particles once they are generated.
 * 
 * @code
 *     std::vector<std::vector<double>> properties(
 *       particle_triangulation.n_locally_owned_active_cells(),
 *       std::vector<double>(dim + 1, 0.));
 * 
 * @endcode
 * 
 * We generate the particles at the position of a single
 * point quadrature. Consequently, one particle will be generated
 * at the centroid of each cell.
 * 
 * @code
 *     Particles::Generators::quadrature_points(particle_triangulation,
 *                                              QMidpoint<dim>(),
 *                                              global_bounding_boxes,
 *                                              particle_handler,
 *                                              mapping,
 *                                              properties);
 * 
 *     pcout << "Number of particles inserted: "
 *           << particle_handler.n_global_particles() << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="BackgroundDOFsandinterpolation"></a> 
 * <h4>Background DOFs and interpolation</h4>
 * 

 * 
 * This function sets up the background degrees of freedom used for the
 * velocity interpolation and allocates the field vector where the entire
 * solution of the velocity field is stored.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::setup_background_dofs()
 *   {
 *     fluid_dh.distribute_dofs(fluid_fe);
 *     const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
 *     IndexSet       locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);
 * 
 *     velocity_field.reinit(locally_owned_dofs,
 *                           locally_relevant_dofs,
 *                           mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This function takes care of interpolating the
 * vortex velocity field to the field vector. This is achieved rather easily
 * by using the VectorTools::interpolate() function.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::interpolate_function_to_field()
 *   {
 *     velocity_field.zero_out_ghost_values();
 *     VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
 *     velocity_field.update_ghost_values();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Timeintegrationofthetrajectories"></a> 
 * <h4>Time integration of the trajectories</h4>
 * 

 * 
 * We integrate the particle trajectories
 * using an analytically defined velocity field. This demonstrates a
 * relatively trivial usage of the particles.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::euler_step_analytical(const double dt)
 *   {
 *     const unsigned int this_mpi_rank =
 *       Utilities::MPI::this_mpi_process(mpi_communicator);
 *     Vector<double> particle_velocity(dim);
 * 
 * @endcode
 * 
 * Looping over all particles in the domain using a particle iterator
 * 
 * @code
 *     for (auto &particle : particle_handler)
 *       {
 * @endcode
 * 
 * We calculate the velocity of the particles using their current
 * location.
 * 
 * @code
 *         Point<dim> particle_location = particle.get_location();
 *         velocity.vector_value(particle_location, particle_velocity);
 * 
 * @endcode
 * 
 * This updates the position of the particles and sets the old position
 * equal to the new position of the particle.
 * 
 * @code
 *         for (int d = 0; d < dim; ++d)
 *           particle_location[d] += particle_velocity[d] * dt;
 * 
 *         particle.set_location(particle_location);
 * 
 * @endcode
 * 
 * We store the processor id (a scalar) and the particle velocity (a
 * vector) in the particle properties. In this example, this is done
 * purely for visualization purposes.
 * 
 * @code
 *         ArrayView<double> properties = particle.get_properties();
 *         for (int d = 0; d < dim; ++d)
 *           properties[d] = particle_velocity[d];
 *         properties[dim] = this_mpi_rank;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In contrast to the previous function in this function we
 * integrate the particle trajectories by interpolating the value of
 * the velocity field at the degrees of freedom to the position of
 * the particles.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::euler_step_interpolated(const double dt)
 *   {
 *     Vector<double> local_dof_values(fluid_fe.dofs_per_cell);
 * 
 * @endcode
 * 
 * We loop over all the local particles. Although this could be achieved
 * directly by looping over all the cells, this would force us
 * to loop over numerous cells which do not contain particles.
 * Rather, we loop over all the particles, but, we get the reference
 * of the cell in which the particle lies and then loop over all particles
 * within that cell. This enables us to gather the values of the velocity
 * out of the `velocity_field` vector once and use them for all particles
 * that lie within the cell.
 * 
 * @code
 *     auto particle = particle_handler.begin();
 *     while (particle != particle_handler.end())
 *       {
 *         const auto cell =
 *           particle->get_surrounding_cell(background_triangulation);
 *         const auto dh_cell =
 *           typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);
 * 
 *         dh_cell->get_dof_values(velocity_field, local_dof_values);
 * 
 * @endcode
 * 
 * Next, compute the velocity at the particle locations by evaluating
 * the finite element solution at the position of the particles.
 * This is essentially an optimized version of the particle advection
 * functionality in step 19, but instead of creating quadrature
 * objects and FEValues objects for each cell, we do the
 * evaluation by hand, which is somewhat more efficient and only
 * matters for this tutorial, because the particle work is the
 * dominant cost of the whole program.
 * 
 * @code
 *         const auto pic = particle_handler.particles_in_cell(cell);
 *         Assert(pic.begin() == particle, ExcInternalError());
 *         for (auto &p : pic)
 *           {
 *             const Point<dim> reference_location = p.get_reference_location();
 *             Tensor<1, dim>   particle_velocity;
 *             for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
 *               {
 *                 const auto comp_j = fluid_fe.system_to_component_index(j);
 * 
 *                 particle_velocity[comp_j.first] +=
 *                   fluid_fe.shape_value(j, reference_location) *
 *                   local_dof_values[j];
 *               }
 * 
 *             Point<dim> particle_location = particle->get_location();
 *             for (int d = 0; d < dim; ++d)
 *               particle_location[d] += particle_velocity[d] * dt;
 *             p.set_location(particle_location);
 * 
 * @endcode
 * 
 * Again, we store the particle velocity and the processor id in the
 * particle properties for visualization purposes.
 * 
 * @code
 *             ArrayView<double> properties = p.get_properties();
 *             for (int d = 0; d < dim; ++d)
 *               properties[d] = particle_velocity[d];
 * 
 *             properties[dim] =
 *               Utilities::MPI::this_mpi_process(mpi_communicator);
 * 
 *             ++particle;
 *           }
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Dataoutput"></a> 
 * <h4>Data output</h4>
 * 

 * 
 * The next two functions take care of writing both the particles
 * and the background mesh to vtu with a pvtu record. This ensures
 * that the simulation results can be visualized when the simulation is
 * launched in parallel.
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::output_particles(const unsigned int it)
 *   {
 *     Particles::DataOut<dim, dim> particle_output;
 * 
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     solution_names.push_back("process_id");
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     particle_output.build_patches(particle_handler,
 *                                   solution_names,
 *                                   data_component_interpretation);
 *     const std::string output_folder(par.output_directory);
 *     const std::string file_name(interpolated_velocity ?
 *                                   "interpolated-particles" :
 *                                   "analytical-particles");
 * 
 *     pcout << "Writing particle output file: " << file_name << "-" << it
 *           << std::endl;
 * 
 *     particle_output.write_vtu_with_pvtu_record(
 *       output_folder, file_name, it, mpi_communicator, 6);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void ParticleTracking<dim>::output_background(const unsigned int it)
 *   {
 *     std::vector<std::string> solution_names(dim, "velocity");
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 * 
 *     DataOut<dim> data_out;
 * 
 * @endcode
 * 
 * Attach the solution data to data_out object
 * 
 * @code
 *     data_out.attach_dof_handler(fluid_dh);
 *     data_out.add_data_vector(velocity_field,
 *                              solution_names,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 *     Vector<float> subdomain(background_triangulation.n_active_cells());
 *     for (unsigned int i = 0; i < subdomain.size(); ++i)
 *       subdomain(i) = background_triangulation.locally_owned_subdomain();
 *     data_out.add_data_vector(subdomain, "subdomain");
 * 
 *     data_out.build_patches(mapping);
 * 
 *     const std::string output_folder(par.output_directory);
 *     const std::string file_name("background");
 * 
 *     pcout << "Writing background field file: " << file_name << "-" << it
 *           << std::endl;
 * 
 *     data_out.write_vtu_with_pvtu_record(
 *       output_folder, file_name, it, mpi_communicator, 6);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Runningthesimulation"></a> 
 * <h4>Running the simulation</h4>
 * This function orchestrates the entire simulation. It is very similar
 * to the other time dependent tutorial programs -- take step-21 or step-26 as
 * an example. Note that we use the DiscreteTime class to monitor the time,
 * the time-step and the step-number. This function is relatively
 * straightforward.
 * 

 * 
 * 
 * @code
 *   template <int dim>
 *   void ParticleTracking<dim>::run()
 *   {
 *     DiscreteTime discrete_time(0, par.final_time, par.time_step);
 * 
 *     generate_particles();
 * 
 *     pcout << "Repartitioning triangulation after particle generation"
 *           << std::endl;
 *     background_triangulation.repartition();
 * 
 * @endcode
 * 
 * We set the initial property of the particles by doing an
 * explicit Euler iteration with a time-step of 0 both in the case
 * of the analytical and the interpolated approach.
 * 
 * @code
 *     if (interpolated_velocity)
 *       {
 *         setup_background_dofs();
 *         interpolate_function_to_field();
 *         euler_step_interpolated(0.);
 *       }
 *     else
 *       euler_step_analytical(0.);
 * 
 *     output_particles(discrete_time.get_step_number());
 *     if (interpolated_velocity)
 *       output_background(discrete_time.get_step_number());
 * 
 * @endcode
 * 
 * The particles are advected by looping over time.
 * 
 * @code
 *     while (!discrete_time.is_at_end())
 *       {
 *         discrete_time.advance_time();
 *         velocity.set_time(discrete_time.get_previous_time());
 * 
 *         if ((discrete_time.get_step_number() % par.repartition_frequency) == 0)
 *           {
 *             background_triangulation.repartition();
 *             if (interpolated_velocity)
 *               setup_background_dofs();
 *           }
 * 
 *         if (interpolated_velocity)
 *           {
 *             interpolate_function_to_field();
 *             euler_step_interpolated(discrete_time.get_previous_step_size());
 *           }
 *         else
 *           euler_step_analytical(discrete_time.get_previous_step_size());
 * 
 * @endcode
 * 
 * After the particles have been moved, it is necessary to identify
 * in which cell they now reside. This is achieved by calling
 * <code>sort_particles_into_subdomains_and_cells</code>
 * 
 * @code
 *         particle_handler.sort_particles_into_subdomains_and_cells();
 * 
 *         if ((discrete_time.get_step_number() % par.output_frequency) == 0)
 *           {
 *             output_particles(discrete_time.get_step_number());
 *             if (interpolated_velocity)
 *               output_background(discrete_time.get_step_number());
 *           }
 *       }
 *   }
 * 
 * } // namespace Step68
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * The remainder of the code, the `main()` function, is standard.
 * We note that we run the particle tracking with the analytical velocity
 * and the interpolated velocity and produce both results
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   using namespace Step68;
 *   using namespace dealii;
 *   deallog.depth_console(1);
 * 
 *   try
 *     {
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       std::string prm_file;
 *       if (argc > 1)
 *         prm_file = argv[1];
 *       else
 *         prm_file = "parameters.prm";
 * 
 *       ParticleTrackingParameters par;
 *       ParameterAcceptor::initialize(prm_file);
 *       {
 *         Step68::ParticleTracking<2> particle_tracking(par, false);
 *         particle_tracking.run();
 *       }
 *       {
 *         Step68::ParticleTracking<2> particle_tracking(par, true);
 *         particle_tracking.run();
 *       }
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
examples/step-68/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序的目录中默认包含一个参数文件的例子。如果你没有在命令行中指定参数文件作为参数，程序将默认尝试读取文件 "参数.prm"，并执行代码。

在任何数量的核心上，模拟输出将看起来像。

@code
bash$ mpirun -np 4 ./step-68 parameters.prm
Number of particles inserted: 606
Repartitioning triangulation after particle generation
Writing particle output file: analytical-particles-0
Writing particle output file: analytical-particles-10
Writing particle output file: analytical-particles-20
Writing particle output file: analytical-particles-30
...
Number of particles inserted: 606
Repartitioning triangulation after particle generation
Writing particle output file: interpolated-particles-0
Writing background field file: background-0
Writing particle output file: interpolated-particles-10
Writing background field file: background-10
Writing particle output file: interpolated-particles-20
Writing background field file: background-20
Writing particle output file: interpolated-particles-30
Writing background field file: background-30
...
Writing particle output file: interpolated-particles-1980
Writing background field file: background-1980
Writing particle output file: interpolated-particles-1990
Writing background field file: background-1990
Writing particle output file: interpolated-particles-2000
Writing background field file: background-2000
@endcode



我们注意到，在默认情况下，模拟以分析速度运行粒子跟踪2000次，然后从头开始，以相同的时间运行速度插值的粒子跟踪。每隔10次迭代就会写出结果。

<a name="Motionoftheparticles"></a><h3> Motion of the particles </h3>


下面的动画显示了粒子在被流场吸引时的轨迹。我们看到，在流动的整个过程中，粒子又回到了它们的初始配置，这是预料之中的事。

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/EbgS5Ch35Xs"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly



<a name="Dynamicloadbalancing"></a><h3> Dynamic load balancing </h3>


下面的动画显示了动态负载平衡的影响。我们清楚地看到，子域自我调整以平衡每个子域的粒子数量。然而，完美的负载平衡并没有达到，部分原因是由于背景网格的粗糙性。

@htmlonly
<p align="center">
  <iframe width="560" height="500" src="https://www.youtube.com/embed/ubUcsR4ECj4"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly






<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序强调了在deal.II中处理粒子的一些主要能力，特别是它们用于分布式并行模拟的能力。然而，这一步骤可以以多种方式进行扩展。

- 可以使用高阶时间积分（例如使用Runge-Kutta 4方法）来提高精度，或者在相同精度下允许更大的时间步长。

- 完整的运动方程（含惯性）可以为粒子求解。在这种情况下，粒子将需要有额外的属性，例如它们的质量，如步骤19，如果想考虑与流体的相互作用，还需要考虑它们的直径。

- 耦合到流动求解器。这一步可以直接与任何解决斯托克斯方程（第32步，第70步）或纳维-斯托克斯方程（例如第57步）的并行程序相耦合。

- 计算两个模型之间最终粒子位置的差异，可以量化内插误差对粒子运动的影响。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-68.cc"
*/
