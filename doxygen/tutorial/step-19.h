/**
@page step_19 The step-19 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Spatialdiscretization">Spatial discretization</a>
        <li><a href="#Dealingwithparticlesprogrammatically">Dealing with particles programmatically</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Globaldefinitions">Global definitions</a>
        <li><a href="#Themainclass">The main class</a>
        <li><a href="#ThecodeCathodeRaySimulatorcodeclassimplementation">The <code>CathodeRaySimulator</code> class implementation</a>
      <ul>
        <li><a href="#ThecodeCathodeRaySimulatorcodeconstructor">The <code>CathodeRaySimulator</code> constructor</a>
        <li><a href="#ThecodeCathodeRaySimulatormake_gridcodefunction">The <code>CathodeRaySimulator::make_grid</code> function</a>
        <li><a href="#ThecodeCathodeRaySimulatorsetup_systemcodefunction">The <code>CathodeRaySimulator::setup_system</code> function</a>
        <li><a href="#ThecodeCathodeRaySimulatorassemble_systemcodefunction">The <code>CathodeRaySimulator::assemble_system</code> function</a>
        <li><a href="#CathodeRaySimulatorsolve">CathodeRaySimulator::solve</a>
        <li><a href="#CathodeRaySimulatorrefine_grid">CathodeRaySimulator::refine_grid</a>
        <li><a href="#CathodeRaySimulatorcreate_particles">CathodeRaySimulator::create_particles</a>
        <li><a href="#CathodeRaySimulatormove_particles">CathodeRaySimulator::move_particles</a>
        <li><a href="#CathodeRaySimulatortrack_lost_particle">CathodeRaySimulator::track_lost_particle</a>
        <li><a href="#CathodeRaySimulatorupdate_timestep_size">CathodeRaySimulator::update_timestep_size</a>
        <li><a href="#ThecodeCathodeRaySimulatoroutput_resultscodefunction">The <code>CathodeRaySimulator::output_results()</code> function</a>
        <li><a href="#CathodeRaySimulatorrun">CathodeRaySimulator::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Avoidingaperformancebottleneckwithparticles"> Avoiding a performance bottleneck with particles </a>
        <li><a href="#Morestatisticsaboutelectrons"> More statistics about electrons </a>
        <li><a href="#Abettersynchronizedvisualization"> A better-synchronized visualization </a>
        <li><a href="#Abettertimeintegrator"> A better time integrator </a>
        <li><a href="#Parallelization"> Parallelization </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-19/doc/intro.dox



 <br> 

<i>
This program was contributed by Wolfgang Bangerth, Rene Gassmoeller, and Peter Munch.


Wolfgang Bangerth acknowledges support through NSF
awards DMS-1821210, EAR-1550901, and OAC-1835673.
</i>

 @note  deal.II中存在对粒子的支持，这主要是由于Rene Gassmoeller的最初努力。如果你在自己的工作中使用粒子功能，请引用出版物 @cite GLHPW2018 来确认这项工作。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


一般来说，有限元方法，特别是deal.II，是为了解决偏微分方程而发明的--换句话说，是为了解决[连续体力学](https://en.wikipedia.org/wiki/Continuum_mechanics)问题。另一方面，有时人们想解决的问题是，跟踪单个物体（"粒子"）以及它们的位置如何演变是有用的。如果这只是导致一组常微分方程，例如，如果你想跟踪太阳系中行星随时间变化的位置，那么deal.II显然不是你合适的工具。另一方面，如果这种演变是由于与偏微分方程的解的相互作用，或者有一个网格来确定哪些粒子与其他粒子相互作用（如在[平滑粒子流体力学（SPH）](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)方法中），那么deal.II对你有支持。

我们在这里要考虑的情况是带电粒子如何在电场中移动。作为动力，我们将考虑[阴极射线]（https://en.wikipedia.org/wiki/Cathode_ray）。由一块被加热的带负电的金属（"阴极"）发出的电子，然后被电场加速到带正电的电极（"阳极"）。阳极通常是环形的，这样大部分电子可以以电子束的形式飞过孔。在过去，它们可能会照亮由[阴极射线管](https://en.wikipedia.org/wiki/Cathode-ray_tube)制成的电视的屏幕。今天，电子束反而在[X射线机](https://en.wikipedia.org/wiki/X-ray_tube)、[电子束光刻](https://en.wikipedia.org/wiki/Electron-beam_lithography)、[电子束焊接](https://en.wikipedia.org/wiki/Electron-beam_welding)和其他一些领域发挥了作用。

然后我们要考虑的方程如下。首先，我们需要描述电场。通过注意到电势 $V$ 满足方程，这是最容易完成的。

@f[


  -\epsilon_0 \Delta V = \rho


@f]

其中 $\epsilon_0$ 是真空的介电常数，而 $\rho$ 是电荷密度。这是由我们将选择的边界条件所增强的，如下所示。

@f{align*}{
  V &= -V_0 && \text{on}\; \Gamma_\text{cathode}\subset\partial\Omega \\
  V &= +V_0 && \text{on}\; \Gamma_\text{anode}\subset\partial\Omega \\
  \epsilon\frac{\partial V}{\partial n} &= 0
   && \text{on}\; \partial\Omega\setminus\Gamma_\text{cathode}\setminus\Gamma_\text{anode}.


@f}

换句话说，我们在两个电极上规定电压 $+V_0$ 和 $-V_0$ ，在其他地方规定绝缘（诺伊曼）边界条件。由于粒子的动力学纯粹是由于电场 $\mathbf E=\nabla V$ ，我们也可以在两个电极上规定 $2V_0$ 和 $0$ --所有重要的是两个电极的电压差。

考虑到这个电势 $V$ 和电场 $\mathbf E=\nabla V$ ，我们可以用微分方程来描述 $i$ 这个粒子的轨迹

@f[
  m {\ddot {\mathbf x}}_i = e\mathbf E,


@f]

其中 $m,e$ 是每个粒子的质量和电荷。在实践中，将其写成位置 $\mathbf x$ 和速度 $\mathbf v$ 的一阶微分方程系统很方便。

@f{align*}{
  {\dot {\mathbf v}}_i &= \frac{e\mathbf E}{m}, \\
  {\dot {\mathbf x}}_i &= {\mathbf v}_i.


@f}

我们将用来处理粒子的deal.II类， Particles::ParticleHandler, 以一种方式存储粒子，因此位置 $\mathbf x_i$ 是 Particles::ParticleHandler 数据结构的一部分。它存储的粒子是按它们所在的单元分类的，因此需要知道每个粒子的位置）。另一方面，速度 $\mathbf v_i$ 与 Particles::ParticleHandler 无关，因此我们将把它存储为每个粒子的 "属性"，并在每个时间步长中更新。属性也可以用来存储我们可能关心的关于每个粒子的任何其他数量：它的电荷，或者如果它们大于一个电子，它的颜色、质量、在空间的位置、化学成分等等。

要完成这个模型，还有两件事要讨论。粒子从哪里开始以及电荷密度 $\rho$ 是什么。

首先，在历史上，阴极射线使用非常大的电场将电子从金属中拉出来。这只产生一个相对较小的电流。我们可以通过加热阴极来做得更好：在这种情况下，统计学上的一部分电子有足够的热能来离开金属；然后电场只要足够强，就可以把它们从宿主的吸引中拉出来。我们将以下列方式对此进行建模。如果（i）电场指向远离电极，即如果 $\mathbf E \cdot \mathbf n < 0$ ，其中 $\mathbf n$ 是指向域外（进入电极）的面的法向量，以及（ii）电场超过一个阈值 $|\mathbf E|\ge E_\text{threshold}$ ，我们将创建一个新粒子。这肯定不是真正发生的足够精确的模型，但对于我们目前的教程程序来说已经足够好了。

第二，原则上我们必须通过以下方式建立电荷密度模型

@f[
  \rho(\mathbf x) = \sum_i e\delta(\mathbf x-\mathbf x_i).


@f]



 @note 现在的问题是，在现实中，一台老式电视中的阴极射线管产生的电流大约为几毫安培。在粒子加速器的更高能量的光束中，电流可能只有几纳安培。但一个安培是每秒流动的 $6\times 10^{18}$ 个电子。现在，正如你将在结果部分看到的，我们实际上只模拟了几微秒（ $10^{-5}$ 秒），但这仍然导致非常非常多的电子 -- 远远超过我们希望用像目前这样小的程序来模拟。因此，让我们假设每个粒子代表 $N$ 个电子。那么粒子的质量和电荷也是 $Nm$ 和 $Ne$ ，我们要解决的方程式是

@f[
  (Nm) {\ddot {\mathbf x}}_i = (Ne)\mathbf E,


@f]

当然，这与上述情况完全相同。另一方面，这些电子 "团块 "的电荷密度由以下公式给出

@f[
  \rho(\mathbf x) = \sum_i (Ne)\delta(\mathbf x-\mathbf x_i).


@f]

我们将在程序中实现这种形式，其中 $N$ 在程序中被选得相当大，以确保粒子实际影响电场。这在实践中可能并不现实。在大多数情况下，没有足够的电子来实际影响整个电场。但现实主义不是我们的目标）。)




 @note  人们可能会问，为什么电场（或者说，电势）的方程没有时间导数，而电子位置的方程却有。从本质上讲，这是一个建模假设。我们假设粒子移动得很慢，以至于在任何时候电场都处于平衡状态。这就是说，换句话说，电子的速度远远小于光速。换句话说，我们可以用电极电压来重新表述  $V_0$  ：由于每伏特的电动势都会使电子加速约600公里/秒（忽略相对论效应），要求  $|\mathbf v_i\|\ll c$  等于说  $2V_0 \ll 500 \text{V}$  。在这个假设下（以及电子总数很小的假设），我们也可以忽略移动电荷产生的磁场，否则也会影响电子的运动。




<a name="Timediscretization"></a><h3>Time discretization</h3>


上面概述的方程形成了一组耦合微分方程。让我们再次把它们集中在一起，以明确这一点。

@f{align*}{


  -\epsilon_0 \Delta V &= \sum_i e\delta(\mathbf x-\mathbf x_i)
  \\
  {\dot {\mathbf x}}_i &= {\mathbf v}_i,
  \\
  {\dot {\mathbf v}}_i &= \frac{e\mathbf E}{m} = \frac{e\mathbf \nabla V}{m}.


@f}

由于电势对粒子位置的依赖性很强，我们不想将其作为一个耦合系统来求解，而是采用一种解耦的方法，首先求解每个时间步长的电势，然后再求解粒子的位置。这与我们在第21步、第31步和第32步（仅举几例）所做的工作的精神是一样的，都可以在第58步讨论的算子分割方法的背景下加以理解。

因此，如果我们用大指数 $n$ 表示时间步长，并且如果我们对ODE使用简单的时间离散化，那么这意味着我们必须在每个时间步长中解决以下方程组。

@f{align*}{


  -\epsilon_0 \Delta V^{(n)} &= \sum_i e\delta(\mathbf x-\mathbf x_i^{(n-1)})
  \\
  \frac{{\mathbf v}_i^{(n)}-{\mathbf v}_i^{(n-1)}}{\Delta t} &= \frac{e\nabla V^{(n)}}{m}
  \\
  \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} &= {\mathbf v}_i^{(n)}.


@f}

当然还有许多更好的方法来做时间离散化（例如简单的[跃迁方案](https://en.wikipedia.org/wiki/Leapfrog_integration)），但这不是本教程程序的重点，因此我们将满足于这里的内容。不过，我们将在本程序的<a href="#extensions">possibilities for extensions</a>部分对这个难题的一个部分进行评论）。

还有一个问题是我们应该如何选择时间步长  $\Delta t$  。这里的限制是， Particles::ParticleHandler 类需要跟踪每个粒子在哪个单元中。如果我们平行运行计算（比如，在step-70中），这尤其是一个问题，因为在这种情况下，每个进程只存储它拥有的那些单元，再加上一层 "幽灵单元"。这在这里并不重要，但一般来说，我们应该确保在每个时间步长中，一个粒子只从一个单元移动到它的任何一个近邻（面、边或顶点的邻居）。如果我们能确保这一点，那么 Particles::ParticleHandler 就能保证能够找出粒子最后在哪个单元。为了做到这一点，一个有用的经验法则是，我们应该选择时间步长，使所有粒子的预期移动距离小于一个细胞的直径。

@f[
  \Delta t \le \frac{h_i}{\|\mathbf v_i\|} \qquad\qquad \forall i,


@f]

或等价的

@f[
  \Delta t \le \min_i \frac{h_i}{\|\mathbf v_i\|}.


@f]

这里， $h_i$ 是粒子 $i$ 所在的单元格最短边的长度--本质上是对单元格大小的衡量。

另一方面，一个粒子可能已经在一个单元的边界上，而邻近的单元可能已经进一步细化。因此，那么穿过那个*邻近*单元的时间实际上将是上述数量的一半，这表明

@f[
  \Delta t \le \min_i \frac{\tfrac 12 h_i}{\|\mathbf v_i\|}.


@f]



但即使这样也是不够的。上面的公式在每次更新粒子位置时使用的是

@f[
\frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} = {\mathbf v}_i^{(n)},


@f]

也就是说，使用当前的*速度 ${\mathbf v}_i^{n}$ 。但是当我们需要选择 $\Delta t$ 时，我们还没有当前的速度 -- 也就是在我们更新了潜能 $V^{(n)}$ 之后，但在我们将速度从 ${\mathbf v}_i^{(n-1)}$ 更新到 ${\mathbf v}_i^{(n)}$ 之前。我们有的只是  ${\mathbf v}_i^{(n-1)}$  。所以我们需要一个额外的安全系数来实现我们的最终选择。

@f[
  \Delta t^{(n)} =
  c_\text{safety} \min_i \frac{\tfrac 12 h_i}{\|\mathbf v_i^{(n-1)}\|}.


@f]

 $c_\text{safety}$ 应该有多大？这取决于与 $\|\mathbf v_i^{(n)}\|$ 相比， $\|\mathbf v_i^{(n-1)}\|$ 可能被低估了多少，而这实际上是很容易评估的。如果沿途遇到的电场大致恒定，那么在一个时间步长中产生的速度为零的粒子，在每个连续的时间步长中大致会获得相等的速度增量。因此， $\|\mathbf v_i^{(n-1)}\|$ 和 $\|\mathbf v_i^{(n)}\|$ 之间的最大差异将是一个系数。因此，我们将选择 $c_\text{safety}=0.5$  。

我们应该考虑的只有另外一种情况。在第一个时间步骤中会发生什么？在那里，任何要被移动的粒子刚刚被创造出来，但它们的速度是零。所以我们不知道我们应该为它们选择什么速度。当然，在所有其他时间步骤中，也有刚刚被创造出来的粒子，但一般来说，具有最高速度的粒子限制了时间步骤的大小，因此新创造出来的具有零速度的粒子并不重要。但是如果我们**只有这样的粒子？

在这种情况下，我们可以使用以下近似值。如果一个粒子从 $\mathbf v^{(0)}=0$ 开始，那么更新公式告诉我们

@f[
  {\mathbf v}_i^{(1)} = \frac{e\nabla V^{(1)}}{m} \Delta t,


@f]

因此

@f[
    \frac{{\mathbf x}_i^{(1)}-{\mathbf x}_i^{(0)}}{\Delta t} = {\mathbf v}_i^{(1)},


@f]

我们可以把它写成

@f[
    {\mathbf x}_i^{(1)} - {\mathbf x}_i^{(0)} = \frac{e\nabla V^{(1)}}{m} \Delta t^2.


@f]

不想让一个粒子移动超过 $\frac 12 h_i$ ，那么就意味着我们应该选择时间步长为

@f[
  \Delta t
  \le
  \min_i
  \sqrt{ \frac{h_i m}{e \|\nabla V^{(1)}\| }}.


@f]

使用关于相邻单元可能小2倍的相同论点，然后得出时间步长为0的最终公式。

@f[
  \Delta t
  =
  \min_i
  \sqrt{ \frac{\frac 12 h_i m}{e \|\nabla V^{(1)}\| } }.


@f]



严格来说，我们必须在每个粒子的位置评估电势 $V^{(1)}$ ，但一个足够好的近似值是使用各自单元顶点的最大值。为什么是顶点而不是中点？因为拉普拉斯方程的解的梯度，即电场，在位于单元顶点的角落奇点上是最大的）。)这样做的好处是，我们可以很好地利用FEValues功能，只要各单元的正交点相同，就可以循环使用预计算的材料。

我们总是可以运行这种方案来估计 $\mathbf v_i^{(n-1)}$ 和 $\mathbf v_i^{(n)}$ 之间的差异，但它依赖于评估每个单元的电场 $\mathbf E$ ，这很昂贵。因此，我们将把这种方法限制在第一个时间步骤上。




<a name="Spatialdiscretization"></a><h3>Spatial discretization</h3>


在讨论了时间离散化之后，对空间离散化的讨论将很简短：我们使用二次有限元，即空间 $Q_2$  ，来近似计算电动势 $V$  。在初始时间步骤中，网格被调整了几次。如果你读过第6步，所有这些都是完全标准的，而且实现起来也没有规定任何形式的惊喜。




<a name="Dealingwithparticlesprogrammatically"></a><h3>Dealing with particles programmatically</h3>


实际上，在deal.II中，添加和移动粒子并不十分困难。要添加一个粒子，本程序的`create_particles()`函数只需使用以下形式的代码片段。

@code
  Particles::Particle<dim> new_particle;
  new_particle.set_location(location);
  new_particle.set_reference_location
      (mapping.transform_real_to_unit_cell(cell, location));
  new_particle.set_id(n_current_particles);


  particle_handler.insert_particle(new_particle, cell);
@endcode

换句话说，它与在 `std::set` 或 `std::map`: 中插入一个对象没有什么不同。 创建对象，设置其属性（这里是当前位置、其参考单元位置和其id）并调用`insert_particle`。唯一可能令人惊讶的是参考位置。为了评估诸如  $\nabla V(\mathbf x_i)$  的东西，有必要在位置  $\mathbf x_i$  评估有限元场。但这需要在参考单元 $\hat{\mathbf x}_i$ 上的点评估有限元形状函数。为了使之有效，每个粒子不仅要存储它的位置和它所在的单元，还要存储该点在单元参考坐标系中对应的位置。

这样，更新粒子的位置就不再困难了。我们只需要调用

@code
  particle->set_location(new_location);
@endcode

我们在`move_particles()`函数中这样做。唯一的区别是，我们必须告诉 Particles::ParticleHandler 类也要找到该位置所对应的单元（而且，在并行计算时，哪个进程拥有该单元）。出于效率的考虑，这在更新所有粒子的位置后最容易完成，并通过 Particles::ParticleHandler::sort_particles_into_subdomains_and_cells() 函数实现。

当然，有些时候，粒子可能会离开有关的域。在这种情况下， Particles::ParticleHandler::sort_particles_into_subdomains_and_cells() 不能找到周围的单元，而只是简单地删除该粒子。但是，跟踪以这种方式丢失的粒子的数量往往是有用的，为此， Particles::ParticleHandler 类提供了一个可以附加的 "信号"。我们在主类的构造函数中展示了如何做到这一点，以计算每个时间步骤中损失了多少粒子。具体来说，这种工作方式是， Particles::ParticleHandler 类有一个 "信号"，人们可以附加一个函数，只要信号被触发就会执行。在这里，这看起来如下。

@code
    particle_handler.signals.particle_lost.connect(
      [this](const typename Particles::ParticleIterator<dim> &        particle,
             const typename Triangulation<dim>::active_cell_iterator &cell)
      {
        this->track_lost_particle(particle, cell);
      });
@endcode

这有点拗口，但实际情况是这样的。我们声明了一个 "捕获"`this`指针的lambda函数（这样我们就可以在lambda函数中访问周围对象的成员函数），它需要两个参数。

- 指的是已经 "丢失 "的粒子。

- 它最后所在的单元格的引用。这个lambda函数然后简单地用这些参数调用 `CathodeRaySimulator::track_lost_particle` 函数。当我们把这个lambda函数附加到信号上时， Particles::ParticleHandler::sort_particles_into_subdomains_and_cells() 函数将为每个找不到新家的粒子触发信号。这让我们有机会记录下粒子的位置，并记录下关于它的统计数据。




 @note  在这个教程程序中，我们通过手工插入粒子，并在我们根据包括静电问题的解决的条件专门选择的位置插入粒子。但在其他情况下，人们主要希望将粒子作为被动对象使用，例如，追踪和可视化流体流动问题的流场。在这些情况下， Particles::Generators 命名空间中有许多函数可以自动生成粒子。例如，这个命名空间中的一个函数也被用于step-70教程程序中。




<a name="Thetestcase"></a><h3>The test case</h3>


这里的测试案例并不意味着是对阴极射线管的真实描述，但它具有正确的一般特征，而且在任何情况下，重点只是演示如何实现使用粒子的deal.II代码。

下图显示了我们要使用的几何图形。

<p align="center"> <img src="https://www.dealii.org/images/steps/developer/step-19.geometry.png" alt="本程序中使用的几何图形" width="600">  </p> 

在这幅图中，边界上用红色和蓝色标记的部分是阴极，保持在一个电动势 $V=-V_0$ 。阴极的红色部分是被加热的部分，导致电子离开金属，然后被电场加速（也显示了一些电场线）。边界的绿色部分是阳极，保持在 $V=+V_0$ 。边界的其余部分满足诺伊曼边界条件。

这种设置模仿了真实的设备。重心角导致电势 $V$ ，其导数（电场 $\mathbf E$ ）有一个奇点--换句话说，它在角的附近变得非常大，允许它把电子从金属中扯出来。这些电子然后被加速推向（绿色）阳极，阳极中间有一个孔，电子可以通过这个孔逃离设备并飞到屏幕上，在那里它们激发 "荧光粉"，然后发出我们从这些老式电视屏幕上看到的光。阴极的非加热部分不受电子发射的影响--在代码中，我们将其标记为电子管的 "聚焦元件"，因为它的负电压会排斥电子，并确保它们不只是垂直于边界从阴极的加热部分飞走，而是事实上将它们的路径弯曲到右边的阳极。

图中的电场线也说明了电场分别连接着负极和正极。电子经历的加速力是沿着这些场线的。最后，图片显示了计算中使用的网格，说明在重租角的顶端以及边界条件改变的所有地方都有奇异点；这些奇异点是可见的，因为网格在这些地方被细化。

实际的利益是要弄清楚从阴极发射的电子中有哪一部分真正通过了阳极上的孔--那些只是反弹到阳极本身的电子除了将电转化为热之外，实际上并没有什么用处。因此，在`track_lost_particle()`函数中（为每个离开域的粒子调用，见上文），我们将估计它可能离开域的位置并在输出中报告。




 @note 值得重申的是，这里使用的几何图形，以及事实上这个程序的任何其他方面，都不是为了代表任何半点现实的东西。教程是我们教授deal.II如何工作的工具，我们经常使用我们有某种直觉的情况，因为这有助于我们解释程序的输出，但这就是我们打算让程序除了作为教学工具之外做任何有用的事情的程度。


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
 * The majority of the include files used in this program are
 * well known from step-6 and similar programs:
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/fe/mapping_q.h>
 * #include <deal.II/matrix_free/fe_point_evaluation.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * 
 * @endcode
 * 
 * The ones that are new are only the following three: The first declares the
 * DiscreteTime class that helps us keep track of time in a time-dependent
 * simulation. The latter two provide all of the particle functionality,
 * namely a way to keep track of particles located on a mesh (the
 * Particles::ParticleHandler class) and the ability to output these
 * particles' locations and their properties for the purposes of
 * visualization (the Particles::DataOut class).
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/particles/particle_handler.h>
 * #include <deal.II/particles/data_out.h>
 * 
 * #include <fstream>
 * 
 * using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Globaldefinitions"></a> 
 * <h3>Global definitions</h3>
 * 

 * 
 * As is customary, we put everything that corresponds to the details of the
 * program into a namespace of its own. At the top, we define a few constants
 * for which we would rather use symbolic names than hard-coded numbers.
 * 

 * 
 * Specifically, we define numbers for
 * @ref GlossBoundaryIndicator "boundary indicators"
 * for the various parts of the geometry, as well as the physical properties
 * of electrons and other specifics of the setup we use here.
 * 

 * 
 * For the boundary indicators, let us start enumerating at some
 * random value 101. The principle here is to use numbers that are
 * *uncommon*. If there are pre-defined boundary indicators previously
 * set by the `GridGenerator` functions, they will likely be small
 * integers starting from zero, but not in this rather randomly chosen
 * range. Using numbers such as those below avoids the possibility for
 * conflicts, and also reduces the temptation to just spell these
 * numbers out in the program (because you will probably never
 * remember which is which, whereas you might have been tempted if
 * they had started at 0).
 * 
 * @code
 * namespace Step19
 * {
 *   namespace BoundaryIds
 *   {
 *     constexpr types::boundary_id open          = 101;
 *     constexpr types::boundary_id cathode       = 102;
 *     constexpr types::boundary_id focus_element = 103;
 *     constexpr types::boundary_id anode         = 104;
 *   } // namespace BoundaryIds
 * 
 *   namespace Constants
 *   {
 *     constexpr double electron_mass   = 9.1093837015e-31;
 *     constexpr double electron_charge = 1.602176634e-19;
 * 
 *     constexpr double V0 = 1;
 * 
 *     constexpr double E_threshold = 0.05;
 * 
 *     constexpr double electrons_per_particle = 3e15;
 *   } // namespace Constants
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * The following is then the main class of this program. It has,
 * fundamentally, the same structure as step-6 and many other
 * tutorial programs. This includes the majority of the member
 * functions (with the purpose of the rest probably self-explanatory
 * from their names) as well as only a small number of member
 * variables beyond those of step-6, all of which are related to
 * dealing with particles.
 * 
 * @code
 *   template <int dim>
 *   class CathodeRaySimulator
 *   {
 *   public:
 *     CathodeRaySimulator();
 * 
 *     void run();
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void assemble_system();
 *     void solve_field();
 *     void refine_grid();
 * 
 *     void create_particles();
 *     void move_particles();
 *     void track_lost_particle(
 *       const typename Particles::ParticleIterator<dim> &        particle,
 *       const typename Triangulation<dim>::active_cell_iterator &cell);
 * 
 * 
 *     void update_timestep_size();
 *     void output_results() const;
 * 
 *     Triangulation<dim>        triangulation;
 *     MappingQGeneric<dim>      mapping;
 *     FE_Q<dim>                 fe;
 *     DoFHandler<dim>           dof_handler;
 *     AffineConstraints<double> constraints;
 * 
 *     SparseMatrix<double> system_matrix;
 *     SparsityPattern      sparsity_pattern;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     Particles::ParticleHandler<dim> particle_handler;
 *     types::particle_index           next_unused_particle_id;
 *     types::particle_index           n_recently_lost_particles;
 *     types::particle_index           n_total_lost_particles;
 *     types::particle_index           n_particles_lost_through_anode;
 * 
 *     DiscreteTime time;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorcodeclassimplementation"></a> 
 * <h3>The <code>CathodeRaySimulator</code> class implementation</h3>
 * 

 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorcodeconstructor"></a> 
 * <h4>The <code>CathodeRaySimulator</code> constructor</h4>
 * 

 * 
 * So then let us get started on the implementation. What the constructor
 * does is really only a straight-forward initialization of all of the member
 * variables at the top. The only two worth mentioning are the
 * `particle_handler`, which is handed a reference to the triangulation
 * on which the particles will live (currently of course still empty,
 * but the particle handler stores the reference and will use it once
 * particles are added -- which happens after the triangulation is built).
 * The other piece of information it gets is how many "properties"
 * each particle needs to store. Here, all we need each particle to
 * remember is its current velocity, i.e., a vector with `dim`
 * components. There are, however, other intrinsic properties that
 * each particle has and that the Particles::ParticleHandler class
 * automatically and always makes sure are available; in particular,
 * these are the current location of a particle, the cell it is on,
 * it's reference location within that cell, and the particle's ID.
 *   

 * 
 * The only other variable of interest is `time`, an object of type
 * DiscreteTime. It keeps track of the current time we are in a
 * time-dependent simulation, and is initialized with the start time
 * (zero) and end time ($10^{-4}$). We will later set the time step
 * size in `update_timestep_size()`.
 *   

 * 
 * The body of the constructor consists of a piece of code we have
 * already discussed in the introduction. Namely, we make sure that the
 * `track_lost_particle()` function is called by the `particle_handler`
 * object every time a particle leaves the domain.
 * 
 * @code
 *   template <int dim>
 *   CathodeRaySimulator<dim>::CathodeRaySimulator()
 *     : mapping(1)
 *     , fe(2)
 *     , dof_handler(triangulation)
 *     , particle_handler(triangulation, mapping, /*n_properties=*/dim)
 *     , next_unused_particle_id(0)
 *     , n_recently_lost_particles(0)
 *     , n_total_lost_particles(0)
 *     , n_particles_lost_through_anode(0)
 *     , time(0, 1e-4)
 *   {
 *     particle_handler.signals.particle_lost.connect(
 *       [this](const typename Particles::ParticleIterator<dim> &        particle,
 *              const typename Triangulation<dim>::active_cell_iterator &cell) {
 *         this->track_lost_particle(particle, cell);
 *       });
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatormake_gridcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::make_grid</code> function</h4>
 * 

 * 
 * The next function is then responsible for generating the mesh on which
 * we want to solve. Recall how the domain looks like:
 * <p align="center">
 * <img
 * src="https://www.dealii.org/images/steps/developer/step-19.geometry.png"
 * alt="The geometry used in this program"
 * width="600">
 * </p>
 * We subdivide this geometry into a mesh of $4\times 2$ cells that looks
 * like this:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *   *---*---*---*---*
 *   \   |   |   |   |
 *    *--*---*---*---*
 *   /   |   |   |   |
 *   *---*---*---*---*
 * @endcode
 * </div>
 * The way this is done is by first defining where the $15=5\times 3$
 * vertices are located -- here, we say that they are on integer points
 * with the middle one on the left side moved to the right by a value of
 * `delta=0.5`.
 *   

 * 
 * In the following, we then have to say which vertices together form
 * the 8 cells. The following code is then entirely equivalent to what
 * we also do in step-14:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::make_grid()
 *   {
 *     static_assert(dim == 2,
 *                   "This function is currently only implemented for 2d.");
 * 
 *     const double       delta = 0.5;
 *     const unsigned int nx    = 5;
 *     const unsigned int ny    = 3;
 * 
 *     const std::vector<Point<dim>> vertices 
 *       = {{0, 0},
 *          {1, 0},
 *          {2, 0},
 *          {3, 0},
 *          {4, 0},
 *          {delta, 1},
 *          {1, 1},
 *          {2, 1},
 *          {3, 1},
 *          {4, 1},
 *          {0, 2},
 *          {1, 2},
 *          {2, 2},
 *          {3, 2},
 *          {4, 2}};
 *     AssertDimension(vertices.size(), nx * ny);
 * 
 *     const std::vector<unsigned int> cell_vertices[(nx - 1) * (ny - 1)] = {
 *       {0, 1, nx + 0, nx + 1},
 *       {1, 2, nx + 1, nx + 2},
 *       {2, 3, nx + 2, nx + 3},
 *       {3, 4, nx + 3, nx + 4},
 * 
 *       {5, nx + 1, 2 * nx + 0, 2 * nx + 1},
 *       {nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2},
 *       {nx + 2, nx + 3, 2 * nx + 2, 2 * nx + 3},
 *       {nx + 3, nx + 4, 2 * nx + 3, 2 * nx + 4}};
 * 
 * @endcode
 * 
 * With these arrays out of the way, we can move to slightly higher
 * higher-level data structures. We create a vector of CellData
 * objects that store for each cell to be created the vertices in
 * question as well as the @ref GlossMaterialId "material id" (which
 * we will here simply set to zero since we don't use it in the program).
 *     

 * 
 * This information is then handed to the
 * Triangulation::create_triangulation() function, and the mesh is twice
 * globally refined.
 * 
 * @code
 *     std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>());
 *     for (unsigned int i = 0; i < cells.size(); ++i)
 *       {
 *         cells[i].vertices    = cell_vertices[i];
 *         cells[i].material_id = 0;
 *       }
 * 
 *     triangulation.create_triangulation(
 *       vertices,
 *       cells,
 *       SubCellData()); // No boundary information
 * 
 *     triangulation.refine_global(2);
 * 
 * @endcode
 * 
 * The remaining part of the function loops over all cells and their faces,
 * and if a face is at the boundary determines which boundary indicator
 * should be applied to it. The various conditions should make sense if
 * you compare the code with the picture of the geometry above.
 *     

 * 
 * Once done with this step, we refine the mesh once more globally.
 * 
 * @code
 *     for (auto &cell : triangulation.active_cell_iterators())
 *       for (auto &face : cell->face_iterators())
 *         if (face->at_boundary())
 *           {
 *             if ((face->center()[0] > 0) && (face->center()[0] < 0.5) &&
 *                 (face->center()[1] > 0) && (face->center()[1] < 2))
 *               face->set_boundary_id(BoundaryIds::cathode);
 *             else if ((face->center()[0] > 0) && (face->center()[0] < 2))
 *               face->set_boundary_id(BoundaryIds::focus_element);
 *             else if ((face->center()[0] > 4 - 1e-12) &&
 *                      ((face->center()[1] > 1.5) || (face->center()[1] < 0.5)))
 *               face->set_boundary_id(BoundaryIds::anode);
 *             else
 *               face->set_boundary_id(BoundaryIds::open);
 *           }
 * 
 *     triangulation.refine_global(1);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatorsetup_systemcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::setup_system</code> function</h4>
 * 

 * 
 * The next function in this program deals with setting up the various
 * objects related to solving the partial differential equations. It is
 * in essence a copy of the corresponding function in step-6 and requires
 * no further discussion.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::cathode,
 *                                              Functions::ConstantFunction<dim>(
 *                                                -Constants::V0),
 *                                              constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::focus_element,
 *                                              Functions::ConstantFunction<dim>(
 *                                                -Constants::V0),
 *                                              constraints);
 *     VectorTools::interpolate_boundary_values(dof_handler,
 *                                              BoundaryIds::anode,
 *                                              Functions::ConstantFunction<dim>(
 *                                                +Constants::V0),
 *                                              constraints);
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
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
 * <a name="ThecodeCathodeRaySimulatorassemble_systemcodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::assemble_system</code> function</h4>
 * 

 * 
 * The function that computes
 * the matrix entries is again in essence a copy of the
 * corresponding function in step-6:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::assemble_system()
 *   {
 *     system_matrix = 0;
 *     system_rhs    = 0;
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell = fe.dofs_per_cell;
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0;
 *         cell_rhs    = 0;
 * 
 *         fe_values.reinit(cell);
 * 
 *         for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *           for (const unsigned int i : fe_values.dof_indices())
 *             {
 *               for (const unsigned int j : fe_values.dof_indices())
 *                 cell_matrix(i, j) +=
 *                   (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                    fe_values.JxW(q_index));           // dx
 *             }
 * 
 * @endcode
 * 
 * The only interesting part of this function is how it forms the right
 * hand side of the linear system. Recall that the right hand side
 * of the PDE is
 * @f[
 * \sum_p (N e)\delta(\mathbf x-\mathbf x_p),
 * @f]
 * where we have used $p$ to index the particles here to avoid
 * confusion with the shape function $\varphi_i$; $\mathbf x_p$
 * is the position of the $p$th particle.
 *         

 * 
 * When multiplied by a test function $\varphi_i$ and integrated over
 * the domain results in a right hand side vector
 * @f{align*}{
 * F_i &= \int_\Omega \varphi_i (\mathbf x)\left[
 * \sum_p (N e)\delta(\mathbf x-\mathbf x_p) \right] dx
 * \\  &=  \sum_p (N e) \varphi_i(\mathbf x_p).
 * @f}
 * Note that the final line no longer contains an integral, and
 * consequently also no occurrence of $dx$ which would require the
 * appearance of the `JxW` symbol in our code.
 *         

 * 
 * For a given cell $K$, this cell's contribution to the right hand
 * side is then
 * @f{align*}{
 * F_i^K &= \sum_{p, \mathbf x_p\in K} (N e) \varphi_i(\mathbf x_p),
 * @f}
 * i.e., we only have to worry about those particles that are actually
 * located on the current cell $K$.
 *         

 * 
 * In practice, what we do here is the following: If there are any
 * particles on the current cell, then we first obtain an iterator range
 * pointing to the first particle of that cell as well as the particle
 * past the last one on this cell (or the end iterator) -- i.e., a
 * half-open range as is common for C++ functions. Knowing now the list
 * of particles, we query their reference locations (with respect to
 * the reference cell), evaluate the shape functions in these reference
 * locations, and compute the force according to the formula above
 * (without any FEValues::JxW).
 *         

 * 
 * @note It is worth pointing out that calling the
 * Particles::ParticleHandler::particles_in_cell() and
 * Particles::ParticleHandler::n_particles_in_cell() functions is not
 * very efficient on problems with a large number of particles. But it
 * illustrates the easiest way to write this algorithm, and so we are
 * willing to incur this cost for the moment for expository purposes.
 * We discuss the issue in more detail in the
 * <a href="#extensions">"possibilities for extensions" section</a>
 * below, and use a better approach in step-70, for example.
 * 
 * @code
 *         if (particle_handler.n_particles_in_cell(cell) > 0)
 *           for (const auto &particle : particle_handler.particles_in_cell(cell))
 *             {
 *               const Point<dim> &reference_location =
 *                 particle.get_reference_location();
 *               for (const unsigned int i : fe_values.dof_indices())
 *                 cell_rhs(i) +=
 *                   (fe.shape_value(i, reference_location) * // phi_i(x_p)
 *                    (-Constants::electrons_per_particle *   // N
 *                     Constants::electron_charge));          // e
 *             }
 * 
 * @endcode
 * 
 * Finally, we can copy the contributions of this cell into
 * the global matrix and right hand side vector:
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices);
 *         constraints.distribute_local_to_global(
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorsolve"></a> 
 * <h4>CathodeRaySimulator::solve</h4>
 * 

 * 
 * The function that solves the linear system is then again exactly as in
 * step-6:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::solve_field()
 *   {
 *     SolverControl            solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>> solver(solver_control);
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.2);
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorrefine_grid"></a> 
 * <h4>CathodeRaySimulator::refine_grid</h4>
 * 

 * 
 * The final field-related function is the one that refines the grid. We will
 * call it a number of times in the first time step to obtain a mesh that
 * is well-adapted to the structure of the solution and, in particular,
 * resolves the various singularities in the solution that are due to
 * re-entrant corners and places where the boundary condition type
 * changes. You might want to refer to step-6 again for more details:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::refine_grid()
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
 *                                                     0.1,
 *                                                     0.03);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorcreate_particles"></a> 
 * <h4>CathodeRaySimulator::create_particles</h4>
 * 

 * 
 * Let us now turn to the functions that deal with particles. The first one
 * is about the creation of particles. As mentioned in the introduction,
 * we want to create a particle at points of the cathode if the the electric
 * field $\mathbf E=\nabla V$ exceeds a certain threshold, i.e., if
 * $|\mathbf E| \ge E_\text{threshold}$, and if furthermore the electric field
 * points into the domain (i.e., if $\mathbf E \cdot \mathbf n < 0$). As is
 * common in the finite element method, we evaluate fields (and their
 * derivatives) at specific evaluation points; typically, these are
 * "quadrature points", and so we create a "quadrature formula" that we will
 * use to designate the points at which we want to evaluate the solution.
 * Here, we will simply take QMidpoint implying that we will only check the
 * threshold condition at the midpoints of faces. We then use this to
 * initialize an object of type FEFaceValues to evaluate the solution at these
 * points.
 *   

 * 
 * All of this will then be used in a loop over all cells, their faces, and
 * specifically those faces that are at the boundary and, moreover, the
 * cathode part of the boundary.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::create_particles()
 *   {
 *     FEFaceValues<dim> fe_face_values(fe,
 *                                      QMidpoint<dim - 1>(),
 *                                      update_quadrature_points |
 *                                        update_gradients |
 *                                        update_normal_vectors);
 * 
 *     std::vector<Tensor<1, dim>> solution_gradients(
 *       fe_face_values.n_quadrature_points);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         if (face->at_boundary() &&
 *             (face->boundary_id() == BoundaryIds::cathode))
 *           {
 *             fe_face_values.reinit(cell, face);
 * 
 * @endcode
 * 
 * So we have found a face on the cathode. Next, we let the
 * FEFaceValues object compute the gradient of the solution at each
 * "quadrature" point, and extract the electric field vector from
 * the gradient in the form of a Tensor variable through the methods
 * discussed in the
 * @ref vector_valued "vector-valued problems" documentation module.
 * 
 * @code
 *             const FEValuesExtractors::Scalar electric_potential(0);
 *             fe_face_values[electric_potential].get_function_gradients(
 *               solution, solution_gradients);
 *             for (const unsigned int q_point :
 *                  fe_face_values.quadrature_point_indices())
 *               {
 *                 const Tensor<1, dim> E = solution_gradients[q_point];
 * 
 * @endcode
 * 
 * Electrons can only escape the cathode if the electric field
 * strength exceeds a threshold and,
 * crucially, if the electric field points *into* the domain.
 * Once we have that checked, we create a new
 * Particles::Particle object at this location and insert it
 * into the Particles::ParticleHandler object with a unique ID.
 *                 

 * 
 * The only thing that may be not obvious here is that we also
 * associate with this particle the location in the reference
 * coordinates of the cell we are currently on. This is done
 * because we will in downstream functions compute quantities
 * such as the electric field at the location of the particle
 * (e.g., to compute the forces that act on it when updating its
 * position in each time step). Evaluating a finite element
 * field at arbitrary coordinates is quite an expensive
 * operation because shape functions are really only defined on
 * the reference cell, and so when asking for the electric field
 * at an arbitrary point requires us first to determine what the
 * reference coordinates of that point are. To avoid having to
 * do this over and over, we determine these coordinates once
 * and for all and then store these reference coordinates
 * directly with the particle.
 * 
 * @code
 *                 if ((E * fe_face_values.normal_vector(q_point) < 0) &&
 *                     (E.norm() > Constants::E_threshold))
 *                   {
 *                     const Point<dim> &location =
 *                       fe_face_values.quadrature_point(q_point);
 * 
 *                     Particles::Particle<dim> new_particle;
 *                     new_particle.set_location(location);
 *                     new_particle.set_reference_location(
 *                       mapping.transform_real_to_unit_cell(cell, location));
 *                     new_particle.set_id(next_unused_particle_id);
 *                     particle_handler.insert_particle(new_particle, cell);
 * 
 *                     ++next_unused_particle_id;
 *                   }
 *               }
 *           }
 * 
 * @endcode
 * 
 * At the end of all of these insertions, we let the `particle_handler`
 * update some internal statistics about the particles it stores.
 * 
 * @code
 *     particle_handler.update_cached_numbers();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatormove_particles"></a> 
 * <h4>CathodeRaySimulator::move_particles</h4>
 * 

 * 
 * The second particle-related function is the one that moves the particles
 * in each time step. To do this, we have to loop over all cells, the
 * particles in each cell, and evaluate the electric field at each of the
 * particles' positions.
 *   

 * 
 * The approach used here is conceptually the same used in the
 * `assemble_system()` function: We loop over all cells, find the particles
 * located there (with the same caveat about the inefficiency of the algorithm
 * used here to find these particles), and use FEPointEvaluation object to
 * evaluate the gradient at these positions:
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::move_particles()
 *   {
 *     const double dt = time.get_next_step_size();
 * 
 *     Vector<double>            solution_values(fe.n_dofs_per_cell());
 *     FEPointEvaluation<1, dim> evaluator(mapping, fe, update_gradients);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (particle_handler.n_particles_in_cell(cell) > 0)
 *         {
 *           const typename Particles::ParticleHandler<
 *             dim>::particle_iterator_range particles_in_cell =
 *             particle_handler.particles_in_cell(cell);
 * 
 *           std::vector<Point<dim>> particle_positions;
 *           for (const auto &particle : particles_in_cell)
 *             particle_positions.push_back(particle.get_reference_location());
 * 
 *           cell->get_dof_values(solution, solution_values);
 * 
 * @endcode
 * 
 * Then we can ask the FEPointEvaluation object for the gradients of
 * the solution (i.e., the electric field $\mathbf E$) at these
 * locations and loop over the individual particles:
 * 
 * @code
 *           evaluator.reinit(cell, particle_positions);
 *           evaluator.evaluate(make_array_view(solution_values),
 *                              EvaluationFlags::gradients);
 * 
 *           {
 *             typename Particles::ParticleHandler<dim>::particle_iterator
 *               particle = particles_in_cell.begin();
 *             for (unsigned int particle_index = 0;
 *                  particle != particles_in_cell.end();
 *                  ++particle, ++particle_index)
 *               {
 *                 const Tensor<1, dim> &E =
 *                   evaluator.get_gradient(particle_index);
 * 
 * @endcode
 * 
 * Having now obtained the electric field at the location of one
 * of the particles, we use this to update first the velocity
 * and then the position. To do so, let us first get the old
 * velocity out of the properties stored with the particle,
 * compute the acceleration, update the velocity, and store this
 * new velocity again in the properties of the particle. Recall
 * that this corresponds to the first of the following set of
 * update equations discussed in the introduction:
 * @f{align*}{
 * \frac{{\mathbf v}_i^{(n)}
 * -{\mathbf v}_i^{(n-1)}}{\Delta t}
 * &= \frac{e\nabla V^{(n)}}{m}
 * \\ \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}
 * {\Delta t} &= {\mathbf v}_i^{(n)}.
 * @f}
 * 
 * @code
 *                 const Tensor<1, dim> old_velocity(particle->get_properties());
 * 
 *                 const Tensor<1, dim> acceleration =
 *                   Constants::electron_charge / Constants::electron_mass * E;
 * 
 *                 const Tensor<1, dim> new_velocity =
 *                   old_velocity + acceleration * dt;
 * 
 *                 particle->set_properties(make_array_view(new_velocity));
 * 
 * @endcode
 * 
 * With the new velocity, we can then also update the location
 * of the particle and tell the particle about it.
 * 
 * @code
 *                 const Point<dim> new_location =
 *                   particle->get_location() + dt * new_velocity;
 *                 particle->set_location(new_location);
 *               }
 *           }
 *         }
 * 
 * @endcode
 * 
 * Having updated the locations and properties (i.e., velocities) of all
 * particles, we need to make sure that the `particle_handler` again knows
 * which cells they are in, and what their locations in the coordinate
 * system of the reference cell are. The following function does that. (It
 * also makes sure that, in parallel computations, particles are moved from
 * one processor to another processor if a particle moves from the subdomain
 * owned by the former to the subdomain owned by the latter.)
 * 
 * @code
 *     particle_handler.sort_particles_into_subdomains_and_cells();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatortrack_lost_particle"></a> 
 * <h4>CathodeRaySimulator::track_lost_particle</h4>
 * 

 * 
 * The final particle-related function is the one that is called whenever a
 * particle is lost from the simulation. This typically happens if it leaves
 * the domain. If that happens, this function is called both the cell (which
 * we can ask for its new location) and the cell it was previously on. The
 * function then keeps track of updating the number of particles lost in this
 * time step, the total number of lost particles, and then estimates whether
 * the particle left through the hole in the middle of the anode. We do so by
 * first checking whether the cell it was in last had an $x$ coordinate to the
 * left of the right boundary (located at $x=4$) and the particle now has a
 * position to the right of the right boundary. If that is so, we compute a
 * direction vector of its motion that is normalized so that the $x$ component
 * of the direction vector is equal to $1$. With this direction vector, we can
 * compute where it would have intersected the line $x=4$. If this intersect
 * is between $0.5$ and $1.5$, then we claim that the particle left through
 * the hole and increment a counter.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::track_lost_particle(
 *     const typename Particles::ParticleIterator<dim> &        particle,
 *     const typename Triangulation<dim>::active_cell_iterator &cell)
 *   {
 *     ++n_recently_lost_particles;
 *     ++n_total_lost_particles;
 * 
 *     const Point<dim> current_location              = particle->get_location();
 *     const Point<dim> approximate_previous_location = cell->center();
 * 
 *     if ((approximate_previous_location[0] < 4) && (current_location[0] > 4))
 *       {
 *         const Tensor<1, dim> direction =
 *           (current_location - approximate_previous_location) /
 *           (current_location[0] - approximate_previous_location[0]);
 * 
 *         const double right_boundary_intercept =
 *           approximate_previous_location[1] +
 *           (4 - approximate_previous_location[0]) * direction[1];
 *         if ((right_boundary_intercept > 0.5) &&
 *             (right_boundary_intercept < 1.5))
 *           ++n_particles_lost_through_anode;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorupdate_timestep_size"></a> 
 * <h4>CathodeRaySimulator::update_timestep_size</h4>
 * 

 * 
 * As discussed at length in the introduction, we need to respect a time step
 * condition whereby particles can not move further than one cell in one time
 * step. To ensure that this is the case, we again first compute the maximal
 * speed of all particles on each cell, and divide the cell size by that
 * speed. We then compute the next time step size as the minimum of this
 * quantity over all cells, using the safety factor discussed in the
 * introduction, and set this as the desired time step size using the
 * DiscreteTime::set_desired_time_step_size() function.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::update_timestep_size()
 *   {
 *     if (time.get_step_number() > 0)
 *       {
 *         double min_cell_size_over_velocity = std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : dof_handler.active_cell_iterators())
 *           if (particle_handler.n_particles_in_cell(cell) > 0)
 *             {
 *               const double cell_size = cell->minimum_vertex_distance();
 * 
 *               double max_particle_velocity(0.0);
 * 
 *               for (const auto &particle :
 *                    particle_handler.particles_in_cell(cell))
 *                 {
 *                   const Tensor<1, dim> velocity(particle.get_properties());
 *                   max_particle_velocity =
 *                     std::max(max_particle_velocity, velocity.norm());
 *                 }
 * 
 *               if (max_particle_velocity > 0)
 *                 min_cell_size_over_velocity =
 *                   std::min(min_cell_size_over_velocity,
 *                            cell_size / max_particle_velocity);
 *             }
 * 
 *         constexpr double c_safety = 0.5;
 *         time.set_desired_next_step_size(c_safety * 0.5 *
 *                                         min_cell_size_over_velocity);
 *       }
 * @endcode
 * 
 * As mentioned in the introduction, we have to treat the very first
 * time step differently since there, particles are not available yet or
 * do not yet have the information associated that we need for the
 * computation of a reasonable step length. The formulas below follow the
 * discussion in the introduction.
 * 
 * @code
 *     else
 *       {
 *         const QTrapezoid<dim> vertex_quadrature;
 *         FEValues<dim> fe_values(fe, vertex_quadrature, update_gradients);
 * 
 *         std::vector<Tensor<1, dim>> field_gradients(vertex_quadrature.size());
 * 
 *         double min_timestep = std::numeric_limits<double>::max();
 * 
 *         for (const auto &cell : dof_handler.active_cell_iterators())
 *           if (particle_handler.n_particles_in_cell(cell) > 0)
 *             {
 *               const double cell_size = cell->minimum_vertex_distance();
 * 
 *               fe_values.reinit(cell);
 *               fe_values.get_function_gradients(solution, field_gradients);
 * 
 *               double max_E = 0;
 *               for (const auto q_point : fe_values.quadrature_point_indices())
 *                 max_E = std::max(max_E, field_gradients[q_point].norm());
 * 
 *               if (max_E > 0)
 *                 min_timestep =
 *                   std::min(min_timestep,
 *                            std::sqrt(0.5 * cell_size *
 *                                      Constants::electron_mass /
 *                                      Constants::electron_charge / max_E));
 *             }
 * 
 *         time.set_desired_next_step_size(min_timestep);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeCathodeRaySimulatoroutput_resultscodefunction"></a> 
 * <h4>The <code>CathodeRaySimulator::output_results()</code> function</h4>
 * 

 * 
 * The final function implementing pieces of the overall algorithm is the one
 * that generates graphical output. In the current context, we want to output
 * both the electric potential field as well as the particle locations and
 * velocities. But we also want to output the electric field, i.e., the
 * gradient of the solution.
 *   

 * 
 * deal.II has a general way how one can compute derived quantities from
 * the solution and output those as well. Here, this is the electric
 * field, but it could also be some other quantity -- say, the norm of the
 * electric field, or in fact anything else one could want to compute from
 * the solution $V_h(\mathbf x)$ or its derivatives. This general solution
 * uses the DataPostprocessor class and, in cases like the one here where we
 * want to output a quantity that represents a vector field, the
 * DataPostprocessorVector class.
 *   

 * 
 * Rather than try and explain how this class works, let us simply refer to
 * the documentation of the DataPostprocessorVector class that has essentially
 * this case as a well-documented example.
 * 
 * @code
 *   template <int dim>
 *   class ElectricFieldPostprocessor : public DataPostprocessorVector<dim>
 *   {
 *   public:
 *     ElectricFieldPostprocessor()
 *       : DataPostprocessorVector<dim>("electric_field", update_gradients)
 *     {}
 * 
 *     virtual void evaluate_scalar_field(
 *       const DataPostprocessorInputs::Scalar<dim> &input_data,
 *       std::vector<Vector<double>> &computed_quantities) const override
 *     {
 *       AssertDimension(input_data.solution_gradients.size(),
 *                       computed_quantities.size());
 * 
 *       for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p)
 *         {
 *           AssertDimension(computed_quantities[p].size(), dim);
 *           for (unsigned int d = 0; d < dim; ++d)
 *             computed_quantities[p][d] = input_data.solution_gradients[p][d];
 *         }
 *     }
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * With this, the `output_results()` function becomes relatively
 * straightforward: We use the DataOut class as we have in almost every one of
 * the previous tutorial programs to output the solution (the "electric
 * potential") and we use the postprocessor defined above to also output its
 * gradient (the "electric field"). This all is then written into a file in
 * VTU format after also associating the current time and time step number
 * with this file.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::output_results() const
 *   {
 *     {
 *       ElectricFieldPostprocessor<dim> electric_field;
 *       DataOut<dim>                    data_out;
 *       data_out.attach_dof_handler(dof_handler);
 *       data_out.add_data_vector(solution, "electric_potential");
 *       data_out.add_data_vector(solution, electric_field);
 *       data_out.build_patches();
 * 
 *       data_out.set_flags(
 *         DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));
 * 
 *       std::ofstream output("solution-" +
 *                            Utilities::int_to_string(time.get_step_number(), 4) +
 *                            ".vtu");
 *       data_out.write_vtu(output);
 *     }
 * 
 * @endcode
 * 
 * Output the particle positions and properties is not more complicated. The
 * Particles::DataOut class plays the role of the DataOut class for
 * particles, and all we have to do is tell that class where to take
 * particles from and how to interpret the `dim` components of the
 * properties -- namely, as a single vector indicating the velocity, rather
 * than as `dim` scalar properties. The rest is then the same as above:
 * 
 * @code
 *     {
 *       Particles::DataOut<dim, dim> particle_out;
 *       particle_out.build_patches(
 *         particle_handler,
 *         std::vector<std::string>(dim, "velocity"),
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>(
 *           dim, DataComponentInterpretation::component_is_part_of_vector));
 * 
 *       particle_out.set_flags(
 *         DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number()));
 * 
 *       std::ofstream output("particles-" +
 *                            Utilities::int_to_string(time.get_step_number(), 4) +
 *                            ".vtu");
 *       particle_out.write_vtu(output);
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="CathodeRaySimulatorrun"></a> 
 * <h4>CathodeRaySimulator::run</h4>
 * 

 * 
 * The last member function of the principal class of this program is then the
 * driver. At the top, it refines the mesh a number of times by solving the
 * problem (with not particles yet created) on a sequence of finer and finer
 * meshes.
 * 
 * @code
 *   template <int dim>
 *   void CathodeRaySimulator<dim>::run()
 *   {
 *     make_grid();
 * 
 * @endcode
 * 
 * do a few refinement cycles up front
 * 
 * @code
 *     const unsigned int n_pre_refinement_cycles = 3;
 *     for (unsigned int refinement_cycle = 0;
 *          refinement_cycle < n_pre_refinement_cycles;
 *          ++refinement_cycle)
 *       {
 *         setup_system();
 *         assemble_system();
 *         solve_field();
 *         refine_grid();
 *       }
 * 
 * 
 * @endcode
 * 
 * Now do the loop over time. The sequence of steps follows closely the
 * outline of the algorithm discussed in the introduction. As discussed in
 * great detail in the documentation of the DiscreteTime class, while we
 * move the field and particle information forward by one time step, the
 * time stored in the `time` variable is not consistent with where (some of)
 * these quantities are (in the diction of DiscreteTime, this is the "update
 * stage"). The call to `time.advance_time()` makes everything consistent
 * again by setting the `time` variable to the time at which the field and
 * particles already are, and once we are in this "consistent stage", we can
 * generate graphical output and write information about the current state
 * of the simulation to screen.
 * 
 * @code
 *     setup_system();
 *     do
 *       {
 *         std::cout << "Timestep " << time.get_step_number() + 1 << std::endl;
 *         std::cout << "  Field degrees of freedom:                 "
 *                   << dof_handler.n_dofs() << std::endl;
 * 
 *         assemble_system();
 *         solve_field();
 * 
 *         create_particles();
 *         std::cout << "  Total number of particles in simulation:  "
 *                   << particle_handler.n_global_particles() << std::endl;
 * 
 *         n_recently_lost_particles = 0;
 *         update_timestep_size();
 *         move_particles();
 * 
 *         time.advance_time();
 * 
 *         output_results();
 * 
 *         std::cout << "  Number of particles lost this time step:  "
 *                   << n_recently_lost_particles << std::endl;
 *         if (n_total_lost_particles > 0)
 *           std::cout << "  Fraction of particles lost through anode: "
 *                     << 1. * n_particles_lost_through_anode /
 *                          n_total_lost_particles
 *                     << std::endl;
 * 
 *         std::cout << std::endl
 *                   << "  Now at t=" << time.get_current_time()
 *                   << ", dt=" << time.get_previous_step_size() << '.'
 *                   << std::endl
 *                   << std::endl;
 *       }
 *     while (time.is_at_end() == false);
 *   }
 * } // namespace Step19
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
 * The final function of the program is then again the `main()` function. It is
 * unchanged in all tutorial programs since step-6 and so there is nothing new
 * to discuss:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       Step19::CathodeRaySimulator<2> cathode_ray_simulator_2d;
 *       cathode_ray_simulator_2d.run();
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
 *   return 0;
 * }
 * @endcode
examples/step-19/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当这个程序运行时，它产生的输出看起来如下。``时间步数1 场自由度：4989 仿真中的粒子总数：20 这个时间步数损失的粒子数：0

  现在在t=2.12647e-07，dt=2.12647e-07。

时间步数2 场自由度：4989 仿真中的粒子总数：24 本时间步数损失的粒子数：0

  现在在t=4.14362e-07，dt=2.01715e-07。

时间步数3 场自由度：4989 仿真中的粒子总数：28 本时间步数损失的粒子数：0

  现在在t=5.96019e-07，dt=1.81657e-07。

时间步数4 场自由度：4989 仿真中的粒子总数。  32 这个时间步长损失的粒子数：0

  现在在t=7.42634e-07，dt=1.46614e-07。


...


  时间步数1000场自由度：4989模拟中的粒子总数。  44 这个时间步长损失的粒子数：6 通过阳极损失的粒子的比例。0.0601266

  现在在t=4.93276e-05，dt=4.87463e-08。

时间步数1001场自由度：4989模拟中的粒子总数。  44 这个时间步长损失的粒子数：0 通过阳极损失的粒子的分数。0.0601266

  现在在t=4.93759e-05，dt=4.82873e-08。


...


时间步数2091场自由度：4989模拟中的粒子总数。  44 这个时间步长损失的粒子数：0 通过阳极损失的粒子的比例。0.0503338

  现在在t=9.99237e-05，dt=4.26254e-08。

时间步数2092场自由度：4989模拟中的粒子总数。  44 这个时间步长损失的粒子数：0 通过阳极损失的粒子的分数。0.0503338

  现在在t=9.99661e-05，dt=4.24442e-08。

时间步数2093场自由度：4989模拟中的粒子总数。  44 这个时间步长损失的粒子数：2 通过阳极损失的粒子的比例。0.050308

  现在在t=0.0001，dt=3.38577e-08。```

随机选取几个时间步长，我们可以用电场的流线和电子的点的形式来可视化解决方案。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-19.solution.0000.png" alt="时间步骤0（t=0秒）的解决方案。" width="500">  <br>  时间步骤0（t=0秒）的解决方案。       <br>  </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-19.solution.1400.png" alt="时间步骤1400（t=0.000068秒）的解决方案。" width="500">  <br>  时间步骤1400（t=0.000068秒）的解决方案。       <br>  </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-19.solution.0700.png" alt="时间步骤700（t=0.000035秒）的解决方案。" width="500">  <br>  在时间步骤700（t=0.000035秒）的解决方案。       <br>  </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-19.solution.2092.png" alt="时间步骤2092（t=0.0001秒）的解决方案。" width="500">  <br>  时间步骤2092（t=0.0001秒）的解决方案。       <br>  </div> </div>

也就是说，更合适的方式是通过创建一个视频，展示这些电子是如何运动的，以及电场是如何随着它们的运动而变化的，从而将这个程序的结果可视化。

@htmlonly
<p align="center">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/HwUtE7xuteE"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly



在这里你可以看到边界的 "焦点元素 "是如何用其负电压排斥电子，并确保它们不会垂直于阴极飞走（就像它们在其轨迹的初始部分那样）。它还显示了电场线如何随着时间的推移而移动，以回应飞过的电荷--换句话说，粒子对电场的反馈，而电场本身驱动着电子的运动。

这部电影表明，电子是以 "成串 "或 "爆裂 "的方式移动的。这种表象的一个因素是电影是如何创建的，是一个伪影。电影的每一帧都对应着一个时间步长，但时间步长是不同的。更具体地说，穿过最小单元的最快粒子决定了时间步长（见介绍中的讨论），因此，每当一个（快速）粒子穿过域的右边缘的小单元时，时间步长都很小；一旦粒子离开域，时间步长又会变长。通过绘制屏幕输出中显示的时间步长，可以很容易地看到这种减速-加速的效果。

然而，这其中的第二部分是真实的。模拟在开始时创造了一大群粒子，而在大约第300个时间步长后，粒子数量就减少了。这可能是因为模拟中的粒子带有负电荷。它们降低了（同样带负电的电极）的电场强度，因此减少了阴极上的点的数量，在这些点上，电场强度超过了将电子从电极中吸引出来所需的阈值。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Avoidingaperformancebottleneckwithparticles"></a><h4> Avoiding a performance bottleneck with particles </h4>


`assemble_system()`、`move_particles()`和`update_timestep_size()`函数都调用 Particles::ParticleHandler::particles_in_cell() 和 Particles::ParticleHandler::n_particles_in_cell() ，查询位于当前单元上的粒子信息。虽然这很方便，但也很低效。为了理解为什么会这样，我们需要知道粒子是如何存储在 Particles::ParticleHandler: 中的，即在一个数据结构中，粒子是以某种线性方式按它们所在的单元进行排序的。因此，为了找到与给定单元相关的粒子，这些函数需要搜索给定单元上的第一个（也可能是最后一个）粒子--这种努力需要花费 ${\cal O}(\log N)$ 次操作，其中 $N$ 是粒子的数量。但这是在每个单元上重复的；假设对于大型计算来说，单元和粒子的数量大致成正比，那么这些函数调用的累积成本是 ${\cal O}(N \log N)$ ，因此大于我们应该对程序的所有部分进行的 ${\cal O}(N)$ 成本。

不过，我们可以使之更便宜。首先，我们可以先调用 Particles::ParticleHandler::n_particles_in_cell(), 而不是 Particles::ParticleHandler::particles_in_cell() ，然后通过计算当前单元上的最后一个粒子到第一个粒子的距离来计算单元上的粒子数。

@code
  const typename Particles::ParticleHandler<dim, spacedim>::particle_iterator_range
    particles_in_cell = particle_handler.particles_in_cell(cell);
  const unsigned int
    n_particles_in_cell = std::distance (particles_in_cell.begin(),
                                         particles_in_cell.end());
@endcode

其中第一个调用当然还是 ${\cal O}(\log N)$ ，但至少第二个调用只需要与当前单元上的粒子数成比例的计算时间，因此，当累积到所有单元时，其成本为 ${\cal O}(N)$ 。

但我们甚至可以通过一些适当的算法设计来摆脱这些调用中的第一个。这是因为粒子的排列方式与单元格相同，因此我们可以在单元格上移动时直接走动它们。下面的算法纲要就是这样做的。

@code
  auto begin_particle_on_cell = particle_handler.begin();
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      unsigned int n_particles_on_cell = 0;
      auto end_particle_on_cell = begin_particle_on_cell;
      while (end_particle_on_cell->get_surrounding_cell(triangulation)
             == cell)
        {
          ++n_particles_on_cell;
          ++end_particle_on_cell;
        }


      ...now operate on the range of particles from begin_particle_on_cell
         to end_particle_on_cell, all of which are known to be on the current
         cell...;


      // Move the begin iterator forward so that it points to the first
      // particle on the next cell
      begin_particle_on_cell = end_particle_on_cell;
    }
@endcode



在这段代码中，我们对每个单元都精确地接触了一次，而且我们从来不需要在大数据结构中搜索每个单元上的第一个或最后一个粒子。因此，该算法总共花费了 ${\cal O}(N)$ 来完成对所有粒子和所有单元的扫瞄。

对这个程序中存在这个问题的所有三个函数实施这个方案并不十分困难。




<a name="Morestatisticsaboutelectrons"></a><h4> More statistics about electrons </h4>


该程序已经计算出了通过阳极上的孔离开该领域的电子的比例。但人们可能还对其他数量感兴趣。例如，这些粒子的平均速度。从每个粒子的属性中获得其速度并不是很困难，就像我们在`move_particles()`函数中所做的那样，并从中计算出统计数据。




<a name="Abettersynchronizedvisualization"></a><h4> A better-synchronized visualization </h4>


如上所述，视频的不同帧之间有不同的时间差，因为我们为每个时间步长创建输出。一个更好的创建电影的方法是在固定的时间间隔内生成一个新的输出文件，不管每个这样的点之间有多少时间步长。




<a name="Abettertimeintegrator"></a><h4> A better time integrator </h4>


我们在这个程序中考虑的问题是一个耦合的、多物理学的问题。但是我们解决它的方法是首先计算（电）势场，然后更新粒子位置。这就是所谓的 "算子分割法"，我们将在第58步中更详细地研究这一概念。

虽然要想出一个不涉及将问题分割成PDE部分和粒子部分的方法是很尴尬的，但人们*可以*（而且可能应该！）想出一个更好的方法来更新粒子的位置。具体来说，我们用来更新粒子位置的方程是

@f{align*}{
  \frac{{\mathbf v}_i^{(n)}-{\mathbf v}_i^{(n-1)}}{\Delta t} &= \frac{e\nabla V^{(n)}}{m}
  \\
  \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}{\Delta t} &= {\mathbf v}_i^{(n)}.


@f}

这相当于一个简单的正向欧拉时间离散化--一种在时间步长上具有一阶精度的方法 $\Delta t$ ，我们知道我们应该避免，因为我们可以做得更好。相反，我们可能想考虑一种方案，如[跃迁方案](https://en.wikipedia.org/wiki/Leapfrog_integration)或更普遍的[折衷积分器](https://en.wikipedia.org/wiki/Symplectic_integrator)，如[Verlet方案](https://en.wikipedia.org/wiki/Verlet_integration)。




<a name="Parallelization"></a><h4> Parallelization </h4>


在写这篇文章时，在作者的一台笔记本电脑上，在发布模式下，该程序的运行时间约为3.5分钟。这是可以接受的。但是，如果我们想让模拟变成三维的呢？如果我们想在任何时候都不使用最多约100个粒子（如这里使用的参数），而是使用100,000个？如果我们需要一个更细的网格？

在这些情况下，人们不只是想在单个处理器上运行程序，实际上是在尽可能多的处理器上运行。这就要求对PDE的解决方案以及粒子进行并行化。在实践中，虽然有大量的挑战要使其高效和良好地扩展，但这些挑战都在deal.II本身中得到了解决。例如，第40步显示了如何将有限元部分并行化，第70步显示了如何将粒子部分也并行化。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-19.cc"
*/
