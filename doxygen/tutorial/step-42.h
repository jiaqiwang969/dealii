/**
@page step_42 The step-42 tutorial program
This tutorial depends on step-41, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Introduction">Introduction</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Reformulationasavariationalinequality">Reformulation as a variational inequality</a>
        <li><a href="#ANewtonmethodfortheplasticnonlinearity">A Newton method for the plastic nonlinearity</a>
        <li><a href="#ActiveSetmethodstosolvethesaddlepointproblem">Active Set methods to solve the saddle point problem</a>
        <li><a href="#Overallalgorithm">Overall algorithm</a>
        <li><a href="#Adaptivemeshrefinement">Adaptive mesh refinement</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeConstitutiveLawcodeclasstemplate">The <code>ConstitutiveLaw</code> class template</a>
      <ul>
        <li><a href="#ConstitutiveLawget_stress_strain_tensor">ConstitutiveLaw::get_stress_strain_tensor</a>
        <li><a href="#ConstitutiveLawget_linearized_stress_strain_tensors">ConstitutiveLaw::get_linearized_stress_strain_tensors</a>
        <li><a href="#ThecodeSphereObstaclecodeclass">The <code>SphereObstacle</code> class</a>
        <li><a href="#ThecodeBitmapFilecodeandcodeChineseObstaclecodeclasses">The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes</a>
      </ul>
        <li><a href="#ThecodePlasticityContactProblemcodeclasstemplate">The <code>PlasticityContactProblem</code> class template</a>
        <li><a href="#ImplementationofthecodePlasticityContactProblemcodeclass">Implementation of the <code>PlasticityContactProblem</code> class</a>
      <ul>
        <li><a href="#PlasticityContactProblemdeclare_parameters">PlasticityContactProblem::declare_parameters</a>
        <li><a href="#ThecodePlasticityContactProblemcodeconstructor">The <code>PlasticityContactProblem</code> constructor</a>
        <li><a href="#PlasticityContactProblemmake_grid">PlasticityContactProblem::make_grid</a>
        <li><a href="#PlasticityContactProblemsetup_system">PlasticityContactProblem::setup_system</a>
        <li><a href="#PlasticityContactProblemcompute_dirichlet_constraints">PlasticityContactProblem::compute_dirichlet_constraints</a>
        <li><a href="#PlasticityContactProblemassemble_mass_matrix_diagonal">PlasticityContactProblem::assemble_mass_matrix_diagonal</a>
        <li><a href="#PlasticityContactProblemupdate_solution_and_constraints">PlasticityContactProblem::update_solution_and_constraints</a>
        <li><a href="#PlasticityContactProblemassemble_newton_system">PlasticityContactProblem::assemble_newton_system</a>
        <li><a href="#PlasticityContactProblemcompute_nonlinear_residual">PlasticityContactProblem::compute_nonlinear_residual</a>
        <li><a href="#PlasticityContactProblemsolve_newton_system">PlasticityContactProblem::solve_newton_system</a>
        <li><a href="#PlasticityContactProblemsolve_newton">PlasticityContactProblem::solve_newton</a>
        <li><a href="#PlasticityContactProblemrefine_grid">PlasticityContactProblem::refine_grid</a>
        <li><a href="#PlasticityContactProblemmove_mesh">PlasticityContactProblem::move_mesh</a>
        <li><a href="#PlasticityContactProblemoutput_results">PlasticityContactProblem::output_results</a>
        <li><a href="#PlasticityContactProblemoutput_contact_force">PlasticityContactProblem::output_contact_force</a>
        <li><a href="#PlasticityContactProblemrun">PlasticityContactProblem::run</a>
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
examples/step-42/doc/intro.dox

 <br> 

<i>This program was contributed by Jörg Frohne (University of Siegen,
Germany) while on a long-term visit to Texas A&amp;M University, with significant
contributions by Timo Heister and Wolfgang Bangerth.
<br>
<br>
The code described here provides the basis for the numerical experiments shown
in the following paper:
<br>
  J. Frohne, T. Heister, W. Bangerth: <b>Efficient numerical methods for the large-scale, parallel
                  solution of elastoplastic contact problems</b><b>Efficient numerical methods for the large-scale, parallel
                  solution of elastoplastic contact problems</b>.
  Accepted for publication in International Journal for Numerical Methods in Engineering, 2015.
</i> 。




<a name="Intro"></a>

<a name="Introduction"></a><h3>Introduction</h3>


这个例子是第41步的延伸，考虑的是三维接触问题，具有各向同性硬化的弹塑性材料行为。换句话说，它考虑的是，如果把一个刚性的障碍物推到一个三维体上，它是如何变形的（接触问题），其中的变形受弹塑性材料法则（一种只能容纳一定最大应力的材料）的制约，随着变形的累积，该材料会硬化。为了说明我们打算做什么，在讨论太多细节之前，让我们只展示一张解决方案的图片（可变形体是一个立方体--实际上只显示了一半--，障碍物对应于一个汉字，将在下面讨论）。

 <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionLi2.png" alt=""> 


这个问题的描述意味着，与第41步相比，我们必须照顾到一个额外的非线性因素：材料行为。由于我们在这里考虑的是一个三维问题，我们还必须考虑到一个事实，即现在接触区是在可变形体的边界，而不是在内部。最后，与第41步相比，我们还必须在处理线性系统和不等式约束时处理悬空节点，因为我们希望使用自适应网格；在后一种情况下，我们将不得不处理优先考虑悬空节点的约束还是不等式的约束更重要。

由于在三维空间中很容易达到几百万个自由度，即使使用自适应网格细化，我们决定使用Trilinos和p4est来并行运行我们的代码，在步骤40的框架上进行并行化。并行化的其他指针可以在步骤32中找到。




<a name="Classicalformulation"></a><h3>Classical formulation</h3>


该问题的经典表述具有以下形式。

@f{align*}
 \varepsilon(\mathbf u) &= A\sigma + \varepsilon^p & &\quad\text{in } \Omega,\\


  -\textrm{div}\ \sigma &= \mathbf f & &\quad\text{in } \Omega,\\
  \varepsilon^p:(\tau - \sigma) &\geq 0\quad\forall\tau\text{ with
  }\mathcal{F}(\tau)\leq 0 & &\quad\text{in } \Omega,\\
  \mathbf u &= 0 & &\quad\text{on }\Gamma_D,\\
  \sigma \mathbf n - [\mathbf n \cdot(\sigma \mathbf n)]\mathbf n &= 0,
  \quad \mathbf n \cdot (\sigma
  \mathbf n) \leq 0 & &\quad\text{on }\Gamma_C,\\
  (\mathbf n \cdot (\sigma
  \mathbf n))(\mathbf n \cdot \mathbf u - g) &= 0,\quad \mathbf n
  \cdot \mathbf u - g \leq 0 & &\quad\text{on } \Gamma_C.


@f}

这里，这些方程的第一个定义了应变 $\varepsilon(\mathbf u)=\frac{1}{2}\left(\nabla \mathbf u
  + \nabla \mathbf u^T\right)$ 和应力 $\sigma$ 之间的关系，通过四阶顺应性张量 $A$ ； $\varepsilon^p$ 提供了应变的塑性成分，确保应力不超过屈服应力。我们将只考虑各向同性的材料，对于这些材料， $A$ 可以用Lam&eacute;模量 $\lambda$ 和 $\mu$ 表示，或者用体模量 $\kappa$ 和 $\mu$ 表示。第二个方程是力的平衡；我们在此不考虑任何体力，并假定 $\mathbf f=0$  。第三行的互补条件意味着，如果 $\mathcal{F}(\sigma)< 0$ ，则 $\varepsilon^p=0$ ，但当且仅当 $\mathcal{F}(\sigma) = 0$ ， $\varepsilon^p$ 可能是一个非零张量，特别是在这种情况下， $\varepsilon^p$ 必须指向 $\partial
\mathcal{F}(\sigma)/\partial \sigma$ 的方向。不等式 $\mathcal{F}(\sigma)\le 0$ 是塑性材料只能支持有限的应力；换句话说，如果外力会导致 $\sigma$ 的应力，那么它们就会产生塑性变形 $\varepsilon^p$ 的反应。这种<i>yield function</i>的典型形式是 $\mathcal{F}(\sigma)=|\sigma^D|-\sigma_{\text{yield}}$ ，其中 $\tau^D
= \tau - \dfrac{1}{3}tr(\tau)I$ 是张量的偏离部分， $|\cdot|$ 表示弗罗本尼斯规范。

进一步的方程描述了 $\Gamma_D$ 上固定的零位移，在可能出现接触的表面 $\Gamma_C=\partial\Omega\backslash\Gamma_D$ 上，障碍物施加的法向力 $\sigma_n=\mathbf n \cdot (\sigma(\mathbf u)
  \mathbf n)$ 是向内的（障碍物对我们的身体没有 "拉力"），切向分量为零 $\mathbf \sigma_t= \sigma \mathbf n - \mathbf \sigma_n \mathbf n
= \sigma \mathbf n - [\mathbf n \cdot(\sigma \mathbf n)]\mathbf n$  。最后一个条件又是一个互补条件，意味着在 $\Gamma_C$ 上，只有当身体与障碍物接触时，法向力才能非零；第二部分描述了障碍物和身体的不可穿透性。最后两个方程通常被称为Signorini接触条件。

大多数材料--尤其是金属--都有这样的特性，即它们在变形时表现出一定的硬化。换句话说， $\sigma_{\text{yield}}$ 随着变形而增加。在实践中，导致硬化的不是弹性变形，而是塑性成分。有不同的构成法则来描述这些材料行为。最简单的称为线性各向同性硬化，由流动函数  $\mathcal{F}(\sigma,\varepsilon^p) = \vert\sigma^D\vert - (\sigma_0 +
\gamma^{\text{iso}}|\varepsilon^p|)$  描述。




<a name="Reformulationasavariationalinequality"></a><h3>Reformulation as a variational inequality</h3>


一般来说，处理不等式是相当笨拙的。在这里，我们必须处理两个问题：塑性和接触问题。正如本页顶部提到的论文中详细描述的那样，我们至少可以重新表述塑性，使其看起来像一个非线性，然后我们可以用牛顿方法处理。这在数学上略显棘手，因为非线性不只是一些平滑的函数，而是在应力达到屈服应力的地方有结点；然而，对于这样的<i>semismooth</i>函数，可以证明牛顿方法仍然收敛。

在不涉及细节的情况下，我们也将摆脱作为独立变量的应力，而完全用位移来工作  $\mathbf u$  。最终，这种重构的目标是，我们希望最终得到一个对称的、正定的问题--比如一个线性化的弹性问题，其空间变量系数由塑性行为产生--需要在每个牛顿步骤中解决。我们希望如此，因为有高效和可扩展的方法来解决这样的线性系统，如用代数多重网格的CG预处理。这与我们继续使用包含位移和应力的混合公式所得到的类似于混合拉普拉斯的鞍点问题（见第20步）是相反的，第20步已经提示了构建良好的求解器和预处理器是多么困难。

说到这里，让我们简单陈述一下我们在重构后得到的问题（同样，细节可以在论文中找到）。找到一个位移 $\mathbf u \in
V^+$ ，以便

@f{align*}
\left(P_{\Pi}(C\varepsilon(\mathbf u)),\varepsilon(\varphi) - \varepsilon(\mathbf u)\right) \geq 0,\quad \forall \varphi\in V^+.


@f}

其中投影仪 $P_\Pi$ 被定义为

@f{align*}
 P_{\Pi}(\tau) \dealcoloneq \begin{cases}
    \tau, & \text{if }\vert\tau^D\vert \leq \sigma_0,\\
    \left[
      \dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}} +
      \left(1-\dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}\right)\dfrac{\sigma_0}{\vert\tau^D\vert}
    \right]\tau^D
    + \dfrac{1}{3}\text{trace}(\tau) I, & \text{if }\vert\tau^D\vert >
    \sigma_0,
  \end{cases}


@f}

和空间 $V^+$ 是满足接触条件的所有位移的空间。

@f{align*}
  V
  &=
  \left\{ \mathbf u\in \left[H^1(\Omega)\right]^{d}:
    \mathbf u = 0 \text{ on } \Gamma_D\right\},
  \\
  V^+
  &=
  \left\{ \mathbf u\in V: \mathbf n \cdot \mathbf u\leq g \text{ on } \Gamma_C \right\}.


@f}



在实际代码中，我们将使用缩写  $\gamma=\dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}$  。

鉴于这种表述，我们将应用两种技术。

- 运行牛顿方法来迭代出投影仪的非线性。

- 为接触条件运行一个主动设置方法，方法与我们在步骤41中所做的基本相同。

一个严格的方法是在我们迭代牛顿方法到收敛时保持活动集的固定（或者也许反过来：在进入下一个牛顿迭代之前找到最终的活动集）。在实践中，事实证明，每个活动集迭代只做一个牛顿步骤就足够了，所以我们将同时迭代它们。我们还将每隔一段时间细化一下网格。




<a name="ANewtonmethodfortheplasticnonlinearity"></a><h3>A Newton method for the plastic nonlinearity</h3>


如前所述，我们将通过应用牛顿方法来处理算子 $P_\Pi$ 的非线性，尽管该算子在严格意义上是不可微的。然而，它满足了<i>slant</i>的可微条件，这就足以使牛顿方法发挥作用。由此产生的方法被称为<i>semi-smooth Newton method</i>，听起来令人印象深刻，但实际上只是一个牛顿方法应用于一个具有适当选择的 "导数 "的半光滑函数。

在目前的情况下，我们将通过在每个迭代 $i$ 中求解以下方程来运行我们的迭代（仍然是不等式，但是线性化）。

@f{align*}
  \label{eq:linearization}
  \left(I_{\Pi}\varepsilon(\tilde {\mathbf u}^{i}),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right) \geq
  \left(\left(I_{\Pi}\varepsilon({\mathbf u}^{i-1}),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right) -
  \left(P_{\Pi}(C\varepsilon({\mathbf u}^{i-1})),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right)\right),
  \quad \forall \varphi\in V^+,


@f}

其中，等级4张量 $I_\Pi=I_\Pi(\varepsilon^D(\mathbf u^{i-1}))$ 由以下公式给出

@f{align}
  I_\Pi = \begin{cases}
    C_{\mu} + C_{\kappa}, & \hspace{-8em} \text{if } \vert C\varepsilon^D(\mathbf u^{i-1}) \vert \leq \sigma_0,
    \\
    \frac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}} C_{\mu} + \frac{\left(1-\frac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}\right)\sigma_0}{\vert C\varepsilon^D(\mathbf u^{i-1}) \vert}\left(C_{\mu} -
      2\mu\dfrac{C\varepsilon^D(\mathbf u^{i-1})\otimes C\varepsilon^D(\mathbf
        u^{i-1})}{\vert C\varepsilon^D(\mathbf u^{i-1})\vert^2}\right) + C_{\kappa}, & \text{ else.}
\end{cases}


@f}

这个张量是 $P_\Pi(C\cdot)$ 围绕 $\varepsilon^D(\mathbf u^{i-1})$ 的（形式）线性化。对于我们这里考虑的线性各向同性材料，投影仪的体积和剪切分量由以下公式给出

@f{gather*}
  C_{\kappa} = \kappa I\otimes I,
  \qquad\qquad\qquad\qquad
  C_{\mu} = 2\mu\left(\mathbb{I}  - \dfrac{1}{3} I\otimes
    I\right),


@f}

其中 $I$ 和 $\mathbb{I}$ 分别是等级为2和4的认同张量。

请注意，这个问题对应于线性弹性接触问题，其中 $I_\Pi$ 扮演弹性张量的角色  $C=A^{-1}$  。事实上，如果材料在某一点上没有塑性，那么 $I_\Pi=C$  。然而，在材料具有塑性的地方， $I_\Pi$ 是一个空间变化的函数。在任何情况下，我们必须解决牛顿迭代的系统 $\tilde {\mathbf u}^{i}$ 使我们更接近重写我们问题的目标，使我们能够使用众所周知的椭圆系统的求解器和预处理器。

作为对牛顿方法的最后说明，让我们提一下，正如牛顿方法常见的那样，我们需要通过控制步长来使其全球化。换句话说，虽然上面的系统求解的是 $\tilde {\mathbf u}^{i}$ ，但最后的迭代结果将是

@f{align*}
  {\mathbf u}^{i} = {\mathbf u}^{i-1} + \alpha_i (\tilde {\mathbf u}^{i} - {\mathbf u}^{i-1})


@f}

其中右边括号中的差值扮演了传统牛顿方向的角色， $\delta {\mathbf u}^{i}$  。我们将用标准的直线搜索来确定 $\alpha^i$ 。




<a name="ActiveSetmethodstosolvethesaddlepointproblem"></a><h3>Active Set methods to solve the saddle point problem</h3>


这个要在每个牛顿步骤中解决的线性化问题基本上与步骤41一样。唯一的区别在于接触区是在边界而不是在域中。但这没有进一步的后果，所以我们参考步骤41的文件，唯一的提示是 $\mathcal{S}$ 这次包含了接触边界的所有顶点 $\Gamma_C$ 。和那里一样，我们需要做的是保持一个自由度子集的固定，导致额外的约束，可以写成一个鞍点问题。然而，正如论文中所讨论的，通过以适当的方式写这些约束，消除自由度之间的耦合，我们最终会得到一组节点，这些节点基本上只是附加了Dirichlet值。




<a name="Overallalgorithm"></a><h3>Overall algorithm</h3>


上述算法结合了阻尼半光滑牛顿法（我们用于非线性构成法）和半光滑牛顿法用于接触。它的工作原理如下。<ol>  <li>  初始化活动和非活动集 $\mathcal{A}_i$ 和 $\mathcal{F}_i$ ，使 $\mathcal{S} = \mathcal{A}_i \cup \mathcal{F}_i$ 和 $\mathcal{A}_i \cap
 \mathcal{F}_i = \emptyset$ 和集 $i = 1$  。这里， $\mathcal{S}$ 是位于可能发生接触的域的表面的所有自由度的集合。  起始值 $\hat U^0 \dealcoloneq
 P_{\mathcal{A}_k}(0)$ 满足我们的障碍条件，也就是说，我们将初始零位移投射到可行位移集合上。

   <li>  组装牛顿矩阵  $A_{pq} \dealcoloneq a'(
 U^{i-1};\varphi_p,\varphi_q)$  和右侧  $F(\hat U^{i-1})$  。  这些对应于线性化的牛顿步骤，暂时忽略了接触不等式。

   <li>  找到满足@f{align*}
 A\tilde U^i + B\Lambda^i & = F, &\\
 \left[B^T\tilde U^i\right]_p & = G_p & \forall p\in\mathcal{A}_i,\\
 \Lambda^i_p & = 0 & \forall p\in\mathcal{F}_i.
 @f}的原始-双数对 $(\tilde U^i,\Lambda^i)$  。

如同步骤-41，我们可以通过消除第一个方程中 ${\cal A}_i$ 的那些自由度来获得这个问题的解决方案，并获得一个线性系统 $\hat {\hat A}(U^{i-1}) \tilde U^i = \hat {\hat H}(U^{i-1})$  。




   <li>  通过应用直线搜索和计算 $U^{i-1}$ 和 $\tilde U^i$ 的线性组合来减弱 $i>2$ 的牛顿迭代。这需要找到一个 $\alpha^i_l \dealcoloneq 2^{-l},(l=0,\ldots,10)$ ，以便@f{gather*}U^i \dealcoloneq \alpha^i_l\bar U^i +
 (1-\alpha^i_l)U^{i-1}@f}。

满足@f{gather*}
   \vert {\hat R}\left({\mathbf u}^{i}\right) \vert < \vert {\hat R}\left({\mathbf u}^{i-1}\right) \vert.
 \f}与 ${\hat R}\left({\mathbf u}\right)=\left(P_{Pi}(C\varepsilon(u)),\varepsilon(\varphi^{i}_p\right)$ ，除了(i)元素 $p\in\mathcal{A}_i$ ，我们设置 ${\hat R}\left({\mathbf u}\right)=0$ ，和(ii)对应于悬挂节点的元素，我们以通常方式消除。

   <li>  通过@f{gather*}\mathcal{A}_{i+1} \dealcoloneq \lbrace p\in\mathcal{S}:\Lambda^i_p +
 c\left(\left[B^TU^i\right]_p - G_p\right) > 0\rbrace,@f}定义新的活动和非活动集。

@f{gather*}\mathcal{F}_{i+1} \dealcoloneq \lbrace p\in\mathcal{S}:\Lambda^i_p +
 c\left(\left[B^TU^i\right]_p - G_p\right) \leq 0\rbrace.@f}



   <li> 项目 $U^i$ ，使其满足接触不等式，@f{gather*}\hat U^i \dealcoloneq P_{\mathcal{A}_{i+1}}(U^i).@f} 。

这里， $P_{\mathcal{A}}(U)$ 是 $\mathcal{A}$ 中的活性成分对间隙@f{gather*}P_{\mathcal{A}}(U)_p \dealcoloneq \begin{cases}
 U_p, & \textrm{if}\quad p\notin\mathcal{A}\\
 g_{h,p}, & \textrm{if}\quad
 p\in\mathcal{A},
 \end{cases}@f}的投影。

其中 $g_{h,p}$ 是<i>gap</i>，表示障碍物与身体未位移配置的距离。

   <li>  如果 $\mathcal{A}_{i+1} = \mathcal{A}_k$ 和 $\left\|
 {\hat R}\left({\mathbf u}^{i}\right) \right\|_{\ell_2} < \delta$ 则停止，否则设置 $i=i+1$ 并转到步骤（1）。这一步确保我们只有在找到正确的活动集和塑性已经迭代到足够的精度时才停止迭代。   </ol> 

在这个算法的第3步中，矩阵 $B\in\mathbb{R}^{n\times m}$ ,  $n>m$ 描述了位移和拉格朗日乘数（接触力）的基数的耦合，在我们的情况下它不是二次的，因为 $\Lambda^k$ 只定义在 $\Gamma_C$ ，即可能发生接触的面。如文中所示，我们可以选择 $B$ 是一个每行只有一个条目的矩阵，（另见H&uuml;eber, Wohlmuth:A primal-dual active set strategy for non-linear multibody contact problems, Comput.Method Appl. Mech.Engrg.194, 2005, pp.3147-3166）。)矢量 $G$ 是由间隙 $g_h$ 的合适近似值定义的。

@f{gather*}G_p = \begin{cases}
g_{h,p}, & \text{if}\quad p\in\mathcal{S}\\
0, & \text{if}\quad p\notin\mathcal{S}.
\end{cases}@f}






<a name="Adaptivemeshrefinement"></a><h3>Adaptive mesh refinement</h3>


由于我们的程序是在三维空间中运行的，所以程序执行的计算很昂贵。因此，使用自适应网格细化是在可接受的运行时间内的一个重要步骤。为了使我们的生活更轻松，我们简单地选择已经在deal.II中实现的KellyErrorEstimator。我们把包含位移 $u$ 的解向量交给它。正如我们将在结果中看到的，它产生了一个相当合理的接触区和塑性的自适应网格。




<a name="Implementation"></a><h3>Implementation</h3>


本教程实质上是步骤40和步骤41的混合体，但我们没有使用PETSc，而是让Trilinos库来处理线性代数的并行化问题（就像步骤32一样）。由于我们试图解决一个类似于步骤41的问题，我们将使用同样的方法，但现在是并行的。

一个困难是处理来自Dirichlet条件的约束，悬挂节点和由接触产生的不平等条件。为此，我们创建了三个AffineConstraints类型的对象，它们描述了各种约束条件，我们将在每次迭代中适当地组合它们。

与第41步相比，该计划有一些新的课程。

 <ul>   <li>   <code>ConstitutiveLaw</code>  描述材料的塑性行为。

 <li>   <code>SphereObstacle</code> 描述一个球体，作为被推入可变形弹性体的障碍物。   是用这个还是下一个类来描述障碍物，由输入参数文件决定。

 <li>   <code>ChineseObstacle</code> （和一个辅助类）是一个允许我们从一个文件中读入障碍物的类。在我们将在结果部分展示的例子中，这个文件将是 <code>'obstacle_file.dat'</code> ，并对应于显示力或力量的中文、日文或韩文符号的数据（见http://www.orientaloutpost.com/："这个词可用于激励--它也可以指力量/运动/推进/力。它可以是任何使你继续前进的内部或外部事物。这是用中文表达动机的最安全方式。如果你的听众是日本人，请看另一个关于动机的条目。这是日语和韩语中的一个词，但它的意思是 "动力 "或 "动能"（没有你可能正在寻找的动机的意思）"）。实质上，我们将假装有一个印章（即对应于平底障碍物的面具，没有中间高度的碎片），我们把它压在身体里。有关的符号看起来如下（也可参见本节顶部的图片，了解最终结果是怎样的）。

    <img src="https://www.dealii.org/images/steps/developer/step-42.character.png" alt="" width="25%">   </ul> 。

除此以外，让我们只对以下方面进行评论。   <ul>   <li>  程序允许你通过参数文件从两个不同的粗略网格中进行选择。这些是立方体 $[0,1]^3$ 或半球体，其开放面朝向正 $z$ 方向。

 <li> 在这两种情况下，我们将假设可能与障碍物接触的边界部分具有边界指标一的惯例。对于这两种网格，我们假定这是一个自由表面，即身体要么在那里接触，要么没有力作用在它身上。对于半球体，弯曲部分的边界指标为零，我们在那里施加零位移。对于盒子，我们沿底部施加零位移，但允许沿边的垂直位移（尽管没有水平位移）。   </ul> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * The set of include files is not much of a surprise any more at this time:
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/sparsity_tools.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/block_sparsity_pattern.h>
 * #include <deal.II/lac/solver_bicgstab.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/trilinos_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h>
 * #include <deal.II/lac/trilinos_vector.h>
 * #include <deal.II/lac/trilinos_parallel_block_vector.h>
 * #include <deal.II/lac/trilinos_precondition.h>
 * #include <deal.II/lac/trilinos_solver.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/manifold_lib.h>
 * 
 * #include <deal.II/distributed/tria.h>
 * #include <deal.II/distributed/grid_refinement.h>
 * #include <deal.II/distributed/solution_transfer.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/numerics/vector_tools.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/fe_field_function.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * Finally, we include two system headers that let us create a directory for
 * output files. The first header provides the <code>mkdir</code> function and
 * the second lets us determine what happened if <code>mkdir</code> fails.
 * 
 * @code
 * #include <sys/stat.h>
 * #include <cerrno>
 * 
 * namespace Step42
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeConstitutiveLawcodeclasstemplate"></a> 
 * <h3>The <code>ConstitutiveLaw</code> class template</h3>
 * 

 * 
 * This class provides an interface for a constitutive law, i.e., for the
 * relationship between strain $\varepsilon(\mathbf u)$ and stress
 * $\sigma$. In this example we are using an elastoplastic material behavior
 * with linear, isotropic hardening. Such materials are characterized by
 * Young's modulus $E$, Poisson's ratio $\nu$, the initial yield stress
 * $\sigma_0$ and the isotropic hardening parameter $\gamma$.  For $\gamma =
 * 0$ we obtain perfect elastoplastic behavior.
 *   

 * 
 * As explained in the paper that describes this program, the first Newton
 * steps are solved with a completely elastic material model to avoid having
 * to deal with both nonlinearities (plasticity and contact) at once. To this
 * end, this class has a function <code>set_sigma_0()</code> that we use later
 * on to simply set $\sigma_0$ to a very large value -- essentially
 * guaranteeing that the actual stress will not exceed it, and thereby
 * producing an elastic material. When we are ready to use a plastic model, we
 * set $\sigma_0$ back to its proper value, using the same function.  As a
 * result of this approach, we need to leave <code>sigma_0</code> as the only
 * non-const member variable of this class.
 * 
 * @code
 *   template <int dim>
 *   class ConstitutiveLaw
 *   {
 *   public:
 *     ConstitutiveLaw(const double E,
 *                     const double nu,
 *                     const double sigma_0,
 *                     const double gamma);
 * 
 *     void set_sigma_0(double sigma_zero);
 * 
 *     bool get_stress_strain_tensor(
 *       const SymmetricTensor<2, dim> &strain_tensor,
 *       SymmetricTensor<4, dim> &      stress_strain_tensor) const;
 * 
 *     void get_linearized_stress_strain_tensors(
 *       const SymmetricTensor<2, dim> &strain_tensor,
 *       SymmetricTensor<4, dim> &      stress_strain_tensor_linearized,
 *       SymmetricTensor<4, dim> &      stress_strain_tensor) const;
 * 
 *   private:
 *     const double kappa;
 *     const double mu;
 *     double       sigma_0;
 *     const double gamma;
 * 
 *     const SymmetricTensor<4, dim> stress_strain_tensor_kappa;
 *     const SymmetricTensor<4, dim> stress_strain_tensor_mu;
 *   };
 * 
 * @endcode
 * 
 * The constructor of the ConstitutiveLaw class sets the required material
 * parameter for our deformable body. Material parameters for elastic
 * isotropic media can be defined in a variety of ways, such as the pair $E,
 * \nu$ (elastic modulus and Poisson's number), using the Lam&eacute;
 * parameters
 * $\lambda,mu$ or several other commonly used conventions. Here, the
 * constructor takes a description of material parameters in the form of
 * $E,\nu$, but since this turns out to these are not the coefficients that
 * appear in the equations of the plastic projector, we immediately convert
 * them into the more suitable set $\kappa,\mu$ of bulk and shear moduli.  In
 * addition, the constructor takes $\sigma_0$ (the yield stress absent any
 * plastic strain) and $\gamma$ (the hardening parameter) as arguments. In
 * this constructor, we also compute the two principal components of the
 * stress-strain relation and its linearization.
 * 
 * @code
 *   template <int dim>
 *   ConstitutiveLaw<dim>::ConstitutiveLaw(double E,
 *                                         double nu,
 *                                         double sigma_0,
 *                                         double gamma)
 *     : kappa(E / (3 * (1 - 2 * nu)))
 *     , mu(E / (2 * (1 + nu)))
 *     , sigma_0(sigma_0)
 *     , gamma(gamma)
 *     , stress_strain_tensor_kappa(kappa *
 *                                  outer_product(unit_symmetric_tensor<dim>(),
 *                                                unit_symmetric_tensor<dim>()))
 *     , stress_strain_tensor_mu(
 *         2 * mu *
 *         (identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(),
 *                                                 unit_symmetric_tensor<dim>()) /
 *                                     3.0))
 *   {}
 * 
 * 
 *   template <int dim>
 *   void ConstitutiveLaw<dim>::set_sigma_0(double sigma_zero)
 *   {
 *     sigma_0 = sigma_zero;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstitutiveLawget_stress_strain_tensor"></a> 
 * <h4>ConstitutiveLaw::get_stress_strain_tensor</h4>
 * 

 * 
 * This is the principal component of the constitutive law. It
 * computes the fourth order symmetric tensor that relates the
 * strain to the stress according to the projection given above,
 * when evaluated at a particular strain point. We need this
 * function to calculate the nonlinear residual in
 * <code>PlasticityContactProblem::residual_nl_system()</code> where
 * we multiply this tensor with the strain given in a quadrature
 * point. The computations follow the formulas laid out in the
 * introduction. In comparing the formulas there with the
 * implementation below, recall that $C_\mu : \varepsilon = \tau_D$
 * and that $C_\kappa : \varepsilon = \kappa
 * \text{trace}(\varepsilon) I = \frac 13 \text{trace}(\tau) I$.
 *   

 * 
 * The function returns whether the quadrature point is plastic to allow for
 * some statistics downstream on how many of the quadrature points are
 * plastic and how many are elastic.
 * 
 * @code
 *   template <int dim>
 *   bool ConstitutiveLaw<dim>::get_stress_strain_tensor(
 *     const SymmetricTensor<2, dim> &strain_tensor,
 *     SymmetricTensor<4, dim> &      stress_strain_tensor) const
 *   {
 *     Assert(dim == 3, ExcNotImplemented());
 * 
 *     SymmetricTensor<2, dim> stress_tensor;
 *     stress_tensor =
 *       (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor;
 * 
 *     const SymmetricTensor<2, dim> deviator_stress_tensor =
 *       deviator(stress_tensor);
 *     const double deviator_stress_tensor_norm = deviator_stress_tensor.norm();
 * 
 *     stress_strain_tensor = stress_strain_tensor_mu;
 *     if (deviator_stress_tensor_norm > sigma_0)
 *       {
 *         const double beta = sigma_0 / deviator_stress_tensor_norm;
 *         stress_strain_tensor *= (gamma + (1 - gamma) * beta);
 *       }
 * 
 *     stress_strain_tensor += stress_strain_tensor_kappa;
 * 
 *     return (deviator_stress_tensor_norm > sigma_0);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstitutiveLawget_linearized_stress_strain_tensors"></a> 
 * <h4>ConstitutiveLaw::get_linearized_stress_strain_tensors</h4>
 * 

 * 
 * This function returns the linearized stress strain tensor, linearized
 * around the solution $u^{i-1}$ of the previous Newton step $i-1$.  The
 * parameter <code>strain_tensor</code> (commonly denoted
 * $\varepsilon(u^{i-1})$) must be passed as an argument, and serves as the
 * linearization point. The function returns the derivative of the nonlinear
 * constitutive law in the variable stress_strain_tensor, as well as the
 * stress-strain tensor of the linearized problem in
 * stress_strain_tensor_linearized.  See
 * PlasticityContactProblem::assemble_nl_system where this function is used.
 * 
 * @code
 *   template <int dim>
 *   void ConstitutiveLaw<dim>::get_linearized_stress_strain_tensors(
 *     const SymmetricTensor<2, dim> &strain_tensor,
 *     SymmetricTensor<4, dim> &      stress_strain_tensor_linearized,
 *     SymmetricTensor<4, dim> &      stress_strain_tensor) const
 *   {
 *     Assert(dim == 3, ExcNotImplemented());
 * 
 *     SymmetricTensor<2, dim> stress_tensor;
 *     stress_tensor =
 *       (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor;
 * 
 *     stress_strain_tensor            = stress_strain_tensor_mu;
 *     stress_strain_tensor_linearized = stress_strain_tensor_mu;
 * 
 *     SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor);
 *     const double deviator_stress_tensor_norm = deviator_stress_tensor.norm();
 * 
 *     if (deviator_stress_tensor_norm > sigma_0)
 *       {
 *         const double beta = sigma_0 / deviator_stress_tensor_norm;
 *         stress_strain_tensor *= (gamma + (1 - gamma) * beta);
 *         stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta);
 *         deviator_stress_tensor /= deviator_stress_tensor_norm;
 *         stress_strain_tensor_linearized -=
 *           (1 - gamma) * beta * 2 * mu *
 *           outer_product(deviator_stress_tensor, deviator_stress_tensor);
 *       }
 * 
 *     stress_strain_tensor += stress_strain_tensor_kappa;
 *     stress_strain_tensor_linearized += stress_strain_tensor_kappa;
 *   }
 * 
 * @endcode
 * 
 * <h3>Equation data: boundary forces, boundary values, obstacles</h3>
 *   

 * 
 * The following should be relatively standard. We need classes for
 * the boundary forcing term (which we here choose to be zero)
 * and boundary values on those part of the boundary that are not part
 * of the contact surface (also chosen to be zero here).
 * 
 * @code
 *   namespace EquationData
 *   {
 *     template <int dim>
 *     class BoundaryForce : public Function<dim>
 *     {
 *     public:
 *       BoundaryForce();
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  values) const override;
 *     };
 * 
 *     template <int dim>
 *     BoundaryForce<dim>::BoundaryForce()
 *       : Function<dim>(dim)
 *     {}
 * 
 * 
 *     template <int dim>
 *     double BoundaryForce<dim>::value(const Point<dim> &,
 *                                      const unsigned int) const
 *     {
 *       return 0.;
 *     }
 * 
 *     template <int dim>
 *     void BoundaryForce<dim>::vector_value(const Point<dim> &p,
 *                                           Vector<double> &  values) const
 *     {
 *       for (unsigned int c = 0; c < this->n_components; ++c)
 *         values(c) = BoundaryForce<dim>::value(p, c);
 *     }
 * 
 * 
 * 
 *     template <int dim>
 *     class BoundaryValues : public Function<dim>
 *     {
 *     public:
 *       BoundaryValues();
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 *     };
 * 
 * 
 *     template <int dim>
 *     BoundaryValues<dim>::BoundaryValues()
 *       : Function<dim>(dim)
 *     {}
 * 
 * 
 *     template <int dim>
 *     double BoundaryValues<dim>::value(const Point<dim> &,
 *                                       const unsigned int) const
 *     {
 *       return 0.;
 *     }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeSphereObstaclecodeclass"></a> 
 * <h4>The <code>SphereObstacle</code> class</h4>
 * 

 * 
 * The following class is the first of two obstacles that can be
 * selected from the input file. It describes a sphere centered
 * at position $x=y=0.5, z=z_{\text{surface}}+0.59$ and radius $r=0.6$,
 * where $z_{\text{surface}}$ is the vertical position of the (flat)
 * surface of the deformable body. The function's <code>value</code>
 * returns the location of the obstacle for a given $x,y$ value if the
 * point actually lies below the sphere, or a large positive value that
 * can't possibly interfere with the deformation if it lies outside
 * the "shadow" of the sphere.
 * 
 * @code
 *     template <int dim>
 *     class SphereObstacle : public Function<dim>
 *     {
 *     public:
 *       SphereObstacle(const double z_surface);
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  values) const override;
 * 
 *     private:
 *       const double z_surface;
 *     };
 * 
 * 
 *     template <int dim>
 *     SphereObstacle<dim>::SphereObstacle(const double z_surface)
 *       : Function<dim>(dim)
 *       , z_surface(z_surface)
 *     {}
 * 
 * 
 *     template <int dim>
 *     double SphereObstacle<dim>::value(const Point<dim> & p,
 *                                       const unsigned int component) const
 *     {
 *       if (component == 0)
 *         return p(0);
 *       else if (component == 1)
 *         return p(1);
 *       else if (component == 2)
 *         {
 *           if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5) < 0.36)
 *             return (-std::sqrt(0.36 - (p(0) - 0.5) * (p(0) - 0.5) -
 *                                (p(1) - 0.5) * (p(1) - 0.5)) +
 *                     z_surface + 0.59);
 *           else
 *             return 1000;
 *         }
 * 
 *       Assert(false, ExcNotImplemented());
 *       return 1e9; // an unreasonable value; ignored in debug mode because of the
 * @endcode
 * 
 * preceding Assert
 * 
 * @code
 *     }
 * 
 * 
 *     template <int dim>
 *     void SphereObstacle<dim>::vector_value(const Point<dim> &p,
 *                                            Vector<double> &  values) const
 *     {
 *       for (unsigned int c = 0; c < this->n_components; ++c)
 *         values(c) = SphereObstacle<dim>::value(p, c);
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBitmapFilecodeandcodeChineseObstaclecodeclasses"></a> 
 * <h4>The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes</h4>
 * 

 * 
 * The following two classes describe the obstacle outlined in the
 * introduction, i.e., the Chinese character. The first of the two,
 * <code>BitmapFile</code> is responsible for reading in data from a picture
 * file stored in pbm ascii format. This data will be bilinearly
 * interpolated and thereby provides a function that describes the obstacle.
 * (The code below shows how one can construct a function by interpolating
 * between given data points. One could use the
 * Functions::InterpolatedUniformGridData, introduced after this tutorial
 * program was written, which does exactly what we want here, but it is
 * instructive to see how to do it by hand.)
 *     

 * 
 * The data which we read from the file will be stored in a double
 * std::vector named obstacle_data.  This vector composes the base to
 * calculate a piecewise bilinear function as a polynomial interpolation.
 * The data we will read from a file consists of zeros (white) and ones
 * (black).
 *     

 * 
 * The <code>hx,hy</code> variables denote the spacing between pixels in $x$
 * and $y$ directions. <code>nx,ny</code> are the numbers of pixels in each
 * of these directions.  <code>get_value()</code> returns the value of the
 * image at a given location, interpolated from the adjacent pixel values.
 * 
 * @code
 *     template <int dim>
 *     class BitmapFile
 *     {
 *     public:
 *       BitmapFile(const std::string &name);
 * 
 *       double get_value(const double x, const double y) const;
 * 
 *     private:
 *       std::vector<double> obstacle_data;
 *       double              hx, hy;
 *       int                 nx, ny;
 * 
 *       double get_pixel_value(const int i, const int j) const;
 *     };
 * 
 * @endcode
 * 
 * The constructor of this class reads in the data that describes
 * the obstacle from the given file name.
 * 
 * @code
 *     template <int dim>
 *     BitmapFile<dim>::BitmapFile(const std::string &name)
 *       : obstacle_data(0)
 *       , hx(0)
 *       , hy(0)
 *       , nx(0)
 *       , ny(0)
 *     {
 *       std::ifstream f(name);
 *       AssertThrow(f,
 *                   ExcMessage(std::string("Can't read from file <") + name +
 *                              ">!"));
 * 
 *       std::string temp;
 *       f >> temp >> nx >> ny;
 * 
 *       AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format."));
 * 
 *       for (int k = 0; k < nx * ny; ++k)
 *         {
 *           double val;
 *           f >> val;
 *           obstacle_data.push_back(val);
 *         }
 * 
 *       hx = 1.0 / (nx - 1);
 *       hy = 1.0 / (ny - 1);
 * 
 *       if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *         std::cout << "Read obstacle from file <" << name << ">" << std::endl
 *                   << "Resolution of the scanned obstacle picture: " << nx
 *                   << " x " << ny << std::endl;
 *     }
 * 
 * @endcode
 * 
 * The following two functions return the value of a given pixel with
 * coordinates $i,j$, which we identify with the values of a function
 * defined at positions <code>i*hx, j*hy</code>, and at arbitrary
 * coordinates $x,y$ where we do a bilinear interpolation between
 * point values returned by the first of the two functions. In the
 * second function, for each $x,y$, we first compute the (integer)
 * location of the nearest pixel coordinate to the bottom left of
 * $x,y$, and then compute the coordinates $\xi,\eta$ within this
 * pixel. We truncate both kinds of variables from both below
 * and above to avoid problems when evaluating the function outside
 * of its defined range as may happen due to roundoff errors.
 * 
 * @code
 *     template <int dim>
 *     double BitmapFile<dim>::get_pixel_value(const int i, const int j) const
 *     {
 *       assert(i >= 0 && i < nx);
 *       assert(j >= 0 && j < ny);
 *       return obstacle_data[nx * (ny - 1 - j) + i];
 *     }
 * 
 *     template <int dim>
 *     double BitmapFile<dim>::get_value(const double x, const double y) const
 *     {
 *       const int ix = std::min(std::max(static_cast<int>(x / hx), 0), nx - 2);
 *       const int iy = std::min(std::max(static_cast<int>(y / hy), 0), ny - 2);
 * 
 *       const double xi  = std::min(std::max((x - ix * hx) / hx, 1.), 0.);
 *       const double eta = std::min(std::max((y - iy * hy) / hy, 1.), 0.);
 * 
 *       return ((1 - xi) * (1 - eta) * get_pixel_value(ix, iy) +
 *               xi * (1 - eta) * get_pixel_value(ix + 1, iy) +
 *               (1 - xi) * eta * get_pixel_value(ix, iy + 1) +
 *               xi * eta * get_pixel_value(ix + 1, iy + 1));
 *     }
 * 
 * @endcode
 * 
 * Finally, this is the class that actually uses the class above. It
 * has a BitmapFile object as a member that describes the height of the
 * obstacle. As mentioned above, the BitmapFile class will provide us
 * with a mask, i.e., values that are either zero or one (and, if you
 * ask for locations between pixels, values that are interpolated between
 * zero and one). This class translates this to heights that are either
 * 0.001 below the surface of the deformable body (if the BitmapFile
 * class reports a one at this location) or 0.999 above the obstacle (if
 * the BitmapFile class reports a zero). The following function should then
 * be self-explanatory.
 * 
 * @code
 *     template <int dim>
 *     class ChineseObstacle : public Function<dim>
 *     {
 *     public:
 *       ChineseObstacle(const std::string &filename, const double z_surface);
 * 
 *       virtual double value(const Point<dim> & p,
 *                            const unsigned int component = 0) const override;
 * 
 *       virtual void vector_value(const Point<dim> &p,
 *                                 Vector<double> &  values) const override;
 * 
 *     private:
 *       const BitmapFile<dim> input_obstacle;
 *       double                z_surface;
 *     };
 * 
 * 
 *     template <int dim>
 *     ChineseObstacle<dim>::ChineseObstacle(const std::string &filename,
 *                                           const double       z_surface)
 *       : Function<dim>(dim)
 *       , input_obstacle(filename)
 *       , z_surface(z_surface)
 *     {}
 * 
 * 
 *     template <int dim>
 *     double ChineseObstacle<dim>::value(const Point<dim> & p,
 *                                        const unsigned int component) const
 *     {
 *       if (component == 0)
 *         return p(0);
 *       if (component == 1)
 *         return p(1);
 *       else if (component == 2)
 *         {
 *           if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0)
 *             return z_surface + 0.999 - input_obstacle.get_value(p(0), p(1));
 *         }
 * 
 *       Assert(false, ExcNotImplemented());
 *       return 1e9; // an unreasonable value; ignored in debug mode because of the
 * @endcode
 * 
 * preceding Assert
 * 
 * @code
 *     }
 * 
 *     template <int dim>
 *     void ChineseObstacle<dim>::vector_value(const Point<dim> &p,
 *                                             Vector<double> &  values) const
 *     {
 *       for (unsigned int c = 0; c < this->n_components; ++c)
 *         values(c) = ChineseObstacle<dim>::value(p, c);
 *     }
 *   } // namespace EquationData
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePlasticityContactProblemcodeclasstemplate"></a> 
 * <h3>The <code>PlasticityContactProblem</code> class template</h3>
 * 

 * 
 * This is the main class of this program and supplies all functions and
 * variables needed to describe the nonlinear contact problem. It is close to
 * step-41 but with some additional features like handling hanging nodes, a
 * Newton method, using Trilinos and p4est for parallel distributed computing.
 * To deal with hanging nodes makes life a bit more complicated since we need
 * another AffineConstraints object now. We create a Newton method for the
 * active set method for the contact situation and to handle the nonlinear
 * operator for the constitutive law.
 *   

 * 
 * The general layout of this class is very much like for most other tutorial
 * programs. To make our life a bit easier, this class reads a set of input
 * parameters from an input file. These parameters, using the ParameterHandler
 * class, are declared in the <code>declare_parameters</code> function (which
 * is static so that it can be called before we even create an object of the
 * current type), and a ParameterHandler object that has been used to read an
 * input file will then be passed to the constructor of this class.
 *   

 * 
 * The remaining member functions are by and large as we have seen in several
 * of the other tutorial programs, though with additions for the current
 * nonlinear system. We will comment on their purpose as we get to them
 * further below.
 * 
 * @code
 *   template <int dim>
 *   class PlasticityContactProblem
 *   {
 *   public:
 *     PlasticityContactProblem(const ParameterHandler &prm);
 * 
 *     void run();
 * 
 *     static void declare_parameters(ParameterHandler &prm);
 * 
 *   private:
 *     void make_grid();
 *     void setup_system();
 *     void compute_dirichlet_constraints();
 *     void update_solution_and_constraints();
 *     void
 *          assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);
 *     void assemble_newton_system(
 *       const TrilinosWrappers::MPI::Vector &linearization_point);
 *     void compute_nonlinear_residual(
 *       const TrilinosWrappers::MPI::Vector &linearization_point);
 *     void solve_newton_system();
 *     void solve_newton();
 *     void refine_grid();
 *     void move_mesh(const TrilinosWrappers::MPI::Vector &displacement) const;
 *     void output_results(const unsigned int current_refinement_cycle);
 * 
 *     void output_contact_force() const;
 * 
 * @endcode
 * 
 * As far as member variables are concerned, we start with ones that we use
 * to indicate the MPI universe this program runs on, a stream we use to let
 * exactly one processor produce output to the console (see step-17) and
 * a variable that is used to time the various sections of the program:
 * 
 * @code
 *     MPI_Comm           mpi_communicator;
 *     ConditionalOStream pcout;
 *     TimerOutput        computing_timer;
 * 
 * @endcode
 * 
 * The next group describes the mesh and the finite element space.
 * In particular, for this parallel program, the finite element
 * space has associated with it variables that indicate which degrees
 * of freedom live on the current processor (the index sets, see
 * also step-40 and the @ref distributed documentation module) as
 * well as a variety of constraints: those imposed by hanging nodes,
 * by Dirichlet boundary conditions, and by the active set of
 * contact nodes. Of the three AffineConstraints variables defined
 * here, the first only contains hanging node constraints, the
 * second also those associated with Dirichlet boundary conditions,
 * and the third these plus the contact constraints.
 *     

 * 
 * The variable <code>active_set</code> consists of those degrees
 * of freedom constrained by the contact, and we use
 * <code>fraction_of_plastic_q_points_per_cell</code> to keep
 * track of the fraction of quadrature points on each cell where
 * the stress equals the yield stress. The latter is only used to
 * create graphical output showing the plastic zone, but not for
 * any further computation; the variable is a member variable of
 * this class since the information is computed as a by-product
 * of computing the residual, but is used only much later. (Note
 * that the vector is a vector of length equal to the number of
 * active cells on the <i>local mesh</i>; it is never used to
 * exchange information between processors and can therefore be
 * a regular deal.II vector.)
 * 
 * @code
 *     const unsigned int                        n_initial_global_refinements;
 *     parallel::distributed::Triangulation<dim> triangulation;
 * 
 *     const unsigned int fe_degree;
 *     FESystem<dim>      fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 *     AffineConstraints<double> constraints_hanging_nodes;
 *     AffineConstraints<double> constraints_dirichlet_and_hanging_nodes;
 *     AffineConstraints<double> all_constraints;
 * 
 *     IndexSet      active_set;
 *     Vector<float> fraction_of_plastic_q_points_per_cell;
 * 
 * 
 * @endcode
 * 
 * The next block of variables corresponds to the solution
 * and the linear systems we need to form. In particular, this
 * includes the Newton matrix and right hand side; the vector
 * that corresponds to the residual (i.e., the Newton right hand
 * side) but from which we have not eliminated the various
 * constraints and that is used to determine which degrees of
 * freedom need to be constrained in the next iteration; and
 * a vector that corresponds to the diagonal of the $B$ matrix
 * briefly mentioned in the introduction and discussed in the
 * accompanying paper.
 * 
 * @code
 *     TrilinosWrappers::SparseMatrix newton_matrix;
 * 
 *     TrilinosWrappers::MPI::Vector solution;
 *     TrilinosWrappers::MPI::Vector newton_rhs;
 *     TrilinosWrappers::MPI::Vector newton_rhs_uncondensed;
 *     TrilinosWrappers::MPI::Vector diag_mass_matrix_vector;
 * 
 * @endcode
 * 
 * The next block contains the variables that describe the material
 * response:
 * 
 * @code
 *     const double         e_modulus, nu, gamma, sigma_0;
 *     ConstitutiveLaw<dim> constitutive_law;
 * 
 * @endcode
 * 
 * And then there is an assortment of other variables that are used
 * to identify the mesh we are asked to build as selected by the
 * parameter file, the obstacle that is being pushed into the
 * deformable body, the mesh refinement strategy, whether to transfer
 * the solution from one mesh to the next, and how many mesh
 * refinement cycles to perform. As possible, we mark these kinds
 * of variables as <code>const</code> to help the reader identify
 * which ones may or may not be modified later on (the output directory
 * being an exception -- it is never modified outside the constructor
 * but it is awkward to initialize in the member-initializer-list
 * following the colon in the constructor since there we have only
 * one shot at setting it; the same is true for the mesh refinement
 * criterion):
 * 
 * @code
 *     const std::string                          base_mesh;
 *     const std::shared_ptr<const Function<dim>> obstacle;
 * 
 *     struct RefinementStrategy
 *     {
 *       enum value
 *       {
 *         refine_global,
 *         refine_percentage,
 *         refine_fix_dofs
 *       };
 *     };
 *     typename RefinementStrategy::value refinement_strategy;
 * 
 *     const bool         transfer_solution;
 *     std::string        output_dir;
 *     const unsigned int n_refinement_cycles;
 *     unsigned int       current_refinement_cycle;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodePlasticityContactProblemcodeclass"></a> 
 * <h3>Implementation of the <code>PlasticityContactProblem</code> class</h3>
 * 

 * 
 * 
 * <a name="PlasticityContactProblemdeclare_parameters"></a> 
 * <h4>PlasticityContactProblem::declare_parameters</h4>
 * 

 * 
 * Let us start with the declaration of run-time parameters that can be
 * selected in the input file. These values will be read back in the
 * constructor of this class to initialize the member variables of this
 * class:
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::declare_parameters(ParameterHandler &prm)
 *   {
 *     prm.declare_entry(
 *       "polynomial degree",
 *       "1",
 *       Patterns::Integer(),
 *       "Polynomial degree of the FE_Q finite element space, typically 1 or 2.");
 *     prm.declare_entry("number of initial refinements",
 *                       "2",
 *                       Patterns::Integer(),
 *                       "Number of initial global mesh refinement steps before "
 *                       "the first computation.");
 *     prm.declare_entry(
 *       "refinement strategy",
 *       "percentage",
 *       Patterns::Selection("global|percentage"),
 *       "Mesh refinement strategy:\n"
 *       " global: one global refinement\n"
 *       " percentage: a fixed percentage of cells gets refined using the Kelly estimator.");
 *     prm.declare_entry("number of cycles",
 *                       "5",
 *                       Patterns::Integer(),
 *                       "Number of adaptive mesh refinement cycles to run.");
 *     prm.declare_entry(
 *       "obstacle",
 *       "sphere",
 *       Patterns::Selection("sphere|read from file"),
 *       "The name of the obstacle to use. This may either be 'sphere' if we should "
 *       "use a spherical obstacle, or 'read from file' in which case the obstacle "
 *       "will be read from a file named 'obstacle.pbm' that is supposed to be in "
 *       "ASCII PBM format.");
 *     prm.declare_entry(
 *       "output directory",
 *       "",
 *       Patterns::Anything(),
 *       "Directory for output files (graphical output and benchmark "
 *       "statistics). If empty, use the current directory.");
 *     prm.declare_entry(
 *       "transfer solution",
 *       "false",
 *       Patterns::Bool(),
 *       "Whether the solution should be used as a starting guess "
 *       "for the next finer mesh. If false, then the iteration starts at "
 *       "zero on every mesh.");
 *     prm.declare_entry("base mesh",
 *                       "box",
 *                       Patterns::Selection("box|half sphere"),
 *                       "Select the shape of the domain: 'box' or 'half sphere'");
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodePlasticityContactProblemcodeconstructor"></a> 
 * <h4>The <code>PlasticityContactProblem</code> constructor</h4>
 * 

 * 
 * Given the declarations of member variables as well as the
 * declarations of run-time parameters that are read from the input
 * file, there is nothing surprising in this constructor. In the body
 * we initialize the mesh refinement strategy and the output directory,
 * creating such a directory if necessary.
 * 
 * @code
 *   template <int dim>
 *   PlasticityContactProblem<dim>::PlasticityContactProblem(
 *     const ParameterHandler &prm)
 *     : mpi_communicator(MPI_COMM_WORLD)
 *     , pcout(std::cout,
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *     , computing_timer(MPI_COMM_WORLD,
 *                       pcout,
 *                       TimerOutput::never,
 *                       TimerOutput::wall_times)
 * 
 *     , n_initial_global_refinements(
 *         prm.get_integer("number of initial refinements"))
 *     , triangulation(mpi_communicator)
 *     , fe_degree(prm.get_integer("polynomial degree"))
 *     , fe(FE_Q<dim>(QGaussLobatto<1>(fe_degree + 1)), dim)
 *     , dof_handler(triangulation)
 * 
 *     , e_modulus(200000)
 *     , nu(0.3)
 *     , gamma(0.01)
 *     , sigma_0(400.0)
 *     , constitutive_law(e_modulus, nu, sigma_0, gamma)
 * 
 *     , base_mesh(prm.get("base mesh"))
 *     , obstacle(prm.get("obstacle") == "read from file" ?
 *                  static_cast<const Function<dim> *>(
 *                    new EquationData::ChineseObstacle<dim>(
 *                      "obstacle.pbm",
 *                      (base_mesh == "box" ? 1.0 : 0.5))) :
 *                  static_cast<const Function<dim> *>(
 *                    new EquationData::SphereObstacle<dim>(
 *                      base_mesh == "box" ? 1.0 : 0.5)))
 * 
 *     , transfer_solution(prm.get_bool("transfer solution"))
 *     , n_refinement_cycles(prm.get_integer("number of cycles"))
 *     , current_refinement_cycle(0)
 * 
 *   {
 *     std::string strat = prm.get("refinement strategy");
 *     if (strat == "global")
 *       refinement_strategy = RefinementStrategy::refine_global;
 *     else if (strat == "percentage")
 *       refinement_strategy = RefinementStrategy::refine_percentage;
 *     else
 *       AssertThrow(false, ExcNotImplemented());
 * 
 *     output_dir = prm.get("output directory");
 *     if (output_dir != "" && *(output_dir.rbegin()) != '/')
 *       output_dir += "/";
 * 
 * @endcode
 * 
 * If necessary, create a new directory for the output.
 * 
 * @code
 *     if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *       {
 *         const int ierr = mkdir(output_dir.c_str(), 0777);
 *         AssertThrow(ierr == 0 || errno == EEXIST, ExcIO());
 *       }
 * 
 *     pcout << "    Using output directory '" << output_dir << "'" << std::endl;
 *     pcout << "    FE degree " << fe_degree << std::endl;
 *     pcout << "    transfer solution " << (transfer_solution ? "true" : "false")
 *           << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemmake_grid"></a> 
 * <h4>PlasticityContactProblem::make_grid</h4>
 * 

 * 
 * The next block deals with constructing the starting mesh.
 * We will use the following helper function and the first
 * block of the <code>make_grid()</code> to construct a
 * mesh that corresponds to a half sphere. deal.II has a function
 * that creates such a mesh, but it is in the wrong location
 * and facing the wrong direction, so we need to shift and rotate
 * it a bit before using it.
 *   

 * 
 * For later reference, as described in the documentation of
 * GridGenerator::half_hyper_ball(), the flat surface of the halfsphere
 * has boundary indicator zero, while the remainder has boundary
 * indicator one.
 * 
 * @code
 *   Point<3> rotate_half_sphere(const Point<3> &in)
 *   {
 *     return {in(2), in(1), -in(0)};
 *   }
 * 
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::make_grid()
 *   {
 *     if (base_mesh == "half sphere")
 *       {
 *         const Point<dim> center(0, 0, 0);
 *         const double     radius = 0.8;
 *         GridGenerator::half_hyper_ball(triangulation, center, radius);
 * @endcode
 * 
 * Since we will attach a different manifold below, we immediately
 * clear the default manifold description:
 * 
 * @code
 *         triangulation.reset_all_manifolds();
 * 
 *         GridTools::transform(&rotate_half_sphere, triangulation);
 *         GridTools::shift(Point<dim>(0.5, 0.5, 0.5), triangulation);
 * 
 *         SphericalManifold<dim> manifold_description(Point<dim>(0.5, 0.5, 0.5));
 *         GridTools::copy_boundary_to_manifold_id(triangulation);
 *         triangulation.set_manifold(0, manifold_description);
 *       }
 * @endcode
 * 
 * Alternatively, create a hypercube mesh. After creating it,
 * assign boundary indicators as follows:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * >     _______
 * >    /  1    /|
 * >   /______ / |
 * >  |       | 8|
 * >  |   8   | /
 * >  |_______|/
 * >      6
 * @endcode
 * </div>
 * In other words, the boundary indicators of the sides of the cube are 8.
 * The boundary indicator of the bottom is 6 and the top has indicator 1.
 * We set these by looping over all cells of all faces and looking at
 * coordinate values of the cell center, and will make use of these
 * indicators later when evaluating which boundary will carry Dirichlet
 * boundary conditions or will be subject to potential contact.
 * (In the current case, the mesh contains only a single cell, and all of
 * its faces are on the boundary, so both the loop over all cells and the
 * query whether a face is on the boundary are, strictly speaking,
 * unnecessary; we retain them simply out of habit: this kind of code can
 * be found in many programs in essentially this form.)
 * 
 * @code
 *     else
 *       {
 *         const Point<dim> p1(0, 0, 0);
 *         const Point<dim> p2(1.0, 1.0, 1.0);
 * 
 *         GridGenerator::hyper_rectangle(triangulation, p1, p2);
 * 
 *         for (const auto &cell : triangulation.active_cell_iterators())
 *           for (const auto &face : cell->face_iterators())
 *             if (face->at_boundary())
 *               {
 *                 if (std::fabs(face->center()[2] - p2[2]) < 1e-12)
 *                   face->set_boundary_id(1);
 *                 if (std::fabs(face->center()[0] - p1[0]) < 1e-12 ||
 *                     std::fabs(face->center()[0] - p2[0]) < 1e-12 ||
 *                     std::fabs(face->center()[1] - p1[1]) < 1e-12 ||
 *                     std::fabs(face->center()[1] - p2[1]) < 1e-12)
 *                   face->set_boundary_id(8);
 *                 if (std::fabs(face->center()[2] - p1[2]) < 1e-12)
 *                   face->set_boundary_id(6);
 *               }
 *       }
 * 
 *     triangulation.refine_global(n_initial_global_refinements);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsetup_system"></a> 
 * <h4>PlasticityContactProblem::setup_system</h4>
 * 

 * 
 * The next piece in the puzzle is to set up the DoFHandler, resize
 * vectors and take care of various other status variables such as
 * index sets and constraint matrices.
 *   

 * 
 * In the following, each group of operations is put into a brace-enclosed
 * block that is being timed by the variable declared at the top of the
 * block (the constructor of the TimerOutput::Scope variable starts the
 * timed section, the destructor that is called at the end of the block
 * stops it again).
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::setup_system()
 *   {
 *     /* setup dofs and get index sets for locally owned and relevant dofs */
 *     {
 *       TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs");
 *       dof_handler.distribute_dofs(fe);
 * 
 *       locally_owned_dofs = dof_handler.locally_owned_dofs();
 *       locally_relevant_dofs.clear();
 *       DoFTools::extract_locally_relevant_dofs(dof_handler,
 *                                               locally_relevant_dofs);
 *     }
 * 
 *     /* setup hanging nodes and Dirichlet constraints */
 *     {
 *       TimerOutput::Scope t(computing_timer, "Setup: constraints");
 *       constraints_hanging_nodes.reinit(locally_relevant_dofs);
 *       DoFTools::make_hanging_node_constraints(dof_handler,
 *                                               constraints_hanging_nodes);
 *       constraints_hanging_nodes.close();
 * 
 *       pcout << "   Number of active cells: "
 *             << triangulation.n_global_active_cells() << std::endl
 *             << "   Number of degrees of freedom: " << dof_handler.n_dofs()
 *             << std::endl;
 * 
 *       compute_dirichlet_constraints();
 *     }
 * 
 *     /* initialization of vectors and the active set */
 *     {
 *       TimerOutput::Scope t(computing_timer, "Setup: vectors");
 *       solution.reinit(locally_relevant_dofs, mpi_communicator);
 *       newton_rhs.reinit(locally_owned_dofs, mpi_communicator);
 *       newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator);
 *       diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator);
 *       fraction_of_plastic_q_points_per_cell.reinit(
 *         triangulation.n_active_cells());
 * 
 *       active_set.clear();
 *       active_set.set_size(dof_handler.n_dofs());
 *     }
 * 
 * @endcode
 * 
 * Finally, we set up sparsity patterns and matrices.
 * We temporarily (ab)use the system matrix to also build the (diagonal)
 * matrix that we use in eliminating degrees of freedom that are in contact
 * with the obstacle, but we then immediately set the Newton matrix back
 * to zero.
 * 
 * @code
 *     {
 *       TimerOutput::Scope                t(computing_timer, "Setup: matrix");
 *       TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
 *                                            mpi_communicator);
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler,
 *                                       sp,
 *                                       constraints_dirichlet_and_hanging_nodes,
 *                                       false,
 *                                       Utilities::MPI::this_mpi_process(
 *                                         mpi_communicator));
 *       sp.compress();
 *       newton_matrix.reinit(sp);
 * 
 * 
 *       TrilinosWrappers::SparseMatrix &mass_matrix = newton_matrix;
 * 
 *       assemble_mass_matrix_diagonal(mass_matrix);
 * 
 *       const unsigned int start = (newton_rhs.local_range().first),
 *                          end   = (newton_rhs.local_range().second);
 *       for (unsigned int j = start; j < end; ++j)
 *         diag_mass_matrix_vector(j) = mass_matrix.diag_element(j);
 *       diag_mass_matrix_vector.compress(VectorOperation::insert);
 * 
 *       mass_matrix = 0;
 *     }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemcompute_dirichlet_constraints"></a> 
 * <h4>PlasticityContactProblem::compute_dirichlet_constraints</h4>
 * 

 * 
 * This function, broken out of the preceding one, computes the constraints
 * associated with Dirichlet-type boundary conditions and puts them into the
 * <code>constraints_dirichlet_and_hanging_nodes</code> variable by merging
 * with the constraints that come from hanging nodes.
 *   

 * 
 * As laid out in the introduction, we need to distinguish between two
 * cases:
 * - If the domain is a box, we set the displacement to zero at the bottom,
 * and allow vertical movement in z-direction along the sides. As
 * shown in the <code>make_grid()</code> function, the former corresponds
 * to boundary indicator 6, the latter to 8.
 * - If the domain is a half sphere, then we impose zero displacement along
 * the curved part of the boundary, associated with boundary indicator zero.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::compute_dirichlet_constraints()
 *   {
 *     constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs);
 *     constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes);
 * 
 *     if (base_mesh == "box")
 *       {
 * @endcode
 * 
 * interpolate all components of the solution
 * 
 * @code
 *         VectorTools::interpolate_boundary_values(
 *           dof_handler,
 *           6,
 *           EquationData::BoundaryValues<dim>(),
 *           constraints_dirichlet_and_hanging_nodes,
 *           ComponentMask());
 * 
 * @endcode
 * 
 * interpolate x- and y-components of the
 * solution (this is a bit mask, so apply
 * operator| )
 * 
 * @code
 *         const FEValuesExtractors::Scalar x_displacement(0);
 *         const FEValuesExtractors::Scalar y_displacement(1);
 *         VectorTools::interpolate_boundary_values(
 *           dof_handler,
 *           8,
 *           EquationData::BoundaryValues<dim>(),
 *           constraints_dirichlet_and_hanging_nodes,
 *           (fe.component_mask(x_displacement) |
 *            fe.component_mask(y_displacement)));
 *       }
 *     else
 *       VectorTools::interpolate_boundary_values(
 *         dof_handler,
 *         0,
 *         EquationData::BoundaryValues<dim>(),
 *         constraints_dirichlet_and_hanging_nodes,
 *         ComponentMask());
 * 
 *     constraints_dirichlet_and_hanging_nodes.close();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemassemble_mass_matrix_diagonal"></a> 
 * <h4>PlasticityContactProblem::assemble_mass_matrix_diagonal</h4>
 * 

 * 
 * The next helper function computes the (diagonal) mass matrix that
 * is used to determine the active set of the active set method we use in
 * the contact algorithm. This matrix is of mass matrix type, but unlike
 * the standard mass matrix, we can make it diagonal (even in the case of
 * higher order elements) by using a quadrature formula that has its
 * quadrature points at exactly the same locations as the interpolation points
 * for the finite element are located. We achieve this by using a
 * QGaussLobatto quadrature formula here, along with initializing the finite
 * element with a set of interpolation points derived from the same quadrature
 * formula. The remainder of the function is relatively straightforward: we
 * put the resulting matrix into the given argument; because we know the
 * matrix is diagonal, it is sufficient to have a loop over only $i$ and
 * not over $j$. Strictly speaking, we could even avoid multiplying the
 * shape function's values at quadrature point <code>q_point</code> by itself
 * because we know the shape value to be a vector with exactly one one which
 * when dotted with itself yields one. Since this function is not time
 * critical we add this term for clarity.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal(
 *     TrilinosWrappers::SparseMatrix &mass_matrix)
 *   {
 *     QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1);
 * 
 *     FEFaceValues<dim> fe_values_face(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary() && face->boundary_id() == 1)
 *             {
 *               fe_values_face.reinit(cell, face);
 *               cell_matrix = 0;
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                    ++q_point)
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                   cell_matrix(i, i) +=
 *                     (fe_values_face[displacement].value(i, q_point) *
 *                      fe_values_face[displacement].value(i, q_point) *
 *                      fe_values_face.JxW(q_point));
 * 
 *               cell->get_dof_indices(local_dof_indices);
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 mass_matrix.add(local_dof_indices[i],
 *                                 local_dof_indices[i],
 *                                 cell_matrix(i, i));
 *             }
 *     mass_matrix.compress(VectorOperation::add);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemupdate_solution_and_constraints"></a> 
 * <h4>PlasticityContactProblem::update_solution_and_constraints</h4>
 * 

 * 
 * The following function is the first function we call in each Newton
 * iteration in the <code>solve_newton()</code> function. What it does is
 * to project the solution onto the feasible set and update the active set
 * for the degrees of freedom that touch or penetrate the obstacle.
 *   

 * 
 * In order to function, we first need to do some bookkeeping: We need
 * to write into the solution vector (which we can only do with fully
 * distributed vectors without ghost elements) and we need to read
 * the Lagrange multiplier and the elements of the diagonal mass matrix
 * from their respective vectors (which we can only do with vectors that
 * do have ghost elements), so we create the respective vectors. We then
 * also initialize the constraints object that will contain constraints
 * from contact and all other sources, as well as an object that contains
 * an index set of all locally owned degrees of freedom that are part of
 * the contact:
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::update_solution_and_constraints()
 *   {
 *     std::vector<bool> dof_touched(dof_handler.n_dofs(), false);
 * 
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs,
 *                                                        mpi_communicator);
 *     distributed_solution = solution;
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs,
 *                                          mpi_communicator);
 *     lambda = newton_rhs_uncondensed;
 * 
 *     TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant(
 *       locally_relevant_dofs, mpi_communicator);
 *     diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;
 * 
 * 
 *     all_constraints.reinit(locally_relevant_dofs);
 *     active_set.clear();
 * 
 * @endcode
 * 
 * The second part is a loop over all cells in which we look at each
 * point where a degree of freedom is defined whether the active set
 * condition is true and we need to add this degree of freedom to
 * the active set of contact nodes. As we always do, if we want to
 * evaluate functions at individual points, we do this with an
 * FEValues object (or, here, an FEFaceValues object since we need to
 * check contact at the surface) with an appropriately chosen quadrature
 * object. We create this face quadrature object by choosing the
 * "support points" of the shape functions defined on the faces
 * of cells (for more on support points, see this
 * @ref GlossSupport "glossary entry"). As a consequence, we have as
 * many quadrature points as there are shape functions per face and
 * looping over quadrature points is equivalent to looping over shape
 * functions defined on a face. With this, the code looks as follows:
 * 
 * @code
 *     Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points());
 *     FEFaceValues<dim>   fe_values_face(fe,
 *                                      face_quadrature,
 *                                      update_quadrature_points);
 * 
 *     const unsigned int dofs_per_face   = fe.n_dofs_per_face();
 *     const unsigned int n_face_q_points = face_quadrature.size();
 * 
 *     std::vector<types::global_dof_index> dof_indices(dofs_per_face);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (!cell->is_artificial())
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary() && face->boundary_id() == 1)
 *             {
 *               fe_values_face.reinit(cell, face);
 *               face->get_dof_indices(dof_indices);
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                    ++q_point)
 *                 {
 * @endcode
 * 
 * At each quadrature point (i.e., at each support point of a
 * degree of freedom located on the contact boundary), we then
 * ask whether it is part of the z-displacement degrees of
 * freedom and if we haven't encountered this degree of
 * freedom yet (which can happen for those on the edges
 * between faces), we need to evaluate the gap between the
 * deformed object and the obstacle. If the active set
 * condition is true, then we add a constraint to the
 * AffineConstraints object that the next Newton update needs
 * to satisfy, set the solution vector's corresponding element
 * to the correct value, and add the index to the IndexSet
 * object that stores which degree of freedom is part of the
 * contact:
 * 
 * @code
 *                   const unsigned int component =
 *                     fe.face_system_to_component_index(q_point).first;
 * 
 *                   const unsigned int index_z = dof_indices[q_point];
 * 
 *                   if ((component == 2) && (dof_touched[index_z] == false))
 *                     {
 *                       dof_touched[index_z] = true;
 * 
 *                       const Point<dim> this_support_point =
 *                         fe_values_face.quadrature_point(q_point);
 * 
 *                       const double obstacle_value =
 *                         obstacle->value(this_support_point, 2);
 *                       const double solution_here = solution(index_z);
 *                       const double undeformed_gap =
 *                         obstacle_value - this_support_point(2);
 * 
 *                       const double c = 100.0 * e_modulus;
 *                       if ((lambda(index_z) /
 *                                diag_mass_matrix_vector_relevant(index_z) +
 *                              c * (solution_here - undeformed_gap) >
 *                            0) &&
 *                           !constraints_hanging_nodes.is_constrained(index_z))
 *                         {
 *                           all_constraints.add_line(index_z);
 *                           all_constraints.set_inhomogeneity(index_z,
 *                                                             undeformed_gap);
 *                           distributed_solution(index_z) = undeformed_gap;
 * 
 *                           active_set.add_index(index_z);
 *                         }
 *                     }
 *                 }
 *             }
 * 
 * @endcode
 * 
 * At the end of this function, we exchange data between processors updating
 * those ghost elements in the <code>solution</code> variable that have been
 * written by other processors. We then merge the Dirichlet constraints and
 * those from hanging nodes into the AffineConstraints object that already
 * contains the active set. We finish the function by outputting the total
 * number of actively constrained degrees of freedom for which we sum over
 * the number of actively constrained degrees of freedom owned by each
 * of the processors. This number of locally owned constrained degrees of
 * freedom is of course the number of elements of the intersection of the
 * active set and the set of locally owned degrees of freedom, which
 * we can get by using <code>operator&</code> on two IndexSets:
 * 
 * @code
 *     distributed_solution.compress(VectorOperation::insert);
 *     solution = distributed_solution;
 * 
 *     all_constraints.close();
 *     all_constraints.merge(constraints_dirichlet_and_hanging_nodes);
 * 
 *     pcout << "         Size of active set: "
 *           << Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(),
 *                                  mpi_communicator)
 *           << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemassemble_newton_system"></a> 
 * <h4>PlasticityContactProblem::assemble_newton_system</h4>
 * 

 * 
 * Given the complexity of the problem, it may come as a bit of a surprise
 * that assembling the linear system we have to solve in each Newton iteration
 * is actually fairly straightforward. The following function builds the
 * Newton right hand side and Newton matrix. It looks fairly innocent because
 * the heavy lifting happens in the call to
 * <code>ConstitutiveLaw::get_linearized_stress_strain_tensors()</code> and in
 * particular in AffineConstraints::distribute_local_to_global(), using the
 * constraints we have previously computed.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::assemble_newton_system(
 *     const TrilinosWrappers::MPI::Vector &linearization_point)
 *   {
 *     TimerOutput::Scope t(computing_timer, "Assembling");
 * 
 *     QGauss<dim>     quadrature_formula(fe.degree + 1);
 *     QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_values_face(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     const EquationData::BoundaryForce<dim> boundary_force;
 *     std::vector<Vector<double>> boundary_force_values(n_face_q_points,
 *                                                       Vector<double>(dim));
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<double>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit(cell);
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 *           std::vector<SymmetricTensor<2, dim>> strain_tensor(n_q_points);
 *           fe_values[displacement].get_function_symmetric_gradients(
 *             linearization_point, strain_tensor);
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             {
 *               SymmetricTensor<4, dim> stress_strain_tensor_linearized;
 *               SymmetricTensor<4, dim> stress_strain_tensor;
 *               constitutive_law.get_linearized_stress_strain_tensors(
 *                 strain_tensor[q_point],
 *                 stress_strain_tensor_linearized,
 *                 stress_strain_tensor);
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 * @endcode
 * 
 * Having computed the stress-strain tensor and its
 * linearization, we can now put together the parts of the
 * matrix and right hand side. In both, we need the linearized
 * stress-strain tensor times the symmetric gradient of
 * $\varphi_i$, i.e. the term $I_\Pi\varepsilon(\varphi_i)$,
 * so we introduce an abbreviation of this term. Recall that
 * the matrix corresponds to the bilinear form
 * $A_{ij}=(I_\Pi\varepsilon(\varphi_i),\varepsilon(\varphi_j))$
 * in the notation of the accompanying publication, whereas
 * the right hand side is $F_i=([I_\Pi-P_\Pi
 * C]\varepsilon(\varphi_i),\varepsilon(\mathbf u))$ where $u$
 * is the current linearization points (typically the last
 * solution). This might suggest that the right hand side will
 * be zero if the material is completely elastic (where
 * $I_\Pi=P_\Pi$) but this ignores the fact that the right
 * hand side will also contain contributions from
 * non-homogeneous constraints due to the contact.
 *                   

 * 
 * The code block that follows this adds contributions that
 * are due to boundary forces, should there be any.
 * 
 * @code
 *                   const SymmetricTensor<2, dim> stress_phi_i =
 *                     stress_strain_tensor_linearized *
 *                     fe_values[displacement].symmetric_gradient(i, q_point);
 * 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     cell_matrix(i, j) +=
 *                       (stress_phi_i *
 *                        fe_values[displacement].symmetric_gradient(j, q_point) *
 *                        fe_values.JxW(q_point));
 * 
 *                   cell_rhs(i) +=
 *                     ((stress_phi_i -
 *                       stress_strain_tensor *
 *                         fe_values[displacement].symmetric_gradient(i,
 *                                                                    q_point)) *
 *                      strain_tensor[q_point] * fe_values.JxW(q_point));
 *                 }
 *             }
 * 
 *           for (const auto &face : cell->face_iterators())
 *             if (face->at_boundary() && face->boundary_id() == 1)
 *               {
 *                 fe_values_face.reinit(cell, face);
 * 
 *                 boundary_force.vector_value_list(
 *                   fe_values_face.get_quadrature_points(),
 *                   boundary_force_values);
 * 
 *                 for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                      ++q_point)
 *                   {
 *                     Tensor<1, dim> rhs_values;
 *                     rhs_values[2] = boundary_force_values[q_point][2];
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                       cell_rhs(i) +=
 *                         (fe_values_face[displacement].value(i, q_point) *
 *                          rhs_values * fe_values_face.JxW(q_point));
 *                   }
 *               }
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           all_constraints.distribute_local_to_global(cell_matrix,
 *                                                      cell_rhs,
 *                                                      local_dof_indices,
 *                                                      newton_matrix,
 *                                                      newton_rhs,
 *                                                      true);
 *         }
 * 
 *     newton_matrix.compress(VectorOperation::add);
 *     newton_rhs.compress(VectorOperation::add);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemcompute_nonlinear_residual"></a> 
 * <h4>PlasticityContactProblem::compute_nonlinear_residual</h4>
 * 

 * 
 * The following function computes the nonlinear residual of the equation
 * given the current solution (or any other linearization point). This
 * is needed in the linear search algorithm where we need to try various
 * linear combinations of previous and current (trial) solution to
 * compute the (real, globalized) solution of the current Newton step.
 *   

 * 
 * That said, in a slight abuse of the name of the function, it actually
 * does significantly more. For example, it also computes the vector
 * that corresponds to the Newton residual but without eliminating
 * constrained degrees of freedom. We need this vector to compute contact
 * forces and, ultimately, to compute the next active set. Likewise, by
 * keeping track of how many quadrature points we encounter on each cell
 * that show plastic yielding, we also compute the
 * <code>fraction_of_plastic_q_points_per_cell</code> vector that we
 * can later output to visualize the plastic zone. In both of these cases,
 * the results are not necessary as part of the line search, and so we may
 * be wasting a small amount of time computing them. At the same time, this
 * information appears as a natural by-product of what we need to do here
 * anyway, and we want to collect it once at the end of each Newton
 * step, so we may as well do it here.
 *   

 * 
 * The actual implementation of this function should be rather obvious:
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::compute_nonlinear_residual(
 *     const TrilinosWrappers::MPI::Vector &linearization_point)
 *   {
 *     QGauss<dim>     quadrature_formula(fe.degree + 1);
 *     QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);
 * 
 *     FEValues<dim> fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_JxW_values);
 * 
 *     FEFaceValues<dim> fe_values_face(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_quadrature_points |
 *                                        update_JxW_values);
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points      = quadrature_formula.size();
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     const EquationData::BoundaryForce<dim> boundary_force;
 *     std::vector<Vector<double>> boundary_force_values(n_face_q_points,
 *                                                       Vector<double>(dim));
 * 
 *     Vector<double> cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     newton_rhs             = 0;
 *     newton_rhs_uncondensed = 0;
 * 
 *     fraction_of_plastic_q_points_per_cell = 0;
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit(cell);
 *           cell_rhs = 0;
 * 
 *           std::vector<SymmetricTensor<2, dim>> strain_tensors(n_q_points);
 *           fe_values[displacement].get_function_symmetric_gradients(
 *             linearization_point, strain_tensors);
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *             {
 *               SymmetricTensor<4, dim> stress_strain_tensor;
 *               const bool              q_point_is_plastic =
 *                 constitutive_law.get_stress_strain_tensor(
 *                   strain_tensors[q_point], stress_strain_tensor);
 *               if (q_point_is_plastic)
 *                 ++fraction_of_plastic_q_points_per_cell(
 *                   cell->active_cell_index());
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   cell_rhs(i) -=
 *                     (strain_tensors[q_point] * stress_strain_tensor *
 *                      fe_values[displacement].symmetric_gradient(i, q_point) *
 *                      fe_values.JxW(q_point));
 * 
 *                   Tensor<1, dim> rhs_values;
 *                   rhs_values = 0;
 *                   cell_rhs(i) += (fe_values[displacement].value(i, q_point) *
 *                                   rhs_values * fe_values.JxW(q_point));
 *                 }
 *             }
 * 
 *           for (const auto &face : cell->face_iterators())
 *             if (face->at_boundary() && face->boundary_id() == 1)
 *               {
 *                 fe_values_face.reinit(cell, face);
 * 
 *                 boundary_force.vector_value_list(
 *                   fe_values_face.get_quadrature_points(),
 *                   boundary_force_values);
 * 
 *                 for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                      ++q_point)
 *                   {
 *                     Tensor<1, dim> rhs_values;
 *                     rhs_values[2] = boundary_force_values[q_point][2];
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                       cell_rhs(i) +=
 *                         (fe_values_face[displacement].value(i, q_point) *
 *                          rhs_values * fe_values_face.JxW(q_point));
 *                   }
 *               }
 * 
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints_dirichlet_and_hanging_nodes.distribute_local_to_global(
 *             cell_rhs, local_dof_indices, newton_rhs);
 * 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *             newton_rhs_uncondensed(local_dof_indices[i]) += cell_rhs(i);
 *         }
 * 
 *     fraction_of_plastic_q_points_per_cell /= quadrature_formula.size();
 *     newton_rhs.compress(VectorOperation::add);
 *     newton_rhs_uncondensed.compress(VectorOperation::add);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsolve_newton_system"></a> 
 * <h4>PlasticityContactProblem::solve_newton_system</h4>
 * 

 * 
 * The last piece before we can discuss the actual Newton iteration
 * on a single mesh is the solver for the linear systems. There are
 * a couple of complications that slightly obscure the code, but
 * mostly it is just setup then solve. Among the complications are:
 *   

 * 
 * - For the hanging nodes we have to apply
 * the AffineConstraints::set_zero function to newton_rhs.
 * This is necessary if a hanging node with solution value $x_0$
 * has one neighbor with value $x_1$ which is in contact with the
 * obstacle and one neighbor $x_2$ which is not in contact. Because
 * the update for the former will be prescribed, the hanging node constraint
 * will have an inhomogeneity and will look like $x_0 = x_1/2 +
 * \text{gap}/2$. So the corresponding entries in the right-hand-side are
 * non-zero with a meaningless value. These values we have to set to zero.
 * - Like in step-40, we need to shuffle between vectors that do and do
 * not have ghost elements when solving or using the solution.
 *   

 * 
 * The rest of the function is similar to step-40 and
 * step-41 except that we use a BiCGStab solver
 * instead of CG. This is due to the fact that for very small hardening
 * parameters $\gamma$, the linear system becomes almost semidefinite though
 * still symmetric. BiCGStab appears to have an easier time with such linear
 * systems.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::solve_newton_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "Solve");
 * 
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs,
 *                                                        mpi_communicator);
 *     distributed_solution = solution;
 * 
 *     constraints_hanging_nodes.set_zero(distributed_solution);
 *     constraints_hanging_nodes.set_zero(newton_rhs);
 * 
 *     TrilinosWrappers::PreconditionAMG preconditioner;
 *     {
 *       TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner");
 * 
 *       std::vector<std::vector<bool>> constant_modes;
 *       DoFTools::extract_constant_modes(dof_handler,
 *                                        ComponentMask(),
 *                                        constant_modes);
 * 
 *       TrilinosWrappers::PreconditionAMG::AdditionalData additional_data;
 *       additional_data.constant_modes        = constant_modes;
 *       additional_data.elliptic              = true;
 *       additional_data.n_cycles              = 1;
 *       additional_data.w_cycle               = false;
 *       additional_data.output_details        = false;
 *       additional_data.smoother_sweeps       = 2;
 *       additional_data.aggregation_threshold = 1e-2;
 * 
 *       preconditioner.initialize(newton_matrix, additional_data);
 *     }
 * 
 *     {
 *       TimerOutput::Scope t(computing_timer, "Solve: iterate");
 * 
 *       TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
 * 
 *       const double relative_accuracy = 1e-8;
 *       const double solver_tolerance =
 *         relative_accuracy *
 *         newton_matrix.residual(tmp, distributed_solution, newton_rhs);
 * 
 *       SolverControl solver_control(newton_matrix.m(), solver_tolerance);
 *       SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control);
 *       solver.solve(newton_matrix,
 *                    distributed_solution,
 *                    newton_rhs,
 *                    preconditioner);
 * 
 *       pcout << "         Error: " << solver_control.initial_value() << " -> "
 *             << solver_control.last_value() << " in "
 *             << solver_control.last_step() << " Bicgstab iterations."
 *             << std::endl;
 *     }
 * 
 *     all_constraints.distribute(distributed_solution);
 * 
 *     solution = distributed_solution;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsolve_newton"></a> 
 * <h4>PlasticityContactProblem::solve_newton</h4>
 * 

 * 
 * This is, finally, the function that implements the damped Newton method
 * on the current mesh. There are two nested loops: the outer loop for the
 * Newton iteration and the inner loop for the line search which will be used
 * only if necessary. To obtain a good and reasonable starting value we solve
 * an elastic problem in the very first Newton step on each mesh (or only on
 * the first mesh if we transfer solutions between meshes). We do so by
 * setting the yield stress to an unreasonably large value in these iterations
 * and then setting it back to the correct value in subsequent iterations.
 *   

 * 
 * Other than this, the top part of this function should be
 * reasonably obvious. We initialize the variable
 * <code>previous_residual_norm</code> to the most negative value
 * representable with double precision numbers so that the
 * comparison whether the current residual is less than that of the
 * previous step will always fail in the first step.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::solve_newton()
 *   {
 *     TrilinosWrappers::MPI::Vector old_solution(locally_owned_dofs,
 *                                                mpi_communicator);
 *     TrilinosWrappers::MPI::Vector residual(locally_owned_dofs,
 *                                            mpi_communicator);
 *     TrilinosWrappers::MPI::Vector tmp_vector(locally_owned_dofs,
 *                                              mpi_communicator);
 *     TrilinosWrappers::MPI::Vector locally_relevant_tmp_vector(
 *       locally_relevant_dofs, mpi_communicator);
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs,
 *                                                        mpi_communicator);
 * 
 *     double residual_norm;
 *     double previous_residual_norm = -std::numeric_limits<double>::max();
 * 
 *     const double correct_sigma = sigma_0;
 * 
 *     IndexSet old_active_set(active_set);
 * 
 *     for (unsigned int newton_step = 1; newton_step <= 100; ++newton_step)
 *       {
 *         if (newton_step == 1 &&
 *             ((transfer_solution && current_refinement_cycle == 0) ||
 *              !transfer_solution))
 *           constitutive_law.set_sigma_0(1e+10);
 *         else if (newton_step == 2 || current_refinement_cycle > 0 ||
 *                  !transfer_solution)
 *           constitutive_law.set_sigma_0(correct_sigma);
 * 
 *         pcout << " " << std::endl;
 *         pcout << "   Newton iteration " << newton_step << std::endl;
 *         pcout << "      Updating active set..." << std::endl;
 * 
 *         {
 *           TimerOutput::Scope t(computing_timer, "update active set");
 *           update_solution_and_constraints();
 *         }
 * 
 *         pcout << "      Assembling system... " << std::endl;
 *         newton_matrix = 0;
 *         newton_rhs    = 0;
 *         assemble_newton_system(solution);
 * 
 *         pcout << "      Solving system... " << std::endl;
 *         solve_newton_system();
 * 
 * @endcode
 * 
 * It gets a bit more hairy after we have computed the
 * trial solution $\tilde{\mathbf u}$ of the current Newton step.
 * We handle a highly nonlinear problem so we have to damp
 * Newton's method using a line search. To understand how we do this,
 * recall that in our formulation, we compute a trial solution
 * in each Newton step and not the update between old and new solution.
 * Since the solution set is a convex set, we will use a line
 * search that tries linear combinations of the
 * previous and the trial solution to guarantee that the
 * damped solution is in our solution set again.
 * At most we apply 5 damping steps.
 *         

 * 
 * There are exceptions to when we use a line search. First,
 * if this is the first Newton step on any mesh, then we don't have
 * any point to compare the residual to, so we always accept a full
 * step. Likewise, if this is the second Newton step on the first mesh
 * (or the second on any mesh if we don't transfer solutions from mesh
 * to mesh), then we have computed the first of these steps using just
 * an elastic model (see how we set the yield stress sigma to an
 * unreasonably large value above). In this case, the first Newton
 * solution was a purely elastic one, the second one a plastic one,
 * and any linear combination would not necessarily be expected to
 * lie in the feasible set -- so we just accept the solution we just
 * got.
 *         

 * 
 * In either of these two cases, we bypass the line search and just
 * update residual and other vectors as necessary.
 * 
 * @code
 *         if ((newton_step == 1) ||
 *             (transfer_solution && newton_step == 2 &&
 *              current_refinement_cycle == 0) ||
 *             (!transfer_solution && newton_step == 2))
 *           {
 *             compute_nonlinear_residual(solution);
 *             old_solution = solution;
 * 
 *             residual                     = newton_rhs;
 *             const unsigned int start_res = (residual.local_range().first),
 *                                end_res   = (residual.local_range().second);
 *             for (unsigned int n = start_res; n < end_res; ++n)
 *               if (all_constraints.is_inhomogeneously_constrained(n))
 *                 residual(n) = 0;
 * 
 *             residual.compress(VectorOperation::insert);
 * 
 *             residual_norm = residual.l2_norm();
 * 
 *             pcout << "      Accepting Newton solution with residual: "
 *                   << residual_norm << std::endl;
 *           }
 *         else
 *           {
 *             for (unsigned int i = 0; i < 5; ++i)
 *               {
 *                 distributed_solution = solution;
 * 
 *                 const double alpha = std::pow(0.5, static_cast<double>(i));
 *                 tmp_vector         = old_solution;
 *                 tmp_vector.sadd(1 - alpha, alpha, distributed_solution);
 * 
 *                 TimerOutput::Scope t(computing_timer, "Residual and lambda");
 * 
 *                 locally_relevant_tmp_vector = tmp_vector;
 *                 compute_nonlinear_residual(locally_relevant_tmp_vector);
 *                 residual = newton_rhs;
 * 
 *                 const unsigned int start_res = (residual.local_range().first),
 *                                    end_res   = (residual.local_range().second);
 *                 for (unsigned int n = start_res; n < end_res; ++n)
 *                   if (all_constraints.is_inhomogeneously_constrained(n))
 *                     residual(n) = 0;
 * 
 *                 residual.compress(VectorOperation::insert);
 * 
 *                 residual_norm = residual.l2_norm();
 * 
 *                 pcout
 *                   << "      Residual of the non-contact part of the system: "
 *                   << residual_norm << std::endl
 *                   << "         with a damping parameter alpha = " << alpha
 *                   << std::endl;
 * 
 *                 if (residual_norm < previous_residual_norm)
 *                   break;
 *               }
 * 
 *             solution     = tmp_vector;
 *             old_solution = solution;
 *           }
 * 
 *         previous_residual_norm = residual_norm;
 * 
 * 
 * @endcode
 * 
 * The final step is to check for convergence. If the active set
 * has not changed across all processors and the residual is
 * less than a threshold of $10^{-10}$, then we terminate
 * the iteration on the current mesh:
 * 
 * @code
 *         if (Utilities::MPI::sum((active_set == old_active_set) ? 0 : 1,
 *                                 mpi_communicator) == 0)
 *           {
 *             pcout << "      Active set did not change!" << std::endl;
 *             if (residual_norm < 1e-10)
 *               break;
 *           }
 * 
 *         old_active_set = active_set;
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemrefine_grid"></a> 
 * <h4>PlasticityContactProblem::refine_grid</h4>
 * 

 * 
 * If you've made it this far into the deal.II tutorial, the following
 * function refining the mesh should not pose any challenges to you
 * any more. It refines the mesh, either globally or using the Kelly
 * error estimator, and if so asked also transfers the solution from
 * the previous to the next mesh. In the latter case, we also need
 * to compute the active set and other quantities again, for which we
 * need the information computed by <code>compute_nonlinear_residual()</code>.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::refine_grid()
 *   {
 *     if (refinement_strategy == RefinementStrategy::refine_global)
 *       {
 *         for (typename Triangulation<dim>::active_cell_iterator cell =
 *                triangulation.begin_active();
 *              cell != triangulation.end();
 *              ++cell)
 *           if (cell->is_locally_owned())
 *             cell->set_refine_flag();
 *       }
 *     else
 *       {
 *         Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
 *         KellyErrorEstimator<dim>::estimate(
 *           dof_handler,
 *           QGauss<dim - 1>(fe.degree + 2),
 *           std::map<types::boundary_id, const Function<dim> *>(),
 *           solution,
 *           estimated_error_per_cell);
 * 
 *         parallel::distributed::GridRefinement ::refine_and_coarsen_fixed_number(
 *           triangulation, estimated_error_per_cell, 0.3, 0.03);
 *       }
 * 
 *     triangulation.prepare_coarsening_and_refinement();
 * 
 *     parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
 *       solution_transfer(dof_handler);
 *     if (transfer_solution)
 *       solution_transfer.prepare_for_coarsening_and_refinement(solution);
 * 
 *     triangulation.execute_coarsening_and_refinement();
 * 
 *     setup_system();
 * 
 *     if (transfer_solution)
 *       {
 *         TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs,
 *                                                            mpi_communicator);
 *         solution_transfer.interpolate(distributed_solution);
 * 
 * @endcode
 * 
 * enforce constraints to make the interpolated solution conforming on
 * the new mesh:
 * 
 * @code
 *         constraints_hanging_nodes.distribute(distributed_solution);
 * 
 *         solution = distributed_solution;
 *         compute_nonlinear_residual(solution);
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemmove_mesh"></a> 
 * <h4>PlasticityContactProblem::move_mesh</h4>
 * 

 * 
 * The remaining three functions before we get to <code>run()</code>
 * have to do with generating output. The following one is an attempt
 * at showing the deformed body in its deformed configuration. To this
 * end, this function takes a displacement vector field and moves every
 * vertex of the (local part) of the mesh by the previously computed
 * displacement. We will call this function with the current
 * displacement field before we generate graphical output, and we will
 * call it again after generating graphical output with the negative
 * displacement field to undo the changes to the mesh so made.
 *   

 * 
 * The function itself is pretty straightforward. All we have to do
 * is keep track which vertices we have already touched, as we
 * encounter the same vertices multiple times as we loop over cells.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::move_mesh(
 *     const TrilinosWrappers::MPI::Vector &displacement) const
 *   {
 *     std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         for (const auto v : cell->vertex_indices())
 *           if (vertex_touched[cell->vertex_index(v)] == false)
 *             {
 *               vertex_touched[cell->vertex_index(v)] = true;
 * 
 *               Point<dim> vertex_displacement;
 *               for (unsigned int d = 0; d < dim; ++d)
 *                 vertex_displacement[d] =
 *                   displacement(cell->vertex_dof_index(v, d));
 * 
 *               cell->vertex(v) += vertex_displacement;
 *             }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemoutput_results"></a> 
 * <h4>PlasticityContactProblem::output_results</h4>
 * 

 * 
 * Next is the function we use to actually generate graphical output. The
 * function is a bit tedious, but not actually particularly complicated.
 * It moves the mesh at the top (and moves it back at the end), then
 * computes the contact forces along the contact surface. We can do
 * so (as shown in the accompanying paper) by taking the untreated
 * residual vector and identifying which degrees of freedom
 * correspond to those with contact by asking whether they have an
 * inhomogeneous constraints associated with them. As always, we need
 * to be mindful that we can only write into completely distributed
 * vectors (i.e., vectors without ghost elements) but that when we
 * want to generate output, we need vectors that do indeed have
 * ghost entries for all locally relevant degrees of freedom.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::output_results(
 *     const unsigned int current_refinement_cycle)
 *   {
 *     TimerOutput::Scope t(computing_timer, "Graphical output");
 * 
 *     pcout << "      Writing graphical output... " << std::flush;
 * 
 *     move_mesh(solution);
 * 
 * @endcode
 * 
 * Calculation of the contact forces
 * 
 * @code
 *     TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs,
 *                                                      mpi_communicator);
 *     const unsigned int start_res = (newton_rhs_uncondensed.local_range().first),
 *                        end_res = (newton_rhs_uncondensed.local_range().second);
 *     for (unsigned int n = start_res; n < end_res; ++n)
 *       if (all_constraints.is_inhomogeneously_constrained(n))
 *         distributed_lambda(n) =
 *           newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n);
 *     distributed_lambda.compress(VectorOperation::insert);
 *     constraints_hanging_nodes.distribute(distributed_lambda);
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs,
 *                                          mpi_communicator);
 *     lambda = distributed_lambda;
 * 
 *     TrilinosWrappers::MPI::Vector distributed_active_set_vector(
 *       locally_owned_dofs, mpi_communicator);
 *     distributed_active_set_vector = 0.;
 *     for (const auto index : active_set)
 *       distributed_active_set_vector[index] = 1.;
 *     distributed_lambda.compress(VectorOperation::insert);
 * 
 *     TrilinosWrappers::MPI::Vector active_set_vector(locally_relevant_dofs,
 *                                                     mpi_communicator);
 *     active_set_vector = distributed_active_set_vector;
 * 
 *     DataOut<dim> data_out;
 * 
 *     data_out.attach_dof_handler(dof_handler);
 * 
 *     const std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_out.add_data_vector(solution,
 *                              std::vector<std::string>(dim, "displacement"),
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 *     data_out.add_data_vector(lambda,
 *                              std::vector<std::string>(dim, "contact_force"),
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 *     data_out.add_data_vector(active_set_vector,
 *                              std::vector<std::string>(dim, "active_set"),
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells());
 *     for (unsigned int i = 0; i < subdomain.size(); ++i)
 *       subdomain(i) = triangulation.locally_owned_subdomain();
 *     data_out.add_data_vector(subdomain, "subdomain");
 * 
 *     data_out.add_data_vector(fraction_of_plastic_q_points_per_cell,
 *                              "fraction_of_plastic_q_points");
 * 
 *     data_out.build_patches();
 * 
 * @endcode
 * 
 * In the remainder of the function, we generate one VTU file on
 * every processor, indexed by the subdomain id of this processor.
 * On the first processor, we then also create a <code>.pvtu</code>
 * file that indexes <i>all</i> of the VTU files so that the entire
 * set of output files can be read at once. These <code>.pvtu</code>
 * are used by Paraview to describe an entire parallel computation's
 * output files. We then do the same again for the competitor of
 * Paraview, the VisIt visualization program, by creating a matching
 * <code>.visit</code> file.
 * 
 * @code
 *     const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record(
 *       output_dir, "solution", current_refinement_cycle, mpi_communicator, 2);
 *     pcout << pvtu_filename << std::endl;
 * 
 *     TrilinosWrappers::MPI::Vector tmp(solution);
 *     tmp *= -1;
 *     move_mesh(tmp);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemoutput_contact_force"></a> 
 * <h4>PlasticityContactProblem::output_contact_force</h4>
 * 

 * 
 * This last auxiliary function computes the contact force by
 * calculating an integral over the contact pressure in z-direction
 * over the contact area. For this purpose we set the contact
 * pressure lambda to 0 for all inactive dofs (whether a degree
 * of freedom is part of the contact is determined just as
 * we did in the previous function). For all
 * active dofs, lambda contains the quotient of the nonlinear
 * residual (newton_rhs_uncondensed) and corresponding diagonal entry
 * of the mass matrix (diag_mass_matrix_vector). Because it is
 * not unlikely that hanging nodes show up in the contact area
 * it is important to apply constraints_hanging_nodes.distribute
 * to the distributed_lambda vector.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::output_contact_force() const
 *   {
 *     TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs,
 *                                                      mpi_communicator);
 *     const unsigned int start_res = (newton_rhs_uncondensed.local_range().first),
 *                        end_res = (newton_rhs_uncondensed.local_range().second);
 *     for (unsigned int n = start_res; n < end_res; ++n)
 *       if (all_constraints.is_inhomogeneously_constrained(n))
 *         distributed_lambda(n) =
 *           newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n);
 *       else
 *         distributed_lambda(n) = 0;
 *     distributed_lambda.compress(VectorOperation::insert);
 *     constraints_hanging_nodes.distribute(distributed_lambda);
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs,
 *                                          mpi_communicator);
 *     lambda = distributed_lambda;
 * 
 *     double contact_force = 0.0;
 * 
 *     QGauss<dim - 1>   face_quadrature_formula(fe.degree + 1);
 *     FEFaceValues<dim> fe_values_face(fe,
 *                                      face_quadrature_formula,
 *                                      update_values | update_JxW_values);
 * 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size();
 * 
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         for (const auto &face : cell->face_iterators())
 *           if (face->at_boundary() && face->boundary_id() == 1)
 *             {
 *               fe_values_face.reinit(cell, face);
 * 
 *               std::vector<Tensor<1, dim>> lambda_values(n_face_q_points);
 *               fe_values_face[displacement].get_function_values(lambda,
 *                                                                lambda_values);
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points;
 *                    ++q_point)
 *                 contact_force +=
 *                   lambda_values[q_point][2] * fe_values_face.JxW(q_point);
 *             }
 *     contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD);
 * 
 *     pcout << "Contact force = " << contact_force << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemrun"></a> 
 * <h4>PlasticityContactProblem::run</h4>
 * 

 * 
 * As in all other tutorial programs, the <code>run()</code> function contains
 * the overall logic. There is not very much to it here: in essence, it
 * performs the loops over all mesh refinement cycles, and within each, hands
 * things over to the Newton solver in <code>solve_newton()</code> on the
 * current mesh and calls the function that creates graphical output for
 * the so-computed solution. It then outputs some statistics concerning both
 * run times and memory consumption that has been collected over the course of
 * computations on this mesh.
 * 
 * @code
 *   template <int dim>
 *   void PlasticityContactProblem<dim>::run()
 *   {
 *     computing_timer.reset();
 *     for (; current_refinement_cycle < n_refinement_cycles;
 *          ++current_refinement_cycle)
 *       {
 *         {
 *           TimerOutput::Scope t(computing_timer, "Setup");
 * 
 *           pcout << std::endl;
 *           pcout << "Cycle " << current_refinement_cycle << ':' << std::endl;
 * 
 *           if (current_refinement_cycle == 0)
 *             {
 *               make_grid();
 *               setup_system();
 *             }
 *           else
 *             {
 *               TimerOutput::Scope t(computing_timer, "Setup: refine mesh");
 *               refine_grid();
 *             }
 *         }
 * 
 *         solve_newton();
 * 
 *         output_results(current_refinement_cycle);
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 * 
 *         Utilities::System::MemoryStats stats;
 *         Utilities::System::get_memory_stats(stats);
 *         pcout << "Peak virtual memory used, resident in kB: " << stats.VmSize
 *               << " " << stats.VmRSS << std::endl;
 * 
 *         if (base_mesh == "box")
 *           output_contact_force();
 *       }
 *   }
 * } // namespace Step42
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * There really isn't much to the <code>main()</code> function. It looks
 * like they always do:
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   using namespace dealii;
 *   using namespace Step42;
 * 
 *   try
 *     {
 *       ParameterHandler prm;
 *       PlasticityContactProblem<3>::declare_parameters(prm);
 *       if (argc != 2)
 *         {
 *           std::cerr << "*** Call this program as <./step-42 input.prm>"
 *                     << std::endl;
 *           return 1;
 *         }
 * 
 *       prm.parse_input(argv[1]);
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(
 *         argc, argv, numbers::invalid_unsigned_int);
 *       {
 *         PlasticityContactProblem<3> problem(prm);
 *         problem.run();
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
examples/step-42/doc/results.dox



<a name="Results"></a><h1>Results</h1>


包含这个程序的目录还包含一些输入参数文件，可以用来创建各种不同的模拟。例如，用 <code>p1_adaptive.prm</code> 参数文件（用球作为障碍物，用盒子作为领域）在16个核心上运行该程序会产生这样的输出。

@code
    Using output directory 'p1adaptive/'
    FE degree 1
    transfer solution false


Cycle 0:
   Number of active cells: 512
   Number of degrees of freedom: 2187


  Newton iteration 1
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 173.076 -> 1.64265e-06 in 7 Bicgstab iterations.
      Accepting Newton solution with residual: 1.64265e-06


   Newton iteration 2
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 57.3622 -> 3.23721e-07 in 8 Bicgstab iterations.
      Accepting Newton solution with residual: 24.9028
      Active set did not change!


   Newton iteration 3
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 24.9028 -> 9.94326e-08 in 7 Bicgstab iterations.
      Residual of the non-contact part of the system: 1.63333
         with a damping parameter alpha = 1
      Active set did not change!


...


  Newton iteration 6
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 1.43188e-07 -> 3.56218e-16 in 8 Bicgstab iterations.
      Residual of the non-contact part of the system: 4.298e-14
         with a damping parameter alpha = 1
      Active set did not change!
      Writing graphical output... p1_adaptive/solution-00.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      1.13s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assembling                      |         6 |     0.463s |        41% |
| Graphical output                |         1 |    0.0257s |       2.3% |
| Residual and lambda             |         4 |    0.0754s |       6.7% |
| Setup                           |         1 |     0.227s |        20% |
| Setup: constraints              |         1 |    0.0347s |       3.1% |
| Setup: distribute DoFs          |         1 |    0.0441s |       3.9% |
| Setup: matrix                   |         1 |    0.0119s |       1.1% |
| Setup: vectors                  |         1 |   0.00155s |      0.14% |
| Solve                           |         6 |     0.246s |        22% |
| Solve: iterate                  |         6 |    0.0631s |       5.6% |
| Solve: setup preconditioner     |         6 |     0.167s |        15% |
| update active set               |         6 |    0.0401s |       3.6% |
+---------------------------------+-----------+------------+------------+


Peak virtual memory used, resident in kB: 541884 77464
Contact force = 37.3058


...


Cycle 3:
   Number of active cells: 14652
   Number of degrees of freedom: 52497


   Newton iteration 1
      Updating active set...
         Size of active set: 145
      Assembling system...
      Solving system...
         Error: 296.309 -> 2.72484e-06 in 10 Bicgstab iterations.
      Accepting Newton solution with residual: 2.72484e-06


...


   Newton iteration 10
      Updating active set...
         Size of active set: 145
      Assembling system...
      Solving system...
         Error: 2.71541e-07 -> 1.5428e-15 in 27 Bicgstab iterations.
      Residual of the non-contact part of the system: 1.89261e-13
         with a damping parameter alpha = 1
      Active set did not change!
      Writing graphical output... p1_adaptive/solution-03.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      38.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assembling                      |        10 |      22.5s |        58% |
| Graphical output                |         1 |     0.327s |      0.85% |
| Residual and lambda             |         9 |      3.75s |       9.8% |
| Setup                           |         1 |      4.83s |        13% |
| Setup: constraints              |         1 |     0.578s |       1.5% |
| Setup: distribute DoFs          |         1 |      0.71s |       1.8% |
| Setup: matrix                   |         1 |     0.111s |      0.29% |
| Setup: refine mesh              |         1 |      4.83s |        13% |
| Setup: vectors                  |         1 |   0.00548s |     0.014% |
| Solve                           |        10 |      5.49s |        14% |
| Solve: iterate                  |        10 |       3.5s |       9.1% |
| Solve: setup preconditioner     |        10 |      1.84s |       4.8% |
| update active set               |        10 |     0.662s |       1.7% |
+---------------------------------+-----------+------------+------------+


Peak virtual memory used, resident in kB: 566052 105788
Contact force = 56.794


...
@endcode



每个周期结束时的表格显示了最近一次网格细化周期的计算时间（这些数字当然是针对产生该输出的机器而言的）和程序不同部分的调用次数，如装配或计算残差。上面的一些数字可以通过将解决方案从一个网格转移到下一个网格来改善，我们在这里没有行使这个选项。当然，你也可以通过使用更多的处理器来使程序运行得更快，特别是在后期的细化周期中：附带的论文显示，至少有1000个内核的良好扩展性。

在一个典型的运行中，你可以看到，对于每一个细化步骤，活动集--接触点--首先被迭代出来。之后，牛顿方法只需要解决塑性问题。对于更细的网格，在最后4或5次牛顿迭代中可以看到二次收敛。

我们不会在这里详细讨论每个输入文件的情况。相反，让我们只展示解决方案的图片（如果单元格的正交点为零，塑性不等式处于活动状态，则域的左半部分被省略）。

 <table align="center">
  <tr>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionColorbar.png">
    </td>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionBall2.png" alt="" width="70%">
    </td>
    <td valign="top">
      &nbsp;
    </td>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionLi2.png" alt="" alt="" width="70%">
    </td>
  </tr>
</table> 

图中显示了适应性细化以及细胞在与球接触过程中的塑化程度。请记住，我们考虑每个正交点的应力偏差部分的规范，以查看是否有弹性或塑性行为。蓝色意味着这个单元只包含弹性正交点，与所有正交点都被塑化的红色单元相反。在顶面的中间--网格最细的地方--非常仔细地看可以看到由障碍物引起的凹陷。这是 <code>move_mesh()</code> 函数的结果。然而，由于我们在这里考虑的障碍物的压痕非常小，所以很难辨别这种效果；我们可以玩玩将网格的顶点按计算出的位移的倍数进行位移。

关于使用该程序可以获得的结果的进一步讨论，见本页面最上方提到的出版物。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h1>Possibilities for extensions</h1>


像往常一样，有多种可能性来扩展这个程序。从算法的角度来看，这个程序在写作时已经达到了我们所能达到的程度，使用了接触不等式、塑性非线性和线性求解器的最佳可用算法。然而，就更现实的情况而言，人们希望用这个程序做一些事情。   <ul>   <li>  将程序从静态扩展到准静态情况，也许可以通过选择后向欧拉模式来实现时间离散化。一些理论结果可以在Jörg Frohne的博士论文中找到，<i>FEM-Simulation
der Umformtechnik metallischer Oberfl&auml;chen im Mikrokosmos</i>，德国锡根大学，2011。

 <li> 考虑有摩擦力的接触问题也将是一个有趣的进步。在几乎每个机械过程中，摩擦都有很大的影响。  为了模拟这种情况，我们必须考虑到接触面的切向应力。摩擦也给我们的问题增加了另一个不等式，因为只要切向应力不超过某个极限，身体和障碍物通常会粘在一起，超过这个极限，两个身体就会互相滑过。

 <li>  如果我们已经模拟了摩擦性接触，下一步要考虑的是接触区的发热。由两个物体之间的摩擦引起的热量会提高可变形物体的温度，并导致一些材料参数的变化。

 <li>  对于接触以及塑性，实施更精确的、与问题相适应的误差估计器可能是有意义的。   </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-42.cc"
*/
