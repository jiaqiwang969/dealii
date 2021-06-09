/**
@page step_44 The step-44 tutorial program
This tutorial depends on step-18.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Listofreferences">List of references</a>
        <li><a href="#Notation"> Notation </a>
        <li><a href="#Kinematics">Kinematics</a>
        <li><a href="#Kinetics">Kinetics</a>
        <li><a href="#Pushforwardandpullbackoperators"> Push-forward and pull-back operators </a>
        <li><a href="#Hyperelasticmaterials">Hyperelastic materials</a>
      <ul>
        <li><a href="#NeoHookeanmaterials"> Neo-Hookean materials </a>
      </ul>
        <li><a href="#Elasticitytensors">Elasticity tensors</a>
        <li><a href="#Principleofstationarypotentialenergyandthethreefieldformulation">Principle of stationary potential energy and the three-field formulation</a>
        <li><a href="#Discretizationofgoverningequations"> Discretization of governing equations </a>
        <li><a href="#Thematerialclass"> The material class </a>
        <li><a href="#Numericalexample"> Numerical example </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Runtimeparameters">Run-time parameters</a>
      <ul>
        <li><a href="#FiniteElementsystem">Finite Element system</a>
        <li><a href="#Geometry">Geometry</a>
        <li><a href="#Materials">Materials</a>
        <li><a href="#Linearsolver">Linear solver</a>
        <li><a href="#Nonlinearsolver">Nonlinear solver</a>
        <li><a href="#Time">Time</a>
        <li><a href="#Allparameters">All parameters</a>
      </ul>
        <li><a href="#Timeclass">Time class</a>
        <li><a href="#CompressibleneoHookeanmaterialwithinathreefieldformulation">Compressible neo-Hookean material within a three-field formulation</a>
        <li><a href="#Quadraturepointhistory">Quadrature point history</a>
        <li><a href="#Quasistaticquasiincompressiblefinitestrainsolid">Quasi-static quasi-incompressible finite-strain solid</a>
        <li><a href="#ImplementationofthecodeSolidcodeclass">Implementation of the <code>Solid</code> class</a>
      <ul>
        <li><a href="#Publicinterface">Public interface</a>
      </ul>
        <li><a href="#Privateinterface">Private interface</a>
      <ul>
        <li><a href="#Threadingbuildingblocksstructures">Threading-building-blocks structures</a>
        <li><a href="#Solidmake_grid">Solid::make_grid</a>
        <li><a href="#Solidsystem_setup">Solid::system_setup</a>
        <li><a href="#Soliddetermine_component_extractors">Solid::determine_component_extractors</a>
        <li><a href="#Solidsetup_qph">Solid::setup_qph</a>
        <li><a href="#Solidupdate_qph_incremental">Solid::update_qph_incremental</a>
        <li><a href="#Solidsolve_nonlinear_timestep">Solid::solve_nonlinear_timestep</a>
        <li><a href="#Solidprint_conv_headerandSolidprint_conv_footer">Solid::print_conv_header and Solid::print_conv_footer</a>
        <li><a href="#Solidget_error_dilation">Solid::get_error_dilation</a>
        <li><a href="#Solidget_error_residual">Solid::get_error_residual</a>
        <li><a href="#Solidget_error_update">Solid::get_error_update</a>
        <li><a href="#Solidget_total_solution">Solid::get_total_solution</a>
        <li><a href="#Solidassemble_system">Solid::assemble_system</a>
        <li><a href="#Solidmake_constraints">Solid::make_constraints</a>
        <li><a href="#Solidassemble_sc">Solid::assemble_sc</a>
        <li><a href="#Solidsolve_linear_system">Solid::solve_linear_system</a>
        <li><a href="#Solidoutput_results">Solid::output_results</a>
      </ul>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-44/doc/intro.dox

 <br> 

<i>This program was contributed by Jean-Paul Pelteret and Andrew McBride.
<br>
This material is based upon work supported by  the German Science Foundation (Deutsche
Forschungsgemeinschaft, DFG), grant STE 544/39-1,  and the National Research Foundation of South Africa.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.439772,https://zenodo.org/badge/DOI/10.5281/zenodo.439772.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程的主题是非线性固体力学。经典的单场方法（例如见步骤18）不能正确描述准不可压缩材料的响应。响应过于僵硬；这种现象被称为锁定。锁定问题可以通过各种替代策略来规避。其中一个策略是三场公式。在这里，它被用来模拟各向同性连续体的三维、完全非线性（几何和材料）响应。材料响应被近似为超弹性。此外，所采用的三场公式对准不可压缩和可压缩材料都有效。

本报告的目的是为使用deal.II处理非线性固体力学的问题提供基础。线性问题在步骤8中得到了解决。在第18步中部分考虑了几何非线性问题的非标准的、超弹性的形式：使用了线性化构成关系的速率形式，问题域随着运动的进行而变化。围绕非线性运动学的重要概念在理论和实施中都没有。然而，第18步确实描述了许多关键概念，以便在deal.II的框架内实现弹性。

我们从非线性运动学的速成课程开始。为了简单起见，我们将注意力限制在准静态问题上。此后，我们介绍了各种关键的应力测量，并描述了构成模型。然后，在解释用于管理材料的类的结构之前，我们详细描述了三场公式。然后介绍了例子问题的设置。

 @note  本教程是针对三维空间的弹性问题而开发的（并在介绍中进行了描述）。  虽然空间维度可以在main()例程中改变，但需要注意的是。  一般来说，二维弹性问题只是作为三维问题的理想化而存在。  也就是说，它们要么是平面应变，要么是平面应力。  这些选择中的任何一个的假设都需要被一致地施加。  更多信息请参见步骤8的说明。

<a name="Listofreferences"></a><h3>List of references</h3>


这里实施的三场公式是由Simo等人（1985）开创的，被称为混合雅各布-压力公式。重要的相关贡献包括Simo和Taylor（1991）以及Miehe（1994）的贡献。这里采用的符号在很大程度上借鉴了Holzapfel（2001）对非线性固体力学理论方面的出色概述。Hughes (2000)对与不可压缩弹性（小应变时）有关的问题作了很好的概述。

<ol>  <li>  J.C. Simo, R.L. Taylor and K.S. Pister (1985), Variational and projection methods for the volume constraint in finite deformation elasto-plasticity,  <em>  Computer Methods in Applied Mechanics and Engineering  </em>  , <strong> 51</strong>, 1-3, 177-208. 		DOI: <a href="http://doi.org/10.1016/0045-7825(85)90033-7">10.1016/0045-7825(85)90033-7</a>;  <li>  J.C. Simo and R.L. Taylor (1991), Quasi-incompressible finite elasticity in principal stretches.Continuum basis and numerical algorithms,  <em>  Computer Methods in Applied Mechanics and Engineering  </em>  , <strong> 85 </strong>, 3, 273-310. 		DOI: <a href="http://doi.org/10.1016/0045-7825(91)90100-K">10.1016/0045-7825(91)90100-K</a>;  <li>  C. Miehe (1994), Aspects of the formulation and finite element implementation of large strain isotropic elasticity  <em>  International Journal for Numerical Methods in Engineering  </em>  <strong> 37 /strong>, 12, 1981-2004. 		DOI: <a href="http://doi.org/10.1002/nme.1620371202">10.1002/nme.1620371202</a>;  <li>  G.A. Holzapfel (2001), Nonlinear Solid Mechanics.A Continuum Approach for Engineering, John Wiley & Sons. 		ISBN: 0-471-82304-X;  <li>  T.J.R. Hughes (2000), The Finite Element Method:线性静态和动态有限元分析》，多佛。 		ISBN: 978-0486411811  </ol> .

<ol>  <li>  J-P. V. Pelteret, D. Davydov, A. McBride, D. K. Vu, and P. Steinmann (2016), 在一个耦合问题中使用这种三场公式的例子记录在<ol>  <li>  J-P.V. Pelteret, D. Davydov, A. McBride, D. K. Vu, and P. Steinmann (2016), Computational electro-and magneto-elasticity for quasi-incompressible media immersed in free space,  <em>  International Journal for Numerical Methods in Engineering  </em>  。 		DOI: <a href="http://doi.org/10.1002/nme.5254">10.1002/nme.5254</a>  </ol>  。

<a name="Notation"></a><h3> Notation </h3>


我们可以把四阶张量看作是将二阶张量（矩阵）映射到自己身上的线性算子，其方式与矩阵将向量映射到向量上一样。有各种四阶单位张量，在即将到来的介绍中会用到。四阶单位张量 $\mathcal{I}$ 和 $\overline{\mathcal{I}}$ 定义如下

@f[
	\mathbf{A} = \mathcal{I}:\mathbf{A}
		\qquad \text{and} \qquad
	\mathbf{A}^T = \overline{\mathcal{I}}:\mathbf{A} \, .


@f]

注意  $\mathcal{I} \neq \overline{\mathcal{I}}^T$  。此外，我们通过以下方式定义对称和偏斜对称的四阶单位张量

@f[
	\mathcal{S} \dealcoloneq \dfrac{1}{2}[\mathcal{I} + \overline{\mathcal{I}}]
		\qquad \text{and} \qquad
	\mathcal{W} \dealcoloneq \dfrac{1}{2}[\mathcal{I} - \overline{\mathcal{I}}] \, ,


@f]

以致于

@f[
	\dfrac{1}{2}[\mathbf{A} + \mathbf{A}^T] = \mathcal{S}:\mathbf{A}
		\qquad \text{and} \qquad
	\dfrac{1}{2}[\mathbf{A} - \mathbf{A}^T] = \mathcal{W}:\mathbf{A} \, .


@f]

identity_tensor()返回的四阶  <code>SymmetricTensor</code>  是  $\mathcal{S}$  。




<a name="Kinematics"></a><h3>Kinematics</h3>


让时间域表示为 $\mathbb{T} = [0,T_{\textrm{end}}]$  ，其中 $t \in \mathbb{T}$ 和 $T_{\textrm{end}}$ 是总的问题持续时间。考虑一个连续体，在时间 $t=0$ 占据参考配置 $\Omega_0$ 。参考配置中的%粒子由位置矢量 $\mathbf{X}$ 识别。身体在后来的时间 $t>0$ 的配置被称为当前配置，表示为 $\Omega$ ，粒子由矢量 $\mathbf{x}$ 识别。参考配置和当前配置之间的非线性映射，表示为  $\boldsymbol{\varphi}$  ，作用如下。

@f[
	\mathbf{x} = \boldsymbol{\varphi}(\mathbf{X},t) \, .


@f]

粒子的位移的材料描述被定义为

@f[
	\mathbf{U}(\mathbf{X},t) = \mathbf{x}(\mathbf{X},t) - \mathbf{X} \, .


@f]



变形梯度 $\mathbf{F}$ 被定义为运动的材料梯度。

@f[
	\mathbf{F}(\mathbf{X},t)
		\dealcoloneq \dfrac{\partial \boldsymbol{\varphi}(\mathbf{X},t)}{\partial \mathbf{X}}
		= \textrm{Grad}\ \mathbf{x}(\mathbf{X},t)
		= \mathbf{I} + \textrm{Grad}\ \mathbf{U} \, .


@f]

变形梯度 $J(\mathbf{X},t) \dealcoloneq \textrm{det}\ \mathbf{F}(\mathbf{X},t) > 0$ 的行列式在参考配置和当前配置中映射出相应的体积元素，分别表示为 $\textrm{d}V$ 和 $\textrm{d}v$  ，为

@f[
	\textrm{d}v = J(\mathbf{X},t)\; \textrm{d}V \, .


@f]



就空间和材料坐标而言，变形的两个重要度量是左和右Cauchy-Green张量，分别表示为 $\mathbf{b} \dealcoloneq \mathbf{F}\mathbf{F}^T$ 和 $\mathbf{C} \dealcoloneq \mathbf{F}^T\mathbf{F}$  。它们都是对称的和正定的。

格林-拉格朗日应变张量的定义为

@f[
	\mathbf{E} \dealcoloneq \frac{1}{2}[\mathbf{C} - \mathbf{I} ]
		= \underbrace{\frac{1}{2}[\textrm{Grad}^T \mathbf{U} +	\textrm{Grad}\mathbf{U}]}_{\boldsymbol{\varepsilon}}
			+ \frac{1}{2}[\textrm{Grad}^T\ \mathbf{U}][\textrm{Grad}\ \mathbf{U}] \, .


@f]

如果假定变形为无限小，那么右边的第二项就可以忽略， $\boldsymbol{\varepsilon}$ （线性化的应变张量）是应变张量的唯一组成部分。从问题的设置来看，这个假设在步骤18中是不成立的，这使得在该教程程序中使用线性化的 $\boldsymbol{\varepsilon}$ 作为应变度量值得怀疑。

为了处理材料在受到体积和剪切型变形时表现出的不同响应，我们考虑将变形梯度 $\mathbf{F}$ 和左Cauchy-Green张量 $\mathbf{b}$ 分解为体积变化（体积）和体积保持（等效）部分。

@f[
	\mathbf{F}
		= (J^{1/3}\mathbf{I})\overline{\mathbf{F}}
	\qquad \text{and} \qquad
	\mathbf{b}
        = (J^{2/3}\mathbf{I})\overline{\mathbf{F}}\,\overline{\mathbf{F}}^T
		=  (J^{2/3}\mathbf{I})\overline{\mathbf{b}} \, .


@f]

显然， $\textrm{det}\ \mathbf{F} = \textrm{det}\ (J^{1/3}\mathbf{I}) = J$  。

空间速度场被表示为 $\mathbf{v}(\mathbf{x},t)$  。空间速度场相对于空间坐标的导数给出空间速度梯度  $\mathbf{l}(\mathbf{x},t)$  ，即

@f[
	\mathbf{l}(\mathbf{x},t)
		\dealcoloneq \dfrac{\partial \mathbf{v}(\mathbf{x},t)}{\partial \mathbf{x}}
		= \textrm{grad}\ \mathbf{v}(\mathbf{x},t) \, ,


@f]

其中 $\textrm{grad} \{\bullet \}
= \frac{\partial \{ \bullet \} }{ \partial \mathbf{x}}
= \frac{\partial \{ \bullet \} }{ \partial \mathbf{X}}\frac{\partial \mathbf{X} }{ \partial \mathbf{x}}
= \textrm{Grad} \{ \bullet \} \mathbf{F}^{-1}$  。




<a name="Kinetics"></a><h3>Kinetics</h3>


考奇应力定理将作用在当前构型 $\mathbf{t}$ 的无穷小表面元素上的考奇牵引力 $\mathrm{d}a$ 等同于考奇应力张量 $\boldsymbol{\sigma}$ （一个空间量）与表面的外向单位法线 $\mathbf{n}$ 的积，即

@f[
	\mathbf{t}(\mathbf{x},t, \mathbf{n}) = \boldsymbol{\sigma}\mathbf{n} \, .


@f]

Cauchy应力是对称的。同样，作用于参考构型 $\mathbf{T}$ 中的无穷小表面元素的第一皮奥拉-基尔霍夫牵引力 $\mathrm{d}A$ 是第一皮奥拉-基尔霍夫应力张量 $\mathbf{P}$ （两点张量）与表面的外向单位法线 $\mathbf{N}$ 的乘积，为

@f[
	\mathbf{T}(\mathbf{X},t, \mathbf{N}) = \mathbf{P}\mathbf{N} \, .


@f]

Cauchy牵引力 $\mathbf{t}$ 和第一个Piola-Kirchhoff牵引力 $\mathbf{T}$ 的关系为

@f[
	\mathbf{t}\mathrm{d}a = \mathbf{T}\mathrm{d}A \, .


@f]

这可以用<a href="http://en.wikipedia.org/wiki/Finite_strain_theory">Nanson's formula</a>来证明。

第一个Piola-Kirchhoff应力张量与Cauchy应力的关系为

@f[
	\mathbf{P} = J \boldsymbol{\sigma}\mathbf{F}^{-T} \, .


@f]

进一步的重要应力测量是（空间）基尔霍夫应力  $\boldsymbol{\tau} = J \boldsymbol{\sigma}$  和（参考）第二Piola-Kirchhoff应力  $\mathbf{S} = {\mathbf{F}}^{-1} \boldsymbol{\tau} {\mathbf{F}}^{-T}$  。




<a name="Pushforwardandpullbackoperators"></a><h3> Push-forward and pull-back operators </h3>


前推和后拉运算符允许人们在材料和空间设置之间转换各种措施。这里使用的应力测量是逆变的，而应变测量是协变的。

二阶协变张量 $(\bullet)^{\text{cov}}$ 的前推和后拉操作分别由以下方法给出。

@f[
	\chi_{*}(\bullet)^{\text{cov}} \dealcoloneq \mathbf{F}^{-T} (\bullet)^{\text{cov}} \mathbf{F}^{-1}
	\qquad \text{and} \qquad
	\chi^{-1}_{*}(\bullet)^{\text{cov}} \dealcoloneq \mathbf{F}^{T} (\bullet)^{\text{cov}} \mathbf{F} \, .


@f]



二阶禁忌张量 $(\bullet)^{\text{con}}$ 的前推和后拉操作分别由以下方法给出。

@f[
	\chi_{*}(\bullet)^{\text{con}} \dealcoloneq \mathbf{F} (\bullet)^{\text{con}} \mathbf{F}^T
	\qquad \text{and} \qquad
	\chi^{-1}_{*}(\bullet)^{\text{con}} \dealcoloneq \mathbf{F}^{-1} (\bullet)^{\text{con}} \mathbf{F}^{-T} \, .


@f]

例如  $\boldsymbol{\tau} = \chi_{*}(\mathbf{S})$  。




<a name="Hyperelasticmaterials"></a><h3>Hyperelastic materials</h3>


超弹性材料的响应受亥姆霍兹自由能函数 $\Psi = \Psi(\mathbf{F}) = \Psi(\mathbf{C}) = \Psi(\mathbf{b})$ 的制约，该函数作为应力的势能。例如，如果Helmholtz自由能取决于右Cauchy-Green张量 $\mathbf{C}$ ，那么各向同性的超弹性响应为

@f[
	\mathbf{S}
		= 2 \dfrac{\partial \Psi(\mathbf{C})}{\partial \mathbf{C}} \, .


@f]

如果亥姆霍兹自由能取决于左Cauchy-Green张量 $\mathbf{b}$ ，那么各向同性的超弹性响应为

@f[
	\boldsymbol{\tau}
		= 2 \dfrac{\partial \Psi(\mathbf{b})}{\partial \mathbf{b}} \mathbf{b}
		=  2 \mathbf{b} \dfrac{\partial \Psi(\mathbf{b})}{\partial \mathbf{b}} \, .


@f]



根据变形梯度的乘法分解，亥姆霍兹自由能可以分解为

@f[
	\Psi(\mathbf{b}) = \Psi_{\text{vol}}(J) + \Psi_{\text{iso}}(\overline{\mathbf{b}}) \, .


@f]

同样，基尔霍夫应力可以分解为体积部分和等效部分 $\boldsymbol{\tau} = \boldsymbol{\tau}_{\text{vol}} + \boldsymbol{\tau}_{\text{iso}}$ ，其中。

@f{align*}
	\boldsymbol{\tau}_{\text{vol}} &=
		2 \mathbf{b} \dfrac{\partial \Psi_{\textrm{vol}}(J)}{\partial \mathbf{b}}
		\\
		&= p J\mathbf{I} \, ,
		\\
	\boldsymbol{\tau}_{\text{iso}} &=
		2 \mathbf{b} \dfrac{\partial \Psi_{\textrm{iso}} (\overline{\mathbf{b}})}{\partial \mathbf{b}}
		\\
		&= \underbrace{( \mathcal{I} - \dfrac{1}{3} \mathbf{I} \otimes \mathbf{I})}_{\mathbb{P}} : \overline{\boldsymbol{\tau}} \, ,


@f}

其中 $p \dealcoloneq \dfrac{\partial \Psi_{\text{vol}}(J)}{\partial J}$ 是压力响应。   $\mathbb{P}$ 是投影张量，它提供了欧拉环境下的偏差算子。虚构的基尔霍夫应力张量 $\overline{\boldsymbol{\tau}}$ 被定义为

@f[
	\overline{\boldsymbol{\tau}}
		\dealcoloneq 2 \overline{\mathbf{b}} \dfrac{\partial \Psi_{\textrm{iso}}(\overline{\mathbf{b}})}{\partial \overline{\mathbf{b}}} \, .


@f]






 @note  上述定义的压力响应与固体力学中广泛使用的压力定义不同，即 $p = - 1/3 \textrm{tr} \boldsymbol{\sigma} = - 1/3 J^{-1} \textrm{tr} \boldsymbol{\tau}$  。这里 $p$ 是静水压力。我们在本教程中使用压力响应（尽管我们把它称为压力）。

<a name="NeoHookeanmaterials"></a><h4> Neo-Hookean materials </h4>


与可压缩<a href="http://en.wikipedia.org/wiki/Neo-Hookean_solid">neo-Hookean material</a>相对应的亥姆霍兹自由能由以下公式给出

@f[
    \Psi \equiv
        \underbrace{\kappa [ \mathcal{G}(J) ] }_{\Psi_{\textrm{vol}}(J)}
        + \underbrace{\bigl[c_1 [ \overline{I}_1 - 3] \bigr]}_{\Psi_{\text{iso}}(\overline{\mathbf{b}})} \, ,


@f]

其中 $\kappa \dealcoloneq \lambda + 2/3 \mu$ 是体积模量（ $\lambda$ 和 $\mu$ 是Lam&eacute; 参数）和 $\overline{I}_1 \dealcoloneq \textrm{tr}\ \overline{\mathbf{b}}$  。函数 $\mathcal{G}(J)$ 被要求是严格凸的，并满足 $\mathcal{G}(1) = 0$ 等条件，进一步的细节见Holzapfel（2001）。在这项工作中  $\mathcal{G} \dealcoloneq \frac{1}{4} [ J^2 - 1 - 2\textrm{ln}J ]$  .

不可压缩性对所有运动施加了等效约束  $J=1$  。对应于不可压缩的新胡克材料的亥姆霍兹自由能由以下公式给出

@f[
    \Psi \equiv
        \underbrace{\bigl[ c_1 [ I_1 - 3] \bigr] }_{\Psi_{\textrm{iso}}(\mathbf{b})} \, ,


@f]

其中  $ I_1 \dealcoloneq \textrm{tr}\mathbf{b} $  。因此，通过从可压缩自由能中去除体积分量并执行  $J=1$  得到不可压缩响应。




<a name="Elasticitytensors"></a><h3>Elasticity tensors</h3>


我们将使用Newton-Raphson策略来解决非线性边界值问题。因此，我们将需要将构成关系线性化。

材料描述中的四阶弹性张量定义为

@f[
	\mathfrak{C}
		= 2\dfrac{\partial \mathbf{S}(\mathbf{C})}{\partial \mathbf{C}}
		= 4\dfrac{\partial^2 \Psi(\mathbf{C})}{\partial \mathbf{C} \partial \mathbf{C}} \, .


@f]

空间描述 $\mathfrak{c}$ 中的四阶弹性张量由 $\mathfrak{C}$ 的推演得到，为

@f[
	\mathfrak{c} = J^{-1} \chi_{*}(\mathfrak{C})
		\qquad \text{and thus} \qquad
	J\mathfrak{c} = 4 \mathbf{b} \dfrac{\partial^2 \Psi(\mathbf{b})} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}	\, .


@f]

四阶弹性张量（对于超弹性材料）同时拥有主要和次要的对称性。

四阶空间弹性张量可以写成以下解耦形式。

@f[
	\mathfrak{c} = \mathfrak{c}_{\text{vol}} + \mathfrak{c}_{\text{iso}} \, ,


@f]

其中

@f{align*}
	J \mathfrak{c}_{\text{vol}}
		&= 4 \mathbf{b} \dfrac{\partial^2 \Psi_{\text{vol}}(J)} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}
		\\
		&= J[\widehat{p}\, \mathbf{I} \otimes \mathbf{I} - 2p \mathcal{I}]
			\qquad \text{where} \qquad
		\widehat{p} \dealcoloneq p + \dfrac{\textrm{d} p}{\textrm{d}J} \, ,
		\\
	J \mathfrak{c}_{\text{iso}}
		&=  4 \mathbf{b} \dfrac{\partial^2 \Psi_{\text{iso}}(\overline{\mathbf{b}})} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}
		\\
		&= \mathbb{P} : \mathfrak{\overline{c}} : \mathbb{P}
			+ \dfrac{2}{3}[\overline{\boldsymbol{\tau}}:\mathbf{I}]\mathbb{P}


			- \dfrac{2}{3}[ \mathbf{I}\otimes\boldsymbol{\tau}_{\text{iso}}
				+ \boldsymbol{\tau}_{\text{iso}} \otimes \mathbf{I} ] \, ,


@f}

其中空间描述中的虚构弹性张量 $\overline{\mathfrak{c}}$ 被定义为

@f[
	\overline{\mathfrak{c}}
		= 4 \overline{\mathbf{b}} \dfrac{ \partial^2 \Psi_{\textrm{iso}}(\overline{\mathbf{b}})} {\partial \overline{\mathbf{b}} \partial \overline{\mathbf{b}}} \overline{\mathbf{b}} \, .


@f]



<a name="Principleofstationarypotentialenergyandthethreefieldformulation"></a><h3>Principle of stationary potential energy and the three-field formulation</h3>


系统的总势能 $\Pi$ 是内部和外部势能之和，分别表示为 $\Pi_{\textrm{int}}$ 和 $\Pi_{\textrm{ext}}$  。我们希望通过最小化势能找到平衡配置。

如上所述，我们采用了三场的表述。我们用 $\mathbf{\Xi} \dealcoloneq \{ \mathbf{u}, \widetilde{p}, \widetilde{J} \}$ 表示主要未知数的集合。独立运动学变量 $\widetilde{J}$ 作为对 $J$ 的约束进入公式，由拉格朗日乘数 $\widetilde{p}$ （压力，我们将看到）强制执行。

这里使用的三场变分原理由以下公式给出

@f[
	\Pi(\mathbf{\Xi}) \dealcoloneq \int_\Omega \bigl[
		\Psi_{\textrm{vol}}(\widetilde{J})
		+ \widetilde{p}\,[J(\mathbf{u}) - \widetilde{J}]
		+ \Psi_{\textrm{iso}}(\overline{\mathbf{b}}(\mathbf{u}))
		\bigr] \textrm{d}v
	+ 	\Pi_{\textrm{ext}} \, ,


@f]

其中外部电势的定义为

@f[
	\Pi_{\textrm{ext}}
		= - \int_\Omega \mathbf{b}^\text{p} \cdot \mathbf{u}~\textrm{d}v


			- \int_{\partial \Omega_{\sigma}} \mathbf{t}^\text{p} \cdot \mathbf{u}~\textrm{d}a \, .


@f]

当前配置 $\partial \Omega$ 的边界由两部分组成： $\partial \Omega = \partial \Omega_{\mathbf{u}} \cup \partial \Omega_{\sigma}$  ，其中 $\partial \Omega_{\mathbf{u}} \cap \partial \Omega_{\boldsymbol{\sigma}} = \emptyset$  。规定的Cauchy牵引力，表示为  $\mathbf{t}^\text{p}$  ，被应用于  $ \partial \Omega_{\boldsymbol{\sigma}}$  ，而运动被规定在边界的其余部分  $\partial \Omega_{\mathbf{u}}$  。每单位电流体积的体力表示为  $\mathbf{b}^\text{p}$  。




势的静止性如下

@f{align*}
	R(\mathbf\Xi;\delta \mathbf{\Xi})
		&= D_{\delta \mathbf{\Xi}}\Pi(\mathbf{\Xi})
		\\
		&= \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \mathbf{u}} \cdot \delta \mathbf{u}
			+ \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \widetilde{p}} \delta \widetilde{p}
			+ \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \widetilde{J}} \delta \tilde{J}
			\\
		&= \int_{\Omega_0}  \left[
			\textrm{grad}\ \delta\mathbf{u} : [ \underbrace{[\widetilde{p} J \mathbf{I}]}_{\equiv \boldsymbol{\tau}_{\textrm{vol}}}
            +  \boldsymbol{\tau}_{\textrm{iso}}]
			+ \delta \widetilde{p}\, [ J(\mathbf{u}) - \widetilde{J}]
			+ \delta \widetilde{J}\left[ \dfrac{\textrm{d} \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}}


            -\widetilde{p}\right]
			\right]~\textrm{d}V
			\\
		&\quad - \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\textrm{d}V


			- \int_{\partial \Omega_{0,\boldsymbol{\sigma}}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\textrm{d}A
			\\
		&=0 \, ,


@f}

对于所有虚拟位移 $\delta \mathbf{u} \in H^1(\Omega)$ ，受 $\delta \mathbf{u} = \mathbf{0}$ 对 $\partial \Omega_{\mathbf{u}}$ 的约束，以及所有虚拟压力 $\delta \widetilde{p} \in L^2(\Omega)$ 和虚拟膨胀 $\delta \widetilde{J} \in L^2(\Omega)$ 。

人们应该注意到，在三个场的表述中 $\boldsymbol{\tau}_{\textrm{vol}} \equiv \widetilde{p} J \mathbf{I}$ ，体积基尔霍夫应力的定义和随后的体积正切与超弹性材料一节中给出的一般形式略有不同，其中 $\boldsymbol{\tau}_{\textrm{vol}} \equiv p J\mathbf{I}$ 。这是因为压力 $\widetilde{p}$ 现在是一个主要的场，而不是一个构成性的派生量。我们需要仔细区分主要场和从构成关系中得到的场。

 @note  虽然变量都是用空间量来表示的，但积分的领域是初始配置。这种方法被称为  <em>  总拉格朗日公式  </em>  。在步骤18中给出的方法，其积分域是当前配置，可以称为  <em>  更新的拉格朗日公式  </em>  。这两种方法的各种优点在文献中被广泛讨论。然而，应该指出的是，它们是等同的。


与残留物相对应的欧拉-拉格朗日方程为：。

@f{align*}
	&\textrm{div}\ \boldsymbol{\sigma} + \mathbf{b}^\text{p} = \mathbf{0} && \textrm{[equilibrium]}
		\\
	&J(\mathbf{u}) = \widetilde{J} 		&& \textrm{[dilatation]}
		\\
	&\widetilde{p} = \dfrac{\textrm{d} \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}} && \textrm{[pressure]} \, .


@f}

第一个方程是空间设置中的（准静态）平衡方程。第二个是约束条件  $J(\mathbf{u}) = \widetilde{J}$  。第三个是压力的定义  $\widetilde{p}$  。

 @note 下面的简化单场推导（ $\mathbf{u}$ 是唯一的主变量）使我们清楚地知道如何将积分的极限转化为参考域。

@f{align*}
\int_{\Omega}\delta \mathbf{u} \cdot [\textrm{div}\ \boldsymbol{\sigma} + \mathbf{b}^\text{p}]~\mathrm{d}v
&=
\int_{\Omega} [-\mathrm{grad}\delta \mathbf{u}:\boldsymbol{\sigma} + \delta \mathbf{u} \cdot\mathbf{b}^\text{p}]~\mathrm{d}v
  + \int_{\partial \Omega} \delta \mathbf{u} \cdot \mathbf{t}^\text{p}~\mathrm{d}a \\
&=


- \int_{\Omega_0} \mathrm{grad}\delta \mathbf{u}:\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot J\mathbf{b}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_0} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \\
&=


- \int_{\Omega_0} \mathrm{grad}\delta \mathbf{u}:\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_{0,\sigma}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \\
&=


- \int_{\Omega_0} [\mathrm{grad}\delta\mathbf{u}]^{\text{sym}} :\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_{0,\sigma}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \, ,


@f}

其中 $[\mathrm{grad}\delta\mathbf{u}]^{\text{sym}} = 1/2[ \mathrm{grad}\delta\mathbf{u} + [\mathrm{grad}\delta\mathbf{u}]^T] $  。

我们将使用迭代牛顿-拉弗森方法来解决非线性剩余方程  $R$  。为了简单起见，我们假设死荷载，即荷载不因变形而改变。

在  $t_{\textrm{n}-1}$  的已知状态和  $t_{\textrm{n}}$  的当前未知状态之间的数量变化被表示为  $\varDelta \{ \bullet \} = { \{ \bullet \} }^{\textrm{n}} - { \{ \bullet \} }^{\textrm{n-1}}$  。在当前迭代 $\textrm{i}$ 的数量值表示为  ${ \{ \bullet \} }^{\textrm{n}}_{\textrm{i}} = { \{ \bullet \} }_{\textrm{i}}$  。迭代  $\textrm{i}$  和  $\textrm{i}+1$  之间的增量变化被表示为  $d \{ \bullet \} \dealcoloneq \{ \bullet \}_{\textrm{i}+1} - \{ \bullet \}_{\textrm{i}}$  。

假设系统的状态在某个迭代中是已知的  $\textrm{i}$  。用牛顿-拉弗森方法求解的非线性治理方程的线性化近似值是：找到  $d \mathbf{\Xi}$  ，以便

@f[
	R(\mathbf{\Xi}_{\mathsf{i}+1}) =
		R(\mathbf{\Xi}_{\mathsf{i}})
		+ D^2_{d \mathbf{\Xi}, \delta \mathbf{\Xi}} \Pi(\mathbf{\Xi_{\mathsf{i}}}) \cdot d \mathbf{\Xi} \equiv 0 \, ,


@f]

然后设置  $\mathbf{\Xi}_{\textrm{i}+1} = \mathbf{\Xi}_{\textrm{i}}
+ d \mathbf{\Xi}$  。切线由以下公式给出

@f[
	D^2_{d \mathbf{\Xi}, \delta \mathbf{\Xi}} \Pi( \mathbf{\Xi}_{\mathsf{i}} )
		= D_{d \mathbf{\Xi}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})
		=: K(\mathbf{\Xi}_{\mathsf{i}}; d \mathbf{\Xi}, \delta \mathbf{\Xi}) \, .


@f]

因此。

@f{align*}
 	K(\mathbf{\Xi}_{\mathsf{i}}; d \mathbf{\Xi}, \delta \mathbf{\Xi})
 		&=
 			D_{d \mathbf{u}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi}) \cdot d \mathbf{u}
 			\\
 				&\quad +
 			 	D_{d \widetilde{p}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})  d \widetilde{p}
 			 \\
 			 	&\quad +
 			  D_{d \widetilde{J}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})  d \widetilde{J} \, ,


@f}

其中

@f{align*}
	D_{d \mathbf{u}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=
 	\int_{\Omega_0} \bigl[ \textrm{grad}\ \delta \mathbf{u} :
 			\textrm{grad}\ d \mathbf{u} [\boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}]
 			+ \textrm{grad}\ \delta \mathbf{u} :[
             \underbrace{[\widetilde{p}J[\mathbf{I}\otimes\mathbf{I} - 2 \mathcal{I}]}_{\equiv J\mathfrak{c}_{\textrm{vol}}} +
             J\mathfrak{c}_{\textrm{iso}}] :\textrm{grad} d \mathbf{u}
 		\bigr]~\textrm{d}V \, ,
 		\\
 	&\quad + \int_{\Omega_0} \delta \widetilde{p} J \mathbf{I} : \textrm{grad}\ d \mathbf{u} ~\textrm{d}V
 	\\
 	D_{d \widetilde{p}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=
 	\int_{\Omega_0} \textrm{grad}\ \delta \mathbf{u} : J \mathbf{I} d \widetilde{p} ~\textrm{d}V


 		-  \int_{\Omega_0} \delta \widetilde{J} d \widetilde{p}  ~\textrm{d}V \, ,
 	\\
 	D_{d \widetilde{J}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=  -\int_{\Omega_0} \delta \widetilde{p} d \widetilde{J}~\textrm{d}V
 	 + \int_{\Omega_0} \delta \widetilde{J}  \dfrac{\textrm{d}^2 \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}\textrm{d}\widetilde{J}} d \widetilde{J} ~\textrm{d}V \, .


@f}



注意，以下条款被称为几何应力和材料对切线矩阵的贡献。

@f{align*}
& \int_{\Omega_0} \textrm{grad}\ \delta \mathbf{u} :
 			\textrm{grad}\ d \mathbf{u} [\boldsymbol{\tau}_{\textrm{iso}} +  \boldsymbol{\tau}_{\textrm{vol}}]~\textrm{d}V
 			&& \quad {[\textrm{Geometrical stress}]} \, ,
 		\\
& \int_{\Omega_0} \textrm{grad} \delta \mathbf{u} :
 			[J\mathfrak{c}_{\textrm{vol}} + J\mathfrak{c}_{\textrm{iso}}] :\textrm{grad}\ d \mathbf{u}
 		~\textrm{d}V
 		&& \quad {[\textrm{Material}]} \, .


@f}






<a name="Discretizationofgoverningequations"></a><h3> Discretization of governing equations </h3>


这里使用的三场公式对准不可压缩材料是有效的，即在 $\nu \rightarrow 0.5$ （其中 $\nu$ 是<a
href="http://en.wikipedia.org/wiki/Poisson's_ratio">Poisson's ratio</a>）的地方，要很好地选择 $\mathbf{u},~\widetilde{p}$ 和 $\widetilde{J}$ 的插值场。通常情况下，选择 $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$ 。这里 $DGPM$ 是FE_DGPMonomial类。一个流行的选择是 $Q_1 \times DGPM_0 \times DGPM_0$ ，它被称为平均扩张法（见Hughes（2000）的直观讨论）。这个代码可以容纳 $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$ 的表述。不连续的近似允许 $\widetilde{p}$ 和 $\widetilde{J}$ 被浓缩出来，并恢复了基于位移的经典方法。

对于完全不可压缩的材料 $\nu = 0.5$ 和三场公式仍将表现出锁定行为。这可以通过在自由能中引入一个额外的约束条件来克服，其形式为  $\int_{\Omega_0} \Lambda [ \widetilde{J} - 1]~\textrm{d}V$  。这里 $\Lambda$ 是一个拉格朗日乘数，用于强制执行等时约束条件。进一步的细节见Miehe (1994)。

线性化的问题可以写成

@f[
	\mathbf{\mathsf{K}}( \mathbf{\Xi}_{\textrm{i}}) d\mathbf{\Xi}
	=
	\mathbf{ \mathsf{F}}(\mathbf{\Xi}_{\textrm{i}})


@f]

其中

@f{align*}
		\underbrace{\begin{bmatrix}
			\mathbf{\mathsf{K}}_{uu}	&	\mathbf{\mathsf{K}}_{u\widetilde{p}}	& \mathbf{0}
			\\
			\mathbf{\mathsf{K}}_{\widetilde{p}u}	&	\mathbf{0}	&	\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}
			\\
			\mathbf{0}	& 	\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}		& \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\end{bmatrix}}_{\mathbf{\mathsf{K}}(\mathbf{\Xi}_{\textrm{i}})}
		\underbrace{\begin{bmatrix}
			d \mathbf{\mathsf{u}}\\
            d \widetilde{\mathbf{\mathsf{p}}} \\
            d \widetilde{\mathbf{\mathsf{J}}}
		\end{bmatrix}}_{d \mathbf{\Xi}}
        =
        \underbrace{\begin{bmatrix}


			-\mathbf{\mathsf{R}}_{u}(\mathbf{u}_{\textrm{i}}) \\


            -\mathbf{\mathsf{R}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}}) \\


           -\mathbf{\mathsf{R}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
		\end{bmatrix}}_{ -\mathbf{\mathsf{R}}(\mathbf{\Xi}_{\textrm{i}}) }
=
        \underbrace{\begin{bmatrix}
			\mathbf{\mathsf{F}}_{u}(\mathbf{u}_{\textrm{i}}) \\
            \mathbf{\mathsf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}}) \\
           \mathbf{\mathsf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
		\end{bmatrix}}_{ \mathbf{\mathsf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .


@f}



在配方中没有压力和膨胀（主要）变量的导数存在。因此，压力和膨胀的不连续有限元插值产生了 $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$ 、 $\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}$ 和 $\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}$ 的块对角矩阵。因此，我们可以很容易地表达每个单元上的场 $\widetilde{p}$ 和 $\widetilde{J}$ ，只需倒置一个局部矩阵并乘以局部右手。然后我们可以将结果插入其余的方程中，并恢复一个经典的基于位移的方法。为了在元素水平上凝结出压力和膨胀的贡献，我们需要以下结果。

@f{align*}
		d \widetilde{\mathbf{\mathsf{p}}}
		& = \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \bigl[
			 \mathbf{\mathsf{F}}_{\widetilde{J}}


			 - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}} d \widetilde{\mathbf{\mathsf{J}}} \bigr]
			\\
		d \widetilde{\mathbf{\mathsf{J}}}
		& = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
			\mathbf{\mathsf{F}}_{\widetilde{p}}


			- \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
			\bigr]
		\\
		 \Rightarrow d \widetilde{\mathbf{\mathsf{p}}}
		&=  \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}


		- \underbrace{\bigl[\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathbf{\mathsf{K}}}}\bigl[ \mathbf{\mathsf{F}}_{\widetilde{p}}


 		- \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}} \bigr]


@f}

因此

@f[
		\underbrace{\bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
		}_{\mathbf{\mathsf{K}}_{\textrm{con}}} d \mathbf{\mathsf{u}}
		=
        \underbrace{
		\Bigl[
		\mathbf{\mathsf{F}}_{u}


			- \mathbf{\mathsf{K}}_{u\widetilde{p}} \bigl[ \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}


			- \overline{\mathbf{\mathsf{K}}}\mathbf{\mathsf{F}}_{\widetilde{p}} \bigr]
		\Bigr]}_{\mathbf{\mathsf{F}}_{\textrm{con}}}


@f]

其中

@f[
		\overline{\overline{\mathbf{\mathsf{K}}}} \dealcoloneq
			\mathbf{\mathsf{K}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{K}}} \mathbf{\mathsf{K}}_{\widetilde{p}u} \, .


@f]

请注意，由于 $\widetilde{p}$ 和 $\widetilde{J}$ 选择的是元素层面的不连续，所有需要反转的矩阵都是在元素层面定义的。

构建各种贡献的程序如下。

- 构建  $\mathbf{\mathsf{K}}$  。

- 形成  $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}$  的元素，并存储在  $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$  中的  $\mathbf{\mathsf{K}}$  。

- 形成 $\overline{\overline{\mathbf{\mathsf{K}}}}$ 并添加到 $\mathbf{\mathsf{K}}_{uu}$ ，得到 $\mathbf{\mathsf{K}}_{\textrm{con}}$ 。

- 修改后的系统矩阵被称为  ${\mathbf{\mathsf{K}}}_{\textrm{store}}$  。   也就是@f[
        \mathbf{\mathsf{K}}_{\textrm{store}}
\dealcoloneq
        \begin{bmatrix}
			\mathbf{\mathsf{K}}_{\textrm{con}}	&	\mathbf{\mathsf{K}}_{u\widetilde{p}}	& \mathbf{0}
			\\
			\mathbf{\mathsf{K}}_{\widetilde{p}u}	&	\mathbf{0}	&	\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
			\\
			\mathbf{0}	& 	\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}		& \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\end{bmatrix} \, .
  @f] 。






<a name="Thematerialclass"></a><h3> The material class </h3>


一个好的面向对象的材料类的设计将有利于本教程扩展到广泛的材料类型。在本教程中，我们只有一个名为Material_Compressible_Neo_Hook_Three_Field的材料类。理想情况下，这个类会派生自超弹性材料（HyperelasticMaterial），而超弹性材料会派生自基类Material。这里使用的三场性质的表述也使问题复杂化。

三场公式的亥姆霍兹自由能函数为  $\Psi = \Psi_\text{vol}(\widetilde{J}) + \Psi_\text{iso}(\overline{\mathbf{b}})$  。Kirchhoff应力的等效部分 ${\boldsymbol{\tau}}_{\text{iso}}(\overline{\mathbf{b}})$ 与使用超弹性材料的单场公式得到的相同。然而，自由能的体积部分现在是一个主要变量的函数  $\widetilde{J}$  。因此，对于三场公式来说，基尔霍夫应力 ${\boldsymbol{\tau}}_{\text{vol}}$ 的体积部分的构成反应（和正切）并不像单场公式那样由超弹性构成法给出。我们可以将术语 $\boldsymbol{\tau}_{\textrm{vol}} \equiv \widetilde{p} J \mathbf{I}$ 标记为体积基尔霍夫应力，但压力 $\widetilde{p}$ 不是由自由能得出的；它是一个主场。

为了有一个灵活的方法，我们决定Material_Compressible_Neo_Hook_Three_Field仍然能够计算并返回一个体积Kirchhoff应力和正切。为了做到这一点，我们选择在与正交点相关的Material_Compressible_Neo_Hook_Three_Field类中存储插值的主域 $\widetilde{p}$ 和 $\widetilde{J}$ 。这个决定应该在以后的阶段，当教程扩展到考虑其他材料时，再重新审视。




<a name="Numericalexample"></a><h3> Numerical example </h3>


这里考虑的数值例子是一个压缩下的几乎不可压缩的块。这个基准问题取自

- S. Reese, P. Wriggers, B.D. Reddy (2000), A new locking-free brick element technique for large deformation problems in elasticity,  <em>  Computers and Structures  </em>  , <strong> 75</strong>, 291-304.   DOI:<a href="http://doi.org/10.1016/S0045-7949(99)00137-6">10.1016/S0045-7949(99)00137-6</a>。

   <img src="https://www.dealii.org/images/steps/developer/step-44.setup.png" alt=""> 

该材料是具有<a href="http://en.wikipedia.org/wiki/Shear_modulus">shear modulus</a> $\mu = 80.194e6$ 和 $\nu = 0.4999$ 的准不可压缩的新胡克式。对于这样一个材料特性的选择，传统的单场 $Q_1$ 方法将锁定。也就是说，响应会过于僵硬。初始和最终配置显示在上面的图片中。利用对称性，我们只求解四分之一的几何体（即一个尺寸为 $0.001$ 的立方体）。域的上表面的内四分之一受到 $p_0$ 的载荷。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * We start by including all the necessary deal.II header files and some C++
 * related ones. They have been discussed in detail in previous tutorial
 * programs, so you need only refer to past tutorials for details.
 * 
 * @code
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/parameter_handler.h>
 * #include <deal.II/base/point.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/symmetric_tensor.h>
 * #include <deal.II/base/tensor.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/work_stream.h>
 * #include <deal.II/dofs/dof_renumbering.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * @endcode
 * 
 * This header gives us the functionality to store
 * data at quadrature points
 * 
 * @code
 * #include <deal.II/base/quadrature_point_data.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/grid_in.h>
 * #include <deal.II/grid/tria.h>
 * 
 * #include <deal.II/fe/fe_dgp_monomial.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_tools.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/fe/mapping_q_eulerian.h>
 * 
 * #include <deal.II/lac/block_sparse_matrix.h>
 * #include <deal.II/lac/block_vector.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition_selector.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/solver_selector.h>
 * #include <deal.II/lac/sparse_direct.h>
 * #include <deal.II/lac/affine_constraints.h>
 * 
 * @endcode
 * 
 * Here are the headers necessary to use the LinearOperator class.
 * These are also all conveniently packaged into a single
 * header file, namely <deal.II/lac/linear_operator_tools.h>
 * but we list those specifically required here for the sake
 * of transparency.
 * 
 * @code
 * #include <deal.II/lac/linear_operator.h>
 * #include <deal.II/lac/packaged_operation.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * Defined in these two headers are some operations that are pertinent to
 * finite strain elasticity. The first will help us compute some kinematic
 * quantities, and the second provides some stanard tensor definitions.
 * 
 * @code
 * #include <deal.II/physics/elasticity/kinematics.h>
 * #include <deal.II/physics/elasticity/standard_tensors.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * 
 * @endcode
 * 
 * We then stick everything that relates to this tutorial program into a
 * namespace of its own, and import all the deal.II function and class names
 * into it:
 * 
 * @code
 * namespace Step44
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameters"></a> 
 * <h3>Run-time parameters</h3>
 *   

 * 
 * There are several parameters that can be set in the code so we set up a
 * ParameterHandler object to read in the choices at run-time.
 * 
 * @code
 *   namespace Parameters
 *   {
 * @endcode
 * 
 * 
 * <a name="FiniteElementsystem"></a> 
 * <h4>Finite Element system</h4>
 * 

 * 
 * As mentioned in the introduction, a different order interpolation should
 * be used for the displacement $\mathbf{u}$ than for the pressure
 * $\widetilde{p}$ and the dilatation $\widetilde{J}$.  Choosing
 * $\widetilde{p}$ and $\widetilde{J}$ as discontinuous (constant) functions
 * at the element level leads to the mean-dilatation method. The
 * discontinuous approximation allows $\widetilde{p}$ and $\widetilde{J}$ to
 * be condensed out and a classical displacement based method is recovered.
 * Here we specify the polynomial order used to approximate the solution.
 * The quadrature order should be adjusted accordingly.
 * 
 * @code
 *     struct FESystem
 *     {
 *       unsigned int poly_degree;
 *       unsigned int quad_order;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 * 
 *     void FESystem::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Finite element system");
 *       {
 *         prm.declare_entry("Polynomial degree",
 *                           "2",
 *                           Patterns::Integer(0),
 *                           "Displacement system polynomial order");
 * 
 *         prm.declare_entry("Quadrature order",
 *                           "3",
 *                           Patterns::Integer(0),
 *                           "Gauss quadrature order");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void FESystem::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Finite element system");
 *       {
 *         poly_degree = prm.get_integer("Polynomial degree");
 *         quad_order  = prm.get_integer("Quadrature order");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Geometry"></a> 
 * <h4>Geometry</h4>
 * 

 * 
 * Make adjustments to the problem geometry and the applied load.  Since the
 * problem modelled here is quite specific, the load scale can be altered to
 * specific values to compare with the results given in the literature.
 * 
 * @code
 *     struct Geometry
 *     {
 *       unsigned int global_refinement;
 *       double       scale;
 *       double       p_p0;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     void Geometry::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Geometry");
 *       {
 *         prm.declare_entry("Global refinement",
 *                           "2",
 *                           Patterns::Integer(0),
 *                           "Global refinement level");
 * 
 *         prm.declare_entry("Grid scale",
 *                           "1e-3",
 *                           Patterns::Double(0.0),
 *                           "Global grid scaling factor");
 * 
 *         prm.declare_entry("Pressure ratio p/p0",
 *                           "100",
 *                           Patterns::Selection("20|40|60|80|100"),
 *                           "Ratio of applied pressure to reference pressure");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void Geometry::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Geometry");
 *       {
 *         global_refinement = prm.get_integer("Global refinement");
 *         scale             = prm.get_double("Grid scale");
 *         p_p0              = prm.get_double("Pressure ratio p/p0");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Materials"></a> 
 * <h4>Materials</h4>
 * 

 * 
 * We also need the shear modulus $ \mu $ and Poisson ration $ \nu $ for the
 * neo-Hookean material.
 * 
 * @code
 *     struct Materials
 *     {
 *       double nu;
 *       double mu;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     void Materials::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Material properties");
 *       {
 *         prm.declare_entry("Poisson's ratio",
 *                           "0.4999",
 *                           Patterns::Double(-1.0, 0.5),
 *                           "Poisson's ratio");
 * 
 *         prm.declare_entry("Shear modulus",
 *                           "80.194e6",
 *                           Patterns::Double(),
 *                           "Shear modulus");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void Materials::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Material properties");
 *       {
 *         nu = prm.get_double("Poisson's ratio");
 *         mu = prm.get_double("Shear modulus");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Linearsolver"></a> 
 * <h4>Linear solver</h4>
 * 

 * 
 * Next, we choose both solver and preconditioner settings.  The use of an
 * effective preconditioner is critical to ensure convergence when a large
 * nonlinear motion occurs within a Newton increment.
 * 
 * @code
 *     struct LinearSolver
 *     {
 *       std::string type_lin;
 *       double      tol_lin;
 *       double      max_iterations_lin;
 *       bool        use_static_condensation;
 *       std::string preconditioner_type;
 *       double      preconditioner_relaxation;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     void LinearSolver::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Linear solver");
 *       {
 *         prm.declare_entry("Solver type",
 *                           "CG",
 *                           Patterns::Selection("CG|Direct"),
 *                           "Type of solver used to solve the linear system");
 * 
 *         prm.declare_entry("Residual",
 *                           "1e-6",
 *                           Patterns::Double(0.0),
 *                           "Linear solver residual (scaled by residual norm)");
 * 
 *         prm.declare_entry(
 *           "Max iteration multiplier",
 *           "1",
 *           Patterns::Double(0.0),
 *           "Linear solver iterations (multiples of the system matrix size)");
 * 
 *         prm.declare_entry("Use static condensation",
 *                           "true",
 *                           Patterns::Bool(),
 *                           "Solve the full block system or a reduced problem");
 * 
 *         prm.declare_entry("Preconditioner type",
 *                           "ssor",
 *                           Patterns::Selection("jacobi|ssor"),
 *                           "Type of preconditioner");
 * 
 *         prm.declare_entry("Preconditioner relaxation",
 *                           "0.65",
 *                           Patterns::Double(0.0),
 *                           "Preconditioner relaxation value");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void LinearSolver::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Linear solver");
 *       {
 *         type_lin                  = prm.get("Solver type");
 *         tol_lin                   = prm.get_double("Residual");
 *         max_iterations_lin        = prm.get_double("Max iteration multiplier");
 *         use_static_condensation   = prm.get_bool("Use static condensation");
 *         preconditioner_type       = prm.get("Preconditioner type");
 *         preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Nonlinearsolver"></a> 
 * <h4>Nonlinear solver</h4>
 * 

 * 
 * A Newton-Raphson scheme is used to solve the nonlinear system of
 * governing equations.  We now define the tolerances and the maximum number
 * of iterations for the Newton-Raphson nonlinear solver.
 * 
 * @code
 *     struct NonlinearSolver
 *     {
 *       unsigned int max_iterations_NR;
 *       double       tol_f;
 *       double       tol_u;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     void NonlinearSolver::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Nonlinear solver");
 *       {
 *         prm.declare_entry("Max iterations Newton-Raphson",
 *                           "10",
 *                           Patterns::Integer(0),
 *                           "Number of Newton-Raphson iterations allowed");
 * 
 *         prm.declare_entry("Tolerance force",
 *                           "1.0e-9",
 *                           Patterns::Double(0.0),
 *                           "Force residual tolerance");
 * 
 *         prm.declare_entry("Tolerance displacement",
 *                           "1.0e-6",
 *                           Patterns::Double(0.0),
 *                           "Displacement error tolerance");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void NonlinearSolver::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Nonlinear solver");
 *       {
 *         max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
 *         tol_f             = prm.get_double("Tolerance force");
 *         tol_u             = prm.get_double("Tolerance displacement");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Time"></a> 
 * <h4>Time</h4>
 * 

 * 
 * Set the timestep size $ \varDelta t $ and the simulation end-time.
 * 
 * @code
 *     struct Time
 *     {
 *       double delta_t;
 *       double end_time;
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     void Time::declare_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Time");
 *       {
 *         prm.declare_entry("End time", "1", Patterns::Double(), "End time");
 * 
 *         prm.declare_entry("Time step size",
 *                           "0.1",
 *                           Patterns::Double(),
 *                           "Time step size");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 *     void Time::parse_parameters(ParameterHandler &prm)
 *     {
 *       prm.enter_subsection("Time");
 *       {
 *         end_time = prm.get_double("End time");
 *         delta_t  = prm.get_double("Time step size");
 *       }
 *       prm.leave_subsection();
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Allparameters"></a> 
 * <h4>All parameters</h4>
 * 

 * 
 * Finally we consolidate all of the above structures into a single
 * container that holds all of our run-time selections.
 * 
 * @code
 *     struct AllParameters : public FESystem,
 *                            public Geometry,
 *                            public Materials,
 *                            public LinearSolver,
 *                            public NonlinearSolver,
 *                            public Time
 * 
 *     {
 *       AllParameters(const std::string &input_file);
 * 
 *       static void declare_parameters(ParameterHandler &prm);
 * 
 *       void parse_parameters(ParameterHandler &prm);
 *     };
 * 
 *     AllParameters::AllParameters(const std::string &input_file)
 *     {
 *       ParameterHandler prm;
 *       declare_parameters(prm);
 *       prm.parse_input(input_file);
 *       parse_parameters(prm);
 *     }
 * 
 *     void AllParameters::declare_parameters(ParameterHandler &prm)
 *     {
 *       FESystem::declare_parameters(prm);
 *       Geometry::declare_parameters(prm);
 *       Materials::declare_parameters(prm);
 *       LinearSolver::declare_parameters(prm);
 *       NonlinearSolver::declare_parameters(prm);
 *       Time::declare_parameters(prm);
 *     }
 * 
 *     void AllParameters::parse_parameters(ParameterHandler &prm)
 *     {
 *       FESystem::parse_parameters(prm);
 *       Geometry::parse_parameters(prm);
 *       Materials::parse_parameters(prm);
 *       LinearSolver::parse_parameters(prm);
 *       NonlinearSolver::parse_parameters(prm);
 *       Time::parse_parameters(prm);
 *     }
 *   } // namespace Parameters
 * 
 * @endcode
 * 
 * 
 * <a name="Timeclass"></a> 
 * <h3>Time class</h3>
 * 

 * 
 * A simple class to store time data. Its functioning is transparent so no
 * discussion is necessary. For simplicity we assume a constant time step
 * size.
 * 
 * @code
 *   class Time
 *   {
 *   public:
 *     Time(const double time_end, const double delta_t)
 *       : timestep(0)
 *       , time_current(0.0)
 *       , time_end(time_end)
 *       , delta_t(delta_t)
 *     {}
 * 
 *     virtual ~Time() = default;
 * 
 *     double current() const
 *     {
 *       return time_current;
 *     }
 *     double end() const
 *     {
 *       return time_end;
 *     }
 *     double get_delta_t() const
 *     {
 *       return delta_t;
 *     }
 *     unsigned int get_timestep() const
 *     {
 *       return timestep;
 *     }
 *     void increment()
 *     {
 *       time_current += delta_t;
 *       ++timestep;
 *     }
 * 
 *   private:
 *     unsigned int timestep;
 *     double       time_current;
 *     const double time_end;
 *     const double delta_t;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="CompressibleneoHookeanmaterialwithinathreefieldformulation"></a> 
 * <h3>Compressible neo-Hookean material within a three-field formulation</h3>
 * 

 * 
 * As discussed in the Introduction, Neo-Hookean materials are a type of
 * hyperelastic materials.  The entire domain is assumed to be composed of a
 * compressible neo-Hookean material.  This class defines the behavior of
 * this material within a three-field formulation.  Compressible neo-Hookean
 * materials can be described by a strain-energy function (SEF) $ \Psi =
 * \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(\widetilde{J})
 * $.
 *   

 * 
 * The isochoric response is given by $
 * \Psi_{\text{iso}}(\overline{\mathbf{b}}) = c_{1} [\overline{I}_{1} - 3] $
 * where $ c_{1} = \frac{\mu}{2} $ and $\overline{I}_{1}$ is the first
 * invariant of the left- or right-isochoric Cauchy-Green deformation tensors.
 * That is $\overline{I}_1 \dealcoloneq \textrm{tr}(\overline{\mathbf{b}})$.
 * In this example the SEF that governs the volumetric response is defined as
 * $ \Psi_{\text{vol}}(\widetilde{J}) = \kappa \frac{1}{4} [ \widetilde{J}^2 -
 * 1 - 2\textrm{ln}\; \widetilde{J} ]$, where $\kappa \dealcoloneq \lambda +
 * 2/3 \mu$ is the <a href="http://en.wikipedia.org/wiki/Bulk_modulus">bulk
 * modulus</a> and $\lambda$ is <a
 * href="http://en.wikipedia.org/wiki/Lam%C3%A9_parameters">Lam&eacute;'s
 * first parameter</a>.
 *   

 * 
 * The following class will be used to characterize the material we work with,
 * and provides a central point that one would need to modify if one were to
 * implement a different material model. For it to work, we will store one
 * object of this type per quadrature point, and in each of these objects
 * store the current state (characterized by the values or measures  of the
 * three fields) so that we can compute the elastic coefficients linearized
 * around the current state.
 * 
 * @code
 *   template <int dim>
 *   class Material_Compressible_Neo_Hook_Three_Field
 *   {
 *   public:
 *     Material_Compressible_Neo_Hook_Three_Field(const double mu, const double nu)
 *       : kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu)))
 *       , c_1(mu / 2.0)
 *       , det_F(1.0)
 *       , p_tilde(0.0)
 *       , J_tilde(1.0)
 *       , b_bar(Physics::Elasticity::StandardTensors<dim>::I)
 *     {
 *       Assert(kappa > 0, ExcInternalError());
 *     }
 * 
 * @endcode
 * 
 * We update the material model with various deformation dependent data
 * based on $F$ and the pressure $\widetilde{p}$ and dilatation
 * $\widetilde{J}$, and at the end of the function include a physical
 * check for internal consistency:
 * 
 * @code
 *     void update_material_data(const Tensor<2, dim> &F,
 *                               const double          p_tilde_in,
 *                               const double          J_tilde_in)
 *     {
 *       det_F                      = determinant(F);
 *       const Tensor<2, dim> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
 *       b_bar                      = Physics::Elasticity::Kinematics::b(F_bar);
 *       p_tilde                    = p_tilde_in;
 *       J_tilde                    = J_tilde_in;
 * 
 *       Assert(det_F > 0, ExcInternalError());
 *     }
 * 
 * @endcode
 * 
 * The second function determines the Kirchhoff stress $\boldsymbol{\tau}
 * = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau()
 *     {
 *       return get_tau_iso() + get_tau_vol();
 *     }
 * 
 * @endcode
 * 
 * The fourth-order elasticity tensor in the spatial setting
 * $\mathfrak{c}$ is calculated from the SEF $\Psi$ as $ J
 * \mathfrak{c}_{ijkl} = F_{iA} F_{jB} \mathfrak{C}_{ABCD} F_{kC} F_{lD}$
 * where $ \mathfrak{C} = 4 \frac{\partial^2 \Psi(\mathbf{C})}{\partial
 * \mathbf{C} \partial \mathbf{C}}$
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc() const
 *     {
 *       return get_Jc_vol() + get_Jc_iso();
 *     }
 * 
 * @endcode
 * 
 * Derivative of the volumetric free energy with respect to
 * $\widetilde{J}$ return $\frac{\partial
 * \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$
 * 
 * @code
 *     double get_dPsi_vol_dJ() const
 *     {
 *       return (kappa / 2.0) * (J_tilde - 1.0 / J_tilde);
 *     }
 * 
 * @endcode
 * 
 * Second derivative of the volumetric free energy wrt $\widetilde{J}$. We
 * need the following computation explicitly in the tangent so we make it
 * public.  We calculate $\frac{\partial^2
 * \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial
 * \widetilde{J}}$
 * 
 * @code
 *     double get_d2Psi_vol_dJ2() const
 *     {
 *       return ((kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde)));
 *     }
 * 
 * @endcode
 * 
 * The next few functions return various data that we choose to store with
 * the material:
 * 
 * @code
 *     double get_det_F() const
 *     {
 *       return det_F;
 *     }
 * 
 *     double get_p_tilde() const
 *     {
 *       return p_tilde;
 *     }
 * 
 *     double get_J_tilde() const
 *     {
 *       return J_tilde;
 *     }
 * 
 *   protected:
 * @endcode
 * 
 * Define constitutive model parameters $\kappa$ (bulk modulus) and the
 * neo-Hookean model parameter $c_1$:
 * 
 * @code
 *     const double kappa;
 *     const double c_1;
 * 
 * @endcode
 * 
 * Model specific data that is convenient to store with the material:
 * 
 * @code
 *     double                  det_F;
 *     double                  p_tilde;
 *     double                  J_tilde;
 *     SymmetricTensor<2, dim> b_bar;
 * 
 * @endcode
 * 
 * The following functions are used internally in determining the result
 * of some of the public functions above. The first one determines the
 * volumetric Kirchhoff stress $\boldsymbol{\tau}_{\textrm{vol}}$:
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_vol() const
 *     {
 *       return p_tilde * det_F * Physics::Elasticity::StandardTensors<dim>::I;
 *     }
 * 
 * @endcode
 * 
 * Next, determine the isochoric Kirchhoff stress
 * $\boldsymbol{\tau}_{\textrm{iso}} =
 * \mathcal{P}:\overline{\boldsymbol{\tau}}$:
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_iso() const
 *     {
 *       return Physics::Elasticity::StandardTensors<dim>::dev_P * get_tau_bar();
 *     }
 * 
 * @endcode
 * 
 * Then, determine the fictitious Kirchhoff stress
 * $\overline{\boldsymbol{\tau}}$:
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_bar() const
 *     {
 *       return 2.0 * c_1 * b_bar;
 *     }
 * 
 * @endcode
 * 
 * Calculate the volumetric part of the tangent $J
 * \mathfrak{c}_\textrm{vol}$:
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc_vol() const
 *     {
 *       return p_tilde * det_F *
 *              (Physics::Elasticity::StandardTensors<dim>::IxI -
 *               (2.0 * Physics::Elasticity::StandardTensors<dim>::S));
 *     }
 * 
 * @endcode
 * 
 * Calculate the isochoric part of the tangent $J
 * \mathfrak{c}_\textrm{iso}$:
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc_iso() const
 *     {
 *       const SymmetricTensor<2, dim> tau_bar = get_tau_bar();
 *       const SymmetricTensor<2, dim> tau_iso = get_tau_iso();
 *       const SymmetricTensor<4, dim> tau_iso_x_I =
 *         outer_product(tau_iso, Physics::Elasticity::StandardTensors<dim>::I);
 *       const SymmetricTensor<4, dim> I_x_tau_iso =
 *         outer_product(Physics::Elasticity::StandardTensors<dim>::I, tau_iso);
 *       const SymmetricTensor<4, dim> c_bar = get_c_bar();
 * 
 *       return (2.0 / dim) * trace(tau_bar) *
 *                Physics::Elasticity::StandardTensors<dim>::dev_P -
 *              (2.0 / dim) * (tau_iso_x_I + I_x_tau_iso) +
 *              Physics::Elasticity::StandardTensors<dim>::dev_P * c_bar *
 *                Physics::Elasticity::StandardTensors<dim>::dev_P;
 *     }
 * 
 * @endcode
 * 
 * Calculate the fictitious elasticity tensor $\overline{\mathfrak{c}}$.
 * For the material model chosen this is simply zero:
 * 
 * @code
 *     SymmetricTensor<4, dim> get_c_bar() const
 *     {
 *       return SymmetricTensor<4, dim>();
 *     }
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="Quadraturepointhistory"></a> 
 * <h3>Quadrature point history</h3>
 * 

 * 
 * As seen in step-18, the <code> PointHistory </code> class offers a method
 * for storing data at the quadrature points.  Here each quadrature point
 * holds a pointer to a material description.  Thus, different material models
 * can be used in different regions of the domain.  Among other data, we
 * choose to store the Kirchhoff stress $\boldsymbol{\tau}$ and the tangent
 * $J\mathfrak{c}$ for the quadrature points.
 * 
 * @code
 *   template <int dim>
 *   class PointHistory
 *   {
 *   public:
 *     PointHistory()
 *       : F_inv(Physics::Elasticity::StandardTensors<dim>::I)
 *       , tau(SymmetricTensor<2, dim>())
 *       , d2Psi_vol_dJ2(0.0)
 *       , dPsi_vol_dJ(0.0)
 *       , Jc(SymmetricTensor<4, dim>())
 *     {}
 * 
 *     virtual ~PointHistory() = default;
 * 
 * @endcode
 * 
 * The first function is used to create a material object and to
 * initialize all tensors correctly: The second one updates the stored
 * values and stresses based on the current deformation measure
 * $\textrm{Grad}\mathbf{u}_{\textrm{n}}$, pressure $\widetilde{p}$ and
 * dilation $\widetilde{J}$ field values.
 * 
 * @code
 *     void setup_lqp(const Parameters::AllParameters &parameters)
 *     {
 *       material =
 *         std::make_shared<Material_Compressible_Neo_Hook_Three_Field<dim>>(
 *           parameters.mu, parameters.nu);
 *       update_values(Tensor<2, dim>(), 0.0, 1.0);
 *     }
 * 
 * @endcode
 * 
 * To this end, we calculate the deformation gradient $\mathbf{F}$ from
 * the displacement gradient $\textrm{Grad}\ \mathbf{u}$, i.e.
 * $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad}\ \mathbf{u}$ and
 * then let the material model associated with this quadrature point
 * update itself. When computing the deformation gradient, we have to take
 * care with which data types we compare the sum $\mathbf{I} +
 * \textrm{Grad}\ \mathbf{u}$: Since $I$ has data type SymmetricTensor,
 * just writing <code>I + Grad_u_n</code> would convert the second
 * argument to a symmetric tensor, perform the sum, and then cast the
 * result to a Tensor (i.e., the type of a possibly nonsymmetric
 * tensor). However, since <code>Grad_u_n</code> is nonsymmetric in
 * general, the conversion to SymmetricTensor will fail. We can avoid this
 * back and forth by converting $I$ to Tensor first, and then performing
 * the addition as between nonsymmetric tensors:
 * 
 * @code
 *     void update_values(const Tensor<2, dim> &Grad_u_n,
 *                        const double          p_tilde,
 *                        const double          J_tilde)
 *     {
 *       const Tensor<2, dim> F = Physics::Elasticity::Kinematics::F(Grad_u_n);
 *       material->update_material_data(F, p_tilde, J_tilde);
 * 
 * @endcode
 * 
 * The material has been updated so we now calculate the Kirchhoff
 * stress $\mathbf{\tau}$, the tangent $J\mathfrak{c}$ and the first and
 * second derivatives of the volumetric free energy.
 *       

 * 
 * We also store the inverse of the deformation gradient since we
 * frequently use it:
 * 
 * @code
 *       F_inv         = invert(F);
 *       tau           = material->get_tau();
 *       Jc            = material->get_Jc();
 *       dPsi_vol_dJ   = material->get_dPsi_vol_dJ();
 *       d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();
 *     }
 * 
 * @endcode
 * 
 * We offer an interface to retrieve certain data.  Here are the kinematic
 * variables:
 * 
 * @code
 *     double get_J_tilde() const
 *     {
 *       return material->get_J_tilde();
 *     }
 * 
 *     double get_det_F() const
 *     {
 *       return material->get_det_F();
 *     }
 * 
 *     const Tensor<2, dim> &get_F_inv() const
 *     {
 *       return F_inv;
 *     }
 * 
 * @endcode
 * 
 * ...and the kinetic variables.  These are used in the material and
 * global tangent matrix and residual assembly operations:
 * 
 * @code
 *     double get_p_tilde() const
 *     {
 *       return material->get_p_tilde();
 *     }
 * 
 *     const SymmetricTensor<2, dim> &get_tau() const
 *     {
 *       return tau;
 *     }
 * 
 *     double get_dPsi_vol_dJ() const
 *     {
 *       return dPsi_vol_dJ;
 *     }
 * 
 *     double get_d2Psi_vol_dJ2() const
 *     {
 *       return d2Psi_vol_dJ2;
 *     }
 * 
 * @endcode
 * 
 * And finally the tangent:
 * 
 * @code
 *     const SymmetricTensor<4, dim> &get_Jc() const
 *     {
 *       return Jc;
 *     }
 * 
 * @endcode
 * 
 * In terms of member functions, this class stores for the quadrature
 * point it represents a copy of a material type in case different
 * materials are used in different regions of the domain, as well as the
 * inverse of the deformation gradient...
 * 
 * @code
 *   private:
 *     std::shared_ptr<Material_Compressible_Neo_Hook_Three_Field<dim>> material;
 * 
 *     Tensor<2, dim> F_inv;
 * 
 * @endcode
 * 
 * ... and stress-type variables along with the tangent $J\mathfrak{c}$:
 * 
 * @code
 *     SymmetricTensor<2, dim> tau;
 *     double                  d2Psi_vol_dJ2;
 *     double                  dPsi_vol_dJ;
 * 
 *     SymmetricTensor<4, dim> Jc;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Quasistaticquasiincompressiblefinitestrainsolid"></a> 
 * <h3>Quasi-static quasi-incompressible finite-strain solid</h3>
 * 

 * 
 * The Solid class is the central class in that it represents the problem at
 * hand. It follows the usual scheme in that all it really has is a
 * constructor, destructor and a <code>run()</code> function that dispatches
 * all the work to private functions of this class:
 * 
 * @code
 *   template <int dim>
 *   class Solid
 *   {
 *   public:
 *     Solid(const std::string &input_file);
 * 
 *     void run();
 * 
 *   private:
 * @endcode
 * 
 * In the private section of this class, we first forward declare a number
 * of objects that are used in parallelizing work using the WorkStream
 * object (see the @ref threads module for more information on this).
 *     

 * 
 * We declare such structures for the computation of tangent (stiffness)
 * matrix and right hand side vector, static condensation, and for updating
 * quadrature points:
 * 
 * @code
 *     struct PerTaskData_ASM;
 *     struct ScratchData_ASM;
 * 
 *     struct PerTaskData_SC;
 *     struct ScratchData_SC;
 * 
 *     struct PerTaskData_UQPH;
 *     struct ScratchData_UQPH;
 * 
 * @endcode
 * 
 * We start the collection of member functions with one that builds the
 * grid:
 * 
 * @code
 *     void make_grid();
 * 
 * @endcode
 * 
 * Set up the finite element system to be solved:
 * 
 * @code
 *     void system_setup();
 * 
 *     void determine_component_extractors();
 * 
 * @endcode
 * 
 * Create Dirichlet constraints for the incremental displacement field:
 * 
 * @code
 *     void make_constraints(const int it_nr);
 * 
 * @endcode
 * 
 * Several functions to assemble the system and right hand side matrices
 * using multithreading. Each of them comes as a wrapper function, one
 * that is executed to do the work in the WorkStream model on one cell,
 * and one that copies the work done on this one cell into the global
 * object that represents it:
 * 
 * @code
 *     void assemble_system();
 * 
 *     void assemble_system_one_cell(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       ScratchData_ASM &                                     scratch,
 *       PerTaskData_ASM &                                     data) const;
 * 
 * @endcode
 * 
 * And similar to perform global static condensation:
 * 
 * @code
 *     void assemble_sc();
 * 
 *     void assemble_sc_one_cell(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       ScratchData_SC &                                      scratch,
 *       PerTaskData_SC &                                      data);
 * 
 *     void copy_local_to_global_sc(const PerTaskData_SC &data);
 * 
 * @endcode
 * 
 * Create and update the quadrature points. Here, no data needs to be
 * copied into a global object, so the copy_local_to_global function is
 * empty:
 * 
 * @code
 *     void setup_qph();
 * 
 *     void update_qph_incremental(const BlockVector<double> &solution_delta);
 * 
 *     void update_qph_incremental_one_cell(
 *       const typename DoFHandler<dim>::active_cell_iterator &cell,
 *       ScratchData_UQPH &                                    scratch,
 *       PerTaskData_UQPH &                                    data);
 * 
 *     void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)
 *     {}
 * 
 * @endcode
 * 
 * Solve for the displacement using a Newton-Raphson method. We break this
 * function into the nonlinear loop and the function that solves the
 * linearized Newton-Raphson step:
 * 
 * @code
 *     void solve_nonlinear_timestep(BlockVector<double> &solution_delta);
 * 
 *     std::pair<unsigned int, double>
 *     solve_linear_system(BlockVector<double> &newton_update);
 * 
 * @endcode
 * 
 * Solution retrieval as well as post-processing and writing data to file:
 * 
 * @code
 *     BlockVector<double>
 *     get_total_solution(const BlockVector<double> &solution_delta) const;
 * 
 *     void output_results() const;
 * 
 * @endcode
 * 
 * Finally, some member variables that describe the current state: A
 * collection of the parameters used to describe the problem setup...
 * 
 * @code
 *     Parameters::AllParameters parameters;
 * 
 * @endcode
 * 
 * ...the volume of the reference configuration...
 * 
 * @code
 *     double vol_reference;
 * 
 * @endcode
 * 
 * ...and description of the geometry on which the problem is solved:
 * 
 * @code
 *     Triangulation<dim> triangulation;
 * 
 * @endcode
 * 
 * Also, keep track of the current time and the time spent evaluating
 * certain functions
 * 
 * @code
 *     Time                time;
 *     mutable TimerOutput timer;
 * 
 * @endcode
 * 
 * A storage object for quadrature point information. As opposed to
 * step-18, deal.II's native quadrature point data manager is employed
 * here.
 * 
 * @code
 *     CellDataStorage<typename Triangulation<dim>::cell_iterator,
 *                     PointHistory<dim>>
 *       quadrature_point_history;
 * 
 * @endcode
 * 
 * A description of the finite-element system including the displacement
 * polynomial degree, the degree-of-freedom handler, number of DoFs per
 * cell and the extractor objects used to retrieve information from the
 * solution vectors:
 * 
 * @code
 *     const unsigned int               degree;
 *     const FESystem<dim>              fe;
 *     DoFHandler<dim>                  dof_handler;
 *     const unsigned int               dofs_per_cell;
 *     const FEValuesExtractors::Vector u_fe;
 *     const FEValuesExtractors::Scalar p_fe;
 *     const FEValuesExtractors::Scalar J_fe;
 * 
 * @endcode
 * 
 * Description of how the block-system is arranged. There are 3 blocks,
 * the first contains a vector DOF $\mathbf{u}$ while the other two
 * describe scalar DOFs, $\widetilde{p}$ and $\widetilde{J}$.
 * 
 * @code
 *     static const unsigned int n_blocks          = 3;
 *     static const unsigned int n_components      = dim + 2;
 *     static const unsigned int first_u_component = 0;
 *     static const unsigned int p_component       = dim;
 *     static const unsigned int J_component       = dim + 1;
 * 
 *     enum
 *     {
 *       u_dof = 0,
 *       p_dof = 1,
 *       J_dof = 2
 *     };
 * 
 *     std::vector<types::global_dof_index> dofs_per_block;
 *     std::vector<types::global_dof_index> element_indices_u;
 *     std::vector<types::global_dof_index> element_indices_p;
 *     std::vector<types::global_dof_index> element_indices_J;
 * 
 * @endcode
 * 
 * Rules for Gauss-quadrature on both the cell and faces. The number of
 * quadrature points on both cells and faces is recorded.
 * 
 * @code
 *     const QGauss<dim>     qf_cell;
 *     const QGauss<dim - 1> qf_face;
 *     const unsigned int    n_q_points;
 *     const unsigned int    n_q_points_f;
 * 
 * @endcode
 * 
 * Objects that store the converged solution and right-hand side vectors,
 * as well as the tangent matrix. There is an AffineConstraints object used
 * to keep track of constraints.  We make use of a sparsity pattern
 * designed for a block system.
 * 
 * @code
 *     AffineConstraints<double> constraints;
 *     BlockSparsityPattern      sparsity_pattern;
 *     BlockSparseMatrix<double> tangent_matrix;
 *     BlockVector<double>       system_rhs;
 *     BlockVector<double>       solution_n;
 * 
 * @endcode
 * 
 * Then define a number of variables to store norms and update norms and
 * normalization factors.
 * 
 * @code
 *     struct Errors
 *     {
 *       Errors()
 *         : norm(1.0)
 *         , u(1.0)
 *         , p(1.0)
 *         , J(1.0)
 *       {}
 * 
 *       void reset()
 *       {
 *         norm = 1.0;
 *         u    = 1.0;
 *         p    = 1.0;
 *         J    = 1.0;
 *       }
 *       void normalize(const Errors &rhs)
 *       {
 *         if (rhs.norm != 0.0)
 *           norm /= rhs.norm;
 *         if (rhs.u != 0.0)
 *           u /= rhs.u;
 *         if (rhs.p != 0.0)
 *           p /= rhs.p;
 *         if (rhs.J != 0.0)
 *           J /= rhs.J;
 *       }
 * 
 *       double norm, u, p, J;
 *     };
 * 
 *     Errors error_residual, error_residual_0, error_residual_norm, error_update,
 *       error_update_0, error_update_norm;
 * 
 * @endcode
 * 
 * Methods to calculate error measures
 * 
 * @code
 *     void get_error_residual(Errors &error_residual);
 * 
 *     void get_error_update(const BlockVector<double> &newton_update,
 *                           Errors &                   error_update);
 * 
 *     std::pair<double, double> get_error_dilation() const;
 * 
 * @endcode
 * 
 * Compute the volume in the spatial configuration
 * 
 * @code
 *     double compute_vol_current() const;
 * 
 * @endcode
 * 
 * Print information to screen in a pleasing way...
 * 
 * @code
 *     static void print_conv_header();
 * 
 *     void print_conv_footer();
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeSolidcodeclass"></a> 
 * <h3>Implementation of the <code>Solid</code> class</h3>
 * 

 * 
 * 
 * <a name="Publicinterface"></a> 
 * <h4>Public interface</h4>
 * 

 * 
 * We initialize the Solid class using data extracted from the parameter file.
 * 
 * @code
 *   template <int dim>
 *   Solid<dim>::Solid(const std::string &input_file)
 *     : parameters(input_file)
 *     , vol_reference(0.)
 *     , triangulation(Triangulation<dim>::maximum_smoothing)
 *     , time(parameters.end_time, parameters.delta_t)
 *     , timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
 *     , degree(parameters.poly_degree)
 *     ,
 * @endcode
 * 
 * The Finite Element System is composed of dim continuous displacement
 * DOFs, and discontinuous pressure and dilatation DOFs. In an attempt to
 * satisfy the Babuska-Brezzi or LBB stability conditions (see Hughes
 * (2000)), we setup a $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$
 * system. $Q_2 \times DGPM_1 \times DGPM_1$ elements satisfy this
 * condition, while $Q_1 \times DGPM_0 \times DGPM_0$ elements do
 * not. However, it has been shown that the latter demonstrate good
 * convergence characteristics nonetheless.
 * 
 * @code
 *     fe(FE_Q<dim>(parameters.poly_degree),
 *        dim, // displacement
 *        FE_DGPMonomial<dim>(parameters.poly_degree - 1),
 *        1, // pressure
 *        FE_DGPMonomial<dim>(parameters.poly_degree - 1),
 *        1)
 *     , // dilatation
 *     dof_handler(triangulation)
 *     , dofs_per_cell(fe.n_dofs_per_cell())
 *     , u_fe(first_u_component)
 *     , p_fe(p_component)
 *     , J_fe(J_component)
 *     , dofs_per_block(n_blocks)
 *     , qf_cell(parameters.quad_order)
 *     , qf_face(parameters.quad_order)
 *     , n_q_points(qf_cell.size())
 *     , n_q_points_f(qf_face.size())
 *   {
 *     Assert(dim == 2 || dim == 3,
 *            ExcMessage("This problem only works in 2 or 3 space dimensions."));
 *     determine_component_extractors();
 *   }
 * 
 * 
 * @endcode
 * 
 * In solving the quasi-static problem, the time becomes a loading parameter,
 * i.e. we increasing the loading linearly with time, making the two concepts
 * interchangeable. We choose to increment time linearly using a constant time
 * step size.
 *   

 * 
 * We start the function with preprocessing, setting the initial dilatation
 * values, and then output the initial grid before starting the simulation
 * proper with the first time (and loading)
 * increment.
 *   

 * 
 * Care must be taken (or at least some thought given) when imposing the
 * constraint $\widetilde{J}=1$ on the initial solution field. The constraint
 * corresponds to the determinant of the deformation gradient in the
 * undeformed configuration, which is the identity tensor. We use
 * FE_DGPMonomial bases to interpolate the dilatation field, thus we can't
 * simply set the corresponding dof to unity as they correspond to the
 * monomial coefficients. Thus we use the VectorTools::project function to do
 * the work for us. The VectorTools::project function requires an argument
 * indicating the hanging node constraints. We have none in this program
 * So we have to create a constraint object. In its original state, constraint
 * objects are unsorted, and have to be sorted (using the
 * AffineConstraints::close function) before they can be used. Have a look at
 * step-21 for more information. We only need to enforce the initial condition
 * on the dilatation. In order to do this, we make use of a
 * ComponentSelectFunction which acts as a mask and sets the J_component of
 * n_components to 1. This is exactly what we want. Have a look at its usage
 * in step-20 for more information.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::run()
 *   {
 *     make_grid();
 *     system_setup();
 *     {
 *       AffineConstraints<double> constraints;
 *       constraints.close();
 * 
 *       const ComponentSelectFunction<dim> J_mask(J_component, n_components);
 * 
 *       VectorTools::project(
 *         dof_handler, constraints, QGauss<dim>(degree + 2), J_mask, solution_n);
 *     }
 *     output_results();
 *     time.increment();
 * 
 * @endcode
 * 
 * We then declare the incremental solution update $\varDelta
 * \mathbf{\Xi} \dealcoloneq \{\varDelta \mathbf{u},\varDelta \widetilde{p},
 * \varDelta \widetilde{J} \}$ and start the loop over the time domain.
 *     

 * 
 * At the beginning, we reset the solution update for this time step...
 * 
 * @code
 *     BlockVector<double> solution_delta(dofs_per_block);
 *     while (time.current() < time.end())
 *       {
 *         solution_delta = 0.0;
 * 
 * @endcode
 * 
 * ...solve the current time step and update total solution vector
 * $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} +
 * \varDelta \mathbf{\Xi}$...
 * 
 * @code
 *         solve_nonlinear_timestep(solution_delta);
 *         solution_n += solution_delta;
 * 
 * @endcode
 * 
 * ...and plot the results before moving on happily to the next time
 * step:
 * 
 * @code
 *         output_results();
 *         time.increment();
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Privateinterface"></a> 
 * <h3>Private interface</h3>
 * 

 * 
 * 
 * <a name="Threadingbuildingblocksstructures"></a> 
 * <h4>Threading-building-blocks structures</h4>
 * 

 * 
 * The first group of private member functions is related to parallelization.
 * We use the Threading Building Blocks library (TBB) to perform as many
 * computationally intensive distributed tasks as possible. In particular, we
 * assemble the tangent matrix and right hand side vector, the static
 * condensation contributions, and update data stored at the quadrature points
 * using TBB. Our main tool for this is the WorkStream class (see the @ref
 * threads module for more information).
 * 

 * 
 * Firstly we deal with the tangent matrix and right-hand side assembly
 * structures. The PerTaskData object stores local contributions to the global
 * system.
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::PerTaskData_ASM
 *   {
 *     FullMatrix<double>                   cell_matrix;
 *     Vector<double>                       cell_rhs;
 *     std::vector<types::global_dof_index> local_dof_indices;
 * 
 *     PerTaskData_ASM(const unsigned int dofs_per_cell)
 *       : cell_matrix(dofs_per_cell, dofs_per_cell)
 *       , cell_rhs(dofs_per_cell)
 *       , local_dof_indices(dofs_per_cell)
 *     {}
 * 
 *     void reset()
 *     {
 *       cell_matrix = 0.0;
 *       cell_rhs    = 0.0;
 *     }
 *   };
 * 
 * 
 * @endcode
 * 
 * On the other hand, the ScratchData object stores the larger objects such as
 * the shape-function values array (<code>Nx</code>) and a shape function
 * gradient and symmetric gradient vector which we will use during the
 * assembly.
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::ScratchData_ASM
 *   {
 *     FEValues<dim>     fe_values;
 *     FEFaceValues<dim> fe_face_values;
 * 
 *     std::vector<std::vector<double>>                  Nx;
 *     std::vector<std::vector<Tensor<2, dim>>>          grad_Nx;
 *     std::vector<std::vector<SymmetricTensor<2, dim>>> symm_grad_Nx;
 * 
 *     ScratchData_ASM(const FiniteElement<dim> &fe_cell,
 *                     const QGauss<dim> &       qf_cell,
 *                     const UpdateFlags         uf_cell,
 *                     const QGauss<dim - 1> &   qf_face,
 *                     const UpdateFlags         uf_face)
 *       : fe_values(fe_cell, qf_cell, uf_cell)
 *       , fe_face_values(fe_cell, qf_face, uf_face)
 *       , Nx(qf_cell.size(), std::vector<double>(fe_cell.n_dofs_per_cell()))
 *       , grad_Nx(qf_cell.size(),
 *                 std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
 *       , symm_grad_Nx(qf_cell.size(),
 *                      std::vector<SymmetricTensor<2, dim>>(
 *                        fe_cell.n_dofs_per_cell()))
 *     {}
 * 
 *     ScratchData_ASM(const ScratchData_ASM &rhs)
 *       : fe_values(rhs.fe_values.get_fe(),
 *                   rhs.fe_values.get_quadrature(),
 *                   rhs.fe_values.get_update_flags())
 *       , fe_face_values(rhs.fe_face_values.get_fe(),
 *                        rhs.fe_face_values.get_quadrature(),
 *                        rhs.fe_face_values.get_update_flags())
 *       , Nx(rhs.Nx)
 *       , grad_Nx(rhs.grad_Nx)
 *       , symm_grad_Nx(rhs.symm_grad_Nx)
 *     {}
 * 
 *     void reset()
 *     {
 *       const unsigned int n_q_points      = Nx.size();
 *       const unsigned int n_dofs_per_cell = Nx[0].size();
 *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *         {
 *           Assert(Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
 *           Assert(grad_Nx[q_point].size() == n_dofs_per_cell,
 *                  ExcInternalError());
 *           Assert(symm_grad_Nx[q_point].size() == n_dofs_per_cell,
 *                  ExcInternalError());
 *           for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
 *             {
 *               Nx[q_point][k]           = 0.0;
 *               grad_Nx[q_point][k]      = 0.0;
 *               symm_grad_Nx[q_point][k] = 0.0;
 *             }
 *         }
 *     }
 *   };
 * 
 * 
 * @endcode
 * 
 * Then we define structures to assemble the statically condensed tangent
 * matrix. Recall that we wish to solve for a displacement-based formulation.
 * We do the condensation at the element level as the $\widetilde{p}$ and
 * $\widetilde{J}$ fields are element-wise discontinuous.  As these operations
 * are matrix-based, we need to setup a number of matrices to store the local
 * contributions from a number of the tangent matrix sub-blocks.  We place
 * these in the PerTaskData struct.
 *   

 * 
 * We choose not to reset any data in the <code>reset()</code> function as the
 * matrix extraction and replacement tools will take care of this
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::PerTaskData_SC
 *   {
 *     FullMatrix<double>                   cell_matrix;
 *     std::vector<types::global_dof_index> local_dof_indices;
 * 
 *     FullMatrix<double> k_orig;
 *     FullMatrix<double> k_pu;
 *     FullMatrix<double> k_pJ;
 *     FullMatrix<double> k_JJ;
 *     FullMatrix<double> k_pJ_inv;
 *     FullMatrix<double> k_bbar;
 *     FullMatrix<double> A;
 *     FullMatrix<double> B;
 *     FullMatrix<double> C;
 * 
 *     PerTaskData_SC(const unsigned int dofs_per_cell,
 *                    const unsigned int n_u,
 *                    const unsigned int n_p,
 *                    const unsigned int n_J)
 *       : cell_matrix(dofs_per_cell, dofs_per_cell)
 *       , local_dof_indices(dofs_per_cell)
 *       , k_orig(dofs_per_cell, dofs_per_cell)
 *       , k_pu(n_p, n_u)
 *       , k_pJ(n_p, n_J)
 *       , k_JJ(n_J, n_J)
 *       , k_pJ_inv(n_p, n_J)
 *       , k_bbar(n_u, n_u)
 *       , A(n_J, n_u)
 *       , B(n_J, n_u)
 *       , C(n_p, n_u)
 *     {}
 * 
 *     void reset()
 *     {}
 *   };
 * 
 * 
 * @endcode
 * 
 * The ScratchData object for the operations we wish to perform here is empty
 * since we need no temporary data, but it still needs to be defined for the
 * current implementation of TBB in deal.II.  So we create a dummy struct for
 * this purpose.
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::ScratchData_SC
 *   {
 *     void reset()
 *     {}
 *   };
 * 
 * 
 * @endcode
 * 
 * And finally we define the structures to assist with updating the quadrature
 * point information. Similar to the SC assembly process, we do not need the
 * PerTaskData object (since there is nothing to store here) but must define
 * one nonetheless. Note that this is because for the operation that we have
 * here -- updating the data on quadrature points -- the operation is purely
 * local: the things we do on every cell get consumed on every cell, without
 * any global aggregation operation as is usually the case when using the
 * WorkStream class. The fact that we still have to define a per-task data
 * structure points to the fact that the WorkStream class may be ill-suited to
 * this operation (we could, in principle simply create a new task using
 * Threads::new_task for each cell) but there is not much harm done to doing
 * it this way anyway.
 * Furthermore, should there be different material models associated with a
 * quadrature point, requiring varying levels of computational expense, then
 * the method used here could be advantageous.
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::PerTaskData_UQPH
 *   {
 *     void reset()
 *     {}
 *   };
 * 
 * 
 * @endcode
 * 
 * The ScratchData object will be used to store an alias for the solution
 * vector so that we don't have to copy this large data structure. We then
 * define a number of vectors to extract the solution values and gradients at
 * the quadrature points.
 * 
 * @code
 *   template <int dim>
 *   struct Solid<dim>::ScratchData_UQPH
 *   {
 *     const BlockVector<double> &solution_total;
 * 
 *     std::vector<Tensor<2, dim>> solution_grads_u_total;
 *     std::vector<double>         solution_values_p_total;
 *     std::vector<double>         solution_values_J_total;
 * 
 *     FEValues<dim> fe_values;
 * 
 *     ScratchData_UQPH(const FiniteElement<dim> & fe_cell,
 *                      const QGauss<dim> &        qf_cell,
 *                      const UpdateFlags          uf_cell,
 *                      const BlockVector<double> &solution_total)
 *       : solution_total(solution_total)
 *       , solution_grads_u_total(qf_cell.size())
 *       , solution_values_p_total(qf_cell.size())
 *       , solution_values_J_total(qf_cell.size())
 *       , fe_values(fe_cell, qf_cell, uf_cell)
 *     {}
 * 
 *     ScratchData_UQPH(const ScratchData_UQPH &rhs)
 *       : solution_total(rhs.solution_total)
 *       , solution_grads_u_total(rhs.solution_grads_u_total)
 *       , solution_values_p_total(rhs.solution_values_p_total)
 *       , solution_values_J_total(rhs.solution_values_J_total)
 *       , fe_values(rhs.fe_values.get_fe(),
 *                   rhs.fe_values.get_quadrature(),
 *                   rhs.fe_values.get_update_flags())
 *     {}
 * 
 *     void reset()
 *     {
 *       const unsigned int n_q_points = solution_grads_u_total.size();
 *       for (unsigned int q = 0; q < n_q_points; ++q)
 *         {
 *           solution_grads_u_total[q]  = 0.0;
 *           solution_values_p_total[q] = 0.0;
 *           solution_values_J_total[q] = 0.0;
 *         }
 *     }
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidmake_grid"></a> 
 * <h4>Solid::make_grid</h4>
 * 

 * 
 * On to the first of the private member functions. Here we create the
 * triangulation of the domain, for which we choose the scaled cube with each
 * face given a boundary ID number.  The grid must be refined at least once
 * for the indentation problem.
 *   

 * 
 * We then determine the volume of the reference configuration and print it
 * for comparison:
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::make_grid()
 *   {
 *     GridGenerator::hyper_rectangle(
 *       triangulation,
 *       (dim == 3 ? Point<dim>(0.0, 0.0, 0.0) : Point<dim>(0.0, 0.0)),
 *       (dim == 3 ? Point<dim>(1.0, 1.0, 1.0) : Point<dim>(1.0, 1.0)),
 *       true);
 *     GridTools::scale(parameters.scale, triangulation);
 *     triangulation.refine_global(std::max(1U, parameters.global_refinement));
 * 
 *     vol_reference = GridTools::volume(triangulation);
 *     std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
 * 
 * @endcode
 * 
 * Since we wish to apply a Neumann BC to a patch on the top surface, we
 * must find the cell faces in this part of the domain and mark them with
 * a distinct boundary ID number.  The faces we are looking for are on the
 * +y surface and will get boundary ID 6 (zero through five are already
 * used when creating the six faces of the cube domain):
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       for (const auto &face : cell->face_iterators())
 *         {
 *           if (face->at_boundary() == true &&
 *               face->center()[1] == 1.0 * parameters.scale)
 *             {
 *               if (dim == 3)
 *                 {
 *                   if (face->center()[0] < 0.5 * parameters.scale &&
 *                       face->center()[2] < 0.5 * parameters.scale)
 *                     face->set_boundary_id(6);
 *                 }
 *               else
 *                 {
 *                   if (face->center()[0] < 0.5 * parameters.scale)
 *                     face->set_boundary_id(6);
 *                 }
 *             }
 *         }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidsystem_setup"></a> 
 * <h4>Solid::system_setup</h4>
 * 

 * 
 * Next we describe how the FE system is setup.  We first determine the number
 * of components per block. Since the displacement is a vector component, the
 * first dim components belong to it, while the next two describe scalar
 * pressure and dilatation DOFs.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::system_setup()
 *   {
 *     timer.enter_subsection("Setup system");
 * 
 *     std::vector<unsigned int> block_component(n_components,
 *                                               u_dof); // Displacement
 *     block_component[p_component] = p_dof;             // Pressure
 *     block_component[J_component] = J_dof;             // Dilatation
 * 
 * @endcode
 * 
 * The DOF handler is then initialized and we renumber the grid in an
 * efficient manner. We also record the number of DOFs per block.
 * 
 * @code
 *     dof_handler.distribute_dofs(fe);
 *     DoFRenumbering::Cuthill_McKee(dof_handler);
 *     DoFRenumbering::component_wise(dof_handler, block_component);
 * 
 *     dofs_per_block =
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
 * 
 *     std::cout << "Triangulation:"
 *               << "\n\t Number of active cells: "
 *               << triangulation.n_active_cells()
 *               << "\n\t Number of degrees of freedom: " << dof_handler.n_dofs()
 *               << std::endl;
 * 
 * @endcode
 * 
 * Setup the sparsity pattern and tangent matrix
 * 
 * @code
 *     tangent_matrix.clear();
 *     {
 *       const types::global_dof_index n_dofs_u = dofs_per_block[u_dof];
 *       const types::global_dof_index n_dofs_p = dofs_per_block[p_dof];
 *       const types::global_dof_index n_dofs_J = dofs_per_block[J_dof];
 * 
 *       BlockDynamicSparsityPattern dsp(n_blocks, n_blocks);
 * 
 *       dsp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
 *       dsp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
 *       dsp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J);
 * 
 *       dsp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
 *       dsp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
 *       dsp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J);
 * 
 *       dsp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u);
 *       dsp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p);
 *       dsp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J);
 *       dsp.collect_sizes();
 * 
 * @endcode
 * 
 * The global system matrix initially has the following structure
 * @f{align*}
 * \underbrace{\begin{bmatrix}
 * \mathsf{\mathbf{K}}_{uu}  & \mathsf{\mathbf{K}}_{u\widetilde{p}} &
 * \mathbf{0}
 * \\ \mathsf{\mathbf{K}}_{\widetilde{p}u} & \mathbf{0} &
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}
 * \\ \mathbf{0} & \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}} &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * \end{bmatrix}}_{\mathsf{\mathbf{K}}(\mathbf{\Xi}_{\textrm{i}})}
 * \underbrace{\begin{bmatrix}
 * d \mathsf{u}
 * \\  d \widetilde{\mathsf{\mathbf{p}}}
 * \\  d \widetilde{\mathsf{\mathbf{J}}}
 * \end{bmatrix}}_{d \mathbf{\Xi}}
 * =
 * \underbrace{\begin{bmatrix}
 * \mathsf{\mathbf{F}}_{u}(\mathbf{u}_{\textrm{i}})
 * \\ \mathsf{\mathbf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}})
 * \\ \mathsf{\mathbf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
 * \end{bmatrix}}_{ \mathsf{\mathbf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .
 * @f}
 * We optimize the sparsity pattern to reflect this structure
 * and prevent unnecessary data creation for the right-diagonal
 * block components.
 * 
 * @code
 *       Table<2, DoFTools::Coupling> coupling(n_components, n_components);
 *       for (unsigned int ii = 0; ii < n_components; ++ii)
 *         for (unsigned int jj = 0; jj < n_components; ++jj)
 *           if (((ii < p_component) && (jj == J_component)) ||
 *               ((ii == J_component) && (jj < p_component)) ||
 *               ((ii == p_component) && (jj == p_component)))
 *             coupling[ii][jj] = DoFTools::none;
 *           else
 *             coupling[ii][jj] = DoFTools::always;
 *       DoFTools::make_sparsity_pattern(
 *         dof_handler, coupling, dsp, constraints, false);
 *       sparsity_pattern.copy_from(dsp);
 *     }
 * 
 *     tangent_matrix.reinit(sparsity_pattern);
 * 
 * @endcode
 * 
 * We then set up storage vectors
 * 
 * @code
 *     system_rhs.reinit(dofs_per_block);
 *     system_rhs.collect_sizes();
 * 
 *     solution_n.reinit(dofs_per_block);
 *     solution_n.collect_sizes();
 * 
 * @endcode
 * 
 * ...and finally set up the quadrature
 * point history:
 * 
 * @code
 *     setup_qph();
 * 
 *     timer.leave_subsection();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Soliddetermine_component_extractors"></a> 
 * <h4>Solid::determine_component_extractors</h4>
 * Next we compute some information from the FE system that describes which
 * local element DOFs are attached to which block component.  This is used
 * later to extract sub-blocks from the global matrix.
 *   

 * 
 * In essence, all we need is for the FESystem object to indicate to which
 * block component a DOF on the reference cell is attached to.  Currently, the
 * interpolation fields are setup such that 0 indicates a displacement DOF, 1
 * a pressure DOF and 2 a dilatation DOF.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::determine_component_extractors()
 *   {
 *     element_indices_u.clear();
 *     element_indices_p.clear();
 *     element_indices_J.clear();
 * 
 *     for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k)
 *       {
 *         const unsigned int k_group = fe.system_to_base_index(k).first.first;
 *         if (k_group == u_dof)
 *           element_indices_u.push_back(k);
 *         else if (k_group == p_dof)
 *           element_indices_p.push_back(k);
 *         else if (k_group == J_dof)
 *           element_indices_J.push_back(k);
 *         else
 *           {
 *             Assert(k_group <= J_dof, ExcInternalError());
 *           }
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solidsetup_qph"></a> 
 * <h4>Solid::setup_qph</h4>
 * The method used to store quadrature information is already described in
 * step-18. Here we implement a similar setup for a SMP machine.
 *   

 * 
 * Firstly the actual QPH data objects are created. This must be done only
 * once the grid is refined to its finest level.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::setup_qph()
 *   {
 *     std::cout << "    Setting up quadrature point data..." << std::endl;
 * 
 *     quadrature_point_history.initialize(triangulation.begin_active(),
 *                                         triangulation.end(),
 *                                         n_q_points);
 * 
 * @endcode
 * 
 * Next we setup the initial quadrature point data.
 * Note that when the quadrature point data is retrieved,
 * it is returned as a vector of smart pointers.
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       {
 *         const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
 *           quadrature_point_history.get_data(cell);
 *         Assert(lqph.size() == n_q_points, ExcInternalError());
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
 *           lqph[q_point]->setup_lqp(parameters);
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solidupdate_qph_incremental"></a> 
 * <h4>Solid::update_qph_incremental</h4>
 * As the update of QP information occurs frequently and involves a number of
 * expensive operations, we define a multithreaded approach to distributing
 * the task across a number of CPU cores.
 *   

 * 
 * To start this, we first we need to obtain the total solution as it stands
 * at this Newton increment and then create the initial copy of the scratch
 * and copy data objects:
 * 
 * @code
 *   template <int dim>
 *   void
 *   Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta)
 *   {
 *     timer.enter_subsection("Update QPH data");
 *     std::cout << " UQPH " << std::flush;
 * 
 *     const BlockVector<double> solution_total(
 *       get_total_solution(solution_delta));
 * 
 *     const UpdateFlags uf_UQPH(update_values | update_gradients);
 *     PerTaskData_UQPH  per_task_data_UQPH;
 *     ScratchData_UQPH  scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);
 * 
 * @endcode
 * 
 * We then pass them and the one-cell update function to the WorkStream to
 * be processed:
 * 
 * @code
 *     WorkStream::run(dof_handler.active_cell_iterators(),
 *                     *this,
 *                     &Solid::update_qph_incremental_one_cell,
 *                     &Solid::copy_local_to_global_UQPH,
 *                     scratch_data_UQPH,
 *                     per_task_data_UQPH);
 * 
 *     timer.leave_subsection();
 *   }
 * 
 * 
 * @endcode
 * 
 * Now we describe how we extract data from the solution vector and pass it
 * along to each QP storage object for processing.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::update_qph_incremental_one_cell(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     ScratchData_UQPH &                                    scratch,
 *     PerTaskData_UQPH & /*data*/)
 *   {
 *     const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
 *       quadrature_point_history.get_data(cell);
 *     Assert(lqph.size() == n_q_points, ExcInternalError());
 * 
 *     Assert(scratch.solution_grads_u_total.size() == n_q_points,
 *            ExcInternalError());
 *     Assert(scratch.solution_values_p_total.size() == n_q_points,
 *            ExcInternalError());
 *     Assert(scratch.solution_values_J_total.size() == n_q_points,
 *            ExcInternalError());
 * 
 *     scratch.reset();
 * 
 * @endcode
 * 
 * We first need to find the values and gradients at quadrature points
 * inside the current cell and then we update each local QP using the
 * displacement gradient and total pressure and dilatation solution
 * values:
 * 
 * @code
 *     scratch.fe_values.reinit(cell);
 *     scratch.fe_values[u_fe].get_function_gradients(
 *       scratch.solution_total, scratch.solution_grads_u_total);
 *     scratch.fe_values[p_fe].get_function_values(
 *       scratch.solution_total, scratch.solution_values_p_total);
 *     scratch.fe_values[J_fe].get_function_values(
 *       scratch.solution_total, scratch.solution_values_J_total);
 * 
 *     for (const unsigned int q_point :
 *          scratch.fe_values.quadrature_point_indices())
 *       lqph[q_point]->update_values(scratch.solution_grads_u_total[q_point],
 *                                    scratch.solution_values_p_total[q_point],
 *                                    scratch.solution_values_J_total[q_point]);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidsolve_nonlinear_timestep"></a> 
 * <h4>Solid::solve_nonlinear_timestep</h4>
 * 

 * 
 * The next function is the driver method for the Newton-Raphson scheme. At
 * its top we create a new vector to store the current Newton update step,
 * reset the error storage objects and print solver header.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::solve_nonlinear_timestep(BlockVector<double> &solution_delta)
 *   {
 *     std::cout << std::endl
 *               << "Timestep " << time.get_timestep() << " @ " << time.current()
 *               << "s" << std::endl;
 * 
 *     BlockVector<double> newton_update(dofs_per_block);
 * 
 *     error_residual.reset();
 *     error_residual_0.reset();
 *     error_residual_norm.reset();
 *     error_update.reset();
 *     error_update_0.reset();
 *     error_update_norm.reset();
 * 
 *     print_conv_header();
 * 
 * @endcode
 * 
 * We now perform a number of Newton iterations to iteratively solve the
 * nonlinear problem.  Since the problem is fully nonlinear and we are
 * using a full Newton method, the data stored in the tangent matrix and
 * right-hand side vector is not reusable and must be cleared at each
 * Newton step. We then initially build the linear system and
 * check for convergence (and store this value in the first iteration).
 * The unconstrained DOFs of the rhs vector hold the out-of-balance
 * forces, and collectively determine whether or not the equilibrium
 * solution has been attained.
 *     

 * 
 * Although for this particular problem we could potentially construct the
 * RHS vector before assembling the system matrix, for the sake of
 * extensibility we choose not to do so. The benefit to assembling the RHS
 * vector and system matrix separately is that the latter is an expensive
 * operation and we can potentially avoid an extra assembly process by not
 * assembling the tangent matrix when convergence is attained. However, this
 * makes parallelizing the code using MPI more difficult. Furthermore, when
 * extending the problem to the transient case additional contributions to
 * the RHS may result from the time discretization and application of
 * constraints for the velocity and acceleration fields.
 * 
 * @code
 *     unsigned int newton_iteration = 0;
 *     for (; newton_iteration < parameters.max_iterations_NR; ++newton_iteration)
 *       {
 *         std::cout << " " << std::setw(2) << newton_iteration << " "
 *                   << std::flush;
 * 
 * @endcode
 * 
 * We construct the linear system, but hold off on solving it
 * (a step that should be significantly more expensive than assembly):
 * 
 * @code
 *         make_constraints(newton_iteration);
 *         assemble_system();
 * 
 * @endcode
 * 
 * We can now determine the normalized residual error and check for
 * solution convergence:
 * 
 * @code
 *         get_error_residual(error_residual);
 *         if (newton_iteration == 0)
 *           error_residual_0 = error_residual;
 * 
 *         error_residual_norm = error_residual;
 *         error_residual_norm.normalize(error_residual_0);
 * 
 *         if (newton_iteration > 0 && error_update_norm.u <= parameters.tol_u &&
 *             error_residual_norm.u <= parameters.tol_f)
 *           {
 *             std::cout << " CONVERGED! " << std::endl;
 *             print_conv_footer();
 * 
 *             break;
 *           }
 * 
 * @endcode
 * 
 * If we have decided that we want to continue with the iteration, we
 * solve the linearized system:
 * 
 * @code
 *         const std::pair<unsigned int, double> lin_solver_output =
 *           solve_linear_system(newton_update);
 * 
 * @endcode
 * 
 * We can now determine the normalized Newton update error:
 * 
 * @code
 *         get_error_update(newton_update, error_update);
 *         if (newton_iteration == 0)
 *           error_update_0 = error_update;
 * 
 *         error_update_norm = error_update;
 *         error_update_norm.normalize(error_update_0);
 * 
 * @endcode
 * 
 * Lastly, since we implicitly accept the solution step we can perform
 * the actual update of the solution increment for the current time
 * step, update all quadrature point information pertaining to
 * this new displacement and stress state and continue iterating:
 * 
 * @code
 *         solution_delta += newton_update;
 *         update_qph_incremental(solution_delta);
 * 
 *         std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
 *                   << std::scientific << lin_solver_output.first << "  "
 *                   << lin_solver_output.second << "  "
 *                   << error_residual_norm.norm << "  " << error_residual_norm.u
 *                   << "  " << error_residual_norm.p << "  "
 *                   << error_residual_norm.J << "  " << error_update_norm.norm
 *                   << "  " << error_update_norm.u << "  " << error_update_norm.p
 *                   << "  " << error_update_norm.J << "  " << std::endl;
 *       }
 * 
 * @endcode
 * 
 * At the end, if it turns out that we have in fact done more iterations
 * than the parameter file allowed, we raise an exception that can be
 * caught in the main() function. The call <code>AssertThrow(condition,
 * exc_object)</code> is in essence equivalent to <code>if (!cond) throw
 * exc_object;</code> but the former form fills certain fields in the
 * exception object that identify the location (filename and line number)
 * where the exception was raised to make it simpler to identify where the
 * problem happened.
 * 
 * @code
 *     AssertThrow(newton_iteration < parameters.max_iterations_NR,
 *                 ExcMessage("No convergence in nonlinear solver!"));
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidprint_conv_headerandSolidprint_conv_footer"></a> 
 * <h4>Solid::print_conv_header and Solid::print_conv_footer</h4>
 * 

 * 
 * This program prints out data in a nice table that is updated
 * on a per-iteration basis. The next two functions set up the table
 * header and footer:
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::print_conv_header()
 *   {
 *     static const unsigned int l_width = 150;
 * 
 *     for (unsigned int i = 0; i < l_width; ++i)
 *       std::cout << "_";
 *     std::cout << std::endl;
 * 
 *     std::cout << "               SOLVER STEP               "
 *               << " |  LIN_IT   LIN_RES    RES_NORM    "
 *               << " RES_U     RES_P      RES_J     NU_NORM     "
 *               << " NU_U       NU_P       NU_J " << std::endl;
 * 
 *     for (unsigned int i = 0; i < l_width; ++i)
 *       std::cout << "_";
 *     std::cout << std::endl;
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   void Solid<dim>::print_conv_footer()
 *   {
 *     static const unsigned int l_width = 150;
 * 
 *     for (unsigned int i = 0; i < l_width; ++i)
 *       std::cout << "_";
 *     std::cout << std::endl;
 * 
 *     const std::pair<double, double> error_dil = get_error_dilation();
 * 
 *     std::cout << "Relative errors:" << std::endl
 *               << "Displacement:\t" << error_update.u / error_update_0.u
 *               << std::endl
 *               << "Force: \t\t" << error_residual.u / error_residual_0.u
 *               << std::endl
 *               << "Dilatation:\t" << error_dil.first << std::endl
 *               << "v / V_0:\t" << error_dil.second * vol_reference << " / "
 *               << vol_reference << " = " << error_dil.second << std::endl;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_dilation"></a> 
 * <h4>Solid::get_error_dilation</h4>
 * 

 * 
 * Calculate the volume of the domain in the spatial configuration
 * 
 * @code
 *   template <int dim>
 *   double Solid<dim>::compute_vol_current() const
 *   {
 *     double vol_current = 0.0;
 * 
 *     FEValues<dim> fe_values(fe, qf_cell, update_JxW_values);
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 * 
 * @endcode
 * 
 * In contrast to that which was previously called for,
 * in this instance the quadrature point data is specifically
 * non-modifiable since we will only be accessing data.
 * We ensure that the right get_data function is called by
 * marking this update function as constant.
 * 
 * @code
 *         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
 *           quadrature_point_history.get_data(cell);
 *         Assert(lqph.size() == n_q_points, ExcInternalError());
 * 
 *         for (const unsigned int q_point : fe_values.quadrature_point_indices())
 *           {
 *             const double det_F_qp = lqph[q_point]->get_det_F();
 *             const double JxW      = fe_values.JxW(q_point);
 * 
 *             vol_current += det_F_qp * JxW;
 *           }
 *       }
 *     Assert(vol_current > 0.0, ExcInternalError());
 *     return vol_current;
 *   }
 * 
 * @endcode
 * 
 * Calculate how well the dilatation $\widetilde{J}$ agrees with $J
 * \dealcoloneq \textrm{det}\ \mathbf{F}$ from the $L^2$ error $ \bigl[
 * \int_{\Omega_0} {[ J - \widetilde{J}]}^{2}\textrm{d}V \bigr]^{1/2}$.
 * We also return the ratio of the current volume of the
 * domain to the reference volume. This is of interest for incompressible
 * media where we want to check how well the isochoric constraint has been
 * enforced.
 * 
 * @code
 *   template <int dim>
 *   std::pair<double, double> Solid<dim>::get_error_dilation() const
 *   {
 *     double dil_L2_error = 0.0;
 * 
 *     FEValues<dim> fe_values(fe, qf_cell, update_JxW_values);
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       {
 *         fe_values.reinit(cell);
 * 
 *         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
 *           quadrature_point_history.get_data(cell);
 *         Assert(lqph.size() == n_q_points, ExcInternalError());
 * 
 *         for (const unsigned int q_point : fe_values.quadrature_point_indices())
 *           {
 *             const double det_F_qp   = lqph[q_point]->get_det_F();
 *             const double J_tilde_qp = lqph[q_point]->get_J_tilde();
 *             const double the_error_qp_squared =
 *               std::pow((det_F_qp - J_tilde_qp), 2);
 *             const double JxW = fe_values.JxW(q_point);
 * 
 *             dil_L2_error += the_error_qp_squared * JxW;
 *           }
 *       }
 * 
 *     return std::make_pair(std::sqrt(dil_L2_error),
 *                           compute_vol_current() / vol_reference);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_residual"></a> 
 * <h4>Solid::get_error_residual</h4>
 * 

 * 
 * Determine the true residual error for the problem.  That is, determine the
 * error in the residual for the unconstrained degrees of freedom.  Note that
 * to do so, we need to ignore constrained DOFs by setting the residual in
 * these vector components to zero.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::get_error_residual(Errors &error_residual)
 *   {
 *     BlockVector<double> error_res(dofs_per_block);
 * 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *       if (!constraints.is_constrained(i))
 *         error_res(i) = system_rhs(i);
 * 
 *     error_residual.norm = error_res.l2_norm();
 *     error_residual.u    = error_res.block(u_dof).l2_norm();
 *     error_residual.p    = error_res.block(p_dof).l2_norm();
 *     error_residual.J    = error_res.block(J_dof).l2_norm();
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_update"></a> 
 * <h4>Solid::get_error_update</h4>
 * 

 * 
 * Determine the true Newton update error for the problem
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::get_error_update(const BlockVector<double> &newton_update,
 *                                     Errors &                   error_update)
 *   {
 *     BlockVector<double> error_ud(dofs_per_block);
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
 *       if (!constraints.is_constrained(i))
 *         error_ud(i) = newton_update(i);
 * 
 *     error_update.norm = error_ud.l2_norm();
 *     error_update.u    = error_ud.block(u_dof).l2_norm();
 *     error_update.p    = error_ud.block(p_dof).l2_norm();
 *     error_update.J    = error_ud.block(J_dof).l2_norm();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidget_total_solution"></a> 
 * <h4>Solid::get_total_solution</h4>
 * 

 * 
 * This function provides the total solution, which is valid at any Newton
 * step. This is required as, to reduce computational error, the total
 * solution is only updated at the end of the timestep.
 * 
 * @code
 *   template <int dim>
 *   BlockVector<double> Solid<dim>::get_total_solution(
 *     const BlockVector<double> &solution_delta) const
 *   {
 *     BlockVector<double> solution_total(solution_n);
 *     solution_total += solution_delta;
 *     return solution_total;
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidassemble_system"></a> 
 * <h4>Solid::assemble_system</h4>
 * 

 * 
 * Since we use TBB for assembly, we simply setup a copy of the
 * data structures required for the process and pass them, along
 * with the assembly functions to the WorkStream object for processing. Note
 * that we must ensure that the matrix and RHS vector are reset before any
 * assembly operations can occur. Furthermore, since we are describing a
 * problem with Neumann BCs, we will need the face normals and so must specify
 * this in the face update flags.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::assemble_system()
 *   {
 *     timer.enter_subsection("Assemble system");
 *     std::cout << " ASM_SYS " << std::flush;
 * 
 *     tangent_matrix = 0.0;
 *     system_rhs     = 0.0;
 * 
 *     const UpdateFlags uf_cell(update_values | update_gradients |
 *                               update_JxW_values);
 *     const UpdateFlags uf_face(update_values | update_normal_vectors |
 *                               update_JxW_values);
 * 
 *     PerTaskData_ASM per_task_data(dofs_per_cell);
 *     ScratchData_ASM scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);
 * 
 * @endcode
 * 
 * The syntax used here to pass data to the WorkStream class
 * is discussed in step-13.
 * 
 * @code
 *     WorkStream::run(
 *       dof_handler.active_cell_iterators(),
 *       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *              ScratchData_ASM &                                     scratch,
 *              PerTaskData_ASM &                                     data) {
 *         this->assemble_system_one_cell(cell, scratch, data);
 *       },
 *       [this](const PerTaskData_ASM &data) {
 *         this->constraints.distribute_local_to_global(data.cell_matrix,
 *                                                      data.cell_rhs,
 *                                                      data.local_dof_indices,
 *                                                      tangent_matrix,
 *                                                      system_rhs);
 *       },
 *       scratch_data,
 *       per_task_data);
 * 
 *     timer.leave_subsection();
 *   }
 * 
 * @endcode
 * 
 * Of course, we still have to define how we assemble the tangent matrix
 * contribution for a single cell.  We first need to reset and initialize some
 * of the scratch data structures and retrieve some basic information
 * regarding the DOF numbering on this cell.  We can precalculate the cell
 * shape function values and gradients. Note that the shape function gradients
 * are defined with regard to the current configuration.  That is
 * $\textrm{grad}\ \boldsymbol{\varphi} = \textrm{Grad}\ \boldsymbol{\varphi}
 * \ \mathbf{F}^{-1}$.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::assemble_system_one_cell(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     ScratchData_ASM &                                     scratch,
 *     PerTaskData_ASM &                                     data) const
 *   {
 *     data.reset();
 *     scratch.reset();
 *     scratch.fe_values.reinit(cell);
 *     cell->get_dof_indices(data.local_dof_indices);
 * 
 *     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
 *       quadrature_point_history.get_data(cell);
 *     Assert(lqph.size() == n_q_points, ExcInternalError());
 * 
 *     for (const unsigned int q_point :
 *          scratch.fe_values.quadrature_point_indices())
 *       {
 *         const Tensor<2, dim> F_inv = lqph[q_point]->get_F_inv();
 *         for (const unsigned int k : scratch.fe_values.dof_indices())
 *           {
 *             const unsigned int k_group = fe.system_to_base_index(k).first.first;
 * 
 *             if (k_group == u_dof)
 *               {
 *                 scratch.grad_Nx[q_point][k] =
 *                   scratch.fe_values[u_fe].gradient(k, q_point) * F_inv;
 *                 scratch.symm_grad_Nx[q_point][k] =
 *                   symmetrize(scratch.grad_Nx[q_point][k]);
 *               }
 *             else if (k_group == p_dof)
 *               scratch.Nx[q_point][k] =
 *                 scratch.fe_values[p_fe].value(k, q_point);
 *             else if (k_group == J_dof)
 *               scratch.Nx[q_point][k] =
 *                 scratch.fe_values[J_fe].value(k, q_point);
 *             else
 *               Assert(k_group <= J_dof, ExcInternalError());
 *           }
 *       }
 * 
 * @endcode
 * 
 * Now we build the local cell stiffness matrix and RHS vector. Since the
 * global and local system matrices are symmetric, we can exploit this
 * property by building only the lower half of the local matrix and copying
 * the values to the upper half.  So we only assemble half of the
 * $\mathsf{\mathbf{k}}_{uu}$, $\mathsf{\mathbf{k}}_{\widetilde{p}
 * \widetilde{p}} = \mathbf{0}$, $\mathsf{\mathbf{k}}_{\widetilde{J}
 * \widetilde{J}}$ blocks, while the whole
 * $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$,
 * $\mathsf{\mathbf{k}}_{u \widetilde{J}} = \mathbf{0}$,
 * $\mathsf{\mathbf{k}}_{u \widetilde{p}}$ blocks are built.
 *     

 * 
 * In doing so, we first extract some configuration dependent variables
 * from our quadrature history objects for the current quadrature point.
 * 
 * @code
 *     for (const unsigned int q_point :
 *          scratch.fe_values.quadrature_point_indices())
 *       {
 *         const SymmetricTensor<2, dim> tau     = lqph[q_point]->get_tau();
 *         const Tensor<2, dim>          tau_ns  = lqph[q_point]->get_tau();
 *         const SymmetricTensor<4, dim> Jc      = lqph[q_point]->get_Jc();
 *         const double                  det_F   = lqph[q_point]->get_det_F();
 *         const double                  p_tilde = lqph[q_point]->get_p_tilde();
 *         const double                  J_tilde = lqph[q_point]->get_J_tilde();
 *         const double dPsi_vol_dJ   = lqph[q_point]->get_dPsi_vol_dJ();
 *         const double d2Psi_vol_dJ2 = lqph[q_point]->get_d2Psi_vol_dJ2();
 *         const SymmetricTensor<2, dim> &I =
 *           Physics::Elasticity::StandardTensors<dim>::I;
 * 
 * @endcode
 * 
 * These two tensors store some precomputed data. Their use will
 * explained shortly.
 * 
 * @code
 *         SymmetricTensor<2, dim> symm_grad_Nx_i_x_Jc;
 *         Tensor<1, dim>          grad_Nx_i_comp_i_x_tau;
 * 
 * @endcode
 * 
 * Next we define some aliases to make the assembly process easier to
 * follow.
 * 
 * @code
 *         const std::vector<double> &                 N = scratch.Nx[q_point];
 *         const std::vector<SymmetricTensor<2, dim>> &symm_grad_Nx =
 *           scratch.symm_grad_Nx[q_point];
 *         const std::vector<Tensor<2, dim>> &grad_Nx = scratch.grad_Nx[q_point];
 *         const double                       JxW = scratch.fe_values.JxW(q_point);
 * 
 *         for (const unsigned int i : scratch.fe_values.dof_indices())
 *           {
 *             const unsigned int component_i =
 *               fe.system_to_component_index(i).first;
 *             const unsigned int i_group = fe.system_to_base_index(i).first.first;
 * 
 * @endcode
 * 
 * We first compute the contributions
 * from the internal forces.  Note, by
 * definition of the rhs as the negative
 * of the residual, these contributions
 * are subtracted.
 * 
 * @code
 *             if (i_group == u_dof)
 *               data.cell_rhs(i) -= (symm_grad_Nx[i] * tau) * JxW;
 *             else if (i_group == p_dof)
 *               data.cell_rhs(i) -= N[i] * (det_F - J_tilde) * JxW;
 *             else if (i_group == J_dof)
 *               data.cell_rhs(i) -= N[i] * (dPsi_vol_dJ - p_tilde) * JxW;
 *             else
 *               Assert(i_group <= J_dof, ExcInternalError());
 * 
 * @endcode
 * 
 * Before we go into the inner loop, we have one final chance to
 * introduce some optimizations. We've already taken into account
 * the symmetry of the system, and we can now precompute some
 * common terms that are repeatedly applied in the inner loop.
 * We won't be excessive here, but will rather focus on expensive
 * operations, namely those involving the rank-4 material stiffness
 * tensor and the rank-2 stress tensor.
 *             

 * 
 * What we may observe is that both of these tensors are contracted
 * with shape function gradients indexed on the "i" DoF. This
 * implies that this particular operation remains constant as we
 * loop over the "j" DoF. For that reason, we can extract this from
 * the inner loop and save the many operations that, for each
 * quadrature point and DoF index "i" and repeated over index "j"
 * are required to double contract a rank-2 symmetric tensor with a
 * rank-4 symmetric tensor, and a rank-1 tensor with a rank-2
 * tensor.
 *             

 * 
 * At the loss of some readability, this small change will reduce
 * the assembly time of the symmetrized system by about half when
 * using the simulation default parameters, and becomes more
 * significant as the h-refinement level increases.
 * 
 * @code
 *             if (i_group == u_dof)
 *               {
 *                 symm_grad_Nx_i_x_Jc    = symm_grad_Nx[i] * Jc;
 *                 grad_Nx_i_comp_i_x_tau = grad_Nx[i][component_i] * tau_ns;
 *               }
 * 
 * @endcode
 * 
 * Now we're prepared to compute the tangent matrix contributions:
 * 
 * @code
 *             for (const unsigned int j :
 *                  scratch.fe_values.dof_indices_ending_at(i))
 *               {
 *                 const unsigned int component_j =
 *                   fe.system_to_component_index(j).first;
 *                 const unsigned int j_group =
 *                   fe.system_to_base_index(j).first.first;
 * 
 * @endcode
 * 
 * This is the $\mathsf{\mathbf{k}}_{uu}$
 * contribution. It comprises a material contribution, and a
 * geometrical stress contribution which is only added along
 * the local matrix diagonals:
 * 
 * @code
 *                 if ((i_group == j_group) && (i_group == u_dof))
 *                   {
 * @endcode
 * 
 * The material contribution:
 * 
 * @code
 *                     data.cell_matrix(i, j) += symm_grad_Nx_i_x_Jc *  
 *                                               symm_grad_Nx[j] * JxW; 
 * 
 * @endcode
 * 
 * The geometrical stress contribution:
 * 
 * @code
 *                     if (component_i == component_j)
 *                       data.cell_matrix(i, j) +=
 *                         grad_Nx_i_comp_i_x_tau * grad_Nx[j][component_j] * JxW;
 *                   }
 * @endcode
 * 
 * Next is the $\mathsf{\mathbf{k}}_{ \widetilde{p} u}$
 * contribution
 * 
 * @code
 *                 else if ((i_group == p_dof) && (j_group == u_dof))
 *                   {
 *                     data.cell_matrix(i, j) += N[i] * det_F *               
 *                                               (symm_grad_Nx[j] * I) * JxW; 
 *                   }
 * @endcode
 * 
 * and lastly the $\mathsf{\mathbf{k}}_{ \widetilde{J}
 * \widetilde{p}}$ and $\mathsf{\mathbf{k}}_{ \widetilde{J}
 * \widetilde{J}}$ contributions:
 * 
 * @code
 *                 else if ((i_group == J_dof) && (j_group == p_dof))
 *                   data.cell_matrix(i, j) -= N[i] * N[j] * JxW;
 *                 else if ((i_group == j_group) && (i_group == J_dof))
 *                   data.cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW;
 *                 else
 *                   Assert((i_group <= J_dof) && (j_group <= J_dof),
 *                          ExcInternalError());
 *               }
 *           }
 *       }
 * 
 * @endcode
 * 
 * Next we assemble the Neumann contribution. We first check to see it the
 * cell face exists on a boundary on which a traction is applied and add
 * the contribution if this is the case.
 * 
 * @code
 *     for (const auto &face : cell->face_iterators())
 *       if (face->at_boundary() && face->boundary_id() == 6)
 *         {
 *           scratch.fe_face_values.reinit(cell, face);
 * 
 *           for (const unsigned int f_q_point :
 *                scratch.fe_face_values.quadrature_point_indices())
 *             {
 *               const Tensor<1, dim> &N =
 *                 scratch.fe_face_values.normal_vector(f_q_point);
 * 
 * @endcode
 * 
 * Using the face normal at this quadrature point we specify the
 * traction in reference configuration. For this problem, a
 * defined pressure is applied in the reference configuration.
 * The direction of the applied traction is assumed not to
 * evolve with the deformation of the domain. The traction is
 * defined using the first Piola-Kirchhoff stress is simply
 * $\mathbf{t} = \mathbf{P}\mathbf{N} = [p_0 \mathbf{I}]
 * \mathbf{N} = p_0 \mathbf{N}$ We use the time variable to
 * linearly ramp up the pressure load.
 *               

 * 
 * Note that the contributions to the right hand side vector we
 * compute here only exist in the displacement components of the
 * vector.
 * 
 * @code
 *               static const double p0 =
 *                 -4.0 / (parameters.scale * parameters.scale);
 *               const double         time_ramp = (time.current() / time.end());
 *               const double         pressure  = p0 * parameters.p_p0 * time_ramp;
 *               const Tensor<1, dim> traction  = pressure * N;
 * 
 *               for (const unsigned int i : scratch.fe_values.dof_indices())
 *                 {
 *                   const unsigned int i_group =
 *                     fe.system_to_base_index(i).first.first;
 * 
 *                   if (i_group == u_dof)
 *                     {
 *                       const unsigned int component_i =
 *                         fe.system_to_component_index(i).first;
 *                       const double Ni =
 *                         scratch.fe_face_values.shape_value(i, f_q_point);
 *                       const double JxW = scratch.fe_face_values.JxW(f_q_point);
 * 
 *                       data.cell_rhs(i) += (Ni * traction[component_i]) * JxW;
 *                     }
 *                 }
 *             }
 *         }
 * 
 * @endcode
 * 
 * Finally, we need to copy the lower half of the local matrix into the
 * upper half:
 * 
 * @code
 *     for (const unsigned int i : scratch.fe_values.dof_indices())
 *       for (const unsigned int j :
 *            scratch.fe_values.dof_indices_starting_at(i + 1))
 *         data.cell_matrix(i, j) = data.cell_matrix(j, i);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidmake_constraints"></a> 
 * <h4>Solid::make_constraints</h4>
 * The constraints for this problem are simple to describe.
 * In this particular example, the boundary values will be calculated for
 * the two first iterations of Newton's algorithm. In general, one would
 * build non-homogeneous constraints in the zeroth iteration (that is, when
 * `apply_dirichlet_bc == true` in the code block that follows) and build
 * only the corresponding homogeneous constraints in the following step. While
 * the current example has only homogeneous constraints, previous experiences
 * have shown that a common error is forgetting to add the extra condition
 * when refactoring the code to specific uses. This could lead to errors that
 * are hard to debug. In this spirit, we choose to make the code more verbose
 * in terms of what operations are performed at each Newton step.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::make_constraints(const int it_nr)
 *   {
 * @endcode
 * 
 * Since we (a) are dealing with an iterative Newton method, (b) are using
 * an incremental formulation for the displacement, and (c) apply the
 * constraints to the incremental displacement field, any non-homogeneous
 * constraints on the displacement update should only be specified at the
 * zeroth iteration. No subsequent contributions are to be made since the
 * constraints will be exactly satisfied after that iteration.
 * 
 * @code
 *     const bool apply_dirichlet_bc = (it_nr == 0);
 * 
 * @endcode
 * 
 * Furthermore, after the first Newton iteration within a timestep, the
 * constraints remain the same and we do not need to modify or rebuild them
 * so long as we do not clear the @p constraints object.
 * 
 * @code
 *     if (it_nr > 1)
 *       {
 *         std::cout << " --- " << std::flush;
 *         return;
 *       }
 * 
 *     std::cout << " CST " << std::flush;
 * 
 *     if (apply_dirichlet_bc)
 *       {
 * @endcode
 * 
 * At the zeroth Newton iteration we wish to apply the full set of
 * non-homogeneous and homogeneous constraints that represent the
 * boundary conditions on the displacement increment. Since in general
 * the constraints may be different at each time step, we need to clear
 * the constraints matrix and completely rebuild it. An example case
 * would be if a surface is accelerating; in such a scenario the change
 * in displacement is non-constant between each time step.
 * 
 * @code
 *         constraints.clear();
 * 
 * @endcode
 * 
 * The boundary conditions for the indentation problem in 3D are as
 * follows: On the -x, -y and -z faces (IDs 0,2,4) we set up a symmetry
 * condition to allow only planar movement while the +x and +z faces
 * (IDs 1,5) are traction free. In this contrived problem, part of the
 * +y face (ID 3) is set to have no motion in the x- and z-component.
 * Finally, as described earlier, the other part of the +y face has an
 * the applied pressure but is also constrained in the x- and
 * z-directions.
 *         

 * 
 * In the following, we will have to tell the function interpolation
 * boundary values which components of the solution vector should be
 * constrained (i.e., whether it's the x-, y-, z-displacements or
 * combinations thereof). This is done using ComponentMask objects (see
 * @ref GlossComponentMask) which we can get from the finite element if we
 * provide it with an extractor object for the component we wish to
 * select. To this end we first set up such extractor objects and later
 * use it when generating the relevant component masks:
 * 
 * @code
 *         const FEValuesExtractors::Scalar x_displacement(0);
 *         const FEValuesExtractors::Scalar y_displacement(1);
 * 
 *         {
 *           const int boundary_id = 0;
 * 
 *           VectorTools::interpolate_boundary_values(
 *             dof_handler,
 *             boundary_id,
 *             Functions::ZeroFunction<dim>(n_components),
 *             constraints,
 *             fe.component_mask(x_displacement));
 *         }
 *         {
 *           const int boundary_id = 2;
 * 
 *           VectorTools::interpolate_boundary_values(
 *             dof_handler,
 *             boundary_id,
 *             Functions::ZeroFunction<dim>(n_components),
 *             constraints,
 *             fe.component_mask(y_displacement));
 *         }
 * 
 *         if (dim == 3)
 *           {
 *             const FEValuesExtractors::Scalar z_displacement(2);
 * 
 *             {
 *               const int boundary_id = 3;
 * 
 *               VectorTools::interpolate_boundary_values(
 *                 dof_handler,
 *                 boundary_id,
 *                 Functions::ZeroFunction<dim>(n_components),
 *                 constraints,
 *                 (fe.component_mask(x_displacement) |
 *                  fe.component_mask(z_displacement)));
 *             }
 *             {
 *               const int boundary_id = 4;
 * 
 *               VectorTools::interpolate_boundary_values(
 *                 dof_handler,
 *                 boundary_id,
 *                 Functions::ZeroFunction<dim>(n_components),
 *                 constraints,
 *                 fe.component_mask(z_displacement));
 *             }
 * 
 *             {
 *               const int boundary_id = 6;
 * 
 *               VectorTools::interpolate_boundary_values(
 *                 dof_handler,
 *                 boundary_id,
 *                 Functions::ZeroFunction<dim>(n_components),
 *                 constraints,
 *                 (fe.component_mask(x_displacement) |
 *                  fe.component_mask(z_displacement)));
 *             }
 *           }
 *         else
 *           {
 *             {
 *               const int boundary_id = 3;
 * 
 *               VectorTools::interpolate_boundary_values(
 *                 dof_handler,
 *                 boundary_id,
 *                 Functions::ZeroFunction<dim>(n_components),
 *                 constraints,
 *                 (fe.component_mask(x_displacement)));
 *             }
 *             {
 *               const int boundary_id = 6;
 * 
 *               VectorTools::interpolate_boundary_values(
 *                 dof_handler,
 *                 boundary_id,
 *                 Functions::ZeroFunction<dim>(n_components),
 *                 constraints,
 *                 (fe.component_mask(x_displacement)));
 *             }
 *           }
 *       }
 *     else
 *       {
 * @endcode
 * 
 * As all Dirichlet constraints are fulfilled exactly after the zeroth
 * Newton iteration, we want to ensure that no further modification are
 * made to those entries. This implies that we want to convert
 * all non-homogeneous Dirichlet constraints into homogeneous ones.
 *         

 * 
 * In this example the procedure to do this is quite straightforward,
 * and in fact we can (and will) circumvent any unnecessary operations
 * when only homogeneous boundary conditions are applied.
 * In a more general problem one should be mindful of hanging node
 * and periodic constraints, which may also introduce some
 * inhomogeneities. It might then be advantageous to keep disparate
 * objects for the different types of constraints, and merge them
 * together once the homogeneous Dirichlet constraints have been
 * constructed.
 * 
 * @code
 *         if (constraints.has_inhomogeneities())
 *           {
 * @endcode
 * 
 * Since the affine constraints were finalized at the previous
 * Newton iteration, they may not be modified directly. So
 * we need to copy them to another temporary object and make
 * modification there. Once we're done, we'll transfer them
 * back to the main @p constraints object.
 * 
 * @code
 *             AffineConstraints<double> homogeneous_constraints(constraints);
 *             for (unsigned int dof = 0; dof != dof_handler.n_dofs(); ++dof)
 *               if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
 *                 homogeneous_constraints.set_inhomogeneity(dof, 0.0);
 * 
 *             constraints.clear();
 *             constraints.copy_from(homogeneous_constraints);
 *           }
 *       }
 * 
 *     constraints.close();
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solidassemble_sc"></a> 
 * <h4>Solid::assemble_sc</h4>
 * Solving the entire block system is a bit problematic as there are no
 * contributions to the $\mathsf{\mathbf{K}}_{ \widetilde{J} \widetilde{J}}$
 * block, rendering it noninvertible (when using an iterative solver).
 * Since the pressure and dilatation variables DOFs are discontinuous, we can
 * condense them out to form a smaller displacement-only system which
 * we will then solve and subsequently post-process to retrieve the
 * pressure and dilatation solutions.
 * 

 * 
 * The static condensation process could be performed at a global level but we
 * need the inverse of one of the blocks. However, since the pressure and
 * dilatation variables are discontinuous, the static condensation (SC)
 * operation can also be done on a per-cell basis and we can produce the
 * inverse of the block-diagonal
 * $\mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}$ block by inverting the
 * local blocks. We can again use TBB to do this since each operation will be
 * independent of one another.
 *   

 * 
 * Using the TBB via the WorkStream class, we assemble the contributions to
 * form
 * $
 * \mathsf{\mathbf{K}}_{\textrm{con}}
 * = \bigl[ \mathsf{\mathbf{K}}_{uu} +
 * \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]
 * $
 * from each element's contributions. These contributions are then added to
 * the global stiffness matrix. Given this description, the following two
 * functions should be clear:
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::assemble_sc()
 *   {
 *     timer.enter_subsection("Perform static condensation");
 *     std::cout << " ASM_SC " << std::flush;
 * 
 *     PerTaskData_SC per_task_data(dofs_per_cell,
 *                                  element_indices_u.size(),
 *                                  element_indices_p.size(),
 *                                  element_indices_J.size());
 *     ScratchData_SC scratch_data;
 * 
 *     WorkStream::run(dof_handler.active_cell_iterators(),
 *                     *this,
 *                     &Solid::assemble_sc_one_cell,
 *                     &Solid::copy_local_to_global_sc,
 *                     scratch_data,
 *                     per_task_data);
 * 
 *     timer.leave_subsection();
 *   }
 * 
 * 
 *   template <int dim>
 *   void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC &data)
 *   {
 *     for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *       for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *         tangent_matrix.add(data.local_dof_indices[i],
 *                            data.local_dof_indices[j],
 *                            data.cell_matrix(i, j));
 *   }
 * 
 * 
 * @endcode
 * 
 * Now we describe the static condensation process. As per usual, we must
 * first find out which global numbers the degrees of freedom on this cell
 * have and reset some data structures:
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::assemble_sc_one_cell(
 *     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     ScratchData_SC &                                      scratch,
 *     PerTaskData_SC &                                      data)
 *   {
 *     data.reset();
 *     scratch.reset();
 *     cell->get_dof_indices(data.local_dof_indices);
 * 
 * @endcode
 * 
 * We now extract the contribution of the dofs associated with the current
 * cell to the global stiffness matrix.  The discontinuous nature of the
 * $\widetilde{p}$ and $\widetilde{J}$ interpolations mean that their is
 * no coupling of the local contributions at the global level. This is not
 * the case with the $\mathbf{u}$ dof.  In other words,
 * $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$,
 * $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{p}}$ and
 * $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$, when extracted
 * from the global stiffness matrix are the element contributions.  This
 * is not the case for $\mathsf{\mathbf{k}}_{uu}$.
 *     

 * 
 * Note: A lower-case symbol is used to denote element stiffness matrices.
 * 

 * 
 * Currently the matrix corresponding to
 * the dof associated with the current element
 * (denoted somewhat loosely as $\mathsf{\mathbf{k}}$)
 * is of the form:
 * @f{align*}
 * \begin{bmatrix}
 * \mathsf{\mathbf{k}}_{uu}  &  \mathsf{\mathbf{k}}_{u\widetilde{p}}
 * & \mathbf{0}
 * \\ \mathsf{\mathbf{k}}_{\widetilde{p}u} & \mathbf{0}  &
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}
 * \\ \mathbf{0}  &  \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}}  &
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}} \end{bmatrix}
 * @f}
 *     

 * 
 * We now need to modify it such that it appear as
 * @f{align*}
 * \begin{bmatrix}
 * \mathsf{\mathbf{k}}_{\textrm{con}}   &
 * \mathsf{\mathbf{k}}_{u\widetilde{p}}    & \mathbf{0}
 * \\ \mathsf{\mathbf{k}}_{\widetilde{p}u} & \mathbf{0} &
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \\ \mathbf{0} & \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}} &
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}} \end{bmatrix}
 * @f}
 * with $\mathsf{\mathbf{k}}_{\textrm{con}} = \bigl[
 * \mathsf{\mathbf{k}}_{uu} +\overline{\overline{\mathsf{\mathbf{k}}}}~
 * \bigr]$ where $               \overline{\overline{\mathsf{\mathbf{k}}}}
 * \dealcoloneq \mathsf{\mathbf{k}}_{u\widetilde{p}}
 * \overline{\mathsf{\mathbf{k}}} \mathsf{\mathbf{k}}_{\widetilde{p}u}
 * $
 * and
 * $
 * \overline{\mathsf{\mathbf{k}}} =
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
 * $.
 *     

 * 
 * At this point, we need to take note of
 * the fact that global data already exists
 * in the $\mathsf{\mathbf{K}}_{uu}$,
 * $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$
 * and
 * $\mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{p}}$
 * sub-blocks.  So if we are to modify them, we must account for the data
 * that is already there (i.e. simply add to it or remove it if
 * necessary).  Since the copy_local_to_global operation is a "+="
 * operation, we need to take this into account
 *     

 * 
 * For the $\mathsf{\mathbf{K}}_{uu}$ block in particular, this means that
 * contributions have been added from the surrounding cells, so we need to
 * be careful when we manipulate this block.  We can't just erase the
 * sub-blocks.
 *     

 * 
 * This is the strategy we will employ to get the sub-blocks we want:
 *     

 * 
 * - $ {\mathsf{\mathbf{k}}}_{\textrm{store}}$:
 * Since we don't have access to $\mathsf{\mathbf{k}}_{uu}$,
 * but we know its contribution is added to
 * the global $\mathsf{\mathbf{K}}_{uu}$ matrix, we just want
 * to add the element wise
 * static-condensation $\overline{\overline{\mathsf{\mathbf{k}}}}$.
 *     

 * 
 * - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$:
 * Similarly, $\mathsf{\mathbf{k}}_{\widetilde{p}
 * \widetilde{J}}$ exists in
 * the subblock. Since the copy
 * operation is a += operation, we
 * need to subtract the existing
 * $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
 * submatrix in addition to
 * "adding" that which we wish to
 * replace it with.
 *     

 * 
 * - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}$:
 * Since the global matrix
 * is symmetric, this block is the
 * same as the one above and we
 * can simply use
 * $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$
 * as a substitute for this one.
 *     

 * 
 * We first extract element data from the
 * system matrix. So first we get the
 * entire subblock for the cell, then
 * extract $\mathsf{\mathbf{k}}$
 * for the dofs associated with
 * the current element
 * 
 * @code
 *     data.k_orig.extract_submatrix_from(tangent_matrix,
 *                                        data.local_dof_indices,
 *                                        data.local_dof_indices);
 * @endcode
 * 
 * and next the local matrices for
 * $\mathsf{\mathbf{k}}_{ \widetilde{p} u}$
 * $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
 * and
 * $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$:
 * 
 * @code
 *     data.k_pu.extract_submatrix_from(data.k_orig,
 *                                      element_indices_p,
 *                                      element_indices_u);
 *     data.k_pJ.extract_submatrix_from(data.k_orig,
 *                                      element_indices_p,
 *                                      element_indices_J);
 *     data.k_JJ.extract_submatrix_from(data.k_orig,
 *                                      element_indices_J,
 *                                      element_indices_J);
 * 
 * @endcode
 * 
 * To get the inverse of $\mathsf{\mathbf{k}}_{\widetilde{p}
 * \widetilde{J}}$, we invert it directly.  This operation is relatively
 * inexpensive since $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
 * since block-diagonal.
 * 
 * @code
 *     data.k_pJ_inv.invert(data.k_pJ);
 * 
 * @endcode
 * 
 * Now we can make condensation terms to
 * add to the $\mathsf{\mathbf{k}}_{uu}$
 * block and put them in
 * the cell local matrix
 * $
 * \mathsf{\mathbf{A}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $:
 * 
 * @code
 *     data.k_pJ_inv.mmult(data.A, data.k_pu);
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{B}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $
 * 
 * @code
 *     data.k_JJ.mmult(data.B, data.A);
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{C}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $
 * 
 * @code
 *     data.k_pJ_inv.Tmmult(data.C, data.B);
 * @endcode
 * 
 * $
 * \overline{\overline{\mathsf{\mathbf{k}}}}
 * =
 * \mathsf{\mathbf{k}}_{u \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $
 * 
 * @code
 *     data.k_pu.Tmmult(data.k_bbar, data.C);
 *     data.k_bbar.scatter_matrix_to(element_indices_u,
 *                                   element_indices_u,
 *                                   data.cell_matrix);
 * 
 * @endcode
 * 
 * Next we place
 * $\mathsf{\mathbf{k}}^{-1}_{ \widetilde{p} \widetilde{J}}$
 * in the
 * $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
 * block for post-processing.  Note again
 * that we need to remove the
 * contribution that already exists there.
 * 
 * @code
 *     data.k_pJ_inv.add(-1.0, data.k_pJ);
 *     data.k_pJ_inv.scatter_matrix_to(element_indices_p,
 *                                     element_indices_J,
 *                                     data.cell_matrix);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solidsolve_linear_system"></a> 
 * <h4>Solid::solve_linear_system</h4>
 * We now have all of the necessary components to use one of two possible
 * methods to solve the linearised system. The first is to perform static
 * condensation on an element level, which requires some alterations
 * to the tangent matrix and RHS vector. Alternatively, the full block
 * system can be solved by performing condensation on a global level.
 * Below we implement both approaches.
 * 
 * @code
 *   template <int dim>
 *   std::pair<unsigned int, double>
 *   Solid<dim>::solve_linear_system(BlockVector<double> &newton_update)
 *   {
 *     unsigned int lin_it  = 0;
 *     double       lin_res = 0.0;
 * 
 *     if (parameters.use_static_condensation == true)
 *       {
 * @endcode
 * 
 * Firstly, here is the approach using the (permanent) augmentation of
 * the tangent matrix. For the following, recall that
 * @f{align*}
 * \mathsf{\mathbf{K}}_{\textrm{store}}
 * \dealcoloneq
 * \begin{bmatrix}
 * \mathsf{\mathbf{K}}_{\textrm{con}}      &
 * \mathsf{\mathbf{K}}_{u\widetilde{p}}    & \mathbf{0}
 * \\  \mathsf{\mathbf{K}}_{\widetilde{p}u}    &       \mathbf{0} &
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \\  \mathbf{0}      &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}                &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}} \end{bmatrix} \, .
 * @f}
 * and
 * @f{align*}
 * d \widetilde{\mathsf{\mathbf{p}}}
 * & =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}} \bigr]
 * \\ d \widetilde{\mathsf{\mathbf{J}}}
 * & =
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d
 * \mathsf{\mathbf{u}} \bigr]
 * \\ \Rightarrow d \widetilde{\mathsf{\mathbf{p}}}
 * &= \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \underbrace{\bigl[\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathsf{\mathbf{K}}}}\bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d
 * \mathsf{\mathbf{u}} \bigr]
 * @f}
 * and thus
 * @f[
 * \underbrace{\bigl[ \mathsf{\mathbf{K}}_{uu} +
 * \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]
 * }_{\mathsf{\mathbf{K}}_{\textrm{con}}} d
 * \mathsf{\mathbf{u}}
 * =
 * \underbrace{
 * \Bigl[
 * \mathsf{\mathbf{F}}_{u}
 * - \mathsf{\mathbf{K}}_{u\widetilde{p}} \bigl[
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \overline{\mathsf{\mathbf{K}}}\mathsf{\mathbf{F}}_{\widetilde{p}}
 * \bigr]
 * \Bigr]}_{\mathsf{\mathbf{F}}_{\textrm{con}}}
 * @f]
 * where
 * @f[
 * \overline{\overline{\mathsf{\mathbf{K}}}} \dealcoloneq
 * \mathsf{\mathbf{K}}_{u\widetilde{p}}
 * \overline{\mathsf{\mathbf{K}}}
 * \mathsf{\mathbf{K}}_{\widetilde{p}u} \, .
 * @f]
 * 

 * 
 * At the top, we allocate two temporary vectors to help with the
 * static condensation, and variables to store the number of
 * linear solver iterations and the (hopefully converged) residual.
 * 
 * @code
 *         BlockVector<double> A(dofs_per_block);
 *         BlockVector<double> B(dofs_per_block);
 * 
 * 
 * @endcode
 * 
 * In the first step of this function, we solve for the incremental
 * displacement $d\mathbf{u}$.  To this end, we perform static
 * condensation to make
 * $\mathsf{\mathbf{K}}_{\textrm{con}}
 * = \bigl[ \mathsf{\mathbf{K}}_{uu} +
 * \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]$
 * and put
 * $\mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}$
 * in the original $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$
 * block. That is, we make $\mathsf{\mathbf{K}}_{\textrm{store}}$.
 * 
 * @code
 *         {
 *           assemble_sc();
 * 
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof)
 *             .vmult(A.block(J_dof), system_rhs.block(p_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{B}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $
 * 
 * @code
 *           tangent_matrix.block(J_dof, J_dof)
 *             .vmult(B.block(J_dof), A.block(J_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $
 * 
 * @code
 *           A.block(J_dof) = system_rhs.block(J_dof);
 *           A.block(J_dof) -= B.block(J_dof);
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof)
 *             .Tvmult(A.block(p_dof), A.block(J_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{u}
 * =
 * \mathsf{\mathbf{K}}_{u \widetilde{p}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $
 * 
 * @code
 *           tangent_matrix.block(u_dof, p_dof)
 *             .vmult(A.block(u_dof), A.block(p_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{F}}_{\text{con}}
 * =
 * \mathsf{\mathbf{F}}_{u}
 * -
 * \mathsf{\mathbf{K}}_{u \widetilde{p}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $
 * 
 * @code
 *           system_rhs.block(u_dof) -= A.block(u_dof);
 * 
 *           timer.enter_subsection("Linear solver");
 *           std::cout << " SLV " << std::flush;
 *           if (parameters.type_lin == "CG")
 *             {
 *               const auto solver_its = static_cast<unsigned int>(
 *                 tangent_matrix.block(u_dof, u_dof).m() *
 *                 parameters.max_iterations_lin);
 *               const double tol_sol =
 *                 parameters.tol_lin * system_rhs.block(u_dof).l2_norm();
 * 
 *               SolverControl solver_control(solver_its, tol_sol);
 * 
 *               GrowingVectorMemory<Vector<double>> GVM;
 *               SolverCG<Vector<double>> solver_CG(solver_control, GVM);
 * 
 * @endcode
 * 
 * We've chosen by default a SSOR preconditioner as it appears to
 * provide the fastest solver convergence characteristics for this
 * problem on a single-thread machine.  However, this might not be
 * true for different problem sizes.
 * 
 * @code
 *               PreconditionSelector<SparseMatrix<double>, Vector<double>>
 *                 preconditioner(parameters.preconditioner_type,
 *                                parameters.preconditioner_relaxation);
 *               preconditioner.use_matrix(tangent_matrix.block(u_dof, u_dof));
 * 
 *               solver_CG.solve(tangent_matrix.block(u_dof, u_dof),
 *                               newton_update.block(u_dof),
 *                               system_rhs.block(u_dof),
 *                               preconditioner);
 * 
 *               lin_it  = solver_control.last_step();
 *               lin_res = solver_control.last_value();
 *             }
 *           else if (parameters.type_lin == "Direct")
 *             {
 * @endcode
 * 
 * Otherwise if the problem is small
 * enough, a direct solver can be
 * utilised.
 * 
 * @code
 *               SparseDirectUMFPACK A_direct;
 *               A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
 *               A_direct.vmult(newton_update.block(u_dof),
 *                              system_rhs.block(u_dof));
 * 
 *               lin_it  = 1;
 *               lin_res = 0.0;
 *             }
 *           else
 *             Assert(false, ExcMessage("Linear solver type not implemented"));
 * 
 *           timer.leave_subsection();
 *         }
 * 
 * @endcode
 * 
 * Now that we have the displacement update, distribute the constraints
 * back to the Newton update:
 * 
 * @code
 *         constraints.distribute(newton_update);
 * 
 *         timer.enter_subsection("Linear solver postprocessing");
 *         std::cout << " PP " << std::flush;
 * 
 * @endcode
 * 
 * The next step after solving the displacement
 * problem is to post-process to get the
 * dilatation solution from the
 * substitution:
 * $
 * d \widetilde{\mathsf{\mathbf{J}}}
 * = \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * \bigr]
 * $
 * 
 * @code
 *         {
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $
 * 
 * @code
 *           tangent_matrix.block(p_dof, u_dof)
 *             .vmult(A.block(p_dof), newton_update.block(u_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $
 * 
 * @code
 *           A.block(p_dof) *= -1.0;
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $
 * 
 * @code
 *           A.block(p_dof) += system_rhs.block(p_dof);
 * @endcode
 * 
 * $
 * d\mathsf{\mathbf{\widetilde{J}}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p}\widetilde{J}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * ]
 * $
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof)
 *             .vmult(newton_update.block(J_dof), A.block(p_dof));
 *         }
 * 
 * @endcode
 * 
 * we ensure here that any Dirichlet constraints
 * are distributed on the updated solution:
 * 
 * @code
 *         constraints.distribute(newton_update);
 * 
 * @endcode
 * 
 * Finally we solve for the pressure
 * update with the substitution:
 * $
 * d \widetilde{\mathsf{\mathbf{p}}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * - \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * \bigr]
 * $
 * 
 * @code
 *         {
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $
 * 
 * @code
 *           tangent_matrix.block(J_dof, J_dof)
 *             .vmult(A.block(J_dof), newton_update.block(J_dof));
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * -\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $
 * 
 * @code
 *           A.block(J_dof) *= -1.0;
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $
 * 
 * @code
 *           A.block(J_dof) += system_rhs.block(J_dof);
 * @endcode
 * 
 * and finally....
 * $
 * d \widetilde{\mathsf{\mathbf{p}}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * - \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * \bigr]
 * $
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof)
 *             .Tvmult(newton_update.block(p_dof), A.block(J_dof));
 *         }
 * 
 * @endcode
 * 
 * We are now at the end, so we distribute all
 * constrained dofs back to the Newton
 * update:
 * 
 * @code
 *         constraints.distribute(newton_update);
 * 
 *         timer.leave_subsection();
 *       }
 *     else
 *       {
 *         std::cout << " ------ " << std::flush;
 * 
 *         timer.enter_subsection("Linear solver");
 *         std::cout << " SLV " << std::flush;
 * 
 *         if (parameters.type_lin == "CG")
 *           {
 * @endcode
 * 
 * Manual condensation of the dilatation and pressure fields on
 * a local level, and subsequent post-processing, took quite a
 * bit of effort to achieve. To recap, we had to produce the
 * inverse matrix
 * $\mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}$, which
 * was permanently written into the global tangent matrix. We then
 * permanently modified $\mathsf{\mathbf{K}}_{uu}$ to produce
 * $\mathsf{\mathbf{K}}_{\textrm{con}}$. This involved the
 * extraction and manipulation of local sub-blocks of the tangent
 * matrix. After solving for the displacement, the individual
 * matrix-vector operations required to solve for dilatation and
 * pressure were carefully implemented. Contrast these many sequence
 * of steps to the much simpler and transparent implementation using
 * functionality provided by the LinearOperator class.
 * 

 * 
 * For ease of later use, we define some aliases for
 * blocks in the RHS vector
 * 
 * @code
 *             const Vector<double> &f_u = system_rhs.block(u_dof);
 *             const Vector<double> &f_p = system_rhs.block(p_dof);
 *             const Vector<double> &f_J = system_rhs.block(J_dof);
 * 
 * @endcode
 * 
 * ... and for blocks in the Newton update vector.
 * 
 * @code
 *             Vector<double> &d_u = newton_update.block(u_dof);
 *             Vector<double> &d_p = newton_update.block(p_dof);
 *             Vector<double> &d_J = newton_update.block(J_dof);
 * 
 * @endcode
 * 
 * We next define some linear operators for the tangent matrix
 * sub-blocks We will exploit the symmetry of the system, so not all
 * blocks are required.
 * 
 * @code
 *             const auto K_uu =
 *               linear_operator(tangent_matrix.block(u_dof, u_dof));
 *             const auto K_up =
 *               linear_operator(tangent_matrix.block(u_dof, p_dof));
 *             const auto K_pu =
 *               linear_operator(tangent_matrix.block(p_dof, u_dof));
 *             const auto K_Jp =
 *               linear_operator(tangent_matrix.block(J_dof, p_dof));
 *             const auto K_JJ =
 *               linear_operator(tangent_matrix.block(J_dof, J_dof));
 * 
 * @endcode
 * 
 * We then construct a LinearOperator that represents the inverse of
 * (square block)
 * $\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}$. Since it is
 * diagonal (or, when a higher order ansatz it used, nearly
 * diagonal), a Jacobi preconditioner is suitable.
 * 
 * @code
 *             PreconditionSelector<SparseMatrix<double>, Vector<double>>
 *               preconditioner_K_Jp_inv("jacobi");
 *             preconditioner_K_Jp_inv.use_matrix(
 *               tangent_matrix.block(J_dof, p_dof));
 *             ReductionControl solver_control_K_Jp_inv(
 *               static_cast<unsigned int>(tangent_matrix.block(J_dof, p_dof).m() *
 *                                         parameters.max_iterations_lin),
 *               1.0e-30,
 *               parameters.tol_lin);
 *             SolverSelector<Vector<double>> solver_K_Jp_inv;
 *             solver_K_Jp_inv.select("cg");
 *             solver_K_Jp_inv.set_control(solver_control_K_Jp_inv);
 *             const auto K_Jp_inv =
 *               inverse_operator(K_Jp, solver_K_Jp_inv, preconditioner_K_Jp_inv);
 * 
 * @endcode
 * 
 * Now we can construct that transpose of
 * $\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}$ and a
 * linear operator that represents the condensed operations
 * $\overline{\mathsf{\mathbf{K}}}$ and
 * $\overline{\overline{\mathsf{\mathbf{K}}}}$ and the final
 * augmented matrix
 * $\mathsf{\mathbf{K}}_{\textrm{con}}$.
 * Note that the schur_complement() operator could also be of use
 * here, but for clarity and the purpose of demonstrating the
 * similarities between the formulation and implementation of the
 * linear solution scheme, we will perform these operations
 * manually.
 * 
 * @code
 *             const auto K_pJ_inv     = transpose_operator(K_Jp_inv);
 *             const auto K_pp_bar     = K_Jp_inv * K_JJ * K_pJ_inv;
 *             const auto K_uu_bar_bar = K_up * K_pp_bar * K_pu;
 *             const auto K_uu_con     = K_uu + K_uu_bar_bar;
 * 
 * @endcode
 * 
 * Lastly, we define an operator for inverse of augmented stiffness
 * matrix, namely $\mathsf{\mathbf{K}}_{\textrm{con}}^{-1}$. Note
 * that the preconditioner for the augmented stiffness matrix is
 * different to the case when we use static condensation. In this
 * instance, the preconditioner is based on a non-modified
 * $\mathsf{\mathbf{K}}_{uu}$, while with the first approach we
 * actually modified the entries of this sub-block. However, since
 * $\mathsf{\mathbf{K}}_{\textrm{con}}$ and
 * $\mathsf{\mathbf{K}}_{uu}$ operate on the same space, it remains
 * adequate for this problem.
 * 
 * @code
 *             PreconditionSelector<SparseMatrix<double>, Vector<double>>
 *               preconditioner_K_con_inv(parameters.preconditioner_type,
 *                                        parameters.preconditioner_relaxation);
 *             preconditioner_K_con_inv.use_matrix(
 *               tangent_matrix.block(u_dof, u_dof));
 *             ReductionControl solver_control_K_con_inv(
 *               static_cast<unsigned int>(tangent_matrix.block(u_dof, u_dof).m() *
 *                                         parameters.max_iterations_lin),
 *               1.0e-30,
 *               parameters.tol_lin);
 *             SolverSelector<Vector<double>> solver_K_con_inv;
 *             solver_K_con_inv.select("cg");
 *             solver_K_con_inv.set_control(solver_control_K_con_inv);
 *             const auto K_uu_con_inv =
 *               inverse_operator(K_uu_con,
 *                                solver_K_con_inv,
 *                                preconditioner_K_con_inv);
 * 
 * @endcode
 * 
 * Now we are in a position to solve for the displacement field.
 * We can nest the linear operations, and the result is immediately
 * written to the Newton update vector.
 * It is clear that the implementation closely mimics the derivation
 * stated in the introduction.
 * 
 * @code
 *             d_u =
 *               K_uu_con_inv * (f_u - K_up * (K_Jp_inv * f_J - K_pp_bar * f_p));
 * 
 *             timer.leave_subsection();
 * 
 * @endcode
 * 
 * The operations need to post-process for the dilatation and
 * pressure fields are just as easy to express.
 * 
 * @code
 *             timer.enter_subsection("Linear solver postprocessing");
 *             std::cout << " PP " << std::flush;
 * 
 *             d_J = K_pJ_inv * (f_p - K_pu * d_u);
 *             d_p = K_Jp_inv * (f_J - K_JJ * d_J);
 * 
 *             lin_it  = solver_control_K_con_inv.last_step();
 *             lin_res = solver_control_K_con_inv.last_value();
 *           }
 *         else if (parameters.type_lin == "Direct")
 *           {
 * @endcode
 * 
 * Solve the full block system with
 * a direct solver. As it is relatively
 * robust, it may be immune to problem
 * arising from the presence of the zero
 * $\mathsf{\mathbf{K}}_{ \widetilde{J} \widetilde{J}}$
 * block.
 * 
 * @code
 *             SparseDirectUMFPACK A_direct;
 *             A_direct.initialize(tangent_matrix);
 *             A_direct.vmult(newton_update, system_rhs);
 * 
 *             lin_it  = 1;
 *             lin_res = 0.0;
 * 
 *             std::cout << " -- " << std::flush;
 *           }
 *         else
 *           Assert(false, ExcMessage("Linear solver type not implemented"));
 * 
 *         timer.leave_subsection();
 * 
 * @endcode
 * 
 * Finally, we again ensure here that any Dirichlet
 * constraints are distributed on the updated solution:
 * 
 * @code
 *         constraints.distribute(newton_update);
 *       }
 * 
 *     return std::make_pair(lin_it, lin_res);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="Solidoutput_results"></a> 
 * <h4>Solid::output_results</h4>
 * Here we present how the results are written to file to be viewed
 * using ParaView or VisIt. The method is similar to that shown in previous
 * tutorials so will not be discussed in detail.
 * 
 * @code
 *   template <int dim>
 *   void Solid<dim>::output_results() const
 *   {
 *     DataOut<dim> data_out;
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *       data_component_interpretation(
 *         dim, DataComponentInterpretation::component_is_part_of_vector);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 *     data_component_interpretation.push_back(
 *       DataComponentInterpretation::component_is_scalar);
 * 
 *     std::vector<std::string> solution_name(dim, "displacement");
 *     solution_name.emplace_back("pressure");
 *     solution_name.emplace_back("dilatation");
 * 
 *     DataOutBase::VtkFlags output_flags;
 *     output_flags.write_higher_order_cells = true;
 *     data_out.set_flags(output_flags);
 * 
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution_n,
 *                              solution_name,
 *                              DataOut<dim>::type_dof_data,
 *                              data_component_interpretation);
 * 
 * @endcode
 * 
 * Since we are dealing with a large deformation problem, it would be nice
 * to display the result on a displaced grid!  The MappingQEulerian class
 * linked with the DataOut class provides an interface through which this
 * can be achieved without physically moving the grid points in the
 * Triangulation object ourselves.  We first need to copy the solution to
 * a temporary vector and then create the Eulerian mapping. We also
 * specify the polynomial degree to the DataOut object in order to produce
 * a more refined output data set when higher order polynomials are used.
 * 
 * @code
 *     Vector<double> soln(solution_n.size());
 *     for (unsigned int i = 0; i < soln.size(); ++i)
 *       soln(i) = solution_n(i);
 *     MappingQEulerian<dim> q_mapping(degree, dof_handler, soln);
 *     data_out.build_patches(q_mapping, degree);
 * 
 *     std::ofstream output("solution-" + std::to_string(dim) + "d-" +
 *                          std::to_string(time.get_timestep()) + ".vtu");
 *     data_out.write_vtu(output);
 *   }
 * 
 * } // namespace Step44
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Mainfunction"></a> 
 * <h3>Main function</h3>
 * Lastly we provide the main driver function which appears
 * no different to the other tutorials.
 * 
 * @code
 * int main()
 * {
 *   using namespace Step44;
 * 
 *   try
 *     {
 *       const unsigned int dim = 3;
 *       Solid<dim>         solid("parameters.prm");
 *       solid.run();
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
examples/step-44/doc/results.dox



<a name="Results"></a><h1>Results</h1>


首先，我们提出了一系列3维结果与文献中的结果的比较（见Reese等人(2000)），以证明该程序按预期工作。

我们首先比较了 $Q_1-DGPM_0-DGPM_0$ 和 $Q_2-DGPM_1-DGPM_1$ 公式的网格细化的收敛性，如下图所总结的。块的上表面的中点的垂直位移被用来评估收敛性。对于不同的载荷参数 $p/p_0$ 值，两种方案都表现出良好的收敛特性。这些结果与文献中的结果一致。低阶公式通常高估了低层次细化的位移，而高阶插值方案则低估了位移，但程度较轻。这个基准，以及其他一系列没有在这里显示的基准，使我们相信代码在正常工作。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_convergence.png" alt="">
	<p align="center">
        Convergence of the $Q_1-DGPM_0-DGPM_0$ formulation in 3-d.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_convergence.png" alt="">
	<p align="center">
        Convergence of the $Q_2-DGPM_1-DGPM_1$ formulation in 3-d.
	</p>
    </td>
  </tr>
</table> 


下面是运行该问题产生的典型屏幕输出。所展示的特殊情况是 $Q_2-DGPM_1-DGPM_1$ 公式的情况。很明显，使用Newton-Raphson方法，可以得到二次收敛的解决方案。在所有的时间步长中，解的收敛是在5个牛顿增量内实现的。收敛后的位移的 $L_2$ -norm比几何尺度小几个数量级。

@code
Grid:
	 Reference volume: 1e-09
Triangulation:
	 Number of active cells: 64
	 Number of degrees of freedom: 2699
    Setting up quadrature point data...


Timestep 1 @ 0.1s
___________________________________________________________________________________________________________________________________________________________
                 SOLVER STEP                   |  LIN_IT   LIN_RES    RES_NORM     RES_U     RES_P      RES_J     NU_NORM      NU_U       NU_P       NU_J
___________________________________________________________________________________________________________________________________________________________
  0  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     786  2.118e-06  1.000e+00  1.000e+00  0.000e+00  0.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00
  1  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     552  1.031e-03  8.563e-02  8.563e-02  9.200e-13  3.929e-08  1.060e-01  3.816e-02  1.060e-01  1.060e-01
  2  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     667  5.602e-06  2.482e-03  2.482e-03  3.373e-15  2.982e-10  2.936e-03  2.053e-04  2.936e-03  2.936e-03
  3  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     856  6.469e-10  2.129e-06  2.129e-06  2.245e-19  1.244e-13  1.887e-06  7.289e-07  1.887e-06  1.887e-06
  4  ASM_R  CONVERGED!
___________________________________________________________________________________________________________________________________________________________
Relative errors:
Displacement:	7.289e-07
Force: 		2.451e-10
Dilatation:	1.353e-07
v / V_0:	1.000e-09 / 1.000e-09 = 1.000e+00



[...]


Timestep 10 @ 1.000e+00s
___________________________________________________________________________________________________________________________________________________________
                 SOLVER STEP                   |  LIN_IT   LIN_RES    RES_NORM     RES_U     RES_P      RES_J     NU_NORM      NU_U       NU_P       NU_J
___________________________________________________________________________________________________________________________________________________________
  0  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     874  2.358e-06  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00
  1  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     658  2.942e-04  1.544e-01  1.544e-01  1.208e+13  1.855e+06  6.014e-02  7.398e-02  6.014e-02  6.014e-02
  2  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     790  2.206e-06  2.908e-03  2.908e-03  7.302e+10  2.067e+03  2.716e-03  1.433e-03  2.716e-03  2.717e-03
  3  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     893  2.374e-09  1.919e-06  1.919e-06  4.527e+07  4.100e+00  1.672e-06  6.842e-07  1.672e-06  1.672e-06
  4  ASM_R  CONVERGED!
___________________________________________________________________________________________________________________________________________________________
Relative errors:
Displacement:	6.842e-07
Force: 		8.995e-10
Dilatation:	1.528e-06
v / V_0:	1.000e-09 / 1.000e-09 = 1.000e+00
@endcode






使用定时器类，我们可以分辨出代码的哪些部分需要最高的计算费用。对于一个有大量自由度的案例（即高度精细化），下面给出了定时器的典型输出。本教程中的大部分代码都是基于Step-18和其他文章中描述、讨论和演示的优化而开发的。超过93%的时间花在线性求解器上，很明显，对于大型三维问题，可能有必要投资一个更好的求解器。SSOR预处理程序不是多线程的，但对于这类实体问题是有效的。研究使用另一种求解器，如通过Trilinos库提供的求解器，可能是有益的。




@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    | 9.874e+02s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble system right-hand side |        53 | 1.727e+00s |  1.75e-01% |
| Assemble tangent matrix         |        43 | 2.707e+01s |  2.74e+00% |
| Linear solver                   |        43 | 9.248e+02s |  9.37e+01% |
| Linear solver postprocessing    |        43 | 2.743e-02s |  2.78e-03% |
| Perform static condensation     |        43 | 1.437e+01s |  1.46e+00% |
| Setup system                    |         1 | 3.897e-01s |  3.95e-02% |
| Update QPH data                 |        43 | 5.770e-01s |  5.84e-02% |
+---------------------------------+-----------+------------+------------+
@endcode




然后我们用ParaView对两种情况的结果进行了可视化。第一个是最粗的网格和最低阶插值方法。   $Q_1-DGPM_0-DGPM_0$  .第二种是在细化网格上使用 $Q_2-DGPM_1-DGPM_1$ 公式。位移的垂直分量、压力 $\widetilde{p}$ 和扩张 $\widetilde{J}$ 场显示如下。


对于第一种情况，很明显，粗略的空间离散化加上大位移导致了低质量的解决方案（加载比为 $p/p_0=80$ ）。此外，元素之间的压力差非常大。元素上的恒定压力场意味着大的压力梯度没有被捕获。然而，应该注意的是，即使在这种离散性差的情况下，在标准 $Q_1$ 位移公式中会出现的锁定现象也不会出现。块体顶面的跟踪节点的最终垂直位移仍在收敛解的12.5%以内。压力解决方案是非常粗略的，在相邻的单元之间有很大的跳跃。很明显，离施加的牵引力最近的体积经历了压缩，而域的外延则处于膨胀状态。膨胀解场和压力场明显相关，正的膨胀表示正压区域，负的表示压缩区域。正如介绍中所讨论的，压缩性压力有一个负号，而扩张性压力有一个正号。这源于体积应变能量函数的定义，与压力的物理现实的解释相反。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-displacement.png" alt="">
	<p align="center">
        Z-displacement solution for the 3-d problem.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-pressure.png" alt="">
	<p align="center">
        Discontinuous piece-wise constant pressure field.
	</p>
    </td>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-dilatation.png" alt="">
	<p align="center">
        Discontinuous piece-wise constant dilatation field.
	</p>
    </td>
  </tr>
</table> 

结合空间细化和高阶插值方案，产生了高质量的解决方案。三个网格细化加上 $Q_2-DGPM_1-DGPM_1$ 公式产生的结果清楚地抓住了问题的力学原理。牵引面的变形得到了很好的解决。我们现在可以观察到所施加的牵引力的实际范围，最大的力被施加在表面的中心点，导致最大的压缩。尽管领域中出现了很高的应变，特别是在施加牵引力的区域的边界，但解决方案仍然是准确的。压力场被捕捉到的细节比以前多得多。压缩和膨胀区域之间有明显的区别和过渡，压力场的线性近似允许在子元素尺度上对压力进行精细的可视化。然而，应该注意的是，压力场仍然是不连续的，可以在一个连续的网格上进行平滑处理，以达到后期处理的目的。




 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-displacement.png" alt="">
	<p align="center">
        Z-displacement solution for the 3-d problem.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-pressure.png" alt="">
	<p align="center">
        Discontinuous linear pressure field.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-dilatation.png" alt="">
	<p align="center">
        Discontinuous linear dilatation field.
	</p>
    </td>
  </tr>
</table> 

这一简要的分析结果表明，三场公式能够有效地规避高度不可压缩介质的体积锁定。混合配方能够准确模拟近乎不可压缩的块体在压缩状态下的位移。命令行输出显示，在极度压缩下的体积变化导致泊松比为0.4999的体积变化小于0.01%。

在运行时间方面，对于类似的自由度数量， $Q_2-DGPM_1-DGPM_1$ 公式往往比 $Q_1-DGPM_0-DGPM_0$ 的计算成本更高（通过为低阶插值增加一个额外的网格细化级别产生）。下图显示了在一台4核（8线程）机器上连续运行的一批测试的情况。高阶方法计算时间的增加可能是由于高阶元素所需的带宽增加。如前所述，使用更好的求解器和预处理程序可以减轻使用高阶公式的费用。据观察，对于给定的问题，与单线程的SSOR预处理程序相比，使用多线程的Jacobi预处理程序可以减少72%的计算运行时间（在最坏的情况下是具有大量自由度的高阶公式）。然而，根据作者的经验，雅可比预处理方法可能不适合某些涉及替代构成模型的有限应变问题。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Normalised_runtime.png" alt="">
	<p align="center">
        Runtime on a 4-core machine, normalised against the lowest grid resolution $Q_1-DGPM_0-DGPM_0$ solution that utilised a SSOR preconditioner.
	</p>
    </td>
  </tr>
</table> 


最后，下面展示了两个不同级别的网格细化的2维问题的位移解决方案的结果。很明显，由于二维模拟的额外约束，所产生的位移场虽然在质量上相似，但与三维模拟的情况不同。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.2d-gr_2.png" alt="">
	<p align="center">
        Y-displacement solution in 2-d for 2 global grid refinement levels.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.2d-gr_5.png" alt="">
	<p align="center">
        Y-displacement solution in 2-d for 5 global grid refinement levels.
	</p>
    </td>
  </tr>
</table> 

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这项工作有许多明显的延伸。

- 首先，可以在自由能函数中加入一个额外的约束条件，以便在材料中强制执行高度的不可压缩性。一个额外的拉格朗日乘数将被引入，但这可以最容易地使用增强的拉格朗日乘数的原则来处理。这在  <em>  Simo和Taylor (1991)  </em>  中得到了证明。

- 这个模型中使用的构成关系是比较基本的。将材料类分成两个独立的类，一个处理体积响应，另一个处理等温线响应，并产生一个通用的材料类（即具有抽象的虚拟函数，派生类必须实现），允许增加更复杂的材料模型，这可能是有益的。这些模型可以包括其他超弹性材料、塑性和粘弹性材料以及其他材料。

- 该程序是为解决单节点多核机器上的问题而开发的。只要稍加努力，该程序就可以通过使用Petsc或Trilinos扩展到大规模的计算环境，使用的技术与step-40中演示的类似。这主要涉及对设置、装配、 <code>PointHistory</code> 和线性求解器例程的修改。

- 由于该程序假定为准静态平衡，为了研究惯性效应很重要的问题，例如涉及冲击的问题，有必要进行扩展以包括动态效应。

- 对于高度非线性问题，负载和解的限制程序可能是必要的。可以增加一个线搜索算法，将步长限制在牛顿增量内，以确保最佳收敛性。也可能需要使用负载限制方法，如Riks方法，来解决涉及几何不稳定性的不稳定问题，如屈曲和快穿。

- 许多物理问题涉及接触。有可能将物体间的摩擦或无摩擦接触的影响纳入这个程序。这将涉及到在自由能函数中增加一个额外的项，因此需要增加装配程序。我们还需要管理接触问题（检测和应力计算）本身。在自由能函数中增加惩罚项的一个替代方法是使用主动集方法，如步骤41中使用的方法。

- 使用LinearOperators的完整缩减程序已经被编码到线性求解器例程中。这也可以通过应用schur_complement()操作符来实现，以更自动化的方式缩减一个或多个字段。

- 最后，自适应网格细化，如步骤6和步骤18所示，可以提供额外的求解精度。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-44.cc"
*/
