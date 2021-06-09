/**
@page step_71 The step-71 tutorial program
@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#AmotivationWhywouldIusethesetools">A motivation: Why would I use these tools?</a>
        <li><a href="#Theoryformagnetomechanicalmaterials">Theory for magneto-mechanical materials</a>
      <ul>
        <li><a href="#Thermodynamicprinciples">Thermodynamic principles</a>
        <li><a href="#Constitutivelaws">Constitutive laws</a>
      <ul>
        <li><a href="#Magnetoelasticconstitutivelaw">Magnetoelastic constitutive law</a>
        <li><a href="#Magnetoviscoelasticconstitutivelaw">Magneto-viscoelastic constitutive law</a>
      </ul>
      </ul>
        <li><a href="#Rheologicalexperiment">Rheological experiment</a>
        <li><a href="#Suggestedliterature">Suggested literature</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#AnintroductoryexampleThefundamentalsofautomaticandsymbolicdifferentiation">An introductory example: The fundamentals of automatic and symbolic differentiation</a>
      <ul>
        <li><a href="#Ananalyticalfunction">An analytical function</a>
        <li><a href="#Computingderivativesusingautomaticdifferentiation">Computing derivatives using automatic differentiation</a>
        <li><a href="#Handcalculatedderivativesoftheanalyticalsolution">Hand-calculated derivatives of the analytical solution</a>
        <li><a href="#Computingderivativesusingsymbolicdifferentiation">Computing derivatives using symbolic differentiation</a>
        <li><a href="#TheSimpleExamplerunfunction">The SimpleExample::run() function</a>
      </ul>
        <li><a href="#AmorecomplexexampleUsingautomaticandsymbolicdifferentiationtocomputederivativesatcontinuumpoints">A more complex example: Using automatic and symbolic differentiation to compute derivatives at continuum points</a>
      <ul>
        <li><a href="#Constitutiveparameters">Constitutive parameters</a>
        <li><a href="#ConstitutivelawsBaseclass">Constitutive laws: Base class</a>
        <li><a href="#Magnetoelasticconstitutivelawusingautomaticdifferentiation">Magnetoelastic constitutive law (using automatic differentiation)</a>
        <li><a href="#Magnetoviscoelasticconstitutivelawusingsymbolicalgebraanddifferentiation">Magneto-viscoelastic constitutive law (using symbolic algebra and differentiation)</a>
      </ul>
        <li><a href="#AmorecomplexexamplecontinuedParametersandhandderivedmaterialclasses">A more complex example (continued): Parameters and hand-derived material classes</a>
      <ul>
        <li><a href="#Magnetoelasticconstitutivelawhandderived">Magnetoelastic constitutive law (hand-derived)</a>
        <li><a href="#Magnetoviscoelasticconstitutivelawhandderived">Magneto-viscoelastic constitutive law (hand-derived)</a>
        <li><a href="#Rheologicalexperimentparameters">Rheological experiment parameters</a>
        <li><a href="#RheologicalexperimentParallelplaterotationalrheometer">Rheological experiment: Parallel plate rotational rheometer</a>
        <li><a href="#TheCoupledConstitutiveLawsrunfunction">The CoupledConstitutiveLaws::run() function</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Introductoryexample">Introductory example</a>
        <li><a href="#Constitutivemodelling">Constitutive modelling</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-71/doc/intro.dox

 <br> 

<i>This program was contributed by Jean-Paul Pelteret.
</i>




<a name="Introduction"></a><h1>Introduction</h1>


本教程的目的很简单，就是介绍[自动](https://en.wikipedia.org/wiki/Automatic_differentiation)和[符号微分](https://en.wikipedia.org/wiki/Computer_algebra)（分别缩写为AD和SD）的基本原理。人们可以在源代码中描述一个函数 $\mathbf f(\mathbf x)$ ，并自动获得导数 $\nabla \mathbf f(\mathbf x)$ （"Jacobian"）、 $\nabla^2 \mathbf f(\mathbf x)$ （"Hessian"）等的表示方法，而无需编写额外的代码行。这样做对解决非线性或优化问题很有帮助，因为人们希望在代码中只描述非线性方程或目标函数，而不必同时提供它们的导数（这对解决非线性问题的牛顿方法或寻找最小化器是必要的）。

由于AD和SD工具在某种程度上独立于有限元和边界值问题，本教程将与你之前可能读过的其他教程不同。它将特别关注这些框架是如何工作的，以及它们背后的原理和思想，并放弃在有限元模拟的直接背景下看待它们。

事实上，我们将研究两组不同的问题，它们的复杂程度大不相同，但当框架正确时，有足够的相似性，同样的AD和SD框架可以被利用。通过这些例子，我们的目的是建立起对使用AD和SD工具所需步骤的理解，以及它们之间的区别，并希望能找出它们可以立即用于改进或简化现有代码的地方。

你想知道什么是AD和SD，这是可信的，首先。好吧，这个问题很容易回答，但如果没有上下文，就没有很好的洞察力。因此，我们不打算在这个介绍中涉及这个问题，而是将其推迟到第一个介绍性的例子中，在这个例子的展开过程中，我们将列出关键点。作为补充，我们应该提到，这两个框架的核心理论在 @ref auto_symb_diff 模块中都有广泛的讨论，所以在此不需要重复。

由于我们必须挑选*个足够有趣的课题来研究，并确定AD和SD在哪里可以有效地使用，所以在教程的后半部分实现的主要问题是对一个耦合的构成法进行建模，特别是一个磁活性材料（具有滞后效应）。作为一种介绍的手段，在介绍的后面将介绍该类材料的一些基础理论。自然，这不是一个广泛受众感兴趣的领域（甚至不是一类材料）。因此，作者希望在前面表示，这个理论和任何后续的推导都不能被认为是本教程的重点。相反，请牢记从相对无害的构成法则描述中产生的问题的复杂性，以及我们可能（在边界值问题的背景下）需要从中推导出什么。我们将在一个有代表性的连续体点的水平上用这些构成法则进行一些计算（所以，仍然是在连续体力学的领域），并将产生一些基准结果，我们可以围绕这些结果对计算性能的主题进行最后讨论。

一旦我们有了可以建立进一步概念的基础，我们将看到如何在有限元（而不是连续体）水平上特别利用AD：这是在步骤-72和步骤-33中涉及的一个主题。但在此之前，让我们花点时间思考一下为什么我们可能要考虑使用这些工具，以及它们可能给你带来什么好处。




<a name="AmotivationWhywouldIusethesetools"></a><h3>A motivation: Why would I use these tools?</h3>


使用AD或SD的主要驱动力通常是，有一些情况需要进行区分，而且这样做有足够的挑战性，使得使用外部工具来执行该特定任务的前景具有吸引力。对AD或SD最有用的情况进行广泛分类，包括（但可能不限于）以下情况。

- <b>Rapid prototyping:</b>对于一类新的问题，你试图快速实现一个解决方案，并希望去除一些复杂的细节（在数学以及代码本身的组织结构方面）。你可能愿意证明任何额外的计算成本是合理的，这将被重组你的代码或修改问题中引入一些复杂的非线性的部分的敏捷性所抵消，只需最小的努力。

- <b>Complex problems:</b>很可能有些问题恰好有一个非线性，对线性化或手工制定有极大的挑战。   让一个在大多数情况下稳健、可靠和准确的工具来为你解决这个挑战，可能会减轻实现某些问题的痛苦。这方面的例子包括第15步，我们解决的非线性PDE的导数并不难推导，但足够繁琐，以至于人们在手工操作时必须注意，而且实现牛顿步骤的相应有限元公式所需的时间不仅仅是实现双线性形式一般所需的几行；第33步（我们实际使用AD）是一个更极端的例子。

- <b>Verification:</b> 对于表现出非线性响应的材料和模拟，准确而非近似的材料切线（机械工程师对材料定律的导数使用的术语）可能是收敛和发散行为之间的区别，特别是在高外部（或耦合）载荷下。   随着问题复杂性的增加，引入细微的（或者，也许不是那么细微的）错误的机会也在增加，这些错误会产生可预见的负面结果。   此外，通过验证实现是完全正确的，也有很多好处。例如，某些类别的问题已知会表现出不稳定性，因此，当你在非线性求解器（例如牛顿方法）中开始失去二次收敛时，那么这对研究者来说可能不是一个巨大的惊喜。然而，很难（如果不是不可能）区分以下两种收敛行为：一种是你接近不稳定的解时产生的收敛行为，另一种是你在材料或有限元线性化中出现了错误，并因此开始偏离最佳收敛路径。例如，拥有一种验证构成法线性化实现的正确性的方法，也许是你用来捕捉这种错误的唯一有意义的方法，假设你没有其他人来检查你的代码。   值得庆幸的是，通过一些战术性的编程，可以很直接地将代码结构化以便重复使用，这样你就可以在生产代码中使用相同的类，并直接在例如单元测试框架中验证它们。

这个教程程序将有两个部分。一部分，我们只是用一组简单的例子来介绍deal.II中自动和符号微分支持的基本思想；另一部分，我们将其应用于一个现实的但更复杂的案例。对于这后半部分，下一节将提供一些关于磁性机械材料的背景--如果你想了解的只是AD和SD的实际情况，你可以跳过这一节，但如果你对如何将AD和SD应用于具体的情况感兴趣，你可能想读完这一节。




<a name="Theoryformagnetomechanicalmaterials"></a><h3>Theory for magneto-mechanical materials</h3>


<a name="Thermodynamicprinciples"></a><h4>Thermodynamic principles</h4>


作为介绍我们将用来为磁活性聚合物建模的磁-机械耦合材料法的前奏，我们将首先对这些构成法则必须认同的突出的热力学进行非常简洁的总结。这里总结的理论基础，由Truesdell和Toupin  @cite Truesdell1960a  以及Coleman和Noll  @cite Coleman1963a  详细描述，并遵循Holzapfel  @cite Holzapfel2007a  所提出的逻辑。

从热力学第一定律出发，并遵循一些技术假设，可以证明动能加内能率与外部来源提供给系统的功率之间的平衡是由以下关系给出的，即左边是一个（任意）体积 $V$ 的能量变化率，右边是作用于该体积的力的总和。

@f[
  D_{t} \int\limits_{V} \left[
    \frac{1}{2} \rho_{0} \mathbf{v} \cdot \mathbf{v}
    + U^{*}_{0} \right] dV
= \int\limits_{V} \left[
  \rho_{0} \mathbf{v} \cdot \mathbf{a}
  + \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}
  + \mathbb{E} \cdot \dot{\mathbb{D}}


  - D_{t} M^{*}_{0}


  - \nabla_{0} \cdot \mathbf{Q}
  + R_{0} \right] dV .


@f]

这里 $D_{t}$ 代表总的时间导数， $\rho_{0}$ 是在拉格朗日参考框架下测量的材料密度， $\mathbf{v}$ 是材料速度， $\mathbf{a}$ 是其加速度， $U^{*}_{0}$ 是每单位参考体积的内能， $\mathbf{P}^{\text{tot}}$ 是总皮拉应力张量， $\dot{\mathbf{F}}$ 是变形梯度张量的时间速率， $\boldsymbol{\mathbb{H}}$ 和 $\boldsymbol{\mathbb{B}}$ 分别是磁场向量和磁感应（或磁通密度）向量， $\mathbb{E}$ 和 $\mathbb{D}$ 是电场向量和电位移向量， $\mathbf{Q}$ 和 $R_{0}$ 代表参考热通向量和热源。材料微分算子 $\nabla_{0} (\bullet) \dealcoloneq \frac{d(\bullet)}{d\mathbf{X}}$ ，其中 $\mathbf{X}$ 是材料位置向量。通过一些条款的重排，引用积分体积 $V$ 的任意性，总的内部能量密度率 $\dot{E}_{0}$ 可以被确定为

@f[
  \dot{E}_{0}
= \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}
  + \mathbb{E} \cdot \dot{\mathbb{D}}


  - \nabla_{0} \cdot \mathbf{Q}
  + R_{0} .


@f]

总的内能不仅包括由于机械变形（第一项）、热通量和热源（第四项和第五项）而产生的贡献，还包括由于储存在磁场和电场本身的内在能量（分别为第二项和第三项）。

热力学第二定律，也被称为熵不平等原则，告诉我们某些热力学过程是不可逆的。在考虑了总熵和熵输入的速度后，可以得出克劳修斯-杜姆不等式。在局部形式下（以及在物质配置中），其内容为

@f[
  \theta \dot{\eta}_{0}


  - R_{0}
  + \nabla_{0} \cdot \mathbf{Q}


  - \frac{1}{\theta} \nabla_{0} \theta \cdot \mathbf{Q}
  \geq 0 .


@f]

量 $\theta$ 是绝对温度， $\eta_{0}$ 代表每单位参考体积的熵值。

用它来代替热力学第一定律结果中的 $R_{0} - \nabla_{0} \cdot \mathbf{Q}$ ，我们现在有了这样的关系

@f[
  \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}
  + \mathbb{E} \cdot \dot{\mathbb{D}}
  + \theta \dot{\eta}_{0}


  - \dot{E}_{0}


  - \frac{1}{\theta} \nabla_{0} \theta \cdot \mathbf{Q}
  \geq 0 .


@f]

傅里叶定律告诉我们，热量从高温区域流向低温区域，根据这一定律，最后一项总是正的，可以忽略不计。这使得局部耗散的不等式变成了

@f[
  \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}
  + \mathbb{E} \cdot \dot{\mathbb{D}}


  - \left[ \dot{E}_{0} - \theta \dot{\eta}_{0}  \right]
  \geq 0 .


@f]

据推测 @cite Holzapfel2007a ，Legendre变换

@f[
  \psi^{*}_{0}
= \psi^{*}_{0} \left( \mathbf{F}, \boldsymbol{\mathbb{B}}, \mathbb{D}, \theta \right)
= E_{0} - \theta \eta_{0} ,


@f]

从中我们可以定义具有所述参数化的自由能密度函数 $\psi^{*}_{0}$ ，它存在并且有效。取此方程的材料速率并将其代入局部耗散不等式，结果是通用表达式

@f[
  \mathcal{D}_{\text{int}}
  = \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}
  + \mathbb{E} \cdot \dot{\mathbb{D}}


  - \dot{\theta} \eta_{0}


  - \dot{\psi}^{*}_{0} \left( \mathbf{F}, \boldsymbol{\mathbb{B}}, \mathbb{D}, \theta \right)
  \geq 0 .


@f]

在等温条件的假设下，并且电场不会以一种被认为是不可忽视的方式激发材料，那么这个耗散不等式就会简化为

@f[
  \mathcal{D}_{\text{int}}
  = \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}


  - \dot{\psi}^{*}_{0} \left( \mathbf{F}, \boldsymbol{\mathbb{B}} \right)
  \geq 0 .


@f]



<a name="Constitutivelaws"></a><h4>Constitutive laws</h4>


当考虑到表现出机械耗散行为的材料时，可以证明这可以通过用代表内部变量的额外参数增加材料自由能密度函数的方式在耗散不等式中得到体现  @cite Holzapfel1996a  。因此，我们把它写成

@f[
  \mathcal{D}_{\text{int}}
  = \mathbf{P}^{\text{tot}} : \dot{\mathbf{F}}
  + \boldsymbol{\mathbb{H}} \cdot \dot{\boldsymbol{\mathbb{B}}}


  - \dot{\psi}^{*}_{0} \left( \mathbf{F}, \mathbf{F}_{v}^{i}, \boldsymbol{\mathbb{B}} \right)
  \geq 0 .


@f]

其中 $\mathbf{F}_{v}^{i} = \mathbf{F}_{v}^{i} \left( t \right)$ 代表与第i个机械耗散（粘性）机制相关的内部变量（其作用类似于变形梯度的测量）。从它的参数化可以推断出，这些内部参数中的每一个都被认为是在时间中演变的。目前，自由能密度函数 $\psi^{*}_{0}$ 是以磁感应 $\boldsymbol{\mathbb{B}}$ 为参数的。这是自然的参数化，是所考虑的平衡法的结果。如果这样一类材料被纳入到有限元模型中，将确定需要采用某种磁问题的表述，即磁矢量势表述。这有它自己的一套挑战，所以在可能的情况下，更简单的磁标量势表述可能是首选。在这种情况下，磁性问题需要在磁场方面进行参数化  $\boldsymbol{\mathbb{H}}$  。为了进行这种重新参数化，我们执行最后的Legendre变换

@f[
  \tilde{\psi}_{0} \left( \mathbf{F}, \mathbf{F}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  = \psi^{*}_{0} \left( \mathbf{F}, \mathbf{F}_{v}^{i}, \boldsymbol{\mathbb{B}} \right)


  - \boldsymbol{\mathbb{H}} \cdot \boldsymbol{\mathbb{B}} .


@f]

同时，我们可以利用材料框架无所谓的原则，以便用对称的变形措施来表达能量密度函数。

@f[
  \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  = \tilde{\psi}_{0} \left( \mathbf{F}, \mathbf{F}_{v}^{i}, \boldsymbol{\mathbb{H}} \right) .


@f]

这两个转换的结果（撇开相当多的明确和隐藏的细节）使减少耗散不等式的最终表达式为

@f[
  \mathcal{D}_{\text{int}}
  = \mathbf{S}^{\text{tot}} : \frac{1}{2} \dot{\mathbf{C}}


  - \boldsymbol{\mathbb{B}} \cdot \dot{\boldsymbol{\mathbb{H}}}


  - \dot{\psi}_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  \geq 0 .


@f]

注意右侧第二项的符号变化，以及时间导数向磁感应矢量的转移）。应力量 $\mathbf{S}^{\text{tot}}$ 被称为总Piola-Kirchhoff应力张量，其能量共轭物 $\mathbf{C} = \mathbf{F}^{T} \cdot \mathbf{F}$ 是右Cauchy-Green变形张量， $\mathbf{C}_{v}^{i} = \mathbf{C}_{v}^{i} \left( t \right)$ 是与`i`th机械耗散（粘性）机制相关的重新参数化内部变量。

对能量密度函数的材料速率进行扩展，并对各种项进行重排，得出的表达式是

@f[
  \mathcal{D}_{\text{int}}
  = \left[ \mathbf{S}^{\text{tot}} - 2 \frac{\partial \psi_{0}}{\partial \mathbf{C}} \right] : \frac{1}{2} \dot{\mathbf{C}}


  - \sum\limits_{i}\left[ 2 \frac{\partial \psi_{0}}{\partial \mathbf{C}_{v}^{i}} \right] : \frac{1}{2} \dot{\mathbf{C}}_{v}^{i}
  + \left[ - \boldsymbol{\mathbb{B}} - \frac{\partial \psi_{0}}{\partial \boldsymbol{\mathbb{H}}} \right] \cdot \dot{\boldsymbol{\mathbb{H}}}
  \geq 0 .


@f]

在这一点上，值得注意的是[偏导数](https://en.wikipedia.org/wiki/Partial_derivative)  $\partial \left( \bullet \right)$  的使用。这是一个重要的细节，对于本教程中的某个设计选择是很重要的。简单提醒一下这意味着什么，一个多变量函数的偏导返回该函数相对于其中一个变量的导数，而其他变量保持不变。

@f[
  \frac{\partial f\left(x, y\right)}{\partial x}
  = \frac{d f\left(x, y\right)}{d x} \Big\vert_{y} .


@f]

更具体到耗散不等式所编码的内容（用非常普遍的自由能密度函数 $\psi_{0}$ ，其参数化尚待正式确定），如果输入变量之一是另一个变量的函数，它也被保持不变，链式规则不再传播，而计算总导数将意味着明智地使用链式规则。通过比较以下两个语句可以更好地理解这一点。

@f{align*}
  \frac{\partial f\left(x, y\left(x\right)\right)}{\partial x}
  &= \frac{d f\left(x, y\left(x\right)\right)}{d x} \Big\vert_{y} \\
  \frac{d f\left(x, y\left(x\right)\right)}{d x}
  &= \frac{d f\left(x, y\left(x\right)\right)}{d x} \Big\vert_{y}
   + \frac{d f\left(x, y\left(x\right)\right)}{d y} \Big\vert_{x} \frac{d y\left(x\right)}{x} .


@f}



回到问题的热力学，我们接下来利用数量的任意性  $\dot{\mathbf{C}}$  和  $\dot{\boldsymbol{\mathbb{H}}}$  ，通过应用科尔曼-诺尔程序  @cite Coleman1963a  ，  @cite Coleman1967a  。这导致了对动力学共轭量的识别

@f[
  \mathbf{S}^{\text{tot}}
  = \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  \dealcoloneq 2 \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}} , \\
  \boldsymbol{\mathbb{B}}
  = \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  \dealcoloneq - \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}} .


@f]

(再次注意，在这个广义的设置中，使用偏导数来定义应力和磁感应)。从耗散功率中剩下的条款（即那些与机械耗散机制有关的条款）来看，如果假定它们是相互独立的，那么，对于每个机制`i`。

@f[
  \frac{\partial \psi_{0}}{\partial \mathbf{C}_{v}^{i}} : \dot{\mathbf{C}}_{v}^{i}
  \leq 0 .


@f]

这一约束必须通过适当选择自由能函数以及仔细考虑内部变量的演化规律来满足。

在构成模型中没有耗散机制的情况下（例如，如果要建模的材料是磁超弹性的），那么自由能密度函数 $\psi_{0} = \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)$ 减少到存储能量密度函数，总应力和磁感应可以被简化

@f{align*}{
  \mathbf{S}^{\text{tot}}
  = \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
  &\dealcoloneq 2 \frac{d \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C}} , \\
  \boldsymbol{\mathbb{B}}
  = \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
  &\dealcoloneq - \frac{d \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}}} ,


@f}

其中算子 $d$ 表示总导数操作。

为了完整起见，应力张量和磁感应的线性化在四阶总参考弹性张量 $\mathcal{H}^{\text{tot}} $ 、二阶磁静力张量 $\mathbb{D}$ 和三阶总参考磁弹性耦合张量 $\mathfrak{P}^{\text{tot}}$ 中得到体现。无论 $\mathbf{S}^{\text{tot}}$ 和 $\boldsymbol{\mathbb{B}}$ 的参数化如何，这些量都可以通过以下方式计算出来

@f{align*}{
  \mathcal{H}^{\text{tot}}
  &= 2 \frac{d \mathbf{S}^{\text{tot}}}{d \mathbf{C}} , \\
  \mathbb{D}
  &= \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}} , \\
  \mathfrak{P}^{\text{tot}}
  &= - \frac{d \mathbf{S}^{\text{tot}}}{d \boldsymbol{\mathbb{H}}} , \\
  \left[ \mathfrak{P}^{\text{tot}} \right]^{T}
  &= 2 \frac{d \boldsymbol{\mathbb{B}}}{d \mathbf{C}} .


@f}

对于速率依赖性材料的情况，这扩展为

@f{align*}{
  \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  &= 4 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C} \otimes d \mathbf{C}} , \\
  \mathbb{D} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  &= -\frac{d^{2} \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}} , \\
  \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)
  &= - 2 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}} \otimes d \mathbf{C}} , \\
  \left[ \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)  \right]^{T}
  &= - 2 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}^{i}, \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C} \otimes d \boldsymbol{\mathbb{H}}} ,


@f}

而对于与速率无关的材料，其线性化为

@f{align*}{
  \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
  &= 4 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C} \otimes d \mathbf{C}} , \\
  \mathbb{D} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
  &= -\frac{d^{2} \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}} , \\
  \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
  &= - 2 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}} \otimes d \mathbf{C}} , \\
  \left[ \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)  \right]^{T}
  &= - 2 \frac{d^{2} \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C} \otimes d \boldsymbol{\mathbb{H}}} .


@f}

它们之间的微妙区别是在计算第一个导数时应用了偏导。我们稍后会看到这对这个具体应用中AD与SD的选择有什么影响。现在，我们将简单介绍在本教程中实现的两种具体材料。

<a name="Magnetoelasticconstitutivelaw"></a><h5>Magnetoelastic constitutive law</h5>


我们要考虑的第一种材料是受磁超弹性构成法支配的材料。这种材料对变形和浸泡在磁场中都有反应，但没有表现出时间或历史相关的行为（如通过粘性阻尼或磁滞的耗散，等等）。这种材料的*存储能量密度函数*只以（当前）场变量为参数，而不是它们的时间导数或过去的值。

我们将选择能量密度函数，它既能捕捉到由于变形和磁化而储存在材料中的能量，也能捕捉到储存在磁场本身的能量，它是

@f[
  \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
= \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
    \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    \right]
+ \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)


- \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    \boldsymbol{\mathbb{H}} \right]


@f]

与

@f[
  f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
= 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    \boldsymbol{\mathbb{H}}}
      {\left(h_{e}^{\text{sat}}\right)^{2}} \right)


@f]

其中变量 $d = \text{tr}(\mathbf{I})$ （ $\mathbf{I}$ 是秩-2身份张量）代表空间维度， $\mathbf{F}$ 是变形梯度张量。为了给 $\psi_{0}$ 的各个组成部分提供一些简单的背景，前两个项与（超弹性）Neohookean材料的储能密度函数非常相似。这里使用的东西与Neohookean材料的唯一区别是弹性剪切模量被磁场敏感的饱和函数 $f_{\mu_{e}}
\left( \boldsymbol{\mathbb{H}} \right)$ 缩放（见 @cite Pelteret2018a ，公式29）。这个函数实际上将导致材料在强磁场的存在下变硬。由于它受一个sigmoid型函数的支配，剪切模量将渐进地收敛于指定的饱和剪切模量。还可以证明， $\psi_{0}$ 中的最后一项是磁场的储能密度函数（从第一原理中得出），由相对渗透率常数缩放。这个定义共同意味着材料是线性磁化的，也就是说，磁化矢量和磁场矢量是对齐的。(这在以电流形式陈述的磁能中当然不明显，但当磁感应和磁化从 $\psi_{0}$ 中导出，并且所有磁场都以 <em> 的电流配置 </em> 表示时，这种关联性就变得很清楚了)。至于磁感应、应力张量和各种材料切线的具体内容，我们将把这些内容推迟到教程正文中介绍，在那里定义了构成法的完整、无辅助的实施。

<a name="Magnetoviscoelasticconstitutivelaw"></a><h5>Magneto-viscoelastic constitutive law</h5>


我们将制定的第二个材料是一个具有单一耗散机制`i`的磁-粘弹性材料。我们将考虑的*自由能量密度函数*被定义为

@f{align*}{
  \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
  \right)
&= \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
+ \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
\boldsymbol{\mathbb{H}} \right)
\\ \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
&= \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
\right)
    \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
    \right]
+ \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)


- \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
    \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
    \boldsymbol{\mathbb{H}} \right]
\\ \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
\boldsymbol{\mathbb{H}} \right)
&= \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
\right)
    \left[ \mathbf{C}_{v} : \left[
      \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
      \mathbf{C} \right] - d - \ln\left(
      \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]


@f}

与

@f[
  f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
= 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
    \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    \boldsymbol{\mathbb{H}}}
      {\left(h_{e}^{\text{sat}}\right)^{2}} \right)


@f]



@f[
  f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
= 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
    \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
    \boldsymbol{\mathbb{H}}}
      {\left(h_{v}^{\text{sat}}\right)^{2}} \right)


@f]

和进化法

@f[
  \dot{\mathbf{C}}_{v} \left( \mathbf{C} \right)
= \frac{1}{\tau} \left[
      \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
        \mathbf{C}\right]^{-1}


    - \mathbf{C}_{v} \right]


@f]

为内部粘性变量。我们已经选择了能量的磁弹性部分 $\psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)$ 来匹配我们探索的第一个材料模型，所以这部分不需要进一步解释。至于粘性部分 $\psi_{0}^{MVE}$ ，自由能的这一部分（与粘性变形张量的演化规律一起）取自 @cite Linder2011a （由 @cite Pelteret2018a 中描述的粘性饱和函数进行额外缩放）。它是在一个热力学上一致的框架中得出的，其核心是在微观层面上模拟聚合物链的运动。

要超越这一点，我们还需要考虑进化规律的时间离散化。选择隐式一阶逆向差分方案，那么

@f[
  \dot{\mathbf{C}}_{v}
\approx \frac{\mathbf{C}_{v}^{(t)} - \mathbf{C}_{v}^{(t-1)}}{\Delta t}
= \frac{1}{\tau} \left[
      \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
        \mathbf{C}\right]^{-1}


    - \mathbf{C}_{v}^{(t)} \right]


@f]

其中上标 $(t)$ 表示该数量是在当前时间步长中提取的， $(t-1)$ 表示在前一时间步长中提取的数量（即历史变量）。时间段大小 $\Delta t$ 是当前时间与上一时间段的差。将这些条款重新排列，使当前时间的所有内部变量量都在方程的左侧，我们可以得到

@f[
\mathbf{C}_{v}^{(t)}
= \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
    \mathbf{C}_{v}^{(t-1)}
  + \frac{\Delta t}{\tau_{v}}
    \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
    \mathbf{C} \right]^{-1}
  \right]


@f]

匹配 @cite Linder2011a 公式54。

<a name="Rheologicalexperiment"></a><h3>Rheological experiment</h3>


现在我们已经展示了所有这些关于热力学和磁力学理论以及构成模型的公式，让我们概述一下这个程序将对所有这些做什么。我们希望对我们制定的材料定律做一些*有意义的事情，因此将它们置于一些机械和磁载荷条件下是有意义的，这些条件在某种程度上代表了在应用或实验室环境中可能发现的一些条件。实现这一目标的方法之一是将这些构成法则嵌入到有限元模型中，以模拟一个设备。不过，在这个例子中，我们将保持简单（毕竟我们关注的是自动和符号微分概念），并将找到一种简明的方法，使用加载条件的分析表达式忠实地复制工业标准的流变学实验。

我们将重现的流变学实验，它理想化了一个用于表征磁活性聚合物的实验室实验，详见 @cite Pelteret2018a （以及 @cite Pelteret2019a ，其中与真实世界的实验一起记录）。下面的图片提供了对问题设置的直观描述。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img
        src="https://www.dealii.org/images/steps/developer/step-71.parallel_plate-geometry.png"
        alt="" height="300">
  <p align="center">
        The basic functional geometry of the parallel-plate rotational
        rheometer. The smooth rotor (blue) applies a torque to an
        experimental sample (red) of radius $r$ and height $H$ while an
        axially aligned magnetic field generated by a a
        magneto-rheological device. Although the time-dependent
        deformation profile of the may be varied, one common experiment
        would be to subject the material to a harmonic torsional
        deformation of constant amplitude and frequency $\omega$.
  </p>
    </td>
    <td align="center">
        <img
        src="https://www.dealii.org/images/steps/developer/step-71.parallel_plate-kinematics.png"
        alt="" height="300">
  <p align="center">
        Schematic of the kinematics of the problem, assuming no
        preloading or compression of the sample. A point $\mathbf{P}$
        located at azimuth $\Theta$ is displaced to location $\mathbf{p}$
        at azimuth $\theta = \Theta + \alpha$.
  </p>
    </td>
  </tr>
</table> 

假设正在测试的是不可压缩的介质，并且通过样品厚度的变形曲线是线性的，那么在样品内某个测量点 $\mathbf{X}$ 的位移，用径向坐标表示，就是

@f{align*}
  r(\mathbf{X})
  &= \frac{R(X_{1}, X_{2})}{\sqrt{\lambda_{3}}} , \\
  \theta(\mathbf{X})
  & = \Theta(X_{1}, X_{2}) + \underbrace{\tau(t)
       \lambda_{3} X_{3}}_{\alpha(X_{3}, t)} , \\
  z(\mathbf{X})
  &= \lambda_{3} X_{3}


@f}

其中 $R(X_{1}, X_{2})$ 和 $\Theta(X_{1}, X_{2})$ 是半径在

-- 的角度， $\lambda_{3}$ 是（恒定的）轴向变形， $\tau(t) = \frac{A}{RH} \sin\left(\omega t\right)$ 是每单位长度的随时间变化的扭转角，将使用固定振幅的正弦波重复振荡 $A$ 来规定。磁场是轴向排列的，即在 $X_{3}$ 方向。

这总结了我们在流变样品内任何一点上全面描述理想化载荷所需的一切。我们将以这样的方式设置问题，即我们在这个样品中 "挑选 "一个有代表性的点，并使其在恒定的轴向变形（默认为压缩载荷）和恒定的、轴向施加的磁场中受到谐波剪切变形。我们将记录该点的应力和磁感应强度，并将数据输出到文件中进行后处理。尽管对这个特定的问题来说没有必要，我们也将计算切线。尽管它们没有直接用于这个特定的工作，但这些二阶导数是在有限元模型中嵌入构成法所需要的（这项工作的一个可能的扩展）。因此，我们将利用这个机会，用辅助微分框架来检查我们的手工计算是否正确。

<a name="Suggestedliterature"></a><h3>Suggested literature</h3>


除了已经提到的 @ref auto_symb_diff 模块外，以下是一些更详细讨论的参考资料

- 磁力学，以及自动分化框架的某些方面。   @cite Pao1978a  ,  @cite Pelteret2019a  , 和

- 使用AD和/或SD实现有限元框架的自动化：  @cite Logg2012a  ,  @cite Korelc2016a  。

 <br> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * We start by including all the necessary deal.II header files and some C++
 * related ones.
 * This first header will give us access to a data structure that will allow
 * us to store arbitrary data within it.
 * 
 * @code
 * #include <deal.II/algorithms/general_data_storage.h>
 * 
 * @endcode
 * 
 * Next come some core classes, including one that provides an implementation
 * for time-stepping.
 * 
 * @code
 * #include <deal.II/base/discrete_time.h>
 * #include <deal.II/base/numbers.h>
 * #include <deal.II/base/parameter_acceptor.h>
 * #include <deal.II/base/symmetric_tensor.h>
 * #include <deal.II/base/tensor.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/utilities.h>
 * 
 * @endcode
 * 
 * Then some headers that define some useful coordinate transformations and
 * kinematic relationships that are often found in nonlinear elasticity.
 * 
 * @code
 * #include <deal.II/physics/transformations.h>
 * #include <deal.II/physics/elasticity/kinematics.h>
 * #include <deal.II/physics/elasticity/standard_tensors.h>
 * 
 * @endcode
 * 
 * The following two headers provide all of the functionality that we need
 * to perform automatic differentiation, and use the symbolic computer algebra
 * system that deal.II can utilize. The headers of all automatic
 * differentiation and symbolic differentiation wrapper classes, and any
 * ancillary data structures that are required, are all collected inside these
 * unifying headers.
 * 
 * @code
 * #include <deal.II/differentiation/ad.h>
 * #include <deal.II/differentiation/sd.h>
 * 
 * @endcode
 * 
 * Including this header allows us the capability to write output to a
 * file stream.
 * 
 * @code
 * #include <fstream>
 * 
 * 
 * @endcode
 * 
 * As per usual, the entire tutorial program is defined within its own unique
 * namespace.
 * 
 * @code
 * namespace Step71
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="AnintroductoryexampleThefundamentalsofautomaticandsymbolicdifferentiation"></a> 
 * <h3>An introductory example: The fundamentals of automatic and symbolic differentiation</h3>
 * 

 * 
 * Automatic and symbolic differentiation have some magical and mystical
 * qualities. Although their use in a project can be beneficial for a
 * multitude of reasons, the barrier to understanding how to use these
 * frameworks or how they can be leveraged may exceed the patience of
 * the developer that is trying to (reliably) integrate them into their work.
 *   

 * 
 * Although it is the wish of the author to successfully illustrate how these
 * tools can be integrated into workflows for finite element modelling, it
 * might be best to first take a step back and start right from the basics.
 * So to start off with, we'll first have a look at differentiating a "simple"
 * mathematical function using both frameworks, so that the fundamental
 * operations (both their sequence and function) can be firmly established and
 * understood with minimal complication. In the second part of this tutorial
 * we will put these fundamentals into practice and build on them further.
 *   

 * 
 * Accompanying the description of the algorithmic steps to use the frameworks
 * will be a simplified view as to what they *might* be doing in the
 * background. This description will be very much one designed to aid
 * understanding, and the reader is encouraged to view the @ref auto_symb_diff
 * module documentation for a far more formal description into how these tools
 * actually work.
 *   

 * 
 * 
 * <a name="Ananalyticalfunction"></a> 
 * <h4>An analytical function</h4>
 * 
 * @code
 *   namespace SimpleExample
 *   {
 * @endcode
 * 
 * In order to convince the reader that these tools are indeed useful in
 * practice, let us choose a function for which it is not too difficult to
 * compute the analytical derivatives by hand. It's just sufficiently
 * complicated to make you think about whether or not you truly want to go
 * through with this exercise, and might also make you question whether you
 * are completely sure that your calculations and implementation for its
 * derivatives are correct. The point, of course, is that differentiation of
 * functions is in a sense relatively formulaic and should be something
 * computers are good at -- if we could build on existing software that
 * understands the rules, we wouldn't have to bother with doing it
 * ourselves.
 *     

 * 
 * We choose the two variable trigonometric function
 * $f(x,y) = \cos\left(\frac{y}{x}\right)$ for this purpose. Notice that
 * this function is templated on the number type. This is done because we
 * can often (but not always) use special auto-differentiable and symbolic
 * types as drop-in replacements for real or complex valued types, and these
 * will then perform some elementary calculations, such as evaluate a
 * function value along with its derivatives. We will exploit that property
 * and make sure that we need only define our function once, and then it can
 * be re-used in whichever context we wish to perform differential
 * operations on it.
 * 
 * @code
 *     template <typename NumberType>
 *     NumberType f(const NumberType &x, const NumberType &y)
 *     {
 *       return std::cos(y / x);
 *     }
 * 
 * @endcode
 * 
 * Rather than revealing this function's derivatives immediately, we'll
 * forward declare functions that return them and defer their definition to
 * later. As implied by the function names, they respectively return
 * the derivatives $\frac{df(x,y)}{dx}$:
 * 
 * @code
 *     double df_dx(const double x, const double y);
 * 
 * @endcode
 * 
 * $\frac{df(x,y)}{dy}$:
 * 
 * @code
 *     double df_dy(const double x, const double y);
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dx^{2}}$:
 * 
 * @code
 *     double d2f_dx_dx(const double x, const double y);
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dx dy}$:
 * 
 * @code
 *     double d2f_dx_dy(const double x, const double y);
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dy dx}$:
 * 
 * @code
 *     double d2f_dy_dx(const double x, const double y);
 * 
 * @endcode
 * 
 * and, lastly, $\frac{d^{2}f(x,y)}{dy^{2}}$:
 * 
 * @code
 *     double d2f_dy_dy(const double x, const double y);
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingderivativesusingautomaticdifferentiation"></a> 
 * <h4>Computing derivatives using automatic differentiation</h4>
 * 

 * 
 * To begin, we'll use AD as the tool to automatically
 * compute derivatives for us. We will evaluate the function with the
 * arguments `x` and `y`, and expect the resulting value and all of the
 * derivatives to match to within the given tolerance.
 * 
 * @code
 *     void
 *     run_and_verify_ad(const double x, const double y, const double tol = 1e-12)
 *     {
 * @endcode
 * 
 * Our function $f(x,y)$ is a scalar-valued function, with arguments that
 * represent the typical input variables that one comes across in
 * algebraic calculations or tensor calculus. For this reason, the
 * Differentiation::AD::ScalarFunction class is the appropriate wrapper
 * class to use to do the computations that we require. (As a point of
 * comparison, if the function arguments represented finite element cell
 * degrees-of-freedom, we'd want to treat them differently.) The spatial
 * dimension of the problem is irrelevant since we have no vector- or
 * tensor-valued arguments to accommodate, so the `dim` template argument
 * is arbitrarily assigned a value of 1. The second template argument
 * stipulates which AD framework will be used (deal.II has support for
 * several external AD frameworks), and what the underlying number type
 * provided by this framework is to be used. This number type
 * influences the maximum order of the differential operation, and the
 * underlying algorithms that are used to compute them. Given its template
 * nature, this choice is a compile-time decision because many (but not
 * all) of the AD libraries exploit compile-time meta-programming to
 * implement these special number types in an efficient manner. The third
 * template parameter states what the result type is; in our case, we're
 * working with `double`s.
 * 
 * @code
 *       constexpr unsigned int                     dim = 1;
 *       constexpr Differentiation::AD::NumberTypes ADTypeCode =
 *         Differentiation::AD::NumberTypes::sacado_dfad_dfad;
 *       using ADHelper =
 *         Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;
 * 
 * @endcode
 * 
 * It is necessary that we pre-register with our @p ADHelper class how many
 * arguments (what we will call "independent variables") the function
 * $f(x,y)$ has. Those arguments are `x` and `y`, so obviously there
 * are two of them.
 * 
 * @code
 *       constexpr unsigned int n_independent_variables = 2;
 * 
 * @endcode
 * 
 * We now have sufficient information to create and initialize an
 * instance of the helper class. We can also get the concrete
 * number type that will be used in all subsequent calculations.
 * This is useful, because we can write everything from here on by
 * referencing this type, and if we ever want to change the framework
 * used, or number type (e.g., if we need more differential operations)
 * then we need only adjust the `ADTypeCode` template parameter.
 * 
 * @code
 *       ADHelper ad_helper(n_independent_variables);
 *       using ADNumberType = typename ADHelper::ad_type;
 * 
 * @endcode
 * 
 * The next step is to register the numerical values of the independent
 * variables with the helper class. This is done because the function
 * and its derivatives will be evaluated for exactly these arguments.
 * Since we register them in the order `{x,y}`, the variable `x` will
 * be assigned component number `0`, and `y` will be component `1`
 * -- a detail that will be used in the next few lines.
 * 
 * @code
 *       ad_helper.register_independent_variables({x, y});
 * 
 * @endcode
 * 
 * We now ask for the helper class to give to us the independent variables
 * with their auto-differentiable representation. These are termed
 * "sensitive variables", because from this point on any operations that
 * we do with the components `independent_variables_ad` are tracked and
 * recorded by the AD framework, and will be considered
 * when we ask for the derivatives of something that they're used to
 * compute. What the helper returns is a `vector` of auto-differentiable
 * numbers, but we can be sure that the zeroth element represents `x`
 * and the first element `y`. Just to make completely sure that there's
 * no ambiguity of what number type these variables are, we suffix all of
 * the auto-differentiable variables with `ad`.
 * 
 * @code
 *       const std::vector<ADNumberType> independent_variables_ad =
 *         ad_helper.get_sensitive_variables();
 *       const ADNumberType &x_ad = independent_variables_ad[0];
 *       const ADNumberType &y_ad = independent_variables_ad[1];
 * 
 * @endcode
 * 
 * We can immediately pass in our sensitive representation of the
 * independent variables to our templated function that computes
 * $f(x,y)$.
 * This also returns an auto-differentiable number.
 * 
 * @code
 *       const ADNumberType f_ad = f(x_ad, y_ad);
 * 
 * @endcode
 * 
 * So now the natural question to ask is what we have actually just
 * computed by passing these special `x_ad` and `y_ad` variables to the
 * function `f`, instead of the original `double` variables `x` and `y`?
 * In other words, how is all of this related to the computation of the
 * derivatives that we were wanting to determine? Or, more concisely: What
 * is so special about this returned `ADNumberType` object that gives it
 * the ability to magically return derivatives?
 *       

 * 
 * In essence, how this *could* be done is the following:
 * This special number can be viewed as a data structure that stores the
 * function value, and the prescribed number of derivatives. For a
 * once-differentiable number expecting two arguments, it might look like
 * this:
 *       

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * struct ADNumberType
 * {
 *   double value;          // The value of the object
 *   double derivatives[2]; // Array of derivatives of the object with
 *                          // respect to x and y
 * };
 * @endcode
 * </div>
 *       

 * 
 * For our independent variable `x_ad`, the starting value of `x_ad.value`
 * would simply be its assigned value (i.e., the real value of that this
 * variable represents). The derivative `x_ad.derivatives[0]` would be
 * initialized to `1`, since `x` is the zeroth independent variable and
 * $\frac{d(x)}{dx} = 1$. The derivative `x.derivatives[1]` would be
 * initialized to zero, since the first independent variable is `y` and
 * $\frac{d(x)}{dy} = 0$.
 *       

 * 
 * For the function derivatives to be meaningful, we must assume that not
 * only is this function differentiable in an analytical sense, but that
 * it is also differentiable at the evaluation point `x,y`.
 * We can exploit both of these assumptions: when we use this number type
 * in mathematical operations, the AD framework *could*
 * overload the operations (e.g., `%operator+()`, `%operator*()` as well
 * as `%sin()`, `%exp()`, etc.) such that the returned result has the
 * expected value. At the same time, it would then compute the derivatives
 * through the knowledge of exactly what function is being overloaded and
 * rigorous application of the chain-rule. So, the `%sin()` function
 * (with its argument `a` itself being a function of the independent
 * variables `x` and `y`) *might* be defined as follows:
 *       

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * ADNumberType sin(const ADNumberType &a)
 * {
 *   ADNumberType output;
 *      

 *   // For the input argument "a", "a.value" is simply its value.
 *   output.value = sin(a.value);
 *      

 *   // We know that the derivative of sin(a) is cos(a), but we need
 *   // to also consider the chain rule and that the input argument
 *   // `a` is also differentiable with respect to the original
 *   // independent variables `x` and `y`. So `a.derivatives[0]`
 *   // and `a.derivatives[1]` respectively represent the partial
 *   // derivatives of `a` with respect to its inputs `x` and `y`.
 *   output.derivatives[0] = cos(a.value)*a.derivatives[0];
 *   output.derivatives[1] = cos(a.value)*a.derivatives[1];
 *      

 *   return output;
 * }
 * @endcode
 * </div>
 *       

 * 
 * All of that could of course also be done for second and even higher
 * order derivatives.
 *       

 * 
 * So it is now clear that with the above representation the
 * `ADNumberType` is carrying around some extra data that represents the
 * various derivatives of differentiable functions with respect to the
 * original (sensitive) independent variables. It should therefore be
 * noted that there is computational overhead associated with using them
 * (as we compute extra functions when doing derivative computations) as
 * well as memory overhead in storing these results. So the prescribed
 * number of levels of differential operations should ideally be kept to a
 * minimum to limit computational cost. We could, for instance, have
 * computed the first derivatives ourself and then have used the
 * Differentiation::AD::VectorFunction helper class to determine the
 * gradient of the collection of dependent functions, which would be the
 * second derivatives of the original scalar function.
 *       

 * 
 * It is also worth noting that because the chain rule is indiscriminately
 * applied and we only see the beginning and end-points of the calculation
 * `{x,y}` $\rightarrow$ `f(x,y)`, we will only ever be able to query
 * the total derivatives of `f`; the partial derivatives
 * (`a.derivatives[0]` and `a.derivatives[1]` in the above example) are
 * intermediate values and are hidden from us.
 * 

 * 
 * Okay, since we now at least have some idea as to exactly what `f_ad`
 * represents and what is encoded within it, let's put all of that to
 * some actual use. To gain access to those hidden derivative results,
 * we register the final result with the helper class. After this point,
 * we can no longer change the value of `f_ad` and have those changes
 * reflected in the results returned by the helper class.
 * 
 * @code
 *       ad_helper.register_dependent_variable(f_ad);
 * 
 * @endcode
 * 
 * The next step is to extract the derivatives (specifically, the function
 * gradient and Hessian). To do so we first create some temporary data
 * structures (with the result type `double`) to store the derivatives
 * (noting that all derivatives are returned at once, and not
 * individually)...
 * 
 * @code
 *       Vector<double>     Df(ad_helper.n_dependent_variables());
 *       FullMatrix<double> D2f(ad_helper.n_dependent_variables(),
 *                              ad_helper.n_independent_variables());
 * 
 * @endcode
 * 
 * ... and we then request that the helper class compute these
 * derivatives, and the function value itself. And that's it. We have
 * everything that we were aiming to get.
 * 
 * @code
 *       const double computed_f = ad_helper.compute_value();
 *       ad_helper.compute_gradient(Df);
 *       ad_helper.compute_hessian(D2f);
 * 
 * @endcode
 * 
 * We can convince ourselves that the AD framework is
 * correct by comparing it to the analytical solution. (Or, if you're
 * like the author, you'll be doing the opposite and will rather verify
 * that your implementation of the analytical solution is correct!)
 * 
 * @code
 *       AssertThrow(std::abs(f(x, y) - computed_f) < tol,
 *                   ExcMessage(std::string("Incorrect value computed for f. ") +
 *                              std::string("Hand-calculated value: ") +
 *                              Utilities::to_string(f(x, y)) +
 *                              std::string(" ; ") +
 *                              std::string("Value computed by AD: ") +
 *                              Utilities::to_string(computed_f)));
 * 
 * @endcode
 * 
 * Because we know the ordering of the independent variables, we know
 * which component of the gradient relates to which derivative...
 * 
 * @code
 *       const double computed_df_dx = Df[0];
 *       const double computed_df_dy = Df[1];
 * 
 *       AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for df/dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_df_dx)));
 *       AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for df/dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(df_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_df_dy)));
 * 
 * @endcode
 * 
 * ... and similar for the Hessian.
 * 
 * @code
 *       const double computed_d2f_dx_dx = D2f[0][0];
 *       const double computed_d2f_dx_dy = D2f[0][1];
 *       const double computed_d2f_dy_dx = D2f[1][0];
 *       const double computed_d2f_dy_dy = D2f[1][1];
 * 
 *       AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dx_dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_d2f_dx_dx)));
 *       AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dx_dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_d2f_dx_dy)));
 *       AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dy_dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_d2f_dy_dx)));
 *       AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dy_dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_d2f_dy_dy)));
 *     }
 * 
 * @endcode
 * 
 * That's pretty great. There wasn't too much work involved in computing
 * second-order derivatives of this trigonometric function.
 * 

 * 
 * 
 * <a name="Handcalculatedderivativesoftheanalyticalsolution"></a> 
 * <h4>Hand-calculated derivatives of the analytical solution</h4>
 * 

 * 
 * Since we now know how much "implementation effort" it takes to have the
 * AD framework compute those derivatives for us, let's
 * compare that to the same computed by hand and implemented in several
 * stand-alone functions.
 * 

 * 
 * Here are the two first derivatives of $f(x,y) =
 * \cos\left(\frac{y}{x}\right)$:
 *     

 * 
 * $\frac{df(x,y)}{dx} = \frac{y}{x^2} \sin\left(\frac{y}{x}\right)$
 * 
 * @code
 *     double df_dx(const double x, const double y)
 *     {
 *       Assert(x != 0.0, ExcDivideByZero());
 *       return y * std::sin(y / x) / (x * x);
 *     }
 * 
 * @endcode
 * 
 * $\frac{df(x,y)}{dx} = -\frac{1}{x} \sin\left(\frac{y}{x}\right)$
 * 
 * @code
 *     double df_dy(const double x, const double y)
 *     {
 *       return -std::sin(y / x) / x;
 *     }
 * 
 * @endcode
 * 
 * And here are the four second derivatives of $f(x,y)$:
 *     

 * 
 * $\frac{d^{2}f(x,y)}{dx^{2}} = -\frac{y}{x^4} (2x
 * \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
 * 
 * @code
 *     double d2f_dx_dx(const double x, const double y)
 *     {
 *       return -y * (2 * x * std::sin(y / x) + y * std::cos(y / x)) /
 *              (x * x * x * x);
 *     }
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dx dy} = \frac{1}{x^3} (x
 * \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$
 * 
 * @code
 *     double d2f_dx_dy(const double x, const double y)
 *     {
 *       return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
 *     }
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dy dx} = \frac{1}{x^3} (x
 * \sin\left(\frac{y}{x}\right) + y \cos\left(\frac{y}{x}\right))$ (as
 * expected, on the basis of [Schwarz's
 * theorem](https://en.wikipedia.org/wiki/Symmetry_of_second_derivatives))
 * 
 * @code
 *     double d2f_dy_dx(const double x, const double y)
 *     {
 *       return (x * std::sin(y / x) + y * std::cos(y / x)) / (x * x * x);
 *     }
 * 
 * @endcode
 * 
 * $\frac{d^{2}f(x,y)}{dy^{2}} = -\frac{1}{x^2}
 * \cos\left(\frac{y}{x}\right)$
 * 
 * @code
 *     double d2f_dy_dy(const double x, const double y)
 *     {
 *       return -(std::cos(y / x)) / (x * x);
 *     }
 * 
 * @endcode
 * 
 * Hmm... there's a lot of places in the above where we could have
 * introduced an error in the above, especially when it comes to applying
 * the chain rule. Although they're no silver bullet, at the very least
 * these AD frameworks can serve as a verification tool to make sure that we
 * haven't made any errors (either by calculation or by implementation) that
 * would negatively affect our results.
 * 

 * 
 * The point of this example of course is that we might have
 * chosen a relatively simple function $f(x,y)$ for which we can
 * hand-verify that the derivatives the AD framework computed is
 * correct. But the AD framework didn't care that the function was
 * simple: It could have been a much much more convoluted
 * expression, or could have depended on more than two variables,
 * and it would still have been able to compute the derivatives --
 * the only difference would have been that *we* wouldn't have
 * been able to come up with the derivatives any more to verify
 * correctness of the AD framework.
 * 

 * 
 * 

 * 
 * 

 * 
 * 
 * <a name="Computingderivativesusingsymbolicdifferentiation"></a> 
 * <h4>Computing derivatives using symbolic differentiation</h4>
 * 

 * 
 * We'll now repeat the same exercise using symbolic differentiation. The
 * term "symbolic differentiation" is a little bit misleading because
 * differentiation is just one tool that the Computer Algebra System (CAS)
 * (i.e., the symbolic framework) provides. Nevertheless, in the context
 * of finite element modeling and applications it is the most common use
 * of a CAS and will therefore be the one that we'll focus on.
 * Once more, we'll supply the argument values `x` and `y` with which to
 * evaluate our function $f(x,y) = \cos\left(\frac{y}{x}\right)$ and its
 * derivatives, and a tolerance with which to test the correctness of the
 * returned results.
 * 
 * @code
 *     void
 *     run_and_verify_sd(const double x, const double y, const double tol = 1e-12)
 *     {
 * @endcode
 * 
 * The first step that we need to take is to form the symbolic variables
 * that represent the function arguments that we wish to differentiate
 * with respect to. Again, these will be the independent variables for
 * our problem and as such are, in some sense, primitive variables that
 * have no dependencies on any other variable. We create these types of
 * (independent) variables by initializing a symbolic type
 * Differentiation::SD::Expression, which is a wrapper to a set of classes
 * used by the symbolic framework, with a unique identifier. On this
 * occasion it makes sense that this identifier, a `std::string`, be
 * simply `"x"` for the $x$ argument, and likewise `"y"` for the $y$
 * argument to the dependent function. Like before, we'll suffix symbolic
 * variable names with `sd` so that we can clearly see which variables are
 * symbolic (as opposed to numeric) in nature.
 * 
 * @code
 *       const Differentiation::SD::Expression x_sd("x");
 *       const Differentiation::SD::Expression y_sd("y");
 * 
 * @endcode
 * 
 * Using the templated function that computes $f(x,y)$, we can pass
 * these independent variables as arguments to the function. The returned
 * result will be another symbolic type that represents the sequence of
 * operations used to compute $\cos\left(\frac{y}{x}\right)$.
 * 
 * @code
 *       const Differentiation::SD::Expression f_sd = f(x_sd, y_sd);
 * 
 * @endcode
 * 
 * At this point it is legitimate to print out the expression `f_sd`, and
 * if we did so
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * std::cout << "f(x,y) = " << f_sd << std::endl;
 * @endcode
 * </div>
 * we would see `f(x,y) = cos(y/x)` printed to the console.
 *       

 * 
 * You might notice that we've constructed our symbolic function `f_sd`
 * with no context as to how we might want to use it: In contrast to the
 * AD approach shown above, what we were returned from calling
 * `f(x_sd, y_sd)` is not the evaluation of the function `f` at some
 * specific point, but is in fact a symbolic representation of the
 * evaluation at a generic, as yet undetermined, point. This is one of the
 * key points that makes symbolic frameworks (the CAS) different from
 * automatic differentiation frameworks. Each of the variables `x_sd` and
 * `y_sd`, and even the composite dependent function `f_sd`, are in some
 * sense respectively "placeholders" for numerical values and a
 * composition of operations. In fact, the individual components that are
 * used to compose the function are also placeholders. The sequence of
 * operations are encoded into in a tree-like data structure (conceptually
 * similar to an [abstract syntax
 * tree](https://en.wikipedia.org/wiki/Abstract_syntax_tree)).
 *       

 * 
 * Once we form these data structures we can defer any operations that we
 * might want to do with them until some later time. Each of these
 * placeholders represents something, but we have the opportunity to
 * define or redefine what they represent at any convenient point in time.
 * So for this particular problem it makes sense that we somehow want to
 * associate "x" and "y" with *some* numerical value (with type yet to be
 * determined), but we could conceptually (and if it made sense) assign
 * the ratio "y/x" a value instead of the variables "x" and "y"
 * individually. We could also associate with "x" or "y" some other
 * symbolic function `g(a,b)`. Any of these operations involves
 * manipulating the recorded tree of operations, and substituting the
 * salient nodes on the tree (and that nodes' subtree) with something
 * else. The key word here is "substitution", and indeed there are many
 * functions in the Differentiation::SD namespace that have this word
 * in their names.
 *       

 * 
 * This capability makes the framework entirely generic.
 * In the context of finite element simulations, the types of operations
 * that we would typically perform with our symbolic types are
 * function composition, differentiation, substitution (partial or
 * complete), and evaluation (i.e., conversion of the symbolic type to its
 * numerical counterpart). But should you need it, a CAS is often capable
 * of more than just this: It could be forming anti-derivatives
 * (integrals) of functions, perform simplifications on the expressions
 * that form a function (e.g., replace $(\sin a)^2 + (\cos a)^2$ by
 * $1$; or, more simply: if the function did an operation like `1+2`, a
 * CAS could replace it by `3`), and so forth: The *expression* that a
 * variable represents is obtained from how the function $f$ is
 * implemented, but a CAS can do with it whatever its functionality
 * happens to be.
 *       

 * 
 * Specifically, to compute the symbolic representation of the first
 * derivatives of the dependent function with respect to its individual
 * independent variables, we use the
 * Differentiation::SD::Expression::differentiate() function with the
 * independent variable given as its argument. Each call will cause the
 * CAS to go through the tree of operations that compose `f_sd` and
 * differentiate each node of the expression tree with respect to the
 * given symbolic argument.
 * 
 * @code
 *       const Differentiation::SD::Expression df_dx_sd = f_sd.differentiate(x_sd);
 *       const Differentiation::SD::Expression df_dy_sd = f_sd.differentiate(y_sd);
 * 
 * @endcode
 * 
 * To compute the symbolic representation of the second derivatives, we
 * simply differentiate the first derivatives with respect to the
 * independent variables. So to compute a higher order derivative, we
 * first need to compute the lower order derivative.
 * (As the return type of the call to `differentiate()` is an expression,
 * we could in principal execute double differentiation directly from the
 * scalar by chaining two calls together. But this is unnecessary in this
 * particular case, since we have the intermediate results at hand.)
 * 
 * @code
 *       const Differentiation::SD::Expression d2f_dx_dx_sd =
 *         df_dx_sd.differentiate(x_sd);
 *       const Differentiation::SD::Expression d2f_dx_dy_sd =
 *         df_dx_sd.differentiate(y_sd);
 *       const Differentiation::SD::Expression d2f_dy_dx_sd =
 *         df_dy_sd.differentiate(x_sd);
 *       const Differentiation::SD::Expression d2f_dy_dy_sd =
 *         df_dy_sd.differentiate(y_sd);
 * @endcode
 * 
 * Printing the expressions for the first and second derivatives, as
 * computed by the CAS, using the statements
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * std::cout << "df_dx_sd: " << df_dx_sd << std::endl;
 * std::cout << "df_dy_sd: " << df_dy_sd << std::endl;
 * std::cout << "d2f_dx_dx_sd: " << d2f_dx_dx_sd << std::endl;
 * std::cout << "d2f_dx_dy_sd: " << d2f_dx_dy_sd << std::endl;
 * std::cout << "d2f_dy_dx_sd: " << d2f_dy_dx_sd << std::endl;
 * std::cout << "d2f_dy_dy_sd: " << d2f_dy_dy_sd << std::endl;
 * @endcode
 * </div>
 * renders the following output:
 * <div class=CodeFragmentInTutorialComment>
 * @code{.sh}
 * df_dx_sd: y*sin(y/x)/x**2
 * df_dy_sd: -sin(y/x)/x
 * d2f_dx_dx_sd: -y**2*cos(y/x)/x**4 - 2*y*sin(y/x)/x**3
 * d2f_dx_dy_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
 * d2f_dy_dx_sd: sin(y/x)/x**2 + y*cos(y/x)/x**3
 * d2f_dy_dy_sd: -cos(y/x)/x**2
 * @endcode
 * </div>
 * This compares favorably to the analytical expressions for these
 * derivatives that were presented earlier.
 * 

 * 
 * Now that we have formed the symbolic expressions for the function and
 * its derivatives, we want to evaluate them for the numeric values for
 * the main function arguments `x` and `y`. To accomplish this, we
 * construct a *substitution map*, which maps the symbolic values to their
 * numerical counterparts.
 * 
 * @code
 *       const Differentiation::SD::types::substitution_map substitution_map =
 *         Differentiation::SD::make_substitution_map(
 *           std::pair<Differentiation::SD::Expression, double>{x_sd, x},
 *           std::pair<Differentiation::SD::Expression, double>{y_sd, y});
 * 
 * @endcode
 * 
 * The last step in the process is to convert all symbolic variables and
 * operations into numerical values, and produce the numerical result of
 * this operation. To do this we combine the substitution map with the
 * symbolic variable in the step we have already mentioned above:
 * "substitution".
 *       

 * 
 * Once we pass this substitution map to the CAS, it will
 * substitute each instance of the symbolic variable (or, more generally,
 * sub-expression) with its numerical counterpart and then propagate these
 * results up the operation tree, simplifying each node on the tree if
 * possible. If the tree is reduced to a single value (i.e., we have
 * substituted all of the independent variables with their numerical
 * counterpart) then the evaluation is complete.
 *       

 * 
 * Due to the strongly-typed nature of C++, we need to instruct the CAS to
 * convert its representation of the result into an intrinsic data type
 * (in this case a `double`). This is the "evaluation" step, and through
 * the template type we define the return type of this process.
 * Conveniently, these two steps can be done at once if we are certain
 * that we've performed a full substitution.
 * 
 * @code
 *       const double computed_f =
 *         f_sd.substitute_and_evaluate<double>(substitution_map);
 * 
 *       AssertThrow(std::abs(f(x, y) - computed_f) < tol,
 *                   ExcMessage(std::string("Incorrect value computed for f. ") +
 *                              std::string("Hand-calculated value: ") +
 *                              Utilities::to_string(f(x, y)) +
 *                              std::string(" ; ") +
 *                              std::string("Value computed by AD: ") +
 *                              Utilities::to_string(computed_f)));
 * 
 * @endcode
 * 
 * We can do the same for the first derivatives...
 * 
 * @code
 *       const double computed_df_dx =
 *         df_dx_sd.substitute_and_evaluate<double>(substitution_map);
 *       const double computed_df_dy =
 *         df_dy_sd.substitute_and_evaluate<double>(substitution_map);
 * 
 *       AssertThrow(std::abs(df_dx(x, y) - computed_df_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for df/dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(df_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_df_dx)));
 *       AssertThrow(std::abs(df_dy(x, y) - computed_df_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for df/dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(df_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by AD: ") +
 *                     Utilities::to_string(computed_df_dy)));
 * 
 * @endcode
 * 
 * ... and the second derivatives.
 * Notice that we can reuse the same substitution map for each of these
 * operations because we wish to evaluate all of these functions for the
 * same values of `x` and `y`. Modifying the values in the substitution
 * map renders the result of same symbolic expression evaluated with
 * different values being assigned to the independent variables.
 * We could also happily have each variable represent a real value in
 * one pass, and a complex value in the next.
 * 
 * @code
 *       const double computed_d2f_dx_dx =
 *         d2f_dx_dx_sd.substitute_and_evaluate<double>(substitution_map);
 *       const double computed_d2f_dx_dy =
 *         d2f_dx_dy_sd.substitute_and_evaluate<double>(substitution_map);
 *       const double computed_d2f_dy_dx =
 *         d2f_dy_dx_sd.substitute_and_evaluate<double>(substitution_map);
 *       const double computed_d2f_dy_dy =
 *         d2f_dy_dy_sd.substitute_and_evaluate<double>(substitution_map);
 * 
 *       AssertThrow(std::abs(d2f_dx_dx(x, y) - computed_d2f_dx_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dx_dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dx_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by SD: ") +
 *                     Utilities::to_string(computed_d2f_dx_dx)));
 *       AssertThrow(std::abs(d2f_dx_dy(x, y) - computed_d2f_dx_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dx_dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dx_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by SD: ") +
 *                     Utilities::to_string(computed_d2f_dx_dy)));
 *       AssertThrow(std::abs(d2f_dy_dx(x, y) - computed_d2f_dy_dx) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dy_dx. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dy_dx(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by SD: ") +
 *                     Utilities::to_string(computed_d2f_dy_dx)));
 *       AssertThrow(std::abs(d2f_dy_dy(x, y) - computed_d2f_dy_dy) < tol,
 *                   ExcMessage(
 *                     std::string("Incorrect value computed for d2f/dy_dy. ") +
 *                     std::string("Hand-calculated value: ") +
 *                     Utilities::to_string(d2f_dy_dy(x, y)) + std::string(" ; ") +
 *                     std::string("Value computed by SD: ") +
 *                     Utilities::to_string(computed_d2f_dy_dy)));
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheSimpleExamplerunfunction"></a> 
 * <h4>The SimpleExample::run() function</h4>
 * 

 * 
 * The function used to drive these initial examples is straightforward.
 * We'll arbitrarily choose some values at which to evaluate the function
 * (although knowing that `x = 0` is not permissible), and then pass these
 * values to the functions that use the AD and SD frameworks.
 * 
 * @code
 *     void run()
 *     {
 *       const double x = 1.23;
 *       const double y = 0.91;
 * 
 *       std::cout << "Simple example using automatic differentiation..."
 *                 << std::endl;
 *       run_and_verify_ad(x, y);
 *       std::cout << "... all calculations are correct!" << std::endl;
 * 
 *       std::cout << "Simple example using symbolic differentiation."
 *                 << std::endl;
 *       run_and_verify_sd(x, y);
 *       std::cout << "... all calculations are correct!" << std::endl;
 *     }
 * 
 *   } // namespace SimpleExample
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="AmorecomplexexampleUsingautomaticandsymbolicdifferentiationtocomputederivativesatcontinuumpoints"></a> 
 * <h3>A more complex example: Using automatic and symbolic differentiation to compute derivatives at continuum points</h3>
 * 

 * 
 * Now that we've introduced the principles behind automatic and symbolic
 * differentiation, we'll put them into action by formulating two coupled
 * magneto-mechanical constitutive laws: one that is rate-independent, and
 * another that exhibits rate-dependent behavior.
 *   

 * 
 * As you will recall from the introduction, the material
 * constitutive laws we will consider are far more complicated than
 * the simple example above. This is not just because of the form of
 * the function $\psi_{0}$ that we will consider, but in particular
 * because $\psi_{0}$ doesn't just depend on two scalar variables, but
 * instead on a whole bunch of *tensors*, each with several
 * components. In some cases, these are *symmetric* tensors, for
 * which only a subset of components is in fact independent, and one has
 * to think about what it actually means to compute a derivative
 * such as $\frac{\partial\psi_{0}}{\partial \mathbf{C}}$ where $\mathbf
 * C$ is a symmetric tensor. How all of this will work will,
 * hopefully, become clear below. It will also become clear that
 * doing this by hand is going to be, at the very best, *exceedingly*
 * *tedious* and, at worst, riddled with hard-to-find bugs.
 * 
 * @code
 *   namespace CoupledConstitutiveLaws
 *   {
 * @endcode
 * 
 * 
 * <a name="Constitutiveparameters"></a> 
 * <h4>Constitutive parameters</h4>
 * 

 * 
 * We start with a description of the various material parameters
 * that appear in the description of the energy function $\psi_{0}$.
 *     

 * 
 * The ConstitutiveParameters class is used to hold these values.
 * Values for all parameters (both constitutive and rheological) are taken
 * from @cite Pelteret2018a, and are given values that produce a
 * constitutive response that is broadly representative of a real,
 * laboratory-made magneto-active polymer, though the specific values used
 * here are of no consequence to the purpose of this program of course.
 *     

 * 
 * The first four constitutive parameters respectively represent
 * - the elastic shear modulus $\mu_{e}$,
 * - the elastic shear modulus at magnetic saturation $\mu_{e}^{\infty}$,
 * - the saturation magnetic field strength for the elastic shear
 * modulus $h_{e}^{\text{sat}}$, and
 * - the Poisson ratio $\nu$.
 * 
 * @code
 *     class ConstitutiveParameters : public ParameterAcceptor
 *     {
 *     public:
 *       ConstitutiveParameters();
 * 
 *       double mu_e       = 30.0e3;
 *       double mu_e_inf   = 250.0e3;
 *       double mu_e_h_sat = 212.2e3;
 *       double nu_e       = 0.49;
 * 
 * @endcode
 * 
 * The next four, which only pertain to the rate-dependent material, are
 * parameters for
 * - the viscoelastic shear modulus $\mu_{v}$,
 * - the viscoelastic shear modulus at magnetic saturation
 * $\mu_{v}^{\infty}$,
 * - the saturation magnetic field strength for the viscoelastic
 * shear modulus $h_{v}^{\text{sat}}$, and
 * - the characteristic relaxation time $\tau$.
 * 
 * @code
 *       double mu_v       = 20.0e3;
 *       double mu_v_inf   = 35.0e3;
 *       double mu_v_h_sat = 92.84e3;
 *       double tau_v      = 0.6;
 * 
 * @endcode
 * 
 * The last parameter is the relative magnetic permeability $\mu_{r}$.
 * 
 * @code
 *       double mu_r = 6.0;
 * 
 *       bool initialized = false;
 *     };
 * 
 * @endcode
 * 
 * The parameters are initialized through the ParameterAcceptor
 * framework, which is discussed in detail in step-60.
 * 
 * @code
 *     ConstitutiveParameters::ConstitutiveParameters()
 *       : ParameterAcceptor("/Coupled Constitutive Laws/Constitutive Parameters/")
 *     {
 *       add_parameter("Elastic shear modulus", mu_e);
 *       add_parameter("Elastic shear modulus at magnetic saturation", mu_e_inf);
 *       add_parameter(
 *         "Saturation magnetic field strength for elastic shear modulus",
 *         mu_e_h_sat);
 *       add_parameter("Poisson ratio", nu_e);
 * 
 *       add_parameter("Viscoelastic shear modulus", mu_v);
 *       add_parameter("Viscoelastic shear modulus at magnetic saturation",
 *                     mu_v_inf);
 *       add_parameter(
 *         "Saturation magnetic field strength for viscoelastic shear modulus",
 *         mu_v_h_sat);
 *       add_parameter("Characteristic relaxation time", tau_v);
 * 
 *       add_parameter("Relative magnetic permeability", mu_r);
 * 
 *       parse_parameters_call_back.connect([&]() { initialized = true; });
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ConstitutivelawsBaseclass"></a> 
 * <h4>Constitutive laws: Base class</h4>
 * 

 * 
 * Since we'll be formulating two constitutive laws for the same class of
 * materials, it makes sense to define a base class that ensures a unified
 * interface to them.
 *     

 * 
 * The class declaration starts with the constructor that will
 * accept the set of constitutive parameters that, in conjunction
 * with the material law itself, dictate the material response.
 * 
 * @code
 *     template <int dim>
 *     class Coupled_Magnetomechanical_Constitutive_Law_Base
 *     {
 *     public:
 *       Coupled_Magnetomechanical_Constitutive_Law_Base(
 *         const ConstitutiveParameters &constitutive_parameters);
 * 
 * @endcode
 * 
 * Instead of computing and returning the kinetic variables or their
 * linearization at will, we'll calculate and store these values within a
 * single method. These cached results will then be returned upon request.
 * We'll defer the precise explanation as to why we'd want to do this to
 * a later stage. What is important for now is to see that this function
 * accepts all of the field variables, namely the magnetic field vector
 * $\boldsymbol{\mathbb{H}}$ and right Cauchy-Green deformation tensor
 * $\mathbf{C}$, as well as the time discretizer. These, in addition to
 * the @p constitutive_parameters, are all the fundamental quantities that
 * are required to compute the material response.
 * 
 * @code
 *       virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
 *                                         const Tensor<1, dim> &         H,
 *                                         const DiscreteTime &time) = 0;
 * 
 * @endcode
 * 
 * The next few functions provide the interface to probe the material
 * response due subject to the applied deformation and magnetic loading.
 *       

 * 
 * Since the class of materials can be expressed in terms of a free energy
 * $\psi_{0}$, we can compute that...
 * 
 * @code
 *       virtual double get_psi() const = 0;
 * 
 * @endcode
 * 
 * ... as well as the two kinetic quantities:
 * - the magnetic induction vector $\boldsymbol{\mathbb{B}}$, and
 * - the total Piola-Kirchhoff stress tensor $\mathbf{S}^{\text{tot}}$
 * 
 * @code
 *       virtual Tensor<1, dim> get_B() const = 0;
 * 
 *       virtual SymmetricTensor<2, dim> get_S() const = 0;
 * 
 * @endcode
 * 
 * ... and the linearization of the kinetic quantities, which are:
 * - the magnetostatic tangent tensor $\mathbb{D}$,
 * - the total referential magnetoelastic coupling tensor
 * $\mathfrak{P}^{\text{tot}}$, and
 * - the total referential elastic tangent tensor
 * $\mathcal{H}^{\text{tot}}$.
 * 
 * @code
 *       virtual SymmetricTensor<2, dim> get_DD() const = 0;
 * 
 *       virtual Tensor<3, dim> get_PP() const = 0;
 * 
 *       virtual SymmetricTensor<4, dim> get_HH() const = 0;
 * 
 * @endcode
 * 
 * We'll also define a method that provides a mechanism for this class
 * instance to do any additional tasks before moving on to the next
 * timestep. Again, the reason for doing this will become clear a little
 * later.
 * 
 * @code
 *       virtual void update_end_of_timestep()
 *       {}
 * 
 * @endcode
 * 
 * In the `protected` part of the class,
 * we store a reference to an instance of the constitutive parameters
 * that govern the material response.
 * For convenience, we also define some functions that return
 * various constitutive parameters (both explicitly defined, as well
 * as calculated).
 *       

 * 
 * The parameters related to the elastic response of the material are,
 * in order:
 * - the elastic shear modulus,
 * - the elastic shear modulus at saturation magnetic field,
 * - the saturation magnetic field strength for the elastic shear
 * modulus,
 * - the Poisson ratio,
 * - the Lam&eacute; parameter, and
 * - the bulk modulus.
 * 
 * @code
 *     protected:
 *       const ConstitutiveParameters &constitutive_parameters;
 * 
 *       double get_mu_e() const;
 * 
 *       double get_mu_e_inf() const;
 * 
 *       double get_mu_e_h_sat() const;
 * 
 *       double get_nu_e() const;
 * 
 *       double get_lambda_e() const;
 * 
 *       double get_kappa_e() const;
 * 
 * @endcode
 * 
 * The parameters related to the elastic response of the material are,
 * in order:
 * - the viscoelastic shear modulus,
 * - the viscoelastic shear modulus at magnetic saturation,
 * - the saturation magnetic field strength for the viscoelastic
 * shear modulus, and
 * - the characteristic relaxation time.
 * 
 * @code
 *       double get_mu_v() const;
 * 
 *       double get_mu_v_inf() const;
 * 
 *       double get_mu_v_h_sat() const;
 * 
 *       double get_tau_v() const;
 * 
 * @endcode
 * 
 * The parameters related to the magnetic response of the material are,
 * in order:
 * - the relative magnetic permeability, and
 * - the magnetic permeability constant $\mu_{0}$ (not really a material
 * constant,
 * but rather a universal constant that we'll group here for
 * simplicity).
 *       

 * 
 * We'll also implement a function that returns the
 * timestep size from the time discretizion.
 * 
 * @code
 *       double get_mu_r() const;
 * 
 *       constexpr double get_mu_0() const;
 *       double           get_delta_t(const DiscreteTime &time) const;
 *     };
 * 
 * 
 * 
 * @endcode
 * 
 * In the following, let us start by implementing the several
 * relatively trivial member functions of the class just defined:
 * 
 * @code
 *     template <int dim>
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::
 *       Coupled_Magnetomechanical_Constitutive_Law_Base(
 *         const ConstitutiveParameters &constitutive_parameters)
 *       : constitutive_parameters(constitutive_parameters)
 *     {
 *       Assert(get_kappa_e() > 0, ExcInternalError());
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e() const
 *     {
 *       return constitutive_parameters.mu_e;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_inf() const
 *     {
 *       return constitutive_parameters.mu_e_inf;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_e_h_sat() const
 *     {
 *       return constitutive_parameters.mu_e_h_sat;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_nu_e() const
 *     {
 *       return constitutive_parameters.nu_e;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_lambda_e() const
 *     {
 *       return 2.0 * get_mu_e() * get_nu_e() / (1.0 - 2.0 * get_nu_e());
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_kappa_e() const
 *     {
 *       return (2.0 * get_mu_e() * (1.0 + get_nu_e())) /
 *              (3.0 * (1.0 - 2.0 * get_nu_e()));
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v() const
 *     {
 *       return constitutive_parameters.mu_v;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_inf() const
 *     {
 *       return constitutive_parameters.mu_v_inf;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_v_h_sat() const
 *     {
 *       return constitutive_parameters.mu_v_h_sat;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_tau_v() const
 *     {
 *       return constitutive_parameters.tau_v;
 *     }
 * 
 * 
 *     template <int dim>
 *     double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_r() const
 *     {
 *       return constitutive_parameters.mu_r;
 *     }
 * 
 * 
 *     template <int dim>
 *     constexpr double
 *     Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_mu_0() const
 *     {
 *       return 4.0 * numbers::PI * 1e-7;
 *     }
 * 
 * 
 *     template <int dim>
 *     double Coupled_Magnetomechanical_Constitutive_Law_Base<dim>::get_delta_t(
 *       const DiscreteTime &time) const
 *     {
 *       return time.get_previous_step_size();
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Magnetoelasticconstitutivelawusingautomaticdifferentiation"></a> 
 * <h4>Magnetoelastic constitutive law (using automatic differentiation)</h4>
 * 

 * 
 * We'll begin by considering a non-dissipative material, namely one that
 * is governed by a magneto-hyperelastic constitutive law that exhibits
 * stiffening when immersed in a magnetic field. As described in
 * the introduction, the stored energy density function for such a material
 * might be given by
 * @f[
 * \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
 * \right]
 * + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
 * - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
 * \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]
 * @f]
 * with
 * @f[
 * f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right) .
 * @f]
 *     

 * 
 * Now on to the class that implements this behavior.
 * Since we expect that this class fully describes a single material, we'll
 * mark it as "final" so that the inheritance tree terminated here.
 * At the top of the class, we define the helper type that we will use in
 * the AD computations for our scalar energy density function. Note that we
 * expect it to return values of type `double`. We also have to specify the
 * number of spatial dimensions, `dim`, so that the link between vector,
 * tensor and symmetric tensor fields and the number of components that they
 * contain may be established. The concrete `ADTypeCode` used for the
 * ADHelper class will be provided as a template argument at the point where
 * this class is actually used.
 * 
 * @code
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     class Magnetoelastic_Constitutive_Law_AD final
 *       : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *     {
 *       using ADHelper =
 *         Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;
 *       using ADNumberType = typename ADHelper::ad_type;
 * 
 *     public:
 *       Magnetoelastic_Constitutive_Law_AD(
 *         const ConstitutiveParameters &constitutive_parameters);
 * 
 * @endcode
 * 
 * Since the public interface to the base class is pure-`virtual`, here
 * we'll declare that this class will override all of these base class
 * methods.
 * 
 * @code
 *       virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
 *                                         const Tensor<1, dim> &         H,
 *                                         const DiscreteTime &) override;
 * 
 *       virtual double get_psi() const override;
 * 
 *       virtual Tensor<1, dim> get_B() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_S() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_DD() const override;
 * 
 *       virtual Tensor<3, dim> get_PP() const override;
 * 
 *       virtual SymmetricTensor<4, dim> get_HH() const override;
 * 
 * @endcode
 * 
 * In the `private` part of the class,
 * we need to define some extractors that will help us set independent
 * variables and later get the computed values related to the dependent
 * variables. If this class were to be used in the context of a finite
 * element problem, then each of these extractors is (most likely) related
 * to the gradient of a component of the solution field (in this case,
 * displacement and magnetic scalar potential). As you can probably infer
 * by now, here "C" denotes the right Cauchy-Green tensor and "H" denotes
 * the magnetic field vector.
 * 
 * @code
 *     private:
 *       const FEValuesExtractors::Vector             H_components;
 *       const FEValuesExtractors::SymmetricTensor<2> C_components;
 * 
 * @endcode
 * 
 * This is an instance of the automatic differentiation helper that
 * we'll set up to do all of the differential calculations related to
 * the constitutive law...
 * 
 * @code
 *       ADHelper ad_helper;
 * 
 * @endcode
 * 
 * ... and the following three member variables will store the output from
 * the
 * @p ad_helper. The @p ad_helper returns the derivatives with respect
 * to all field variables at once, so we'll retain the full gradient
 * vector and Hessian matrix. From that, we'll extract the individual
 * entries that we're actually interested in.
 * 
 * @code
 *       double             psi;
 *       Vector<double>     Dpsi;
 *       FullMatrix<double> D2psi;
 *     };
 * 
 * @endcode
 * 
 * When setting up the field component extractors, it is completely
 * arbitrary as to how they are ordered. But it is important that the
 * extractors do not have overlapping indices. The total number of
 * components of these extractors defines the number of independent
 * variables that the
 * @p ad_helper needs to track, and with respect to which we'll be taking
 * derivatives. The resulting data structures @p Dpsi and @p D2psi must also
 * be sized accordingly. Once the @p ad_helper is configured (its input
 * argument being the total number of components of $\mathbf{C}$ and
 * $\boldsymbol{\mathbb{H}}$), we can directly interrogate it as to how many
 * independent variables it uses.
 * 
 * @code
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::
 *       Magnetoelastic_Constitutive_Law_AD(
 *         const ConstitutiveParameters &constitutive_parameters)
 *       : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
 *           constitutive_parameters)
 *       , H_components(0)
 *       , C_components(Tensor<1, dim>::n_independent_components)
 *       , ad_helper(Tensor<1, dim>::n_independent_components +
 *                   SymmetricTensor<2, dim>::n_independent_components)
 *       , psi(0.0)
 *       , Dpsi(ad_helper.n_independent_variables())
 *       , D2psi(ad_helper.n_independent_variables(),
 *               ad_helper.n_independent_variables())
 *     {}
 * 
 * @endcode
 * 
 * As stated before, due to the way that the automatic differentiation
 * libraries
 * work, the @p ad_helper will always returns the derivatives of the energy
 * density function with respect to all field variables simultaneously.
 * For this reason, it does not make sense to compute the derivatives in
 * the functions `get_B()`, `get_S()`, etc. because we'd be doing a lot of
 * extra computations that are then simply discarded. So, the best way to
 * deal with that is to have a single function call that does all of the
 * calculations up-front, and then we extract the stored data as its needed.
 * That's what we'll do in the `update_internal_data()` method. As the
 * material is rate-independent, we can ignore the DiscreteTime argument.
 * 
 * @code
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     void
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::update_internal_data(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H,
 *       const DiscreteTime &)
 *     {
 *       Assert(determinant(C) > 0, ExcInternalError());
 * 
 * @endcode
 * 
 * Since we reuse the @p ad_helper data structure at each time step,
 * we need to clear it of all stale information before use.
 * 
 * @code
 *       ad_helper.reset();
 * 
 * @endcode
 * 
 * The next step is to set the values for all field components.
 * These define the "point" around which we'll be computing the function
 * gradients and their linearization.
 * The extractors that we created before provide the association between
 * the fields and the registry within the @p ad_helper -- they'll be used
 * repeatedly to ensure that we have the correct interpretation of which
 * variable corresponds to which component of `H` or `C`.
 * 
 * @code
 *       ad_helper.register_independent_variable(H, H_components);
 *       ad_helper.register_independent_variable(C, C_components);
 * 
 * @endcode
 * 
 * Now that we've done the initial setup, we can retrieve the AD
 * counterparts of our fields. These are truly the independent variables
 * for the energy function, and are "sensitive" to the calculations that
 * are performed with them. Notice that the AD number are treated as a
 * special number type, and can be used in many templated classes (in this
 * example, as the scalar type for the Tensor and SymmetricTensor class).
 * 
 * @code
 *       const Tensor<1, dim, ADNumberType> H_ad =
 *         ad_helper.get_sensitive_variables(H_components);
 *       const SymmetricTensor<2, dim, ADNumberType> C_ad =
 *         ad_helper.get_sensitive_variables(C_components);
 * 
 * @endcode
 * 
 * We can also use them in many functions that are templated on the
 * scalar type. So, for these intermediate values that we require,
 * we can perform tensor operations and some mathematical functions.
 * The resulting type will also be an automatically differentiable
 * number, which encodes the operations performed in these functions.
 * 
 * @code
 *       const ADNumberType det_F_ad = std::sqrt(determinant(C_ad));
 *       const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
 *       AssertThrow(det_F_ad > ADNumberType(0.0),
 *                   ExcMessage("Volumetric Jacobian must be positive."));
 * 
 * @endcode
 * 
 * Next we'll compute the scaling function that will cause the shear
 * modulus to change (increase) under the influence of a magnetic field...
 * 
 * @code
 *       const ADNumberType f_mu_e_ad =
 *         1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
 *                 std::tanh((2.0 * H_ad * H_ad) /
 *                           (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()));
 * 
 * @endcode
 * 
 * ... and then we can define the material stored energy density function.
 * We'll see later that this example is sufficiently complex to warrant
 * the use of AD to, at the very least, verify an unassisted
 * implementation.
 * 
 * @code
 *       const ADNumberType psi_ad =
 *         0.5 * this->get_mu_e() * f_mu_e_ad *
 *           (trace(C_ad) - dim - 2.0 * std::log(det_F_ad))                 
 *         + this->get_lambda_e() * std::log(det_F_ad) * std::log(det_F_ad) 
 *         - 0.5 * this->get_mu_0() * this->get_mu_r() * det_F_ad *
 *             (H_ad * C_inv_ad * H_ad); 
 * 
 * @endcode
 * 
 * The stored energy density function is, in fact, the dependent variable
 * for this problem, so as a final step in the  "configuration" phase,
 * we register its definition with the @p ad_helper.
 * 
 * @code
 *       ad_helper.register_dependent_variable(psi_ad);
 * 
 * @endcode
 * 
 * Finally, we can retrieve the resulting value of the stored energy
 * density function, as well as its gradient and Hessian with respect
 * to the input fields, and cache them.
 * 
 * @code
 *       psi = ad_helper.compute_value();
 *       ad_helper.compute_gradient(Dpsi);
 *       ad_helper.compute_hessian(D2psi);
 *     }
 * 
 * @endcode
 * 
 * The following few functions then allow for querying the so-stored value
 * of $\psi_{0}$, and to extract the desired components of the gradient
 * vector and Hessian matrix. We again make use of the extractors to express
 * which parts of the total gradient vector and Hessian matrix we wish to
 * retrieve. They only return the derivatives of the energy function, so
 * for our definitions of the kinetic variables and their linearization a
 * few more manipulations are required to form the desired result.
 * 
 * @code
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     double Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_psi() const
 *     {
 *       return psi;
 *     }
 * 
 * 
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     Tensor<1, dim>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_B() const
 *     {
 *       const Tensor<1, dim> dpsi_dH =
 *         ad_helper.extract_gradient_component(Dpsi, H_components);
 *       return -dpsi_dH;
 *     }
 * 
 * 
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     SymmetricTensor<2, dim>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_S() const
 *     {
 *       const SymmetricTensor<2, dim> dpsi_dC =
 *         ad_helper.extract_gradient_component(Dpsi, C_components);
 *       return 2.0 * dpsi_dC;
 *     }
 * 
 * 
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     SymmetricTensor<2, dim>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_DD() const
 *     {
 *       const Tensor<2, dim> dpsi_dH_dH =
 *         ad_helper.extract_hessian_component(D2psi, H_components, H_components);
 *       return -symmetrize(dpsi_dH_dH);
 *     }
 * 
 * @endcode
 * 
 * Note that for coupled terms the order of the extractor
 * arguments is especially important, as it dictates the order in which
 * the directional derivatives are taken. So, if we'd reversed the order
 * of the extractors in the call to `extract_hessian_component()` then we'd
 * actually have been retrieving part of $\left[ \mathfrak{P}^{\text{tot}}
 * \right]^{T}$.
 * 
 * @code
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     Tensor<3, dim>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_PP() const
 *     {
 *       const Tensor<3, dim> dpsi_dC_dH =
 *         ad_helper.extract_hessian_component(D2psi, C_components, H_components);
 *       return -2.0 * dpsi_dC_dH;
 *     }
 * 
 * 
 *     template <int dim, Differentiation::AD::NumberTypes ADTypeCode>
 *     SymmetricTensor<4, dim>
 *     Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode>::get_HH() const
 *     {
 *       const SymmetricTensor<4, dim> dpsi_dC_dC =
 *         ad_helper.extract_hessian_component(D2psi, C_components, C_components);
 *       return 4.0 * dpsi_dC_dC;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Magnetoviscoelasticconstitutivelawusingsymbolicalgebraanddifferentiation"></a> 
 * <h4>Magneto-viscoelastic constitutive law (using symbolic algebra and differentiation)</h4>
 * 

 * 
 * The second material law that we'll consider will be one that represents
 * a magneto-viscoelastic material with a single dissipative mechanism.
 * We'll consider the free energy density function for such a material to
 * be defined as
 * @f{align*}{
 * \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
 * \right)
 * &= \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * \\ \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * &= \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
 * \right)
 * \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
 * \right]
 * + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
 * - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
 * \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]
 * \\ \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * &= \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)
 * \left[ \mathbf{C}_{v} : \left[
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right] - d - \ln\left(
 * \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
 * @f}
 * with
 * @f[
 * f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
 * @f]
 * @f[
 * f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{v}^{\text{sat}}\right)^{2}} \right),
 * @f]
 * in conjunction with the evolution law for the internal viscous variable
 * @f[
 * \mathbf{C}_{v}^{(t)}
 * = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
 * \mathbf{C}_{v}^{(t-1)}
 * + \frac{\Delta t}{\tau_{v}}
 * \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right]^{-1}
 * \right]
 * @f]
 * that was discretized using a first-order backward difference
 * approximation.
 *     

 * 
 * Again, let us see how this is implemented in a concrete class. Instead of
 * the AD framework used in the previous class, we will now utilize the SD
 * approach. To support this,
 * the class constructor accepts not only the @p constitutive_parameters,
 * but also two additional variables that will be used to initialize
 * a Differentiation::SD::BatchOptimizer. We'll give more context to this
 * later.
 * 
 * @code
 *     template <int dim>
 *     class Magnetoviscoelastic_Constitutive_Law_SD final
 *       : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *     {
 *     public:
 *       Magnetoviscoelastic_Constitutive_Law_SD(
 *         const ConstitutiveParameters &               constitutive_parameters,
 *         const Differentiation::SD::OptimizerType     optimizer_type,
 *         const Differentiation::SD::OptimizationFlags optimization_flags);
 * 
 * @endcode
 * 
 * Like for the automatic differentiation helper, the
 * Differentiation::SD::BatchOptimizer will return a collection of
 * results all at once. So, in order to do that just once, we'll utilize
 * a similar approach to before and do all of the expensive calculations
 * within the `update_internal_data()` function, and cache the results
 * for layer extraction.
 * 
 * @code
 *       virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
 *                                         const Tensor<1, dim> &         H,
 *                                         const DiscreteTime &time) override;
 * 
 *       virtual double get_psi() const override;
 * 
 *       virtual Tensor<1, dim> get_B() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_S() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_DD() const override;
 * 
 *       virtual Tensor<3, dim> get_PP() const override;
 * 
 *       virtual SymmetricTensor<4, dim> get_HH() const override;
 * 
 * @endcode
 * 
 * Since we're dealing with a rate dependent material, we'll have to
 * update the history variable at the appropriate time. That will be the
 * purpose of this function.
 * 
 * @code
 *       virtual void update_end_of_timestep() override;
 * 
 * @endcode
 * 
 * In the `private` part of the class, we will want to
 * keep track of the internal viscous deformation, so the following
 * two (real-valued, non-symbolic) member variables respectively hold
 * - the value of internal variable time step (and, if embedded within a
 * nonlinear solver framework, Newton step), and
 * - the value of internal variable at the previous timestep.
 *       

 * 
 * (We've labeled these variables "Q" so that they're easy to identify;
 * in a sea of calculations it is not necessarily easy to distinguish
 * `Cv` or `C_v` from `C`.)
 * 
 * @code
 *     private:
 *       SymmetricTensor<2, dim> Q_t;
 *       SymmetricTensor<2, dim> Q_t1;
 * 
 * @endcode
 * 
 * As we'll be using symbolic types, we'll need to define some symbolic
 * variables to use with the framework. (They are all suffixed with "SD"
 * to make it easy to distinguish the symbolic types or expressions from
 * real-valued types or scalars.) This can be done once up front
 * (potentially even as `static` variables) to minimize the overhead
 * associated with creating these variables. For the ultimate in generic
 * programming, we can even describe the constitutive parameters
 * symbolically, *potentially* allowing a single class instance to be
 * reused with different inputs for these values too.
 *       

 * 
 * These are the symbolic scalars that represent the elastic, viscous, and
 * magnetic material parameters
 * (defined mostly in the same order as they appear in the @p ConstitutiveParameters
 * class). We also store a symbolic expression, @p delta_t_sd, that represents the
 * time step size):
 * 
 * @code
 *       const Differentiation::SD::Expression mu_e_sd;
 *       const Differentiation::SD::Expression mu_e_inf_sd;
 *       const Differentiation::SD::Expression mu_e_h_sat_sd;
 *       const Differentiation::SD::Expression lambda_e_sd;
 *       const Differentiation::SD::Expression mu_v_sd;
 *       const Differentiation::SD::Expression mu_v_inf_sd;
 *       const Differentiation::SD::Expression mu_v_h_sat_sd;
 *       const Differentiation::SD::Expression tau_v_sd;
 *       const Differentiation::SD::Expression delta_t_sd;
 *       const Differentiation::SD::Expression mu_r_sd;
 * 
 * @endcode
 * 
 * Next we define some tensorial symbolic variables that represent the
 * independent field variables, upon which the energy density function
 * is parameterized:
 * 
 * @code
 *       const Tensor<1, dim, Differentiation::SD::Expression>          H_sd;
 *       const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_sd;
 * 
 * @endcode
 * 
 * And similarly we have the symbolic representation of the internal
 * viscous variables (both its current value and its value at the
 * previous timestep):
 * 
 * @code
 *       const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t_sd;
 *       const SymmetricTensor<2, dim, Differentiation::SD::Expression> Q_t1_sd;
 * 
 * @endcode
 * 
 * We should also store the definitions of the dependent expressions:
 * Although we'll only compute them once, we require them to retrieve
 * data from the @p optimizer that is declared below.
 * Furthermore, when serializing a material class like this one (not done
 * as a part of this tutorial) we'd either need to serialize these
 * expressions as well or we'd need to reconstruct them upon reloading.
 * 
 * @code
 *       Differentiation::SD::Expression                          psi_sd;
 *       Tensor<1, dim, Differentiation::SD::Expression>          B_sd;
 *       SymmetricTensor<2, dim, Differentiation::SD::Expression> S_sd;
 *       SymmetricTensor<2, dim, Differentiation::SD::Expression> BB_sd;
 *       Tensor<3, dim, Differentiation::SD::Expression>          PP_sd;
 *       SymmetricTensor<4, dim, Differentiation::SD::Expression> HH_sd;
 * 
 * @endcode
 * 
 * The next variable is then the optimizer that is used to evaluate the
 * dependent functions. More specifically, it provides the possibility to
 * accelerate the evaluation of the symbolic dependent expressions. This
 * is a vital tool, because the native evaluation of lengthy expressions
 * (using no method of acceleration, but rather direct evaluation directly
 * of the symbolic expressions) can be very slow. The
 * Differentiation::SD::BatchOptimizer class provides a mechanism by which
 * to transform the symbolic expression tree into another code path that,
 * for example, shares intermediate results between the various dependent
 * expressions (meaning that these intermediate values only get calculated
 * once per evaluation) and/or compiling the code using a just-in-time
 * compiler (thereby retrieving near-native performance for the evaluation
 * step).
 *       

 * 
 * Performing this code transformation is very computationally expensive,
 * so we store the optimizer so that it is done just once per class
 * instance. This also further motivates the decision to make the
 * constitutive parameters themselves symbolic. We could then reuse a
 * single instance
 * of this @p optimizer across several materials (with the same energy
 * function, of course) and potentially multiple continuum points (if
 * embedded within a finite element simulation).
 *       

 * 
 * As specified by the template parameter, the numerical result will be of
 * type <tt>double</tt>.
 * 
 * @code
 *       Differentiation::SD::BatchOptimizer<double> optimizer;
 * 
 * @endcode
 * 
 * During the evaluation phase, we must map the symbolic variables to
 * their real-valued counterparts. The next method will provide this
 * functionality.
 *       

 * 
 * The final method of this class will configure the @p optimizer.
 * 
 * @code
 *       Differentiation::SD::types::substitution_map
 *       make_substitution_map(const SymmetricTensor<2, dim> &C,
 *                             const Tensor<1, dim> &         H,
 *                             const double                   delta_t) const;
 * 
 *       void initialize_optimizer();
 *     };
 * 
 * @endcode
 * 
 * As the resting deformation state is one at which the material is
 * considered to be completely relaxed, the internal viscous variables are
 * initialized with the identity tensor, i.e. $\mathbf{C}_{v} = \mathbf{I}$.
 * The various symbolic variables representing the constitutive parameters,
 * time step size, and field and internal variables all get a unique
 * identifier. The optimizer is passed the two parameters that declare which
 * optimization (acceleration) technique should be applied, as well as
 * which additional steps should be taken by the CAS to help improve
 * performance during evaluation.
 * 
 * @code
 *     template <int dim>
 *     Magnetoviscoelastic_Constitutive_Law_SD<dim>::
 *       Magnetoviscoelastic_Constitutive_Law_SD(
 *         const ConstitutiveParameters &               constitutive_parameters,
 *         const Differentiation::SD::OptimizerType     optimizer_type,
 *         const Differentiation::SD::OptimizationFlags optimization_flags)
 *       : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
 *           constitutive_parameters)
 *       , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
 *       , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
 *       , mu_e_sd("mu_e")
 *       , mu_e_inf_sd("mu_e_inf")
 *       , mu_e_h_sat_sd("mu_e_h_sat")
 *       , lambda_e_sd("lambda_e")
 *       , mu_v_sd("mu_v")
 *       , mu_v_inf_sd("mu_v_inf")
 *       , mu_v_h_sat_sd("mu_v_h_sat")
 *       , tau_v_sd("tau_v")
 *       , delta_t_sd("delta_t")
 *       , mu_r_sd("mu_r")
 *       , H_sd(Differentiation::SD::make_vector_of_symbols<dim>("H"))
 *       , C_sd(Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("C"))
 *       , Q_t_sd(
 *           Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t"))
 *       , Q_t1_sd(
 *           Differentiation::SD::make_symmetric_tensor_of_symbols<2, dim>("Q_t1"))
 *       , optimizer(optimizer_type, optimization_flags)
 *     {
 *       initialize_optimizer();
 *     }
 * 
 * @endcode
 * 
 * The substitution map simply pairs all of the following data together:
 * - the constitutive parameters (with values retrieved from the base
 * class),
 * - the time step size (with its value retrieved from the time
 * discretizer),
 * - the field values (with their values being prescribed by an external
 * function that is calling into this @p Magnetoviscoelastic_Constitutive_Law_SD instance), and
 * - the current and previous internal viscous deformation (with their
 * values
 * stored within this class instance).
 * 
 * @code
 *     template <int dim>
 *     Differentiation::SD::types::substitution_map
 *     Magnetoviscoelastic_Constitutive_Law_SD<dim>::make_substitution_map(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H,
 *       const double                   delta_t) const
 *     {
 *       return Differentiation::SD::make_substitution_map(
 *         std::make_pair(mu_e_sd, this->get_mu_e()),
 *         std::make_pair(mu_e_inf_sd, this->get_mu_e_inf()),
 *         std::make_pair(mu_e_h_sat_sd, this->get_mu_e_h_sat()),
 *         std::make_pair(lambda_e_sd, this->get_lambda_e()),
 *         std::make_pair(mu_v_sd, this->get_mu_v()),
 *         std::make_pair(mu_v_inf_sd, this->get_mu_v_inf()),
 *         std::make_pair(mu_v_h_sat_sd, this->get_mu_v_h_sat()),
 *         std::make_pair(tau_v_sd, this->get_tau_v()),
 *         std::make_pair(delta_t_sd, delta_t),
 *         std::make_pair(mu_r_sd, this->get_mu_r()),
 *         std::make_pair(H_sd, H),
 *         std::make_pair(C_sd, C),
 *         std::make_pair(Q_t_sd, Q_t),
 *         std::make_pair(Q_t1_sd, Q_t1));
 *     }
 * 
 * @endcode
 * 
 * Due to the "natural" use of the symbolic expressions, much of the
 * procedure to configure the @p optimizer looks very similar to that which
 * is used to construct the automatic differentiation helper.
 * Nevertheless, we'll detail these steps again to highlight the differences
 * that underlie the two frameworks.
 *     

 * 
 * The function starts with expressions that symbolically encode the
 * determinant of the deformation gradient (as expressed in terms of the
 * right Cauchy-Green deformation tensor, our primary field variable), as
 * well as the inverse of $\mathbf{C}$ itself:
 * 
 * @code
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law_SD<dim>::initialize_optimizer()
 *     {
 *       const Differentiation::SD::Expression det_F_sd =
 *         std::sqrt(determinant(C_sd));
 *       const SymmetricTensor<2, dim, Differentiation::SD::Expression> C_inv_sd =
 *         invert(C_sd);
 * 
 * @endcode
 * 
 * Next is the symbolic representation of the saturation function for
 * the elastic part of the free energy density function, followed by the
 * magnetoelastic contribution to the free energy density function.
 * This all has the same structure as we'd seen previously.
 * 
 * @code
 *       const Differentiation::SD::Expression f_mu_e_sd =
 *         1.0 +
 *         (mu_e_inf_sd / mu_e_sd - 1.0) *
 *           std::tanh((2.0 * H_sd * H_sd) / (mu_e_h_sat_sd * mu_e_h_sat_sd));
 * 
 *       const Differentiation::SD::Expression psi_ME_sd =
 *         0.5 * mu_e_sd * f_mu_e_sd *
 *           (trace(C_sd) - dim - 2.0 * std::log(det_F_sd)) +
 *         lambda_e_sd * std::log(det_F_sd) * std::log(det_F_sd) -
 *         0.5 * this->get_mu_0() * mu_r_sd * det_F_sd * (H_sd * C_inv_sd * H_sd);
 * 
 * @endcode
 * 
 * In addition, we define the magneto-viscoelastic contribution to the
 * free energy density function. The first component required to implement
 * this is a scaling function that will cause the viscous shear modulus to
 * change (increase) under the influence of a magnetic field (see
 * @cite Pelteret2018a, equation 29). Thereafter we can compute the
 * dissipative component of the energy density function; its expression
 * is stated in @cite Pelteret2018a (equation 28), which is a
 * straight-forward extension of an energy density function formulated in
 * @cite Linder2011a (equation 46).
 * 
 * @code
 *       const Differentiation::SD::Expression f_mu_v_sd =
 *         1.0 +
 *         (mu_v_inf_sd / mu_v_sd - 1.0) *
 *           std::tanh((2.0 * H_sd * H_sd) / (mu_v_h_sat_sd * mu_v_h_sat_sd));
 * 
 *       const Differentiation::SD::Expression psi_MVE_sd =
 *         0.5 * mu_v_sd * f_mu_v_sd *
 *         (Q_t_sd * (std::pow(det_F_sd, -2.0 / dim) * C_sd) - dim -
 *          std::log(determinant(Q_t_sd)));
 * 
 * @endcode
 * 
 * From these building blocks, we can then define the material's total
 * free energy density function:
 * 
 * @code
 *       psi_sd = psi_ME_sd + psi_MVE_sd;
 * 
 * @endcode
 * 
 * As it stands, to the CAS the variable @p Q_t_sd appears
 * to be independent of @p C_sd. Our tensorial symbolic expression
 * @p Q_t_sd just has an identifier associated with it, and there is
 * nothing that links it to the other tensorial symbolic expression
 * @p C_sd. So any derivatives taken with respect to @p C_sd will ignore
 * this inherent dependence which, as we can see from the evolution law,
 * is in fact
 * $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t \right)$.
 * This means that deriving any function $f = f(\mathbf{C}, \mathbf{Q})$
 * with respect to  $\mathbf{C}$ will return partial derivatives
 * $\frac{\partial f(\mathbf{C}, \mathbf{Q})}{\partial \mathbf{C}}
 * \Big\vert_{\mathbf{Q}}$ as opposed to the total derivative
 * $\frac{d f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))}{d \mathbf{C}} =
 * \frac{\partial f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))}{\partial
 * \mathbf{C}} \Big\vert_{\mathbf{Q}} + \frac{\partial f(\mathbf{C},
 * \mathbf{Q}(\mathbf{C}))}{\partial \mathbf{Q}}
 * \Big\vert_{\mathbf{C}} : \frac{d \mathbf{Q}(\mathbf{C}))}{d
 * \mathbf{C}}$.
 *       

 * 
 * By contrast, with the current AD libraries the total derivative would
 * always be returned. This implies that the computed kinetic variables
 * would be incorrect for this class of material model, making AD the
 * incorrect tool from which to derive (at the continuum point level) the
 * constitutive law for this dissipative material from an energy density
 * function.
 *       

 * 
 * It is this specific level of control that characterizes a defining
 * difference difference between the SD and AD frameworks. In a few lines
 * we'll be manipulating the expression for the internal variable
 * @p Q_t_sd such that it produces the correct linearization.
 * 

 * 
 * But, first, we'll compute the symbolic expressions for the kinetic
 * variables, i.e., the magnetic induction vector and the Piola-Kirchhoff
 * stress tensor. The code that performs the differentiation quite closely
 * mimics the definition stated in the theory.
 * 
 * @code
 *       B_sd = -Differentiation::SD::differentiate(psi_sd, H_sd);
 *       S_sd = 2.0 * Differentiation::SD::differentiate(psi_sd, C_sd);
 * 
 * @endcode
 * 
 * Since the next step is to linearize the above, it is the appropriate
 * time to inform the CAS of the explicit dependency of @p Q_t_sd on @p C_sd,
 * i.e., state that $\mathbf{C}_{v} = \mathbf{C}_{v} \left( \mathbf{C}, t
 * \right)$. This means that all future differential operations made with
 * respect
 * to @p C_sd will take into account this dependence (i.e., compute total derivatives).
 * In other words, we will transform some expression such that their
 * intrinsic parameterization changes from $f(\mathbf{C}, \mathbf{Q})$
 * to $f(\mathbf{C}, \mathbf{Q}(\mathbf{C}))$.
 *       

 * 
 * To do this, we consider the time-discrete evolution law.
 * From that, we have the explicit expression for the internal
 * variable in terms of its history as well as the primary
 * field variable. That is what it described in this expression:
 * 
 * @code
 *       const SymmetricTensor<2, dim, Differentiation::SD::Expression>
 *         Q_t_sd_explicit =
 *           (1.0 / (1.0 + delta_t_sd / tau_v_sd)) *
 *           (Q_t1_sd +
 *            (delta_t_sd / tau_v_sd * std::pow(det_F_sd, 2.0 / dim) * C_inv_sd));
 * 
 * @endcode
 * 
 * Next we produce an intermediate substitution map, which will take
 * every instance of @p Q_t_sd (our identifier) found in an expression
 * and replace it with the full expression held in @p Q_t_sd_explicit.
 * 
 * @code
 *       const Differentiation::SD::types::substitution_map
 *         substitution_map_explicit = Differentiation::SD::make_substitution_map(
 *           std::make_pair(Q_t_sd, Q_t_sd_explicit));
 * 
 * @endcode
 * 
 * We can the perform this substitution on the two kinetic variables
 * and immediately differentiate the result that appears after that
 * substitution with the field variables. (If you'd like, this could
 * be split up into two steps with the intermediate results stored in
 * a temporary variable.) Again, if you overlook the "complexity"
 * generated by the substitution, these calls that linearize the kinetic
 * variables and produce the three tangent tensors quite closely resembles
 * what's stated in the theory.
 * 
 * @code
 *       BB_sd = symmetrize(Differentiation::SD::differentiate(
 *         Differentiation::SD::substitute(B_sd, substitution_map_explicit),
 *         H_sd));
 *       PP_sd = -Differentiation::SD::differentiate(
 *         Differentiation::SD::substitute(S_sd, substitution_map_explicit), H_sd);
 *       HH_sd =
 *         2.0 *
 *         Differentiation::SD::differentiate(
 *           Differentiation::SD::substitute(S_sd, substitution_map_explicit),
 *           C_sd);
 * 
 * @endcode
 * 
 * Now we need to tell the @p optimizer what entries we need to provide
 * numerical values for in order for it to successfully perform its
 * calculations. These essentially act as the input arguments to
 * all dependent functions that the @p optimizer must evaluate.
 * They are, collectively, the independent variables
 * for the problem, the history variables, the time step sie and the
 * constitutive parameters (since we've not hard encoded them in the
 * energy density function).
 *       

 * 
 * So what we really want is to provide it a collection of
 * symbols, which one could accomplish in this way:
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * optimizer.register_symbols(Differentiation::SD::make_symbol_map(
 *   mu_e_sd, mu_e_inf_sd, mu_e_h_sat_sd, lambda_e_sd,
 *   mu_v_sd, mu_v_inf_sd, mu_v_h_sat_sd, tau_v_sd,
 *   delta_t_sd, mu_r_sd,
 *   H_sd, C_sd,
 *   Q_t_sd, Q_t1_sd));
 * @endcode
 * </div>
 * But this is all actually already encoded as the keys of the
 * substitution map. Doing the above would also mean that we
 * need to manage the symbols in two places (here and when constructing
 * the substitution map), which is annoying and a potential source of
 * error if this material class is modified or extended.
 * Since we're not interested in the values at this point,
 * it is alright if the substitution map is filled with invalid data
 * for the values associated with each key entry.
 * So we'll simply create a fake substitution map, and extract the
 * symbols from that. Note that any substitution map passed to the
 * @p optimizer will have to, at the very least, contain entries for
 * these symbols.
 * 
 * @code
 *       optimizer.register_symbols(
 *         Differentiation::SD::Utilities::extract_symbols(
 *           make_substitution_map({}, {}, 0)));
 * 
 * @endcode
 * 
 * We then inform the optimizer of what values we want calculated, which
 * in our situation encompasses all of the dependent variables (namely
 * the energy density function and its various derivatives).
 * 
 * @code
 *       optimizer.register_functions(psi_sd, B_sd, S_sd, BB_sd, PP_sd, HH_sd);
 * 
 * @endcode
 * 
 * The last step is to finalize the optimizer. With this call it will
 * determine an equivalent code path that will evaluate all of the
 * dependent functions at once, but with less computational
 * cost than when evaluating the symbolic expression directly.
 * Note: This is an expensive call, so we want execute it as few times
 * as possible. We've done it in the constructor of our class, which
 * achieves the goal of being called only once per class instance.
 * 
 * @code
 *       optimizer.optimize();
 *     }
 * 
 * @endcode
 * 
 * Since the configuration of the @p optimizer was done up front, there's
 * very little to do each time we want to compute kinetic variables or
 * their linearization (derivatives).
 * 
 * @code
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_internal_data(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H,
 *       const DiscreteTime &           time)
 *     {
 * @endcode
 * 
 * To update the internal history variable, we first need to compute
 * a few fundamental quantities, which we've seen before.
 * We can also ask the time discretizer for the time step size that
 * was used to iterate from the previous time step to the current one.
 * 
 * @code
 *       const double delta_t = this->get_delta_t(time);
 * 
 *       const double                  det_F = std::sqrt(determinant(C));
 *       const SymmetricTensor<2, dim> C_inv = invert(C);
 *       AssertThrow(det_F > 0.0,
 *                   ExcMessage("Volumetric Jacobian must be positive."));
 * 
 * @endcode
 * 
 * Now we can update the (real valued) internal viscous deformation
 * tensor, as per the definition given by the evolution law in conjunction
 * with the chosen time discretization scheme.
 * 
 * @code
 *       Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
 *             (Q_t1 + (delta_t / this->get_tau_v()) * std::pow(det_F, 2.0 / dim) *
 *                       C_inv);
 * 
 * @endcode
 * 
 * Next we pass the optimizer the numeric values that we wish the
 * independent variables, time step size and (implicit to this call),
 * the constitutive parameters to represent.
 * 
 * @code
 *       const auto substitution_map = make_substitution_map(C, H, delta_t);
 * 
 * @endcode
 * 
 * When making this next call, the call path used to (numerically)
 * evaluate the dependent functions is quicker than dictionary
 * substitution.
 * 
 * @code
 *       optimizer.substitute(substitution_map);
 *     }
 * 
 * @endcode
 * 
 * Having called `update_internal_data()`, it is then valid to
 * extract data from the optimizer.
 * When doing the evaluation, we need the exact symbolic expressions of
 * the data to extracted from the optimizer. The implication of this
 * is that we needed to store the symbolic expressions of all dependent
 * variables for the lifetime of the optimizer (naturally, the same
 * is implied for the input variables).
 * 
 * @code
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_psi() const
 *     {
 *       return optimizer.evaluate(psi_sd);
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_B() const
 *     {
 *       return optimizer.evaluate(B_sd);
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim>
 *     Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_S() const
 *     {
 *       return optimizer.evaluate(S_sd);
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim>
 *     Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_DD() const
 *     {
 *       return optimizer.evaluate(BB_sd);
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_PP() const
 *     {
 *       return optimizer.evaluate(PP_sd);
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<4, dim>
 *     Magnetoviscoelastic_Constitutive_Law_SD<dim>::get_HH() const
 *     {
 *       return optimizer.evaluate(HH_sd);
 *     }
 * 
 * @endcode
 * 
 * When moving forward in time, the "current" state of the internal variable
 * instantaneously defines the state at the "previous" timestep. As such, we
 * record value of history variable for use as the "past value" at the next
 * time step.
 * 
 * @code
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law_SD<dim>::update_end_of_timestep()
 *     {
 *       Q_t1 = Q_t;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="AmorecomplexexamplecontinuedParametersandhandderivedmaterialclasses"></a> 
 * <h3>A more complex example (continued): Parameters and hand-derived material classes</h3>
 * 

 * 
 * Now that we've seen how the AD and SD frameworks can make light(er) work
 * of defining these constitutive laws, we'll implement the equivalent
 * classes by hand for the purpose of verification and to do some
 * preliminary benchmarking of the frameworks versus a native
 * implementation.
 *     

 * 
 * At the expense of the author's sanity, what is documented below
 * (hopefully accurately) are the full definitions for the kinetic variables
 * and their tangents, as well as some intermediate computations. Since the
 * structure and design of the constitutive law classes has been outlined
 * earlier, we'll gloss over it and simply delineate between the various
 * stages of calculations in the `update_internal_data()` method definition.
 * It should be easy enough to link the derivative calculations (with their
 * moderately expressive variable names) to their documented definitions
 * that appear in the class descriptions.
 * We will, however, take the opportunity to present two different paradigms
 * for implementing constitutive law classes. The second will provide more
 * flexibility than the first (thereby making it more easily extensible,
 * in the author's opinion) at the expense of some performance.
 * 

 * 
 * 
 * <a name="Magnetoelasticconstitutivelawhandderived"></a> 
 * <h4>Magnetoelastic constitutive law (hand-derived)</h4>
 * 

 * 
 * From the stored energy that, as mentioned earlier, is defined as
 * @f[
 * \psi_{0} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * = \frac{1}{2} \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
 * \right]
 * + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
 * - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
 * \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]
 * @f]
 * with
 * @f[
 * f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right) ,
 * \\ \text{det}(\mathbf{F}) = \sqrt{\text{det}(\mathbf{C})}
 * @f]
 * for this magnetoelastic material, the first derivatives that correspond
 * to the magnetic induction vector and total Piola-Kirchhoff stress
 * tensor are
 * @f[
 * \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right)
 * \dealcoloneq - \frac{d \psi_{0}}{d \boldsymbol{\mathbb{H}}}
 * = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
 * (\text{det}(\mathbf{F}))
 * \right] \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{d \boldsymbol{\mathbb{H}}}
 * + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \left[ \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}
 * \right]
 * @f]
 * @f{align}
 * \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right)
 * \dealcoloneq 2 \frac{d \psi_{0} \left( \mathbf{C},
 * \boldsymbol{\mathbb{H}} \right)}{d \mathbf{C}}
 * &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
 * - 2 \frac{1}{\text{det}(\mathbf{F})}
 * \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
 * + 4 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right)
 * \frac{1}{\text{det}(\mathbf{F})} \frac{d\,\text{det}(\mathbf{F})}{d
 * \mathbf{C}}
 * - \mu_{0} \mu_{r} \left[
 * \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right] \frac{d\,\text{det}(\mathbf{F})}{d
 * \mathbf{C}} + \text{det}(\mathbf{F}) \frac{d \left[
 * \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}
 * \right]}{d \mathbf{C}} \right]
 * \\ &= \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ \mathbf{I} - \mathbf{C}^{-1} \right]
 * + 2 \lambda_{e} \ln \left(\text{det}(\mathbf{F}) \right) \mathbf{C}^{-1}
 * - \mu_{0} \mu_{r} \left[
 * \frac{1}{2}  \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
 * \cdot \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})
 * \mathbf{C}^{-1}
 * - \text{det}(\mathbf{F})
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
 * \right]
 * @f}
 * with
 * @f[
 * \frac{d f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
 * \boldsymbol{\mathbb{H}}}
 * = \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
 * \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}}
 * \boldsymbol{\mathbb{H}} \right]
 * @f]
 * @f[
 * \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
 * = \mathbf{I}
 * \quad \text{(the second-order identity tensor)}
 * @f]
 * @f[
 * \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}}
 * = \frac{1}{2} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
 * @f]
 * @f[
 * \frac{d C^{-1}_{ab}}{d C_{cd}}
 * = - \text{sym} \left( C^{-1}_{ac} C^{-1}_{bd} \right)
 * = -\frac{1}{2} \left[ C^{-1}_{ac} C^{-1}_{bd} + C^{-1}_{ad} C^{-1}_{bc}
 * \right]
 * @f]
 * @f[
 * \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
 * = - \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
 * @f]
 * The use of the symmetry operator $\text{sym} \left( \bullet \right)$ in
 * the one derivation above helps to ensure that the resulting rank-4
 * tensor, which holds minor symmetries due to the symmetry of $\mathbf{C}$,
 * still maps rank-2 symmetric tensors to rank-2 symmetric tensors. See the
 * SymmetricTensor class documentation and the introduction to step-44 and
 * for further explanation as to what symmetry means in the context of
 * fourth-order tensors.
 *     

 * 
 * The linearization of each of the kinematic variables with respect to
 * their arguments are
 * @f[
 * \mathbb{D} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
 * = - \frac{1}{2} \mu_{e} \left[ \text{tr}(\mathbf{C}) - d - 2 \ln
 * (\text{det}(\mathbf{F}))
 * \right] \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{d \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}}
 * + \mu_{0} \mu_{r} \text{det}(\mathbf{F}) \mathbf{C}^{-1}
 * @f]
 * @f{align}
 * \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d \boldsymbol{\mathbb{H}}}
 * &= - \mu_{e}
 * \left[ \frac{d\,\text{tr}(\mathbf{C})}{d \mathbf{C}}
 * - 2 \frac{1}{\text{det}(\mathbf{F})}
 * \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \right]
 * \otimes \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}}
 * \right)}}{d \boldsymbol{\mathbb{H}}}
 * + \mu_{0} \mu_{r} \left[
 * \frac{d\,\text{det}(\mathbf{F})}{d \mathbf{C}} \otimes
 * \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}
 * \right]}{d \boldsymbol{\mathbb{H}}} \right]
 * + \text{det}(\mathbf{F})
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
 * \cdot \boldsymbol{\mathbb{H}}
 * \right]}{d \mathbf{C} \otimes d \boldsymbol{\mathbb{H}}}
 * \\ &= - \mu_{e}
 * \left[ \mathbf{I} - \mathbf{C}^{-1} \right] \otimes
 * \frac{d f_{\mu_{e} \left( \boldsymbol{\mathbb{H}} \right)}}{d
 * \boldsymbol{\mathbb{H}}}
 * + \mu_{0} \mu_{r} \left[
 * \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
 * \right]
 * + \text{det}(\mathbf{F})
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
 * \cdot \boldsymbol{\mathbb{H}}
 * \right]}{d \mathbf{C} \otimes \mathbf{C} \boldsymbol{\mathbb{H}}}
 * @f}
 * @f{align}
 * \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d \mathbf{C}}
 * &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
 * + 4 \lambda_{e} \left[ \mathbf{C}^{-1} \otimes \left[
 * \frac{1}{\text{det}(\mathbf{F})} \frac{d \, \text{det}(\mathbf{F})}{d
 * \mathbf{C}} \right] + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}} \right]
 * \\ &- \mu_{0} \mu_{r}  \left[
 * \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \frac{d \left[
 * \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]}{d \mathbf{C}}
 * + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right] \mathbf{C}^{-1} \otimes \frac{d \,
 * \text{det}(\mathbf{F})}{d \mathbf{C}}
 * + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}}
 * \right]
 * \\ &+ 2 \mu_{0} \mu_{r} \left[ \left[
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
 * \right] \otimes \frac{d \, \text{det}(\mathbf{F})}{d \mathbf{C}}
 * - \text{det}(\mathbf{F})
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
 * \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
 * \mathbf{C}}
 * \right]
 * \\ &= 2 \mu_{e} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ - \frac{d \mathbf{C}^{-1}}{d \mathbf{C}} \right]
 * + 4 \lambda_{e} \left[ \frac{1}{2} \mathbf{C}^{-1} \otimes
 * \mathbf{C}^{-1} + \ln \left(\text{det}(\mathbf{F}) \right) \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}} \right]
 * \\ &- \mu_{0} \mu_{r}  \left[
 * - \text{det}(\mathbf{F}) \mathbf{C}^{-1} \otimes \left[ \left[
 * \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \right]
 * + \frac{1}{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F})  \mathbf{C}^{-1}
 * \otimes \mathbf{C}^{-1}
 * + \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right] \text{det}(\mathbf{F}) \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}}
 * \right]
 * \\ &+ 2 \mu_{0} \mu_{r} \left[ \frac{1}{2} \text{det}(\mathbf{F}) \left[
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}} \right]
 * \right] \otimes \mathbf{C}^{-1}
 * - \text{det}(\mathbf{F})
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1}
 * \cdot \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
 * \mathbf{C}}
 * \right]
 * @f}
 * with
 * @f[
 * \frac{d^{2} f_{\mu_{e}} \left( \boldsymbol{\mathbb{H}} \right)}{d
 * \boldsymbol{\mathbb{H}} \otimes d \boldsymbol{\mathbb{H}}}
 * = -2 \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
 * \text{sech}^{2} \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
 * \left[ \frac{4} {\left(h_{e}^{\text{sat}}\right)^{2}} \mathbf{I}
 * \right]
 * @f]
 * @f[
 * \frac{d \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}
 * \right]}{d \boldsymbol{\mathbb{H}}}
 * = 2 \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
 * @f]
 * @f[
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d
 * \boldsymbol{\mathbb{H}}} \Rightarrow \frac{d^{2} \left[ \mathbb{H}_{e}
 * C^{-1}_{ef} \mathbb{H}_{f}
 * \right]}{d C_{ab} d \mathbb{H}_{c}}
 * = - C^{-1}_{ac} C^{-1}_{be} \mathbb{H}_{e} - C^{-1}_{ae} \mathbb{H}_{e}
 * C^{-1}_{bc}
 * @f]
 * @f{align}
 * \frac{d^{2} \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}}\right]}{d \mathbf{C} \otimes d \mathbf{C}}
 * &= -\frac{d \left[\left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
 * \right] \otimes
 * \left[ \mathbf{C}^{-1} \cdot \boldsymbol{\mathbb{H}}
 * \right]\right]}{d \mathbf{C}}
 * \\ \Rightarrow
 * \frac{d^{2} \left[ \mathbb{H}_{e} C^{-1}_{ef} \mathbb{H}_{f}
 * \right]}{d C_{ab} d C_{cd}}
 * &= \text{sym} \left( C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf}
 * \mathbb{H}_{f} C^{-1}_{bd}
 * + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f}
 * C^{-1}_{ad} \right)
 * \\ &= \frac{1}{2} \left[
 * C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{cf} \mathbb{H}_{f} C^{-1}_{bd}
 * + C^{-1}_{ae} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{bc}
 * + C^{-1}_{ce} \mathbb{H}_{e} C^{-1}_{bf} \mathbb{H}_{f} C^{-1}_{ad}
 * + C^{-1}_{be} \mathbb{H}_{e} C^{-1}_{df} \mathbb{H}_{f} C^{-1}_{ac}
 * \right]
 * @f}
 *     

 * 
 * Well, that escalated quickly -- although the the definition of $\psi_{0}$
 * and $f_{\mu_e}$ might have given some hints that the calculating
 * the kinetic fields and their linearization would take some effort, it is
 * likely that there's a little more complexity to the final definitions
 * that perhaps initially thought.
 * Knowing what we now do, it's probably fair to say that we really do not
 * want to compute first and second derivatives of these functions with
 * respect to their arguments -- regardless of well we did in calculus
 * classes, or how good a programmer we may be.
 *     

 * 
 * In the class method definition where these are ultimately implemented,
 * we've composed these calculations slightly differently. Some intermediate
 * steps are also retained to give another perspective of how to
 * systematically compute the derivatives. Additionally, some calculations
 * are decomposed less or further to reuse some of the intermediate values
 * and, hopefully, aid the reader to follow the derivative operations.
 * 
 * @code
 *     template <int dim>
 *     class Magnetoelastic_Constitutive_Law final
 *       : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *     {
 *     public:
 *       Magnetoelastic_Constitutive_Law(
 *         const ConstitutiveParameters &constitutive_parameters);
 * 
 *       virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
 *                                         const Tensor<1, dim> &         H,
 *                                         const DiscreteTime &) override;
 * 
 *       virtual double get_psi() const override;
 * 
 *       virtual Tensor<1, dim> get_B() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_S() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_DD() const override;
 * 
 *       virtual Tensor<3, dim> get_PP() const override;
 * 
 *       virtual SymmetricTensor<4, dim> get_HH() const override;
 * 
 *     private:
 *       double                  psi;
 *       Tensor<1, dim>          B;
 *       SymmetricTensor<2, dim> S;
 *       SymmetricTensor<2, dim> BB;
 *       Tensor<3, dim>          PP;
 *       SymmetricTensor<4, dim> HH;
 *     };
 * 
 * 
 *     template <int dim>
 *     Magnetoelastic_Constitutive_Law<dim>::Magnetoelastic_Constitutive_Law(
 *       const ConstitutiveParameters &constitutive_parameters)
 *       : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
 *           constitutive_parameters)
 *       , psi(0.0)
 *     {}
 * 
 * @endcode
 * 
 * For this class's update method, we'll simply precompute a collection of
 * intermediate values (for function evaluations, derivative calculations,
 * and the like) and "manually" arrange them in the order that's required
 * to maximize their reuse. This means that we have to manage this
 * ourselves, and decide what values must be compute before others, all
 * while keeping some semblance of order or structure in the code itself.
 * It's effective, but perhaps a little tedious. It also doesn't do too much
 * to help future extension of the class, because all of these values remain
 * local to this single method.
 *     

 * 
 * Interestingly, this basic technique of precomputing intermediate
 * expressions that are used in more than one place has a name:
 * [common subexpression elimination
 * (CSE)](https://en.wikipedia.org/wiki/Common_subexpression_elimination).
 * It is a strategy used by Computer Algebra Systems to reduce the
 * computational expense when they are tasked with evaluating similar
 * expressions.
 * 
 * @code
 *     template <int dim>
 *     void Magnetoelastic_Constitutive_Law<dim>::update_internal_data(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H,
 *       const DiscreteTime &)
 *     {
 *       const double                  det_F = std::sqrt(determinant(C));
 *       const SymmetricTensor<2, dim> C_inv = invert(C);
 *       AssertThrow(det_F > 0.0,
 *                   ExcMessage("Volumetric Jacobian must be positive."));
 * 
 * @endcode
 * 
 * The saturation function for the magneto-elastic energy.
 * 
 * @code
 *       const double two_h_dot_h_div_h_sat_squ =
 *         (2.0 * H * H) / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat());
 *       const double tanh_two_h_dot_h_div_h_sat_squ =
 *         std::tanh(two_h_dot_h_div_h_sat_squ);
 * 
 *       const double f_mu_e =
 *         1.0 + (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
 *                 tanh_two_h_dot_h_div_h_sat_squ;
 * 
 * @endcode
 * 
 * The first derivative of the saturation function, noting that
 * $\frac{d \tanh(x)}{dx} = \text{sech}^{2}(x)$.
 * 
 * @code
 *       const double dtanh_two_h_dot_h_div_h_sat_squ =
 *         std::pow(1.0 / std::cosh(two_h_dot_h_div_h_sat_squ), 2.0);
 *       const Tensor<1, dim> dtwo_h_dot_h_div_h_sat_squ_dH =
 *         2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) * H;
 * 
 *       const Tensor<1, dim> df_mu_e_dH =
 *         (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
 *         (dtanh_two_h_dot_h_div_h_sat_squ * dtwo_h_dot_h_div_h_sat_squ_dH);
 * 
 * @endcode
 * 
 * The second derivative of saturation function, noting that
 * $\frac{d \text{sech}^{2}(x)}{dx} = -2 \tanh(x) \text{sech}^{2}(x)$.
 * 
 * @code
 *       const double d2tanh_two_h_dot_h_div_h_sat_squ =
 *         -2.0 * tanh_two_h_dot_h_div_h_sat_squ * dtanh_two_h_dot_h_div_h_sat_squ;
 *       const SymmetricTensor<2, dim> d2two_h_dot_h_div_h_sat_squ_dH_dH =
 *         2.0 * 2.0 / (this->get_mu_e_h_sat() * this->get_mu_e_h_sat()) *
 *         Physics::Elasticity::StandardTensors<dim>::I;
 * 
 *       const SymmetricTensor<2, dim> d2f_mu_e_dH_dH =
 *         (this->get_mu_e_inf() / this->get_mu_e() - 1.0) *
 *         (d2tanh_two_h_dot_h_div_h_sat_squ *
 *            symmetrize(outer_product(dtwo_h_dot_h_div_h_sat_squ_dH,
 *                                     dtwo_h_dot_h_div_h_sat_squ_dH)) +
 *          dtanh_two_h_dot_h_div_h_sat_squ * d2two_h_dot_h_div_h_sat_squ_dH_dH);
 * 
 * @endcode
 * 
 * Some intermediate quantities attained directly from the
 * field / kinematic variables.
 * 
 * @code
 *       const double         log_det_F         = std::log(det_F);
 *       const double         tr_C              = trace(C);
 *       const Tensor<1, dim> C_inv_dot_H       = C_inv * H;
 *       const double         H_dot_C_inv_dot_H = H * C_inv_dot_H;
 * 
 * @endcode
 * 
 * First derivatives of the intermediate quantities.
 * 
 * @code
 *       const SymmetricTensor<2, dim> d_tr_C_dC =
 *         Physics::Elasticity::StandardTensors<dim>::I;
 *       const SymmetricTensor<2, dim> ddet_F_dC     = 0.5 * det_F * C_inv;
 *       const SymmetricTensor<2, dim> dlog_det_F_dC = 0.5 * C_inv;
 * 
 *       const Tensor<1, dim> dH_dot_C_inv_dot_H_dH = 2.0 * C_inv_dot_H;
 * 
 *       SymmetricTensor<4, dim> dC_inv_dC;
 *       for (unsigned int A = 0; A < dim; ++A)
 *         for (unsigned int B = A; B < dim; ++B)
 *           for (unsigned int C = 0; C < dim; ++C)
 *             for (unsigned int D = C; D < dim; ++D)
 *               dC_inv_dC[A][B][C][D] -=               
 *                 0.5 * (C_inv[A][C] * C_inv[B][D]     
 *                        + C_inv[A][D] * C_inv[B][C]); 
 * 
 *       const SymmetricTensor<2, dim> dH_dot_C_inv_dot_H_dC =
 *         -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H));
 * 
 * @endcode
 * 
 * Second derivatives of the intermediate quantities.
 * 
 * @code
 *       const SymmetricTensor<4, dim> d2log_det_F_dC_dC = 0.5 * dC_inv_dC;
 * 
 *       const SymmetricTensor<4, dim> d2det_F_dC_dC =
 *         0.5 * (outer_product(C_inv, ddet_F_dC) + det_F * dC_inv_dC);
 * 
 *       const SymmetricTensor<2, dim> d2H_dot_C_inv_dot_H_dH_dH = 2.0 * C_inv;
 * 
 *       Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
 *       for (unsigned int A = 0; A < dim; ++A)
 *         for (unsigned int B = 0; B < dim; ++B)
 *           for (unsigned int C = 0; C < dim; ++C)
 *             d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
 *               C_inv[A][C] * C_inv_dot_H[B] + 
 *               C_inv_dot_H[A] * C_inv[B][C];  
 * 
 *       SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
 *       for (unsigned int A = 0; A < dim; ++A)
 *         for (unsigned int B = A; B < dim; ++B)
 *           for (unsigned int C = 0; C < dim; ++C)
 *             for (unsigned int D = C; D < dim; ++D)
 *               d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
 *                 0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
 *                        C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
 *                        C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
 *                        C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);
 * 
 * @endcode
 * 
 * The stored energy density function.
 * 
 * @code
 *       psi =
 *         (0.5 * this->get_mu_e() * f_mu_e) *
 *           (tr_C - dim - 2.0 * std::log(det_F)) +
 *         this->get_lambda_e() * (std::log(det_F) * std::log(det_F)) -
 *         (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H);
 * 
 * @endcode
 * 
 * The kinetic quantities.
 * 
 * @code
 *       B = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) *
 *             df_mu_e_dH 
 *           + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
 *               dH_dot_C_inv_dot_H_dH; 
 * 
 *       S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                        
 *             (d_tr_C_dC - 2.0 * dlog_det_F_dC)                              
 *           + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC) 
 *           - 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *            
 *               (H_dot_C_inv_dot_H * ddet_F_dC                               
 *                + det_F * dH_dot_C_inv_dot_H_dC);                           
 * 
 * @endcode
 * 
 * The linearization of the kinetic quantities.
 * 
 * @code
 *       BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * 
 *              d2f_mu_e_dH_dH                                             
 *            + 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
 *                d2H_dot_C_inv_dot_H_dH_dH; 
 * 
 *       PP = -2.0 * (0.5 * this->get_mu_e()) *                                  
 *              outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),   
 *                            df_mu_e_dH)                                        
 *            +                                                                  
 *            2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                
 *              (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) 
 *               + det_F * d2H_dot_C_inv_dot_H_dC_dH);                           
 * 
 *       HH =
 *         4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) 
 *         + 4.0 * this->get_lambda_e() *                                       
 *             (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               
 *              + 2.0 * log_det_F * d2log_det_F_dC_dC)                          
 *         - 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *                
 *             (H_dot_C_inv_dot_H * d2det_F_dC_dC                               
 *              + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)               
 *              + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)               
 *              + det_F * d2H_dot_C_inv_dot_H_dC_dC);                           
 *     }
 * 
 *     template <int dim>
 *     double Magnetoelastic_Constitutive_Law<dim>::get_psi() const
 *     {
 *       return psi;
 *     }
 * 
 *     template <int dim>
 *     Tensor<1, dim> Magnetoelastic_Constitutive_Law<dim>::get_B() const
 *     {
 *       return B;
 *     }
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_S() const
 *     {
 *       return S;
 *     }
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim> Magnetoelastic_Constitutive_Law<dim>::get_DD() const
 *     {
 *       return BB;
 *     }
 * 
 *     template <int dim>
 *     Tensor<3, dim> Magnetoelastic_Constitutive_Law<dim>::get_PP() const
 *     {
 *       return PP;
 *     }
 * 
 *     template <int dim>
 *     SymmetricTensor<4, dim> Magnetoelastic_Constitutive_Law<dim>::get_HH() const
 *     {
 *       return HH;
 *     }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Magnetoviscoelasticconstitutivelawhandderived"></a> 
 * <h4>Magneto-viscoelastic constitutive law (hand-derived)</h4>
 * 

 * 
 * As mentioned before, the free energy density function for the
 * magneto-viscoelastic material with one dissipative mechanism that we'll
 * be considering is defined as
 * @f[
 * \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
 * \right)
 * = \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * + \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * @f]
 * @f[
 * \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * = \frac{1}{2} \mu_{e} f_{\mu_{e}^{ME}} \left( \boldsymbol{\mathbb{H}}
 * \right)
 * \left[ \text{tr}(\mathbf{C}) - d - 2 \ln (\text{det}(\mathbf{F}))
 * \right]
 * + \lambda_{e} \ln^{2} \left(\text{det}(\mathbf{F}) \right)
 * - \frac{1}{2} \mu_{0} \mu_{r} \text{det}(\mathbf{F})
 * \left[ \boldsymbol{\mathbb{H}} \cdot \mathbf{C}^{-1} \cdot
 * \boldsymbol{\mathbb{H}} \right]
 * @f]
 * @f[
 * \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * = \frac{1}{2} \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)
 * \left[ \mathbf{C}_{v} : \left[
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right] - d - \ln\left(
 * \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
 * @f]
 * with
 * @f[
 * f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{e}^{\infty}}{\mu_{e}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{e}^{\text{sat}}\right)^{2}} \right)
 * @f]
 * @f[
 * f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)
 * = 1 + \left[ \frac{\mu_{v}^{\infty}}{\mu_{v}} - 1 \right]
 * \tanh \left( 2 \frac{\boldsymbol{\mathbb{H}} \cdot
 * \boldsymbol{\mathbb{H}}}
 * {\left(h_{v}^{\text{sat}}\right)^{2}} \right)
 * @f]
 * and the evolution law
 * @f[
 * \dot{\mathbf{C}}_{v} \left( \mathbf{C} \right)
 * = \frac{1}{\tau} \left[
 * \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}\right]^{-1}
 * - \mathbf{C}_{v} \right]
 * @f]
 * that itself is parameterized in terms of $\mathbf{C}$.
 * By design, the magnetoelastic part of the energy
 * $\psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)$
 * is identical to that of the magnetoelastic material presented earlier.
 * So, for the derivatives of the various contributions stemming from this
 * part of the energy, please refer to the previous section. We'll continue
 * to highlight the specific contributions from those terms by
 * superscripting the salient terms with $ME$, while contributions from the
 * magneto-viscoelastic component are superscripted with $MVE$.
 * Furthermore, the magnetic saturation function
 * $f_{\mu_{v}}^{MVE} \left( \boldsymbol{\mathbb{H}} \right)$
 * for the damping term has the identical form as that of the elastic
 * term (i.e.,
 * $f_{\mu_{e}}^{ME} \left( \boldsymbol{\mathbb{H}} \right)$
 * ), and so the structure of its derivatives are identical to that
 * seen before; the only change is for the three constitutive parameters
 * that are now associated with the viscous shear modulus $\mu_{v}$ rather
 * than the elastic shear modulus $\mu_{e}$.
 *     

 * 
 * For this magneto-viscoelastic material, the first derivatives that
 * correspond to the magnetic induction vector and total Piola-Kirchhoff
 * stress tensor are
 * @f[
 * \boldsymbol{\mathbb{B}} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * \dealcoloneq - \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
 * \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} \equiv
 * \boldsymbol{\mathbb{B}}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right)
 * + \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) =  - \frac{d \psi_{0}^{ME} \left(
 * \mathbf{C}, \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}}}
 * - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
 * @f]
 * @f[
 * \mathbf{S}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * \dealcoloneq 2 \frac{\partial \psi_{0} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
 * \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} \equiv
 * \mathbf{S}^{\text{tot}, ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right)
 * + \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}}
 * \right)
 * =  2 \frac{d \psi_{0}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}}
 * \right)}{d \mathbf{C}}
 * + 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
 * @f]
 * with the viscous contributions being
 * @f[
 * \boldsymbol{\mathbb{B}}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * = - \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \boldsymbol{\mathbb{H}}}
 * \Big\vert_{\mathbf{C}, \mathbf{C}_{v}} = - \frac{1}{2} \mu_{v}
 * \left[ \mathbf{C}_{v} : \left[
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right] - d - \ln\left(
 * \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
 * \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{\partial \boldsymbol{\mathbb{H}}}
 * @f]
 * @f[
 * \mathbf{S}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}}
 * \right)
 * = 2 \frac{\partial \psi_{0}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)}{\partial \mathbf{C}}
 * \Big\vert_{\mathbf{C}_{v}, \boldsymbol{\mathbb{H}}} = \mu_{v}
 * f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
 * \frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}^{-1} \right]
 * + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}_{v}
 * \right]
 * @f]
 * and with
 * @f[
 * \frac{\partial f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{\partial \boldsymbol{\mathbb{H}}} \equiv \frac{d
 * f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)}{d
 * \boldsymbol{\mathbb{H}}} .
 * @f]
 * The time-discretized evolution law,
 * @f[
 * \mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)
 * = \frac{1}{1 + \frac{\Delta t}{\tau_{v}}} \left[
 * \mathbf{C}_{v}^{(t-1)}
 * + \frac{\Delta t}{\tau_{v}}
 * \left[\left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right]^{-1}
 * \right]
 * @f]
 * will also dictate how the linearization of the internal
 * variable with respect to the field variables is composed.
 *     

 * 
 * Observe that in order to attain the *correct* expressions for the
 * magnetic induction vector and total Piola-Kirchhoff stress tensor for
 * this dissipative material, we must adhere strictly to the outcome of
 * applying the Coleman-Noll procedure: we must take *partial derivatives*
 * of the free energy density function with respect to the field variables.
 * (For our non-dissipative magnetoelastic material, taking either partial
 * or total derivatives would have had the same result, so there was no
 * need to draw your attention to this before.)
 * The crucial part of the operation is to freeze the internal variable
 * $\mathbf{C}_{v}^{(t)} \left( \mathbf{C} \right)$ while computing the
 * derivatives of $\psi_{0}^{MVE} \left( \mathbf{C},
 * \mathbf{C}_{v} \left( \mathbf{C} \right), \boldsymbol{\mathbb{H}}
 * \right)$ with respect to $\mathbf{C}$ -- the dependence of
 * $\mathbf{C}_{v}^{(t)}$ on $\mathbf{C}$ is not to be taken into account.
 * When deciding whether to use AD or SD to perform this task
 * the choice is clear -- only the symbolic framework provides a mechanism
 * to do this; as was mentioned before, AD can only return total derivatives
 * so it is unsuitable for the task.
 *     

 * 
 * To wrap things up, we'll present the material tangents for this
 * rate-dependent coupled material. The linearization of both kinetic
 * variables with respect to their arguments are
 * @f[
 * \mathbb{D} \left( \mathbf{C}, \mathbf{C}_{v}, \boldsymbol{\mathbb{H}}
 * \right) = \frac{d \boldsymbol{\mathbb{B}}}{d \boldsymbol{\mathbb{H}}}
 * \equiv \mathbb{D}^{ME} \left( \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * + \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = \frac{d \boldsymbol{\mathbb{B}}^{ME}}{d
 * \boldsymbol{\mathbb{H}}}
 * + \frac{d \boldsymbol{\mathbb{B}}^{MVE}}{d \boldsymbol{\mathbb{H}}}
 * @f]
 * @f[
 * \mathfrak{P}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot}}}{d
 * \boldsymbol{\mathbb{H}}} \equiv \mathfrak{P}^{\text{tot}, ME} \left(
 * \mathbf{C}, \boldsymbol{\mathbb{H}} \right)
 * + \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = - \frac{d \mathbf{S}^{\text{tot},
 * ME}}{d \boldsymbol{\mathbb{H}}}
 * - \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \boldsymbol{\mathbb{H}}}
 * @f]
 * @f[
 * \mathcal{H}^{\text{tot}} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot}}}{d
 * \mathbf{C}} \equiv \mathcal{H}^{\text{tot}, ME} \left( \mathbf{C},
 * \boldsymbol{\mathbb{H}} \right)
 * + \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = 2 \frac{d \mathbf{S}^{\text{tot},
 * ME}}{d \mathbf{C}}
 * + 2 \frac{d \mathbf{S}^{\text{tot}, MVE}}{d \mathbf{C}}
 * @f]
 * where the tangents for the viscous contributions are
 * @f[
 * \mathbb{D}^{MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = - \frac{1}{2} \mu_{v}
 * \left[ \mathbf{C}_{v} : \left[
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C} \right] - d - \ln\left(
 * \text{det}\left(\mathbf{C}_{v}\right) \right)  \right]
 * \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{\partial \boldsymbol{\mathbb{H}} \otimes
 * d \boldsymbol{\mathbb{H}}}
 * @f]
 * @f[
 * \mathfrak{P}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right) = - \mu_{v}
 * \left[  \left[ \mathbf{C}_{v} : \mathbf{C} \right] \left[ -
 * \frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}^{-1} \right]
 * + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}_{v}
 * \right] \otimes \frac{d f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{d \boldsymbol{\mathbb{H}}}
 * @f]
 * @f{align}
 * \mathcal{H}^{\text{tot}, MVE} \left( \mathbf{C}, \mathbf{C}_{v},
 * \boldsymbol{\mathbb{H}} \right)
 * &= 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ - \frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}^{-1} \right] \otimes
 * \left[ \mathbf{C}_{v} + \mathbf{C} : \frac{d \mathbf{C}_{v}}{d
 * \mathbf{C}} \right]
 * \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[ \mathbf{C}_{v} : \mathbf{C} \right]
 * \left[
 * \frac{1}{d^{2}}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
 * - \frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}} \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}}
 * \right]
 * \\ &+ 2 \mu_{v} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}} \right)
 * \left[
 * -\frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \mathbf{C}_{v} \otimes \mathbf{C}^{-1}
 * + \left[\text{det}\left(\mathbf{F}\right)\right]^{-\frac{2}{d}}
 * \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
 * \right]
 * @f}
 * with
 * @f[
 * \frac{\partial^{2} f_{\mu_{v}^{MVE}} \left( \boldsymbol{\mathbb{H}}
 * \right)}{\partial \boldsymbol{\mathbb{H}} \otimes
 * d \boldsymbol{\mathbb{H}}} \equiv \frac{d^{2} f_{\mu_{v}^{MVE}} \left(
 * \boldsymbol{\mathbb{H}} \right)}{d \boldsymbol{\mathbb{H}} \otimes d
 * \boldsymbol{\mathbb{H}}}
 * @f]
 * and, from the evolution law,
 * @f[
 * \frac{d \mathbf{C}_{v}}{d \mathbf{C}}
 * \equiv \frac{d \mathbf{C}_{v}^{(t)}}{d \mathbf{C}}
 * = \frac{\frac{\Delta t}{\tau_{v}} }{1 + \frac{\Delta t}{\tau_{v}}}
 * \left[
 * \frac{1}{d}
 * \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}}
 * \mathbf{C}^{-1} \otimes \mathbf{C}^{-1}
 * + \left[\text{det}\left(\mathbf{F}\right)\right]^{\frac{2}{d}} \frac{d
 * \mathbf{C}^{-1}}{d \mathbf{C}}
 * \right] .
 * @f]
 * Notice that just the last term of $\mathcal{H}^{\text{tot}, MVE}$
 * contains the tangent of the internal variable. The linearization of this
 * particular evolution law is linear. For an example of a nonlinear
 * evolution law, for which this linearization must be solved for in an
 * iterative manner, see @cite Koprowski-Theiss2011a.
 * 
 * @code
 *     template <int dim>
 *     class Magnetoviscoelastic_Constitutive_Law final
 *       : public Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *     {
 *     public:
 *       Magnetoviscoelastic_Constitutive_Law(
 *         const ConstitutiveParameters &constitutive_parameters);
 * 
 *       virtual void update_internal_data(const SymmetricTensor<2, dim> &C,
 *                                         const Tensor<1, dim> &         H,
 *                                         const DiscreteTime &time) override;
 * 
 *       virtual double get_psi() const override;
 * 
 *       virtual Tensor<1, dim> get_B() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_S() const override;
 * 
 *       virtual SymmetricTensor<2, dim> get_DD() const override;
 * 
 *       virtual Tensor<3, dim> get_PP() const override;
 * 
 *       virtual SymmetricTensor<4, dim> get_HH() const override;
 * 
 *       virtual void update_end_of_timestep() override;
 * 
 *     private:
 *       SymmetricTensor<2, dim> Q_t;
 *       SymmetricTensor<2, dim> Q_t1;
 * 
 *       double                  psi;
 *       Tensor<1, dim>          B;
 *       SymmetricTensor<2, dim> S;
 *       SymmetricTensor<2, dim> BB;
 *       Tensor<3, dim>          PP;
 *       SymmetricTensor<4, dim> HH;
 * 
 * @endcode
 * 
 * A data structure that is used to store all intermediate calculations.
 * We'll see shortly precisely how this can be leveraged to make the part
 * of the code where we actually perform calculations clean and easy
 * (well, at least easier) to follow and maintain. But for now, we can say
 * that it will allow us to move the parts of the code where we compute
 * the derivatives of intermediate quantities away from where they are
 * used.
 * 
 * @code
 *       mutable GeneralDataStorage cache;
 * 
 * @endcode
 * 
 * The next two functions are used to update the state of the field and
 * internal variables, and will be called before we perform any
 * detailed calculations.
 * 
 * @code
 *       void set_primary_variables(const SymmetricTensor<2, dim> &C,
 *                                  const Tensor<1, dim> &         H) const;
 * 
 *       void update_internal_variable(const DiscreteTime &time);
 * 
 * @endcode
 * 
 * The remainder of the class interface is dedicated to methods that
 * are used to compute the components required to calculate the free
 * energy density function, and all of its derivatives:
 * 

 * 
 * The kinematic, or field, variables.
 * 
 * @code
 *       const Tensor<1, dim> &get_H() const;
 * 
 *       const SymmetricTensor<2, dim> &get_C() const;
 * 
 * @endcode
 * 
 * A generalized formulation for the saturation function, with the
 * required constitutive parameters passed as arguments to each function.
 * 
 * @code
 *       double get_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;
 * 
 *       double get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;
 * 
 *       double get_f_mu(const double mu,
 *                       const double mu_inf,
 *                       const double mu_h_sat) const;
 * 
 * @endcode
 * 
 * A generalized formulation for the first derivative of saturation
 * function, with the required constitutive parameters passed as arguments
 * to each function.
 * 
 * @code
 *       double get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;
 * 
 *       Tensor<1, dim>
 *       get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const;
 * 
 *       Tensor<1, dim> get_df_mu_dH(const double mu,
 *                                   const double mu_inf,
 *                                   const double mu_h_sat) const;
 * 
 * @endcode
 * 
 * A generalized formulation for the second derivative of saturation
 * function, with the required constitutive parameters passed as arguments
 * to each function.
 * 
 * @code
 *       double get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const;
 * 
 *       SymmetricTensor<2, dim>
 *       get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const;
 * 
 *       SymmetricTensor<2, dim> get_d2f_mu_dH_dH(const double mu,
 *                                                const double mu_inf,
 *                                                const double mu_h_sat) const;
 * 
 * @endcode
 * 
 * Intermediate quantities attained directly from the
 * field / kinematic variables.
 * 
 * @code
 *       const double &get_det_F() const;
 * 
 *       const SymmetricTensor<2, dim> &get_C_inv() const;
 * 
 *       const double &get_log_det_F() const;
 * 
 *       const double &get_trace_C() const;
 * 
 *       const Tensor<1, dim> &get_C_inv_dot_H() const;
 * 
 *       const double &get_H_dot_C_inv_dot_H() const;
 * 
 * @endcode
 * 
 * First derivatives of the intermediate quantities.
 * 
 * @code
 *       const SymmetricTensor<4, dim> &get_dC_inv_dC() const;
 * 
 *       const SymmetricTensor<2, dim> &get_d_tr_C_dC() const;
 * 
 *       const SymmetricTensor<2, dim> &get_ddet_F_dC() const;
 * 
 *       const SymmetricTensor<2, dim> &get_dlog_det_F_dC() const;
 * 
 *       const Tensor<1, dim> &get_dH_dot_C_inv_dot_H_dH() const;
 * 
 *       const SymmetricTensor<2, dim> &get_dH_dot_C_inv_dot_H_dC() const;
 * 
 * @endcode
 * 
 * Derivative of internal variable with respect to field variables.
 * Notice that we only need this one derivative of the internal variable,
 * as this variable is only differentiated as part of the linearization
 * of the kinetic variables.
 * 
 * @code
 *       const SymmetricTensor<4, dim> &
 *       get_dQ_t_dC(const DiscreteTime &time) const;
 * 
 * @endcode
 * 
 * Second derivatives of the intermediate quantities.
 * 
 * @code
 *       const SymmetricTensor<4, dim> &get_d2log_det_F_dC_dC() const;
 * 
 *       const SymmetricTensor<4, dim> &get_d2det_F_dC_dC() const;
 * 
 *       const SymmetricTensor<2, dim> &get_d2H_dot_C_inv_dot_H_dH_dH() const;
 * 
 *       const Tensor<3, dim> &get_d2H_dot_C_inv_dot_H_dC_dH() const;
 * 
 *       const SymmetricTensor<4, dim> &get_d2H_dot_C_inv_dot_H_dC_dC() const;
 *     };
 * 
 * 
 *     template <int dim>
 *     Magnetoviscoelastic_Constitutive_Law<
 *       dim>::Magnetoviscoelastic_Constitutive_Law(const ConstitutiveParameters
 *                                                    &constitutive_parameters)
 *       : Coupled_Magnetomechanical_Constitutive_Law_Base<dim>(
 *           constitutive_parameters)
 *       , Q_t(Physics::Elasticity::StandardTensors<dim>::I)
 *       , Q_t1(Physics::Elasticity::StandardTensors<dim>::I)
 *       , psi(0.0)
 *     {}
 * 
 * 
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_data(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H,
 *       const DiscreteTime &           time)
 *     {
 * @endcode
 * 
 * Record the applied deformation state as well as the magnetic load.
 * Thereafter, update internal (viscous) variable based on new deformation
 * state.
 * 
 * @code
 *       set_primary_variables(C, H);
 *       update_internal_variable(time);
 * 
 * @endcode
 * 
 * Get the values for the elastic and viscous saturation function based
 * on the current magnetic field...
 * 
 * @code
 *       const double f_mu_e = get_f_mu(this->get_mu_e(),
 *                                      this->get_mu_e_inf(),
 *                                      this->get_mu_e_h_sat());
 * 
 *       const double f_mu_v = get_f_mu(this->get_mu_v(),
 *                                      this->get_mu_v_inf(),
 *                                      this->get_mu_v_h_sat());
 * 
 * @endcode
 * 
 * ... as well as their first derivatives...
 * 
 * @code
 *       const Tensor<1, dim> df_mu_e_dH = get_df_mu_dH(this->get_mu_e(),
 *                                                      this->get_mu_e_inf(),
 *                                                      this->get_mu_e_h_sat());
 * 
 *       const Tensor<1, dim> df_mu_v_dH = get_df_mu_dH(this->get_mu_v(),
 *                                                      this->get_mu_v_inf(),
 *                                                      this->get_mu_v_h_sat());
 * 
 * 
 * @endcode
 * 
 * ... and their second derivatives.
 * 
 * @code
 *       const SymmetricTensor<2, dim> d2f_mu_e_dH_dH =
 *         get_d2f_mu_dH_dH(this->get_mu_e(),
 *                          this->get_mu_e_inf(),
 *                          this->get_mu_e_h_sat());
 * 
 *       const SymmetricTensor<2, dim> d2f_mu_v_dH_dH =
 *         get_d2f_mu_dH_dH(this->get_mu_v(),
 *                          this->get_mu_v_inf(),
 *                          this->get_mu_v_h_sat());
 * 
 * @endcode
 * 
 * Intermediate quantities. Note that, since we're fetching these values
 * from a cache that has a lifetime that outlasts this function call, we
 * can alias the result rather than copying the value from the cache.
 * 
 * @code
 *       const double &                 det_F = get_det_F();
 *       const SymmetricTensor<2, dim> &C_inv = get_C_inv();
 * 
 *       const double &log_det_F         = get_log_det_F();
 *       const double &tr_C              = get_trace_C();
 *       const double &H_dot_C_inv_dot_H = get_H_dot_C_inv_dot_H();
 * 
 * @endcode
 * 
 * First derivatives of intermediate values, as well as the that of the
 * internal variable with respect to the right Cauchy-Green deformation
 * tensor.
 * 
 * @code
 *       const SymmetricTensor<2, dim> &d_tr_C_dC     = get_d_tr_C_dC();
 *       const SymmetricTensor<2, dim> &ddet_F_dC     = get_ddet_F_dC();
 *       const SymmetricTensor<2, dim> &dlog_det_F_dC = get_dlog_det_F_dC();
 * 
 *       const SymmetricTensor<4, dim> &dQ_t_dC = get_dQ_t_dC(time);
 * 
 *       const Tensor<1, dim> &dH_dot_C_inv_dot_H_dH = get_dH_dot_C_inv_dot_H_dH();
 * 
 *       const SymmetricTensor<2, dim> &dH_dot_C_inv_dot_H_dC =
 *         get_dH_dot_C_inv_dot_H_dC();
 * 
 * @endcode
 * 
 * Second derivatives of intermediate values.
 * 
 * @code
 *       const SymmetricTensor<4, dim> &d2log_det_F_dC_dC =
 *         get_d2log_det_F_dC_dC();
 * 
 *       const SymmetricTensor<4, dim> &d2det_F_dC_dC = get_d2det_F_dC_dC();
 * 
 *       const SymmetricTensor<2, dim> &d2H_dot_C_inv_dot_H_dH_dH =
 *         get_d2H_dot_C_inv_dot_H_dH_dH();
 * 
 *       const Tensor<3, dim> &d2H_dot_C_inv_dot_H_dC_dH =
 *         get_d2H_dot_C_inv_dot_H_dC_dH();
 * 
 *       const SymmetricTensor<4, dim> &d2H_dot_C_inv_dot_H_dC_dC =
 *         get_d2H_dot_C_inv_dot_H_dC_dC();
 * 
 * @endcode
 * 
 * Since the definitions of the linearizations become particularly
 * lengthy, we'll decompose the free energy density function into three
 * additive components:
 * - the "Neo-Hookean"-like term,
 * - the rate-dependent term, and
 * - the term that resembles that of the energy stored in the magnetic
 * field.
 *       

 * 
 * To remain consistent, each of these contributions will be individually
 * added to the variables that we want to compute in that same order.
 *       

 * 
 * So, first of all this is the energy density function itself:
 * 
 * @code
 *       psi = (0.5 * this->get_mu_e() * f_mu_e) *
 *               (tr_C - dim - 2.0 * std::log(det_F)) +
 *             this->get_lambda_e() * (std::log(det_F) * std::log(det_F));
 *       psi += (0.5 * this->get_mu_v() * f_mu_v) *
 *              (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
 *               std::log(determinant(Q_t)));
 *       psi -=
 *         (0.5 * this->get_mu_0() * this->get_mu_r()) * det_F * (H * C_inv * H);
 * 
 * @endcode
 * 
 * ... followed by the magnetic induction vector and Piola-Kirchhoff
 * stress:
 * 
 * @code
 *       B =
 *         -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) * df_mu_e_dH;
 *       B -= (0.5 * this->get_mu_v()) *
 *            (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
 *             std::log(determinant(Q_t))) *
 *            df_mu_v_dH;
 *       B += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
 *            dH_dot_C_inv_dot_H_dH;
 * 
 *       S = 2.0 * (0.5 * this->get_mu_e() * f_mu_e) *                         
 *             (d_tr_C_dC - 2.0 * dlog_det_F_dC)                               
 *           + 2.0 * this->get_lambda_e() * (2.0 * log_det_F * dlog_det_F_dC); 
 *       S += 2.0 * (0.5 * this->get_mu_v() * f_mu_v) *
 *            ((Q_t * C) *
 *               ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) * ddet_F_dC) +
 *             std::pow(det_F, -2.0 / dim) * Q_t);                // dC/dC = II
 *       S -= 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * 
 *            (H_dot_C_inv_dot_H * ddet_F_dC                      
 *             + det_F * dH_dot_C_inv_dot_H_dC);                  
 * 
 * @endcode
 * 
 * ... and lastly the tangents due to the linearization of the kinetic
 * variables.
 * 
 * @code
 *       BB = -(0.5 * this->get_mu_e() * (tr_C - dim - 2.0 * log_det_F)) *
 *            d2f_mu_e_dH_dH;
 *       BB -= (0.5 * this->get_mu_v()) *
 *             (Q_t * (std::pow(det_F, -2.0 / dim) * C) - dim -
 *              std::log(determinant(Q_t))) *
 *             d2f_mu_v_dH_dH;
 *       BB += 0.5 * this->get_mu_0() * this->get_mu_r() * det_F *
 *             d2H_dot_C_inv_dot_H_dH_dH;
 * 
 *       PP = -2.0 * (0.5 * this->get_mu_e()) *
 *            outer_product(Tensor<2, dim>(d_tr_C_dC - 2.0 * dlog_det_F_dC),
 *                          df_mu_e_dH);
 *       PP -= 2.0 * (0.5 * this->get_mu_v()) *
 *             outer_product(Tensor<2, dim>((Q_t * C) *
 *                                            ((-2.0 / dim) *
 *                                             std::pow(det_F, -2.0 / dim - 1.0) *
 *                                             ddet_F_dC) +
 *                                          std::pow(det_F, -2.0 / dim) * Q_t),
 *                           df_mu_v_dH);
 *       PP += 2.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) *
 *             (outer_product(Tensor<2, dim>(ddet_F_dC), dH_dot_C_inv_dot_H_dH) +
 *              det_F * d2H_dot_C_inv_dot_H_dC_dH);
 * 
 *       HH =
 *         4.0 * (0.5 * this->get_mu_e() * f_mu_e) * (-2.0 * d2log_det_F_dC_dC) 
 *         + 4.0 * this->get_lambda_e() *                                       
 *             (2.0 * outer_product(dlog_det_F_dC, dlog_det_F_dC)               
 *              + 2.0 * log_det_F * d2log_det_F_dC_dC);                         
 *       HH += 4.0 * (0.5 * this->get_mu_v() * f_mu_v) *
 *             (outer_product((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
 *                              ddet_F_dC,
 *                            C * dQ_t_dC + Q_t) +
 *              (Q_t * C) *
 *                (outer_product(ddet_F_dC,
 *                               (-2.0 / dim) * (-2.0 / dim - 1.0) *
 *                                 std::pow(det_F, -2.0 / dim - 2.0) * ddet_F_dC) +
 *                 ((-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
 *                  d2det_F_dC_dC)) +
 *              outer_product(Q_t,
 *                            (-2.0 / dim) * std::pow(det_F, -2.0 / dim - 1.0) *
 *                              ddet_F_dC) +
 *              std::pow(det_F, -2.0 / dim) * dQ_t_dC);
 *       HH -= 4.0 * (0.5 * this->get_mu_0() * this->get_mu_r()) * 
 *             (H_dot_C_inv_dot_H * d2det_F_dC_dC                  
 *              + outer_product(ddet_F_dC, dH_dot_C_inv_dot_H_dC)  
 *              + outer_product(dH_dot_C_inv_dot_H_dC, ddet_F_dC)  
 *              + det_F * d2H_dot_C_inv_dot_H_dC_dC);              
 * 
 * 
 * @endcode
 * 
 * Now that we're done using all of those temporary variables stored
 * in our cache, we can clear it out to free up some memory.
 * 
 * @code
 *       cache.reset();
 *     }
 * 
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law<dim>::get_psi() const
 *     {
 *       return psi;
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_B() const
 *     {
 *       return B;
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim>
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_S() const
 *     {
 *       return S;
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim>
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_DD() const
 *     {
 *       return BB;
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<3, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_PP() const
 *     {
 *       return PP;
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<4, dim>
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_HH() const
 *     {
 *       return HH;
 *     }
 * 
 * 
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law<dim>::update_end_of_timestep()
 *     {
 *       Q_t1 = Q_t;
 *     }
 * 
 * 
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law<dim>::update_internal_variable(
 *       const DiscreteTime &time)
 *     {
 *       const double delta_t = this->get_delta_t(time);
 * 
 *       Q_t = (1.0 / (1.0 + delta_t / this->get_tau_v())) *
 *             (Q_t1 + (delta_t / this->get_tau_v()) *
 *                       std::pow(get_det_F(), 2.0 / dim) * get_C_inv());
 *     }
 * 
 * @endcode
 * 
 * The next few functions implement the generalized formulation for the
 * saturation function, as well as its various derivatives.
 * 
 * @code
 *     template <int dim>
 *     double
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_two_h_dot_h_div_h_sat_squ(
 *       const double mu_h_sat) const
 *     {
 *       const Tensor<1, dim> &H = get_H();
 *       return (2.0 * H * H) / (mu_h_sat * mu_h_sat);
 *     }
 * 
 * 
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law<
 *       dim>::get_tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
 *     {
 *       return std::tanh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat));
 *     }
 * 
 * @endcode
 * 
 * A scaling function that will cause the shear modulus
 * to change (increase) under the influence of a magnetic
 * field.
 * 
 * @code
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law<dim>::get_f_mu(
 *       const double mu,
 *       const double mu_inf,
 *       const double mu_h_sat) const
 *     {
 *       return 1.0 +
 *              (mu_inf / mu - 1.0) * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
 *     }
 * 
 * @endcode
 * 
 * First derivative of scaling function
 * 
 * @code
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law<
 *       dim>::get_dtanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
 *     {
 *       return std::pow(1.0 / std::cosh(get_two_h_dot_h_div_h_sat_squ(mu_h_sat)),
 *                       2.0);
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<
 *       dim>::get_dtwo_h_dot_h_div_h_sat_squ_dH(const double mu_h_sat) const
 *     {
 *       return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) * get_H();
 *     }
 * 
 * 
 *     template <int dim>
 *     Tensor<1, dim> Magnetoviscoelastic_Constitutive_Law<dim>::get_df_mu_dH(
 *       const double mu,
 *       const double mu_inf,
 *       const double mu_h_sat) const
 *     {
 *       return (mu_inf / mu - 1.0) *
 *              (get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
 *               get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat));
 *     }
 * 
 * 
 *     template <int dim>
 *     double Magnetoviscoelastic_Constitutive_Law<
 *       dim>::get_d2tanh_two_h_dot_h_div_h_sat_squ(const double mu_h_sat) const
 *     {
 *       return -2.0 * get_tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
 *              get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat);
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim> Magnetoviscoelastic_Constitutive_Law<
 *       dim>::get_d2two_h_dot_h_div_h_sat_squ_dH_dH(const double mu_h_sat) const
 *     {
 *       return 2.0 * 2.0 / (mu_h_sat * mu_h_sat) *
 *              Physics::Elasticity::StandardTensors<dim>::I;
 *     }
 * 
 * 
 *     template <int dim>
 *     SymmetricTensor<2, dim>
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2f_mu_dH_dH(
 *       const double mu,
 *       const double mu_inf,
 *       const double mu_h_sat) const
 *     {
 *       return (mu_inf / mu - 1.0) *
 *              (get_d2tanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
 *                 symmetrize(
 *                   outer_product(get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat),
 *                                 get_dtwo_h_dot_h_div_h_sat_squ_dH(mu_h_sat))) +
 *               get_dtanh_two_h_dot_h_div_h_sat_squ(mu_h_sat) *
 *                 get_d2two_h_dot_h_div_h_sat_squ_dH_dH(mu_h_sat));
 *     }
 * 
 * @endcode
 * 
 * For the cached calculation approach that we've adopted for this material
 * class, the root of all calculations are the field variables, and the
 * immutable ancillary data such as the constitutive parameters and time
 * step size. As such, we need to enter them into the cache in a different
 * manner to the other variables, since they are inputs that are prescribed
 * from outside the class itself. This function simply adds them to the
 * cache directly from the input arguments, checking that there is no
 * equivalent data there in the first place (we expect to call the
 * `update_internal_data()` method only once per time step, or Newton
 * iteration).
 * 
 * @code
 *     template <int dim>
 *     void Magnetoviscoelastic_Constitutive_Law<dim>::set_primary_variables(
 *       const SymmetricTensor<2, dim> &C,
 *       const Tensor<1, dim> &         H) const
 *     {
 * @endcode
 * 
 * Set value for $\boldsymbol{\mathbb{H}}$.
 * 
 * @code
 *       const std::string name_H("H");
 *       Assert(!cache.stores_object_with_name(name_H),
 *              ExcMessage(
 *                "The primary variable has already been added to the cache."));
 *       cache.add_unique_copy(name_H, H);
 * 
 * @endcode
 * 
 * Set value for $\mathbf{C}$.
 * 
 * @code
 *       const std::string name_C("C");
 *       Assert(!cache.stores_object_with_name(name_C),
 *              ExcMessage(
 *                "The primary variable has already been added to the cache."));
 *       cache.add_unique_copy(name_C, C);
 *     }
 * 
 * @endcode
 * 
 * After that, we can fetch them from the cache at any point in time.
 * 
 * @code
 *     template <int dim>
 *     const Tensor<1, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_H() const
 *     {
 *       const std::string name("H");
 *       Assert(cache.stores_object_with_name(name),
 *              ExcMessage("Primary variables must be added to the cache."));
 *       return cache.template get_object_with_name<Tensor<1, dim>>(name);
 *     }
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_C() const
 *     {
 *       const std::string name("C");
 *       Assert(cache.stores_object_with_name(name),
 *              ExcMessage("Primary variables must be added to the cache."));
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * @endcode
 * 
 * With the primary variables guaranteed to be in the cache when we need
 * them, we can not compute all intermediate values (either directly, or
 * indirectly) from them.
 *     

 * 
 * If the cache does not already store the value that we're looking for,
 * then we quickly calculate it, store it in the cache and return the value
 * just stored in the cache. That way we can return it as a reference and
 * avoid copying the object. The same goes for any values that a compound
 * function might depend on. Said another way, if there is a dependency
 * chain of calculations that come before the one that we're currently
 * interested in doing, then we're guaranteed to resolve the dependencies
 * before we proceed with using any of those values. Although there is a
 * cost to fetching data from the cache, the "resolved dependency" concept
 * might be sufficiently convenient to make it worth looking past the extra
 * cost. If these material laws are embedded within a finite element
 * framework, then the added cost might not even be noticeable.
 * 
 * @code
 *     template <int dim>
 *     const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_det_F() const
 *     {
 *       const std::string name("det_F");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const double det_F = std::sqrt(determinant(get_C()));
 *           AssertThrow(det_F > 0.0,
 *                       ExcMessage("Volumetric Jacobian must be positive."));
 *           cache.add_unique_copy(name, det_F);
 *         }
 * 
 *       return cache.template get_object_with_name<double>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv() const
 *     {
 *       const std::string name("C_inv");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           cache.add_unique_copy(name, invert(get_C()));
 *         }
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const double &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_log_det_F() const
 *     {
 *       const std::string name("log(det_F)");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, std::log(get_det_F()));
 * 
 *       return cache.template get_object_with_name<double>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const double &Magnetoviscoelastic_Constitutive_Law<dim>::get_trace_C() const
 *     {
 *       const std::string name("trace(C)");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, trace(get_C()));
 * 
 *       return cache.template get_object_with_name<double>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const Tensor<1, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_C_inv_dot_H() const
 *     {
 *       const std::string name("C_inv_dot_H");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, get_C_inv() * get_H());
 * 
 *       return cache.template get_object_with_name<Tensor<1, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const double &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_H_dot_C_inv_dot_H() const
 *     {
 *       const std::string name("H_dot_C_inv_dot_H");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, get_H() * get_C_inv_dot_H());
 * 
 *       return cache.template get_object_with_name<double>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<4, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_dQ_t_dC(
 *       const DiscreteTime &time) const
 *     {
 *       const std::string name("dQ_t_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const double  delta_t = this->get_delta_t(time);
 *           const double &det_F   = get_det_F();
 * 
 *           const SymmetricTensor<4, dim> dQ_t_dC =
 *             (1.0 / (1.0 + delta_t / this->get_tau_v())) *
 *             (delta_t / this->get_tau_v()) *
 *             ((2.0 / dim) * std::pow(det_F, 2.0 / dim - 1.0) *
 *                outer_product(get_C_inv(), get_ddet_F_dC()) +
 *              std::pow(det_F, 2.0 / dim) * get_dC_inv_dC());
 * 
 *           cache.add_unique_copy(name, dQ_t_dC);
 *         }
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<4, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_dC_inv_dC() const
 *     {
 *       const std::string name("dC_inv_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const SymmetricTensor<2, dim> &C_inv = get_C_inv();
 *           SymmetricTensor<4, dim>        dC_inv_dC;
 * 
 *           for (unsigned int A = 0; A < dim; ++A)
 *             for (unsigned int B = A; B < dim; ++B)
 *               for (unsigned int C = 0; C < dim; ++C)
 *                 for (unsigned int D = C; D < dim; ++D)
 *                   dC_inv_dC[A][B][C][D] -=               
 *                     0.5 * (C_inv[A][C] * C_inv[B][D]     
 *                            + C_inv[A][D] * C_inv[B][C]); 
 * 
 *           cache.add_unique_copy(name, dC_inv_dC);
 *         }
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d_tr_C_dC() const
 *     {
 *       const std::string name("d_tr_C_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name,
 *                               Physics::Elasticity::StandardTensors<dim>::I);
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_ddet_F_dC() const
 *     {
 *       const std::string name("ddet_F_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, 0.5 * get_det_F() * get_C_inv());
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_dlog_det_F_dC() const
 *     {
 *       const std::string name("dlog_det_F_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, 0.5 * get_C_inv());
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const Tensor<1, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dH() const
 *     {
 *       const std::string name("dH_dot_C_inv_dot_H_dH");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, 2.0 * get_C_inv_dot_H());
 * 
 *       return cache.template get_object_with_name<Tensor<1, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_dH_dot_C_inv_dot_H_dC() const
 *     {
 *       const std::string name("dH_dot_C_inv_dot_H_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const Tensor<1, dim> C_inv_dot_H = get_C_inv_dot_H();
 *           cache.add_unique_copy(
 *             name, -symmetrize(outer_product(C_inv_dot_H, C_inv_dot_H)));
 *         }
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<4, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2log_det_F_dC_dC() const
 *     {
 *       const std::string name("d2log_det_F_dC_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, 0.5 * get_dC_inv_dC());
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<4, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2det_F_dC_dC() const
 *     {
 *       const std::string name("d2det_F_dC_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name,
 *                               0.5 *
 *                                 (outer_product(get_C_inv(), get_ddet_F_dC()) +
 *                                  get_det_F() * get_dC_inv_dC()));
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<2, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dH_dH()
 *       const
 *     {
 *       const std::string name("d2H_dot_C_inv_dot_H_dH_dH");
 *       if (cache.stores_object_with_name(name) == false)
 *         cache.add_unique_copy(name, 2.0 * get_C_inv());
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<2, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const Tensor<3, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dH()
 *       const
 *     {
 *       const std::string name("d2H_dot_C_inv_dot_H_dC_dH");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
 *           const SymmetricTensor<2, dim> &C_inv       = get_C_inv();
 * 
 *           Tensor<3, dim> d2H_dot_C_inv_dot_H_dC_dH;
 *           for (unsigned int A = 0; A < dim; ++A)
 *             for (unsigned int B = 0; B < dim; ++B)
 *               for (unsigned int C = 0; C < dim; ++C)
 *                 d2H_dot_C_inv_dot_H_dC_dH[A][B][C] -=
 *                   C_inv[A][C] * C_inv_dot_H[B] + 
 *                   C_inv_dot_H[A] * C_inv[B][C];  
 * 
 *           cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dH);
 *         }
 * 
 *       return cache.template get_object_with_name<Tensor<3, dim>>(name);
 *     }
 * 
 * 
 *     template <int dim>
 *     const SymmetricTensor<4, dim> &
 *     Magnetoviscoelastic_Constitutive_Law<dim>::get_d2H_dot_C_inv_dot_H_dC_dC()
 *       const
 *     {
 *       const std::string name("d2H_dot_C_inv_dot_H_dC_dC");
 *       if (cache.stores_object_with_name(name) == false)
 *         {
 *           const Tensor<1, dim> &         C_inv_dot_H = get_C_inv_dot_H();
 *           const SymmetricTensor<2, dim> &C_inv       = get_C_inv();
 * 
 *           SymmetricTensor<4, dim> d2H_dot_C_inv_dot_H_dC_dC;
 *           for (unsigned int A = 0; A < dim; ++A)
 *             for (unsigned int B = A; B < dim; ++B)
 *               for (unsigned int C = 0; C < dim; ++C)
 *                 for (unsigned int D = C; D < dim; ++D)
 *                   d2H_dot_C_inv_dot_H_dC_dC[A][B][C][D] +=
 *                     0.5 * (C_inv_dot_H[A] * C_inv_dot_H[C] * C_inv[B][D] +
 *                            C_inv_dot_H[A] * C_inv_dot_H[D] * C_inv[B][C] +
 *                            C_inv_dot_H[B] * C_inv_dot_H[C] * C_inv[A][D] +
 *                            C_inv_dot_H[B] * C_inv_dot_H[D] * C_inv[A][C]);
 * 
 *           cache.add_unique_copy(name, d2H_dot_C_inv_dot_H_dC_dC);
 *         }
 * 
 *       return cache.template get_object_with_name<SymmetricTensor<4, dim>>(name);
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="Rheologicalexperimentparameters"></a> 
 * <h4>Rheological experiment parameters</h4>
 * 

 * 
 * The @p RheologicalExperimentParameters class is used to drive the
 * numerical experiments that are to be conducted on the coupled materials
 * that we've implemented constitutive laws for.
 * 
 * @code
 *     class RheologicalExperimentParameters : public ParameterAcceptor
 *     {
 *     public:
 *       RheologicalExperimentParameters();
 * 
 * @endcode
 * 
 * These are  dimensions of the rheological specimen that is to be
 * simulated. They, effectively, define the measurement point for our
 * virtual experiment.
 * 
 * @code
 *       double sample_radius = 0.01;
 *       double sample_height = 0.001;
 * 
 * @endcode
 * 
 * The three steady-state loading parameters are respectively
 * - the axial stretch,
 * - the shear strain amplitude, and
 * - the axial magnetic field strength.
 * 
 * @code
 *       double lambda_2 = 0.95;
 *       double gamma_12 = 0.05;
 *       double H_2      = 60.0e3;
 * 
 * @endcode
 * 
 * Moreover, the parameters for the time-dependent rheological
 * loading conditions are
 * - the loading cycle frequency,
 * - the number of load cycles, and
 * - the number of discrete timesteps per cycle.
 * 
 * @code
 *       double       frequency         = 1.0 / (2.0 * numbers::PI);
 *       unsigned int n_cycles          = 5;
 *       unsigned int n_steps_per_cycle = 2500;
 * 
 * @endcode
 * 
 * We also declare some self-explanatory parameters related to output
 * data generated for the experiments conducted with rate-dependent and
 * rate-independent materials.
 * 
 * @code
 *       bool        output_data_to_file = true;
 *       std::string output_filename_rd =
 *         "experimental_results-rate_dependent.csv";
 *       std::string output_filename_ri =
 *         "experimental_results-rate_independent.csv";
 * 
 * @endcode
 * 
 * The next few functions compute time-related parameters for the
 * experiment...
 * 
 * @code
 *       double start_time() const;
 * 
 *       double end_time() const;
 * 
 *       double delta_t() const;
 * 
 * @endcode
 * 
 * ... while the following two prescribe the mechanical and magnetic
 * loading at any given time...
 * 
 * @code
 *       Tensor<1, 3> get_H(const double time) const;
 * 
 *       Tensor<2, 3> get_F(const double time) const;
 * 
 * @endcode
 * 
 * ... and this last one outputs the status of the experiment to the
 * console.
 * 
 * @code
 *       bool print_status(const int step_number) const;
 * 
 *       bool initialized = false;
 *     };
 * 
 * 
 * 
 *     RheologicalExperimentParameters::RheologicalExperimentParameters()
 *       : ParameterAcceptor("/Coupled Constitutive Laws/Rheological Experiment/")
 *     {
 *       add_parameter("Experimental sample radius", sample_radius);
 *       add_parameter("Experimental sample radius", sample_height);
 * 
 *       add_parameter("Axial stretch", lambda_2);
 *       add_parameter("Shear strain amplitude", gamma_12);
 *       add_parameter("Axial magnetic field strength", H_2);
 * 
 *       add_parameter("Frequency", frequency);
 *       add_parameter("Number of loading cycles", n_cycles);
 *       add_parameter("Discretisation for each cycle", n_steps_per_cycle);
 * 
 *       add_parameter("Output experimental results to file", output_data_to_file);
 *       add_parameter("Output file name (rate dependent constitutive law)",
 *                     output_filename_rd);
 *       add_parameter("Output file name (rate independent constitutive law)",
 *                     output_filename_ri);
 * 
 *       parse_parameters_call_back.connect([&]() -> void { initialized = true; });
 *     }
 * 
 * 
 *     double RheologicalExperimentParameters::start_time() const
 *     {
 *       return 0.0;
 *     }
 * 
 * 
 *     double RheologicalExperimentParameters::end_time() const
 *     {
 *       return n_cycles / frequency;
 *     }
 * 
 * 
 *     double RheologicalExperimentParameters::delta_t() const
 *     {
 *       return (end_time() - start_time()) / (n_steps_per_cycle * n_cycles);
 *     }
 * 
 * 
 *     bool
 *     RheologicalExperimentParameters::print_status(const int step_number) const
 *     {
 *       return (step_number % (n_cycles * n_steps_per_cycle / 100)) == 0;
 *     }
 * 
 * @endcode
 * 
 * The applied magnetic field is always aligned with the axis of rotation
 * of the rheometer's rotor.
 * 
 * @code
 *     Tensor<1, 3> RheologicalExperimentParameters::get_H(const double) const
 *     {
 *       return Tensor<1, 3>({0.0, 0.0, H_2});
 *     }
 * 
 * @endcode
 * 
 * The applied deformation (gradient) is computed based on the geometry
 * of the rheometer and the sample, the sampling point, and the experimental
 * parameters. From the displacement profile documented in the introduction,
 * the deformation gradient may be expressed in Cartesian coordinates as
 * @f[
 * \mathbf{F} = \begin{bmatrix}
 * \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
 * & -\frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
 * & -\tau R \sqrt{\lambda_{3}} \sin\left(\Theta + \alpha\right)
 * \\  \frac{\sin\left(\alpha\right)}{\sqrt{\lambda_{3}}}
 * & \frac{\cos\left(\alpha\right)}{\sqrt{\lambda_{3}}}
 * & -\tau R \sqrt{\lambda_{3}} \cos\left(\Theta + \alpha\right)
 * \\  0 & 0 & \lambda_{3}
 * \end{bmatrix}
 * @f]
 * 
 * @code
 *     Tensor<2, 3> RheologicalExperimentParameters::get_F(const double time) const
 *     {
 *       AssertThrow((sample_radius > 0.0 && sample_height > 0.0),
 *                   ExcMessage("Non-physical sample dimensions"));
 *       AssertThrow(lambda_2 > 0.0,
 *                   ExcMessage("Non-physical applied axial stretch"));
 * 
 *       const double sqrt_lambda_2     = std::sqrt(lambda_2);
 *       const double inv_sqrt_lambda_2 = 1.0 / sqrt_lambda_2;
 * 
 *       const double alpha_max =
 *         std::atan(std::tan(gamma_12) * sample_height /
 *                   sample_radius); // Small strain approximation
 *       const double A       = sample_radius * alpha_max;
 *       const double w       = 2.0 * numbers::PI * frequency; // in rad /s
 *       const double gamma_t = A * std::sin(w * time);
 *       const double tau_t =
 *         gamma_t /
 *         (sample_radius * sample_height); // Torsion angle per unit length
 *       const double alpha_t = tau_t * lambda_2 * sample_height;
 * 
 *       Tensor<2, 3> F;
 *       F[0][0] = inv_sqrt_lambda_2 * std::cos(alpha_t);
 *       F[0][1] = -inv_sqrt_lambda_2 * std::sin(alpha_t);
 *       F[0][2] = -tau_t * sample_radius * sqrt_lambda_2 * std::sin(alpha_t);
 *       F[1][0] = inv_sqrt_lambda_2 * std::sin(alpha_t);
 *       F[1][1] = inv_sqrt_lambda_2 * std::cos(alpha_t);
 *       F[1][2] = tau_t * sample_radius * sqrt_lambda_2 * std::cos(alpha_t);
 *       F[2][0] = 0.0;
 *       F[2][1] = 0.0;
 *       F[2][2] = lambda_2;
 * 
 *       AssertThrow((F[0][0] > 0) && (F[1][1] > 0) && (F[2][2] > 0),
 *                   ExcMessage("Non-physical deformation gradient component."));
 *       AssertThrow(std::abs(determinant(F) - 1.0) < 1e-6,
 *                   ExcMessage("Volumetric Jacobian is not equal to unity."));
 * 
 *       return F;
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="RheologicalexperimentParallelplaterotationalrheometer"></a> 
 * <h4>Rheological experiment: Parallel plate rotational rheometer</h4>
 * 

 * 
 * This is the function that will drive the numerical experiments.
 * 
 * @code
 *     template <int dim>
 *     void run_rheological_experiment(
 *       const RheologicalExperimentParameters &experimental_parameters,
 *       Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *         &material_hand_calculated,
 *       Coupled_Magnetomechanical_Constitutive_Law_Base<dim>
 *         &               material_assisted_computation,
 *       TimerOutput &     timer,
 *       const std::string filename)
 *     {
 * @endcode
 * 
 * We can take the hand-implemented constitutive law and compare the
 * results that we attain with it to those that we get using AD or SD.
 * In this way, we can verify that they produce identical results (which
 * indicates that either both implementations have a high probability of
 * being correct, or that they're incorrect with identical flaws being
 * present in both). Either way, it is a decent sanity check for the
 * fully self-implemented variants and can certainly be used as a
 * debugging strategy when differences between the results are
 * detected).
 * 
 * @code
 *       const auto check_material_class_results =
 *         [](
 *           const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &to_verify,
 *           const Coupled_Magnetomechanical_Constitutive_Law_Base<dim> &blessed,
 *           const double tol = 1e-6) {
 *           (void)to_verify;
 *           (void)blessed;
 *           (void)tol;
 * 
 *           Assert(std::abs(blessed.get_psi() - to_verify.get_psi()) < tol,
 *                  ExcMessage("No match for psi. Error: " +
 *                             Utilities::to_string(std::abs(
 *                               blessed.get_psi() - to_verify.get_psi()))));
 * 
 *           Assert((blessed.get_B() - to_verify.get_B()).norm() < tol,
 *                  ExcMessage("No match for B. Error: " +
 *                             Utilities::to_string(
 *                               (blessed.get_B() - to_verify.get_B()).norm())));
 *           Assert((blessed.get_S() - to_verify.get_S()).norm() < tol,
 *                  ExcMessage("No match for S. Error: " +
 *                             Utilities::to_string(
 *                               (blessed.get_S() - to_verify.get_S()).norm())));
 * 
 *           Assert((blessed.get_DD() - to_verify.get_DD()).norm() < tol,
 *                  ExcMessage("No match for BB. Error: " +
 *                             Utilities::to_string(
 *                               (blessed.get_DD() - to_verify.get_DD()).norm())));
 *           Assert((blessed.get_PP() - to_verify.get_PP()).norm() < tol,
 *                  ExcMessage("No match for PP. Error: " +
 *                             Utilities::to_string(
 *                               (blessed.get_PP() - to_verify.get_PP()).norm())));
 *           Assert((blessed.get_HH() - to_verify.get_HH()).norm() < tol,
 *                  ExcMessage("No match for HH. Error: " +
 *                             Utilities::to_string(
 *                               (blessed.get_HH() - to_verify.get_HH()).norm())));
 *         };
 * 
 * @endcode
 * 
 * We'll be outputting the constitutive response of the material to file
 * for post-processing, so here we declare a `stream` that will act as
 * a buffer for this output. We'll use a simple CSV format for the
 * outputted results.
 * 
 * @code
 *       std::ostringstream stream;
 *       stream
 *         << "Time;Axial magnetic field strength [A/m];Axial magnetic induction [T];Shear strain [%];Shear stress [Pa]\n";
 * 
 * @endcode
 * 
 * Using the DiscreteTime class, we iterate through each timestep using
 * a fixed time step size.
 * 
 * @code
 *       for (DiscreteTime time(experimental_parameters.start_time(),
 *                              experimental_parameters.end_time() +
 *                                experimental_parameters.delta_t(),
 *                              experimental_parameters.delta_t());
 *            time.is_at_end() == false;
 *            time.advance_time())
 *         {
 *           if (experimental_parameters.print_status(time.get_step_number()))
 *             std::cout << "Timestep = " << time.get_step_number()
 *                       << " @ time = " << time.get_current_time() << "s."
 *                       << std::endl;
 * 
 * @endcode
 * 
 * We fetch and compute the loading to be applied to the material
 * at this time step...
 * 
 * @code
 *           const Tensor<1, dim> H =
 *             experimental_parameters.get_H(time.get_current_time());
 *           const Tensor<2, dim> F =
 *             experimental_parameters.get_F(time.get_current_time());
 *           const SymmetricTensor<2, dim> C =
 *             Physics::Elasticity::Kinematics::C(F);
 * 
 * @endcode
 * 
 * ... then we update the state of the materials...
 * 
 * @code
 *           {
 *             TimerOutput::Scope timer_section(timer, "Hand calculated");
 *             material_hand_calculated.update_internal_data(C, H, time);
 *             material_hand_calculated.update_end_of_timestep();
 *           }
 * 
 *           {
 *             TimerOutput::Scope timer_section(timer, "Assisted computation");
 *             material_assisted_computation.update_internal_data(C, H, time);
 *             material_assisted_computation.update_end_of_timestep();
 *           }
 * 
 * @endcode
 * 
 * ... and test for discrepancies between the two.
 * 
 * @code
 *           check_material_class_results(material_hand_calculated,
 *                                        material_assisted_computation);
 * 
 *           if (experimental_parameters.output_data_to_file)
 *             {
 * @endcode
 * 
 * The next thing that we will do is collect some results to
 * post-process. All quantities are in the "current configuration"
 * (rather than the "reference configuration", in which all
 * quantities computed by the constitutive laws are framed).
 * 
 * @code
 *               const Tensor<1, dim> h =
 *                 Physics::Transformations::Covariant::push_forward(H, F);
 *               const Tensor<1, dim> b =
 *                 Physics::Transformations::Piola::push_forward(
 *                   material_hand_calculated.get_B(), F);
 *               const SymmetricTensor<2, dim> sigma =
 *                 Physics::Transformations::Piola::push_forward(
 *                   material_hand_calculated.get_S(), F);
 *               stream << time.get_current_time() << ";" << h[2] << ";" << b[2]
 *                      << ";" << F[1][2] * 100.0 << ";" << sigma[1][2] << "\n";
 *             }
 *         }
 * 
 * @endcode
 * 
 * Finally, we output the strain-stress and magnetic loading history to
 * file.
 * 
 * @code
 *       if (experimental_parameters.output_data_to_file)
 *         {
 *           std::ofstream output(filename);
 *           output << stream.str();
 *         }
 *     }
 * 
 * @endcode
 * 
 * 
 * <a name="TheCoupledConstitutiveLawsrunfunction"></a> 
 * <h4>The CoupledConstitutiveLaws::run() function</h4>
 * 

 * 
 * The purpose of this driver function is to read in all of the parameters
 * from file and, based off of that, create a representative instance of
 * each constitutive law and invoke the function that conducts a rheological
 * experiment with it.
 * 
 * @code
 *     void run(int argc, char *argv[])
 *     {
 *       using namespace dealii;
 * 
 *       constexpr unsigned int dim = 3;
 * 
 *       const ConstitutiveParameters          constitutive_parameters;
 *       const RheologicalExperimentParameters experimental_parameters;
 * 
 *       std::string parameter_file;
 *       if (argc > 1)
 *         parameter_file = argv[1];
 *       else
 *         parameter_file = "parameters.prm";
 *       ParameterAcceptor::initialize(parameter_file, "used_parameters.prm");
 * 
 * @endcode
 * 
 * We start the actual work by configuring and running the experiment
 * using our rate-independent constitutive law. The automatically
 * differentiable number type is hard-coded here, but with some clever
 * templating it is possible to select which framework to use at run time
 * (e.g., as selected through the parameter file). We'll simultaneously
 * perform the experiments with the counterpary material law that was
 * fully implemented by hand, and check what it computes against our
 * assisted implementation.
 * 
 * @code
 *       {
 *         TimerOutput timer(std::cout,
 *                           TimerOutput::summary,
 *                           TimerOutput::wall_times);
 *         std::cout
 *           << "Coupled magnetoelastic constitutive law using automatic differentiation."
 *           << std::endl;
 * 
 *         constexpr Differentiation::AD::NumberTypes ADTypeCode =
 *           Differentiation::AD::NumberTypes::sacado_dfad_dfad;
 * 
 *         Magnetoelastic_Constitutive_Law<dim> material(constitutive_parameters);
 *         Magnetoelastic_Constitutive_Law_AD<dim, ADTypeCode> material_ad(
 *           constitutive_parameters);
 * 
 *         run_rheological_experiment(experimental_parameters,
 *                                    material,
 *                                    material_ad,
 *                                    timer,
 *                                    experimental_parameters.output_filename_ri);
 * 
 *         std::cout << "... all calculations are correct!" << std::endl;
 *       }
 * 
 * @endcode
 * 
 * Next we do the same for the rate-dependent constitutive law.
 * The highest performance option is selected as default if SymEngine
 * is set up to use the LLVM just-in-time compiler which (in conjunction
 * with some aggressive compilation flags) produces the fastest code
 * evaluation path of all of the available option. As a fall-back, the
 * so called "lambda" optimizer (which only requires a C++11 compliant
 * compiler) will be selected. At the same time, we'll ask the CAS to
 * perform common subexpression elimination to minimize the number of
 * intermediate calculations used during evaluation.
 * We'll record how long it takes to execute the "initialization" step
 * inside the constructor for the SD implementation, as this is where the
 * abovementioned transformations occur.
 * 
 * @code
 *       {
 *         TimerOutput timer(std::cout,
 *                           TimerOutput::summary,
 *                           TimerOutput::wall_times);
 *         std::cout
 *           << "Coupled magneto-viscoelastic constitutive law using symbolic differentiation."
 *           << std::endl;
 * 
 * #ifdef DEAL_II_SYMENGINE_WITH_LLVM
 *         std::cout << "Using LLVM optimizer." << std::endl;
 *         constexpr Differentiation::SD::OptimizerType optimizer_type =
 *           Differentiation::SD::OptimizerType::llvm;
 *         constexpr Differentiation::SD::OptimizationFlags optimization_flags =
 *           Differentiation::SD::OptimizationFlags::optimize_all;
 * #else
 *         std::cout << "Using lambda optimizer." << std::endl;
 *         constexpr Differentiation::SD::OptimizerType optimizer_type =
 *           Differentiation::SD::OptimizerType::lambda;
 *         constexpr Differentiation::SD::OptimizationFlags optimization_flags =
 *           Differentiation::SD::OptimizationFlags::optimize_cse;
 * #endif
 * 
 *         Magnetoviscoelastic_Constitutive_Law<dim> material(
 *           constitutive_parameters);
 * 
 *         timer.enter_subsection("Initialize symbolic CL");
 *         Magnetoviscoelastic_Constitutive_Law_SD<dim> material_sd(
 *           constitutive_parameters, optimizer_type, optimization_flags);
 *         timer.leave_subsection();
 * 
 *         run_rheological_experiment(experimental_parameters,
 *                                    material,
 *                                    material_sd,
 *                                    timer,
 *                                    experimental_parameters.output_filename_rd);
 * 
 *         std::cout << "... all calculations are correct!" << std::endl;
 *       }
 *     }
 * 
 *   } // namespace CoupledConstitutiveLaws
 * 
 * } // namespace Step71
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * The main function only calls the driver functions for the two sets of
 * examples that are to be executed.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   Step71::SimpleExample::run();
 *   Step71::CoupledConstitutiveLaws::run(argc, argv);
 * 
 *   return 0;
 * }
 * @endcode
examples/step-71/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Introductoryexample"></a><h3>Introductory example</h3>


第一个探索性的例子产生了以下输出。经核实，这三种实现方式产生的结果是相同的。

@code
> ./step-71
Simple example using automatic differentiation...
... all calculations are correct!
Simple example using symbolic differentiation.
... all calculations are correct!
@endcode



<a name="Constitutivemodelling"></a><h3>Constitutive modelling</h3>


为了帮助总结虚拟实验本身的结果，下面是一些图表，显示了材料样品内选定位置的剪切应力，与剪切应变的关系。这些图表显示了在三种不同的磁载荷下的应力-应变曲线，以及（机械）载荷曲线的最后一个周期，当速率依赖型材料达到可重复（"稳态"）响应时。这些类型的图表通常被称为[Lissajous plots](https://en.wikipedia.org/wiki/Lissajous_curve)。粘弹性材料的曲线所呈现的椭圆面积提供了某种衡量材料耗散能量多少的方法，其椭圆度表明粘性反应相对于弹性反应的相位变化。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-71.lissajous_plot-me.png" alt="" width="400">
	<p align="center">
        Lissajous plot for the magneto-elastic material.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-71.lissajous_plot-mve.png" alt="" width="400">
	<p align="center">
        Lissajous plot for the magneto-viscoelastic material.
	</p>
    </td>
  </tr>
</table> 

看到磁弹性材料的反应有一个与加载曲线相匹配的卸载曲线并不奇怪--毕竟该材料是非耗散性的。但在这里可以清楚地注意到，随着施加磁场的增加，曲线的梯度是如何增加的。沿着这条曲线任何一点的切线都与瞬时剪切模量有关，由于能量密度函数的定义方式，我们预计剪切模量会随着磁场强度的增加而增加。我们观察到磁-粘弹性材料的行为大致相同。由加载-卸载曲线追踪的椭圆的主轴有一个斜率，随着施加更大的磁载荷而增加。同时，材料耗散的能量也越多。

至于代码输出，这是打印到控制台的与磁弹性材料进行的流变学实验有关的部分的内容。

@code
Coupled magnetoelastic constitutive law using automatic differentiation.
Timestep = 0 @ time = 0s.
Timestep = 125 @ time = 0.314159s.
Timestep = 250 @ time = 0.628318s.
Timestep = 375 @ time = 0.942477s.
...
Timestep = 12250 @ time = 30.7876s.
Timestep = 12375 @ time = 31.1018s.
Timestep = 12500 @ time = 31.4159s.
... all calculations are correct!
@endcode



而这部分输出与用磁涡流材料进行的实验有关。

@code
Coupled magneto-viscoelastic constitutive law using symbolic differentiation.
Using LLVM optimizer.
Timestep = 0 @ time = 0s.
Timestep = 125 @ time = 0.314159s.
Timestep = 250 @ time = 0.628318s.
Timestep = 375 @ time = 0.942477s.
...
Timestep = 12250 @ time = 30.7876s.
Timestep = 12375 @ time = 31.1018s.
Timestep = 12500 @ time = 31.4159s.
... all calculations are correct!
@endcode



计时器的输出也被发射到控制台，因此我们可以比较进行手工计算和辅助计算所需的时间，并对使用AD和SD框架的开销有一些了解。下面是使用AD框架的磁弹性实验的时间，基于Trilinos库的Sacado组件。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |       3.2s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assisted computation            |     12501 |      3.02s |        95% |
| Hand calculated                 |     12501 |    0.0464s |       1.5% |
+---------------------------------+-----------+------------+------------+
@endcode

关于使用自动微分进行的计算（作为提醒，这是使用Sacado库结合动态前向自动微分类型进行的两级微分），我们观察到辅助计算需要大约 $65 \times$ 的时间来计算所需的数量。这看起来确实是一个相当大的开销，但是，正如介绍中提到的，这是否可以接受，完全是主观的，取决于环境的。在对导数进行必要的手工计算、验证其正确性、实现它们以及验证实现的正确性方面，你是否更看重计算机时间而不是人的时间？如果你开发的研究代码只在相对较少的实验中运行，你可能更看重自己的时间。如果你开发的是一个将在万核集群上反复运行数小时的生产代码，你的考虑可能就不同了。在任何情况下，AD方法的一个很好的特点是，当函数和类在标量类型上被模板化时，有 "滴入 "能力。这意味着开始使用它需要付出最小的努力。

相比之下，使用准时制（JIT）编译的符号代数实现的磁涡弹材料的时间表明，在初始化过程中付出一些不可忽视的代价，计算本身的执行效率要高得多。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      1.34s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assisted computation            |     12501 |     0.376s |        28% |
| Hand calculated                 |     12501 |     0.368s |        27% |
| Initialize symbolic CL          |         1 |     0.466s |        35% |
+---------------------------------+-----------+------------+------------+
@endcode

由于初始化阶段很可能只需要在每个线程中执行一次，这个初始的昂贵阶段可以通过重复使用一个 Differentiation::SD::BatchOptimizer 实例来抵消。尽管与磁弹性构成法相比，磁弹性构成法有更多的条款需要计算，但它在执行动能变量和切线的计算方面仍然快了一个数量级。而与使用缓存方案的手工计算变量相比，计算时间几乎相等。因此，尽管使用符号框架需要在如何实现和操作符号表达方面进行范式转变，但它可以提供AD框架所缺乏的良好性能和灵活性。

在数据缓存这一点上，事实上，在用这种材料进行的数值实验中，与使用中间值的实现相比，磁涡流材料实现的数值缓存所增加的成本大约是 $6\times$ ，在`update_internal_data()`中花费的时间增加。下面是删除缓存数据结构时为 "手工计算 "变体提取的时间比较样本输出。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      1.01s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assisted computation            |     12501 |     0.361s |        36% |
| Hand calculated                 |     12501 |    0.0562s |       5.6% |
| Initialize symbolic CL          |         1 |     0.469s |        47% |
+---------------------------------+-----------+------------+------------+
@endcode



通过一些小的调整，我们可以很容易地测试批量优化器的不同优化方案。因此，让我们比较一下与 "LLVM "批处理优化器设置相关的计算费用和其他方案。下面是 "lambda "优化方法的时间报告（保留了CSE的使用）。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      3.87s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assisted computation            |     12501 |      3.12s |        81% |
| Hand calculated                 |     12501 |     0.394s |        10% |
| Initialize symbolic CL          |         1 |     0.209s |       5.4% |
+---------------------------------+-----------+------------+------------+
@endcode

这里的主要观察结果是，与 "LLVM "方法相比，在 "辅助计算 "部分花费的时间要多一个数量级。

最后，我们将测试 "字典 "替换与CSE的结合情况。字典替换只是在本地CAS框架内进行了所有的评估，没有对底层数据结构进行任何转换。在这种情况下，只有使用缓存中间结果的CSE才能提供任何 "加速"。考虑到这一点，下面是这个选择的结果。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |  1.54e+03s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assisted computation            |     12501 |  1.54e+03s |     1e+02% |
| Hand calculated                 |     12501 |     0.563s |         0% |
| Initialize symbolic CL          |         1 |     0.184s |         0% |
+---------------------------------+-----------+------------+------------+
@endcode

不用说，与其他两种方法相比，这些结果花了相当多的时间来产生......字典 "替换方法也许只适用于简单的表达式，或者当调用的数量足够少的时候。

<a name="SowhichframeworkshouldIuse"></a><h1>So, which framework should I use?</h1>


也许你已经相信这些工具有一些优点，并能对你有直接的帮助或用途。现在明显的问题是要使用哪一个。特别是在连续点水平上，你将使用这些框架来计算构成法的导数，我们可以说以下几点。

- 自动区分可能提供了进入辅助区分世界的最简单的切入点。

- 考虑到一个构成框架的足够通用的实现，AD通常可以被用作内在标量类型的替代品，然后可以利用辅助类来计算一阶（以及可能的高阶）导数，只需付出最小的努力。

- 作为对上述观点的限定，"直接替换 "并不意味着你必须对这些数字所通过的算法不持争议态度。有可能在不经意间进行的操作，在进行区分时，会返回一个错误的结果。   所以这绝对是一个人应该注意的事情。   一个具体的例子。当计算一个张量的特征值时，如果该张量是对角线的，那么得到结果的捷径就是直接返回对角线条目（从输入张量中提取的）。就计算特征值本身而言，这是完全正确的，但是不通过算法来计算非对角线张量的特征值会产生意想不到的副作用，即特征值看起来（对AD框架而言）是完全相互脱钩的，它们的交叉敏感度没有被编码在返回的结果中。在进行微分时，导数张量的许多条目将被丢失。为了解决这个问题，我们必须确保使用标准的特征值求解算法，这样返回的特征值对彼此的敏感性就会在结果中得到编码。

- 涉及AD数字类型的计算可能很昂贵。随着微分运算顺序的增加，费用也会增加（有时相当可观）。这可能会被周围操作的计算复杂性所缓解（例如线性求解），但最终还是要看具体问题。

- AD被限制在只需要总导数的情况下。如果一个微分运算需要相对于自变量的偏导，那么使用它是不合适的。

- 每个AD库都有自己的怪癖（说起来很悲哀，但根据作者的经验，是真的），所以可能需要一些试验和错误来找到合适的库和选择AD号来满足你的目的。这些 "怪癖 "的原因往往归结于库背后的整体理念（数据结构、模板元编程的使用等）以及导数计算的数学实现（例如，使用对数函数改变基础的结果操作可能会限制输入值的域--当然，细节都是对用户隐藏的）。   此外，一个库可能比另一个库能更快地计算出所需的结果，所以在这方面进行一些初步探索可能是有益的。

- 符号微分（好吧，一般来说，使用CAS）提供了最灵活的框架，可以进行辅助计算。

- SD框架可以做AD框架能做的所有事情，还有一个好处是可以对何时进行某些操纵和操作进行低层次控制。

- 加速表达式的评估是可能的，与一些手工实现相比，有可能导致SD框架接近原生的性能（当然，这种比较取决于整个程序设计），但代价是初始优化调用。

- 巧妙地使用 Differentiation::SD::BatchOptimizer 可以将优化依赖表达式的昂贵调用的费用降到最低。   对 Differentiation::SD::BatchOptimizer 进行序列化的可能性，往往（但不总是）这种昂贵的调用可以做一次，然后在以后的模拟中重复使用。

- 例如，如果两个或更多的材料法只因其材料参数而不同，那么只要这些材料参数被认为是象征性的，就可以在它们之间共享一个批次优化器。这意味着你可以 "区分一次，在许多情况下评估"。

- SD框架可以部分地被用作标量类型的 "直接替换"，但人们（至少）必须在它周围增加一些框架来执行值替换步骤，将符号类型转换为它们的数字对应物。

- 在一些专门的算法中可能无法使用SD数字。   例如，如果一个算法的退出点或代码分支是基于（符号）输入参数应该采取的一些具体的数值，那么显然这是不可能的。我们要么重新实现专门针对SD数字类型的算法（有点不方便，但经常是可能的，因为 Differentiation::SD::Expression 类支持条件反射），要么必须使用创造性的手段来解决这个具体问题（例如，引入一个符号表达式来表示这个算法返回的结果，如果在要使用它的环境中是有意义的，也许可以将它声明为一个[符号函数]（https://dealii.org/developer/doxygen/deal.II/namespaceDifferentiation_1_1SD.html#a876041f6048705c7a8ad0855cdb1bd7a）。这以后可以用它的数值来替代，如果宣布为符号函数，那么它的递延导数也可以作为替代的结果纳入计算中。)

- 使用SD的最大缺点是，使用它需要一个范式的转变，人们必须以不同的方式来构建大多数问题，以便充分利用它。仔细考虑如何使用和重用数据结构也是让它有效工作的关键）。这可能意味着，人们需要对它进行一番玩耍，并建立起对典型操作顺序的理解，以及每一步在操作基础数据方面的具体作用。如果人们有时间和意愿这样做，那么使用这个工具的好处可能是巨大的。

<a name="Possibilitiesforextension"></a><h1>Possibilities for extension</h1>


有几个合乎逻辑的方法可以扩展这个计划。

- 也许最明显的扩展是实施和测试其他构成模型。   这可能仍然属于磁-机械耦合问题的范畴，也许可以考虑替代能量函数的 "Neo-Hookean "型弹性部分，改变耗散能量的构成法则（及其相关的演化法则），或者包括磁滞效应或这些材料试图模拟的复合聚合物的损坏模型。

- 当然，所实现的模型可以被修改或完全替换为专注于物理学其他方面的模型，如电活性聚合物、生物力学材料、弹塑性介质等。

- 对粘弹性演化法实施不同的时间微调方案。

- 与其直接从能量密度函数推导出一切，不如使用 Differentiation::AD::VectorFunction 直接线性化动力学量。   这将意味着只需要一个可微分的自动微分的数字类型，并且肯定会大大改善性能。   这种方法也为耗散材料提供了机会，比如这里考虑的磁涡弹材料，可以与AD结合起来实现。这是因为线性化调用了因变量相对于场变量的总导数，这正是AD框架所能提供的。

- 调查使用其他可自动微分的数字类型和框架（如ADOL-C）。由于每个AD库都有自己的实现，选择使用哪个库可能会导致性能的提高，在最不幸的情况下，计算也会更加稳定。至少可以说，对于deal.II支持的AD库，结果的准确性应该基本不受这个决定的影响。

- 在有限元模拟中嵌入这些构成法则中的一个。

如果不费吹灰之力，人们可以考虑重新编写非线性问题求解器，比如在步骤15中实现的使用AD或SD方法来计算牛顿矩阵的求解器。事实上，这在第72步中已经完成。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-71.cc"
*/
