/**
@page step_78 The step-78 tutorial program
This tutorial depends on step-26.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Particularitiesoftheequationsystem">Particularities of the equation system</a>
        <li><a href="#Schemeforthenumericalsolution">Scheme for the numerical solution</a>
        <li><a href="#TestCase">Test Case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#SolutionClass">Solution Class</a>
        <li><a href="#EquationData">Equation Data</a>
        <li><a href="#ThecodeBlackScholescodeClass">The <code>BlackScholes</code> Class</a>
        <li><a href="#ThecodeBlackScholescodeImplementation">The <code>BlackScholes</code> Implementation</a>
      <ul>
        <li><a href="#codeBlackScholessetup_systemcode"><code>BlackScholes::setup_system</code></a>
        <li><a href="#codeBlackScholessolve_time_stepcode"><code>BlackScholes::solve_time_step</code></a>
        <li><a href="#codeBlackScholesadd_results_for_outputcode"><code>BlackScholes::add_results_for_output</code></a>
        <li><a href="#codeBlackScholesrefine_gridcode"><code>BlackScholes::refine_grid</code></a>
        <li><a href="#codeBlackScholesprocess_solutioncode"><code>BlackScholes::process_solution</code></a>
        <li><a href="#codeBlackScholeswrite_convergence_tablecode"><code>BlackScholes::write_convergence_table</code> </a>
        <li><a href="#codeBlackScholesruncode"><code>BlackScholes::run</code></a>
      </ul>
        <li><a href="#ThecodemaincodeFunction">The <code>main</code> Function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-78/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


布莱克-斯科尔斯方程是一个偏微分方程，它有点不符合普通的方案。它描述了一个 "欧洲看涨 "股票期权的公平价格是多少。不用说得太详细，股票 "期权 "是一个人可以从银行购买的合同，它允许我，但不要求我，在未来某个固定的时间 $T$ 以固定的价格 $K$ 购买特定股票。那么，作为这种期权的买方要回答的问题是："我认为这样的合约值多少钱？"，或者作为卖方，"我需要为这个合约收取多少钱？"，这既是合约在时间 $t<T$ 前的函数，也是股票价格 $S$ 的函数。费舍尔-布莱克和迈伦-斯科尔斯在假设股票价格表现出随机的价格波动，具有一定的 "波动性"，再加上背景的指数价格上涨（可以认为是通货膨胀率，随着时间的推移，所有货币都会贬值）的情况下，得出了这种期权的公平价格 $V(S,t)$ 的偏微分方程。由于他们的工作，布莱克和斯科尔斯在1997年获得了诺贝尔经济科学奖，这使得这是第一个处理有人获得诺贝尔奖的问题的教程程序  @cite black1973pricing  。

该公式如下。

@f{align*}{
    &\frac{\partial V}{\partial t} + \frac{\sigma^2S^2}{2} \
    \frac{\partial^2 V}{\partial S^2} + \
    rS\frac{\partial V}{\partial S} - rV = 0, \
    \quad\quad &&\forall S\in \Omega, t \in (0,T)
    \\
    &V(0,t) = 0, \
    &&\forall t \in (0,T)
    \\
    &V(S,t) \rightarrow S, \
    && \text{as } S \rightarrow \infty, \forall t \in (0,T)
    \\
    &V(S,T) = \max(S-K,0) \
    &&\forall S \in \Omega


@f}

其中

@f{align*}{
    V(S,t): && \text{Value of call option at time t and asset price S} \\
    \sigma: && \text{Volatility of the underlying asset} \\
    r: && \text{Risk free interest rate} \\
    K : && \text{Strike price for purchasing asset}


@f}



我们应该这样解释这个方程，它是一个时间依赖的偏微分方程，一个 "空间 "变量  $S$  是股票的价格， $V(S,t)$  是时间  $t$  的期权价格，如果当时的股票价格是  $S$  。

<a name="Particularitiesoftheequationsystem"></a><h3>Particularities of the equation system</h3>


在继续讨论其数值解法之前，这个方程中存在一些值得讨论的怪异现象。首先，"空间 "域 $\Omega\subset\mathbb{R}$ 是无界的，因此 $S$ 的值也可以是无界的。这是因为股票价格可能有一个实际的上限，但没有一个概念上的上限。那么边界条件 $V(S,t)\rightarrow S$ 作为 $S\rightarrow \infty$ 可以解释如下。如果股票价格（今天或 $t=T$ 时）是 $S\gg K$ ，那么允许我以价格 $K$ 购买股票的期权的价值是多少？人们期望它是 $V\approx S-K$ 加上一些通货膨胀的调整，或者，如果我们真的真正考虑到 $S$ 的巨大价值，我们可以忽略 $K$ ，得出无限边界的边界值应该是上面所说的 $V\rightarrow S$ 形式。

在实践中，对于我们使用有限元方法来解决这个问题，我们需要对 $\Omega$ 进行约束。由于这个方程描述的是价格，而谈论价格为负值是没有意义的，我们将设定 $\Omega$ 的下限为0。然后，对于上限，我们将选择一个非常大的数字，一个 $S$ 不太可能达到的数字。我们将称其为 $S_\text{max}$ 。所以， $\Omega=[0,S_\text{max}]$  。

第二，在截断域之后，我们需要问在这个现在的有限边界上我们应该摆出什么边界值。为了解决这个问题，我们使用 "看跌-看涨 "平价  @cite stoll1969relationship  。一个 "拉动期权 "是指我们被允许，但不是必须，在未来的某个时间以价格  $K$  向某人出售*股票  $T$  。这说明

@f{align*}{
    V(S,t)+Ke^{-r(T-t)}=P(S,t)+S


@f}

其中 $V(S,t)$ 是看涨期权的价值，而 $P(S,t)$ 是看跌期权的价值。由于我们预计 $P(S,t) \rightarrow 0$ 为 $S \rightarrow \infty$ ，这说明

@f{align*}{
    V(S,t) \rightarrow S-Ke^{-r(T-t)},


@f}

而我们可以用这个作为我们的有限点的合理边界条件  $S_\text{max}$  。

Block-Scholes 方程的第二个复杂之处在于，我们得到的是一个最终条件，而不是初始条件。这是因为我们知道期权在 $t=T$ 时的价值：如果 $T$ 时的股票价格是 $S<K$ ，那么我们就没有动力使用我们的期权买入价格 $K$ ，因为我们可以在公开市场上以更低的价格买入该股票。所以 $V(S,T)=0$ 为 $S<K$  。另一方面，如果在 $T$ 时我们有 $S>K$ ，那么我们可以通过期权以 $K$ 的价格买入股票，并立即在市场上以 $S$ 的价格再次卖出，给我带来 $S-K$ 的利润。换句话说， $V(S,T)=S-K$ 换取 $S>K$  。因此，我们只知道在*结束时间*的 $V$ 的值，而不知道初始时间--事实上，找出当前时间的公平价格（传统上被认为是 $t=0$ ）是解决这些方程的关键所在。

这意味着这不是一个在时间上向前推进的方程，实际上是在时间上*向后推进的。因此，通过改变变量 $\tau=T-t$ 来反向解决这个问题是有意义的，现在 $\tau$ 表示打击时间 $T$ 之前的时间。

有了这一切，我们最终得到了以下问题。

@f{align*}{
    &-\frac{\partial V}{\partial \tau} + \frac{\sigma^2S^2}{2} \
    \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV=0\
    , \quad\quad &&\forall S\in [0,S_\text{max}], \tau \in [0,T]
    \\
    &V(0,\tau) = 0, \
    &&\forall \tau \in [0,T]
    \\
    &V(S_\text{max},\tau)=S_\text{max}-Ke^{-r\tau}, \
    &&\forall \tau \in [0,T]
    \\
    &V(S,0) = \max(S-K,0) \
    &&\forall S \in [0,S_\text{max}]


@f}



从概念上讲，这是一个变量 $V$ 的平流-扩散-反应问题：有一个二阶导数的扩散项，一个一阶导数的平流项，以及一个零阶反应项。我们可以预期这个问题在实践中会有一些宽容，因为对于现实中的系数值，它是扩散主导的。但是，由于问题中的平流项，我们将不得不小心地进行网格细化和时间步长的选择。还有一个问题是，扩散项是以非保守的形式写的，因此按部分积分并不明显。这将在下一节中讨论。

<a name="Schemeforthenumericalsolution"></a><h3>Scheme for the numerical solution</h3>


我们将使用IMEX方法解决这个问题。特别是，我们首先用theta方法进行时间离散，随后将为平流和扩散项选择不同的theta值。让 $V^n(S)$ 近似于 $V(S,\tau_n)$  。

@f{align*}{
    0=&-\frac{V^n(S)-V^{n-1}(S)}{k_n} \\
    &+\frac{\sigma^2S^2}{2}\left[(1-\theta)\frac{d^2V^{n-1}(S)}{dS^2} + \
    \theta \frac{d^2V^{n}(S)}{dS^2}\right] \\
    &+rS\left[(1-\theta)\frac{dV^{n-1}(S)}{dS} + \
    \theta\frac{dV^{n}(S)}{dS}\right]  \\
    &-r\left[(1-\theta)V^{n-1}(S) + \theta V^n(S)\right]


@f}

这里， $k_n=\tau_n-\tau_{n-1}$ 是时间步长。鉴于这种时间离散化，我们可以通过与测试函数相乘，然后通过部分积分来进行空间离散化。因为其中有一些有趣的细节，由于这个方程中的平流和非平流项，这个过程将被详细解释。

因此，我们首先用测试函数相乘， $\{\phi_i(S)\}_{i\in\mathbb{N}}$  。

@f{align*}{
    0=&-\int_0^{S_\text{max}}\phi_i(S)\left[V^n(S)-V^{n-1}(S)\right]dS \\
    &+k_n\int_0^{S_\text{max}}\phi_i(S)\left[\frac{\sigma^2S^2}{2} \
    \left[(1-\theta)\frac{d^2V^{n-1}(S)}{dS^2} + \
     \theta \frac{d^2V^{n}(S)}{dS^2}\right]\right]dS \\
    &+k_n\int_0^{S_\text{max}}\phi_i(S)\left[rS\left[(1-\theta)
     \frac{dV^{n-1}(S)}{dS}\
     + \theta\frac{dV^{n}(S)}{dS}\right]\right]dS  \\
    &-k_n\int_0^{S_\text{max}}\phi_i(S)\left[r\left[(1-\theta)V^{n-1}(S)\
     + \theta V^n(S)\right]\right]dS


@f}




像往常一样，（1）变成 $-\textbf{M}V^n+\textbf{M}V^{n-1}$ ，（4）变成 $k_n\left[-r(1-\theta)\textbf{M}V^{n-1} - \theta r\textbf{M}V^n\right]$ ，其中 $\textbf{M}_{i,j}=\left(\phi_i(S),\phi_j(S)\right)$ ，我们不仅用 $V$ 表示函数 $V(S)$ ，而且用离散化后的节点值向量表示。

有趣的部分来自于（2）和（3）。


对于（2），我们有。

@f{align*}{
    &k_n\int_0^{S_\text{max}}\phi_i(S)\left[\frac{\sigma^2S^2}{2} \
     \left[(1-\theta)\frac{d^2V^{n-1}(S)}{dS^2} + \
     \theta \frac{d^2V^{n}(S)}{dS^2}\right]\right]dS \\
    &=k_n(1-\theta)\int_0^{S_\text{max}}\phi_i(S)\frac{\sigma^2S^2}{2} \
     \frac{d^2V^{n-1}(S)}{dS^2} \
    +k_n\theta\int_0^{S_\text{max}}\phi_i(S)\frac{\sigma^2S^2}{2} \
     \frac{d^2V^{n}(S)}{dS^2}


@f}



这里有两个积分，或多或少都是一样的，区别在于积分前面的系数略有不同，以及V的时间步骤不同。因此，考虑一般的积分，我们将用部分积分的方法来解决这个问题。

@f{align*}{
    &\int_{0}^{S_\text{max}} \phi_i(S)\frac{\sigma^2S^2}{2}
        \frac{d^2V^n(S)}{dS^2}dS \\
    &= \phi_i(S)\frac{1}{2}\sigma^2S^2\frac{dV^n(S)}{dS}\Bigg|_0^{S_{max}} - \
    \int_0^{S_\text{max}} \phi_i(S)\sigma^2S\frac{dV^n(S)}{dS}dS - \
    \int_0^{S_\text{max}} \frac{d\phi_i(S)}{dS}\frac{1}{2}\sigma^2S^2 \
    \frac{dV^n(S)}{dS}dS \\
    &= -\int_0^{S_\text{max}} \phi_i(S)\sigma^2S\frac{dV^n(S)}{dS}dS - \
    \int_0^{S_\text{max}} \frac{d\phi_i(S)}{dS}\frac{1}{2}\sigma^2S^2 \
    \frac{dV^n(S)}{dS}dS \\
    &= -\int_0^{S_\text{max}} \phi_i(S)\sigma^2S \sum_j V_j^n
        \frac{d\phi_j(S)}{dS}dS \


    -\int_0^{S_\text{max}} \frac{d\phi_i(S)}{dS}\frac{1}{2} \
    \sigma^2S^2  \sum_k V_k^n \frac{d\phi_k(S)}{dS}dS \\
    &= -\sum_j \sigma^2 \int_0^{S_\text{max}} \phi_i(S)S
        \frac{d\phi_j(S)}{dS}dS V_j^n\


    - \sum_k \frac{1}{2}\sigma^2 \int_0^{S_\text{max}} \frac{d\phi_i(S)}{dS}S^2\
    \frac{d\phi_k}{dS}dS V_k^n \\
    &= -\sum_j \sigma^2 \left(\phi_i(S)S, \frac{d\phi_j(S)}{dS}\right) V_j^n \


    - \sum_k \frac{1}{2}\sigma^2 \left(\frac{d\phi_i(S)}{dS}S^2,\
    \frac{d\phi_k(S)}{dS}\right) V_k^n \\
    &= -\sigma^2\textbf{B}V^n - \frac{1}{2}\sigma^2\textbf{D}V^n, \quad\quad \
    \textbf{B}_{i,j} = \left(\phi_i(S)S, \frac{d\phi_j(S)}{dS}\right),\
    \textbf{D}_{i,j} = \left(\frac{d\phi_i(S)}{dS}S^2,\frac{d\phi_j(S)}{dS}\right)


@f}



因此，在加入常数并将 $V^n$ 换成 $V^{n-1}$ （如适用）后，我们对（2）得出如下结论。

@f{align*}{
    &k_n\int_0^{S_\text{max}}\phi_i(S)\left[\frac{\sigma^2S^2}{2}
        \left[(1-\theta)\
    \frac{d^2V^{n-1}(S)}{dS^2} + \
    \theta \frac{d^2V^{n}(S)}{dS^2}\right]\right]dS \\
    &= k_n\left[-(1-\theta)\sigma^2\textbf{B}V^{n-1}\


     -(1-\theta)\frac{1}{2}\sigma^2\textbf{D}V^{n-1} \


    -\theta\sigma^2\textbf{B}V^{n}


     -\theta\frac{1}{2}\sigma^2\textbf{D}V^{n}\right]


@f}

但是，由于矩阵 $\textbf{B}$ 涉及一个平流项，我们将在那里选择 $\theta=0$ --换句话说，我们使用显式欧拉方法来处理平流问题。相反，由于矩阵 $\textbf{D}$ 涉及扩散项，我们将在那里选择 $\theta=1/2$ --即我们用二阶Crank-Nicolson方法处理扩散。

因此，我们得出以下结论。

@f{align*}{
    k_n\left[-\frac{1}{4}\sigma^2\textbf{D}V^{n-1} \


    -\frac{1}{4}\sigma^2\textbf{D}V^n \


    - \sigma^2\textbf{B}V^{n-1}\right]


@f}



现在，要处理（3）。为此，我们将再次通过考虑上述的一般情况来进行。

@f{align*}{
    &\int_{0}^{S_\text{max}} \phi_i(S)rS\frac{dV^n}{dS}dS \\
    &= \phi_i(S)rSV^n\Bigg|_0^{S_\text{max}} - \int_0^{S_\text{max}}
        \left[r\phi_i(S) \
    + r\frac{d\phi_i(S)}{dS}S \right]V^ndS \\
    &= -\int_0^{S_\text{max}} r\phi_i(S)V^ndS - \
    \int_0^{S_\text{max}} r\frac{d\phi_i(S)}{dS}SV^ndS \\
    &= -\int_0^{S_\text{max}} r\phi_i(S) \sum_j V_j^n\phi_j(S)dS \


    -\int_0^{S_\text{max}} rS\frac{d\phi_i(S)}{dS} \sum_k V_k\phi_k(S)dS \\
    &= -\sum_j r\left(\phi_i(S), \phi_j(S)\right) V_j^n -\
     \sum_k r\left(S\frac{d\phi_i(S)}{dS}, \phi_k(S)\right)V_k^n \\
    &= -r\textbf{M}V^n -r\textbf{A}V^n, \quad\quad\
    \textbf{M}_{i,j} = \left(\phi_i(S), \phi_j(S)\right),\
    \textbf{A}_{i,j} = \left(S\frac{d\phi_i(S)}{dS}, \phi_j(S)\right)


@f}



因此，在加入常数并将 $V^n$ 换成 $V^{n-1}$ （如适用）后，我们对（3）得出如下结论。

@f{align*}{
    &k_n\int_0^{S_\text{max}}\phi_i(S)\left[rS\left[(1-\theta)
        \frac{dV^{n-1}(S)}{dS} +\
     \theta\frac{dV^{n}(S)}{dS}\right]\right]dS \\
    &= k_n\left[-(1-\theta)r\textbf{M}V^{n-1} -(1-\theta)r\textbf{A}V^{n-1}\


    -\theta r\textbf{M}V^n -\theta r\textbf{A}V^n\right]


@f}

就像以前一样，我们将用 $\theta=0$ 来表示矩阵 $\textbf{A}$ ，用 $\theta=\frac{1}{2}$ 表示矩阵 $\textbf{M}$  。因此，我们对（3）得出以下结果。

@f{align*}{
    k_n\left[-\frac{1}{2}r\textbf{M}V^{n-1} - \frac{1}{2}r\textbf{M}V^n \


    -r\textbf{A}V^{n-1}\right]


@f}



现在，把所有的东西放在一起，我们得到以下布莱克-斯科尔斯方程的离散形式。

@f{align*}{
    0&= \\
    &-\textbf{M}V^n+\textbf{M}V^{n-1} \\
    & +k_n\left[-\frac{1}{4}\sigma^2\textbf{D}V^{n-1} \


    -\frac{1}{4}\sigma^2\textbf{D}V^n \


    - \sigma^2\textbf{B}V^n \


     -\frac{1}{2}r\textbf{M}V^{n-1} - \frac{1}{2}r\textbf{M}V^n \


    -r\textbf{A}V^n \


     -r\frac{1}{2}\textbf{M}V^{n-1} - \frac{1}{2} r\textbf{M}V^n\right] \\
    &= -\textbf{M}V^n + \textbf{M}V^{n-1} +\
    k_n\left[- \frac{1}{4}\sigma^2\textbf{D}V^{n-1} -\
    \frac{1}{4}\sigma^2\textbf{D}V^n - r\textbf{M}V^{n-1} -\
    r\textbf{M}V^n  - \sigma^2\textbf{B}V^{n-1} - r\textbf{A}V^{n-1}\right]


@f}

因此，我们总共有。

@f{equation}{
    0 = \textbf{M}V^n - \textbf{M}V^{n-1} +\
    k_n\left[ \frac{1}{4}\sigma^2\textbf{D}V^{n-1} +\
    \frac{1}{4}\sigma^2\textbf{D}V^n + r\textbf{M}V^{n-1} + r\textbf{M}V^n  +\
    \sigma^2\textbf{B}V^{n-1} + r\textbf{A}V^{n-1}\right]\tag{*}


@f}



像往常一样，我们可以把未知量写在左边，把已知量写在右边。这就导致了在每个时间步骤中必须解决的以下线性系统。

@f{align*}{
    \left[\textbf{M}+\frac{1}{4}k_n\sigma^2\textbf{D}+k_nr\textbf{M}\right]V^n\
     =\
    \left[-\frac{1}{4}k_n\sigma^2\textbf{D}-\
    k_nr\textbf{M}+k_n\sigma^2\textbf{B}-\
    k_nr\textbf{A}+\textbf{M}\right]V^{n-1}


@f}









<a name="TestCase"></a><h3>Test Case</h3> 对于这个程序，我们将使用制造解决方案的方法（MMS）来测试它是否正确工作。这意味着，我们将选择我们的解决方案是类似于步骤7的某个函数。对于我们的案例，我们将使用。


@f{align*}{
    V(S,\tau) = -\tau^2 - S^2 + 6\tag{**}


@f}

这意味着，利用我们的PDE，我们得出了以下问题。

@f{align*}{
    &-\frac{\partial V}{\partial \tau} +\
    \frac{\sigma^2S^2}{2}\frac{\partial^2 V}{\partial S^2} +\
    rS\frac{\partial V}{\partial S} - rV = f(S,\tau) \\
    &V(0,\tau) = -\tau^2 + 6 \\
    &V(S_\text{max}, \tau) = -S_\text{max}^2 - \tau^2 + 6 \\
    &V(S, 0) = -S^2 + 6


@f}

其中， $f(S,\tau) = 2\tau - \sigma^2S^2 - 2rS^2 - r(-\tau^2 - S^2 + 6)$  。这个设置现在有方程本身和 $S=0$ 处的边界条件的右手边，这是我们以前没有的，同时还有不符合实际情况的 "最终 "条件（或者，用 $\tau$ -时间 "初始条件"）。我们将以这样的方式在代码中实现这一点，以便于交换 -- 上述变化的引入只是为了能够使用制造的解决方案。

如果程序工作正常，那么它应该产生（**）作为解决方案。这确实意味着我们需要在一定程度上修改我们的变异形式，以考虑到非零的右手边。

首先，我们定义如下。

@f{align*}{
    F^n_i = \left(\phi_i(S), f^n(S)\right), && \text{where } f^n(S) =\
     f(S,\tau_n)


@f}

因此，我们得出了新的方程式。

@f{align*}{
    \left[\textbf{M}+\frac{1}{4}k_n\sigma^2\textbf{D}+k_nr\textbf{M}\right]V^n\
     =\
     \left[-\frac{1}{4}k_n\sigma^2\textbf{D}-\
     k_nr\textbf{M}+k_n\sigma^2\textbf{B}-\
     k_nr\textbf{A}+\textbf{M}\right]V^{n-1} -\
      k_n\left[\frac{1}{2}F^{n-1}+\frac{1}{2}F^n\right]


@f}



然后我们按上述方法解决这个方程。


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
 * The program starts with the usual include files, all of which you should have
 * seen before by now:
 * 
 * @code
 * #include <deal.II/base/convergence_table.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/logstream.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/utilities.h>
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_accessor.h>
 * #include <deal.II/dofs/dof_tools.h>
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_values.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/tria_accessor.h>
 * #include <deal.II/grid/tria_iterator.h>
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/precondition.h>
 * #include <deal.II/lac/sparse_matrix.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/vector.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/data_out_stack.h>
 * #include <deal.II/numerics/error_estimator.h>
 * #include <deal.II/numerics/matrix_tools.h>
 * #include <deal.II/numerics/solution_transfer.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * Then the usual placing of all content of this program into a namespace and
 * the importation of the deal.II namespace into the one we will work in. We
 * also define an identifier to allow for the MMS code to be run when
 * <code>MMS</code> is defined. Otherwise, the program solves the original
 * problem:
 * 
 * @code
 * namespace BlackScholesSolver
 * {
 *   using namespace dealii;
 * 
 * #define MMS
 * 
 * @endcode
 * 
 * 
 * <a name="SolutionClass"></a> 
 * <h3>Solution Class</h3>
 * 

 * 
 * This section creates a class for the known solution when testing using the
 * MMS. Here we are using $v(\tau,S) = -\tau^2 -S^2 + 6$ for the solution. We
 * need to include the solution equation and the gradient for the H1 seminorm
 * calculation.
 * 
 * @code
 *   template <int dim>
 *   class Solution : public Function<dim>
 *   {
 *   public:
 *     Solution(const double maturity_time);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     virtual Tensor<1, dim>
 *     gradient(const Point<dim> & p,
 *              const unsigned int component = 0) const override;
 * 
 *   private:
 *     const double maturity_time;
 *   };
 * 
 * 
 *   template <int dim>
 *   Solution<dim>::Solution(const double maturity_time)
 *     : maturity_time(maturity_time)
 *   {
 *     Assert(dim == 1, ExcNotImplemented());
 *   }
 * 
 * 
 *   template <int dim>
 *   double Solution<dim>::value(const Point<dim> & p,
 *                               const unsigned int component) const
 *   {
 *     return -Utilities::fixed_power<2, double>(p(component)) -
 *            Utilities::fixed_power<2, double>(this->get_time()) + 6;
 *   }
 * 
 * 
 *   template <int dim>
 *   Tensor<1, dim> Solution<dim>::gradient(const Point<dim> & p,
 *                                          const unsigned int component) const
 *   {
 *     return Point<dim>(-2 * p(component));
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="EquationData"></a> 
 * <h3>Equation Data</h3>
 * 

 * 
 * In the following classes and functions, we implement the right hand side
 * and boundary values that define this problem and for which we need function
 * objects. The right hand side is chosen as discussed at the end of the
 * introduction.
 *   

 * 
 * First, we handle the initial condition.
 * 
 * @code
 *   template <int dim>
 *   class InitialConditions : public Function<dim>
 *   {
 *   public:
 *     InitialConditions(const double strike_price);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 *     const double strike_price;
 *   };
 * 
 * 
 *   template <int dim>
 *   InitialConditions<dim>::InitialConditions(const double strike_price)
 *     : strike_price(strike_price)
 *   {}
 * 
 * 
 *   template <int dim>
 *   double InitialConditions<dim>::value(const Point<dim> & p,
 *                                        const unsigned int component) const
 *   {
 * #ifdef MMS
 *     return -Utilities::fixed_power<2, double>(p(component)) + 6;
 * #else
 *     return std::max(p(component) - strike_price, 0.);
 * #endif
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Next, we handle the left boundary condition.
 * 
 * @code
 *   template <int dim>
 *   class LeftBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 *   };
 * 
 * 
 *   template <int dim>
 *   double LeftBoundaryValues<dim>::value(const Point<dim> &,
 *                                         const unsigned int /*component*/) const
 *   {
 * #ifdef MMS
 *     return -Utilities::fixed_power<2, double>(this->get_time()) + 6;
 * #else
 *     return 0.;
 * #endif
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Then, we handle the right boundary condition.
 * 
 * @code
 *   template <int dim>
 *   class RightBoundaryValues : public Function<dim>
 *   {
 *   public:
 *     RightBoundaryValues(const double strike_price, const double interest_rate);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 *     const double strike_price;
 *     const double interest_rate;
 *   };
 * 
 * 
 *   template <int dim>
 *   RightBoundaryValues<dim>::RightBoundaryValues(const double strike_price,
 *                                                 const double interest_rate)
 *     : strike_price(strike_price)
 *     , interest_rate(interest_rate)
 *   {}
 * 
 * 
 *   template <int dim>
 *   double RightBoundaryValues<dim>::value(const Point<dim> & p,
 *                                          const unsigned int component) const
 *   {
 * #ifdef MMS
 *     return -Utilities::fixed_power<2, double>(p(component)) -
 *            Utilities::fixed_power<2, double>(this->get_time()) + 6;
 * #else
 *     return (p(component) - strike_price) *
 *            exp((-interest_rate) * (this->get_time()));
 * #endif
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * Finally, we handle the right hand side.
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     RightHandSide(const double asset_volatility, const double interest_rate);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 *     const double asset_volatility;
 *     const double interest_rate;
 *   };
 * 
 * 
 *   template <int dim>
 *   RightHandSide<dim>::RightHandSide(const double asset_volatility,
 *                                     const double interest_rate)
 *     : asset_volatility(asset_volatility)
 *     , interest_rate(interest_rate)
 *   {}
 * 
 * 
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 * #ifdef MMS
 *     return 2 * (this->get_time()) -
 *            Utilities::fixed_power<2, double>(asset_volatility * p(component)) -
 *            2 * interest_rate * Utilities::fixed_power<2, double>(p(component)) -
 *            interest_rate *
 *              (-Utilities::fixed_power<2, double>(p(component)) -
 *               Utilities::fixed_power<2, double>(this->get_time()) + 6);
 * #else
 *     (void)p;
 *     (void)component;
 *     return 0.0;
 * #endif
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBlackScholescodeClass"></a> 
 * <h3>The <code>BlackScholes</code> Class</h3>
 * 

 * 
 * The next piece is the declaration of the main class of this program. This
 * is very similar to the Step-26 tutorial, with some modifications. New
 * matrices had to be added to calculate the A and B matrices, as well as the
 * $V_{diff}$ vector mentioned in the introduction. We also define the
 * parameters used in the problem.
 *   

 * 
 * - <code>maximum_stock_price</code>: The imposed upper bound on the spatial
 * domain. This is the maximum allowed stock price.
 * - <code>maturity_time</code>: The upper bound on the time domain. This is
 * when the option expires.\n
 * - <code>asset_volatility</code>: The volatility of the stock price.\n
 * - <code>interest_rate</code>: The risk free interest rate.\n
 * - <code>strike_price</code>: The agreed upon price that the buyer will
 * have the option of purchasing  the stocks at the expiration time.
 *   

 * 
 * Some slight differences between this program and step-26 are the creation
 * of the <code>a_matrix</code> and the <code>b_matrix</code>, which is
 * described in the introduction. We then also need to store the current time,
 * the size of the time step, and the number of the current time step.
 * Next, we will store the output into a <code>DataOutStack</code>
 * variable because we will be layering the solution at each time on top of
 * one another to create the solution manifold. Then, we have a variable that
 * stores the current cycle and number of cycles that we will run when
 * calculating the solution. The cycle is one full solution calculation given
 * a mesh. We refine the mesh once in between each cycle to exhibit the
 * convergence properties of our program. Finally, we store the convergence
 * data into a convergence table.
 *   

 * 
 * As far as member functions are concerned, we have a function that
 * calculates the convergence information for each cycle, called
 * <code>process_solution</code>. This is just like what is done in step-7.
 * 
 * @code
 *   template <int dim>
 *   class BlackScholes
 *   {
 *   public:
 *     BlackScholes();
 * 
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void solve_time_step();
 *     void refine_grid();
 *     void process_solution();
 *     void add_results_for_output();
 *     void write_convergence_table();
 * 
 *     const double maximum_stock_price;
 *     const double maturity_time;
 *     const double asset_volatility;
 *     const double interest_rate;
 *     const double strike_price;
 * 
 *     Triangulation<dim> triangulation;
 *     FE_Q<dim>          fe;
 *     DoFHandler<dim>    dof_handler;
 * 
 *     AffineConstraints<double> constraints;
 * 
 *     SparsityPattern      sparsity_pattern;
 *     SparseMatrix<double> mass_matrix;
 *     SparseMatrix<double> laplace_matrix;
 *     SparseMatrix<double> a_matrix;
 *     SparseMatrix<double> b_matrix;
 *     SparseMatrix<double> system_matrix;
 * 
 *     Vector<double> solution;
 *     Vector<double> system_rhs;
 * 
 *     double time;
 *     double time_step;
 * 
 *     const double       theta;
 *     const unsigned int n_cycles;
 *     const unsigned int n_time_steps;
 * 
 *     DataOutStack<dim>        data_out_stack;
 *     std::vector<std::string> solution_names;
 * 
 *     ConvergenceTable convergence_table;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBlackScholescodeImplementation"></a> 
 * <h3>The <code>BlackScholes</code> Implementation</h3>
 * 

 * 
 * Now, we get to the implementation of the main class. We will set the values
 * for the various parameters used in the problem. These were chosen because
 * they are fairly normal values for these parameters. Although the stock
 * price has no upper bound in reality (it is in fact infinite), we impose
 * an upper bound that is twice the strike price. This is a somewhat arbitrary
 * choice to be twice the strike price, but it is large enough to see the
 * interesting parts of the solution.
 * 
 * @code
 *   template <int dim>
 *   BlackScholes<dim>::BlackScholes()
 *     : maximum_stock_price(1.)
 *     , maturity_time(1.)
 *     , asset_volatility(.2)
 *     , interest_rate(0.05)
 *     , strike_price(0.5)
 *     , fe(1)
 *     , dof_handler(triangulation)
 *     , time(0.0)
 *     , theta(0.5)
 *     , n_cycles(4)
 *     , n_time_steps(5000)
 *   {
 *     Assert(dim == 1, ExcNotImplemented());
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholessetup_systemcode"></a> 
 * <h4><code>BlackScholes::setup_system</code></h4>
 * 

 * 
 * The next function sets up the DoFHandler object, computes
 * the constraints, and sets the linear algebra objects to their correct
 * sizes. We also compute the mass matrix here by calling a function from the
 * library. We will compute the other 3 matrices next, because these need to
 * be computed 'by hand'.
 *   

 * 
 * Note, that the time step is initialized here because the maturity time was
 * needed to compute the time step.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::setup_system()
 *   {
 *     dof_handler.distribute_dofs(fe);
 * 
 *     time_step = maturity_time / n_time_steps;
 * 
 *     constraints.clear();
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     constraints.close();
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs());
 *     DoFTools::make_sparsity_pattern(dof_handler,
 *                                     dsp,
 *                                     constraints,
 *                                     /*keep_constrained_dofs = */ true);
 *     sparsity_pattern.copy_from(dsp);
 * 
 *     mass_matrix.reinit(sparsity_pattern);
 *     laplace_matrix.reinit(sparsity_pattern);
 *     a_matrix.reinit(sparsity_pattern);
 *     b_matrix.reinit(sparsity_pattern);
 *     system_matrix.reinit(sparsity_pattern);
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler,
 *                                       QGauss<dim>(fe.degree + 1),
 *                                       mass_matrix);
 * 
 * @endcode
 * 
 * Below is the code to create the Laplace matrix with non-constant
 * coefficients. This corresponds to the matrix D in the introduction. This
 * non-constant coefficient is represented in the
 * <code>current_coefficient</code> variable.
 * 
 * @code
 *     const unsigned int dofs_per_cell = fe.dofs_per_cell;
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     QGauss<dim>        quadrature_formula(fe.degree + 1);
 *     FEValues<dim>      fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0.;
 *         fe_values.reinit(cell);
 *         for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *           {
 *             const double current_coefficient =
 *               fe_values.quadrature_point(q_index).square();
 *             for (const unsigned int i : fe_values.dof_indices())
 *               {
 *                 for (const unsigned int j : fe_values.dof_indices())
 *                   cell_matrix(i, j) +=
 *                     (current_coefficient *              // (x_q)^2
 *                      fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
 *                      fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
 *                      fe_values.JxW(q_index));           // dx
 *               }
 *           }
 *         cell->get_dof_indices(local_dof_indices);
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             for (const unsigned int j : fe_values.dof_indices())
 *               laplace_matrix.add(local_dof_indices[i],
 *                                  local_dof_indices[j],
 *                                  cell_matrix(i, j));
 *           }
 *       }
 * 
 * @endcode
 * 
 * Now we will create the A matrix. Below is the code to create the matrix A
 * as discussed in the introduction. The non constant coefficient is again
 * represented in  the <code>current_coefficient</code> variable.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0.;
 *         fe_values.reinit(cell);
 *         for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *           {
 *             const Tensor<1, dim> current_coefficient =
 *               fe_values.quadrature_point(q_index);
 *             for (const unsigned int i : fe_values.dof_indices())
 *               {
 *                 for (const unsigned int j : fe_values.dof_indices())
 *                   {
 *                     cell_matrix(i, j) +=
 *                       (current_coefficient *               // x_q
 *                        fe_values.shape_grad(i, q_index) *  // grad phi_i(x_q)
 *                        fe_values.shape_value(j, q_index) * // phi_j(x_q)
 *                        fe_values.JxW(q_index));            // dx
 *                   }
 *               }
 *           }
 *         cell->get_dof_indices(local_dof_indices);
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             for (const unsigned int j : fe_values.dof_indices())
 *               a_matrix.add(local_dof_indices[i],
 *                            local_dof_indices[j],
 *                            cell_matrix(i, j));
 *           }
 *       }
 * 
 * @endcode
 * 
 * Finally we will create the matrix B. Below is the code to create the
 * matrix B as discussed in the introduction. The non constant coefficient
 * is again represented in the <code>current_coefficient</code> variable.
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       {
 *         cell_matrix = 0.;
 *         fe_values.reinit(cell);
 *         for (const unsigned int q_index : fe_values.quadrature_point_indices())
 *           {
 *             const Tensor<1, dim> current_coefficient =
 *               fe_values.quadrature_point(q_index);
 *             for (const unsigned int i : fe_values.dof_indices())
 *               {
 *                 for (const unsigned int j : fe_values.dof_indices())
 *                   cell_matrix(i, j) +=
 *                     (current_coefficient *               // x_q
 *                      fe_values.shape_value(i, q_index) * // phi_i(x_q)
 *                      fe_values.shape_grad(j, q_index) *  // grad phi_j(x_q)
 *                      fe_values.JxW(q_index));            // dx
 *               }
 *           }
 *         cell->get_dof_indices(local_dof_indices);
 *         for (const unsigned int i : fe_values.dof_indices())
 *           {
 *             for (const unsigned int j : fe_values.dof_indices())
 *               b_matrix.add(local_dof_indices[i],
 *                            local_dof_indices[j],
 *                            cell_matrix(i, j));
 *           }
 *       }
 * 
 *     solution.reinit(dof_handler.n_dofs());
 *     system_rhs.reinit(dof_handler.n_dofs());
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholessolve_time_stepcode"></a> 
 * <h4><code>BlackScholes::solve_time_step</code></h4>
 * 

 * 
 * The next function is the one that solves the actual linear system for a
 * single time step. The only interesting thing here is that the matrices
 * we have built are symmetric positive definite, so we can use the
 * conjugate gradient method.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::solve_time_step()
 *   {
 *     SolverControl                          solver_control(1000, 1e-12);
 *     SolverCG<Vector<double>>               cg(solver_control);
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner;
 *     preconditioner.initialize(system_matrix, 1.0);
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 *     constraints.distribute(solution);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholesadd_results_for_outputcode"></a> 
 * <h4><code>BlackScholes::add_results_for_output</code></h4>
 * 

 * 
 * This is simply the function to stitch the solution pieces together. For
 * this, we create a new layer at each time, and then add the solution vector
 * for that timestep. The function then stitches this together with the old
 * solutions using 'build_patches'.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::add_results_for_output()
 *   {
 *     data_out_stack.new_parameter_value(time, time_step);
 *     data_out_stack.attach_dof_handler(dof_handler);
 *     data_out_stack.add_data_vector(solution, solution_names);
 *     data_out_stack.build_patches(2);
 *     data_out_stack.finish_parameter_value();
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholesrefine_gridcode"></a> 
 * <h4><code>BlackScholes::refine_grid</code></h4>
 * 

 * 
 * It is somewhat unnecessary to have a function for the global refinement
 * that we do. The reason for the function is to allow for the possibility of
 * an adaptive refinement later.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::refine_grid()
 *   {
 *     triangulation.refine_global(1);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholesprocess_solutioncode"></a> 
 * <h4><code>BlackScholes::process_solution</code></h4>
 * 

 * 
 * This is where we calculate the convergence and error data to evaluate the
 * effectiveness of the program. Here, we calculate the $L^2$, $H^1$ and
 * $L^{\infty}$ norms.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::process_solution()
 *   {
 *     Solution<dim> sol(maturity_time);
 *     sol.set_time(time);
 *     Vector<float> difference_per_cell(triangulation.n_active_cells());
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       sol,
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe.degree + 1),
 *                                       VectorTools::L2_norm);
 *     const double L2_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::L2_norm);
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       sol,
 *                                       difference_per_cell,
 *                                       QGauss<dim>(fe.degree + 1),
 *                                       VectorTools::H1_seminorm);
 *     const double H1_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::H1_seminorm);
 *     const QTrapezoid<1>  q_trapezoid;
 *     const QIterated<dim> q_iterated(q_trapezoid, fe.degree * 2 + 1);
 *     VectorTools::integrate_difference(dof_handler,
 *                                       solution,
 *                                       sol,
 *                                       difference_per_cell,
 *                                       q_iterated,
 *                                       VectorTools::Linfty_norm);
 *     const double Linfty_error =
 *       VectorTools::compute_global_error(triangulation,
 *                                         difference_per_cell,
 *                                         VectorTools::Linfty_norm);
 *     const unsigned int n_active_cells = triangulation.n_active_cells();
 *     const unsigned int n_dofs         = dof_handler.n_dofs();
 *     convergence_table.add_value("cells", n_active_cells);
 *     convergence_table.add_value("dofs", n_dofs);
 *     convergence_table.add_value("L2", L2_error);
 *     convergence_table.add_value("H1", H1_error);
 *     convergence_table.add_value("Linfty", Linfty_error);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholeswrite_convergence_tablecode"></a> 
 * <h4><code>BlackScholes::write_convergence_table</code> </h4>
 * 

 * 
 * This next part is building the convergence and error tables. By this, we
 * need to set the settings for how to output the data that was calculated
 * during <code>BlackScholes::process_solution</code>. First, we will create
 * the headings and set up the cells properly. During this, we will also
 * prescribe the precision of our results. Then we will write the calculated
 * errors based on the $L^2$, $H^1$, and $L^{\infty}$ norms to the console and
 * to the error LaTeX file.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::write_convergence_table()
 *   {
 *     convergence_table.set_precision("L2", 3);
 *     convergence_table.set_precision("H1", 3);
 *     convergence_table.set_precision("Linfty", 3);
 *     convergence_table.set_scientific("L2", true);
 *     convergence_table.set_scientific("H1", true);
 *     convergence_table.set_scientific("Linfty", true);
 *     convergence_table.set_tex_caption("cells", "\\# cells");
 *     convergence_table.set_tex_caption("dofs", "\\# dofs");
 *     convergence_table.set_tex_caption("L2", "@fL^2@f-error");
 *     convergence_table.set_tex_caption("H1", "@fH^1@f-error");
 *     convergence_table.set_tex_caption("Linfty", "@fL^\\infty@f-error");
 *     convergence_table.set_tex_format("cells", "r");
 *     convergence_table.set_tex_format("dofs", "r");
 *     std::cout << std::endl;
 *     convergence_table.write_text(std::cout);
 *     std::string error_filename = "error";
 *     error_filename += "-global";
 *     error_filename += ".tex";
 *     std::ofstream error_table_file(error_filename);
 *     convergence_table.write_tex(error_table_file);
 * 
 * @endcode
 * 
 * Next, we will make the convergence table. We will again write this to
 * the console and to the convergence LaTeX file.
 * 
 * @code
 *     convergence_table.add_column_to_supercolumn("cells", "n cells");
 *     std::vector<std::string> new_order;
 *     new_order.emplace_back("n cells");
 *     new_order.emplace_back("H1");
 *     new_order.emplace_back("L2");
 *     convergence_table.set_column_order(new_order);
 *     convergence_table.evaluate_convergence_rates(
 *       "L2", ConvergenceTable::reduction_rate);
 *     convergence_table.evaluate_convergence_rates(
 *       "L2", ConvergenceTable::reduction_rate_log2);
 *     convergence_table.evaluate_convergence_rates(
 *       "H1", ConvergenceTable::reduction_rate);
 *     convergence_table.evaluate_convergence_rates(
 *       "H1", ConvergenceTable::reduction_rate_log2);
 *     std::cout << std::endl;
 *     convergence_table.write_text(std::cout);
 *     std::string conv_filename = "convergence";
 *     conv_filename += "-global";
 *     switch (fe.degree)
 *       {
 *         case 1:
 *           conv_filename += "-q1";
 *           break;
 *         case 2:
 *           conv_filename += "-q2";
 *           break;
 *         default:
 *           Assert(false, ExcNotImplemented());
 *       }
 *     conv_filename += ".tex";
 *     std::ofstream table_file(conv_filename);
 *     convergence_table.write_tex(table_file);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="codeBlackScholesruncode"></a> 
 * <h4><code>BlackScholes::run</code></h4>
 * 

 * 
 * Now we get to the main driver of the program. This is where we do all the
 * work of looping through the time steps and calculating the solution vector
 * each time. Here at the top, we set the initial refinement value and then
 * create a mesh. Then we refine this mesh once. Next, we set up the
 * data_out_stack object to store our solution. Finally, we start a for loop
 * to loop through the cycles. This lets us recalculate a solution for each
 * successive mesh refinement. At the beginning of each iteration, we need to
 * reset the time and time step number. We introduce an if statement to
 * accomplish this because we don't want to do this on the first iteration.
 * 
 * @code
 *   template <int dim>
 *   void BlackScholes<dim>::run()
 *   {
 *     GridGenerator::hyper_cube(triangulation, 0.0, maximum_stock_price, true);
 *     triangulation.refine_global(0);
 * 
 *     solution_names.emplace_back("u");
 *     data_out_stack.declare_data_vector(solution_names,
 *                                        DataOutStack<dim>::dof_vector);
 * 
 *     Vector<double> vmult_result;
 *     Vector<double> forcing_terms;
 * 
 *     for (unsigned int cycle = 0; cycle < n_cycles; cycle++)
 *       {
 *         if (cycle != 0)
 *           {
 *             refine_grid();
 *             time = 0.0;
 *           }
 * 
 *         setup_system();
 * 
 *         std::cout << std::endl
 *                   << "===========================================" << std::endl
 *                   << "Cycle " << cycle << ':' << std::endl
 *                   << "Number of active cells: "
 *                   << triangulation.n_active_cells() << std::endl
 *                   << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *                   << std::endl
 *                   << std::endl;
 * 
 *         VectorTools::interpolate(dof_handler,
 *                                  InitialConditions<dim>(strike_price),
 *                                  solution);
 * 
 *         if (cycle == (n_cycles - 1))
 *           {
 *             add_results_for_output();
 *           }
 * 
 * @endcode
 * 
 * Next, we run the main loop which runs until we exceed the maturity
 * time. We first compute the right hand side of the equation, which is
 * described in the introduction. Recall that it contains the term
 * $\left[-\frac{1}{4}k_n\sigma^2\mathbf{D}-k_nr\mathbf{M}+k_n\sigma^2
 * \mathbf{B}-k_nr\mathbf{A}+\mathbf{M}\right]V^{n-1}$. We put these
 * terms into the variable system_rhs, with the help of a temporary
 * vector:
 * 
 * @code
 *         vmult_result.reinit(dof_handler.n_dofs());
 *         forcing_terms.reinit(dof_handler.n_dofs());
 *         for (unsigned int timestep_number = 0; timestep_number < n_time_steps;
 *              ++timestep_number)
 *           {
 *             time += time_step;
 * 
 *             if (timestep_number % 1000 == 0)
 *               std::cout << "Time step " << timestep_number << " at t=" << time
 *                         << std::endl;
 * 
 *             mass_matrix.vmult(system_rhs, solution);
 * 
 *             laplace_matrix.vmult(vmult_result, solution);
 *             system_rhs.add(
 *               (-1) * (1 - theta) * time_step *
 *                 Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
 *               vmult_result);
 *             mass_matrix.vmult(vmult_result, solution);
 * 
 *             system_rhs.add((-1) * (1 - theta) * time_step * interest_rate * 2,
 *                            vmult_result);
 * 
 *             a_matrix.vmult(vmult_result, solution);
 *             system_rhs.add((-1) * time_step * interest_rate, vmult_result);
 * 
 *             b_matrix.vmult(vmult_result, solution);
 *             system_rhs.add(
 *               (-1) * Utilities::fixed_power<2, double>(asset_volatility) *
 *                 time_step * 1,
 *               vmult_result);
 * 
 * @endcode
 * 
 * The second piece is to compute the contributions of the source
 * terms. This corresponds to the term $-k_n\left[\frac{1}{2}F^{n-1}
 * +\frac{1}{2}F^n\right]$. The following code calls
 * VectorTools::create_right_hand_side to compute the vectors $F$,
 * where we set the time of the right hand side (source) function
 * before we evaluate it. The result of this all ends up in the
 * forcing_terms variable:
 * 
 * @code
 *             RightHandSide<dim> rhs_function(asset_volatility, interest_rate);
 *             rhs_function.set_time(time);
 *             VectorTools::create_right_hand_side(dof_handler,
 *                                                 QGauss<dim>(fe.degree + 1),
 *                                                 rhs_function,
 *                                                 forcing_terms);
 *             forcing_terms *= time_step * theta;
 *             system_rhs -= forcing_terms;
 * 
 *             rhs_function.set_time(time - time_step);
 *             VectorTools::create_right_hand_side(dof_handler,
 *                                                 QGauss<dim>(fe.degree + 1),
 *                                                 rhs_function,
 *                                                 forcing_terms);
 *             forcing_terms *= time_step * (1 - theta);
 *             system_rhs -= forcing_terms;
 * 
 * @endcode
 * 
 * Next, we add the forcing terms to the ones that come from the
 * time stepping, and also build the matrix $\left[\mathbf{M}+
 * \frac{1}{4}k_n\sigma^2\mathbf{D}+k_nr\mathbf{M}\right]$ that we
 * have to invert in each time step. The final piece of these
 * operations is to eliminate hanging node constrained degrees of
 * freedom from the linear system:
 * 
 * @code
 *             system_matrix.copy_from(mass_matrix);
 *             system_matrix.add(
 *               (theta)*time_step *
 *                 Utilities::fixed_power<2, double>(asset_volatility) * 0.5,
 *               laplace_matrix);
 *             system_matrix.add((time_step)*interest_rate * theta * (1 + 1),
 *                               mass_matrix);
 * 
 *             constraints.condense(system_matrix, system_rhs);
 * 
 * @endcode
 * 
 * There is one more operation we need to do before we can solve it:
 * boundary values. To this end, we create a boundary value object,
 * set the proper time to the one of the current time step, and
 * evaluate it as we have done many times before. The result is used
 * to also set the correct boundary values in the linear system:
 * 
 * @code
 *             {
 *               RightBoundaryValues<dim> right_boundary_function(strike_price,
 *                                                                interest_rate);
 *               LeftBoundaryValues<dim>  left_boundary_function;
 *               right_boundary_function.set_time(time);
 *               left_boundary_function.set_time(time);
 *               std::map<types::global_dof_index, double> boundary_values;
 *               VectorTools::interpolate_boundary_values(dof_handler,
 *                                                        0,
 *                                                        left_boundary_function,
 *                                                        boundary_values);
 *               VectorTools::interpolate_boundary_values(dof_handler,
 *                                                        1,
 *                                                        right_boundary_function,
 *                                                        boundary_values);
 *               MatrixTools::apply_boundary_values(boundary_values,
 *                                                  system_matrix,
 *                                                  solution,
 *                                                  system_rhs);
 *             }
 * 
 * @endcode
 * 
 * With this out of the way, all we have to do is solve the system,
 * generate graphical data on the last cycle, and create the
 * convergence table data.
 * 
 * @code
 *             solve_time_step();
 * 
 *             if (cycle == (n_cycles - 1))
 *               {
 *                 add_results_for_output();
 *               }
 *           }
 * #ifdef MMS
 *         process_solution();
 * #endif
 *       }
 * 
 *     const std::string filename = "solution.vtk";
 *     std::ofstream     output(filename);
 *     data_out_stack.write_vtk(output);
 * 
 * #ifdef MMS
 *     write_convergence_table();
 * #endif
 *   }
 * 
 * } // namespace BlackScholesSolver
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodemaincodeFunction"></a> 
 * <h3>The <code>main</code> Function</h3>
 * 

 * 
 * Having made it this far, there is, again, nothing much to discuss for the
 * main function of this program: it looks like all such functions since step-6.
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace BlackScholesSolver;
 * 
 *       BlackScholes<1> black_scholes_solver;
 *       black_scholes_solver.run();
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
 *     }
 *   return 0;
 * }
 * @endcode
examples/step-78/doc/results.dox



<a name="Results"></a><h1>Results</h1>



下面是该程序的输出。

@code
===========================================
Number of active cells: 1
Number of degrees of freedom: 2


Time step 0 at t=0.0002
[...]


Cycle 7:
Number of active cells: 128
Number of degrees of freedom: 129


Time step 0 at t=0.0002
Time step 1000 at t=0.2002
Time step 2000 at t=0.4002
Time step 3000 at t=0.6002
Time step 4000 at t=0.8002


cells dofs    L2        H1      Linfty
    1    2 1.667e-01 5.774e-01 2.222e-01
    2    3 3.906e-02 2.889e-01 5.380e-02
    4    5 9.679e-03 1.444e-01 1.357e-02
    8    9 2.405e-03 7.218e-02 3.419e-03
   16   17 5.967e-04 3.609e-02 8.597e-04
   32   33 1.457e-04 1.804e-02 2.155e-04
   64   65 3.307e-05 9.022e-03 5.388e-05
  128  129 5.016e-06 4.511e-03 1.342e-05


n cells         H1                  L2
      1 5.774e-01    -    - 1.667e-01    -    -
      2 2.889e-01 2.00 1.00 3.906e-02 4.27 2.09
      4 1.444e-01 2.00 1.00 9.679e-03 4.04 2.01
      8 7.218e-02 2.00 1.00 2.405e-03 4.02 2.01
     16 3.609e-02 2.00 1.00 5.967e-04 4.03 2.01
     32 1.804e-02 2.00 1.00 1.457e-04 4.10 2.03
     64 9.022e-03 2.00 1.00 3.307e-05 4.41 2.14
    128 4.511e-03 2.00 1.00 5.016e-06 6.59 2.72
@endcode



更有趣的是收敛表的输出。它们被输出到控制台，以及一个LaTeX文件中。收敛表如上所示。在这里，你可以看到，相对于 $H^1$ -norm，解决方案的收敛率为 $\mathcal{O}(h)$ ，相对于 $L^2$ -norm，解决方案的收敛率为 $\mathcal{O}(h^2)$ 。


下面是解决方案的可视化。

<div style="text-align:center;"> <img src="https://www.dealii.org/images/steps/developer/step-78.mms-solution.png" alt="MMS问题的解决方案。"> </div>


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-78.cc"
*/
