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
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


The Black-Scholes equation is a partial differential equation that falls a bit
out of the ordinary scheme. It describes what the fair price of a "European
call" stock option is. Without going into too much detail, a stock "option" is
a contract one can buy from a bank that allows me, but not requires me, to buy
a specific stock at a fixed price $K$ at a fixed future time $T$ in the
future. The question one would then want to answer as a buyer of such an
option is "How much do I think such a contract is worth?", or as the seller
"How much do I need to charge for this contract?", both as a function of the
time $t<T$ before the contract is up at time $T$ and as a function of the stock
price $S$. Fischer Black and Myron Scholes derived a partial differential
equation for the fair price $V(S,t)$ for such options under the assumption that
stock prices exhibit random price fluctuations with a given level of
"volatility" plus a background exponential price increase (which one can think
of as the inflation rate that simply devalues all money over time). For their
work, Black and Scholes received the Nobel Prize in Economic Sciences in 1997,
making this the first tutorial program dealing with a problem for which someone
has gotten a Nobel Prize @cite black1973pricing.

The equation reads as follows:
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
where
@f{align*}{
    V(S,t): && \text{Value of call option at time t and asset price S} \\
    \sigma: && \text{Volatility of the underlying asset} \\
    r: && \text{Risk free interest rate} \\
    K : && \text{Strike price for purchasing asset}
@f}

The way we should interpret this equation is that it is a time-dependent partial
differential equation of one "space" variable
$S$ as the price of the stock, and $V(S,t)$ is the price of the option at time
$t$ if the stock price at that time were $S$.

<a name="Particularitiesoftheequationsystem"></a><h3>Particularities of the equation system</h3>


There are a number of oddities in this equation that are worth discussing before
moving on to its numerical solution. First, the "spatial" domain
$\Omega\subset\mathbb{R}$ is unbounded, and thus $S$ can be unbounded in value.
This is because there may be a practical upper bound for stock prices, but not a
conceptual one. The boundary conditions $V(S,t)\rightarrow S$ as
$S\rightarrow \infty$ can then be interpreted as follows: What is the value of
an option that allows me to buy a stock at price $K$ if the stock price (today
or at time $t=T$) is $S\gg K$? One would expect that it is $V\approx S-K$ plus
some adjustment for inflation, or, if we really truly consider huge values of
$S$, we can neglect $K$ and arrive at the statement that the boundary values at
the infinite boundary should be of the form $V\rightarrow S$ as stated above.

In practice, for us to use a finite element method to solve this, we are going
to need to bound $\Omega$. Since this equation describes prices, and it doesn't
make sense to talk about prices being negative, we will set the lower bound of
$\Omega$ to be 0. Then, for an upper bound, we will choose a very large number,
one that $S$ is not very likely to ever get to. We will call this $S_\text{max}$.
So, $\Omega=[0,S_\text{max}]$.

Second, after truncating the domain, we need to ask what boundary values we
should pose at this now finite boundary. To take care of this, we use "put-call"
parity @cite stoll1969relationship. A "pull option" is one in which we are
allowed, but not required, to *sell* a stock at price $K$ to someone at a future
time $T$. This says
@f{align*}{
    V(S,t)+Ke^{-r(T-t)}=P(S,t)+S
@f}
where $V(S,t)$ is the value of the call option, and $P(S,t)$ is the value of the
put option. Since we expect $P(S,t) \rightarrow 0$ as $S \rightarrow \infty$,
this says
@f{align*}{
    V(S,t) \rightarrow S-Ke^{-r(T-t)},
@f}
and we can use this as a reasonable boundary condition at our finite point
$S_\text{max}$.

The second complication of the Block-Scholes equation is that we are given a
final condition, and not an initial condition. This is because we know what the
option is worth at time $t=T$: If the stock price at $T$ is $S<K$, then we have
no incentive to use our option of buying a price $K$ because we can buy that stock
for cheaper on the open market. So $V(S,T)=0$ for $S<K$. On the other hand, if
at time $T$ we have $S>K$, then we can buy the stock at price $K$ via the option
and immediately sell it again on the market for price $S$, giving me a profit of
$S-K$. In other words, $V(S,T)=S-K$ for $S>K$. So, we only know
values for $V$ at the *end time* but not the initial time -- in fact, finding
out what a fair price at the current time (conventionally taken to be $t=0$) is
what solving these equations is all about.

This means that this is not an equation that is posed going forward in
time, but in fact going *backward* in time. Thus it makes sense to solve this
problem in reverse by making the change of variables $\tau=T-t$ where now $\tau$
denotes the time before the strike time $T$.

With all of this, we finally end up with the following problem:
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

Conceptually, this is an advection-diffusion-reaction problem for the variable
$V$: There is both a second-order derivative diffusion term, a first-order
derivative advection term, and a zeroth-order reaction term.
We can expect this problem to be a little bit forgiving in practice because for
realistic values of the coefficients, it is diffusive dominated. But, because of
the advective terms in the problem, we will have to be careful with mesh
refinement and time step choice. There is also the issue that the diffusion term
 is written in a non-conservative form and so integration by parts is not
 immediately obvious. This will be discussed in the next section.

<a name="Schemeforthenumericalsolution"></a><h3>Scheme for the numerical solution</h3>


We will solve this problem using an IMEX method. In particular, we first discretize
in time with the theta method and will later pick different values of theta for
the advective and diffusive terms.
Let $V^n(S)$ approximate $V(S,\tau_n)$:
@f{align*}{
    0=&-\frac{V^n(S)-V^{n-1}(S)}{k_n} \\
    &+\frac{\sigma^2S^2}{2}\left[(1-\theta)\frac{d^2V^{n-1}(S)}{dS^2} + \
    \theta \frac{d^2V^{n}(S)}{dS^2}\right] \\
    &+rS\left[(1-\theta)\frac{dV^{n-1}(S)}{dS} + \
    \theta\frac{dV^{n}(S)}{dS}\right]  \\
    &-r\left[(1-\theta)V^{n-1}(S) + \theta V^n(S)\right]
@f}
Here, $k_n=\tau_n-\tau_{n-1}$ is the time step size. Given this time
discretization, we can proceed to discretize space by multiplying with test
functions and then integrating by parts. Because there are some interesting
details in this due to the advective and non-advective terms in this equation,
this process will be explained in detail.

So, we begin by multiplying by test functions, $\{\phi_i(S)\}_{i\in\mathbb{N}}$:
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


As usual, (1) becomes $-\textbf{M}V^n+\textbf{M}V^{n-1}$ and (4) becomes
$k_n\left[-r(1-\theta)\textbf{M}V^{n-1} - \theta r\textbf{M}V^n\right]$, where
$\textbf{M}_{i,j}=\left(\phi_i(S),\phi_j(S)\right)$, and where we have taken the
liberty of denoting by $V$ not only the function $V(S)$ but also the vector of
nodal values after discretization.

The interesting parts come from (2) and (3).


For (2), we have:
@f{align*}{
    &k_n\int_0^{S_\text{max}}\phi_i(S)\left[\frac{\sigma^2S^2}{2} \
     \left[(1-\theta)\frac{d^2V^{n-1}(S)}{dS^2} + \
     \theta \frac{d^2V^{n}(S)}{dS^2}\right]\right]dS \\
    &=k_n(1-\theta)\int_0^{S_\text{max}}\phi_i(S)\frac{\sigma^2S^2}{2} \
     \frac{d^2V^{n-1}(S)}{dS^2} \
    +k_n\theta\int_0^{S_\text{max}}\phi_i(S)\frac{\sigma^2S^2}{2} \
     \frac{d^2V^{n}(S)}{dS^2}
@f}

There are two integrals here, that are more or less the same, with the
differences being a slightly different coefficient in front of the integral,
and a different time step for V. Therefore, we will outline this integral in the
general case, and account for the differences at the end. So, consider the
general integral, which we will solve using integration by parts:
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

So, after adding in the constants and exchanging $V^n$ for $V^{n-1}$ where
applicable, we arrive at the following for (2):
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
But, because the matrix $\textbf{B}$ involves an advective term, we will choose
$\theta=0$ there -- in other words, we use an explicit Euler method to treat
advection. Conversely, since the matrix $\textbf{D}$ involves the diffusive term,
we will choose $\theta=1/2$ there -- i.e., we treat diffusion using the second
order Crank-Nicolson method.

So, we arrive at the following:
@f{align*}{
    k_n\left[-\frac{1}{4}\sigma^2\textbf{D}V^{n-1} \
    -\frac{1}{4}\sigma^2\textbf{D}V^n \
    - \sigma^2\textbf{B}V^{n-1}\right]
@f}

Now, to handle (3). For this, we will again proceed by considering the general
case like above.

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

So, again after adding in the constants and exchanging $V^n$ for $V^{n-1}$ where
applicable, we arrive at the following for (3):
@f{align*}{
    &k_n\int_0^{S_\text{max}}\phi_i(S)\left[rS\left[(1-\theta)
        \frac{dV^{n-1}(S)}{dS} +\
     \theta\frac{dV^{n}(S)}{dS}\right]\right]dS \\
    &= k_n\left[-(1-\theta)r\textbf{M}V^{n-1} -(1-\theta)r\textbf{A}V^{n-1}\
    -\theta r\textbf{M}V^n -\theta r\textbf{A}V^n\right]
@f}
Just as before, we will use $\theta=0$ for the matrix $\textbf{A}$ and
$\theta=\frac{1}{2}$ for the matrix $\textbf{M}$. So, we arrive at the
following for (3):
@f{align*}{
    k_n\left[-\frac{1}{2}r\textbf{M}V^{n-1} - \frac{1}{2}r\textbf{M}V^n \
    -r\textbf{A}V^{n-1}\right]
@f}

Now, putting everything together, we obtain the following discrete form for the
Black-Scholes Equation:
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
So, altogether we have:

@f{equation}{
    0 = \textbf{M}V^n - \textbf{M}V^{n-1} +\
    k_n\left[ \frac{1}{4}\sigma^2\textbf{D}V^{n-1} +\
    \frac{1}{4}\sigma^2\textbf{D}V^n + r\textbf{M}V^{n-1} + r\textbf{M}V^n  +\
    \sigma^2\textbf{B}V^{n-1} + r\textbf{A}V^{n-1}\right]\tag{*}
@f}

As usual, we can write this with the unknown quantities on the left and the
known ones on the right. This leads to the following linear system that would
have to be solved in each time step:

@f{align*}{
    \left[\textbf{M}+\frac{1}{4}k_n\sigma^2\textbf{D}+k_nr\textbf{M}\right]V^n\
     =\
    \left[-\frac{1}{4}k_n\sigma^2\textbf{D}-\
    k_nr\textbf{M}+k_n\sigma^2\textbf{B}-\
    k_nr\textbf{A}+\textbf{M}\right]V^{n-1}
@f}




<a name="TestCase"></a><h3>Test Case</h3>

For this program, we will use the Method of Manufactured Solutions (MMS) to test
 that it is working correctly. This means that we will choose our solution to be
  a certain function similar to step-7. For our case, we will use:
@f{align*}{
    V(S,\tau) = -\tau^2 - S^2 + 6\tag{**}
@f}
This means that, using our PDE, we arrive at the following problem:
@f{align*}{
    &-\frac{\partial V}{\partial \tau} +\
    \frac{\sigma^2S^2}{2}\frac{\partial^2 V}{\partial S^2} +\
    rS\frac{\partial V}{\partial S} - rV = f(S,\tau) \\
    &V(0,\tau) = -\tau^2 + 6 \\
    &V(S_\text{max}, \tau) = -S_\text{max}^2 - \tau^2 + 6 \\
    &V(S, 0) = -S^2 + 6
@f}
Where, $f(S,\tau) = 2\tau - \sigma^2S^2 - 2rS^2 - r(-\tau^2 - S^2 + 6)$.
This set-up now has right hand sides for the equation itself and for the
boundary conditions at $S=0$ that we did not have before, along with "final"
conditions (or, with $\tau$-time "initial conditions") that do not match the
real situation. We will implement this in such a way in the code that it is easy
to exchange -- the introduction of the changes above is just meant to enable the
 use of a manufactured solution.

If the program is working correctly, then it should produce (**) as the
solution. This does mean that we need to modify our variational form somewhat to
account for the non-zero right hand side.

First, we define the following:
@f{align*}{
    F^n_i = \left(\phi_i(S), f^n(S)\right), && \text{where } f^n(S) =\
     f(S,\tau_n)
@f}
So, we arrive at the new equation:

@f{align*}{
    \left[\textbf{M}+\frac{1}{4}k_n\sigma^2\textbf{D}+k_nr\textbf{M}\right]V^n\
     =\
     \left[-\frac{1}{4}k_n\sigma^2\textbf{D}-\
     k_nr\textbf{M}+k_n\sigma^2\textbf{B}-\
     k_nr\textbf{A}+\textbf{M}\right]V^{n-1} -\
      k_n\left[\frac{1}{2}F^{n-1}+\frac{1}{2}F^n\right]
@f}

We then solve this equation as outlined above.
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
<a name="Results"></a><h1>Results</h1>



Below is the output of the program:
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

What is more interesting is the output of the convergence tables. They are
outputted into the console, as well into a LaTeX file. The convergence tables
are shown above. Here, you can see that the the solution has a convergence rate
of $\mathcal{O}(h)$ with respect to the $H^1$-norm, and the solution has a convergence rate
of $\mathcal{O}(h^2)$ with respect to the $L^2$-norm.


Below is the visualization of the solution.

<div style="text-align:center;">
  <img src="https://www.dealii.org/images/steps/developer/step-78.mms-solution.png"
       alt="Solution of the MMS problem.">
</div>
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-78.cc"
*/
