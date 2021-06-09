/**
@page step_37 The step-37 tutorial program
This tutorial depends on step-16, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The test case</a>
        <li><a href="#Matrixvectorproductimplementation">Matrix-vector product implementation</a>
        <li><a href="#Combinationwithmultigrid">Combination with multigrid</a>
        <li><a href="#UsingCPUdependentinstructionsvectorization">Using CPU-dependent instructions (vectorization)</a>
        <li><a href="#Runningmultigridonlargescaleparallelcomputers">Running multigrid on large-scale parallel computers</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Matrixfreeimplementation">Matrix-free implementation</a>
      <ul>
        <li><a href="#Computationofcoefficient">Computation of coefficient</a>
        <li><a href="#LocalevaluationofLaplaceoperator">Local evaluation of Laplace operator</a>
      </ul>
        <li><a href="#LaplaceProblemclass">LaplaceProblem class</a>
      <ul>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_rhs">LaplaceProblem::assemble_rhs</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Comparisonwithasparsematrix">Comparison with a sparse matrix</a>
        <li><a href="#ResultsforlargescaleparallelcomputationsonSuperMUC"> Results for large-scale parallel computations on SuperMUC</a>
        <li><a href="#Adaptivity"> Adaptivity</a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions</a>
      <ul>
        <li><a href="#Kellyerrorestimator"> Kelly error estimator </a>
        <li><a href="#Sharedmemoryparallelization"> Shared-memory parallelization</a>
        <li><a href="#InhomogeneousDirichletboundaryconditions"> Inhomogeneous Dirichlet boundary conditions </a>
      <ul>
        <li><a href="#UseFEEvaluationread_dof_values_plaintoavoidresolvingconstraints"> Use FEEvaluation::read_dof_values_plain() to avoid resolving constraints </a>
        <li><a href="#UseLaplaceOperatorwithasecondAffineConstraintsobjectwithoutDirichletconditions"> Use LaplaceOperator with a second AffineConstraints object without Dirichlet conditions </a>
    </ul>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-37/doc/intro.dox

 <br> 

<i>
This program was contributed by Katharina Kormann and Martin
Kronbichler.


The algorithm for the matrix-vector product is based on the article <a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic interface
for parallel cell-based finite element operator application</a><a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">A generic interface
for parallel cell-based finite element operator application</a> by Martin
Kronbichler and Katharina Kormann, Computers and Fluids 63:135&ndash;147,
2012, and the paper &quot;Parallel finite element operator application: Graph
partitioning and coloring&quot; by Katharina Kormann and Martin Kronbichler
in: Proceedings of the 7th IEEE International Conference on e-Science, 2011.


This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA). The
large-scale computations shown in the results section of this tutorial program
were supported by Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu)
by providing computing time on the GCS Supercomputer SuperMUC at Leibniz
Supercomputing Centre (LRZ, www.lrz.de) through project id pr83te. </i> 。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个例子展示了如何在超立方体上实现一个无矩阵的方法，即不明确存储矩阵元素的方法，用于具有可变系数的二阶泊松方程。该线性系统将用多网格方法求解，并使用MPI的大规模并行性。

无矩阵方法的主要动机是，在今天的处理器上，对主内存的访问（即对不适合缓存的对象）已经成为许多偏微分方程求解器的瓶颈。为了执行基于矩阵的矩阵-向量乘积，现代CPU花在等待数据从内存到达的时间远远多于实际进行浮点乘法和加法的时间。因此，如果我们可以通过重新计算矩阵元素来代替在内存中查找矩阵元素，或者更确切地说，这些条目所代表的运算符&mdash;，我们可能会在整体运行时间方面获胜，即使这需要大量的额外浮点运算。也就是说，用一个微不足道的实现来实现这一点是不够的，我们需要真正关注细节来获得性能。这个教程程序和上面提到的论文展示了如何实现这样一个方案，并演示了可以获得的速度提升。




<a name="Thetestcase"></a><h3>The test case</h3>


在这个例子中，我们考虑泊松问题@f{eqnarray*} -
\nabla \cdot a(\mathbf x) \nabla u &=& 1, \\ u &=& 0 \quad \text{on}\
\partial \Omega @f}，其中 $a(\mathbf x)$ 是一个可变系数。下面，我们将解释如何在不明确形成矩阵的情况下实现这个问题的矩阵-向量乘积。当然，对于其他方程也可以用类似的方法进行构造。

我们选择 $\Omega=[0,1]^3$ 和 $a(\mathbf x)=\frac{1}{0.05 +
2\|\mathbf x\|^2}$ 作为域。由于系数是围绕原点对称的，但域却不是，我们最终会得到一个非对称的解决方案。




<a name="Matrixvectorproductimplementation"></a><h3>Matrix-vector product implementation</h3>


为了找出我们如何编写一个执行矩阵-向量乘积的代码，但不需要存储矩阵元素，让我们先看看一个有限元矩阵<i>A</i>是如何组装起来的。

@f{eqnarray*}
A = \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}}
P_{\mathrm{cell,{loc-glob}}}^T A_{\mathrm{cell}} P_{\mathrm{cell,{loc-glob}}}.


@f}

在这个公式中，矩阵<i>P</i><sub>cell,loc-glob</sub>是一个矩形矩阵，定义了从当前单元的局部自由度到全局自由度的索引映射。可以建立这个算子的信息通常被编码在 <code>local_dof_indices</code> 变量中，并在deal.II中用于汇编调用填充矩阵。这里，<i>A</i><sub>cell</sub>表示与<i>A</i>相关的单元矩阵。

如果我们要进行矩阵-向量乘积，因此我们可以使用

@f{eqnarray*}
y &=& A\cdot u = \left(\sum_{\text{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
A_\mathrm{cell} P_\mathrm{cell,{loc-glob}}\right) \cdot u
\\
&=& \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
A_\mathrm{cell} u_\mathrm{cell}
\\
&=& \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
v_\mathrm{cell},


@f}

其中<i>u</i><sub>cell</sub>是<i>u</i>在各单元自由度处的值，而<i>v</i><sub>cell</sub>=<i>A</i><sub>cell</sub><i>u</i><sub>cell</sub>相应为结果。因此，实现拉普拉斯的局部作用的一个天真尝试是使用以下代码。

@code
Matrixfree<dim>::vmult (Vector<double>       &dst,
                        const Vector<double> &src) const
{
  dst = 0;


  QGauss<dim>  quadrature_formula(fe.degree+1);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_gradients | update_JxW_values|
                           update_quadrature_points);


  const unsigned int   dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int   n_q_points    = quadrature_formula.size();


  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_src (dofs_per_cell),
                       cell_dst (dofs_per_cell);
  const Coefficient<dim> coefficient;
  std::vector<double> coefficient_values(n_q_points);


  std::vector<unsigned int> local_dof_indices (dofs_per_cell);


  for (const auto & cell : dof_handler.active_cell_iterators())
    {
      cell_matrix = 0;
      fe_values.reinit (cell);
      coefficient.value_list(fe_values.get_quadrature_points(),
                             coefficient_values);


      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (fe_values.shape_grad(i,q) *
                                 fe_values.shape_grad(j,q) *
                                 fe_values.JxW(q)*
                                 coefficient_values[q]);


      cell->get_dof_indices (local_dof_indices);


      for (unsigned int i=0; i<dofs_per_cell; ++i)
        cell_src(i) = src(local_dof_indices(i));


      cell_matrix.vmult (cell_dst, cell_src);


      for (unsigned int i=0; i<dofs_per_cell; ++i)
        dst(local_dof_indices(i)) += cell_dst;
    }
}
@endcode



在这里，我们忽略了边界条件以及我们可能有的任何悬空节点，尽管使用AffineConstraints类来包括这两者都不是很困难。请注意，我们首先以通常的方式生成局部矩阵，作为每个局部矩阵项的所有正交点的总和。为了形成上述公式中表达的实际乘积，我们提取细胞相关自由度的 <code>src</code> 的值（<i>P</i><sub>cell,loc-glob</sub>的作用），乘以局部矩阵（<i>A</i><sub>cell</sub>），最后把结果加到目标向量 <code>dst</code> （<i>P</i><sub>cell,loc-glob</sub><sup>T</sup>的动作，加在所有元素上）。原则上不会比这更难。

虽然这段代码是完全正确的，但它非常慢。对于每个单元，我们生成一个局部矩阵，这需要三个嵌套循环，循环长度等于局部自由度的数量来计算。然后，乘法本身是由两个嵌套循环完成的，这意味着它要便宜得多。

改善这一点的一个方法是认识到，从概念上讲，局部矩阵可以被认为是三个矩阵的乘积。

@f{eqnarray*}
A_\mathrm{cell} = B_\mathrm{cell}^T D_\mathrm{cell} B_\mathrm{cell},


@f}

对于拉普拉斯算子的例子，<i>q</i>*dim+<i>d,i</i>的第1个元素<sub>cell</sub>是由 <code>fe_values.shape_grad(i,q)[d]</code> 给出。这个矩阵由 <code>dim*n_q_points</code> 行和 @p dofs_per_cell 列组成。矩阵<i>D</i><sub>cell</sub>是对角线，包含了 <code>fe_values.JxW(q) * coefficient_values[q]</code> 的值（或者说， @p 这些值中每一个的dim副本）。这种有限元矩阵的表示方法经常可以在工程文献中找到。

当单元格矩阵被应用于一个矢量时。

@f{eqnarray*}
A_\mathrm{cell}\cdot u_\mathrm{cell} = B_\mathrm{cell}^T
D_\mathrm{cell} B_\mathrm{cell} \cdot u_\mathrm{cell},


@f}

这样就不会形成矩阵-矩阵乘积，而是每次用一个矩阵与一个矢量从右到左相乘，这样就只形成三个连续的矩阵-矢量乘积。这种方法去掉了局部矩阵计算中的三个嵌套循环，从而将一个单元格的工作复杂度从类似 $\mathcal
{O}(\mathrm{dofs\_per\_cell}^3)$ 降低到 $\mathcal
{O}(\mathrm{dofs\_per\_cell}^2)$  。对这种算法的解释是，我们首先将本地DoF上的值向量转换为正交点上的梯度向量。在第二个循环中，我们把这些梯度乘以积分权重和系数。第三次循环应用第二个梯度（转置形式），这样我们就得到了单元斗室上的（拉普拉斯）值矢量。

上述代码的瓶颈是对每一个 FEValues::reinit 的调用所做的操作，其花费的时间和其他步骤加起来差不多（至少如果网格是非结构化的；deal.II可以识别结构化网格上的梯度往往是不变的）。这当然不理想，我们希望能做得更好。reinit函数所做的是计算实空间的梯度，使用从实空间到参考单元的转换的Jacobian来转换参考单元上的梯度。这是为单元格上的每个基函数和每个正交点进行的。雅各布系数并不取决于基函数，但它在不同的正交点上通常是不同的。如果你只建立一次矩阵，就像我们在以前所有的教程程序中所做的那样，没有什么需要优化的，因为 FEValues::reinit 需要在每个单元上调用。在这个过程中，转换是在计算局部矩阵元素时应用的。

然而，在一个无矩阵的实现中，我们会经常计算这些积分，因为迭代求解器在求解过程中会多次应用矩阵。因此，我们需要考虑是否可以缓存一些在运算器应用中被重用的数据，也就是积分计算。另一方面，我们意识到我们不能缓存太多的数据，否则我们又回到了内存访问成为主导因素的情况。因此，我们不会在矩阵<i>B</i>中存储转换后的梯度，因为一般来说，对于曲线网格的每个基函数和每个元素上的正交点，它们都是不同的。

诀窍是去掉雅各布变换的因素，首先只在参考单元上应用梯度。这个操作将本地道夫上的值向量插值到正交点上的（单位坐标）梯度向量。在这里，我们首先应用我们从梯度中分解出来的雅各布，然后应用正交点的权重，最后应用转置的雅各布来准备第三个循环，通过单元格上的梯度测试并对正交点求和。

让我们再次用矩阵的方式来写。让矩阵<i>B</i><sub>cell</sub>表示与单元有关的梯度矩阵，每一行包含正交点上的值。它由矩阵与矩阵的乘积构成@f{eqnarray*} B_\mathrm{cell} =
J_\mathrm{cell}^{-\mathrm T} B_\mathrm{ref\_cell}, @f}，其中<i>B</i><sub>ref_cell</sub>表示参考单元的梯度，<i>J</i><sup>-T</sup><sub>cell</sub>表示从单位到实数单元的变换的反转置Jacobian（在变换的语言中，由<i>J</i><sup>-T</sup><sub>cell</sub>表示协变变换的操作）。<i>J</i><sup>-T</sup><sub>cell</sub>是块对角线的，块的大小等于问题的维度。每个对角线块都是雅各布变换，从参考单元到实际单元。

把事情放在一起，我们发现

@f{eqnarray*}
A_\mathrm{cell} = B_\mathrm{cell}^T D B_\mathrm{cell}
                = B_\mathrm{ref\_cell}^T J_\mathrm{cell}^{-1}
                  D_\mathrm{cell}
                  J_\mathrm{cell}^{-\mathrm T} B_\mathrm{ref\_cell},


@f}

所以我们要计算积（从右边开始计算局部积）。

@f{eqnarray*}
v_\mathrm{cell} = B_\mathrm{ref\_cell}^T J_\mathrm{cell}^{-1} D J_\mathrm{cell}^{-\mathrm T}
B_\mathrm{ref\_cell} u_\mathrm{cell}, \quad
v = \sum_{\mathrm{cell}=1}^{\mathrm{n\_cells}} P_\mathrm{cell,{loc-glob}}^T
v_\mathrm{cell}.


@f}



@code
  FEValues<dim> fe_values_reference (fe, quadrature_formula,
                                     update_gradients);
  Triangulation<dim> reference_cell;
  GridGenerator::hyper_cube(reference_cell, 0., 1.);
  fe_values_reference.reinit (reference_cell.begin());


  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_inverse_jacobians | update_JxW_values |
                           update_quadrature_points);


  for (const auto & cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit (cell);
      coefficient.value_list(fe_values.get_quadrature_points(),
                             coefficient_values);


      cell->get_dof_indices (local_dof_indices);


      for (unsigned int i=0; i<dofs_per_cell; ++i)
        cell_src(i) = src(local_dof_indices(i));


      temp_vector = 0;
      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            temp_vector(q*dim+d) +=
              fe_values_reference.shape_grad(i,q)[d] * cell_src(i);


      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // apply the transposed inverse Jacobian of the mapping
          Tensor<1,dim> temp;
          for (unsigned int d=0; d<dim; ++d)
            temp[d] = temp_vector(q*dim+d);
          for (unsigned int d=0; d<dim; ++d)
            {
              double sum = 0;
              for (unsigned int e=0; e<dim; ++e)
                sum += fe_values.inverse_jacobian(q)[e][d] *
                               temp[e];
              temp_vector(q*dim+d) = sum;
            }


          // multiply by coefficient and integration weight
          for (unsigned int d=0; d<dim; ++d)
            temp_vector(q*dim+d) *= fe_values.JxW(q) * coefficient_values[q];


          // apply the inverse Jacobian of the mapping
          for (unsigned int d=0; d<dim; ++d)
            temp[d] = temp_vector(q*dim+d);
          for (unsigned int d=0; d<dim; ++d)
            {
              double sum = 0;
              for (unsigned int e=0; e<dim; ++e)
                sum += fe_values.inverse_jacobian(q)[d][e] *
                       temp[e];
              temp_vector(q*dim+d) = sum;
            }
        }


      cell_dst = 0;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int d=0; d<dim; ++d)
            cell_dst(i) += fe_values_reference.shape_grad(i,q)[d] *
                                   temp_vector(q*dim+d);


      for (unsigned int i=0; i<dofs_per_cell; ++i)
        dst(local_dof_indices(i)) += cell_dst(i);
    }
}
@endcode



注意我们如何为参考单元梯度创建一个额外的FEValues对象，以及如何将其初始化为参考单元。然后，实际的导数数据是由反的、转置的Jacobian（deal.II将Jacobian矩阵从实单元到单位单元称为inverse_jacobian，因为正向转换是从单位单元到实单元）应用的。因子 $J_\mathrm{cell}^{-1} D_\mathrm{cell} J_\mathrm{cell}^{-\mathrm T}$ 是块对角线超过正交的。在这种形式下，人们意识到可变系数（可能通过张量表示）和一般网格拓扑结构的雅各布变换对变换单元格导数的系数有类似的影响。

在这一点上，人们可能会想，为什么我们要分别存储矩阵 $J_\mathrm{cell}^{-\mathrm T}$ 和系数，而不是只存储完整的因子 $J_\mathrm{cell}^{-1} D_\mathrm{cell}
J_\mathrm{cell}^{-\mathrm T}$  。后者会使用更少的内存，因为张量是对称的，在三维中具有六个独立的值，而对于前者，我们需要九个条目用于反转雅各布系数，一个用于正交权重和雅各布行列式，一个用于系数，总共是11个双数。原因是前者允许通过一个共同的缓存数据框架来实现通用的微分算子，而后者则专门存储拉普拉斯的系数。如果应用需要，这种专门化可能会得到回报，值得考虑。请注意，deal.II中的实现足够聪明，可以检测笛卡尔或仿生几何，其中雅各布系数在整个单元中是恒定的，不需要为每个单元存储（实际上在不同的单元中也常常是相同的）。

从操作数的角度来看，最后的优化是利用基函数中的张量积结构，这是最为关键的。这是可能的，因为我们已经从<i>B</i><sub>ref_cell</sub>描述的参考单元操作中剔除了梯度，即对参考单元的完全规则的数据域进行插值操作。我们举例说明在两个空间维度上降低复杂度的过程，但是同样的技术也可以用在更高的维度上。在参考单元上，基函数是张量积形式的  $\phi(x,y,z) = \varphi_i(x) \varphi_j(y)$  。矩阵<i>B</i><sub>ref_cell</sub>计算第一分量的部分具有 $B_\mathrm{sub\_cell}^x = B_\mathrm{grad,x} \otimes B_\mathrm{val,y}$ 的形式，其中<i>B</i><sub>grad,x</sub>和<i>B</i><sub>val,y</sub>包含所有一维正交点上所有一维基函数的评价。用含有属于基函数 $\varphi_i(x) \varphi_j(y)$ 的系数的<i>U</i>组成矩阵<i>U(j,i)</i>，我们得到 $(B_\mathrm{grad,x} \otimes
B_\mathrm{val,y})u_\mathrm{cell} = B_\mathrm{val,y} U B_\mathrm{grad,x}$ 。这就把计算这个乘积的复杂度从 $p^4$ 降低到 $2 p^3$ ，其中<i>p</i>-1是有限元的度数（即，等价地，<i>p</i>是每个坐标方向上的形状函数的数量），或者一般来说 $p^{2d}$ 到 $d p^{d+1}$ 。我们之所以用多项式度数来看复杂度，是因为我们希望能够到高度数，可能会增加多项式度数<i>p</i>而不是网格分辨率。像这里使用的中等度数的好算法是独立于维度的多项式度数的线性算法，而不是基于矩阵的方案或通过FEValues的天真评价。在deal.II的实现中所使用的技术自20世纪80年代以来就已经在谱元界建立起来。

实现一个无矩阵和基于单元的有限元算子，与以前的教程程序中显示的通常的矩阵装配代码相比，需要一个有点不同的程序设计。做到这一点的数据结构是MatrixFree类和FEEvaluation类，前者收集所有数据并在所有单元上发出一个（并行）循环，后者利用张量积结构评估有限元基函数。

本教程中展示的无矩阵的矩阵-向量乘积的实现比使用稀疏矩阵的线性元素的矩阵-向量乘积要慢，但由于张量乘积结构降低了复杂度，并且在计算过程中减少了内存传输，所以对所有高阶元素来说速度更快。当在一个多核处理器上工作时，减少内存传输的影响特别有利，因为在这个处理器上有几个处理单元共享内存的访问。在这种情况下，一个受计算约束的算法将显示出几乎完美的并行加速（除了可能通过涡轮模式改变处理器的时钟频率，这取决于有多少个核心在工作），而一个受内存传输约束的算法可能无法实现类似的加速（即使工作是完全并行的，我们可以期待像稀疏矩阵-向量产品那样的完美缩放）。这种实现方式的另一个好处是，我们不必建立稀疏矩阵本身，这也可能是相当昂贵的，这取决于基础微分方程。此外，上述框架可以简单地推广到非线性运算，正如我们在步骤48中所展示的那样。




<a name="Combinationwithmultigrid"></a><h3>Combination with multigrid</h3>


上面，我们花了很大的力气来实现一个不实际存储矩阵元素的矩阵-向量积。然而，在许多用户代码中，人们想要的不仅仅是做一些矩阵-向量乘积&mdash；在求解线性系统时，人们希望尽可能少做这些操作。理论上，我们可以使用CG方法，而不需要预处理；然而，这对拉普拉斯的效率并不高。相反，预调节器是用来提高收敛速度的。不幸的是，大多数比较常用的预处理方法，如SSOR、ILU或代数多网格（AMG）不能在这里使用，因为它们的实现需要了解系统矩阵的元素。

一个解决方案是使用几何多网格方法，如步骤16所示。众所周知，它们的速度非常快，而且适合我们的目的，因为所有的成分，包括不同网格层之间的转移，都可以用与单元格集合相关的矩阵-向量产品来表示。我们需要做的就是找到一个基于矩阵-向量乘积而不是所有矩阵条目的平滑器。一个这样的候选方法是阻尼雅可比迭代，它需要访问矩阵对角线，但它在阻尼所有高频误差方面往往不够好。雅可比方法的特性可以通过所谓的切比雪夫迭代进行几次改进。切比雪夫迭代由矩阵-向量乘积的多项式表达式描述，其中的系数可以被选择来实现某些特性，在这种情况下，可以平滑误差的高频成分，这些误差与雅可比预处理矩阵的特征值相关。在零度时，具有最佳阻尼参数的雅可比方法被检索出来，而高阶修正被用来改善平滑特性。切比雪夫平滑法在多网格中的有效性已经被证明，例如在文章<a href="http://www.sciencedirect.com/science/article/pii/S0021999103001943">
<i>M. Adams, M. Brezina, J. Hu, R. Tuminaro. Parallel multigrid smoothers:
polynomial versus Gauss&ndash;Seidel, J. Comput. Phys. 188:593&ndash;610,
2003</i><i>M. Adams, M. Brezina, J. Hu, R. Tuminaro. Parallel multigrid smoothers:
polynomial versus Gauss&ndash;Seidel, J. Comput. Phys. 188:593&ndash;610,
2003</i></a>中。这篇文章还指出了我们在这里利用的切比雪夫平滑器的另一个优势，即它们很容易并行化，而SOR/Gauss&ndash;Seidel平滑依赖于替换，对于这种替换，天真的并行化在矩阵的对角线子块上工作，从而降低了效率（更多细节见例如Y. Saad, Iterative Methods for Sparse Linear Systems, SIAM, 2nd edition, 2003, chapters 11 & 12）。

然后，在多网格框架中的实现就很简单了。本程序中的多网格实现与step-16类似，包括自适应性。




<a name="UsingCPUdependentinstructionsvectorization"></a><h3>Using CPU-dependent instructions (vectorization)</h3>


FEEvaluation中的计算内核是以优化使用计算资源的方式来编写的。为了达到这个目的，他们不对双倍数据类型进行操作，而是对我们称之为VectorizedArray的东西进行操作（例如，查看 FEEvaluationBase::get_value, 的返回类型，对于标量元素是VectorizedArray，对于矢量有限元素是Tensor of VectorizedArray）。VectorizedArray是一个双数或浮点数的短阵列，其长度取决于使用的特定计算机系统。例如，基于x86-64的系统支持流式SIMD扩展（SSE），处理器的矢量单元可以通过一条CPU指令处理两个双数（或四个单精度浮点数）。较新的处理器（大约从2012年起）支持所谓的高级向量扩展（AVX），有256位操作数，可以分别使用四个双数和八个浮点数。矢量化是一个单指令/多数据（SIMD）的概念，也就是说，一条CPU指令被用来同时处理多个数据值。通常情况下，有限元程序不会明确使用矢量化，因为这个概念的好处只体现在算术密集型操作中。大部分典型的有限元工作负载都受到内存带宽的限制（对稀疏矩阵和向量的操作），在这种情况下，额外的计算能力是无用的。

不过，在幕后，优化的BLAS包可能严重依赖矢量化。另外，优化的编译器可能会自动将涉及标准代码的循环转化为更有效的矢量化形式（deal.II在矢量更新的常规循环中使用OpenMP SIMD pragmas）。然而，数据流必须非常有规律，才能让编译器产生高效的代码。例如，受益于矢量化的原型操作（矩阵-矩阵乘积）的自动矢量化在大多数编译器上都失败了（截至2012年初编写本教程并在2016年底更新时，gcc和英特尔编译器都无法为 FullMatrix::mmult 函数，甚至在更简单的情况下也不行，即矩阵边界是编译时常量而不是 FullMatrix::mmult). 中的运行时常量。此外，有可能被一起处理的数据可能没有以连续的方式布置在内存中，或者没有对处理器需要的地址边界进行必要的对齐。或者编译器可能无法证明数据阵列在一次加载几个元素时不会重叠。

因此，在deal.II的无矩阵实现中，我们选择在最适合于有限元计算的层次上应用矢量化。所有单元的计算通常是完全相同的（除了从向量读写时使用的间接寻址中的索引），因此SIMD可以用来一次处理几个单元。在下面的所有内容中，你可以考虑用一个向量数组来保存几个单元的数据。记住，它与空间维度和元素数量无关，例如在Tensor或Point中。

请注意，矢量化取决于代码运行的CPU，以及代码的编译对象。为了给你的计算机生成最快的FEEvaluation内核，你应该用所谓的<i>native</i>处理器变体编译deal.II。当使用gcc编译器时，可以通过在cmake构建设置中设置变量<tt>CMAKE_CXX_FLAGS</tt>为<tt>"-march=native"</tt>来启用它（在命令行中，指定<tt>-DCMAKE_CXX_FLAGS="-march=native"</tt>，更多信息见deal.II阅读手册）。其他编译器也有类似的选项。我们在本例的run()函数中输出当前的矢量化长度。




<a name="Runningmultigridonlargescaleparallelcomputers"></a><h3>Running multigrid on large-scale parallel computers</h3>


如上所述，无矩阵框架中的所有组件都可以通过MPI使用领域分解轻松实现并行化。由于在deal.II中通过p4est（详见step-40）可以很容易地访问大规模的并行网格，而且基于单元格的循环与无矩阵评估<i>only</i>需要在每个处理器上将网格分解成大小基本相同的块，因此编写一个使用分布式内存工作的并行程序所需的工作相对较少。虽然其他使用MPI的教程程序依赖于PETSc或Trilinos，但这个程序使用deal.II自己的并行向量设施。

deal.II并行向量类， LinearAlgebra::distributed::Vector, 持有解决方案的处理器本地部分以及重影自由度的数据字段，即由远程处理器拥有的自由度，但由当前处理器拥有的单元访问。在 @ref GlossLocallyActiveDof "术语表 "中，这些自由度被称为本地活动自由度。函数 MatrixFree::initialize_dof_vector() 提供了一个设置这种设计的方法。请注意，悬挂节点可以与额外的重影自由度有关，这些自由度必须包括在分布式矢量中，但不属于 @ref
GlossLocallyActiveDof "词汇表 "意义上的本地活动自由度。此外，分布式向量持有本地拥有但其他处理器需要的DoF的MPI元数据。这个向量类设计的一个好处是对重影项的访问方式。在向量的存储方案中，数据阵列延伸到解决方案的处理器本地部分之外，有更多的向量条目可用于重影自由度。这为所有本地活动自由度提供了一个连续的索引范围。(注意，索引范围取决于网格的具体配置。)由于无矩阵操作可以被认为是在做性能关键的线性代数，而性能关键的代码不能把时间浪费在做MPI全局到MPI局部的索引转换上，一个MPI等级的局部索引空间的可用性是很重要的。这里访问事物的方式是直接数组访问。这是通过 LinearAlgebra::distributed::Vector::local_element(), 提供的，但实际上很少需要，因为所有这些都发生在FEEvaluation的内部。

 LinearAlgebra::distributed::Vector 的设计与我们之前在step-40和step-32中使用的 PETScWrappers::MPI::Vector 和 TrilinosWrappers::MPI::Vector 数据类型类似，但由于我们不需要这些库的其他并行功能，所以我们使用deal.II的 LinearAlgebra::distributed::Vector 类来代替在本教程程序中链接另一个大型库。还要注意的是，PETSc和Trilinos向量不提供对直接数组访问的幽灵条目的细粒度控制，因为它们抽象出了必要的实现细节。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * First include the necessary files from the deal.II library.
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/function.h>
 * #include <deal.II/base/timer.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/solver_cg.h>
 * #include <deal.II/lac/la_parallel_vector.h>
 * #include <deal.II/lac/precondition.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * 
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * 
 * #include <deal.II/multigrid/multigrid.h>
 * #include <deal.II/multigrid/mg_transfer_matrix_free.h>
 * #include <deal.II/multigrid/mg_tools.h>
 * #include <deal.II/multigrid/mg_coarse.h>
 * #include <deal.II/multigrid/mg_smoother.h>
 * #include <deal.II/multigrid/mg_matrix.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * This includes the data structures for the efficient implementation of
 * matrix-free methods or more generic finite element operators with the class
 * MatrixFree.
 * 
 * @code
 * #include <deal.II/matrix_free/matrix_free.h>
 * #include <deal.II/matrix_free/operators.h>
 * #include <deal.II/matrix_free/fe_evaluation.h>
 * 
 * #include <iostream>
 * #include <fstream>
 * 
 * 
 * namespace Step37
 * {
 *   using namespace dealii;
 * 
 * 
 * @endcode
 * 
 * To be efficient, the operations performed in the matrix-free
 * implementation require knowledge of loop lengths at compile time, which
 * are given by the degree of the finite element. Hence, we collect the
 * values of the two template parameters that can be changed at one place in
 * the code. Of course, one could make the degree of the finite element a
 * run-time parameter by compiling the computational kernels for all degrees
 * that are likely (say, between 1 and 6) and selecting the appropriate
 * kernel at run time. Here, we simply choose second order $Q_2$ elements
 * and choose dimension 3 as standard.
 * 
 * @code
 *   const unsigned int degree_finite_element = 2;
 *   const unsigned int dimension             = 3;
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * We define a variable coefficient function for the Poisson problem. It is
 * similar to the function in step-5 but we use the form $a(\mathbf
 * x)=\frac{1}{0.05 + 2\|\bf x\|^2}$ instead of a discontinuous one. It is
 * merely to demonstrate the possibilities of this implementation, rather
 * than making much sense physically. We define the coefficient in the same
 * way as functions in earlier tutorial programs. There is one new function,
 * namely a @p value method with template argument @p number.
 * 
 * @code
 *   template <int dim>
 *   class Coefficient : public Function<dim>
 *   {
 *   public:
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *     template <typename number>
 *     number value(const Point<dim, number> &p,
 *                  const unsigned int        component = 0) const;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * This is the new function mentioned above: Evaluate the coefficient for
 * abstract type @p number. It might be just a usual double, but it can also
 * be a somewhat more complicated type that we call VectorizedArray. This
 * data type is essentially a short array of doubles as discussed in the
 * introduction that holds data from several cells. For example, we evaluate
 * the coefficient shown here not on a simple point as usually done, but we
 * hand it a Point<dim,VectorizedArray<double> > point, which is actually a
 * collection of four points in the case of AVX. Do not confuse the entries
 * in VectorizedArray with the different coordinates of the point. Indeed,
 * the data is laid out such that <code>p[0]</code> returns a
 * VectorizedArray, which in turn contains the x-coordinate for the first
 * point and the second point. You may access the coordinates individually
 * using e.g. <code>p[0][j]</code>, j=0,1,2,3, but it is recommended to
 * define operations on a VectorizedArray as much as possible in order to
 * make use of vectorized operations.
 *   

 * 
 * In the function implementation, we assume that the number type overloads
 * basic arithmetic operations, so we just write the code as usual. The base
 * class function @p value is then computed from the templated function with
 * double type, in order to avoid duplicating code.
 * 
 * @code
 *   template <int dim>
 *   template <typename number>
 *   number Coefficient<dim>::value(const Point<dim, number> &p,
 *                                  const unsigned int /*component*/) const
 *   {
 *     return 1. / (0.05 + 2. * p.square());
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   double Coefficient<dim>::value(const Point<dim> & p,
 *                                  const unsigned int component) const
 *   {
 *     return value<double>(p, component);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Matrixfreeimplementation"></a> 
 * <h3>Matrix-free implementation</h3>
 * 

 * 
 * The following class, called <code>LaplaceOperator</code>, implements the
 * differential operator. For all practical purposes, it is a matrix, i.e.,
 * you can ask it for its size (member functions <code>m(), n()</code>) and
 * you can apply it to a vector (the <code>vmult()</code> function). The
 * difference to a real matrix of course lies in the fact that this class
 * does not actually store the <i>elements</i> of the matrix, but only knows
 * how to compute the action of the operator when applied to a vector.
 *   

 * 
 * The infrastructure describing the matrix size, the initialization from a
 * MatrixFree object, and the various interfaces to matrix-vector products
 * through vmult() and Tvmult() methods, is provided by the class
 * MatrixFreeOperator::Base from which this class derives. The
 * LaplaceOperator class defined here only has to provide a few interfaces,
 * namely the actual action of the operator through the apply_add() method
 * that gets used in the vmult() functions, and a method to compute the
 * diagonal entries of the underlying matrix. We need the diagonal for the
 * definition of the multigrid smoother. Since we consider a problem with
 * variable coefficient, we further implement a method that can fill the
 * coefficient values.
 *   

 * 
 * Note that the file <code>include/deal.II/matrix_free/operators.h</code>
 * already contains an implementation of the Laplacian through the class
 * MatrixFreeOperators::LaplaceOperator. For educational purposes, the
 * operator is re-implemented in this tutorial program, explaining the
 * ingredients and concepts used there.
 *   

 * 
 * This program makes use of the data cache for finite element operator
 * application that is integrated in deal.II. This data cache class is
 * called MatrixFree. It contains mapping information (Jacobians) and index
 * relations between local and global degrees of freedom. It also contains
 * constraints like the ones from hanging nodes or Dirichlet boundary
 * conditions. Moreover, it can issue a loop over all cells in %parallel,
 * making sure that only cells are worked on that do not share any degree of
 * freedom (this makes the loop thread-safe when writing into destination
 * vectors). This is a more advanced strategy compared to the WorkStream
 * class described in the @ref threads module. Of course, to not destroy
 * thread-safety, we have to be careful when writing into class-global
 * structures.
 *   

 * 
 * The class implementing the Laplace operator has three template arguments,
 * one for the dimension (as many deal.II classes carry), one for the degree
 * of the finite element (which we need to enable efficient computations
 * through the FEEvaluation class), and one for the underlying scalar
 * type. We want to use <code>double</code> numbers (i.e., double precision,
 * 64-bit floating point) for the final matrix, but floats (single
 * precision, 32-bit floating point numbers) for the multigrid level
 * matrices (as that is only a preconditioner, and floats can be processed
 * twice as fast). The class FEEvaluation also takes a template argument for
 * the number of quadrature points in one dimension. In the code below, we
 * hard-code it to <code>fe_degree+1</code>. If we wanted to change it
 * independently of the polynomial degree, we would need to add a template
 * parameter as is done in the MatrixFreeOperators::LaplaceOperator class.
 *   

 * 
 * As a sidenote, if we implemented several different operations on the same
 * grid and degrees of freedom (like a mass matrix and a Laplace matrix), we
 * would define two classes like the current one for each of the operators
 * (derived from the MatrixFreeOperators::Base class), and let both of them
 * refer to the same MatrixFree data cache from the general problem
 * class. The interface through MatrixFreeOperators::Base requires us to
 * only provide a minimal set of functions. This concept allows for writing
 * complex application codes with many matrix-free operations.
 *   

 * 
 * @note Storing values of type <code>VectorizedArray<number></code>
 * requires care: Here, we use the deal.II table class which is prepared to
 * hold the data with correct alignment. However, storing e.g. an
 * <code>std::vector<VectorizedArray<number> ></code> is not possible with
 * vectorization: A certain alignment of the data with the memory address
 * boundaries is required (essentially, a VectorizedArray that is 32 bytes
 * long in case of AVX needs to start at a memory address that is divisible
 * by 32). The table class (as well as the AlignedVector class it is based
 * on) makes sure that this alignment is respected, whereas std::vector does
 * not in general, which may lead to segmentation faults at strange places
 * for some systems or suboptimal performance for other systems.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   class LaplaceOperator
 *     : public MatrixFreeOperators::
 *         Base<dim, LinearAlgebra::distributed::Vector<number>>
 *   {
 *   public:
 *     using value_type = number;
 * 
 *     LaplaceOperator();
 * 
 *     void clear() override;
 * 
 *     void evaluate_coefficient(const Coefficient<dim> &coefficient_function);
 * 
 *     virtual void compute_diagonal() override;
 * 
 *   private:
 *     virtual void apply_add(
 *       LinearAlgebra::distributed::Vector<number> &      dst,
 *       const LinearAlgebra::distributed::Vector<number> &src) const override;
 * 
 *     void
 *     local_apply(const MatrixFree<dim, number> &                   data,
 *                 LinearAlgebra::distributed::Vector<number> &      dst,
 *                 const LinearAlgebra::distributed::Vector<number> &src,
 *                 const std::pair<unsigned int, unsigned int> &cell_range) const;
 * 
 *     void local_compute_diagonal(
 *       const MatrixFree<dim, number> &              data,
 *       LinearAlgebra::distributed::Vector<number> & dst,
 *       const unsigned int &                         dummy,
 *       const std::pair<unsigned int, unsigned int> &cell_range) const;
 * 
 *     Table<2, VectorizedArray<number>> coefficient;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * This is the constructor of the @p LaplaceOperator class. All it does is
 * to call the default constructor of the base class
 * MatrixFreeOperators::Base, which in turn is based on the Subscriptor
 * class that asserts that this class is not accessed after going out of scope
 * e.g. in a preconditioner.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
 *     : MatrixFreeOperators::Base<dim,
 *                                 LinearAlgebra::distributed::Vector<number>>()
 *   {}
 * 
 * 
 * 
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::clear()
 *   {
 *     coefficient.reinit(0, 0);
 *     MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
 *       clear();
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Computationofcoefficient"></a> 
 * <h4>Computation of coefficient</h4>
 * 

 * 
 * To initialize the coefficient, we directly give it the Coefficient class
 * defined above and then select the method
 * <code>coefficient_function.value</code> with vectorized number (which the
 * compiler can deduce from the point data type). The use of the
 * FEEvaluation class (and its template arguments) will be explained below.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::evaluate_coefficient(
 *     const Coefficient<dim> &coefficient_function)
 *   {
 *     const unsigned int n_cells = this->data->n_cell_batches();
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*this->data);
 * 
 *     coefficient.reinit(n_cells, phi.n_q_points);
 *     for (unsigned int cell = 0; cell < n_cells; ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           coefficient(cell, q) =
 *             coefficient_function.value(phi.quadrature_point(q));
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LocalevaluationofLaplaceoperator"></a> 
 * <h4>Local evaluation of Laplace operator</h4>
 * 

 * 
 * Here comes the main function of this class, the evaluation of the
 * matrix-vector product (or, in general, a finite element operator
 * evaluation). This is done in a function that takes exactly four
 * arguments, the MatrixFree object, the destination and source vectors, and
 * a range of cells that are to be worked on. The method
 * <code>cell_loop</code> in the MatrixFree class will internally call this
 * function with some range of cells that is obtained by checking which
 * cells are possible to work on simultaneously so that write operations do
 * not cause any race condition. Note that the cell range used in the loop
 * is not directly the number of (active) cells in the current mesh, but
 * rather a collection of batches of cells.  In other word, "cell" may be
 * the wrong term to begin with, since FEEvaluation groups data from several
 * cells together. This means that in the loop over quadrature points we are
 * actually seeing a group of quadrature points of several cells as one
 * block. This is done to enable a higher degree of vectorization.  The
 * number of such "cells" or "cell batches" is stored in MatrixFree and can
 * be queried through MatrixFree::n_cell_batches(). Compared to the deal.II
 * cell iterators, in this class all cells are laid out in a plain array
 * with no direct knowledge of level or neighborship relations, which makes
 * it possible to index the cells by unsigned integers.
 *   

 * 
 * The implementation of the Laplace operator is quite simple: First, we
 * need to create an object FEEvaluation that contains the computational
 * kernels and has data fields to store temporary results (e.g. gradients
 * evaluated on all quadrature points on a collection of a few cells). Note
 * that temporary results do not use a lot of memory, and since we specify
 * template arguments with the element order, the data is stored on the
 * stack (without expensive memory allocation). Usually, one only needs to
 * set two template arguments, the dimension as a first argument and the
 * degree of the finite element as the second argument (this is equal to the
 * number of degrees of freedom per dimension minus one for FE_Q
 * elements). However, here we also want to be able to use float numbers for
 * the multigrid preconditioner, which is the last (fifth) template
 * argument. Therefore, we cannot rely on the default template arguments and
 * must also fill the third and fourth field, consequently. The third
 * argument specifies the number of quadrature points per direction and has
 * a default value equal to the degree of the element plus one. The fourth
 * argument sets the number of components (one can also evaluate
 * vector-valued functions in systems of PDEs, but the default is a scalar
 * element), and finally the last argument sets the number type.
 *   

 * 
 * Next, we loop over the given cell range and then we continue with the
 * actual implementation: <ol> <li>Tell the FEEvaluation object the (macro)
 * cell we want to work on.  <li>Read in the values of the source vectors
 * (@p read_dof_values), including the resolution of constraints. This
 * stores $u_\mathrm{cell}$ as described in the introduction.  <li>Compute
 * the unit-cell gradient (the evaluation of finite element
 * functions). Since FEEvaluation can combine value computations with
 * gradient computations, it uses a unified interface to all kinds of
 * derivatives of order between zero and two. We only want gradients, no
 * values and no second derivatives, so we set the function arguments to
 * true in the gradient slot (second slot), and to false in the values slot
 * (first slot). There is also a third slot for the Hessian which is
 * false by default, so it needs not be given. Note that the FEEvaluation
 * class internally evaluates shape functions in an efficient way where one
 * dimension is worked on at a time (using the tensor product form of shape
 * functions and quadrature points as mentioned in the introduction). This
 * gives complexity equal to $\mathcal O(d^2 (p+1)^{d+1})$ for polynomial
 * degree $p$ in $d$ dimensions, compared to the naive approach with loops
 * over all local degrees of freedom and quadrature points that is used in
 * FEValues and costs $\mathcal O(d (p+1)^{2d})$.  <li>Next comes the
 * application of the Jacobian transformation, the multiplication by the
 * variable coefficient and the quadrature weight. FEEvaluation has an
 * access function @p get_gradient that applies the Jacobian and returns the
 * gradient in real space. Then, we just need to multiply by the (scalar)
 * coefficient, and let the function @p submit_gradient apply the second
 * Jacobian (for the test function) and the quadrature weight and Jacobian
 * determinant (JxW). Note that the submitted gradient is stored in the same
 * data field as where it is read from in @p get_gradient. Therefore, you
 * need to make sure to not read from the same quadrature point again after
 * having called @p submit_gradient on that particular quadrature point. In
 * general, it is a good idea to copy the result of @p get_gradient when it
 * is used more often than once.  <li>Next follows the summation over
 * quadrature points for all test functions that corresponds to the actual
 * integration step. For the Laplace operator, we just multiply by the
 * gradient, so we call the integrate function with the respective argument
 * set. If you have an equation where you test by both the values of the
 * test functions and the gradients, both template arguments need to be set
 * to true. Calling first the integrate function for values and then
 * gradients in a separate call leads to wrong results, since the second
 * call will internally overwrite the results from the first call. Note that
 * there is no function argument for the second derivative for integrate
 * step.  <li>Eventually, the local contributions in the vector
 * $v_\mathrm{cell}$ as mentioned in the introduction need to be added into
 * the result vector (and constraints are applied). This is done with a call
 * to @p distribute_local_to_global, the same name as the corresponding
 * function in the AffineConstraints (only that we now store the local vector
 * in the FEEvaluation object, as are the indices between local and global
 * degrees of freedom).  </ol>
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::local_apply(
 *     const MatrixFree<dim, number> &                   data,
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src,
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const
 *   {
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         AssertDimension(coefficient.size(0), data.n_cell_batches());
 *         AssertDimension(coefficient.size(1), phi.n_q_points);
 * 
 *         phi.reinit(cell);
 *         phi.read_dof_values(src);
 *         phi.evaluate(EvaluationFlags::gradients);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q), q);
 *         phi.integrate(EvaluationFlags::gradients);
 *         phi.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * This function implements the loop over all cells for the
 * Base::apply_add() interface. This is done with the @p cell_loop of the
 * MatrixFree class, which takes the operator() of this class with arguments
 * MatrixFree, OutVector, InVector, cell_range. When working with MPI
 * parallelization (but no threading) as is done in this tutorial program,
 * the cell loop corresponds to the following three lines of code:
 *   

 * 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 * src.update_ghost_values();
 * local_apply(*this->data, dst, src, std::make_pair(0U,
 *                                                   data.n_cell_batches()));
 * dst.compress(VectorOperation::add);
 * @endcode
 * </div>
 *   

 * 
 * Here, the two calls update_ghost_values() and compress() perform the data
 * exchange on processor boundaries for MPI, once for the source vector
 * where we need to read from entries owned by remote processors, and once
 * for the destination vector where we have accumulated parts of the
 * residuals that need to be added to the respective entry of the owner
 * processor. However, MatrixFree::cell_loop does not only abstract away
 * those two calls, but also performs some additional optimizations. On the
 * one hand, it will split the update_ghost_values() and compress() calls in
 * a way to allow for overlapping communication and computation. The
 * local_apply function is then called with three cell ranges representing
 * partitions of the cell range from 0 to MatrixFree::n_cell_batches(). On
 * the other hand, cell_loop also supports thread parallelism in which case
 * the cell ranges are split into smaller chunks and scheduled in an
 * advanced way that avoids access to the same vector entry by several
 * threads. That feature is explained in step-48.
 *   

 * 
 * Note that after the cell loop, the constrained degrees of freedom need to
 * be touched once more for sensible vmult() operators: Since the assembly
 * loop automatically resolves constraints (just as the
 * AffineConstraints::distribute_local_to_global() call does), it does not
 * compute any contribution for constrained degrees of freedom, leaving the
 * respective entries zero. This would represent a matrix that had empty
 * rows and columns for constrained degrees of freedom. However, iterative
 * solvers like CG only work for non-singular matrices. The easiest way to
 * do that is to set the sub-block of the matrix that corresponds to
 * constrained DoFs to an identity matrix, in which case application of the
 * matrix would simply copy the elements of the right hand side vector into
 * the left hand side. Fortunately, the vmult() implementations
 * MatrixFreeOperators::Base do this automatically for us outside the
 * apply_add() function, so we do not need to take further action here.
 *   

 * 
 * When using the combination of MatrixFree and FEEvaluation in parallel
 * with MPI, there is one aspect to be careful about &mdash; the indexing
 * used for accessing the vector. For performance reasons, MatrixFree and
 * FEEvaluation are designed to access vectors in MPI-local index space also
 * when working with multiple processors. Working in local index space means
 * that no index translation needs to be performed at the place the vector
 * access happens, apart from the unavoidable indirect addressing. However,
 * local index spaces are ambiguous: While it is standard convention to
 * access the locally owned range of a vector with indices between 0 and the
 * local size, the numbering is not so clear for the ghosted entries and
 * somewhat arbitrary. For the matrix-vector product, only the indices
 * appearing on locally owned cells (plus those referenced via hanging node
 * constraints) are necessary. However, in deal.II we often set all the
 * degrees of freedom on ghosted elements as ghosted vector entries, called
 * the @ref GlossLocallyRelevantDof "locally relevant DoFs described in the
 * glossary". In that case, the MPI-local index of a ghosted vector entry
 * can in general be different in the two possible ghost sets, despite
 * referring to the same global index. To avoid problems, FEEvaluation
 * checks that the partitioning of the vector used for the matrix-vector
 * product does indeed match with the partitioning of the indices in
 * MatrixFree by a check called
 * LinearAlgebra::distributed::Vector::partitioners_are_compatible. To
 * facilitate things, the MatrixFreeOperators::Base class includes a
 * mechanism to fit the ghost set to the correct layout. This happens in the
 * ghost region of the vector, so keep in mind that the ghost region might
 * be modified in both the destination and source vector after a call to a
 * vmult() method. This is legitimate because the ghost region of a
 * distributed deal.II vector is a mutable section and filled on
 * demand. Vectors used in matrix-vector products must not be ghosted upon
 * entry of vmult() functions, so no information gets lost.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::apply_add(
 *     LinearAlgebra::distributed::Vector<number> &      dst,
 *     const LinearAlgebra::distributed::Vector<number> &src) const
 *   {
 *     this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * The following function implements the computation of the diagonal of the
 * operator. Computing matrix entries of a matrix-free operator evaluation
 * turns out to be more complicated than evaluating the
 * operator. Fundamentally, we could obtain a matrix representation of the
 * operator by applying the operator on <i>all</i> unit vectors. Of course,
 * that would be very inefficient since we would need to perform <i>n</i>
 * operator evaluations to retrieve the whole matrix. Furthermore, this
 * approach would completely ignore the matrix sparsity. On an individual
 * cell, however, this is the way to go and actually not that inefficient as
 * there usually is a coupling between all degrees of freedom inside the
 * cell.
 *   

 * 
 * We first initialize the diagonal vector to the correct parallel
 * layout. This vector is encapsulated in a member called
 * inverse_diagonal_entries of type DiagonalMatrix in the base class
 * MatrixFreeOperators::Base. This member is a shared pointer that we first
 * need to initialize and then get the vector representing the diagonal
 * entries in the matrix. As to the actual diagonal computation, we again
 * use the cell_loop infrastructure of MatrixFree to invoke a local worker
 * routine called local_compute_diagonal(). Since we will only write into a
 * vector but not have any source vector, we put a dummy argument of type
 * <tt>unsigned int</tt> in place of the source vector to confirm with the
 * cell_loop interface. After the loop, we need to set the vector entries
 * subject to Dirichlet boundary conditions to one (either those on the
 * boundary described by the AffineConstraints object inside MatrixFree or
 * the indices at the interface between different grid levels in adaptive
 * multigrid). This is done through the function
 * MatrixFreeOperators::Base::set_constrained_entries_to_one() and matches
 * with the setting in the matrix-vector product provided by the Base
 * operator. Finally, we need to invert the diagonal entries which is the
 * form required by the Chebyshev smoother based on the Jacobi iteration. In
 * the loop, we assert that all entries are non-zero, because they should
 * either have obtained a positive contribution from integrals or be
 * constrained and treated by @p set_constrained_entries_to_one() following
 * cell_loop.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
 *   {
 *     this->inverse_diagonal_entries.reset(
 *       new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
 *     LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
 *       this->inverse_diagonal_entries->get_vector();
 *     this->data->initialize_dof_vector(inverse_diagonal);
 *     unsigned int dummy = 0;
 *     this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
 *                           this,
 *                           inverse_diagonal,
 *                           dummy);
 * 
 *     this->set_constrained_entries_to_one(inverse_diagonal);
 * 
 *     for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
 *       {
 *         Assert(inverse_diagonal.local_element(i) > 0.,
 *                ExcMessage("No diagonal entry in a positive definite operator "
 *                           "should be zero"));
 *         inverse_diagonal.local_element(i) =
 *           1. / inverse_diagonal.local_element(i);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * In the local compute loop, we compute the diagonal by a loop over all
 * columns in the local matrix and putting the entry 1 in the <i>i</i>th
 * slot and a zero entry in all other slots, i.e., we apply the cell-wise
 * differential operator on one unit vector at a time. The inner part
 * invoking FEEvaluation::evaluate, the loop over quadrature points, and
 * FEEvalution::integrate, is exactly the same as in the local_apply
 * function. Afterwards, we pick out the <i>i</i>th entry of the local
 * result and put it to a temporary storage (as we overwrite all entries in
 * the array behind FEEvaluation::get_dof_value() with the next loop
 * iteration). Finally, the temporary storage is written to the destination
 * vector. Note how we use FEEvaluation::get_dof_value() and
 * FEEvaluation::submit_dof_value() to read and write to the data field that
 * FEEvaluation uses for the integration on the one hand and writes into the
 * global vector on the other hand.
 *   

 * 
 * Given that we are only interested in the matrix diagonal, we simply throw
 * away all other entries of the local matrix that have been computed along
 * the way. While it might seem wasteful to compute the complete cell matrix
 * and then throw away everything but the diagonal, the integration are so
 * efficient that the computation does not take too much time. Note that the
 * complexity of operator evaluation per element is $\mathcal
 * O((p+1)^{d+1})$ for polynomial degree $k$, so computing the whole matrix
 * costs us $\mathcal O((p+1)^{2d+1})$ operations, not too far away from
 * $\mathcal O((p+1)^{2d})$ complexity for computing the diagonal with
 * FEValues. Since FEEvaluation is also considerably faster due to
 * vectorization and other optimizations, the diagonal computation with this
 * function is actually the fastest (simple) variant. (It would be possible
 * to compute the diagonal with sum factorization techniques in $\mathcal
 * O((p+1)^{d+1})$ operations involving specifically adapted
 * kernels&mdash;but since such kernels are only useful in that particular
 * context and the diagonal computation is typically not on the critical
 * path, they have not been implemented in deal.II.)
 *   

 * 
 * Note that the code that calls distribute_local_to_global on the vector to
 * accumulate the diagonal entries into the global matrix has some
 * limitations. For operators with hanging node constraints that distribute
 * an integral contribution of a constrained DoF to several other entries
 * inside the distribute_local_to_global call, the vector interface used
 * here does not exactly compute the diagonal entries, but lumps some
 * contributions located on the diagonal of the local matrix that would end
 * up in a off-diagonal position of the global matrix to the diagonal. The
 * result is correct up to discretization accuracy as explained in <a
 * href="http://dx.doi.org/10.4208/cicp.101214.021015a">Kormann (2016),
 * section 5.3</a>, but not mathematically equal. In this tutorial program,
 * no harm can happen because the diagonal is only used for the multigrid
 * level matrices where no hanging node constraints appear.
 * 
 * @code
 *   template <int dim, int fe_degree, typename number>
 *   void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
 *     const MatrixFree<dim, number> &             data,
 *     LinearAlgebra::distributed::Vector<number> &dst,
 *     const unsigned int &,
 *     const std::pair<unsigned int, unsigned int> &cell_range) const
 *   {
 *     FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);
 * 
 *     AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
 *       {
 *         AssertDimension(coefficient.size(0), data.n_cell_batches());
 *         AssertDimension(coefficient.size(1), phi.n_q_points);
 * 
 *         phi.reinit(cell);
 *         for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
 *           {
 *             for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
 *               phi.submit_dof_value(VectorizedArray<number>(), j);
 *             phi.submit_dof_value(make_vectorized_array<number>(1.), i);
 * 
 *             phi.evaluate(EvaluationFlags::gradients);
 *             for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *               phi.submit_gradient(coefficient(cell, q) * phi.get_gradient(q),
 *                                   q);
 *             phi.integrate(EvaluationFlags::gradients);
 *             diagonal[i] = phi.get_dof_value(i);
 *           }
 *         for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
 *           phi.submit_dof_value(diagonal[i], i);
 *         phi.distribute_local_to_global(dst);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemclass"></a> 
 * <h3>LaplaceProblem class</h3>
 * 

 * 
 * This class is based on the one in step-16. However, we replaced the
 * SparseMatrix<double> class by our matrix-free implementation, which means
 * that we can also skip the sparsity patterns. Notice that we define the
 * LaplaceOperator class with the degree of finite element as template
 * argument (the value is defined at the top of the file), and that we use
 * float numbers for the multigrid level matrices.
 *   

 * 
 * The class also has a member variable to keep track of all the detailed
 * timings for setting up the entire chain of data before we actually go
 * about solving the problem. In addition, there is an output stream (that
 * is disabled by default) that can be used to output details for the
 * individual setup operations instead of the summary only that is printed
 * out by default.
 *   

 * 
 * Since this program is designed to be used with MPI, we also provide the
 * usual @p pcout output stream that only prints the information of the
 * processor with MPI rank 0. The grid used for this programs can either be
 * a distributed triangulation based on p4est (in case deal.II is configured
 * to use p4est), otherwise it is a serial grid that only runs without MPI.
 * 
 * @code
 *   template <int dim>
 *   class LaplaceProblem
 *   {
 *   public:
 *     LaplaceProblem();
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_rhs();
 *     void solve();
 *     void output_results(const unsigned int cycle) const;
 * 
 * #ifdef DEAL_II_WITH_P4EST
 *     parallel::distributed::Triangulation<dim> triangulation;
 * #else
 *     Triangulation<dim> triangulation;
 * #endif
 * 
 *     FE_Q<dim>       fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     MappingQ1<dim> mapping;
 * 
 *     AffineConstraints<double> constraints;
 *     using SystemMatrixType =
 *       LaplaceOperator<dim, degree_finite_element, double>;
 *     SystemMatrixType system_matrix;
 * 
 *     MGConstrainedDoFs mg_constrained_dofs;
 *     using LevelMatrixType = LaplaceOperator<dim, degree_finite_element, float>;
 *     MGLevelObject<LevelMatrixType> mg_matrices;
 * 
 *     LinearAlgebra::distributed::Vector<double> solution;
 *     LinearAlgebra::distributed::Vector<double> system_rhs;
 * 
 *     double             setup_time;
 *     ConditionalOStream pcout;
 *     ConditionalOStream time_details;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * When we initialize the finite element, we of course have to use the
 * degree specified at the top of the file as well (otherwise, an exception
 * will be thrown at some point, since the computational kernel defined in
 * the templated LaplaceOperator class and the information from the finite
 * element read out by MatrixFree will not match). The constructor of the
 * triangulation needs to set an additional flag that tells the grid to
 * conform to the 2:1 cell balance over vertices, which is needed for the
 * convergence of the geometric multigrid routines. For the distributed
 * grid, we also need to specifically enable the multigrid hierarchy.
 * 
 * @code
 *   template <int dim>
 *   LaplaceProblem<dim>::LaplaceProblem()
 *     :
 * #ifdef DEAL_II_WITH_P4EST
 *     triangulation(
 *       MPI_COMM_WORLD,
 *       Triangulation<dim>::limit_level_difference_at_vertices,
 *       parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
 *     ,
 * #else
 *     triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
 *     ,
 * #endif
 *     fe(degree_finite_element)
 *     , dof_handler(triangulation)
 *     , setup_time(0.)
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *     ,
 * @endcode
 * 
 * The LaplaceProblem class holds an additional output stream that
 * collects detailed timings about the setup phase. This stream, called
 * time_details, is disabled by default through the @p false argument
 * specified here. For detailed timings, removing the @p false argument
 * prints all the details.
 * 
 * @code
 *     time_details(std::cout,
 *                  false && Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * The setup stage is in analogy to step-16 with relevant changes due to the
 * LaplaceOperator class. The first thing to do is to set up the DoFHandler,
 * including the degrees of freedom for the multigrid levels, and to
 * initialize constraints from hanging nodes and homogeneous Dirichlet
 * conditions. Since we intend to use this programs in %parallel with MPI,
 * we need to make sure that the constraints get to know the locally
 * relevant degrees of freedom, otherwise the storage would explode when
 * using more than a few hundred millions of degrees of freedom, see
 * step-40.
 * 

 * 
 * Once we have created the multigrid dof_handler and the constraints, we
 * can call the reinit function for the global matrix operator as well as
 * each level of the multigrid scheme. The main action is to set up the
 * <code> MatrixFree </code> instance for the problem. The base class of the
 * <code>LaplaceOperator</code> class, MatrixFreeOperators::Base, is
 * initialized with a shared pointer to MatrixFree object. This way, we can
 * simply create it here and then pass it on to the system matrix and level
 * matrices, respectively. For setting up MatrixFree, we need to activate
 * the update flag in the AdditionalData field of MatrixFree that enables
 * the storage of quadrature point coordinates in real space (by default, it
 * only caches data for gradients (inverse transposed Jacobians) and JxW
 * values). Note that if we call the reinit function without specifying the
 * level (i.e., giving <code>level = numbers::invalid_unsigned_int</code>),
 * MatrixFree constructs a loop over the active cells. In this tutorial, we
 * do not use threads in addition to MPI, which is why we explicitly disable
 * it by setting the MatrixFree::AdditionalData::tasks_parallel_scheme to
 * MatrixFree::AdditionalData::none. Finally, the coefficient is evaluated
 * and vectors are initialized as explained above.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::setup_system()
 *   {
 *     Timer time;
 *     setup_time = 0;
 * 
 *     system_matrix.clear();
 *     mg_matrices.clear_elements();
 * 
 *     dof_handler.distribute_dofs(fe);
 *     dof_handler.distribute_mg_dofs();
 * 
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
 *           << std::endl;
 * 
 *     IndexSet locally_relevant_dofs;
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 *     VectorTools::interpolate_boundary_values(
 *       mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
 *     constraints.close();
 *     setup_time += time.wall_time();
 *     time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *     time.restart();
 * 
 *     {
 *       typename MatrixFree<dim, double>::AdditionalData additional_data;
 *       additional_data.tasks_parallel_scheme =
 *         MatrixFree<dim, double>::AdditionalData::none;
 *       additional_data.mapping_update_flags =
 *         (update_gradients | update_JxW_values | update_quadrature_points);
 *       std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
 *         new MatrixFree<dim, double>());
 *       system_mf_storage->reinit(mapping,
 *                                 dof_handler,
 *                                 constraints,
 *                                 QGauss<1>(fe.degree + 1),
 *                                 additional_data);
 *       system_matrix.initialize(system_mf_storage);
 *     }
 * 
 *     system_matrix.evaluate_coefficient(Coefficient<dim>());
 * 
 *     system_matrix.initialize_dof_vector(solution);
 *     system_matrix.initialize_dof_vector(system_rhs);
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *     time.restart();
 * 
 * @endcode
 * 
 * Next, initialize the matrices for the multigrid method on all the
 * levels. The data structure MGConstrainedDoFs keeps information about
 * the indices subject to boundary conditions as well as the indices on
 * edges between different refinement levels as described in the step-16
 * tutorial program. We then go through the levels of the mesh and
 * construct the constraints and matrices on each level. These follow
 * closely the construction of the system matrix on the original mesh,
 * except the slight difference in naming when accessing information on
 * the levels rather than the active cells.
 * 
 * @code
 *     const unsigned int nlevels = triangulation.n_global_levels();
 *     mg_matrices.resize(0, nlevels - 1);
 * 
 *     std::set<types::boundary_id> dirichlet_boundary;
 *     dirichlet_boundary.insert(0);
 *     mg_constrained_dofs.initialize(dof_handler);
 *     mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
 *                                                        dirichlet_boundary);
 * 
 *     for (unsigned int level = 0; level < nlevels; ++level)
 *       {
 *         IndexSet relevant_dofs;
 *         DoFTools::extract_locally_relevant_level_dofs(dof_handler,
 *                                                       level,
 *                                                       relevant_dofs);
 *         AffineConstraints<double> level_constraints;
 *         level_constraints.reinit(relevant_dofs);
 *         level_constraints.add_lines(
 *           mg_constrained_dofs.get_boundary_indices(level));
 *         level_constraints.close();
 * 
 *         typename MatrixFree<dim, float>::AdditionalData additional_data;
 *         additional_data.tasks_parallel_scheme =
 *           MatrixFree<dim, float>::AdditionalData::none;
 *         additional_data.mapping_update_flags =
 *           (update_gradients | update_JxW_values | update_quadrature_points);
 *         additional_data.mg_level = level;
 *         std::shared_ptr<MatrixFree<dim, float>> mg_mf_storage_level(
 *           new MatrixFree<dim, float>());
 *         mg_mf_storage_level->reinit(mapping,
 *                                     dof_handler,
 *                                     level_constraints,
 *                                     QGauss<1>(fe.degree + 1),
 *                                     additional_data);
 * 
 *         mg_matrices[level].initialize(mg_mf_storage_level,
 *                                       mg_constrained_dofs,
 *                                       level);
 *         mg_matrices[level].evaluate_coefficient(Coefficient<dim>());
 *       }
 *     setup_time += time.wall_time();
 *     time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_rhs"></a> 
 * <h4>LaplaceProblem::assemble_rhs</h4>
 * 

 * 
 * The assemble function is very simple since all we have to do is to
 * assemble the right hand side. Thanks to FEEvaluation and all the data
 * cached in the MatrixFree class, which we query from
 * MatrixFreeOperators::Base, this can be done in a few lines. Since this
 * call is not wrapped into a MatrixFree::cell_loop (which would be an
 * alternative), we must not forget to call compress() at the end of the
 * assembly to send all the contributions of the right hand side to the
 * owner of the respective degree of freedom.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::assemble_rhs()
 *   {
 *     Timer time;
 * 
 *     system_rhs = 0;
 *     FEEvaluation<dim, degree_finite_element> phi(
 *       *system_matrix.get_matrix_free());
 *     for (unsigned int cell = 0;
 *          cell < system_matrix.get_matrix_free()->n_cell_batches();
 *          ++cell)
 *       {
 *         phi.reinit(cell);
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q)
 *           phi.submit_value(make_vectorized_array<double>(1.0), q);
 *         phi.integrate(EvaluationFlags::values);
 *         phi.distribute_local_to_global(system_rhs);
 *       }
 *     system_rhs.compress(VectorOperation::add);
 * 
 *     setup_time += time.wall_time();
 *     time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s" << std::endl;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * The solution process is similar as in step-16. We start with the setup of
 * the transfer. For LinearAlgebra::distributed::Vector, there is a very
 * fast transfer class called MGTransferMatrixFree that does the
 * interpolation between the grid levels with the same fast sum
 * factorization kernels that get also used in FEEvaluation.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::solve()
 *   {
 *     Timer                            time;
 *     MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs);
 *     mg_transfer.build(dof_handler);
 *     setup_time += time.wall_time();
 *     time_details << "MG build transfer time     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *     time.restart();
 * 
 * @endcode
 * 
 * As a smoother, this tutorial program uses a Chebyshev iteration instead
 * of SOR in step-16. (SOR would be very difficult to implement because we
 * do not have the matrix elements available explicitly, and it is
 * difficult to make it work efficiently in %parallel.)  The smoother is
 * initialized with our level matrices and the mandatory additional data
 * for the Chebyshev smoother. We use a relatively high degree here (5),
 * since matrix-vector products are comparably cheap. We choose to smooth
 * out a range of $[1.2 \hat{\lambda}_{\max}/15,1.2 \hat{\lambda}_{\max}]$
 * in the smoother where $\hat{\lambda}_{\max}$ is an estimate of the
 * largest eigenvalue (the factor 1.2 is applied inside
 * PreconditionChebyshev). In order to compute that eigenvalue, the
 * Chebyshev initialization performs a few steps of a CG algorithm
 * without preconditioner. Since the highest eigenvalue is usually the
 * easiest one to find and a rough estimate is enough, we choose 10
 * iterations. Finally, we also set the inner preconditioner type in the
 * Chebyshev method which is a Jacobi iteration. This is represented by
 * the DiagonalMatrix class that gets the inverse diagonal entry provided
 * by our LaplaceOperator class.
 *     

 * 
 * On level zero, we initialize the smoother differently because we want
 * to use the Chebyshev iteration as a solver. PreconditionChebyshev
 * allows the user to switch to solver mode where the number of iterations
 * is internally chosen to the correct value. In the additional data
 * object, this setting is activated by choosing the polynomial degree to
 * @p numbers::invalid_unsigned_int. The algorithm will then attack all
 * eigenvalues between the smallest and largest one in the coarse level
 * matrix. The number of steps in the Chebyshev smoother are chosen such
 * that the Chebyshev convergence estimates guarantee to reduce the
 * residual by the number specified in the variable @p
 * smoothing_range. Note that for solving, @p smoothing_range is a
 * relative tolerance and chosen smaller than one, in this case, we select
 * three orders of magnitude, whereas it is a number larger than 1 when
 * only selected eigenvalues are smoothed.
 *     

 * 
 * From a computational point of view, the Chebyshev iteration is a very
 * attractive coarse grid solver as long as the coarse size is
 * moderate. This is because the Chebyshev method performs only
 * matrix-vector products and vector updates, which typically parallelize
 * better to the largest cluster size with more than a few tens of
 * thousands of cores than inner product involved in other iterative
 * methods. The former involves only local communication between neighbors
 * in the (coarse) mesh, whereas the latter requires global communication
 * over all processors.
 * 
 * @code
 *     using SmootherType =
 *       PreconditionChebyshev<LevelMatrixType,
 *                             LinearAlgebra::distributed::Vector<float>>;
 *     mg::SmootherRelaxation<SmootherType,
 *                            LinearAlgebra::distributed::Vector<float>>
 *                                                          mg_smoother;
 *     MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
 *     smoother_data.resize(0, triangulation.n_global_levels() - 1);
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       {
 *         if (level > 0)
 *           {
 *             smoother_data[level].smoothing_range     = 15.;
 *             smoother_data[level].degree              = 5;
 *             smoother_data[level].eig_cg_n_iterations = 10;
 *           }
 *         else
 *           {
 *             smoother_data[0].smoothing_range = 1e-3;
 *             smoother_data[0].degree          = numbers::invalid_unsigned_int;
 *             smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
 *           }
 *         mg_matrices[level].compute_diagonal();
 *         smoother_data[level].preconditioner =
 *           mg_matrices[level].get_matrix_diagonal_inverse();
 *       }
 *     mg_smoother.initialize(mg_matrices, smoother_data);
 * 
 *     MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>>
 *       mg_coarse;
 *     mg_coarse.initialize(mg_smoother);
 * 
 * @endcode
 * 
 * The next step is to set up the interface matrices that are needed for the
 * case with hanging nodes. The adaptive multigrid realization in deal.II
 * implements an approach called local smoothing. This means that the
 * smoothing on the finest level only covers the local part of the mesh
 * defined by the fixed (finest) grid level and ignores parts of the
 * computational domain where the terminal cells are coarser than this
 * level. As the method progresses to coarser levels, more and more of the
 * global mesh will be covered. At some coarser level, the whole mesh will
 * be covered. Since all level matrices in the multigrid method cover a
 * single level in the mesh, no hanging nodes appear on the level matrices.
 * At the interface between multigrid levels, homogeneous Dirichlet boundary
 * conditions are set while smoothing. When the residual is transferred to
 * the next coarser level, however, the coupling over the multigrid
 * interface needs to be taken into account. This is done by the so-called
 * interface (or edge) matrices that compute the part of the residual that
 * is missed by the level matrix with
 * homogeneous Dirichlet conditions. We refer to the @ref mg_paper
 * "Multigrid paper by Janssen and Kanschat" for more details.
 *     

 * 
 * For the implementation of those interface matrices, there is already a
 * pre-defined class MatrixFreeOperators::MGInterfaceOperator that wraps
 * the routines MatrixFreeOperators::Base::vmult_interface_down() and
 * MatrixFreeOperators::Base::vmult_interface_up() in a new class with @p
 * vmult() and @p Tvmult() operations (that were originally written for
 * matrices, hence expecting those names). Note that vmult_interface_down
 * is used during the restriction phase of the multigrid V-cycle, whereas
 * vmult_interface_up is used during the prolongation phase.
 *     

 * 
 * Once the interface matrix is created, we set up the remaining Multigrid
 * preconditioner infrastructure in complete analogy to step-16 to obtain
 * a @p preconditioner object that can be applied to a matrix.
 * 
 * @code
 *     mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix(
 *       mg_matrices);
 * 
 *     MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
 *       mg_interface_matrices;
 *     mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
 *     for (unsigned int level = 0; level < triangulation.n_global_levels();
 *          ++level)
 *       mg_interface_matrices[level].initialize(mg_matrices[level]);
 *     mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface(
 *       mg_interface_matrices);
 * 
 *     Multigrid<LinearAlgebra::distributed::Vector<float>> mg(
 *       mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
 *     mg.set_edge_matrices(mg_interface, mg_interface);
 * 
 *     PreconditionMG<dim,
 *                    LinearAlgebra::distributed::Vector<float>,
 *                    MGTransferMatrixFree<dim, float>>
 *       preconditioner(dof_handler, mg, mg_transfer);
 * 
 * @endcode
 * 
 * The setup of the multigrid routines is quite easy and one cannot see
 * any difference in the solve process as compared to step-16. All the
 * magic is hidden behind the implementation of the LaplaceOperator::vmult
 * operation. Note that we print out the solve time and the accumulated
 * setup time through standard out, i.e., in any case, whereas detailed
 * times for the setup operations are only printed in case the flag for
 * detail_times in the constructor is changed.
 * 

 * 
 * 
 * @code
 *     SolverControl solver_control(100, 1e-12 * system_rhs.l2_norm());
 *     SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control);
 *     setup_time += time.wall_time();
 *     time_details << "MG build smoother time     (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *     pcout << "Total setup time               (wall) " << setup_time << "s\n";
 * 
 *     time.reset();
 *     time.start();
 *     constraints.set_zero(solution);
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner);
 * 
 *     constraints.distribute(solution);
 * 
 *     pcout << "Time solve (" << solver_control.last_step() << " iterations)"
 *           << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) "
 *           << time.cpu_time() << "s/" << time.wall_time() << "s\n";
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemoutput_results"></a> 
 * <h4>LaplaceProblem::output_results</h4>
 * 

 * 
 * Here is the data output, which is a simplified version of step-5. We use
 * the standard VTU (= compressed VTK) output for each grid produced in the
 * refinement process. In addition, we use a compression algorithm that is
 * optimized for speed rather than disk usage. The default setting (which
 * optimizes for disk usage) makes saving the output take about 4 times as
 * long as running the linear solver, while setting
 * DataOutBase::VtkFlags::compression_level to
 * DataOutBase::VtkFlags::best_speed lowers this to only one fourth the time
 * of the linear solve.
 *   

 * 
 * We disable the output when the mesh gets too large. A variant of this
 * program has been run on hundreds of thousands MPI ranks with as many as
 * 100 billion grid cells, which is not directly accessible to classical
 * visualization tools.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const
 *   {
 *     Timer time;
 *     if (triangulation.n_global_active_cells() > 1000000)
 *       return;
 * 
 *     DataOut<dim> data_out;
 * 
 *     solution.update_ghost_values();
 *     data_out.attach_dof_handler(dof_handler);
 *     data_out.add_data_vector(solution, "solution");
 *     data_out.build_patches(mapping);
 * 
 *     DataOutBase::VtkFlags flags;
 *     flags.compression_level = DataOutBase::VtkFlags::best_speed;
 *     data_out.set_flags(flags);
 *     data_out.write_vtu_with_pvtu_record(
 *       "./", "solution", cycle, MPI_COMM_WORLD, 3);
 * 
 *     time_details << "Time write output          (CPU/wall) " << time.cpu_time()
 *                  << "s/" << time.wall_time() << "s\n";
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * The function that runs the program is very similar to the one in
 * step-16. We do few refinement steps in 3D compared to 2D, but that's
 * it.
 *   

 * 
 * Before we run the program, we output some information about the detected
 * vectorization level as discussed in the introduction.
 * 
 * @code
 *   template <int dim>
 *   void LaplaceProblem<dim>::run()
 *   {
 *     {
 *       const unsigned int n_vect_doubles = VectorizedArray<double>::size();
 *       const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;
 * 
 *       pcout << "Vectorization over " << n_vect_doubles
 *             << " doubles = " << n_vect_bits << " bits ("
 *             << Utilities::System::get_current_vectorization_level() << ")"
 *             << std::endl;
 *     }
 * 
 *     for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle)
 *       {
 *         pcout << "Cycle " << cycle << std::endl;
 * 
 *         if (cycle == 0)
 *           {
 *             GridGenerator::hyper_cube(triangulation, 0., 1.);
 *             triangulation.refine_global(3 - dim);
 *           }
 *         triangulation.refine_global(1);
 *         setup_system();
 *         assemble_rhs();
 *         solve();
 *         output_results(cycle);
 *         pcout << std::endl;
 *       };
 *   }
 * } // namespace Step37
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
 * Apart from the fact that we set up the MPI framework according to step-40,
 * there are no surprises in the main function.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace Step37;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
 * 
 *       LaplaceProblem<dimension> laplace_problem;
 *       laplace_problem.run();
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
 * 
 *   return 0;
 * }
 * @endcode
examples/step-37/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


由于这个例子解决的是与步骤5相同的问题（除了不同的系数），所以对解决方案没有什么可说的。我们还是展示了一张图片，通过等高线和体积渲染来说明解决方案的大小。

 <img src="https://www.dealii.org/images/steps/developer/step-37.solution.png" alt=""> 

更有趣的是评估多网格求解器的某些方面。当我们在二维运行这个程序时，对于二次（ $Q_2$ ）元素，我们得到以下输出（当在一个核心上以释放模式运行时）。

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 81
Total setup time               (wall) 0.00159788s
Time solve (6 iterations)  (CPU/wall) 0.000951s/0.000951052s


Cycle 1
Number of degrees of freedom: 289
Total setup time               (wall) 0.00114608s
Time solve (6 iterations)  (CPU/wall) 0.000935s/0.000934839s


Cycle 2
Number of degrees of freedom: 1089
Total setup time               (wall) 0.00244665s
Time solve (6 iterations)  (CPU/wall) 0.00207s/0.002069s


Cycle 3
Number of degrees of freedom: 4225
Total setup time               (wall) 0.00678205s
Time solve (6 iterations)  (CPU/wall) 0.005616s/0.00561595s


Cycle 4
Number of degrees of freedom: 16641
Total setup time               (wall) 0.0241671s
Time solve (6 iterations)  (CPU/wall) 0.019543s/0.0195441s


Cycle 5
Number of degrees of freedom: 66049
Total setup time               (wall) 0.0967851s
Time solve (6 iterations)  (CPU/wall) 0.07457s/0.0745709s


Cycle 6
Number of degrees of freedom: 263169
Total setup time               (wall) 0.346374s
Time solve (6 iterations)  (CPU/wall) 0.260042s/0.265033s
@endcode



如同步骤16，我们看到随着自由度的增加，CG的迭代次数保持不变。恒定的迭代次数（加上最佳的计算特性）意味着当问题大小在一个周期内翻两番时，计算时间大约翻了四倍。该代码在存储方面也非常有效。大约200-400万个自由度适合于1GB的内存，也见下面的MPI结果。一个有趣的事实是，尽管没有建立矩阵，但解决一个线性系统比设置要便宜（大约一半的时间花在 DoFHandler::distribute_dofs() 和 DoFHandler::distribute_mg_dofs() 的调用上）。这表明这种方法的效率很高，但也表明deal.II数据结构的设置相当昂贵，设置成本必须在几个系统求解中摊销。

如果我们在三个空间维度上运行程序，就不会有太大变化。由于我们使用了均匀的网格细化，我们得到的元素数量是八倍，每个周期的自由度大约是八倍。

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 125
Total setup time               (wall) 0.00231099s
Time solve (6 iterations)  (CPU/wall) 0.000692s/0.000922918s


Cycle 1
Number of degrees of freedom: 729
Total setup time               (wall) 0.00289083s
Time solve (6 iterations)  (CPU/wall) 0.001534s/0.0024128s


Cycle 2
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0143182s
Time solve (6 iterations)  (CPU/wall) 0.010785s/0.0107841s


Cycle 3
Number of degrees of freedom: 35937
Total setup time               (wall) 0.087064s
Time solve (6 iterations)  (CPU/wall) 0.063522s/0.06545s


Cycle 4
Number of degrees of freedom: 274625
Total setup time               (wall) 0.596306s
Time solve (6 iterations)  (CPU/wall) 0.427757s/0.431765s


Cycle 5
Number of degrees of freedom: 2146689
Total setup time               (wall) 4.96491s
Time solve (6 iterations)  (CPU/wall) 3.53126s/3.56142s
@endcode



既然如此简单，我们看看如果我们增加多项式的度数会发生什么。当在三维中选择度数为4，即在 $\mathcal Q_4$ 元素上，通过改变程序顶部的一行<code>const unsigned int degree_finite_element=4;</code>，我们得到以下程序输出。

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 729
Total setup time               (wall) 0.00633097s
Time solve (6 iterations)  (CPU/wall) 0.002829s/0.00379395s


Cycle 1
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0174279s
Time solve (6 iterations)  (CPU/wall) 0.012255s/0.012254s


Cycle 2
Number of degrees of freedom: 35937
Total setup time               (wall) 0.082655s
Time solve (6 iterations)  (CPU/wall) 0.052362s/0.0523629s


Cycle 3
Number of degrees of freedom: 274625
Total setup time               (wall) 0.507943s
Time solve (6 iterations)  (CPU/wall) 0.341811s/0.345788s


Cycle 4
Number of degrees of freedom: 2146689
Total setup time               (wall) 3.46251s
Time solve (7 iterations)  (CPU/wall) 3.29638s/3.3265s


Cycle 5
Number of degrees of freedom: 16974593
Total setup time               (wall) 27.8989s
Time solve (7 iterations)  (CPU/wall) 26.3705s/27.1077s
@endcode



由于一定网格上的 $\mathcal Q_4$ 元素对应于一半网格大小的 $\mathcal Q_2$ 元素，我们可以比较第四周期使用四度多项式和第五周期使用二次多项式的运行时间，两者都是210万自由度。令人惊讶的效果是，尽管多用了一次线性迭代， $\mathcal Q_4$ 元素的求解器实际上比四次方的情况略快。高阶多项式的速度与低阶多项式类似，甚至比低阶多项式更快，这是通过和分解进行无矩阵算子评估的主要优势之一，见<a
href="http://dx.doi.org/10.1016/j.compfluid.2012.04.012">matrix-free
paper</a>。这与基于矩阵的方法有根本的不同，后者随着多项式度数的增加和耦合的密集，每个未知数的成本会越来越高。

此外，对于更高的订单，设置也变得更便宜，这是因为需要设置的元素更少。

最后，让我们看一下度数为8的时间，这相当于低阶方法的另一轮网格细化。

@code
Vectorization over 2 doubles = 128 bits (SSE2)
Cycle 0
Number of degrees of freedom: 4913
Total setup time               (wall) 0.0842004s
Time solve (8 iterations)  (CPU/wall) 0.019296s/0.0192959s


Cycle 1
Number of degrees of freedom: 35937
Total setup time               (wall) 0.327048s
Time solve (8 iterations)  (CPU/wall) 0.07517s/0.075999s


Cycle 2
Number of degrees of freedom: 274625
Total setup time               (wall) 2.12335s
Time solve (8 iterations)  (CPU/wall) 0.448739s/0.453698s


Cycle 3
Number of degrees of freedom: 2146689
Total setup time               (wall) 16.1743s
Time solve (8 iterations)  (CPU/wall) 3.95003s/3.97717s


Cycle 4
Number of degrees of freedom: 16974593
Total setup time               (wall) 130.8s
Time solve (8 iterations)  (CPU/wall) 31.0316s/31.767s
@endcode



在这里，初始化似乎比以前慢得多，这主要是由于矩阵对角线的计算，它实际上是在每个单元格上计算一个729 x 729的矩阵，扔掉除对角线以外的所有东西。然而，解算时间再次非常接近四次方的情况，这表明理论上预期的随着多项式程度的增加而出现的线性增长几乎完全被更好的计算特性和高阶方法在几个单元上的自由度份额较小而增加了评估的复杂性所抵消。

<a name="Comparisonwithasparsematrix"></a><h3>Comparison with a sparse matrix</h3>


为了了解无矩阵实现的能力，我们通过测量问题初始化的计算时间（分配DoF、设置和装配矩阵、设置多网格结构）以及无矩阵变体和基于稀疏矩阵的变体的实际求解时间，将上面的3D例子与基于稀疏矩阵的变体的性能进行比较。如上图所示，我们将预处理程序建立在浮点数上，将实际的矩阵和向量建立在双数上。测试在英特尔酷睿i7-5500U笔记本处理器（两个核心，支持<a
href="http://en.wikipedia.org/wiki/Advanced_Vector_Extensions">AVX</a>，即用一条CPU指令就可以完成对双数的四次操作，这在FEEvaluation中被大量使用）、优化模式和两个MPI行列上运行。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="2">Sparse matrix</th>
    <th colspan="2">Matrix-free implementation</th>
  </tr>
  <tr>
    <th>n_dofs</th>
    <th>Setup + assemble</th>
    <th>&nbsp;Solve&nbsp;</th>
    <th>Setup + assemble</th>
    <th>&nbsp;Solve&nbsp;</th>
  </tr>
  <tr>
    <td align="right">125</td>
    <td align="center">0.0042s</td>
    <td align="center">0.0012s</td>
    <td align="center">0.0022s</td>
    <td align="center">0.00095s</td>
  </tr>
  <tr>
    <td align="right">729</td>
    <td align="center">0.012s</td>
    <td align="center">0.0040s</td>
    <td align="center">0.0027s</td>
    <td align="center">0.0021s</td>
  </tr>
  <tr>
    <td align="right">4,913</td>
    <td align="center">0.082s</td>
    <td align="center">0.012s</td>
    <td align="center">0.011s</td>
    <td align="center">0.0057s</td>
  </tr>
  <tr>
    <td align="right">35,937</td>
    <td align="center">0.73s</td>
    <td align="center">0.13s</td>
    <td align="center">0.048s</td>
    <td align="center">0.040s</td>
  </tr>
  <tr>
    <td align="right">274,625</td>
    <td align="center">5.43s</td>
    <td align="center">1.01s</td>
    <td align="center">0.33s</td>
    <td align="center">0.25s</td>
  </tr>
  <tr>
    <td align="right">2,146,689</td>
    <td align="center">43.8s</td>
    <td align="center">8.24s</td>
    <td align="center">2.42s</td>
    <td align="center">2.06s</td>
  </tr>
</table> 

该表清楚地显示，无矩阵实现的求解速度是两倍以上，而在初始化成本方面，则是六倍以上。随着问题大小被放大8倍，我们注意到，时间通常也会上升8倍（因为求解器的迭代次数恒定为6次）。主要的偏差是在5k到36k自由度的稀疏矩阵中，时间增加了12倍。这是处理器中的（L3）缓存不能再容纳矩阵-向量乘积所需的所有数据的阈值，所有的矩阵元素必须从主内存中获取。

当然，这种情况不一定适用于所有情况，因为在有些问题上，对矩阵项的了解可以使解算器的效果好得多（如当系数的变化比上面的例子更强烈时）。此外，这也取决于计算机系统。目前的系统具有良好的内存性能，因此稀疏矩阵的性能相当好。尽管如此，对于本例中使用的<i>Q</i><sub>2</sub>元素，无矩阵的实现已经给出了一个不错的速度。这一点对于时间依赖性或非线性问题尤其明显，在这些问题中，稀疏矩阵需要一次又一次地被重新组合，有了这个类，这就变得容易多了。当然，由于产品的复杂性更好，当元素的阶数增加时，该方法获得了越来越大的优势（无矩阵实现每个自由度的成本为4<i>d</i><sup>2</sup><i>p</i>，而稀疏矩阵为2<i>p<sup>d</sup></i>，所以无论如何它在4阶以上的3d中会获胜）。

<a name="ResultsforlargescaleparallelcomputationsonSuperMUC"></a><h3> Results for large-scale parallel computations on SuperMUC</h3>


正如介绍和代码中的注释所解释的，这个程序可以用MPI并行运行。事实证明，几何多栅方案工作得非常好，可以扩展到非常大的机器。据作者所知，这里显示的几何多网格结果是截至2016年底用deal.II完成的最大计算，在<a
href="https://www.lrz.de/services/compute/supermuc/systemdescription/">complete
SuperMUC Phase 1</a>的多达147456个核心上运行。超过1000个核心的可扩展性的要素是，没有任何依赖于全局问题大小的数据结构被完整地保存在一个处理器上，并且通信不是太频繁，以避免遇到网络的延迟问题。  对于用迭代求解器求解的PDEs，通信延迟往往是限制因素，而不是网络的吞吐量。以SuperMUC系统为例，两个处理器之间的点对点延迟在1e-6到1e-5秒之间，取决于MPI网络中的距离。这一类的矩阵-向量产品与 @p LaplaceOperator 涉及几个点对点通信步骤，与每个核心上的计算交错进行。由此产生的矩阵-向量乘积的延迟约为1e-4秒。全局通信，例如一个 @p MPI_Allreduce 操作，在MPI网络中的所有等级上累积每个等级的单一数字之和，其延迟为1e-4秒。这个程序中使用的多网格V型循环也是全局通信的一种形式。想一想发生在单个处理器上的粗略网格求解。在开始之前，它积累了来自所有处理器的贡献。当完成后，粗网格解决方案被转移到更细的层次，在那里越来越多的处理器帮助平滑，直到细网格。从本质上讲，这是在网络中的处理器上的一个树状模式，并由网格控制。相对于 @p MPI_Allreduce 的操作，在还原中的树被优化为MPI网络中的实际链接，多网格V-cycle是根据网格的划分来做的。因此，我们不能期望有同样的优化效果。此外，多网格循环并不是简单地在细化树上走来走去，而是在做平滑的时候在每一层上进行通信。换句话说，多网格中的全局通信更具挑战性，与提供较少优化机会的网格有关。测得的V型周期的延迟在6e-3和2e-2秒之间，即与60至200次MPI_Allreduce操作相同。

下图显示了在 $\mathcal Q_3$ 元素上进行的缩放实验。沿着这条线，问题的大小保持不变，因为核的数量在增加。当内核数量增加一倍时，人们期望计算时间减少一半，灰色虚线表示。结果显示，在达到0.1秒左右的绝对时间之前，该实现显示了几乎理想的行为。解算器的公差已经被设定为解算器执行五次迭代。这种绘制数据的方式是该算法的<b>strong scaling</b>。当我们走到非常大的核心数时，曲线会提前变平，这是因为SuperMUC中的通信网络，距离较远的处理器之间的通信会稍慢一些。

 <img src="https://www.dealii.org/images/steps/developer/step-37.scaling_strong.png" alt=""> 

此外，该图还包含了<b>weak scaling</b>的结果，列出了当处理器内核和元素的数量都以相同的速度增加时，算法的表现。在这种情况下，我们期望计算时间保持不变。在算法上，CG的迭代次数恒定在5次，所以我们从这一点来看是好的。图中的线条是这样排列的：每个数据系列中的左上角点代表每个处理器的相同大小，即131,072个元素（或每个核心大约350万个自由度）。表示理想的强缩放的灰色线条相隔8个相同的系数。结果再次表明，缩放比例几乎是理想的。当从288个核到147456个核时，并行效率在75%左右，每个核的局部问题大小为75万自由度，在288个核上需要1.0秒，在2304个核上需要1.03秒，在18000个核上需要1.19秒，在147000个核上需要1.35秒。这些算法对处理器的利用率也达到了很高。在147k核心上最大的计算在SuperMUC上达到约1.7 PFLOPs/s，其中算术峰值为3.2 PFLOPs/s。对于一个迭代式PDE求解器来说，这是一个非常高的数字，而且通常只有密集线性代数才会达到显著的数字。稀疏线性代数被限制在这个数值的十分之一。

正如介绍中提到的，无矩阵方法减少了数据结构的内存消耗。除了由于更少的内存传输而带来的更高的性能外，该算法还允许非常大的问题被装入内存。下图显示了随着我们增加问题的大小，直到计算耗尽内存的上限时的计算时间。我们对1k核、8k核和65k核进行了计算，发现问题的大小几乎可以在两个数量级上进行理想的扩展。这张图中显示的最大的计算涉及2920亿（ $2.92 \cdot 10^{11}$ ）个自由度。在147k核心的DG计算中，上述算法也被运行，涉及多达5490亿（2^39）个自由度。

 <img src="https://www.dealii.org/images/steps/developer/step-37.scaling_size.png" alt=""> 

最后，我们注意到，在对上述大规模系统进行测试的同时，deal.II中的多网格算法也得到了改进。原始版本包含了基于MGSmootherPrecondition的次优代码，其中一些MPI_Allreduce命令（检查所有向量条目是否为零）在每一级的平滑操作上都要进行，这在65k核以上的系统中才变得明显。然而，下面的图片显示，改进已经在较小的规模上得到了回报，这里显示的是对 $\mathcal Q_5$ 元素在多达14336个内核上的计算。

 <img src="https://www.dealii.org/images/steps/developer/step-37.scaling_oldnew.png" alt=""> 




<a name="Adaptivity"></a><h3> Adaptivity</h3>


正如代码中所解释的，这里介绍的算法是为运行在自适应细化的网格上准备的。如果只有部分网格被细化，多网格循环将以局部平滑的方式运行，并通过 MatrixFreeOperators::Base 类对细化程度不同的界面施加迪里切条件进行平滑。由于自由度在层次上的分布方式，将层次单元的所有者与第一个下级活动单元的所有者联系起来，在MPI中不同的处理器之间可能存在不平衡，这限制了可扩展性，约为2到5倍。

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions</h3>


<a name="Kellyerrorestimator"></a><h4> Kelly error estimator </h4>


如上所述，代码已经准备好用于局部自适应h-精简。对于泊松方程，可以采用KellyErrorEstimator类中实现的Kelly误差指标。然而，我们需要小心处理平行向量的鬼魂指数。为了评估误差指标中的跳跃项，每个MPI进程需要知道本地相关的DoF。然而 MatrixFree::initialize_dof_vector() 函数只用一些本地相关的DoF来初始化向量。在向量中提供的鬼魂指数是一个严格的集合，只有那些在单元积分（包括约束解决）中被触及的指数。这种选择有性能上的原因，因为与矩阵-向量乘积相比，发送所有本地相关的自由度会过于昂贵。因此，原样的解决方案向量不适合KellyErrorEstimator类。诀窍是改变分区的幽灵部分，例如使用一个临时向量和 LinearAlgebra::distributed::Vector::copy_locally_owned_data_from() ，如下所示。

@code
IndexSet locally_relevant_dofs;
DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
LinearAlgebra::distributed::Vector<double> copy_vec(solution);
solution.reinit(dof_handler.locally_owned_dofs(),
                locally_relevant_dofs,
                triangulation.get_communicator());
solution.copy_locally_owned_data_from(copy_vec);
constraints.distribute(solution);
solution.update_ghost_values();
@endcode



<a name="Sharedmemoryparallelization"></a><h4> Shared-memory parallelization</h4>


这个程序只用MPI来并行化。作为一种选择，MatrixFree循环也可以在混合模式下发出，例如通过在集群的节点上使用MPI并行化，在一个节点的共享内存区域内通过Intel TBB使用线程。要使用这一点，就需要在主函数的MPI_InitFinalize数据结构中同时设置线程数，并将 MatrixFree::AdditionalData::tasks_parallel_scheme 设置为partition_color，以便真正并行地进行循环。这个用例将在步骤-48中讨论。

<a name="InhomogeneousDirichletboundaryconditions"></a><h4> Inhomogeneous Dirichlet boundary conditions </h4>


所提出的程序假定了同质的Dirichlet边界条件。当进入非均质条件时，情况就有点复杂了。为了理解如何实现这样的设置，让我们首先回顾一下这些条件是如何在数学公式中出现的，以及它们是如何在基于矩阵的变体中实现的。从本质上讲，非均质Dirichlet条件将解决方案中的一些节点值设定为给定值，而不是通过变分原理来确定它们。

@f{eqnarray*}
u_h(\mathbf{x}) = \sum_{i\in \mathcal N} \varphi_i(\mathbf{x}) u_i =
\sum_{i\in \mathcal N \setminus \mathcal N_D} \varphi_i(\mathbf{x}) u_i +
\sum_{i\in \mathcal N_D} \varphi_i(\mathbf{x}) g_i,


@f}

其中 $u_i$ 表示解决方案的节点值， $\mathcal N$ 表示所有节点的集合。集合 $\mathcal N_D\subset \mathcal N$ 是受迪里希特边界条件约束的节点子集，其中解被强制等于 $u_i = g_i = g(\mathbf{x}_i)$ 作为迪里希特约束的节点点上的边界值插值 $i\in \mathcal
N_D$  。然后我们把这个解的表示插入到弱的形式中，例如上面所示的拉普拉斯，并把已知量移到右边。

@f{eqnarray*}
(\nabla \varphi_i, \nabla u_h)_\Omega &=& (\varphi_i, f)_\Omega \quad \Rightarrow \\
\sum_{j\in \mathcal N \setminus \mathcal N_D}(\nabla \varphi_i,\nabla \varphi_j)_\Omega \, u_j &=&
(\varphi_i, f)_\Omega


-\sum_{j\in \mathcal N_D} (\nabla \varphi_i,\nabla\varphi_j)_\Omega\, g_j.


@f}

在这个公式中，对所有的基函数 $\varphi_i$ 与 $i\in N \setminus \mathcal N_D$ 进行测试，这些基函数与迪里希特条件约束的节点没有关系。

在deal.II的实现中，右手边的积分 $(\nabla \varphi_i,\nabla \varphi_j)_\Omega$ 已经包含在我们在每个单元格上组装的局部矩阵贡献中。当使用 AffineConstraints::distribute_local_to_global() 时，正如在步骤6和步骤7的教程程序中首次描述的那样，我们可以通过将本地矩阵的列<i>j</i>和行<i>i</i>相乘来说明不均匀约束的贡献<i>j</i> 的局部矩阵根据积分 $(\varphi_i,
\varphi_j)_\Omega$ 乘以不均匀性，然后从全局右侧向量中的位置<i>i</i>中减去所得，也见 @ref
constraints  模块。实质上，我们使用一些从方程左侧被消除的积分来最终确定右侧的贡献。当首先将所有条目写进左侧矩阵，然后通过 MatrixTools::apply_boundary_values(). 消除矩阵的行和列时，也会涉及类似的数学。

原则上，属于受限自由度的成分可以从线性系统中剔除，因为它们不携带任何信息。实际上，在deal.II中，我们总是保持线性系统的大小不变，以避免处理两种不同的编号系统，并避免对两种不同的索引集产生混淆。为了确保在不向受限行添加任何东西时，线性系统不会变得奇异，我们再向矩阵对角线添加假条目，否则与真实条目无关。

在无矩阵方法中，我们需要采取不同的方法，因为 @p LaplaceOperator类代表了<b>homogeneous</b>算子的矩阵-向量乘积（最后一个公式的左手边）。  传递给 MatrixFree::reinit() 的AffineConstraints对象是否包含不均匀约束并不重要，只要它代表一个<b>linear</b>算子， MatrixFree::cell_loop() 调用将只解决约束的同质部分。

在我们的无矩阵代码中，非均质条件的贡献最终会在右侧计算中与矩阵算子完全脱钩，并由上述不同的函数处理。因此，我们需要明确地生成进入右手边的数据，而不是使用矩阵装配的副产品。由于我们已经知道如何在一个向量上应用算子，我们可以尝试对一个向量使用这些设施，我们只设置Dirichlet值。

@code
  // interpolate boundary values on vector solution
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(mapping,
                                           dof_handler,
                                           0,
                                           BoundaryValueFunction<dim>(),
                                           boundary_values);
  for (const std::pair<const types::global_dof_index, double> &pair : boundary_values)
    if (solution.locally_owned_elements().is_element(pair.first))
      solution(pair.first) = pair.second;
@endcode

或者说，如果我们已经将不均匀约束填充到AffineConstraints对象中。

@code
  solution = 0;
  constraints.distribute(solution);
@endcode



然后我们可以将向量 @p solution 传递给 @p  LaplaceOperator::vmult_add() 函数，并将新的贡献添加到 @p system_rhs向量中，在 @p LaplaceProblem::assemble_rhs() 函数中被填充。然而，这个想法并不奏效，因为vmult()函数中使用的 FEEvaluation::read_dof_values() 调用假定所有约束条件的值都是同质的（否则运算符就不是线性运算符，而是仿射运算符）。为了同时检索不均匀性的值，我们可以选择以下两种策略中的一种。

<a name="UseFEEvaluationread_dof_values_plaintoavoidresolvingconstraints"></a><h5> Use FEEvaluation::read_dof_values_plain() to avoid resolving constraints </h5>


FEEvaluation类有一个设施，正是为了解决这个要求。对于非均质的Dirichlet值，我们确实希望在从向量 @p solution. 中读取数据时跳过隐含的均质（Dirichlet）约束。例如，我们可以扩展 @p  LaplaceProblem::assemble_rhs() 函数来处理非均质的Dirichlet值，如下所示，假设Dirichlet值已经被插值到对象 @p constraints:  中

@code
template <int dim>
void LaplaceProblem<dim>::assemble_rhs()
{
  solution = 0;
  constraints.distribute(solution);
  solution.update_ghost_values();
  system_rhs = 0;


  const Table<2, VectorizedArray<double>> &coefficient = system_matrix.get_coefficient();
  FEEvaluation<dim, degree_finite_element> phi(*system_matrix.get_matrix_free());
  for (unsigned int cell = 0;
       cell < system_matrix.get_matrix_free()->n_cell_batches();
       ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(solution);
      phi.evaluate(EvaluationFlags::gradients);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          phi.submit_gradient(-coefficient(cell, q) * phi.get_gradient(q), q);
          phi.submit_value(make_vectorized_array<double>(1.0), q);
        }
      phi.integrate(EvaluationFlags::values|EvaluationFlags::gradients);
      phi.distribute_local_to_global(system_rhs);
    }
  system_rhs.compress(VectorOperation::add);
}
@endcode



在这段代码中，我们用忽略所有约束的 FEEvaluation::read_dof_values_plain() 代替了用于暂定解向量的 FEEvaluation::read_dof_values() 函数。由于这种设置，我们必须确保其他约束条件，例如通过悬挂节点，已经正确地分布到输入向量中，因为它们没有像 FEEvaluation::read_dof_values_plain(). 那样被解决。 在循环内部，我们然后评估拉普拉斯，并用 @p LaplaceOperator 类中的 FEEvaluation::submit_gradient() 重复二次导数调用，但符号调换，因为我们想根据上述公式减去右侧向量的迪里希条件的贡献。当我们调用 FEEvaluation::integrate() 时，我们将关于值槽和第一导数槽的两个参数设置为真，以说明在正交点的循环中加入的两个项。一旦右手边集合完毕，我们就继续求解同质问题的线性系统，比如说涉及到一个变量 @p solution_update. 在求解之后，我们可以将 @p solution_update 加入到包含最终（非同质）解决方案的 @p solution 向量中。

请注意，拉普拉斯的负号与我们需要用来建立右手边的强制力的正号是一个更普遍的概念。我们所实施的只不过是牛顿的非线性方程方法，但应用于线性系统。我们在迪里切特边界条件方面使用了对变量 @p solution 的初始猜测，并计算了残差 $r = f - Au_0$  。然后线性系统被解为  $\Delta u = A^{-1} (f-Au)$  ，我们最后计算出  $u = u_0 + \Delta u$  。对于一个线性系统，我们显然在一次迭代后就能达到精确解。如果我们想将代码扩展到非线性问题，我们会将 @p assemble_rhs() 函数重新命名为一个更具描述性的名字，如 @p  assemble_residual()，计算残差的（弱）形式，而 @p LaplaceOperator::apply_add() 函数将得到残差相对于解变量的线性化。

<a name="UseLaplaceOperatorwithasecondAffineConstraintsobjectwithoutDirichletconditions"></a><h5> Use LaplaceOperator with a second AffineConstraints object without Dirichlet conditions </h5>


获得重新使用 @p  LaplaceOperator::apply_add() 函数的第二个替代方法是添加第二个LaplaceOperator，跳过Dirichlet约束。为了做到这一点，我们初始化第二个MatrixFree对象，它没有任何边界值约束。这个 @p matrix_free 对象然后被传递给一个 @p LaplaceOperator 类实例 @p inhomogeneous_operator，它只用于创建右手边。

@code
template <int dim>
void LaplaceProblem<dim>::assemble_rhs()
{
  system_rhs = 0;
  AffineConstraints<double> no_constraints;
  no_constraints.close();
  LaplaceOperator<dim, degree_finite_element, double> inhomogeneous_operator;


  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_gradients | update_JxW_values | update_quadrature_points);
  std::shared_ptr<MatrixFree<dim, double>> matrix_free(
    new MatrixFree<dim, double>());
  matrix_free->reinit(dof_handler,
                      no_constraints,
                      QGauss<1>(fe.degree + 1),
                      additional_data);
  inhomogeneous_operator.initialize(matrix_free);


  solution = 0.0;
  constraints.distribute(solution);
  inhomogeneous_operator.evaluate_coefficient(Coefficient<dim>());
  inhomogeneous_operator.vmult(system_rhs, solution);
  system_rhs *= -1.0;


  FEEvaluation<dim, degree_finite_element> phi(
    *inhomogeneous_operator.get_matrix_free());
  for (unsigned int cell = 0;
       cell < inhomogeneous_operator.get_matrix_free()->n_cell_batches();
       ++cell)
    {
      phi.reinit(cell);
      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        phi.submit_value(make_vectorized_array<double>(1.0), q);
      phi.integrate(EvaluationFlags::values);
      phi.distribute_local_to_global(system_rhs);
    }
  system_rhs.compress(VectorOperation::add);
}
@endcode



这种技术的更复杂的实现可以重新使用原始的MatrixFree对象。这可以通过用多个块初始化MatrixFree对象来实现，其中每个块对应于不同的AffineConstraints对象。这样做需要对LaplaceOperator类进行大量的修改，但是库中的 MatrixFreeOperators::LaplaceOperator 类可以做到这一点。关于如何设置块的更多信息，请参见 MatrixFreeOperators::Base 中关于块的讨论。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-37.cc"
*/
