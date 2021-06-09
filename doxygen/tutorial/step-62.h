/**
@page step_62 The step-62 tutorial program
This tutorial depends on step-8, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Elasticequations">Elastic equations</a>
        <li><a href="#Simulationparameters">Simulation parameters</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Auxiliaryclassesandfunctions">Auxiliary classes and functions</a>
      <ul>
        <li><a href="#TheRightHandSideclass">The `RightHandSide` class</a>
        <li><a href="#ThePMLclass">The `PML` class</a>
        <li><a href="#TheRhoclass">The `Rho` class</a>
        <li><a href="#TheParametersclass">The `Parameters` class</a>
        <li><a href="#TheQuadratureCacheclass">The `QuadratureCache` class</a>
        <li><a href="#Theget_stiffness_tensorfunction">The `get_stiffness_tensor()` function</a>
      </ul>
        <li><a href="#TheElasticWaveclass">The `ElasticWave` class</a>
        <li><a href="#Implementationoftheauxiliaryclasses">Implementation of the auxiliary classes</a>
      <ul>
        <li><a href="#TheRightHandSideclassimplementation">The `RightHandSide` class implementation</a>
        <li><a href="#ThePMLclassimplementation">The `PML` class implementation</a>
        <li><a href="#TheRhoclassimplementation">The `Rho` class implementation</a>
        <li><a href="#TheParametersclassimplementation">The `Parameters` class implementation</a>
        <li><a href="#TheQuadratureCacheclassimplementation">The `QuadratureCache` class implementation</a>
      </ul>
        <li><a href="#ImplementationoftheElasticWaveclass">Implementation of the `ElasticWave` class</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#ElasticWavesetup_system">ElasticWave::setup_system</a>
        <li><a href="#ElasticWaveassemble_system">ElasticWave::assemble_system</a>
        <li><a href="#ElasticWavesolve">ElasticWave::solve</a>
        <li><a href="#ElasticWaveinitialize_position_vector">ElasticWave::initialize_position_vector</a>
        <li><a href="#ElasticWavestore_frequency_step_data">ElasticWave::store_frequency_step_data</a>
        <li><a href="#ElasticWaveoutput_results">ElasticWave::output_results</a>
        <li><a href="#ElasticWavesetup_quadrature_cache">ElasticWave::setup_quadrature_cache</a>
        <li><a href="#ElasticWavefrequency_sweep">ElasticWave::frequency_sweep</a>
        <li><a href="#ElasticWaverun">ElasticWave::run</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Resonancefrequencyandbandgap">Resonance frequency and bandgap</a>
        <li><a href="#Modeprofile">Mode profile</a>
        <li><a href="#Experimentalapplications">Experimental applications</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-62/doc/intro.dox

 <br> 

<i>This program was contributed by Daniel Garcia-Sanchez.</i> <br>  。




 @note  作为这个程序的前提条件，你需要安装HDF5、复杂的PETSc和p4est库。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II与这些附加库的安装情况。

<a name="Introduction"></a><h1>Introduction</h1> 声子晶体是一种周期性的纳米结构，可以改变机械振动或[声子]的运动（https://en.wikipedia.org/wiki/Phonon）。声子结构可用于分散、引导和限制机械振动。这些结构在[量子信息](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.86.1391)方面有潜在的应用，并已被用于研究[宏观量子现象](https://science.sciencemag.org/content/358/6360/203)。声子晶体通常是在[洁净室]中制造的(https://en.wikipedia.org/wiki/Cleanroom)。


在本教程中，我们展示了如何设计一个[声子超晶格空腔](https://doi.org/10.1103/PhysRevA.94.033813)，这是一种特殊类型的声子晶体，可用于限制机械振动。声子超晶格空腔是由两个[分布式布拉格反射器](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)、镜子和一个 $\lambda/2$ 空腔组成，其中 $\lambda$ 是声学波长。声学DBRs是周期性结构，其中一组具有对比物理特性（声速指数）的双层堆栈被重复 $N$ 次。超晶格空腔通常通过[分子束外延](https://en.wikipedia.org/wiki/Molecular-beam_epitaxy)在[砷化镓](https://en.wikipedia.org/wiki/Gallium_arsenide)晶片上生长。双层对应于砷化镓/砷化铝镜像对。如下图所示，镜像层（棕色和绿色）的厚度为 $\lambda/4$ ，空腔（蓝色）的厚度为 $\lambda/2$  。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.01.svg" height="200" /> 

在本教程中，我们计算了[带隙](https://en.wikipedia.org/wiki/Band_gap)和声子超晶格空腔的机械共振，但这里介绍的代码可以很容易地用于设计和计算其他类型的[声子晶体](https://science.sciencemag.org/content/358/6360/203)。

该装置是一个波导，其中的波从左到右。本教程的模拟是在二维进行的，但代码是独立于维度的，可以很容易地用于三维模拟。波导的宽度等于域的 $y$ 维，波导的长度等于域的 $x$ 维。有两个取决于波导宽度的制度。

- 单一模式。在这种情况下，结构的宽度要比波长小得多。   这种情况可以用有限元法（我们在这里采取的方法）或用简单的半分析法[一维转移矩阵形式]（https://en.wikipedia.org/wiki/Transfer_matrix）来解决。

- 多模。在这种情况下，结构的宽度比波长大。   这种情况可以用有限元法或[散射矩阵形式主义]（https://doi.org/10.1103/PhysRevA.94.033813）来解决。   尽管我们在本教程中没有研究这种情况，但通过增加波导宽度参数（jupyter笔记本中的`dimension_y'），很容易达到多模制度。

本教程的模拟是在频域进行的。为了计算传输频谱，我们使用了时域[FDTD](https://meep.readthedocs.io/en/latest/Python_Tutorials/Resonant_Modes_and_Transmission_in_a_Waveguide_Cavity/)模拟中常用的一个[程序]。在结构的左侧产生一个特定频率的脉冲，在结构的右侧测量传输的能量。仿真运行了两次。首先，我们运行声子结构的模拟并测量传输能量。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.02.svg" height="200" /> 

然后，我们运行没有声子结构的模拟，并测量传输的能量。我们使用没有结构的模拟来进行校准。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.03.svg" height="200" /> 

传输系数相当于第一次模拟的能量除以校准能量。我们对每个频率步骤重复这一程序。




<a name="Elasticequations"></a><h3>Elastic equations</h3> 我们在这里要模拟的是弹性波的传输。因此，对问题的正确描述使用了弹性方程，在时域中，弹性方程由以下几项给出


@f[
\rho\partial_{tt} u_i - \partial_j (c_{ijkl} \varepsilon_{kl}) = f_i,
\qquad i=0,1,2


@f]

其中刚度张量 $c_{ijkl}$ 取决于空间坐标，应变是位移的对称梯度，由以下公式给出

@f[
\varepsilon_{kl} =\frac{1}{2}(\partial_k u_l + \partial_l u_k)


@f]



[完美匹配层（PML）](https://en.wikipedia.org/wiki/Perfectly_matched_layer)可以用来在边界处截断解决方案。PML是一种导致复杂坐标拉伸的变换。

本教程程序没有采用时域方法，而是通过对时间变量进行傅里叶变换，将上述方程转换为频域。频域中的弹性方程的内容如下

@f{eqnarray*}
\nabla\cdot(\boldsymbol{\bar\sigma} \xi \boldsymbol{\Lambda})&=&-\omega^2\rho\xi\mathbf{\bar u}\\
\boldsymbol{\bar \sigma} &=&\mathbf{C}\boldsymbol{\bar\varepsilon}\\
\boldsymbol{\bar\varepsilon}&=&\frac{1}{2}[(\nabla\mathbf{\bar{u}}\boldsymbol{\Lambda}+\boldsymbol{\Lambda}^\mathrm{T}(\nabla\mathbf{\bar{u}})^\mathrm{T})]\\
\xi &=&\prod_i^\textrm{dim}s_i\\
\boldsymbol{\Lambda} &=& \operatorname{diag}(1/s_0,1/s_1,1/s_2)\qquad\textrm{for 3D}\\
\boldsymbol{\Lambda} &=& \operatorname{diag}(1/s_0,1/s_1)\qquad\textrm{for 2D}


@f}

其中系数 $s_i = 1+is_i'(x,y,z)$ 说明了吸收情况。3D中有3个 $s_i$ 系数，2D中有2个。 $s_i$ 的虚部在PML外等于零。PML仅对精确的波浪方程是无反射的。当方程组被离散化时，PML就不再是无反射的了。只要介质是缓慢变化的，反射就可以变得任意小，见[绝热定理](https://doi.org/10.1103/PhysRevE.66.066608)。在代码中，已经使用了PML的二次开启。线性和立方开启也是[已知可行的](https://doi.org/10.1364/OE.16.011376)。这些方程可以扩展为

@f[


-\omega^2\rho \xi  u_m - \partial_n \left(\frac{\xi}{s_n}c_{mnkl}
\varepsilon_{kl}\right) = f_m


@f]



@f[
\varepsilon_{kl} =\frac{1}{2}\left(\frac{1}{s_k}\partial_k u_l
+ \frac{1}{s_l}\partial_l u_k\right)


@f]

其中对重复指数（这里是 $n$ ，以及 $k$ 和 $l$ ）的求和一如既往地隐含着。请注意，应用PML的复数坐标拉伸后，应变不再是对称的。这组方程可以写成

@f[


-\omega^2\rho \xi  u_m - \partial_n \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k u_l
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l u_k\right) = f_m


@f]



与应变一样，应力张量在PML内也不是对称的（ $s_j\neq 0$ ）。事实上，PML内部的场不是物理的。介绍张量 $\alpha_{mnkl}$ 和 $\beta_{mnkl}$ 是有用的。

@f[


-\omega^2\rho \xi  u_m - \partial_n \left(\alpha_{mnkl}\partial_k u_l
+  \beta_{mnkl}\partial_l u_k\right) = f_m


@f]



我们可以乘以 $\varphi_m$ 并在 $\Omega$ 域上进行积分，并进行部分积分。

@f{eqnarray*}


-\omega^2\int_\Omega\rho\xi\varphi_m u_m + \int_\Omega\partial_n\varphi_m \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k u_l
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l u_k\right) = \int_\Omega\varphi_m f_m


@f}

正是这组方程，我们要解决一组频率 $\omega$ ，以计算传输系数与频率的关系。这个线性系统变成

@f{eqnarray*}
AU&=&F\\
A_{ij} &=& -\omega^2\int_\Omega\rho \xi\varphi_m^i \varphi_m^j + \int_\Omega\partial_n\varphi_m^i \left(\frac{\xi c_{mnkl}}{2s_n s_k} \partial_k \varphi_l^j
+ \frac{\xi c_{mnkl}}{2s_n s_l} \partial_l \varphi_k^j\right)\\
F_i &=& \int_\Omega\varphi_m^i f_m


@f}



<a name="Simulationparameters"></a><h3>Simulation parameters</h3> 在本教程中，我们使用python [jupyter notebook](https://github.com/dealii/dealii/blob/master/example/step-62/step-62.ipynb)来设置参数和运行模拟。首先，我们创建一个HDF5文件，在其中存储参数和模拟的结果。


每个模拟（位移和校准）都存储在一个单独的HDF5组中。

@code{.py}
import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import scipy.constants as constants
import scipy.optimize


# This considerably reduces the size of the svg data
plt.rcParams['svg.fonttype'] = 'none'


h5_file = h5py.File('results.h5', 'w')
data = h5_file.create_group('data')
displacement = data.create_group('displacement')
calibration = data.create_group('calibration')


# Set the parameters
for group in [displacement, calibration]:
    # Dimensions of the domain
    # The waveguide length is equal to dimension_x
    group.attrs['dimension_x'] = 2e-5
    # The waveguide width is equal to dimension_y
    group.attrs['dimension_y'] = 2e-8


    # Position of the probe that we use to measure the flux
    group.attrs['probe_pos_x']   = 8e-6
    group.attrs['probe_pos_y']   = 0
    group.attrs['probe_width_y'] = 2e-08


    # Number of points in the probe
    group.attrs['nb_probe_points'] = 5


    # Global refinement
    group.attrs['grid_level'] = 1


    # Cavity
    group.attrs['cavity_resonance_frequency'] = 20e9
    group.attrs['nb_mirror_pairs']            = 15


    # Material
    group.attrs['poissons_ratio'] = 0.27
    group.attrs['youngs_modulus'] = 270000000000.0
    group.attrs['material_a_rho'] = 3200
    if group == displacement:
        group.attrs['material_b_rho'] = 2000
    else:
        group.attrs['material_b_rho'] = 3200
    group.attrs['lambda'] = (group.attrs['youngs_modulus'] * group.attrs['poissons_ratio'] /
                           ((1 + group.attrs['poissons_ratio']) *
                           (1 - 2 * group.attrs['poissons_ratio'])))
    group.attrs['mu']= (group.attrs['youngs_modulus'] / (2 * (1 + group.attrs['poissons_ratio'])))


    # Force
    group.attrs['max_force_amplitude'] = 1e26
    group.attrs['force_sigma_x']       = 1e-7
    group.attrs['force_sigma_y']       = 1
    group.attrs['max_force_width_x']   = 3e-7
    group.attrs['max_force_width_y']   = 2e-8
    group.attrs['force_x_pos']         = -8e-6
    group.attrs['force_y_pos']         = 0


    # PML
    group.attrs['pml_x']            = True
    group.attrs['pml_y']            = False
    group.attrs['pml_width_x']      = 1.8e-6
    group.attrs['pml_width_y']      = 5e-7
    group.attrs['pml_coeff']        = 1.6
    group.attrs['pml_coeff_degree'] = 2


    # Frequency sweep
    group.attrs['center_frequency']    = 20e9
    group.attrs['frequency_range']     = 0.5e9
    group.attrs['start_frequency']     = group.attrs['center_frequency'] - group.attrs['frequency_range'] / 2
    group.attrs['stop_frequency']      = group.attrs['center_frequency'] + group.attrs['frequency_range'] / 2
    group.attrs['nb_frequency_points'] = 400


    # Other parameters
    if group == displacement:
        group.attrs['simulation_name'] = 'phononic_cavity_displacement'
    else:
        group.attrs['simulation_name'] = 'phononic_cavity_calibration'
    group.attrs['save_vtu_files'] = False


h5_file.close()
@endcode




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
 * Most of the include files we need for this program have already been
 * discussed in previous programs, in particular in step-40.
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h>
 * #include <deal.II/base/function.h>
 * 
 * #include <deal.II/base/index_set.h>
 * #include <deal.II/base/quadrature_lib.h>
 * #include <deal.II/base/timer.h>
 * #include <deal.II/base/utilities.h>
 * 
 * #include <deal.II/dofs/dof_handler.h>
 * #include <deal.II/dofs/dof_tools.h>
 * 
 * #include <deal.II/fe/fe_q.h>
 * #include <deal.II/fe/fe_system.h>
 * #include <deal.II/fe/fe_values.h>
 * 
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_refinement.h>
 * 
 * #include <deal.II/lac/affine_constraints.h>
 * #include <deal.II/lac/dynamic_sparsity_pattern.h>
 * #include <deal.II/lac/full_matrix.h>
 * #include <deal.II/lac/generic_linear_algebra.h>
 * #include <deal.II/lac/petsc_solver.h>
 * #include <deal.II/lac/vector.h>
 * 
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/error_estimator.h>
 * 
 * #include <fstream>
 * #include <iostream>
 * 
 * @endcode
 * 
 * The following header provides the Tensor class that we use to represent the
 * material properties.
 * 
 * @code
 * #include <deal.II/base/tensor.h>
 * 
 * 
 * @endcode
 * 
 * The following header is necessary for the HDF5 interface of deal.II.
 * 
 * @code
 * #include <deal.II/base/hdf5.h>
 * 
 * @endcode
 * 
 * This header is required for the function VectorTools::point_value that we use
 * to evaluate the result of the simulation.
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * We need these headers for the function
 * GridTools::find_active_cell_around_point that we use in the function
 * `ElasticWave::store_frequency_step_data()`
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h>
 * #include <deal.II/grid/grid_tools_cache.h>
 * 
 * namespace step62
 * {
 *   using namespace dealii;
 * 
 * @endcode
 * 
 * 
 * <a name="Auxiliaryclassesandfunctions"></a> 
 * <h3>Auxiliary classes and functions</h3>
 * The following classes are used to store the parameters of the simulation.
 * 

 * 
 * 
 * <a name="TheRightHandSideclass"></a> 
 * <h4>The `RightHandSide` class</h4>
 * This class is used to define the force pulse on the left side of the
 * structure:
 * 
 * @code
 *   template <int dim>
 *   class RightHandSide : public Function<dim>
 *   {
 *   public:
 *     RightHandSide(HDF5::Group &data);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component) const override;
 * 
 *   private:
 * @endcode
 * 
 * The variable `data` is the HDF5::Group in which all the simulation
 * results will be stored. Note that the variables `RightHandSide::data`,
 * `PML::data`, `Rho::data` and `Parameters::data` point to the same group
 * of the HDF5 file. When a HDF5::Group is copied, it will point to the same
 * group of the HDF5 file.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The simulation parameters are stored in `data` as HDF5 attributes. The
 * following attributes are defined in the jupyter notebook, stored in
 * `data` as HDF5 attributes and then read by the constructor.
 * 
 * @code
 *     const double     max_force_amplitude;
 *     const double     force_sigma_x;
 *     const double     force_sigma_y;
 *     const double     max_force_width_x;
 *     const double     max_force_width_y;
 *     const Point<dim> force_center;
 * 
 *   public:
 * @endcode
 * 
 * In this particular simulation the force has only a $x$ component,
 * $F_y=0$.
 * 
 * @code
 *     const unsigned int force_component = 0;
 *   };
 * 
 * @endcode
 * 
 * 
 * <a name="ThePMLclass"></a> 
 * <h4>The `PML` class</h4>
 * This class is used to define the shape of the Perfectly Matches
 * Layer (PML) to absorb waves traveling towards the boundary:
 * 
 * @code
 *   template <int dim>
 *   class PML : public Function<dim, std::complex<double>>
 *   {
 *   public:
 *     PML(HDF5::Group &data);
 * 
 *     virtual std::complex<double>
 *     value(const Point<dim> &p, const unsigned int component) const override;
 * 
 *   private:
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const double pml_coeff;
 *     const int    pml_coeff_degree;
 *     const double dimension_x;
 *     const double dimension_y;
 *     const bool   pml_x;
 *     const bool   pml_y;
 *     const double pml_width_x;
 *     const double pml_width_y;
 *     const double a_coeff_x;
 *     const double a_coeff_y;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRhoclass"></a> 
 * <h4>The `Rho` class</h4>
 * This class is used to define the mass density.
 * 
 * @code
 *   template <int dim>
 *   class Rho : public Function<dim>
 *   {
 *   public:
 *     Rho(HDF5::Group &data);
 * 
 *     virtual double value(const Point<dim> & p,
 *                          const unsigned int component = 0) const override;
 * 
 *   private:
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const double       lambda;
 *     const double       mu;
 *     const double       material_a_rho;
 *     const double       material_b_rho;
 *     const double       cavity_resonance_frequency;
 *     const unsigned int nb_mirror_pairs;
 *     const double       dimension_y;
 *     const unsigned int grid_level;
 *     double             average_rho_width;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheParametersclass"></a> 
 * <h4>The `Parameters` class</h4>
 * This class contains all the parameters that will be used in the simulation.
 * 
 * @code
 *   template <int dim>
 *   class Parameters
 *   {
 *   public:
 *     Parameters(HDF5::Group &data);
 * 
 * @endcode
 * 
 * HDF5::Group in which all the simulation results will be stored.
 * 
 * @code
 *     HDF5::Group data;
 * 
 * @endcode
 * 
 * The same as before, the following attributes are defined in the jupyter
 * notebook, stored in `data` as HDF5 attributes and then read by the
 * constructor.
 * 
 * @code
 *     const std::string        simulation_name;
 *     const bool               save_vtu_files;
 *     const double             start_frequency;
 *     const double             stop_frequency;
 *     const unsigned int       nb_frequency_points;
 *     const double             lambda;
 *     const double             mu;
 *     const double             dimension_x;
 *     const double             dimension_y;
 *     const unsigned int       nb_probe_points;
 *     const unsigned int       grid_level;
 *     const Point<dim>         probe_start_point;
 *     const Point<dim>         probe_stop_point;
 *     const RightHandSide<dim> right_hand_side;
 *     const PML<dim>           pml;
 *     const Rho<dim>           rho;
 * 
 *   private:
 *     const double comparison_float_constant = 1e-12;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheQuadratureCacheclass"></a> 
 * <h4>The `QuadratureCache` class</h4>
 * The calculation of the mass and stiffness matrices is very expensive. These
 * matrices are the same for all the frequency steps. The right hand side
 * vector is also the same for all the frequency steps. We use this class to
 * store these objects and re-use them at each frequency step. Note that here
 * we don't store the assembled mass and stiffness matrices and right hand
 * sides, but instead the data for a single cell. `QuadratureCache` class is
 * very similar to the `PointHistory` class that has been used in step-18.
 * 
 * @code
 *   template <int dim>
 *   class QuadratureCache
 *   {
 *   public:
 *     QuadratureCache(const unsigned int dofs_per_cell);
 * 
 *   private:
 *     unsigned int dofs_per_cell;
 * 
 *   public:
 * @endcode
 * 
 * We store the mass and stiffness matrices in the variables
 * mass_coefficient and stiffness_coefficient. We store as well the
 * right_hand_side and JxW values which are going to be the same for all the
 * frequency steps.
 * 
 * @code
 *     FullMatrix<std::complex<double>>  mass_coefficient;
 *     FullMatrix<std::complex<double>>  stiffness_coefficient;
 *     std::vector<std::complex<double>> right_hand_side;
 *     double                            JxW;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Theget_stiffness_tensorfunction"></a> 
 * <h4>The `get_stiffness_tensor()` function</h4>
 * 

 * 
 * This function returns the stiffness tensor of the material. For the sake of
 * simplicity we consider the stiffness to be isotropic and homogeneous; only
 * the density $\rho$ depends on the position. As we have previously shown in
 * step-8, if the stiffness is isotropic and homogeneous, the stiffness
 * coefficients $c_{ijkl}$ can be expressed as a function of the two
 * coefficients $\lambda$ and $\mu$. The coefficient tensor reduces to
 * @f[
 * c_{ijkl}
 * =
 * \lambda \delta_{ij} \delta_{kl} +
 * \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).
 * @f]
 * 
 * @code
 *   template <int dim>
 *   SymmetricTensor<4, dim> get_stiffness_tensor(const double lambda,
 *                                                const double mu)
 *   {
 *     SymmetricTensor<4, dim> stiffness_tensor;
 *     for (unsigned int i = 0; i < dim; ++i)
 *       for (unsigned int j = 0; j < dim; ++j)
 *         for (unsigned int k = 0; k < dim; ++k)
 *           for (unsigned int l = 0; l < dim; ++l)
 *             stiffness_tensor[i][j][k][l] =
 *               (((i == k) && (j == l) ? mu : 0.0) +
 *                ((i == l) && (j == k) ? mu : 0.0) +
 *                ((i == j) && (k == l) ? lambda : 0.0));
 *     return stiffness_tensor;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheElasticWaveclass"></a> 
 * <h3>The `ElasticWave` class</h3>
 * 

 * 
 * Next let's declare the main class of this program. Its structure is very
 * similar to the step-40 tutorial program. The main differences are:
 * - The sweep over the frequency values.
 * - We save the stiffness and mass matrices in `quadrature_cache` and
 * use them for each frequency step.
 * - We store the measured energy by the probe for each frequency step in the
 * HDF5 file.
 * 
 * @code
 *   template <int dim>
 *   class ElasticWave
 *   {
 *   public:
 *     ElasticWave(const Parameters<dim> &parameters);
 *     void run();
 * 
 *   private:
 *     void setup_system();
 *     void assemble_system(const double omega,
 *                          const bool   calculate_quadrature_data);
 *     void solve();
 *     void initialize_probe_positions_vector();
 *     void store_frequency_step_data(const unsigned int frequency_idx);
 *     void output_results();
 * 
 * @endcode
 * 
 * This is called before every frequency step to set up a pristine state
 * for the cache variables.
 * 
 * @code
 *     void setup_quadrature_cache();
 * 
 * @endcode
 * 
 * This function loops over the frequency vector and runs the simulation for
 * each frequency step.
 * 
 * @code
 *     void frequency_sweep();
 * 
 * @endcode
 * 
 * The parameters are stored in this variable.
 * 
 * @code
 *     Parameters<dim> parameters;
 * 
 *     MPI_Comm mpi_communicator;
 * 
 *     parallel::distributed::Triangulation<dim> triangulation;
 * 
 *     QGauss<dim> quadrature_formula;
 * 
 * @endcode
 * 
 * We store the mass and stiffness matrices for each cell this vector.
 * 
 * @code
 *     std::vector<QuadratureCache<dim>> quadrature_cache;
 * 
 * 
 *     FESystem<dim>   fe;
 *     DoFHandler<dim> dof_handler;
 * 
 *     IndexSet locally_owned_dofs;
 *     IndexSet locally_relevant_dofs;
 * 
 *     AffineConstraints<std::complex<double>> constraints;
 * 
 *     LinearAlgebraPETSc::MPI::SparseMatrix system_matrix;
 *     LinearAlgebraPETSc::MPI::Vector       locally_relevant_solution;
 *     LinearAlgebraPETSc::MPI::Vector       system_rhs;
 * 
 * 
 * @endcode
 * 
 * This vector contains the range of frequencies that we are going to
 * simulate.
 * 
 * @code
 *     std::vector<double> frequency;
 * 
 * @endcode
 * 
 * This vector contains the coordinates $(x,y)$ of the points of the
 * measurement probe.
 * 
 * @code
 *     FullMatrix<double> probe_positions;
 * 
 * @endcode
 * 
 * HDF5 datasets to store the frequency and `probe_positions` vectors.
 * 
 * @code
 *     HDF5::DataSet frequency_dataset;
 *     HDF5::DataSet probe_positions_dataset;
 * 
 * @endcode
 * 
 * HDF5 dataset that stores the values of the energy measured by the probe.
 * 
 * @code
 *     HDF5::DataSet displacement;
 * 
 * 
 *     ConditionalOStream pcout;
 *     TimerOutput        computing_timer;
 *   };
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Implementationoftheauxiliaryclasses"></a> 
 * <h3>Implementation of the auxiliary classes</h3>
 * 

 * 
 * 
 * <a name="TheRightHandSideclassimplementation"></a> 
 * <h4>The `RightHandSide` class implementation</h4>
 * 

 * 
 * The constructor reads all the parameters from the HDF5::Group `data` using
 * the HDF5::Group::get_attribute() function.
 * 
 * @code
 *   template <int dim>
 *   RightHandSide<dim>::RightHandSide(HDF5::Group &data)
 *     : Function<dim>(dim)
 *     , data(data)
 *     , max_force_amplitude(data.get_attribute<double>("max_force_amplitude"))
 *     , force_sigma_x(data.get_attribute<double>("force_sigma_x"))
 *     , force_sigma_y(data.get_attribute<double>("force_sigma_y"))
 *     , max_force_width_x(data.get_attribute<double>("max_force_width_x"))
 *     , max_force_width_y(data.get_attribute<double>("max_force_width_y"))
 *     , force_center(Point<dim>(data.get_attribute<double>("force_x_pos"),
 *                               data.get_attribute<double>("force_y_pos")))
 *   {}
 * 
 * @endcode
 * 
 * This function defines the spatial shape of the force vector pulse which
 * takes the form of a Gaussian function
 * @f{align*}
 * F_x &=
 * \left\{
 * \begin{array}{ll}
 * a \exp(- (\frac{(x-b_x)^2 }{ 2 \sigma_x^2}+\frac{(y-b_y)^2 }{ 2
 * \sigma_y^2}))
 * & \text{if}\, x_\textrm{min} <x<x_\textrm{max}\, \text{and}\,
 * y_\textrm{min} <y<y_\textrm{max}  \\ 0 & \text{otherwise},
 * \end{array}
 * \right.\\ F_y &= 0
 * @f}
 * where $a$ is the maximum amplitude that takes the force and $\sigma_x$ and
 * $\sigma_y$ are the standard deviations for the $x$ and $y$ components. Note
 * that the pulse has been cropped to $x_\textrm{min}<x<x_\textrm{max}$ and
 * $y_\textrm{min} <y<y_\textrm{max}$.
 * 
 * @code
 *   template <int dim>
 *   double RightHandSide<dim>::value(const Point<dim> & p,
 *                                    const unsigned int component) const
 *   {
 *     if (component == force_component)
 *       {
 *         if (std::abs(p[0] - force_center[0]) < max_force_width_x / 2 &&
 *             std::abs(p[1] - force_center[1]) < max_force_width_y / 2)
 *           {
 *             return max_force_amplitude *
 *                    std::exp(-(std::pow(p[0] - force_center[0], 2) /
 *                                 (2 * std::pow(force_sigma_x, 2)) +
 *                               std::pow(p[1] - force_center[1], 2) /
 *                                 (2 * std::pow(force_sigma_y, 2))));
 *           }
 *         else
 *           {
 *             return 0;
 *           }
 *       }
 *     else
 *       {
 *         return 0;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ThePMLclassimplementation"></a> 
 * <h4>The `PML` class implementation</h4>
 * 

 * 
 * As before, the constructor reads all the parameters from the HDF5::Group
 * `data` using the HDF5::Group::get_attribute() function. As we have
 * discussed, a quadratic turn-on of the PML has been defined in the jupyter
 * notebook. It is possible to use a linear, cubic or another power degree by
 * changing the parameter `pml_coeff_degree`. The parameters `pml_x` and
 * `pml_y` can be used to turn on and off the `x` and `y` PMLs.
 * 
 * @code
 *   template <int dim>
 *   PML<dim>::PML(HDF5::Group &data)
 *     : Function<dim, std::complex<double>>(dim)
 *     , data(data)
 *     , pml_coeff(data.get_attribute<double>("pml_coeff"))
 *     , pml_coeff_degree(data.get_attribute<int>("pml_coeff_degree"))
 *     , dimension_x(data.get_attribute<double>("dimension_x"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , pml_x(data.get_attribute<bool>("pml_x"))
 *     , pml_y(data.get_attribute<bool>("pml_y"))
 *     , pml_width_x(data.get_attribute<double>("pml_width_x"))
 *     , pml_width_y(data.get_attribute<double>("pml_width_y"))
 *     , a_coeff_x(pml_coeff / std::pow(pml_width_x, pml_coeff_degree))
 *     , a_coeff_y(pml_coeff / std::pow(pml_width_y, pml_coeff_degree))
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * The PML coefficient for the `x` component takes the form
 * $s'_x = a_x x^{\textrm{degree}}$
 * 
 * @code
 *   template <int dim>
 *   std::complex<double> PML<dim>::value(const Point<dim> & p,
 *                                        const unsigned int component) const
 *   {
 *     double calculated_pml_x_coeff = 0;
 *     double calculated_pml_y_coeff = 0;
 * 
 *     if ((component == 0) && pml_x)
 *       {
 *         const double pml_x_start_position = dimension_x / 2 - pml_width_x;
 *         if (std::abs(p[0]) > pml_x_start_position)
 *           {
 *             const double x_prime = std::abs(p[0]) - pml_x_start_position;
 *             calculated_pml_x_coeff =
 *               a_coeff_x * std::pow(x_prime, pml_coeff_degree);
 *           }
 *       }
 * 
 *     if ((component == 1) && pml_y)
 *       {
 *         const double pml_y_start_position = dimension_y / 2 - pml_width_y;
 *         if (std::abs(p[1]) > pml_y_start_position)
 *           {
 *             const double y_prime = std::abs(p[1]) - pml_y_start_position;
 *             calculated_pml_y_coeff =
 *               a_coeff_y * std::pow(y_prime, pml_coeff_degree);
 *           }
 *       }
 * 
 *     return 1. + std::max(calculated_pml_x_coeff, calculated_pml_y_coeff) *
 *                   std::complex<double>(0., 1.);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheRhoclassimplementation"></a> 
 * <h4>The `Rho` class implementation</h4>
 * 

 * 
 * This class is used to define the mass density. As we have explaine before,
 * a phononic superlattice cavity is formed by two
 * [Distributed Reflector](https://en.wikipedia.org/wiki/Band_gap),
 * mirrors and a $\lambda/2$ cavity where $\lambda$ is the acoustic
 * wavelength. Acoustic DBRs are periodic structures where a set of bilayer
 * stacks with contrasting physical properties (sound velocity index) is
 * repeated $N$ times. The change of in the wave velocity is generated by
 * alternating layers with different density.
 * 
 * @code
 *   template <int dim>
 *   Rho<dim>::Rho(HDF5::Group &data)
 *     : Function<dim>(1)
 *     , data(data)
 *     , lambda(data.get_attribute<double>("lambda"))
 *     , mu(data.get_attribute<double>("mu"))
 *     , material_a_rho(data.get_attribute<double>("material_a_rho"))
 *     , material_b_rho(data.get_attribute<double>("material_b_rho"))
 *     , cavity_resonance_frequency(
 *         data.get_attribute<double>("cavity_resonance_frequency"))
 *     , nb_mirror_pairs(data.get_attribute<int>("nb_mirror_pairs"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , grid_level(data.get_attribute<int>("grid_level"))
 *   {
 * @endcode
 * 
 * In order to increase the precision we use
 * [subpixel
 * smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/).
 * 
 * @code
 *     average_rho_width = dimension_y / (std::pow(2.0, grid_level));
 *     data.set_attribute("average_rho_width", average_rho_width);
 *   }
 * 
 * 
 * 
 *   template <int dim>
 *   double Rho<dim>::value(const Point<dim> &p,
 *                          const unsigned int /*component*/) const
 *   {
 * @endcode
 * 
 * The speed of sound is defined by
 * @f[
 * c = \frac{K_e}{\rho}
 * @f]
 * where $K_e$ is the effective elastic constant and $\rho$ the density.
 * Here we consider the case in which the waveguide width is much smaller
 * than the wavelength. In this case it can be shown that for the two
 * dimensional case
 * @f[
 * K_e = 4\mu\frac{\lambda +\mu}{\lambda+2\mu}
 * @f]
 * and for the three dimensional case $K_e$ is equal to the Young's modulus.
 * @f[
 * K_e = \mu\frac{3\lambda +2\mu}{\lambda+\mu}
 * @f]
 * 
 * @code
 *     double elastic_constant;
 *     if (dim == 2)
 *       {
 *         elastic_constant = 4 * mu * (lambda + mu) / (lambda + 2 * mu);
 *       }
 *     else if (dim == 3)
 *       {
 *         elastic_constant = mu * (3 * lambda + 2 * mu) / (lambda + mu);
 *       }
 *     else
 *       {
 *         Assert(false, ExcInternalError());
 *       }
 *     const double material_a_speed_of_sound =
 *       std::sqrt(elastic_constant / material_a_rho);
 *     const double material_a_wavelength =
 *       material_a_speed_of_sound / cavity_resonance_frequency;
 *     const double material_b_speed_of_sound =
 *       std::sqrt(elastic_constant / material_b_rho);
 *     const double material_b_wavelength =
 *       material_b_speed_of_sound / cavity_resonance_frequency;
 * 
 * @endcode
 * 
 * The density $\rho$ takes the following form
 * <img alt="Phononic superlattice cavity"
 * src="https://www.dealii.org/images/steps/developer/step-62.04.svg"
 * height="200" />
 * where the brown color represents material_a and the green color
 * represents material_b.
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_transition_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4);
 *         if (std::abs(p[0]) >=
 *               (layer_transition_center - average_rho_width / 2) &&
 *             std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2))
 *           {
 *             const double coefficient =
 *               (std::abs(p[0]) -
 *                (layer_transition_center - average_rho_width / 2)) /
 *               average_rho_width;
 *             return (1 - coefficient) * material_a_rho +
 *                    coefficient * material_b_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * Here we define the
 * [subpixel
 * smoothing](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/)
 * which improves the precision of the simulation.
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_transition_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 4;
 *         if (std::abs(p[0]) >=
 *               (layer_transition_center - average_rho_width / 2) &&
 *             std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2))
 *           {
 *             const double coefficient =
 *               (std::abs(p[0]) -
 *                (layer_transition_center - average_rho_width / 2)) /
 *               average_rho_width;
 *             return (1 - coefficient) * material_b_rho +
 *                    coefficient * material_a_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * then the cavity
 * 
 * @code
 *     if (std::abs(p[0]) <= material_a_wavelength / 2)
 *       {
 *         return material_a_rho;
 *       }
 * 
 * @endcode
 * 
 * the material_a layers
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 4 + material_a_wavelength / 8;
 *         const double layer_width = material_a_wavelength / 4;
 *         if (std::abs(p[0]) >= (layer_center - layer_width / 2) &&
 *             std::abs(p[0]) <= (layer_center + layer_width / 2))
 *           {
 *             return material_a_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * the material_b layers
 * 
 * @code
 *     for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++)
 *       {
 *         const double layer_center =
 *           material_a_wavelength / 2 +
 *           idx * (material_b_wavelength / 4 + material_a_wavelength / 4) +
 *           material_b_wavelength / 8;
 *         const double layer_width = material_b_wavelength / 4;
 *         if (std::abs(p[0]) >= (layer_center - layer_width / 2) &&
 *             std::abs(p[0]) <= (layer_center + layer_width / 2))
 *           {
 *             return material_b_rho;
 *           }
 *       }
 * 
 * @endcode
 * 
 * and finally the default is material_a.
 * 
 * @code
 *     return material_a_rho;
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheParametersclassimplementation"></a> 
 * <h4>The `Parameters` class implementation</h4>
 * 

 * 
 * The constructor reads all the parameters from the HDF5::Group `data` using
 * the HDF5::Group::get_attribute() function.
 * 
 * @code
 *   template <int dim>
 *   Parameters<dim>::Parameters(HDF5::Group &data)
 *     : data(data)
 *     , simulation_name(data.get_attribute<std::string>("simulation_name"))
 *     , save_vtu_files(data.get_attribute<bool>("save_vtu_files"))
 *     , start_frequency(data.get_attribute<double>("start_frequency"))
 *     , stop_frequency(data.get_attribute<double>("stop_frequency"))
 *     , nb_frequency_points(data.get_attribute<int>("nb_frequency_points"))
 *     , lambda(data.get_attribute<double>("lambda"))
 *     , mu(data.get_attribute<double>("mu"))
 *     , dimension_x(data.get_attribute<double>("dimension_x"))
 *     , dimension_y(data.get_attribute<double>("dimension_y"))
 *     , nb_probe_points(data.get_attribute<int>("nb_probe_points"))
 *     , grid_level(data.get_attribute<int>("grid_level"))
 *     , probe_start_point(data.get_attribute<double>("probe_pos_x"),
 *                         data.get_attribute<double>("probe_pos_y") -
 *                           data.get_attribute<double>("probe_width_y") / 2)
 *     , probe_stop_point(data.get_attribute<double>("probe_pos_x"),
 *                        data.get_attribute<double>("probe_pos_y") +
 *                          data.get_attribute<double>("probe_width_y") / 2)
 *     , right_hand_side(data)
 *     , pml(data)
 *     , rho(data)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheQuadratureCacheclassimplementation"></a> 
 * <h4>The `QuadratureCache` class implementation</h4>
 * 

 * 
 * We need to reserve enough space for the mass and stiffness matrices and the
 * right hand side vector.
 * 
 * @code
 *   template <int dim>
 *   QuadratureCache<dim>::QuadratureCache(const unsigned int dofs_per_cell)
 *     : dofs_per_cell(dofs_per_cell)
 *     , mass_coefficient(dofs_per_cell, dofs_per_cell)
 *     , stiffness_coefficient(dofs_per_cell, dofs_per_cell)
 *     , right_hand_side(dofs_per_cell)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationoftheElasticWaveclass"></a> 
 * <h3>Implementation of the `ElasticWave` class</h3>
 * 

 * 
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * This is very similar to the constructor of step-40. In addition we create
 * the HDF5 datasets `frequency_dataset`, `position_dataset` and
 * `displacement`. Note the use of the `template` keyword for the creation of
 * the HDF5 datasets. It is a C++ requirement to use the `template` keyword in
 * order to treat `create_dataset` as a dependent template name.
 * 
 * @code
 *   template <int dim>
 *   ElasticWave<dim>::ElasticWave(const Parameters<dim> &parameters)
 *     : parameters(parameters)
 *     , mpi_communicator(MPI_COMM_WORLD)
 *     , triangulation(mpi_communicator,
 *                     typename Triangulation<dim>::MeshSmoothing(
 *                       Triangulation<dim>::smoothing_on_refinement |
 *                       Triangulation<dim>::smoothing_on_coarsening))
 *     , quadrature_formula(2)
 *     , fe(FE_Q<dim>(1), dim)
 *     , dof_handler(triangulation)
 *     , frequency(parameters.nb_frequency_points)
 *     , probe_positions(parameters.nb_probe_points, dim)
 *     , frequency_dataset(parameters.data.template create_dataset<double>(
 *         "frequency",
 *         std::vector<hsize_t>{parameters.nb_frequency_points}))
 *     , probe_positions_dataset(parameters.data.template create_dataset<double>(
 *         "position",
 *         std::vector<hsize_t>{parameters.nb_probe_points, dim}))
 *     , displacement(
 *         parameters.data.template create_dataset<std::complex<double>>(
 *           "displacement",
 *           std::vector<hsize_t>{parameters.nb_probe_points,
 *                                parameters.nb_frequency_points}))
 *     , pcout(std::cout,
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
 *     , computing_timer(mpi_communicator,
 *                       pcout,
 *                       TimerOutput::summary,
 *                       TimerOutput::wall_times)
 *   {}
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesetup_system"></a> 
 * <h4>ElasticWave::setup_system</h4>
 * 

 * 
 * There is nothing new in this function, the only difference with step-40 is
 * that we don't have to apply boundary conditions because we use the PMLs to
 * truncate the domain.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::setup_system()
 *   {
 *     TimerOutput::Scope t(computing_timer, "setup");
 * 
 *     dof_handler.distribute_dofs(fe);
 * 
 *     locally_owned_dofs = dof_handler.locally_owned_dofs();
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
 * 
 *     locally_relevant_solution.reinit(locally_owned_dofs,
 *                                      locally_relevant_dofs,
 *                                      mpi_communicator);
 * 
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator);
 * 
 *     constraints.clear();
 *     constraints.reinit(locally_relevant_dofs);
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints);
 * 
 *     constraints.close();
 * 
 *     DynamicSparsityPattern dsp(locally_relevant_dofs);
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
 *     SparsityTools::distribute_sparsity_pattern(dsp,
 *                                                locally_owned_dofs,
 *                                                mpi_communicator,
 *                                                locally_relevant_dofs);
 * 
 *     system_matrix.reinit(locally_owned_dofs,
 *                          locally_owned_dofs,
 *                          dsp,
 *                          mpi_communicator);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveassemble_system"></a> 
 * <h4>ElasticWave::assemble_system</h4>
 * 

 * 
 * This function is also very similar to step-40, though there are notable
 * differences. We assemble the system for each frequency/omega step. In the
 * first step we set `calculate_quadrature_data = True` and we calculate the
 * mass and stiffness matrices and the right hand side vector. In the
 * subsequent steps we will use that data to accelerate the calculation.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::assemble_system(const double omega,
 *                                          const bool   calculate_quadrature_data)
 *   {
 *     TimerOutput::Scope t(computing_timer, "assembly");
 * 
 *     FEValues<dim>      fe_values(fe,
 *                             quadrature_formula,
 *                             update_values | update_gradients |
 *                               update_quadrature_points | update_JxW_values);
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
 *     const unsigned int n_q_points    = quadrature_formula.size();
 * 
 *     FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell, dofs_per_cell);
 *     Vector<std::complex<double>>     cell_rhs(dofs_per_cell);
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 * 
 * @endcode
 * 
 * Here we store the value of the right hand side, rho and the PML.
 * 
 * @code
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));
 *     std::vector<double>         rho_values(n_q_points);
 *     std::vector<Vector<std::complex<double>>> pml_values(
 *       n_q_points, Vector<std::complex<double>>(dim));
 * 
 * @endcode
 * 
 * We calculate the stiffness tensor for the $\lambda$ and $\mu$ that have
 * been defined in the jupyter notebook. Note that contrary to $\rho$ the
 * stiffness is constant among for the whole domain.
 * 
 * @code
 *     const SymmetricTensor<4, dim> stiffness_tensor =
 *       get_stiffness_tensor<dim>(parameters.lambda, parameters.mu);
 * 
 * @endcode
 * 
 * We use the same method of step-20 for vector-valued problems.
 * 
 * @code
 *     const FEValuesExtractors::Vector displacement(0);
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell_matrix = 0;
 *           cell_rhs    = 0;
 * 
 * @endcode
 * 
 * We have to calculate the values of the right hand side, rho and
 * the PML only if we are going to calculate the mass and the
 * stiffness matrices. Otherwise we can skip this calculation which
 * considerably reduces the total calculation time.
 * 
 * @code
 *           if (calculate_quadrature_data)
 *             {
 *               fe_values.reinit(cell);
 * 
 *               parameters.right_hand_side.vector_value_list(
 *                 fe_values.get_quadrature_points(), rhs_values);
 *               parameters.rho.value_list(fe_values.get_quadrature_points(),
 *                                         rho_values);
 *               parameters.pml.vector_value_list(
 *                 fe_values.get_quadrature_points(), pml_values);
 *             }
 * 
 * @endcode
 * 
 * We have done this in step-18. Get a pointer to the quadrature
 * cache data local to the present cell, and, as a defensive
 * measure, make sure that this pointer is within the bounds of the
 * global array:
 * 
 * @code
 *           QuadratureCache<dim> *local_quadrature_points_data =
 *             reinterpret_cast<QuadratureCache<dim> *>(cell->user_pointer());
 *           Assert(local_quadrature_points_data >= &quadrature_cache.front(),
 *                  ExcInternalError());
 *           Assert(local_quadrature_points_data <= &quadrature_cache.back(),
 *                  ExcInternalError());
 *           for (unsigned int q = 0; q < n_q_points; ++q)
 *             {
 * @endcode
 * 
 * The quadrature_data variable is used to store the mass and
 * stiffness matrices, the right hand side vector and the value
 * of `JxW`.
 * 
 * @code
 *               QuadratureCache<dim> &quadrature_data =
 *                 local_quadrature_points_data[q];
 * 
 * @endcode
 * 
 * Below we declare the force vector and the parameters of the
 * PML $s$ and $\xi$.
 * 
 * @code
 *               Tensor<1, dim>                       force;
 *               Tensor<1, dim, std::complex<double>> s;
 *               std::complex<double>                 xi(1, 0);
 * 
 * @endcode
 * 
 * The following block is calculated only in the first frequency
 * step.
 * 
 * @code
 *               if (calculate_quadrature_data)
 *                 {
 * @endcode
 * 
 * Store the value of `JxW`.
 * 
 * @code
 *                   quadrature_data.JxW = fe_values.JxW(q);
 * 
 *                   for (unsigned int component = 0; component < dim; ++component)
 *                     {
 * @endcode
 * 
 * Convert vectors to tensors and calculate xi
 * 
 * @code
 *                       force[component] = rhs_values[q][component];
 *                       s[component]     = pml_values[q][component];
 *                       xi *= s[component];
 *                     }
 * 
 * @endcode
 * 
 * Here we calculate the $\alpha_{mnkl}$ and $\beta_{mnkl}$
 * tensors.
 * 
 * @code
 *                   Tensor<4, dim, std::complex<double>> alpha;
 *                   Tensor<4, dim, std::complex<double>> beta;
 *                   for (unsigned int m = 0; m < dim; ++m)
 *                     for (unsigned int n = 0; n < dim; ++n)
 *                       for (unsigned int k = 0; k < dim; ++k)
 *                         for (unsigned int l = 0; l < dim; ++l)
 *                           {
 *                             alpha[m][n][k][l] = xi *
 *                                                 stiffness_tensor[m][n][k][l] /
 *                                                 (2.0 * s[n] * s[k]);
 *                             beta[m][n][k][l] = xi *
 *                                                stiffness_tensor[m][n][k][l] /
 *                                                (2.0 * s[n] * s[l]);
 *                           }
 * 
 *                   for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                     {
 *                       const Tensor<1, dim> phi_i =
 *                         fe_values[displacement].value(i, q);
 *                       const Tensor<2, dim> grad_phi_i =
 *                         fe_values[displacement].gradient(i, q);
 * 
 *                       for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                         {
 *                           const Tensor<1, dim> phi_j =
 *                             fe_values[displacement].value(j, q);
 *                           const Tensor<2, dim> grad_phi_j =
 *                             fe_values[displacement].gradient(j, q);
 * 
 * @endcode
 * 
 * calculate the values of the mass matrix.
 * 
 * @code
 *                           quadrature_data.mass_coefficient[i][j] =
 *                             rho_values[q] * xi * phi_i * phi_j;
 * 
 * @endcode
 * 
 * Loop over the $mnkl$ indices of the stiffness
 * tensor.
 * 
 * @code
 *                           std::complex<double> stiffness_coefficient = 0;
 *                           for (unsigned int m = 0; m < dim; ++m)
 *                             for (unsigned int n = 0; n < dim; ++n)
 *                               for (unsigned int k = 0; k < dim; ++k)
 *                                 for (unsigned int l = 0; l < dim; ++l)
 *                                   {
 * @endcode
 * 
 * Here we calculate the stiffness matrix.
 * Note that the stiffness matrix is not
 * symmetric because of the PMLs. We use the
 * gradient function (see the
 * [documentation](https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html))
 * which is a <code>Tensor@<2,dim@></code>.
 * The matrix $G_{ij}$ consists of entries
 * @f[
 * G_{ij}=
 * \frac{\partial\phi_i}{\partial x_j}
 * =\partial_j \phi_i
 * @f]
 * Note the position of the indices $i$ and
 * $j$ and the notation that we use in this
 * tutorial: $\partial_j\phi_i$. As the
 * stiffness tensor is not symmetric, it is
 * very easy to make a mistake.
 * 
 * @code
 *                                     stiffness_coefficient +=
 *                                       grad_phi_i[m][n] *
 *                                       (alpha[m][n][k][l] * grad_phi_j[l][k] +
 *                                        beta[m][n][k][l] * grad_phi_j[k][l]);
 *                                   }
 * 
 * @endcode
 * 
 * We save the value of the stiffness matrix in
 * quadrature_data
 * 
 * @code
 *                           quadrature_data.stiffness_coefficient[i][j] =
 *                             stiffness_coefficient;
 *                         }
 * 
 * @endcode
 * 
 * and the value of the right hand side in
 * quadrature_data.
 * 
 * @code
 *                       quadrature_data.right_hand_side[i] =
 *                         phi_i * force * fe_values.JxW(q);
 *                     }
 *                 }
 * 
 * @endcode
 * 
 * We loop again over the degrees of freedom of the cells to
 * calculate the system matrix. These loops are really quick
 * because we have already calculated the stiffness and mass
 * matrices, only the value of $\omega$ changes.
 * 
 * @code
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i)
 *                 {
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j)
 *                     {
 *                       std::complex<double> matrix_sum = 0;
 *                       matrix_sum += -std::pow(omega, 2) *
 *                                     quadrature_data.mass_coefficient[i][j];
 *                       matrix_sum += quadrature_data.stiffness_coefficient[i][j];
 *                       cell_matrix(i, j) += matrix_sum * quadrature_data.JxW;
 *                     }
 *                   cell_rhs(i) += quadrature_data.right_hand_side[i];
 *                 }
 *             }
 *           cell->get_dof_indices(local_dof_indices);
 *           constraints.distribute_local_to_global(cell_matrix,
 *                                                  cell_rhs,
 *                                                  local_dof_indices,
 *                                                  system_matrix,
 *                                                  system_rhs);
 *         }
 * 
 *     system_matrix.compress(VectorOperation::add);
 *     system_rhs.compress(VectorOperation::add);
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesolve"></a> 
 * <h4>ElasticWave::solve</h4>
 * 

 * 
 * This is even more simple than in step-40. We use the parallel direct solver
 * MUMPS which requires less options than an iterative solver. The drawback is
 * that it does not scale very well. It is not straightforward to solve the
 * Helmholtz equation with an iterative solver. The shifted Laplacian
 * multigrid method is a well known approach to precondition this system, but
 * this is beyond the scope of this tutorial.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::solve()
 *   {
 *     TimerOutput::Scope              t(computing_timer, "solve");
 *     LinearAlgebraPETSc::MPI::Vector completely_distributed_solution(
 *       locally_owned_dofs, mpi_communicator);
 * 
 *     SolverControl                    solver_control;
 *     PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator);
 *     solver.solve(system_matrix, completely_distributed_solution, system_rhs);
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations."
 *           << std::endl;
 *     constraints.distribute(completely_distributed_solution);
 *     locally_relevant_solution = completely_distributed_solution;
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveinitialize_position_vector"></a> 
 * <h4>ElasticWave::initialize_position_vector</h4>
 * 

 * 
 * We use this function to calculate the values of the position vector.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::initialize_probe_positions_vector()
 *   {
 *     for (unsigned int position_idx = 0;
 *          position_idx < parameters.nb_probe_points;
 *          ++position_idx)
 *       {
 * @endcode
 * 
 * Because of the way the operator + and - are overloaded to subtract
 * two points, the following has to be done:
 * `Point_b<dim> + (-Point_a<dim>)`
 * 
 * @code
 *         const Point<dim> p =
 *           (position_idx / ((double)(parameters.nb_probe_points - 1))) *
 *             (parameters.probe_stop_point + (-parameters.probe_start_point)) +
 *           parameters.probe_start_point;
 *         probe_positions[position_idx][0] = p[0];
 *         probe_positions[position_idx][1] = p[1];
 *         if (dim == 3)
 *           {
 *             probe_positions[position_idx][2] = p[2];
 *           }
 *       }
 *   }
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavestore_frequency_step_data"></a> 
 * <h4>ElasticWave::store_frequency_step_data</h4>
 * 

 * 
 * This function stores in the HDF5 file the measured energy by the probe.
 * 
 * @code
 *   template <int dim>
 *   void
 *   ElasticWave<dim>::store_frequency_step_data(const unsigned int frequency_idx)
 *   {
 *     TimerOutput::Scope t(computing_timer, "store_frequency_step_data");
 * 
 * @endcode
 * 
 * We store the displacement in the $x$ direction; the displacement in the
 * $y$ direction is negligible.
 * 
 * @code
 *     const unsigned int probe_displacement_component = 0;
 * 
 * @endcode
 * 
 * The vector coordinates contains the coordinates in the HDF5 file of the
 * points of the probe that are located in locally owned cells. The vector
 * displacement_data contains the value of the displacement at these points.
 * 
 * @code
 *     std::vector<hsize_t>              coordinates;
 *     std::vector<std::complex<double>> displacement_data;
 * 
 *     const auto &mapping = get_default_linear_mapping(triangulation);
 *     GridTools::Cache<dim, dim> cache(triangulation, mapping);
 *     typename Triangulation<dim, dim>::active_cell_iterator cell_hint{};
 *     std::vector<bool>                                      marked_vertices = {};
 *     const double                                           tolerance = 1.e-10;
 * 
 *     for (unsigned int position_idx = 0;
 *          position_idx < parameters.nb_probe_points;
 *          ++position_idx)
 *       {
 *         Point<dim> point;
 *         for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *           {
 *             point[dim_idx] = probe_positions[position_idx][dim_idx];
 *           }
 *         bool point_in_locally_owned_cell = false;
 *         {
 *           auto cell_and_ref_point = GridTools::find_active_cell_around_point(
 *             cache, point, cell_hint, marked_vertices, tolerance);
 *           if (cell_and_ref_point.first.state() == IteratorState::valid)
 *             {
 *               cell_hint = cell_and_ref_point.first;
 *               point_in_locally_owned_cell =
 *                 cell_and_ref_point.first->is_locally_owned();
 *             }
 *         }
 *         if (point_in_locally_owned_cell)
 *           {
 * @endcode
 * 
 * Then we can store the values of the displacement in the points of
 * the probe in `displacement_data`.
 * 
 * @code
 *             Vector<std::complex<double>> tmp_vector(dim);
 *             VectorTools::point_value(dof_handler,
 *                                      locally_relevant_solution,
 *                                      point,
 *                                      tmp_vector);
 *             coordinates.emplace_back(position_idx);
 *             coordinates.emplace_back(frequency_idx);
 *             displacement_data.emplace_back(
 *               tmp_vector(probe_displacement_component));
 *           }
 *       }
 * 
 * @endcode
 * 
 * We write the displacement data in the HDF5 file. The call
 * HDF5::DataSet::write_selection() is MPI collective which means that all
 * the processes have to participate.
 * 
 * @code
 *     if (coordinates.size() > 0)
 *       {
 *         displacement.write_selection(displacement_data, coordinates);
 *       }
 * @endcode
 * 
 * Therefore even if the process has no data to write it has to participate
 * in the collective call. For this we can use HDF5::DataSet::write_none().
 * Note that we have to specify the data type, in this case
 * `std::complex<double>`.
 * 
 * @code
 *     else
 *       {
 *         displacement.write_none<std::complex<double>>();
 *       }
 * 
 * @endcode
 * 
 * If the variable `save_vtu_files` in the input file equals `True` then all
 * the data will be saved as vtu. The procedure to write `vtu` files has
 * been described in step-40.
 * 
 * @code
 *     if (parameters.save_vtu_files)
 *       {
 *         std::vector<std::string> solution_names(dim, "displacement");
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *           interpretation(
 *             dim, DataComponentInterpretation::component_is_part_of_vector);
 * 
 *         DataOut<dim> data_out;
 *         data_out.add_data_vector(dof_handler,
 *                                  locally_relevant_solution,
 *                                  solution_names,
 *                                  interpretation);
 *         Vector<float> subdomain(triangulation.n_active_cells());
 *         for (unsigned int i = 0; i < subdomain.size(); ++i)
 *           subdomain(i) = triangulation.locally_owned_subdomain();
 *         data_out.add_data_vector(subdomain, "subdomain");
 * 
 *         std::vector<Vector<double>> force(
 *           dim, Vector<double>(triangulation.n_active_cells()));
 *         std::vector<Vector<double>> pml(
 *           dim, Vector<double>(triangulation.n_active_cells()));
 *         Vector<double> rho(triangulation.n_active_cells());
 * 
 *         for (auto &cell : triangulation.active_cell_iterators())
 *           {
 *             if (cell->is_locally_owned())
 *               {
 *                 for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *                   {
 *                     force[dim_idx](cell->active_cell_index()) =
 *                       parameters.right_hand_side.value(cell->center(), dim_idx);
 *                     pml[dim_idx](cell->active_cell_index()) =
 *                       parameters.pml.value(cell->center(), dim_idx).imag();
 *                   }
 *                 rho(cell->active_cell_index()) =
 *                   parameters.rho.value(cell->center());
 *               }
 * @endcode
 * 
 * And on the cells that we are not interested in, set the
 * respective value to a bogus value in order to make sure that if
 * we were somehow wrong about our assumption we would find out by
 * looking at the graphical output:
 * 
 * @code
 *             else
 *               {
 *                 for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *                   {
 *                     force[dim_idx](cell->active_cell_index()) = -1e+20;
 *                     pml[dim_idx](cell->active_cell_index())   = -1e+20;
 *                   }
 *                 rho(cell->active_cell_index()) = -1e+20;
 *               }
 *           }
 * 
 *         for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx)
 *           {
 *             data_out.add_data_vector(force[dim_idx],
 *                                      "force_" + std::to_string(dim_idx));
 *             data_out.add_data_vector(pml[dim_idx],
 *                                      "pml_" + std::to_string(dim_idx));
 *           }
 *         data_out.add_data_vector(rho, "rho");
 * 
 *         data_out.build_patches();
 * 
 *         std::stringstream  frequency_idx_stream;
 *         const unsigned int nb_number_positions =
 *           ((unsigned int)std::log10(parameters.nb_frequency_points)) + 1;
 *         frequency_idx_stream << std::setw(nb_number_positions)
 *                              << std::setfill('0') << frequency_idx;
 *         std::string filename = (parameters.simulation_name + "_" +
 *                                 frequency_idx_stream.str() + ".vtu");
 *         data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator);
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaveoutput_results"></a> 
 * <h4>ElasticWave::output_results</h4>
 * 

 * 
 * This function writes the datasets that have not already been written.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::output_results()
 *   {
 * @endcode
 * 
 * The vectors `frequency` and `position` are the same for all the
 * processes. Therefore any of the processes can write the corresponding
 * `datasets`. Because the call HDF5::DataSet::write is MPI collective, the
 * rest of the processes will have to call HDF5::DataSet::write_none.
 * 
 * @code
 *     if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *       {
 *         frequency_dataset.write(frequency);
 *         probe_positions_dataset.write(probe_positions);
 *       }
 *     else
 *       {
 *         frequency_dataset.write_none<double>();
 *         probe_positions_dataset.write_none<double>();
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavesetup_quadrature_cache"></a> 
 * <h4>ElasticWave::setup_quadrature_cache</h4>
 * 

 * 
 * We use this function at the beginning of our computations to set up initial
 * values of the cache variables. This function has been described in step-18.
 * There are no differences with the function of step-18.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::setup_quadrature_cache()
 *   {
 *     triangulation.clear_user_data();
 * 
 *     {
 *       std::vector<QuadratureCache<dim>> tmp;
 *       quadrature_cache.swap(tmp);
 *     }
 * 
 *     quadrature_cache.resize(triangulation.n_locally_owned_active_cells() *
 *                               quadrature_formula.size(),
 *                             QuadratureCache<dim>(fe.n_dofs_per_cell()));
 *     unsigned int cache_index = 0;
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       if (cell->is_locally_owned())
 *         {
 *           cell->set_user_pointer(&quadrature_cache[cache_index]);
 *           cache_index += quadrature_formula.size();
 *         }
 *     Assert(cache_index == quadrature_cache.size(), ExcInternalError());
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWavefrequency_sweep"></a> 
 * <h4>ElasticWave::frequency_sweep</h4>
 * 

 * 
 * For clarity we divide the function `run` of step-40 into the functions
 * `run` and `frequency_sweep`. In the function `frequency_sweep` we place the
 * iteration over the frequency vector.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::frequency_sweep()
 *   {
 *     for (unsigned int frequency_idx = 0;
 *          frequency_idx < parameters.nb_frequency_points;
 *          ++frequency_idx)
 *       {
 *         pcout << parameters.simulation_name + " frequency idx: "
 *               << frequency_idx << '/' << parameters.nb_frequency_points - 1
 *               << std::endl;
 * 
 * 
 * 
 *         setup_system();
 *         if (frequency_idx == 0)
 *           {
 *             pcout << "   Number of active cells :       "
 *                   << triangulation.n_active_cells() << std::endl;
 *             pcout << "   Number of degrees of freedom : "
 *                   << dof_handler.n_dofs() << std::endl;
 *           }
 * 
 *         if (frequency_idx == 0)
 *           {
 * @endcode
 * 
 * Write the simulation parameters only once
 * 
 * @code
 *             parameters.data.set_attribute("active_cells",
 *                                           triangulation.n_active_cells());
 *             parameters.data.set_attribute("degrees_of_freedom",
 *                                           dof_handler.n_dofs());
 *           }
 * 
 * @endcode
 * 
 * We calculate the frequency and omega values for this particular step.
 * 
 * @code
 *         const double current_loop_frequency =
 *           (parameters.start_frequency +
 *            frequency_idx *
 *              (parameters.stop_frequency - parameters.start_frequency) /
 *              (parameters.nb_frequency_points - 1));
 *         const double current_loop_omega =
 *           2 * numbers::PI * current_loop_frequency;
 * 
 * @endcode
 * 
 * In the first frequency step we calculate the mass and stiffness
 * matrices and the right hand side. In the subsequent frequency steps
 * we will use those values. This improves considerably the calculation
 * time.
 * 
 * @code
 *         assemble_system(current_loop_omega,
 *                         (frequency_idx == 0) ? true : false);
 *         solve();
 * 
 *         frequency[frequency_idx] = current_loop_frequency;
 *         store_frequency_step_data(frequency_idx);
 * 
 *         computing_timer.print_summary();
 *         computing_timer.reset();
 *         pcout << std::endl;
 *       }
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticWaverun"></a> 
 * <h4>ElasticWave::run</h4>
 * 

 * 
 * This function is very similar to the one in step-40.
 * 
 * @code
 *   template <int dim>
 *   void ElasticWave<dim>::run()
 *   {
 * #ifdef DEBUG
 *     pcout << "Debug mode" << std::endl;
 * #else
 *     pcout << "Release mode" << std::endl;
 * #endif
 * 
 *     {
 *       Point<dim> p1;
 *       p1(0) = -parameters.dimension_x / 2;
 *       p1(1) = -parameters.dimension_y / 2;
 *       if (dim == 3)
 *         {
 *           p1(2) = -parameters.dimension_y / 2;
 *         }
 *       Point<dim> p2;
 *       p2(0) = parameters.dimension_x / 2;
 *       p2(1) = parameters.dimension_y / 2;
 *       if (dim == 3)
 *         {
 *           p2(2) = parameters.dimension_y / 2;
 *         }
 *       std::vector<unsigned int> divisions(dim);
 *       divisions[0] = int(parameters.dimension_x / parameters.dimension_y);
 *       divisions[1] = 1;
 *       if (dim == 3)
 *         {
 *           divisions[2] = 1;
 *         }
 *       GridGenerator::subdivided_hyper_rectangle(triangulation,
 *                                                 divisions,
 *                                                 p1,
 *                                                 p2);
 *     }
 * 
 *     triangulation.refine_global(parameters.grid_level);
 * 
 *     setup_quadrature_cache();
 * 
 *     initialize_probe_positions_vector();
 * 
 *     frequency_sweep();
 * 
 *     output_results();
 *   }
 * } // namespace step62
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h4>The main function</h4>
 * 

 * 
 * The main function is very similar to the one in step-40.
 * 
 * @code
 * int main(int argc, char *argv[])
 * {
 *   try
 *     {
 *       using namespace dealii;
 *       const unsigned int dim = 2;
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
 * 
 *       HDF5::File data_file("results.h5",
 *                            HDF5::File::FileAccessMode::create,
 *                            MPI_COMM_WORLD);
 *       auto       data = data_file.create_group("data");
 * 
 * @endcode
 * 
 * Each of the simulations (displacement and calibration) is stored in a
 * separate HDF5 group:
 * 
 * @code
 *       const std::vector<std::string> group_names = {"displacement",
 *                                                     "calibration"};
 *       for (auto group_name : group_names)
 *         {
 * @endcode
 * 
 * For each of these two group names, we now create the group and put
 * attributes into these groups.
 * Specifically, these are:
 * - The dimensions of the waveguide (in $x$ and $y$ directions)
 * - The position of the probe (in $x$ and $y$ directions)
 * - The number of points in the probe
 * - The global refinement level
 * - The cavity resonance frequency
 * - The number of mirror pairs
 * - The material properties
 * - The force parameters
 * - The PML parameters
 * - The frequency parameters
 * 

 * 
 * 
 * @code
 *           auto group = data.create_group(group_name);
 * 
 *           group.set_attribute<double>("dimension_x", 2e-5);
 *           group.set_attribute<double>("dimension_y", 2e-8);
 *           group.set_attribute<double>("probe_pos_x", 8e-6);
 *           group.set_attribute<double>("probe_pos_y", 0);
 *           group.set_attribute<double>("probe_width_y", 2e-08);
 *           group.set_attribute<unsigned int>("nb_probe_points", 5);
 *           group.set_attribute<unsigned int>("grid_level", 1);
 *           group.set_attribute<double>("cavity_resonance_frequency", 20e9);
 *           group.set_attribute<unsigned int>("nb_mirror_pairs", 15);
 * 
 *           group.set_attribute<double>("poissons_ratio", 0.27);
 *           group.set_attribute<double>("youngs_modulus", 270000000000.0);
 *           group.set_attribute<double>("material_a_rho", 3200);
 * 
 *           if (group_name == std::string("displacement"))
 *             group.set_attribute<double>("material_b_rho", 2000);
 *           else
 *             group.set_attribute<double>("material_b_rho", 3200);
 * 
 *           group.set_attribute(
 *             "lambda",
 *             group.get_attribute<double>("youngs_modulus") *
 *               group.get_attribute<double>("poissons_ratio") /
 *               ((1 + group.get_attribute<double>("poissons_ratio")) *
 *                (1 - 2 * group.get_attribute<double>("poissons_ratio"))));
 *           group.set_attribute("mu",
 *                               group.get_attribute<double>("youngs_modulus") /
 *                                 (2 * (1 + group.get_attribute<double>(
 *                                             "poissons_ratio"))));
 * 
 *           group.set_attribute<double>("max_force_amplitude", 1e26);
 *           group.set_attribute<double>("force_sigma_x", 1e-7);
 *           group.set_attribute<double>("force_sigma_y", 1);
 *           group.set_attribute<double>("max_force_width_x", 3e-7);
 *           group.set_attribute<double>("max_force_width_y", 2e-8);
 *           group.set_attribute<double>("force_x_pos", -8e-6);
 *           group.set_attribute<double>("force_y_pos", 0);
 * 
 *           group.set_attribute<bool>("pml_x", true);
 *           group.set_attribute<bool>("pml_y", false);
 *           group.set_attribute<double>("pml_width_x", 1.8e-6);
 *           group.set_attribute<double>("pml_width_y", 5e-7);
 *           group.set_attribute<double>("pml_coeff", 1.6);
 *           group.set_attribute<unsigned int>("pml_coeff_degree", 2);
 * 
 *           group.set_attribute<double>("center_frequency", 20e9);
 *           group.set_attribute<double>("frequency_range", 0.5e9);
 *           group.set_attribute<double>(
 *             "start_frequency",
 *             group.get_attribute<double>("center_frequency") -
 *               group.get_attribute<double>("frequency_range") / 2);
 *           group.set_attribute<double>(
 *             "stop_frequency",
 *             group.get_attribute<double>("center_frequency") +
 *               group.get_attribute<double>("frequency_range") / 2);
 *           group.set_attribute<unsigned int>("nb_frequency_points", 400);
 * 
 *           if (group_name == std::string("displacement"))
 *             group.set_attribute<std::string>(
 *               "simulation_name", std::string("phononic_cavity_displacement"));
 *           else
 *             group.set_attribute<std::string>(
 *               "simulation_name", std::string("phononic_cavity_calibration"));
 * 
 *           group.set_attribute<bool>("save_vtu_files", false);
 *         }
 * 
 *       {
 * @endcode
 * 
 * Displacement simulation. The parameters are read from the
 * displacement HDF5 group and the results are saved in the same HDF5
 * group.
 * 
 * @code
 *         auto                    displacement = data.open_group("displacement");
 *         step62::Parameters<dim> parameters(displacement);
 * 
 *         step62::ElasticWave<dim> elastic_problem(parameters);
 *         elastic_problem.run();
 *       }
 * 
 *       {
 * @endcode
 * 
 * Calibration simulation. The parameters are read from the calibration
 * HDF5 group and the results are saved in the same HDF5 group.
 * 
 * @code
 *         auto                    calibration = data.open_group("calibration");
 *         step62::Parameters<dim> parameters(calibration);
 * 
 *         step62::ElasticWave<dim> elastic_problem(parameters);
 *         elastic_problem.run();
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
examples/step-62/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Resonancefrequencyandbandgap"></a><h3>Resonance frequency and bandgap</h3>


在[jupyter notebook](https://github.com/dealii/dealii/blob/master/example/step-62/step-62.ipynb)中用以下代码分析了结果

@code{.py}
h5_file = h5py.File('results.h5', 'r')
data = h5_file['data']


# Gaussian function that we use to fit the resonance
def resonance_f(freq, freq_m, quality_factor, max_amplitude):
    omega = 2 * constants.pi * freq
    omega_m = 2 * constants.pi * freq_m
    gamma = omega_m / quality_factor
    return max_amplitude * omega_m**2 * gamma**2 / (((omega_m**2 - omega**2)**2 + gamma**2 * omega**2))


frequency = data['displacement']['frequency'][...]
# Average the probe points
displacement = np.mean(data['displacement']['displacement'], axis=0)
calibration_displacement = np.mean(data['calibration']['displacement'], axis=0)
reflection_coefficient = displacement / calibration_displacement
reflectivity = (np.abs(np.mean(data['displacement']['displacement'][...]**2, axis=0))/
                np.abs(np.mean(data['calibration']['displacement'][...]**2, axis=0)))


try:
    x_data = frequency
    y_data = reflectivity
    quality_factor_guess = 1e3
    freq_guess = x_data[np.argmax(y_data)]
    amplitude_guess = np.max(y_data)
    fit_result, covariance = scipy.optimize.curve_fit(resonance_f, x_data, y_data,
                                                      [freq_guess, quality_factor_guess, amplitude_guess])
    freq_m = fit_result[0]
    quality_factor = np.abs(fit_result[1])
    max_amplitude = fit_result[2]
    y_data_fit = resonance_f(x_data, freq_m, quality_factor, max_amplitude)


    fig = plt.figure()
    plt.plot(frequency / 1e9, reflectivity, frequency / 1e9, y_data_fit)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('amplitude^2 (a.u.)')
    plt.title('Transmission\n' + 'freq = ' + "%.7g" % (freq_guess / 1e9) + 'GHz Q = ' + "%.6g" % quality_factor)
except:
    fig = plt.figure()
    plt.plot(frequency / 1e9, reflectivity)
    plt.xlabel('frequency (GHz)')
    plt.ylabel('amplitude^2 (a.u.)')
    plt.title('Transmission')


fig = plt.figure()
plt.plot(frequency / 1e9, np.angle(reflection_coefficient))
plt.xlabel('frequency (GHz)')
plt.ylabel('phase (rad)')
plt.title('Phase (transmission coefficient)\n')


plt.show()
h5_file.close()
@endcode



一个声腔的特点是[共振频率](https://en.wikipedia.org/wiki/Resonance)和[品质因子](https://en.wikipedia.org/wiki/Q_factor)。质量因子等于谐振器中储存的能量与每周期耗散的能量之间的比率，这大约相当于谐振频率与[半满宽度（FWHM）](https://en.wikipedia.org/wiki/Full_width_at_half_maximum)之间的比率。FWHM等于振动功率大于谐振频率的一半的带宽。

@f[
Q = \frac{f_r}{\Delta f} = \frac{\omega_r}{\Delta \omega} =
2 \pi \times \frac{\text{energy stored}}{\text{energy dissipated per cycle}}


@f]



机械共振 $a^2$ 的振幅的平方作为频率的函数有一个高斯形状

@f[
a^2 = a_\textrm{max}^2\frac{\omega^2\Gamma^2}{(\omega_r^2-\omega^2)^2+\Gamma^2\omega^2}


@f]

其中 $f_r = \frac{\omega_r}{2\pi}$ 是共振频率， $\Gamma=\frac{\omega_r}{Q}$ 是耗损率。我们在jupyter笔记本中使用前面的方程式来拟合机械共振。

鉴于我们所选择的参数值，人们可以通过分析来估计谐振频率。事实上，我们在这个程序中得到的结果证实了这一点：声子超晶格空腔在20GHz时表现出机械共振，质量系数为5046。下面的图片显示了在共振频率附近的传输振幅和相位与频率的关系。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.05.png" height="400" />  <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.06.png" height="400" />  。

上面的图片表明，周期性结构有其预期的效果：它实际上只让一个非常特定频率的波通过，而所有其他的波都被反射。当然，这正是人们建造这类设备的目的。但这并不十分容易。在实践中，实际上只有一个 "带隙"，也就是说，该设备只在一定的频率范围内阻止20GHz频率以外的波。事实上，要想知道这个被阻挡的 "间隙 "有多大，我们可以通过输入文件中的适当参数将频率范围扩大到16GHz。然后我们得到以下图像。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.07.png" height="400" /> 

这张图片表明的是，在18到22GHz左右的范围内，确实只有频率为20GHz的波被允许通过，但在这个范围之外，还有很多其他频率可以通过该设备。

<a name="Modeprofile"></a><h3>Mode profile</h3>


我们可以用Paraview或VisIt检查模态轮廓。正如我们所讨论的，在共振时，所有的机械能都被传递，运动的振幅在腔内被放大。可以看出，PML对于截断解决方案是非常有效的。下图显示了共振时的模式轮廓。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.08.png" height="400" /> 

另一方面，在共振之外，所有的机械能都被反射。下面的图片显示了19.75GHz时的轮廓。注意力脉冲和反射波在位置 $x=-8\mu\textrm{m}$ 的干扰。

 <img alt="Phononic superlattice cavity" src="https://www.dealii.org/images/steps/developer/step-62.09.png" height="400" /> 

<a name="Experimentalapplications"></a><h3>Experimental applications</h3>


声波超晶格空腔在[量子光学机械学](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.86.1391)中找到了应用。这里我们介绍了二维超晶格空腔的模拟，但这个代码也可以用来模拟 "现实世界 "的三维设备，如[微柱超晶格空腔](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.060101)，它是研究宏观量子现象的有希望的候选者。微柱超晶格空腔的20GHz模式本质上是一个机械谐波振荡器，与环境隔离得非常好。如果该装置在稀释冰箱中被冷却到20mK，那么该模式就会成为一个宏观的量子谐波振荡器。




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


我们可以不在C++文件中设置参数，而是用一个Python脚本来设置参数，并将其保存在HDF5文件中，我们将使用该文件进行模拟。然后deal.II程序将从HDF5文件中读取参数。

@code{.py}
import numpy as np
import h5py
import matplotlib.pyplot as plt
import subprocess
import scipy.constants as constants
import scipy.optimize


# This considerably reduces the size of the svg data
plt.rcParams['svg.fonttype'] = 'none'


h5_file = h5py.File('results.h5', 'w')
data = h5_file.create_group('data')
displacement = data.create_group('displacement')
calibration = data.create_group('calibration')


# Set the parameters
for group in [displacement, calibration]:
    # Dimensions of the domain
    # The waveguide length is equal to dimension_x
    group.attrs['dimension_x'] = 2e-5
    # The waveguide width is equal to dimension_y
    group.attrs['dimension_y'] = 2e-8


    # Position of the probe that we use to measure the flux
    group.attrs['probe_pos_x']   = 8e-6
    group.attrs['probe_pos_y']   = 0
    group.attrs['probe_width_y'] = 2e-08


    # Number of points in the probe
    group.attrs['nb_probe_points'] = 5


    # Global refinement
    group.attrs['grid_level'] = 1


    # Cavity
    group.attrs['cavity_resonance_frequency'] = 20e9
    group.attrs['nb_mirror_pairs']            = 15


    # Material
    group.attrs['poissons_ratio'] = 0.27
    group.attrs['youngs_modulus'] = 270000000000.0
    group.attrs['material_a_rho'] = 3200
    if group == displacement:
        group.attrs['material_b_rho'] = 2000
    else:
        group.attrs['material_b_rho'] = 3200
    group.attrs['lambda'] = (group.attrs['youngs_modulus'] * group.attrs['poissons_ratio'] /
                           ((1 + group.attrs['poissons_ratio']) *
                           (1 - 2 * group.attrs['poissons_ratio'])))
    group.attrs['mu']= (group.attrs['youngs_modulus'] / (2 * (1 + group.attrs['poissons_ratio'])))


    # Force
    group.attrs['max_force_amplitude'] = 1e26
    group.attrs['force_sigma_x']       = 1e-7
    group.attrs['force_sigma_y']       = 1
    group.attrs['max_force_width_x']   = 3e-7
    group.attrs['max_force_width_y']   = 2e-8
    group.attrs['force_x_pos']         = -8e-6
    group.attrs['force_y_pos']         = 0


    # PML
    group.attrs['pml_x']            = True
    group.attrs['pml_y']            = False
    group.attrs['pml_width_x']      = 1.8e-6
    group.attrs['pml_width_y']      = 5e-7
    group.attrs['pml_coeff']        = 1.6
    group.attrs['pml_coeff_degree'] = 2


    # Frequency sweep
    group.attrs['center_frequency']    = 20e9
    group.attrs['frequency_range']     = 0.5e9
    group.attrs['start_frequency']     = group.attrs['center_frequency'] - group.attrs['frequency_range'] / 2
    group.attrs['stop_frequency']      = group.attrs['center_frequency'] + group.attrs['frequency_range'] / 2
    group.attrs['nb_frequency_points'] = 400


    # Other parameters
    if group == displacement:
        group.attrs['simulation_name'] = 'phononic_cavity_displacement'
    else:
        group.attrs['simulation_name'] = 'phononic_cavity_calibration'
    group.attrs['save_vtu_files'] = False


h5_file.close()
@endcode



为了读取HDF5参数，我们必须使用 HDF5::File::FileAccessMode::open 标志。

@code{.py}
      HDF5::File data_file("results.h5",
                           HDF5::File::FileAccessMode::open,
                           MPI_COMM_WORLD);
      auto       data = data_file.open_group("data");
@endcode




 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-62.cc"
*/
