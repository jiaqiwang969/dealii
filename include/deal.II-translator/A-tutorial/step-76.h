/**
  @page step_76 The step-76 tutorial program  
* This tutorial depends on   step-67  .
* @htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a><a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#SharedmemoryandhybridparallelizationwithMPI30">Shared-memory and hybrid parallelization with MPI-3.0</a><a href="#SharedmemoryandhybridparallelizationwithMPI30">Shared-memory and hybrid parallelization with MPI-3.0</a>
      <ul>
        <li><a href="#Motivation">Motivation</a><a href="#Motivation">Motivation</a>
        <li><a href="#BasicMPI30commands">Basic MPI-3.0 commands</a><a href="#BasicMPI30commands">Basic MPI-3.0 commands</a>
        <li><a href="#MPI30andLinearAlgebradistributedVector">MPI-3.0 and LinearAlgebra::distributed::Vector</a><a href="#MPI30andLinearAlgebradistributedVector">MPI-3.0 and LinearAlgebra::distributed::Vector</a>
        <li><a href="#MPI30andMatrixFree">MPI-3.0 and MatrixFree</a><a href="#MPI30andMatrixFree">MPI-3.0 and MatrixFree</a>
      </ul>
        <li><a href="#Cellcentricloops">Cell-centric loops</a><a href="#Cellcentricloops">Cell-centric loops</a>
      <ul>
        <li><a href="#MotivationFCLvsCCL">Motivation: FCL vs. CCL</a><a href="#MotivationFCLvsCCL">Motivation: FCL vs. CCL</a>
        <li><a href="#CellcentricloopsandMatrixFree">Cell-centric loops and MatrixFree</a><a href="#CellcentricloopsandMatrixFree">Cell-centric loops and MatrixFree</a>
      </ul>
        <li><a href="#ProvidinglambdastoMatrixFreeloops">Providing lambdas to MatrixFree loops</a><a href="#ProvidinglambdastoMatrixFreeloops">Providing lambdas to MatrixFree loops</a>
        <li><a href="#VectorizedArrayType">VectorizedArrayType</a><a href="#VectorizedArrayType">VectorizedArrayType</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a><a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Parametersandutilityfunctions">Parameters and utility functions</a><a href="#Parametersandutilityfunctions">Parameters and utility functions</a>
        <li><a href="#EuleroperatorusingacellcentricloopandMPI30sharedmemory">Euler operator using a cell-centric loop and MPI-3.0 shared memory</a><a href="#EuleroperatorusingacellcentricloopandMPI30sharedmemory">Euler operator using a cell-centric loop and MPI-3.0 shared memory</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a><a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a>
        <li><a href="#BlockGaussSeidellikepreconditioners">Block Gauss-Seidel-like preconditioners</a><a href="#BlockGaussSeidellikepreconditioners">Block Gauss-Seidel-like preconditioners</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a><a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
* 
*   <br>  
* <i>
This program was contributed by Martin Kronbichler, Peter Munch, and David
Schneider. Many of the features shown here have been added to deal.II during
and for the development of the deal.II-based, efficient, matrix-free
finite-element library for high-dimensional partial differential equations
hyper.deal (see https://github.com/hyperdeal/hyperdeal). For more details and
for applications of the presented features in slightly different contexts
(high-dimensional advection equation and Vlasov-Poisson equations) see the release
paper @cite munch2020hyperdeal.
* 
This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA) and
by the Bavarian government through the project "High-order matrix-free finite
element implementations with hybrid parallelization and improved data locality"
within the KONWIHR program.
</i>
* <a name="Intro"></a><a name="Introduction"></a><h1>Introduction</h1>
* 

* This tutorial program solves the Euler equations of fluid dynamics, using anexplicit time integrator with the matrix-free framework applied to ahigh-order discontinuous Galerkin discretization in space. The numericalapproach used here is identical to that used in   step-67  , however, we utilizedifferent advanced MatrixFree techniques to reach even a higher throughput.
* The two main features of this tutorial are:
* 
*  - the usage of shared-memory features from MPI-3.0 and
* 
*  - the usage of cell-centric loops, which allow to write to the global vector only  once and, therefore, are ideal for the usage of shared memory.
* Further topics we discuss in this tutorial are the usage and benefits of thetemplate argument VectorizedArrayType (instead of simply usingVectorizedArray<Number>) as well as the possibility to pass lambdas toMatrixFree loops.
* For details on the numerics, we refer to the documentation of   step-67  . Weconcentrate here only on the key differences.
* <a name="SharedmemoryandhybridparallelizationwithMPI30"></a><h3>Shared-memory and hybrid parallelization with MPI-3.0</h3>
* 

* <a name="Motivation"></a><h4>Motivation</h4>
* 

* There exist many shared-memory libraries that are based on threads like TBB,OpenMP, or TaskFlow. Integrating such libraries into existing MPI programsallows one to use shared memory. However, these libraries come with an overheadfor the programmer, since all parallelizable code sections have to be found andtransformed according to the library used, including the difficulty when somethird-party numerical library, like an iterative solver package, only relies onMPI.
* Considering a purely MPI-parallelized FEM application, one can identify thatthe major time and memory benefit of using shared memory would come fromaccessing the part of the solution vector owned by the processes on the samecompute node without the need to make explicit copies and buffering them.Fur this propose, MPI-3.0 provides shared-memory features based on so-calledwindows, where processes can directly access the data of the neighbors on the sameshared-memory domain.
* <a name="BasicMPI30commands"></a><h4>Basic MPI-3.0 commands</h4>
* 

* A few relevant MPI-3.0 commands are worth discussing in detail.A new MPI communicator   <code>comm_sm</code>  , which consists of processes fromthe communicator   <code>comm</code>   that have access to the same shared memory,can be created via:
* @code
MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &comm_sm);
@endcode
* 
* The following code snippet shows the simplified allocation routines ofshared memory for the value type   <code>T</code>   and the size  <code>local_size</code>  , as well as, how to query pointers to the data belongingto processes in the same shared-memory domain:
* @code
MPI_Win          win;         // window
T              data_this;   // pointer to locally-owned data
std::vector<T> data_others; // pointers to shared data
* 
// configure shared memory
MPI_Info info;
MPI_Info_create(&info);
MPI_Info_set(info, "alloc_shared_noncontig", "true");
* 
// allocate shared memory
MPI_Win_allocate_shared(local_size sizeof(T), sizeof(T), info, comm_sm, &data_this, &win);
* 
// get pointers to the shared data owned by the processes in the same sm domain
data_others.resize(size_sm);
int disp_unit = 0; // displacement size
* 
- an output parameter we don't need right now
MPI_Aint ssize = 0; // window size
* 
- an output parameter we don't need right now
for (int i = 0; i < size_sm; ++i)
  MPI_Win_shared_query(win, i, &ssize, &disp_unit, &data_others[i]);
* 
Assert(data_this == data_others[rank_sm], ExcMessage("Something went wrong!"));
@endcode
* 
* Once the data is not needed anymore, the window has to be freed, which alsofrees the locally-owned data:
* @code
MPI_Win_free(&win);
@endcode
* 
* <a name="MPI30andLinearAlgebradistributedVector"></a><h4>MPI-3.0 and LinearAlgebra::distributed::Vector</h4>
* 

* The commands mentioned in the last section are integrated into  LinearAlgebra::distributed::Vector   and are used to allocate shared memory ifan optional (second) communicator is provided to the reinit()-functions.
* For example, a vector can be set up with a partitioner (containing the globalcommunicator) and a sub-communicator (containing the processes on the samecompute node):
* @code
vec.reinit(partitioner, comm_sm);
@endcode
* 
* Locally owned values and ghost values can be processed as usual. However, nowusers also have read access to the values of the shared-memory neighbors viathe function:
* @code
const std::vector<ArrayView<const Number>> &
LinearAlgebra::distributed::Vector::shared_vector_data() const;
@endcode
* 
* <a name="MPI30andMatrixFree"></a><h4>MPI-3.0 and MatrixFree</h4>
* 

* While   LinearAlgebra::distributed::Vector   provides the option to allocateshared memory and to access the values of shared memory of neighboring processesin a coordinated way, it does not actually exploit the benefits of theusage of shared memory itself.
* The MatrixFree infrastructure, however, does:
* 
*  - On the one hand, within the matrix-free loops   MatrixFree::loop(),      MatrixFree::cell_loop(),   and   MatrixFree::loop_cell_centric(),   only ghost  values that need to be updated   <em>  are  </em>   updated. Ghost values from  shared-memory neighbors can be accessed directly, making buffering, i.e.,  copying of the values into the ghost region of a vector possibly redundant.  To deal with possible race conditions, necessary synchronizations are  performed within MatrixFree. In the case that values have to be buffered,  values are copied directly from the neighboring shared-memory process,  bypassing more expensive MPI operations based on   <code>MPI_ISend</code>   and    <code>MPI_IRecv</code>  .
* 
*  - On the other hand, classes like FEEvaluation and FEFaceEvaluation can read  directly from the shared memory, so buffering the values is indeed  not necessary in certain cases.
* To be able to use the shared-memory capabilities of MatrixFree, MatrixFreehas to be appropriately configured by providing the user-created sub-communicator:
* @code
typename MatrixFree<dim, Number>::AdditionalData additional_data;
* 
// set flags as usual (not shown)
* 
additional_data.communicator_sm = comm_sm;
* 
data.reinit(mapping, dof_handler, constraint, quadrature, additional_data);
@endcode
* 
* 

* <a name="Cellcentricloops"></a><h3>Cell-centric loops</h3>
* 

* <a name="MotivationFCLvsCCL"></a><h4>Motivation: FCL vs. CCL</h4>
* 

* "Face-centric loops" (short FCL) visit cells and faces (inner and boundary ones) inseparate loops. As a consequence, each entity is visited only once and fluxesbetween cells are evaluated only once. How to perform face-centric loopswith the help of   MatrixFree::loop()   by providing three functions (one forthe cell integrals, one for the inner, and one for the boundary faces) hasbeen presented in   step-59   and   step-67  .
* "Cell-centric loops" (short CCL or ECL (for element-centric loops)in the hyper.deal release paper), incontrast, process a cell and in direct succession process all itsfaces (i.e., visit all faces twice). Their benefit has become clear formodern CPU processor architecture in the literature   @cite KronbichlerKormann2019  ,although this kind of loop implies that fluxes have to be computed twice (foreach side of an interior face). CCL has two primary advantages:
* 
*  - On the one hand, entries in the solution vector are written exactly once  back to main memory in the case of CCL, while in the case of FCL at least once  despite of cache-efficient scheduling of cell and face loops-due to cache  capacity misses.
* 
*  - On the other hand, since each entry of the solution vector is accessed exactly  once, no synchronization between threads is needed while accessing the solution  vector in the case of CCL. This absence of race conditions during writing into  the destination vector makes CCL particularly suitable for shared-memory  parallelization.
* One should also note that although fluxes are computed twice in the case of CCL,this does not automatically translate into doubling of the computation, sincevalues already interpolated to the cell quadrature points can be interpolatedto a face with a simple 1D interpolation.
* <a name="CellcentricloopsandMatrixFree"></a><h4>Cell-centric loops and MatrixFree</h4>
* 

* For cell-centric loop implementations, the function   MatrixFree::loop_cell_centric()  can be used, to which the user can pass a function that should be performed oneach cell.
* To derive an appropriate function, which can be passed in   MatrixFree::loop_cell_centric(),  one might, in principle, transform/merge the following three functions, which canbe passed to a   MatrixFree::loop():  
* @code
matrix_free.template loop<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed on cells
* 
    FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);
* 
        // some operations on the cell quadrature points
* 
        phi.integrate_scatter(cell_evaluation_flags, dst);
      }
  },
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed inner faces
* 
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data,  /*is_interior_face=*/ true);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_p(data,  /*is_interior_face=*/ false);
* 
    for (unsigned int face = range.first; face < range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, face_evaluation_flags);
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, face_evaluation_flags);
* 
        // some operations on the face quadrature points
* 
        phi_m.integrate_scatter(face_evaluation_flags, dst);
        phi_p.integrate_scatter(face_evaluation_flags, dst);
      }
  },
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    // operation performed boundary faces
* 
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data,  /*is_interior_face=*/ true);
* 
    for (unsigned int face = range.first; face < range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src, face_evaluation_flags);
* 
        // some operations on the face quadrature points
* 
        phi_m.integrate_scatter(face_evaluation_flags, dst);
      }
  },
  dst,
  src);
@endcode
* 
* in the following way:
* @code
matrix_free.template loop_cell_centric<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    FEEvaluation<dim, degree, degree + 1, 1, Number>     phi(data);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_m(data,  /*is_interior_face=*/ true);
    FEFaceEvaluation<dim, degree, degree + 1, 1, Number> phi_p(data,  /*is_interior_face=*/ false);
* 
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);
* 
        // some operations on the cell quadrature points
* 
        phi.integrate_scatter(cell_evaluation_flags, dst);
* 
        // loop over all faces of cell
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          {
            if (data.get_faces_by_cells_boundary_id(cell, face)[0] ==
                numbers::internal_face_boundary_id)
              {
                // internal face
                phi_m.reinit(cell, face);
                phi_m.gather_evaluate(src, face_evaluation_flags);
                phi_p.reinit(cell, face);
                phi_p.gather_evaluate(src, face_evaluation_flags);
* 
                // some operations on the face quadrature points
* 
                phi_m.integrate_scatter(face_evaluation_flags, dst);
              }
            else
              {
                // boundary face
                phi_m.reinit(cell, face);
                phi_m.gather_evaluate(src, face_evaluation_flags);
* 
                // some operations on the face quadrature points
* 
                phi_m.integrate_scatter(face_evaluation_flags, dst);
              }
          }
      }
  },
  dst,
  src);
@endcode
* 
* It should be noted that FEFaceEvaluation is initialized now with two numbers,the cell number and the local face number. The given example onlyhighlights how to transform face-centric loops into cell-centric loops andis by no means efficient, since data is read and written multiple timesfrom and to the global vector as well as computations are performedredundantly. Below, we will discuss advanced techniques that target these issues.
* To be able to use   MatrixFree::loop_cell_centric(),   following flags of   MatrixFree::AdditionalData  have to be enabled:
* @code
typename MatrixFree<dim, Number>::AdditionalData additional_data;
* 
// set flags as usual (not shown)
* 
additional_data.hold_all_faces_to_owned_cells       = true;
additional_data.mapping_update_flags_faces_by_cells =
  additional_data.mapping_update_flags_inner_faces |
  additional_data.mapping_update_flags_boundary_faces;
* 
data.reinit(mapping, dof_handler, constraint, quadrature, additional_data);
@endcode
* 
* In particular, these flags enable that the internal data structures are set upfor all faces of the cells.
* Currently, cell-centric loops in deal.II only work for uniformly refined meshesand if no constraints are applied (which is the standard case DG is normallyused).
* 

* <a name="ProvidinglambdastoMatrixFreeloops"></a><h3>Providing lambdas to MatrixFree loops</h3>
* 

* The examples given above have already used lambdas, which have been provided tomatrix-free loops. The following short examples present how to transform functions betweena version where a class and a pointer to one of its methods are used and avariant where lambdas are utilized.
* In the following code, a class and a pointer to one of its methods, which shouldbe interpreted as cell integral, are passed to   MatrixFree::loop():  
* @code
void
local_apply_cell(const MatrixFree<dim, Number> &              data,
                 VectorType &                                 dst,
                 const VectorType &                           src,
                 const std::pair<unsigned int, unsigned int> &range) const
{
  FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, cell_evaluation_flags);
* 
      // some operations on the quadrature points
* 
      phi.integrate_scatter(cell_evaluation_flags, dst);
    }
}
@endcode
* 
* @code
matrix_free.cell_loop(&Operator::local_apply_cell, this, dst, src);
@endcode
* 
* However, it is also possible to pass an anonymous function via a lambda functionwith the same result:
* @code
matrix_free.template cell_loop<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto range) {
    FEEvaluation<dim, degree, degree + 1, 1, Number> phi(data);
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, cell_evaluation_flags);
* 
        // some operations on the quadrature points
* 
        phi.integrate_scatter(cell_evaluation_flags, dst);
      }
  },
  dst,
  src);
@endcode
* 
* <a name="VectorizedArrayType"></a><h3>VectorizedArrayType</h3>
* 

* The class VectorizedArray<Number> is a key component to achieve the highnode-level performance of the matrix-free algorithms in deal.II.It is a wrapper class around a short vector of   $n$   entries of type Number andmaps arithmetic operations to appropriate single-instruction/multiple-data(SIMD) concepts by intrinsic functions. The length of the vector can bequeried by   VectorizedArray::size()   and its underlying number type by  VectorizedArray::value_type.  
* In the default case (  <code>VectorizedArray<Number></code>  ), the vector length isset at compile time of the library tomatch the highest value supported by the given processor architecture.However, also a second optional template argument can bespecified as   <code>VectorizedArray<Number, size></code>, where <code>size</code>   explicitlycontrols the  vector length within the capabilities of a particular instructionset. A full list of supported vector lengths is presented in the following table:
*   <table align="center" class="doxtable">
  <tr>
   <th>double</th>
   <th>float</th>
   <th>ISA</th>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 1></code></td>
   <td><code>VectorizedArray<float, 1></code></td>
   <td>(auto-vectorization)</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 2></code></td>
   <td><code>VectorizedArray<float, 4></code></td>
   <td>SSE2/AltiVec</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 4></code></td>
   <td><code>VectorizedArray<float, 8></code></td>
   <td>AVX/AVX2</td>
  </tr>
  <tr>
   <td><code>VectorizedArray<double, 8></code></td>
   <td><code>VectorizedArray<float, 16></code></td>
   <td>AVX-512</td>
  </tr>
</table>  
* This allows users to select the vector length/ISA and, as a consequence, thenumber of cells to be processed at once in matrix-free operator evaluations,possibly reducing the pressure on the caches, an severe issue for very highdegrees (and dimensions).
* A possible further reason to reduce the number of filled lanesis to simplify debugging: instead of having to look at, e.g., 8cells, one can concentrate on a single cell.
* The interface of VectorizedArray also enables the replacement by any type witha matching interface. Specifically, this prepares deal.II for the   <code>std::simd</code>  class that is planned to become part of the C++23 standard. The following tablecompares the deal.II-specific SIMD classes and the equivalent C++23 classes:
* 

*   <table align="center" class="doxtable">
  <tr>
   <th>VectorizedArray (deal.II)</th>
   <th>std::simd (C++23)</th>
  </tr>
  <tr>
   <td><code>VectorizedArray<Number></code></td>
   <td><code>std::experimental::native_simd<Number></code></td>
  </tr>
  <tr>
   <td><code>VectorizedArray<Number, size></code></td>
   <td><code>std::experimental::fixed_size_simd<Number, size></code></td>
  </tr>
</table>  
* 

*  <a name="CommProg"></a> <h1> The commented program</h1>
*   <a name="Parametersandutilityfunctions"></a>  <h3>Parameters and utility functions</h3>
* 

* 
*  The same includes as in   step-67  :
* 

* 
* @code
 #include <deal.II/base/conditional_ostream.h>
 #include <deal.II/base/function.h>
 #include <deal.II/base/logstream.h>
 #include <deal.II/base/time_stepping.h>
 #include <deal.II/base/timer.h>
 #include <deal.II/base/utilities.h>
 #include <deal.II/base/vectorization.h>
* 
 #include <deal.II/distributed/tria.h>
* 
 #include <deal.II/dofs/dof_handler.h>
* 
 #include <deal.II/fe/fe_dgq.h>
 #include <deal.II/fe/fe_system.h>
* 
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/tria.h>
 #include <deal.II/grid/tria_accessor.h>
 #include <deal.II/grid/tria_iterator.h>
* 
 #include <deal.II/lac/affine_constraints.h>
 #include <deal.II/lac/la_parallel_vector.h>
* 
 #include <deal.II/matrix_free/fe_evaluation.h>
 #include <deal.II/matrix_free/matrix_free.h>
 #include <deal.II/matrix_free/operators.h>
* 
 #include <deal.II/numerics/data_out.h>
* 
 #include <fstream>
 #include <iomanip>
 #include <iostream>
* 
 @endcode
* 
*  A new include for categorizing of cells according to their boundary IDs:
* 

* 
* @code
 #include <deal.II/matrix_free/tools.h>
* 
 
* 
 namespace Euler_DG
 {
   using namespace dealii;
* 
 @endcode
* 
*  The same input parameters as in   step-67  :
* 

* 
* @code
   constexpr unsigned int testcase             = 1;
   constexpr unsigned int dimension            = 2;
   constexpr unsigned int n_global_refinements = 2;
   constexpr unsigned int fe_degree            = 5;
   constexpr unsigned int n_q_points_1d        = fe_degree + 2;
* 
 @endcode
* 
*  This parameter specifies the size of the shared-memory group. Currently, only the values 1 and   numbers::invalid_unsigned_int   is possible, leading to the options that the memory features can be turned off or all processes having access to the same shared-memory domain are grouped together.
* 

* 
* @code
   constexpr unsigned int group_size = numbers::invalid_unsigned_int;
* 
   using Number = double;
* 
 @endcode
* 
*  Here, the type of the data structure is chosen for vectorization. In the default case, VectorizedArray<Number> is used, i.e., the highest instruction-set-architecture extension available on the given hardware with the maximum number of vector lanes is used. However, one might reduce the number of filled lanes, e.g., by writing   <code>using VectorizedArrayType = VectorizedArray<Number, 4></code>   to only process 4 cells.
* 

* 
* @code
   using VectorizedArrayType = VectorizedArray<Number>;
* 
 @endcode
* 
*  The following parameters have not changed:
* 

* 
* @code
   constexpr double gamma       = 1.4;
   constexpr double final_time  = testcase == 0 ? 10 : 2.0;
   constexpr double output_tick = testcase == 0 ? 1 : 0.05;
* 
   const double courant_number = 0.15 / std::pow(fe_degree, 1.5);
* 
 @endcode
* 
*  Specify max number of time steps useful for performance studies.
* 

* 
* @code
   constexpr unsigned int max_time_steps = numbers::invalid_unsigned_int;
* 
 @endcode
* 
*  Runge-Kutta-related functions copied from   step-67   and slightly modified with the purpose to minimize global vector access:
* 

* 
* @code
   enum LowStorageRungeKuttaScheme
   {
     stage_3_order_3,
     stage_5_order_4,
     stage_7_order_4,
     stage_9_order_5,
   };
   constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4;
* 
 
* 
   class LowStorageRungeKuttaIntegrator
   {
   public:
     LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme)
     {
       TimeStepping::runge_kutta_method lsrk;
       switch (scheme)
         {
           case stage_3_order_3:
             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3;
             break;
           case stage_5_order_4:
             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4;
             break;
           case stage_7_order_4:
             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4;
             break;
           case stage_9_order_5:
             lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5;
             break;
* 
           default:
             AssertThrow(false, ExcNotImplemented());
         }
       TimeStepping::LowStorageRungeKutta<
         LinearAlgebra::distributed::Vector<Number>>
                           rk_integrator(lsrk);
       std::vector<double> ci; // not used
       rk_integrator.get_coefficients(ai, bi, ci);
     }
* 
     unsigned int n_stages() const
     {
       return bi.size();
     }
* 
     template <typename VectorType, typename Operator>
     void perform_time_step(const Operator &pde_operator,
                            const double    current_time,
                            const double    time_step,
                            VectorType &    solution,
                            VectorType &    vec_ri,
                            VectorType &    vec_ki) const
     {
       AssertDimension(ai.size() + 1, bi.size());
* 
       vec_ki.swap(solution);
* 
       double sum_previous_bi = 0;
       for (unsigned int stage = 0; stage < bi.size(); ++stage)
         {
           const double c_i = stage == 0 ? 0 : sum_previous_bi + ai[stage
* 
- 1];
* 
           pde_operator.perform_stage(stage,
                                      current_time + c_i time_step,
                                      bi[stage] time_step,
                                      (stage == bi.size()
* 
- 1 ?
                                         0 :
                                         ai[stage] time_step),
                                      (stage % 2 == 0 ? vec_ki : vec_ri),
                                      (stage % 2 == 0 ? vec_ri : vec_ki),
                                      solution);
* 
           if (stage > 0)
             sum_previous_bi += bi[stage
* 
- 1];
         }
     }
* 
   private:
     std::vector<double> bi;
     std::vector<double> ai;
   };
* 
 
 @endcode
* 
*  Euler-specific utility functions from   step-67  :
* 

* 
* @code
   enum EulerNumericalFlux
   {
     lax_friedrichs_modified,
     harten_lax_vanleer,
   };
   constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified;
* 
 
* 
   template <int dim>
   class ExactSolution : public Function<dim>
   {
   public:
     ExactSolution(const double time)
       : Function<dim>(dim + 2, time)
     {}
* 
     virtual double value(const Point<dim> & p,
                          const unsigned int component = 0) const override;
   };
* 
 
* 
   template <int dim>
   double ExactSolution<dim>::value(const Point<dim> & x,
                                    const unsigned int component) const
   {
     const double t = this->get_time();
* 
     switch (testcase)
       {
         case 0:
           {
             Assert(dim == 2, ExcNotImplemented());
             const double beta = 5;
* 
             Point<dim> x0;
             x0[0] = 5.;
             const double radius_sqr =
               (x
* 
- x0).norm_square()
* 
- 2. (x[0]
* 
- x0[0]) t + t t;
             const double factor =
               beta / (numbers::PI 2) std::exp(1.
* 
- radius_sqr);
             const double density_log = std::log2(
               std::abs(1.
* 
- (gamma
* 
- 1.) / gamma 0.25 factor factor));
             const double density = std::exp2(density_log (1. / (gamma
* 
- 1.)));
             const double u       = 1.
* 
- factor (x[1]
* 
- x0[1]);
             const double v       = factor (x[0]
* 
- t
* 
- x0[0]);
* 
             if (component == 0)
               return density;
             else if (component == 1)
               return density u;
             else if (component == 2)
               return density v;
             else
               {
                 const double pressure =
                   std::exp2(density_log (gamma / (gamma
* 
- 1.)));
                 return pressure / (gamma
* 
- 1.) +
                        0.5 (density u u + density v v);
               }
           }
* 
         case 1:
           {
             if (component == 0)
               return 1.;
             else if (component == 1)
               return 0.4;
             else if (component == dim + 1)
               return 3.097857142857143;
             else
               return 0.;
           }
* 
         default:
           Assert(false, ExcNotImplemented());
           return 0.;
       }
   }
* 
 
* 
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE 
     Tensor<1, dim, Number>
     euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
   {
     const Number inverse_density = Number(1.) / conserved_variables[0];
* 
     Tensor<1, dim, Number> velocity;
     for (unsigned int d = 0; d < dim; ++d)
       velocity[d] = conserved_variables[1 + d] inverse_density;
* 
     return velocity;
   }
* 
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE 
     Number
     euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
   {
     const Tensor<1, dim, Number> velocity =
       euler_velocity<dim>(conserved_variables);
* 
     Number rho_u_dot_u = conserved_variables[1] velocity[0];
     for (unsigned int d = 1; d < dim; ++d)
       rho_u_dot_u += conserved_variables[1 + d] velocity[d];
* 
     return (gamma
* 
- 1.) (conserved_variables[dim + 1]
* 
- 0.5 rho_u_dot_u);
   }
* 
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE 
     Tensor<1, dim + 2, Tensor<1, dim, Number>>
     euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
   {
     const Tensor<1, dim, Number> velocity =
       euler_velocity<dim>(conserved_variables);
     const Number pressure = euler_pressure<dim>(conserved_variables);
* 
     Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
     for (unsigned int d = 0; d < dim; ++d)
       {
         flux[0][d] = conserved_variables[1 + d];
         for (unsigned int e = 0; e < dim; ++e)
           flux[e + 1][d] = conserved_variables[e + 1] velocity[d];
         flux[d + 1][d] += pressure;
         flux[dim + 1][d] =
           velocity[d] (conserved_variables[dim + 1] + pressure);
       }
* 
     return flux;
   }
* 
   template <int n_components, int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE 
     Tensor<1, n_components, Number>
     operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
               const Tensor<1, dim, Number> &                         vector)
   {
     Tensor<1, n_components, Number> result;
     for (unsigned int d = 0; d < n_components; ++d)
       result[d] = matrix[d] vector;
     return result;
   }
* 
   template <int dim, typename Number>
   inline DEAL_II_ALWAYS_INLINE 
     Tensor<1, dim + 2, Number>
     euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
                          const Tensor<1, dim + 2, Number> &u_p,
                          const Tensor<1, dim, Number> &    normal)
   {
     const auto velocity_m = euler_velocity<dim>(u_m);
     const auto velocity_p = euler_velocity<dim>(u_p);
* 
     const auto pressure_m = euler_pressure<dim>(u_m);
     const auto pressure_p = euler_pressure<dim>(u_p);
* 
     const auto flux_m = euler_flux<dim>(u_m);
     const auto flux_p = euler_flux<dim>(u_p);
* 
     switch (numerical_flux_type)
       {
         case lax_friedrichs_modified:
           {
             const auto lambda =
               0.5 std::sqrt(std::max(velocity_p.norm_square() +
                                          gamma pressure_p (1. / u_p[0]),
                                        velocity_m.norm_square() +
                                          gamma pressure_m (1. / u_m[0])));
* 
             return 0.5 (flux_m normal + flux_p normal) +
                    0.5 lambda (u_m
* 
- u_p);
           }
* 
         case harten_lax_vanleer:
           {
             const auto avg_velocity_normal =
               0.5 ((velocity_m + velocity_p) normal);
             const auto   avg_c = std::sqrt(std::abs(
               0.5 gamma
               (pressure_p (1. / u_p[0]) + pressure_m (1. / u_m[0]))));
             const Number s_pos =
               std::max(Number(), avg_velocity_normal + avg_c);
             const Number s_neg =
               std::min(Number(), avg_velocity_normal
* 
- avg_c);
             const Number inverse_s = Number(1.) / (s_pos
* 
- s_neg);
* 
             return inverse_s
                    ((s_pos (flux_m normal)
* 
- s_neg (flux_p normal))
* 
-
                     s_pos s_neg (u_m
* 
- u_p));
           }
* 
         default:
           {
             Assert(false, ExcNotImplemented());
             return {};
           }
       }
   }
* 
 
* 
 @endcode
* 
*  General-purpose utility functions from   step-67  :
* 

* 
* @code
   template <int dim, typename VectorizedArrayType>
   VectorizedArrayType
   evaluate_function(const Function<dim> &                  function,
                     const Point<dim, VectorizedArrayType> &p_vectorized,
                     const unsigned int                     component)
   {
     VectorizedArrayType result;
     for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
       {
         Point<dim> p;
         for (unsigned int d = 0; d < dim; ++d)
           p[d] = p_vectorized[d][v];
         result[v] = function.value(p, component);
       }
     return result;
   }
* 
 
   template <int dim, typename VectorizedArrayType, int n_components = dim + 2>
   Tensor<1, n_components, VectorizedArrayType>
   evaluate_function(const Function<dim> &                  function,
                     const Point<dim, VectorizedArrayType> &p_vectorized)
   {
     AssertDimension(function.n_components, n_components);
     Tensor<1, n_components, VectorizedArrayType> result;
     for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
       {
         Point<dim> p;
         for (unsigned int d = 0; d < dim; ++d)
           p[d] = p_vectorized[d][v];
         for (unsigned int d = 0; d < n_components; ++d)
           result[d][v] = function.value(p, d);
       }
     return result;
   }
* 
 
 @endcode
* 
*   <a name="EuleroperatorusingacellcentricloopandMPI30sharedmemory"></a>  <h3>Euler operator using a cell-centric loop and MPI-3.0 shared memory</h3>
* 

* 
*  Euler operator from   step-67   with some changes as detailed below:
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   class EulerOperator
   {
   public:
     static constexpr unsigned int n_quadrature_points_1d = n_points_1d;
* 
     EulerOperator(TimerOutput &timer_output);
* 
     ~EulerOperator();
* 
     void reinit(const Mapping<dim> &   mapping,
                 const DoFHandler<dim> &dof_handler);
* 
     void set_inflow_boundary(const types::boundary_id       boundary_id,
                              std::unique_ptr<Function<dim>> inflow_function);
* 
     void set_subsonic_outflow_boundary(
       const types::boundary_id       boundary_id,
       std::unique_ptr<Function<dim>> outflow_energy);
* 
     void set_wall_boundary(const types::boundary_id boundary_id);
* 
     void set_body_force(std::unique_ptr<Function<dim>> body_force);
* 
     void
     perform_stage(const unsigned int                                stage,
                   const Number                                      cur_time,
                   const Number                                      bi,
                   const Number                                      ai,
                   const LinearAlgebra::distributed::Vector<Number> &current_ri,
                   LinearAlgebra::distributed::Vector<Number> &      vec_ki,
                   LinearAlgebra::distributed::Vector<Number> &solution) const;
* 
     void project(const Function<dim> &                       function,
                  LinearAlgebra::distributed::Vector<Number> &solution) const;
* 
     std::array<double, 3> compute_errors(
       const Function<dim> &                             function,
       const LinearAlgebra::distributed::Vector<Number> &solution) const;
* 
     double compute_cell_transport_speed(
       const LinearAlgebra::distributed::Vector<Number> &solution) const;
* 
     void
     initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const;
* 
   private:
 @endcode
* 
*  Instance of SubCommunicatorWrapper containing the sub-communicator, which we need to pass to   MatrixFree::reinit()   to be able to exploit MPI-3.0 shared-memory capabilities:
* 

* 
* @code
     MPI_Comm subcommunicator;
* 
     MatrixFree<dim, Number, VectorizedArrayType> data;
* 
     TimerOutput &timer;
* 
     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
       inflow_boundaries;
     std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
                                    subsonic_outflow_boundaries;
     std::set<types::boundary_id>   wall_boundaries;
     std::unique_ptr<Function<dim>> body_force;
   };
* 
 
* 
 @endcode
* 
*  New constructor, which creates a sub-communicator. The user can specify the size of the sub-communicator via the global parameter group_size. If the size is set to
* 
*  - , all MPI processes of a shared-memory domain are combined to a group. The specified size is decisive for the benefit of the shared-memory capabilities of MatrixFree and, therefore, setting the   <code>size</code> to <code>-1</code>   is a reasonable choice. By setting, the size to   <code>1</code>   users explicitly disable the MPI-3.0 shared-memory features of MatrixFree and rely completely on MPI-2.0 features, like   <code>MPI_Isend</code>   and   <code>MPI_Irecv</code>  .
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer)
     : timer(timer)
   {
 #if DEAL_II_MPI_VERSION_GTE(3, 0)
     if (group_size == 1)
       {
         this->subcommunicator = MPI_COMM_SELF;
       }
     else if (group_size == numbers::invalid_unsigned_int)
       {
         const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
* 
         MPI_Comm_split_type(MPI_COMM_WORLD,
                             MPI_COMM_TYPE_SHARED,
                             rank,
                             MPI_INFO_NULL,
                             &subcommunicator);
       }
     else
       {
         Assert(false, ExcNotImplemented());
       }
 #else
     (void)subcommunicator;
     (void)group_size;
     this->subcommunicator = MPI_COMM_SELF;
 #endif
   }
* 
 
 @endcode
* 
*  New destructor responsible for freeing of the sub-communicator.
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   EulerOperator<dim, degree, n_points_1d>::~EulerOperator()
   {
 #ifdef DEAL_II_WITH_MPI
     if (this->subcommunicator != MPI_COMM_SELF)
       MPI_Comm_free(&subcommunicator);
 #endif
   }
* 
 
 @endcode
* 
*  Modified reinit() function to setup the internal data structures in MatrixFree in a way that it is usable by the cell-centric loops and the MPI-3.0 shared-memory capabilities are used:
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::reinit(
     const Mapping<dim> &   mapping,
     const DoFHandler<dim> &dof_handler)
   {
     const std::vector<const DoFHandler<dim>> dof_handlers = {&dof_handler};
     const AffineConstraints<double>            dummy;
     const std::vector<const AffineConstraints<double>> constraints = {&dummy};
     const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
                                                     QGauss<1>(fe_degree + 1)};
* 
     typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
       additional_data;
     additional_data.mapping_update_flags =
       (update_gradients | update_JxW_values | update_quadrature_points |
        update_values);
     additional_data.mapping_update_flags_inner_faces =
       (update_JxW_values | update_quadrature_points | update_normal_vectors |
        update_values);
     additional_data.mapping_update_flags_boundary_faces =
       (update_JxW_values | update_quadrature_points | update_normal_vectors |
        update_values);
     additional_data.tasks_parallel_scheme =
       MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none;
* 
 @endcode
* 
*  Categorize cells so that all lanes have the same boundary IDs for each face. This is strictly not necessary, however, allows to write simpler code in   EulerOperator::perform_stage()   without masking, since it is guaranteed that all cells grouped together (in a VectorizedArray) have to perform exactly the same operation also on the faces.
* 

* 
* @code
     MatrixFreeTools::categorize_by_boundary_ids(dof_handler.get_triangulation(),
                                                 additional_data);
* 
 @endcode
* 
*  Enable MPI-3.0 shared-memory capabilities within MatrixFree by providing the sub-communicator:
* 

* 
* @code
     additional_data.communicator_sm = subcommunicator;
* 
     data.reinit(
       mapping, dof_handlers, constraints, quadratures, additional_data);
   }
* 
 
 @endcode
* 
*  The following function does an entire stage of a Runge--Kutta update and is
* 

* 
* 
*  - alongside the slightly modified setup
* 
*  - the heart of this tutorial compared to   step-67  .   
*   In contrast to   step-67  , we are not executing the advection step (using   MatrixFree::loop())   and the inverse mass-matrix step (using   MatrixFree::cell_loop())   in sequence, but evaluate everything in one go inside of   MatrixFree::loop_cell_centric().   This function expects a single function that is executed on each locally-owned (macro) cell as parameter so that we need to loop over all faces of that cell and perform needed integration steps on our own.   
*   The following function contains to a large extent copies of the following functions from   step-67   so that comments related the evaluation of the weak form are skipped here:
* 

* 
* 
*  -   <code>EulerDG::EulerOperator::local_apply_cell</code>  
* 

* 
* 
*  -   <code>EulerDG::EulerOperator::local_apply_face</code>  
* 

* 
* 
*  -   <code>EulerDG::EulerOperator::local_apply_boundary_face</code>  
* 

* 
* 
*  -   <code>EulerDG::EulerOperator::local_apply_inverse_mass_matrix</code>  
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::perform_stage(
     const unsigned int                                stage,
     const Number                                      current_time,
     const Number                                      bi,
     const Number                                      ai,
     const LinearAlgebra::distributed::Vector<Number> &current_ri,
     LinearAlgebra::distributed::Vector<Number> &      vec_ki,
     LinearAlgebra::distributed::Vector<Number> &      solution) const
   {
     for (auto &i : inflow_boundaries)
       i.second->set_time(current_time);
     for (auto &i : subsonic_outflow_boundaries)
       i.second->set_time(current_time);
* 
 @endcode
* 
*  Run a cell-centric loop by calling   MatrixFree::loop_cell_centric()   and providing a lambda containing the effects of the cell, face and boundary-face integrals:
* 

* 
* @code
     data.template loop_cell_centric<LinearAlgebra::distributed::Vector<Number>,
                                     LinearAlgebra::distributed::Vector<Number>>(
       [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
         using FECellIntegral = FEEvaluation<dim,
                                             degree,
                                             n_points_1d,
                                             dim + 2,
                                             Number,
                                             VectorizedArrayType>;
         using FEFaceIntegral = FEFaceEvaluation<dim,
                                                 degree,
                                                 n_points_1d,
                                                 dim + 2,
                                                 Number,
                                                 VectorizedArrayType>;
* 
         FECellIntegral phi(data);
         FECellIntegral phi_temp(data);
         FEFaceIntegral phi_m(data, true);
         FEFaceIntegral phi_p(data, false);
* 
         Tensor<1, dim, VectorizedArrayType>     constant_body_force;
         const Functions::ConstantFunction<dim>constant_function =
           dynamic_cast<Functions::ConstantFunction<dim>>(body_force.get());
* 
         if (constant_function)
           constant_body_force =
             evaluate_function<dim, VectorizedArrayType, dim>(
              constant_function, Point<dim, VectorizedArrayType>());
* 
         const dealii::internal::EvaluatorTensorProduct<
           dealii::internal::EvaluatorVariant::evaluate_evenodd,
           dim,
           n_points_1d,
           n_points_1d,
           VectorizedArrayType>
           eval(AlignedVector<VectorizedArrayType>(),
                data.get_shape_info().data[0].shape_gradients_collocation_eo,
                AlignedVector<VectorizedArrayType>());
* 
         AlignedVector<VectorizedArrayType> buffer(phi.static_n_q_points
                                                   phi.n_components);
* 
 @endcode
* 
*  Loop over all cell batches:
* 

* 
* @code
         for (unsigned int cell = cell_range.first; cell < cell_range.second;
              ++cell)
           {
             phi.reinit(cell);
* 
             if (ai != Number())
               phi_temp.reinit(cell);
* 
 @endcode
* 
*  Read values from global vector and compute the values at the quadrature points:
* 

* 
* @code
             if (ai != Number() && stage == 0)
               {
                 phi.read_dof_values(src);
* 
                 for (unsigned int i = 0;
                      i < phi.static_dofs_per_component (dim + 2);
                      ++i)
                   phi_temp.begin_dof_values()[i] = phi.begin_dof_values()[i];
* 
                 phi.evaluate(EvaluationFlags::values);
               }
             else
               {
                 phi.gather_evaluate(src, EvaluationFlags::values);
               }
* 
 @endcode
* 
*  Buffer the computed values at the quadrature points, since these are overridden by   FEEvaluation::submit_value()   in the next step, however, are needed later on for the face integrals:
* 

* 
* @code
             for (unsigned int i = 0; i < phi.static_n_q_points (dim + 2); ++i)
               buffer[i] = phi.begin_values()[i];
* 
 @endcode
* 
*  Apply the cell integral at the cell quadrature points. See also the function   <code>EulerOperator::local_apply_cell()</code>   from   step-67  :
* 

* 
* @code
             for (unsigned int q = 0; q < phi.n_q_points; ++q)
               {
                 const auto w_q = phi.get_value(q);
                 phi.submit_gradient(euler_flux<dim>(w_q), q);
                 if (body_force.get() != nullptr)
                   {
                     const Tensor<1, dim, VectorizedArrayType> force =
                       constant_function ?
                         constant_body_force :
                         evaluate_function<dim, VectorizedArrayType, dim>(
                          body_force, phi.quadrature_point(q));
* 
                     Tensor<1, dim + 2, VectorizedArrayType> forcing;
                     for (unsigned int d = 0; d < dim; ++d)
                       forcing[d + 1] = w_q[0] force[d];
                     for (unsigned int d = 0; d < dim; ++d)
                       forcing[dim + 1] += force[d] w_q[d + 1];
* 
                     phi.submit_value(forcing, q);
                   }
               }
* 
 @endcode
* 
*  Test with the gradient of the test functions in the quadrature points. We skip the interpolation back to the support points of the element, since we first collect all contributions in the cell quadrature points and only perform the interpolation back as the final step.
* 

* 
* @code
             {
               autovalues_ptr   = phi.begin_values();
               autogradient_ptr = phi.begin_gradients();
* 
               for (unsigned int c = 0; c < dim + 2; ++c)
                 {
                   if (dim >= 1 && body_force.get() == nullptr)
                     eval.template gradients<0, false, false>(
                       gradient_ptr + phi.static_n_q_points 0, values_ptr);
                   else if (dim >= 1)
                     eval.template gradients<0, false, true>(
                       gradient_ptr + phi.static_n_q_points 0, values_ptr);
                   if (dim >= 2)
                     eval.template gradients<1, false, true>(
                       gradient_ptr + phi.static_n_q_points 1, values_ptr);
                   if (dim >= 3)
                     eval.template gradients<2, false, true>(
                       gradient_ptr + phi.static_n_q_points 2, values_ptr);
* 
                   values_ptr += phi.static_n_q_points;
                   gradient_ptr += phi.static_n_q_points dim;
                 }
             }
* 
 @endcode
* 
*  Loop over all faces of the current cell:
* 

* 
* @code
             for (unsigned int face = 0;
                  face < GeometryInfo<dim>::faces_per_cell;
                  ++face)
               {
 @endcode
* 
*  Determine the boundary ID of the current face. Since we have set up MatrixFree in a way that all filled lanes have guaranteed the same boundary ID, we can select the boundary ID of the first lane.
* 

* 
* @code
                 const auto boundary_ids =
                   data.get_faces_by_cells_boundary_id(cell, face);
* 
                 Assert(std::equal(boundary_ids.begin(),
                                   boundary_ids.begin() +
                                     data.n_active_entries_per_cell_batch(cell),
                                   boundary_ids.begin()),
                        ExcMessage("Boundary IDs of lanes differ."));
* 
                 const auto boundary_id = boundary_ids[0];
* 
                 phi_m.reinit(cell, face);
* 
 @endcode
* 
*  Interpolate the values from the cell quadrature points to the quadrature points of the current face via a simple 1D interpolation:
* 

* 
* @code
                 internal::FEFaceNormalEvaluationImpl<dim,
                                                      n_points_1d
* 
- 1,
                                                      VectorizedArrayType>::
                   template interpolate_quadrature<true, false>(
                     dim + 2,
                     data.get_shape_info(),
                     buffer.data(),
                     phi_m.begin_values(),
                     false,
                     face);
* 
 @endcode
* 
*  Check if the face is an internal or a boundary face and select a different code path based on this information:
* 

* 
* @code
                 if (boundary_id == numbers::internal_face_boundary_id)
                   {
 @endcode
* 
*  Process and internal face. The following lines of code are a copy of the function   <code>EulerDG::EulerOperator::local_apply_face</code>   from   step-67  :
* 

* 
* @code
                     phi_p.reinit(cell, face);
                     phi_p.gather_evaluate(src, EvaluationFlags::values);
* 
                     for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                       {
                         const auto numerical_flux =
                           euler_numerical_flux<dim>(phi_m.get_value(q),
                                                     phi_p.get_value(q),
                                                     phi_m.get_normal_vector(q));
                         phi_m.submit_value(-numerical_flux, q);
                       }
                   }
                 else
                   {
 @endcode
* 
*  Process a boundary face. These following lines of code are a copy of the function   <code>EulerDG::EulerOperator::local_apply_boundary_face</code>   from   step-67  :
* 

* 
* @code
                     for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                       {
                         const auto w_m    = phi_m.get_value(q);
                         const auto normal = phi_m.get_normal_vector(q);
* 
                         auto rho_u_dot_n = w_m[1] normal[0];
                         for (unsigned int d = 1; d < dim; ++d)
                           rho_u_dot_n += w_m[1 + d] normal[d];
* 
                         bool at_outflow = false;
* 
                         Tensor<1, dim + 2, VectorizedArrayType> w_p;
* 
                         if (wall_boundaries.find(boundary_id) !=
                             wall_boundaries.end())
                           {
                             w_p[0] = w_m[0];
                             for (unsigned int d = 0; d < dim; ++d)
                               w_p[d + 1] =
                                 w_m[d + 1]
* 
- 2. rho_u_dot_n normal[d];
                             w_p[dim + 1] = w_m[dim + 1];
                           }
                         else if (inflow_boundaries.find(boundary_id) !=
                                  inflow_boundaries.end())
                           w_p = evaluate_function(
                            inflow_boundaries.find(boundary_id)->second,
                             phi_m.quadrature_point(q));
                         else if (subsonic_outflow_boundaries.find(
                                    boundary_id) !=
                                  subsonic_outflow_boundaries.end())
                           {
                             w_p = w_m;
                             w_p[dim + 1] =
                               evaluate_function(*subsonic_outflow_boundaries
                                                    .find(boundary_id)
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 


* 
* 
->second,
                                                 phi_m.quadrature_point(q),
                                                 dim + 1);
                             at_outflow = true;
                           }
                         else
                           AssertThrow(false,
                                       ExcMessage(
                                         "Unknown boundary id, did "
                                         "you set a boundary condition for "
                                         "this part of the domain boundary?"));
* 
                         auto flux = euler_numerical_flux<dim>(w_m, w_p, normal);
* 
                         if (at_outflow)
                           for (unsigned int v = 0;
                                v < VectorizedArrayType::size();
                                ++v)
                             {
                               if (rho_u_dot_n[v] <
* 
-1e-12)
                                 for (unsigned int d = 0; d < dim; ++d)
                                   flux[d + 1][v] = 0.;
                             }
* 
                         phi_m.submit_value(-flux, q);
                       }
                   }
* 
 @endcode
* 
*  Evaluate local integrals related to cell by quadrature and add into cell contribution via a simple 1D interpolation:
* 

* 
* @code
                 internal::FEFaceNormalEvaluationImpl<dim,
                                                      n_points_1d
* 
- 1,
                                                      VectorizedArrayType>::
                   template interpolate_quadrature<false, true>(
                     dim + 2,
                     data.get_shape_info(),
                     phi_m.begin_values(),
                     phi.begin_values(),
                     false,
                     face);
               }
* 
 @endcode
* 
*  Apply inverse mass matrix in the cell quadrature points. See also the function   <code>EulerDG::EulerOperator::local_apply_inverse_mass_matrix()</code>   from   step-67  :
* 

* 
* @code
             for (unsigned int q = 0; q < phi.static_n_q_points; ++q)
               {
                 const auto factor = VectorizedArrayType(1.0) / phi.JxW(q);
                 for (unsigned int c = 0; c < dim + 2; ++c)
                   phi.begin_values()[c phi.static_n_q_points + q] =
                     phi.begin_values()[c phi.static_n_q_points + q] factor;
               }
* 
 @endcode
* 
*  Transform values from collocation space to the original Gauss-Lobatto space:
* 

* 
* @code
             internal::FEEvaluationImplBasisChange<
               dealii::internal::EvaluatorVariant::evaluate_evenodd,
               internal::EvaluatorQuantity::hessian,
               dim,
               degree + 1,
               n_points_1d,
               VectorizedArrayType,
               VectorizedArrayType>::do_backward(dim + 2,
                                                 data.get_shape_info()
                                                   .data[0]
                                                   .inverse_shape_values_eo,
                                                 false,
                                                 phi.begin_values(),
                                                 phi.begin_dof_values());
* 
 @endcode
* 
*  Perform Runge-Kutta update and write results back to global vectors:
* 

* 
* @code
             if (ai == Number())
               {
                 for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q)
                   phi.begin_dof_values()[q] = bi phi.begin_dof_values()[q];
                 phi.distribute_local_to_global(solution);
               }
             else
               {
                 if (stage != 0)
                   phi_temp.read_dof_values(solution);
* 
                 for (unsigned int q = 0; q < phi.static_dofs_per_cell; ++q)
                   {
                     const auto K_i = phi.begin_dof_values()[q];
* 
                     phi.begin_dof_values()[q] =
                       phi_temp.begin_dof_values()[q] + (ai K_i);
* 
                     phi_temp.begin_dof_values()[q] += bi K_i;
                   }
                 phi.set_dof_values(dst);
                 phi_temp.set_dof_values(solution);
               }
           }
       },
       vec_ki,
       current_ri,
       true,
       MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values);
   }
* 
 
* 
 @endcode
* 
*  From here, the code of   step-67   has not changed.
* 

* 
* @code
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::initialize_vector(
     LinearAlgebra::distributed::Vector<Number> &vector) const
   {
     data.initialize_dof_vector(vector);
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary(
     const types::boundary_id       boundary_id,
     std::unique_ptr<Function<dim>> inflow_function)
   {
     AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
                     subsonic_outflow_boundaries.end() &&
                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                 ExcMessage("You already set the boundary with id " +
                            std::to_string(static_cast<int>(boundary_id)) +
                            " to another type of boundary before now setting " +
                            "it as inflow"));
     AssertThrow(inflow_function->n_components == dim + 2,
                 ExcMessage("Expected function with dim+2 components"));
* 
     inflow_boundaries[boundary_id] = std::move(inflow_function);
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary(
     const types::boundary_id       boundary_id,
     std::unique_ptr<Function<dim>> outflow_function)
   {
     AssertThrow(inflow_boundaries.find(boundary_id) ==
                     inflow_boundaries.end() &&
                   wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                 ExcMessage("You already set the boundary with id " +
                            std::to_string(static_cast<int>(boundary_id)) +
                            " to another type of boundary before now setting " +
                            "it as subsonic outflow"));
     AssertThrow(outflow_function->n_components == dim + 2,
                 ExcMessage("Expected function with dim+2 components"));
* 
     subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary(
     const types::boundary_id boundary_id)
   {
     AssertThrow(inflow_boundaries.find(boundary_id) ==
                     inflow_boundaries.end() &&
                   subsonic_outflow_boundaries.find(boundary_id) ==
                     subsonic_outflow_boundaries.end(),
                 ExcMessage("You already set the boundary with id " +
                            std::to_string(static_cast<int>(boundary_id)) +
                            " to another type of boundary before now setting " +
                            "it as wall boundary"));
* 
     wall_boundaries.insert(boundary_id);
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::set_body_force(
     std::unique_ptr<Function<dim>> body_force)
   {
     AssertDimension(body_force->n_components, dim);
* 
     this->body_force = std::move(body_force);
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   void EulerOperator<dim, degree, n_points_1d>::project(
     const Function<dim> &                       function,
     LinearAlgebra::distributed::Vector<Number> &solution) const
   {
     FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType>
       phi(data, 0, 1);
     MatrixFreeOperators::CellwiseInverseMassMatrix<dim,
                                                    degree,
                                                    dim + 2,
                                                    Number,
                                                    VectorizedArrayType>
       inverse(phi);
     solution.zero_out_ghost_values();
     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
       {
         phi.reinit(cell);
         for (unsigned int q = 0; q < phi.n_q_points; ++q)
           phi.submit_dof_value(evaluate_function(function,
                                                  phi.quadrature_point(q)),
                                q);
         inverse.transform_from_q_points_to_basis(dim + 2,
                                                  phi.begin_dof_values(),
                                                  phi.begin_dof_values());
         phi.set_dof_values(solution);
       }
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors(
     const Function<dim> &                             function,
     const LinearAlgebra::distributed::Vector<Number> &solution) const
   {
     TimerOutput::Scope t(timer, "compute errors");
     double             errors_squared[3] = {};
     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number, VectorizedArrayType>
       phi(data, 0, 0);
* 
     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
       {
         phi.reinit(cell);
         phi.gather_evaluate(solution, EvaluationFlags::values);
         VectorizedArrayType local_errors_squared[3] = {};
         for (unsigned int q = 0; q < phi.n_q_points; ++q)
           {
             const auto error =
               evaluate_function(function, phi.quadrature_point(q))
* 
-
               phi.get_value(q);
             const auto JxW = phi.JxW(q);
* 
             local_errors_squared[0] += error[0] error[0] JxW;
             for (unsigned int d = 0; d < dim; ++d)
               local_errors_squared[1] += (error[d + 1] error[d + 1]) JxW;
             local_errors_squared[2] += (error[dim + 1] error[dim + 1]) JxW;
           }
         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
              ++v)
           for (unsigned int d = 0; d < 3; ++d)
             errors_squared[d] += local_errors_squared[d][v];
       }
* 
     Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);
* 
     std::array<double, 3> errors;
     for (unsigned int d = 0; d < 3; ++d)
       errors[d] = std::sqrt(errors_squared[d]);
* 
     return errors;
   }
* 
 
* 
   template <int dim, int degree, int n_points_1d>
   double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed(
     const LinearAlgebra::distributed::Vector<Number> &solution) const
   {
     TimerOutput::Scope t(timer, "compute transport speed");
     Number             max_transport = 0;
     FEEvaluation<dim, degree, degree + 1, dim + 2, Number, VectorizedArrayType>
       phi(data, 0, 1);
* 
     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
       {
         phi.reinit(cell);
         phi.gather_evaluate(solution, EvaluationFlags::values);
         VectorizedArrayType local_max = 0.;
         for (unsigned int q = 0; q < phi.n_q_points; ++q)
           {
             const auto solution = phi.get_value(q);
             const auto velocity = euler_velocity<dim>(solution);
             const auto pressure = euler_pressure<dim>(solution);
* 
             const auto          inverse_jacobian = phi.inverse_jacobian(q);
             const auto          convective_speed = inverse_jacobian velocity;
             VectorizedArrayType convective_limit = 0.;
             for (unsigned int d = 0; d < dim; ++d)
               convective_limit =
                 std::max(convective_limit, std::abs(convective_speed[d]));
* 
             const auto speed_of_sound =
               std::sqrt(gamma pressure (1. / solution[0]));
* 
             Tensor<1, dim, VectorizedArrayType> eigenvector;
             for (unsigned int d = 0; d < dim; ++d)
               eigenvector[d] = 1.;
             for (unsigned int i = 0; i < 5; ++i)
               {
                 eigenvector = transpose(inverse_jacobian)
                               (inverse_jacobian eigenvector);
                 VectorizedArrayType eigenvector_norm = 0.;
                 for (unsigned int d = 0; d < dim; ++d)
                   eigenvector_norm =
                     std::max(eigenvector_norm, std::abs(eigenvector[d]));
                 eigenvector /= eigenvector_norm;
               }
             const auto jac_times_ev   = inverse_jacobian eigenvector;
             const auto max_eigenvalue = std::sqrt(
               (jac_times_ev jac_times_ev) / (eigenvector eigenvector));
             local_max =
               std::max(local_max,
                        max_eigenvalue speed_of_sound + convective_limit);
           }
* 
         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
              ++v)
           for (unsigned int d = 0; d < 3; ++d)
             max_transport = std::max(max_transport, local_max[v]);
       }
* 
     max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);
* 
     return max_transport;
   }
* 
 
* 
   template <int dim>
   class EulerProblem
   {
   public:
     EulerProblem();
* 
     void run();
* 
   private:
     void make_grid_and_dofs();
* 
     void output_results(const unsigned int result_number);
* 
     LinearAlgebra::distributed::Vector<Number> solution;
* 
     ConditionalOStream pcout;
* 
 #ifdef DEAL_II_WITH_P4EST
     parallel::distributed::Triangulation<dim> triangulation;
 #else
     Triangulation<dim> triangulation;
 #endif
* 
     FESystem<dim>        fe;
     MappingQGeneric<dim> mapping;
     DoFHandler<dim>      dof_handler;
* 
     TimerOutput timer;
* 
     EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator;
* 
     double time, time_step;
* 
     class Postprocessor : public DataPostprocessor<dim>
     {
     public:
       Postprocessor();
* 
       virtual void evaluate_vector_field(
         const DataPostprocessorInputs::Vector<dim> &inputs,
         std::vector<Vector<double>> &computed_quantities) const override;
* 
       virtual std::vector<std::string> get_names() const override;
* 
       virtual std::vector<
         DataComponentInterpretation::DataComponentInterpretation>
       get_data_component_interpretation() const override;
* 
       virtual UpdateFlags get_needed_update_flags() const override;
* 
     private:
       const bool do_schlieren_plot;
     };
   };
* 
 
* 
   template <int dim>
   EulerProblem<dim>::Postprocessor::Postprocessor()
     : do_schlieren_plot(dim == 2)
   {}
* 
 
* 
   template <int dim>
   void EulerProblem<dim>::Postprocessor::evaluate_vector_field(
     const DataPostprocessorInputs::Vector<dim> &inputs,
     std::vector<Vector<double>> &               computed_quantities) const
   {
     const unsigned int n_evaluation_points = inputs.solution_values.size();
* 
     if (do_schlieren_plot == true)
       Assert(inputs.solution_gradients.size() == n_evaluation_points,
              ExcInternalError());
* 
     Assert(computed_quantities.size() == n_evaluation_points,
            ExcInternalError());
     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
     Assert(computed_quantities[0].size() ==
              dim + 2 + (do_schlieren_plot == true ? 1 : 0),
            ExcInternalError());
* 
     for (unsigned int q = 0; q < n_evaluation_points; ++q)
       {
         Tensor<1, dim + 2> solution;
         for (unsigned int d = 0; d < dim + 2; ++d)
           solution[d] = inputs.solution_values[q](d);
* 
         const double         density  = solution[0];
         const Tensor<1, dim> velocity = euler_velocity<dim>(solution);
         const double         pressure = euler_pressure<dim>(solution);
* 
         for (unsigned int d = 0; d < dim; ++d)
           computed_quantities[q](d) = velocity[d];
         computed_quantities[q](dim)     = pressure;
         computed_quantities[q](dim + 1) = std::sqrt(gamma pressure / density);
* 
         if (do_schlieren_plot == true)
           computed_quantities[q](dim + 2) =
             inputs.solution_gradients[q][0] inputs.solution_gradients[q][0];
       }
   }
* 
 
* 
   template <int dim>
   std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const
   {
     std::vector<std::string> names;
     for (unsigned int d = 0; d < dim; ++d)
       names.emplace_back("velocity");
     names.emplace_back("pressure");
     names.emplace_back("speed_of_sound");
* 
     if (do_schlieren_plot == true)
       names.emplace_back("schlieren_plot");
* 
     return names;
   }
* 
 
* 
   template <int dim>
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
   EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const
   {
     std::vector<DataComponentInterpretation::DataComponentInterpretation>
       interpretation;
     for (unsigned int d = 0; d < dim; ++d)
       interpretation.push_back(
         DataComponentInterpretation::component_is_part_of_vector);
     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
     interpretation.push_back(DataComponentInterpretation::component_is_scalar);
* 
     if (do_schlieren_plot == true)
       interpretation.push_back(
         DataComponentInterpretation::component_is_scalar);
* 
     return interpretation;
   }
* 
 
* 
   template <int dim>
   UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const
   {
     if (do_schlieren_plot == true)
       return update_values | update_gradients;
     else
       return update_values;
   }
* 
 
* 
   template <int dim>
   EulerProblem<dim>::EulerProblem()
     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
 #ifdef DEAL_II_WITH_P4EST
     , triangulation(MPI_COMM_WORLD)
 #endif
     , fe(FE_DGQ<dim>(fe_degree), dim + 2)
     , mapping(fe_degree)
     , dof_handler(triangulation)
     , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
     , euler_operator(timer)
     , time(0)
     , time_step(0)
   {}
* 
 
* 
   template <int dim>
   void EulerProblem<dim>::make_grid_and_dofs()
   {
     switch (testcase)
       {
         case 0:
           {
             Point<dim> lower_left;
             for (unsigned int d = 1; d < dim; ++d)
               lower_left[d] =
* 
-5;
* 
             Point<dim> upper_right;
             upper_right[0] = 10;
             for (unsigned int d = 1; d < dim; ++d)
               upper_right[d] = 5;
* 
             GridGenerator::hyper_rectangle(triangulation,
                                            lower_left,
                                            upper_right);
             triangulation.refine_global(2);
* 
             euler_operator.set_inflow_boundary(
               0, std::make_unique<ExactSolution<dim>>(0));
* 
             break;
           }
* 
         case 1:
           {
             GridGenerator::channel_with_cylinder(
               triangulation, 0.03, 1, 0, true);
* 
             euler_operator.set_inflow_boundary(
               0, std::make_unique<ExactSolution<dim>>(0));
             euler_operator.set_subsonic_outflow_boundary(
               1, std::make_unique<ExactSolution<dim>>(0));
* 
             euler_operator.set_wall_boundary(2);
             euler_operator.set_wall_boundary(3);
* 
             if (dim == 3)
               euler_operator.set_body_force(
                 std::make_unique<Functions::ConstantFunction<dim>>(
                   std::vector<double>({0., 0.,
* 
-0.2})));
* 
             break;
           }
* 
         default:
           Assert(false, ExcNotImplemented());
       }
* 
     triangulation.refine_global(n_global_refinements);
* 
     dof_handler.distribute_dofs(fe);
* 
     euler_operator.reinit(mapping, dof_handler);
     euler_operator.initialize_vector(solution);
* 
     std::locale s = pcout.get_stream().getloc();
     pcout.get_stream().imbue(std::locale(""));
     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
           << " ( = " << (dim + 2) << " [vars] x "
           << triangulation.n_global_active_cells() << " [cells] x "
           << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )"
           << std::endl;
     pcout.get_stream().imbue(s);
   }
* 
 
* 
   template <int dim>
   void EulerProblem<dim>::output_results(const unsigned int result_number)
   {
     const std::array<double, 3> errors =
       euler_operator.compute_errors(ExactSolution<dim>(time), solution);
     const std::string quantity_name = testcase == 0 ? "error" : "norm";
* 
     pcout << "Time:" << std::setw(8) << std::setprecision(3) << time
           << ", dt: " << std::setw(8) << std::setprecision(2) << time_step
           << ", " << quantity_name << " rho: " << std::setprecision(4)
           << std::setw(10) << errors[0] << ", rho u: " << std::setprecision(4)
           << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4)
           << std::setw(10) << errors[2] << std::endl;
* 
     {
       TimerOutput::Scope t(timer, "output");
* 
       Postprocessor postprocessor;
       DataOut<dim>  data_out;
* 
       DataOutBase::VtkFlags flags;
       flags.write_higher_order_cells = true;
       data_out.set_flags(flags);
* 
       data_out.attach_dof_handler(dof_handler);
       {
         std::vector<std::string> names;
         names.emplace_back("density");
         for (unsigned int d = 0; d < dim; ++d)
           names.emplace_back("momentum");
         names.emplace_back("energy");
* 
         std::vector<DataComponentInterpretation::DataComponentInterpretation>
           interpretation;
         interpretation.push_back(
           DataComponentInterpretation::component_is_scalar);
         for (unsigned int d = 0; d < dim; ++d)
           interpretation.push_back(
             DataComponentInterpretation::component_is_part_of_vector);
         interpretation.push_back(
           DataComponentInterpretation::component_is_scalar);
* 
         data_out.add_data_vector(dof_handler, solution, names, interpretation);
       }
       data_out.add_data_vector(solution, postprocessor);
* 
       LinearAlgebra::distributed::Vector<Number> reference;
       if (testcase == 0 && dim == 2)
         {
           reference.reinit(solution);
           euler_operator.project(ExactSolution<dim>(time), reference);
           reference.sadd(-1., 1, solution);
           std::vector<std::string> names;
           names.emplace_back("error_density");
           for (unsigned int d = 0; d < dim; ++d)
             names.emplace_back("error_momentum");
           names.emplace_back("error_energy");
* 
           std::vector<DataComponentInterpretation::DataComponentInterpretation>
             interpretation;
           interpretation.push_back(
             DataComponentInterpretation::component_is_scalar);
           for (unsigned int d = 0; d < dim; ++d)
             interpretation.push_back(
               DataComponentInterpretation::component_is_part_of_vector);
           interpretation.push_back(
             DataComponentInterpretation::component_is_scalar);
* 
           data_out.add_data_vector(dof_handler,
                                    reference,
                                    names,
                                    interpretation);
         }
* 
       Vector<double> mpi_owner(triangulation.n_active_cells());
       mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
       data_out.add_data_vector(mpi_owner, "owner");
* 
       data_out.build_patches(mapping,
                              fe.degree,
                              DataOut<dim>::curved_inner_cells);
* 
       const std::string filename =
         "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu";
       data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
     }
   }
* 
 
* 
   template <int dim>
   void EulerProblem<dim>::run()
   {
     {
       const unsigned int n_vect_number = VectorizedArrayType::size();
       const unsigned int n_vect_bits   = 8 sizeof(Number) n_vect_number;
* 
       pcout << "Running with "
             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
             << " MPI processes" << std::endl;
       pcout << "Vectorization over " << n_vect_number << " "
             << (std::is_same<Number, double>::value ? "doubles" : "floats")
             << " = " << n_vect_bits << " bits ("
             << Utilities::System::get_current_vectorization_level() << ")"
             << std::endl;
     }
* 
     make_grid_and_dofs();
* 
     const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme);
* 
     LinearAlgebra::distributed::Vector<Number> rk_register_1;
     LinearAlgebra::distributed::Vector<Number> rk_register_2;
     rk_register_1.reinit(solution);
     rk_register_2.reinit(solution);
* 
     euler_operator.project(ExactSolution<dim>(time), solution);
* 
 
     double min_vertex_distance = std::numeric_limits<double>::max();
     for (const auto &cell : triangulation.active_cell_iterators())
       if (cell->is_locally_owned())
         min_vertex_distance =
           std::min(min_vertex_distance, cell->minimum_vertex_distance());
     min_vertex_distance =
       Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);
* 
     time_step = courant_number integrator.n_stages() /
                 euler_operator.compute_cell_transport_speed(solution);
     pcout << "Time step size: " << time_step
           << ", minimal h: " << min_vertex_distance
           << ", initial transport scaling: "
           << 1. / euler_operator.compute_cell_transport_speed(solution)
           << std::endl
           << std::endl;
* 
     output_results(0);
* 
     unsigned int timestep_number = 0;
* 
     while (time < final_time
* 
- 1e-12 && timestep_number < max_time_steps)
       {
         ++timestep_number;
         if (timestep_number % 5 == 0)
           time_step =
             courant_number integrator.n_stages() /
             Utilities::truncate_to_n_digits(
               euler_operator.compute_cell_transport_speed(solution), 3);
* 
         {
           TimerOutput::Scope t(timer, "rk time stepping total");
           integrator.perform_time_step(euler_operator,
                                        time,
                                        time_step,
                                        solution,
                                        rk_register_1,
                                        rk_register_2);
         }
* 
         time += time_step;
* 
         if (static_cast<int>(time / output_tick) !=
               static_cast<int>((time
* 
- time_step) / output_tick) ||
             time >= final_time
* 
- 1e-12)
           output_results(
             static_cast<unsigned int>(std::round(time / output_tick)));
       }
* 
     timer.print_wall_time_statistics(MPI_COMM_WORLD);
     pcout << std::endl;
   }
* 
 } // namespace Euler_DG
* 
 
 int main(int argc, char*argv)
 {
   using namespace Euler_DG;
   using namespace dealii;
* 
   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
* 
   try
     {
       deallog.depth_console(0);
* 
       EulerProblem<dimension> euler_problem;
       euler_problem.run();
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
* 
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
* 
   return 0;
 }
 @endcode
* <a name="Results"></a><h1>Results</h1>
* 

* Running the program with the default settings on a machine with 40 processesproduces the following output:
* @code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 27.648.000 ( = 5 [vars] x 25.600 [cells] x 216 [dofs/cell/var] )
Time step size: 0.000295952, minimal h: 0.0075, initial transport scaling: 0.00441179
Time:       0, dt:   0.0003, norm rho:  5.385e-16, rho u:  1.916e-16, energy: 1.547e-15
+--------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed         |     17.52s    10 |     17.52s |     17.52s    11 |
|                                      |                  |                               |
| Section                  | no. calls |   min time  rank |   avg time |   max time  rank |
+--------------------------------------+------------------+------------+------------------+
| compute errors           |         1 |  0.009594s    16 |  0.009705s |  0.009819s     8 |
| compute transport speed  |        22 |    0.1366s     0 |    0.1367s |    0.1368s    18 |
| output                   |         1 |     1.233s     0 |     1.233s |     1.233s    32 |
| rk time stepping total   |       100 |     8.746s    35 |     8.746s |     8.746s     0 |
| rk_stage
* 
- integrals L_h |       500 |     8.742s    36 |     8.742s |     8.743s     2 |
+--------------------------------------+------------------+------------+------------------+
@endcode
* 
* and the following visual output:
*   <table align="center" class="doxtable" style="width:85%">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_010.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_025.png" alt="" width="100%">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_050.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_100.png" alt="" width="100%">
    </td>
  </tr>
</table>  
* As a reference, the results of   step-67   using FCL are:
* @code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 27.648.000 ( = 5 [vars] x 25.600 [cells] x 216 [dofs/cell/var] )
Time step size: 0.000295952, minimal h: 0.0075, initial transport scaling: 0.00441179
Time:       0, dt:   0.0003, norm rho:  5.385e-16, rho u:  1.916e-16, energy: 1.547e-15
+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     13.33s     0 |     13.34s |     13.35s    34 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |         1 |  0.007977s    10 |  0.008053s |  0.008161s    30 |
| compute transport speed       |        22 |    0.1228s    34 |    0.2227s |    0.3845s     0 |
| output                        |         1 |     1.255s     3 |     1.257s |     1.259s    27 |
| rk time stepping total        |       100 |     11.15s     0 |     11.32s |     11.42s    34 |
| rk_stage
* 
- integrals L_h      |       500 |     8.719s    10 |     8.932s |     9.196s     0 |
| rk_stage
* 
- inv mass + vec upd |       500 |     1.944s     0 |     2.377s |      2.55s    10 |
+-------------------------------------------+------------------+------------+------------------+
@endcode
* 
* By the modifications shown in this tutorial, we were able to achieve a speedup of27% for the Runge-Kutta stages.
* <a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>
* 

* The algorithms are easily extendable to higher dimensions: a high-dimensional<a href="https://github.com/hyperdeal/hyperdeal/blob/a9e67b4e625ff1dde2fed93ad91cdfacfaa3acdf/include/hyper.deal/operators/advection/advection_operation.h#L219-L569">advection operator based on cell-centric loops</a>is part of the hyper.deal library. An extension of cell-centric loopsto locally-refined meshes is more involved.
* <a name="ExtensiontothecompressibleNavierStokesequations"></a><h4>Extension to the compressible Navier-Stokes equations</h4>
* 

* The solver presented in this tutorial program can also be extended to thecompressible Navier–Stokes equations by adding viscous terms, as alsosuggested in   step-67  . To keep as much of the performance obtained here despitethe additional cost of elliptic terms, e.g. via an interior penalty method, thattutorial has proposed to switch the basis from FE_DGQ to FE_DGQHermite like inthe   step-59   tutorial program. The reasoning behind this switch is that in thecase of FE_DGQ all values of neighboring cells (i.e.,   $k+1$   layers) are needed,whilst in the case of FE_DGQHermite only 2 layers, making the lattersignificantly more suitable for higher degrees. The additional layers have to be,on the one hand, loaded from main memory during flux computation and, one theother hand, have to be communicated. Using the shared-memory capabilitiesintroduced in this tutorial, the second point can be eliminated on a singlecompute node or its influence can be reduced in a hybrid context.
* <a name="BlockGaussSeidellikepreconditioners"></a><h4>Block Gauss-Seidel-like preconditioners</h4>
* 

* Cell-centric loops could be used to create block Gauss-Seidel preconditionersthat are multiplicative within one process and additive over processes. Thesetype of preconditioners use during flux computation, in contrast to Jacobi-typepreconditioners, already updated values from neighboring cells. The followingpseudo-code visualizes how this could in principal be achieved:
* @code
// vector monitor if cells have been updated or not
Vector<Number> visit_flags(data.n_cell_batches () + data.n_ghost_cell_batches ());
* 
// element centric loop with a modified kernel
data.template loop_cell_centric<VectorType, VectorType>(
  [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
* 
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        // cell integral as usual (not shown)
* 
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          {
            const auto boundary_id = data.get_faces_by_cells_boundary_id(cell, face)[0];
* 
            if (boundary_id == numbers::internal_face_boundary_id)
              {
                phi_p.reinit(cell, face);
* 
                const auto flags = phi_p.read_cell_data(visit_flags);
                const auto all_neighbors_have_been_updated =
                  std::min(flags.begin(),
                           flags().begin() + data.n_active_entries_per_cell_batch(cell) == 1;
* 
                if(all_neighbors_have_been_updated)
                  phi_p.gather_evaluate(dst, EvaluationFlags::values);
                else
                  phi_p.gather_evaluate(src, EvaluationFlags::values);
* 
                // continue as usual (not shown)
              }
            else
              {
                // boundary integral as usual (not shown)
              }
          }
* 
        // continue as above and apply your favorite algorithm to invert
        // the cell-local operator (not shown)
* 
        // make cells as updated
        phi.set_cell_data(visit_flags, VectorizedArrayType(1.0));
      }
  },
  dst,
  src,
  true,
  MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values);
@endcode
* 
* For this purpose, one can exploit the cell-data vector capabilities ofMatrixFree and the range-based iteration capabilities of VectorizedArray.
* Please note that in the given example we process   <code>VectorizedArrayType::size()</code>  number of blocks, since each lane corresponds to one block. We consider blocksas updated if all blocks processed by a vector register have been updated. Inthe case of Cartesian meshes this is a reasonable approach, however, forgeneral unstructured meshes this conservative approach might lead to a decrease in theefficiency of the preconditioner. A reduction of cells processed in parallelby explicitly reducing the number of lanes used by   <code>VectorizedArrayType</code>  might increase the quality of the preconditioner, but with the cost that eachiteration might be more expensive. This dilemma leads us to a further"possibility for extension": vectorization within an element.
* 

* <a name="PlainProg"></a><h1> The plain program</h1>  @include "step-76.cc"  
* */
