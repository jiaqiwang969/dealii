// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


/**
 * @defgroup geodynamics The geodynamics demonstration suite
 * deal.II's   @ref Tutorial   "tutorial"
 * contains a set of programs that togetherform the geodynamics demonstration
 * suite. The idea of these programs is todemonstrate techniques for advanced
 * finite element software usingapplications from geodynamics, i.e. the
 * investigation of processes in thesolid earth. By doing so, these programs
 * are supposed to provide a basisfor more specialized, dedicated programs
 * that can solve actual geodynamicsproblems, for example as part of the work
 * of graduate students orpostdocs. A more thorough discussion of the
 * motivation for these programsfollows below. Currently, the geodynamics
 * testsuite contains the followingprograms:
 *
 *  -   step-8  : Elasticity
 *
 *  -   step-18  : A %parallel elasticity solver
 *
 *  -   step-20  : Porous media flow
 *
 *  -   step-21  : Multiphase flow through porous media
 *
 *  -   step-22  : Stokes flow
 *
 *  -   step-31  : Thermal convection (Boussinesq flow)
 *
 *  -   step-32  : A %parallel Boussinesq solver for mantle convection
 * Some of these programs were developed under contract from the
 * CaliforniaInstitute of Technology with support by the National Science
 * Foundationunder Award No. EAR-0426271, the first of the grants that
 * fundedthe <a target="_top" href="http://www.geodynamics.org">Computational
 * Infrastructure in Geodynamics</a> initiative. The recipient, Wolfgang
 * Bangerth, gratefullyacknowledges this source of support.
 *
 *  <h3>Rationale</h3> Adaptive mesh refinement (AMR) has long been identified
 * as a key technologythat would aid in the accurate and efficient numerical
 * solution of a number ofgeodynamics applications. It has been discussed in
 * the geodynamics communityfor several years and has been a continuous topic
 * on the task list of CIGsince its inception. Yet, relatively little has
 * happened in this direction sofar. Only recently have there been attempts to
 * use AMR in geodynamics: CIGsponsored a workshop on AMR technique in Boulder
 * in October 2007; acollaboration between George Biros, Omar Ghattas, Mike
 * Gurnis, and ShijieZhong's groups is currently developing a %parallel
 * adaptive mantle convectionsolver; and some of the principal developers of
 * deal.II eventually developedthe <a
 * href="https://aspect.geodynamics.org">ASPECT code</a> for the simulationof
 * mantle convection that is by now a rather established and widely usedcode.
 * One of the reasons for the slow adoption of AMR techniques in geodynamics
 * isthe relatively steep initial hurdle: codes have to provide the data
 * structuresand algorithms to deal with adaptive meshes, finite elements have
 * to be ableto deal with hanging nodes, etc. To do so efficiently and in
 * sufficientgenerality adds several 10,000 lines of code to finite element
 * programs, toomuch for the average student to do within the time frame of a
 * dissertation. Onthe other hand, there are libraries that provide the
 * infrastructure code onwhich applications supporting AMR can rapidly be
 * built. deal.IIof course provides exactly this infrastructure. The goal of
 * the geodynamics testsuite is to write programs for a variety oftopics
 * relevant to geodynamics. Continuing in the style of the existing
 * tutorialprograms
 *
 *  -  an extensive introduction explaining the background andformulation of an application as well as the concepts of the numerical schemeused in its solution; detailed comments throughout the code explainingimplementation details; and a section showing numerical results
 *
 * -  we intend toprovide the resulting programs as well-documented applications solving modelproblems. In particular, they are aimed at the following goals:  <ul>    <li>   <i>Starting points:</i> The existing tutorial of deal.II has proven to  be an excellent starting point for graduate students and researchers to  jump-start developing their own applications. By providing programs that are  already close to the targeted application, first results can often be  obtained very quickly, both maintaining the initial enthusiasm during  development as well as allowing to spend research time on implementing  application specific behavior rather than using months of work on basic  infrastructure code supporting AMR.
 * Supporting this point is the fact that although there are  <a
 * href="https://www.dealii.org/publications.html">more than 1,000
 * publications</a> presenting results obtained with deal.II, we are aware of
 * only a relatively small number of applications that have been built with
 * deal.II from  scratch; all others have started as modifications of one of
 * the tutorial  programs.
 * <li>   <i>Training:</i> The tutorial programs we propose to write will
 * provide students and researchers with a reference implementation of current
 * numerical technology such as AMR, higher order elements, sophisticated
 * linear and nonlinear solvers, stabilization techniques, etc. Providing
 * these  as starting points for further development by others will also serve
 * the  goal of training a new generation of geodynamicists in modern
 * numerical  algorithms.
 * <li>   <i>Extending equations and formulations:</i> In deal.II, it is
 * fairly  simple to extend a set of equations by another equation, for
 * example an  additional advected quantity that enters the existing equations
 * as a right  hand side or in one of the coefficients. Since applications
 * typically use  blocked matrices rather than the
 * one-big-matrix-for-everything approach, it  is also not complicated to find
 * suitable linear solvers for augmented  equations. Consequently, deal.II is
 * a good tool for trying out more complex  formulations of problems, or more
 * complete models and their effects on the  accuracy of solutions.
 * <li>   <i>Rapid prototyping and benchmarking:</i> deal.II provides many
 * interchangeable components that allow rapid prototyping of finite element
 * kinds and orders, stabilization techniques, or linear solvers. For example,
 * typically only a few lines of code have to be changed to replace low-order
 * by high-order elements. Through this, it becomes relatively simple to try
 * out higher order elements, a different block elimination solver, or a
 * different stabilization technique. In turn, this may help in benchmarking
 * applications both regarding computing times to solve as well as concerning
 * the accuracy of numerical solutions. The applications in this module will
 * already have been benchmarked for  correctness. Existing tutorial programs
 * typically employ simpler rather than  more complicated solver schemes for
 * exposition but frequently suggest more  complicated schemes including hints
 * on how they might be implemented in an  appendix.
 * <li>   <i>Try algorithms:</i> The rapid prototyping abilities of deal.II may  also help in determining best algorithms on the scale of programs to which  deal.II is applicable, and then to implement this particular algorithm  (without the ability to change it easily) in a dedicated program that can  run on larger scale machines. For example, a small mantle convection code  built on deal.II may be used to determine whether second order elements are  useful for this purpose (see, for example, the results shown in    step-31  ). If so, then one may use this kind of knowledge in larger codes,  such as the ASPECT code mentioned above.  </ul>
 *
 *
 *
 */
