## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------


########################################################################
##                                                                    ##
##               The deal.II project configuration file               ##
##                                                                    ##
########################################################################


#
# General information
#

SET(DEAL_II_PACKAGE_NAME "deal.II")
SET(DEAL_II_PACKAGE_VERSION "10.0.0-pre")
SET(DEAL_II_PACKAGE_VENDOR "The deal.II Authors <http://www.dealii.org/>")
SET(DEAL_II_PACKAGE_DESCRIPTION "Library for solving partial differential equations with the finite element method")

SET(DEAL_II_VERSION_MAJOR "10")
SET(DEAL_II_VERSION_MINOR "0")
SET(DEAL_II_VERSION_SUBMINOR "0")
SET(DEAL_II_VERSION "10.0.0")

SET(DEAL_II_GIT_BRANCH "master")
SET(DEAL_II_GIT_REVISION "8311ff0c68b3b71d2fe3cfb8950e17e823a18e2d")
SET(DEAL_II_GIT_SHORTREV "8311ff0c68")
SET(DEAL_II_GIT_TAG "v9.3.0")

SET(DEAL_II_PROJECT_CONFIG_NAME "deal.II")

SET(DEAL_II_BUILD_TYPE "DebugRelease")
SET(DEAL_II_BUILD_TYPES "DEBUG;RELEASE")


#
# Information about the project location
#

SET(DEAL_II_DOCHTML_RELDIR "doc")
SET(DEAL_II_DOCREADME_RELDIR "./")
SET(DEAL_II_EXAMPLES_RELDIR "examples")
SET(DEAL_II_EXECUTABLE_RELDIR "bin")
SET(DEAL_II_INCLUDE_RELDIR "include")
SET(DEAL_II_LIBRARY_RELDIR "lib")
SET(DEAL_II_PROJECT_CONFIG_RELDIR "lib/cmake/deal.II")
SET(DEAL_II_SHARE_RELDIR "share/deal.II")

#
# Determine DEAL_II_PATH from CMAKE_CURRENT_LIST_DIR:
#

SET(DEAL_II_PATH "${CMAKE_CURRENT_LIST_DIR}")
SET(_path "${DEAL_II_PROJECT_CONFIG_RELDIR}")
WHILE(NOT "${_path}" STREQUAL "")
  GET_FILENAME_COMPONENT(DEAL_II_PATH "${DEAL_II_PATH}" PATH)
  GET_FILENAME_COMPONENT(_path "${_path}" PATH)
ENDWHILE()

#
# Print a message after inclusion of this file:
#

SET(DEAL_II_PROJECT_CONFIG_INCLUDED TRUE)

SET(DEAL_II_BUILD_DIR TRUE)

IF(NOT ${DEAL_II_PACKAGE_NAME}_FIND_QUIETLY)
  IF(DEAL_II_BUILD_DIR)
    MESSAGE(STATUS
      "Using the ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION} build directory found at ${DEAL_II_PATH}"
      )
  ELSE()
    MESSAGE(STATUS
      "Using the ${DEAL_II_PACKAGE_NAME}-${DEAL_II_PACKAGE_VERSION} installation found at ${DEAL_II_PATH}"
      )
  ENDIF()
ENDIF()


#
# Include all convenience macros:
#

FILE(GLOB _macro_files
  "${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/macros/*.cmake"
  )
FOREACH(file ${_macro_files})
  IF(NOT ${DEAL_II_PACKAGE_NAME}_FIND_QUIETLY)
    MESSAGE(STATUS "Include macro ${file}")
  ENDIF()
  INCLUDE(${file})
ENDFOREACH()


#
# Compiler and linker configuration
#

SET(DEAL_II_CXX_COMPILER "/usr/bin/c++")
SET(DEAL_II_C_COMPILER "/usr/bin/cc")

# used for all targets:
SET(DEAL_II_CXX_FLAGS "-pedantic -fPIC -Wall -Wextra -Wmissing-braces -Woverloaded-virtual -Wpointer-arith -Wsign-compare  -Wswitch -Wsynth -Wwrite-strings -Wno-placement-new  -Wno-literal-suffix -Wno-psabi -Wno-class-memaccess -fopenmp-simd -Wno-unused-local-typedefs")

# _additionally_ used for debug targets:
SET(DEAL_II_CXX_FLAGS_DEBUG "-O0 -ggdb -Wa,--compress-debug-sections")

# _additionally_ used for release targets:
SET(DEAL_II_CXX_FLAGS_RELEASE "-O2 -funroll-loops -funroll-all-loops -fstrict-aliasing -Wno-unused-local-typedefs")

# used for all targets:
SET(DEAL_II_LINKER_FLAGS "-rdynamic -fuse-ld=gold -lpthread")

# _additionally_ used for debug targets:
SET(DEAL_II_LINKER_FLAGS_DEBUG "-ggdb")

# _additionally_ used for release targets:
SET(DEAL_II_LINKER_FLAGS_RELEASE "")

# used for all targets:
SET(DEAL_II_USER_DEFINITIONS "BOOST_NO_AUTO_PTR")

# _additionally_ used for debug targets:
SET(DEAL_II_USER_DEFINITIONS_DEBUG "DEBUG")

# _additionally_ used for release targets:
SET(DEAL_II_USER_DEFINITIONS_RELEASE "")

#
# MPI runtime:
#

SET(DEAL_II_MPIEXEC "")
SET(DEAL_II_MPIEXEC_NUMPROC_FLAG "")
SET(DEAL_II_MPIEXEC_PREFLAGS "")
SET(DEAL_II_MPIEXEC_POSTFLAGS "")

#
# CUDA specific setup:
#
SET(CMAKE_CUDA_ARCHITECTURES "")

#
# Build a static executable:
#

SET(DEAL_II_STATIC_EXECUTABLE "OFF")


#
# Information about include directories and libraries
#

# Full list of include directories:
SET(DEAL_II_INCLUDE_DIRS "${DEAL_II_PATH}/include;/workspaces/dealii/include/;/workspaces/dealii/bundled/taskflow-2.5.0/include;/usr/include;/usr/include/suitesparse;/usr/include/opencascade")

# Full list of libraries for the debug target:
SET(DEAL_II_LIBRARIES_DEBUG "${DEAL_II_PATH}/lib/libdeal_II.g.so;/usr/lib/aarch64-linux-gnu/libtbb.so;/usr/lib/aarch64-linux-gnu/libz.so;/usr/lib/aarch64-linux-gnu/libboost_iostreams.so;/usr/lib/aarch64-linux-gnu/libboost_serialization.so;/usr/lib/aarch64-linux-gnu/libboost_system.so;/usr/lib/aarch64-linux-gnu/libboost_thread.so;/usr/lib/aarch64-linux-gnu/libboost_regex.so;/usr/lib/aarch64-linux-gnu/libboost_chrono.so;/usr/lib/aarch64-linux-gnu/libboost_date_time.so;/usr/lib/aarch64-linux-gnu/libboost_atomic.so;/usr/lib/aarch64-linux-gnu/libumfpack.so;/usr/lib/aarch64-linux-gnu/libcholmod.so;/usr/lib/aarch64-linux-gnu/libccolamd.so;/usr/lib/aarch64-linux-gnu/libcolamd.so;/usr/lib/aarch64-linux-gnu/libcamd.so;/usr/lib/aarch64-linux-gnu/libsuitesparseconfig.so;/usr/lib/aarch64-linux-gnu/libamd.so;/usr/lib/aarch64-linux-gnu/libmetis.so;rt;/usr/lib/aarch64-linux-gnu/libarpack.so;/usr/lib/aarch64-linux-gnu/liblapack.so;/usr/lib/aarch64-linux-gnu/libblas.so;/usr/lib/aarch64-linux-gnu/libassimp.so;/usr/lib/aarch64-linux-gnu/libgmsh.so;/usr/lib/aarch64-linux-gnu/libgsl.so;/usr/lib/aarch64-linux-gnu/libgslcblas.so;/usr/lib/aarch64-linux-gnu/libmuparser.so;/usr/lib/aarch64-linux-gnu/libTKBO.so;/usr/lib/aarch64-linux-gnu/libTKBool.so;/usr/lib/aarch64-linux-gnu/libTKBRep.so;/usr/lib/aarch64-linux-gnu/libTKernel.so;/usr/lib/aarch64-linux-gnu/libTKFeat.so;/usr/lib/aarch64-linux-gnu/libTKFillet.so;/usr/lib/aarch64-linux-gnu/libTKG2d.so;/usr/lib/aarch64-linux-gnu/libTKG3d.so;/usr/lib/aarch64-linux-gnu/libTKGeomAlgo.so;/usr/lib/aarch64-linux-gnu/libTKGeomBase.so;/usr/lib/aarch64-linux-gnu/libTKHLR.so;/usr/lib/aarch64-linux-gnu/libTKIGES.so;/usr/lib/aarch64-linux-gnu/libTKMath.so;/usr/lib/aarch64-linux-gnu/libTKMesh.so;/usr/lib/aarch64-linux-gnu/libTKOffset.so;/usr/lib/aarch64-linux-gnu/libTKPrim.so;/usr/lib/aarch64-linux-gnu/libTKShHealing.so;/usr/lib/aarch64-linux-gnu/libTKSTEP.so;/usr/lib/aarch64-linux-gnu/libTKSTEPAttr.so;/usr/lib/aarch64-linux-gnu/libTKSTEPBase.so;/usr/lib/aarch64-linux-gnu/libTKSTEP209.so;/usr/lib/aarch64-linux-gnu/libTKSTL.so;/usr/lib/aarch64-linux-gnu/libTKTopAlgo.so;/usr/lib/aarch64-linux-gnu/libTKXSBase.so;/usr/lib/aarch64-linux-gnu/libsundials_idas.so;/usr/lib/aarch64-linux-gnu/libsundials_arkode.so;/usr/lib/aarch64-linux-gnu/libsundials_kinsol.so;/usr/lib/aarch64-linux-gnu/libsundials_nvecserial.so")

# Full list of libraries for the release target:
SET(DEAL_II_LIBRARIES_RELEASE "${DEAL_II_PATH}/lib/libdeal_II.so;/usr/lib/aarch64-linux-gnu/libtbb.so;/usr/lib/aarch64-linux-gnu/libz.so;/usr/lib/aarch64-linux-gnu/libboost_iostreams.so;/usr/lib/aarch64-linux-gnu/libboost_serialization.so;/usr/lib/aarch64-linux-gnu/libboost_system.so;/usr/lib/aarch64-linux-gnu/libboost_thread.so;/usr/lib/aarch64-linux-gnu/libboost_regex.so;/usr/lib/aarch64-linux-gnu/libboost_chrono.so;/usr/lib/aarch64-linux-gnu/libboost_date_time.so;/usr/lib/aarch64-linux-gnu/libboost_atomic.so;/usr/lib/aarch64-linux-gnu/libumfpack.so;/usr/lib/aarch64-linux-gnu/libcholmod.so;/usr/lib/aarch64-linux-gnu/libccolamd.so;/usr/lib/aarch64-linux-gnu/libcolamd.so;/usr/lib/aarch64-linux-gnu/libcamd.so;/usr/lib/aarch64-linux-gnu/libsuitesparseconfig.so;/usr/lib/aarch64-linux-gnu/libamd.so;/usr/lib/aarch64-linux-gnu/libmetis.so;rt;/usr/lib/aarch64-linux-gnu/libarpack.so;/usr/lib/aarch64-linux-gnu/liblapack.so;/usr/lib/aarch64-linux-gnu/libblas.so;/usr/lib/aarch64-linux-gnu/libassimp.so;/usr/lib/aarch64-linux-gnu/libgmsh.so;/usr/lib/aarch64-linux-gnu/libgsl.so;/usr/lib/aarch64-linux-gnu/libgslcblas.so;/usr/lib/aarch64-linux-gnu/libmuparser.so;/usr/lib/aarch64-linux-gnu/libTKBO.so;/usr/lib/aarch64-linux-gnu/libTKBool.so;/usr/lib/aarch64-linux-gnu/libTKBRep.so;/usr/lib/aarch64-linux-gnu/libTKernel.so;/usr/lib/aarch64-linux-gnu/libTKFeat.so;/usr/lib/aarch64-linux-gnu/libTKFillet.so;/usr/lib/aarch64-linux-gnu/libTKG2d.so;/usr/lib/aarch64-linux-gnu/libTKG3d.so;/usr/lib/aarch64-linux-gnu/libTKGeomAlgo.so;/usr/lib/aarch64-linux-gnu/libTKGeomBase.so;/usr/lib/aarch64-linux-gnu/libTKHLR.so;/usr/lib/aarch64-linux-gnu/libTKIGES.so;/usr/lib/aarch64-linux-gnu/libTKMath.so;/usr/lib/aarch64-linux-gnu/libTKMesh.so;/usr/lib/aarch64-linux-gnu/libTKOffset.so;/usr/lib/aarch64-linux-gnu/libTKPrim.so;/usr/lib/aarch64-linux-gnu/libTKShHealing.so;/usr/lib/aarch64-linux-gnu/libTKSTEP.so;/usr/lib/aarch64-linux-gnu/libTKSTEPAttr.so;/usr/lib/aarch64-linux-gnu/libTKSTEPBase.so;/usr/lib/aarch64-linux-gnu/libTKSTEP209.so;/usr/lib/aarch64-linux-gnu/libTKSTL.so;/usr/lib/aarch64-linux-gnu/libTKTopAlgo.so;/usr/lib/aarch64-linux-gnu/libTKXSBase.so;/usr/lib/aarch64-linux-gnu/libsundials_idas.so;/usr/lib/aarch64-linux-gnu/libsundials_arkode.so;/usr/lib/aarch64-linux-gnu/libsundials_kinsol.so;/usr/lib/aarch64-linux-gnu/libsundials_nvecserial.so")

# Full list of libraries with "debug" and "optimized" keywords for easy use with TARGET_LINK_LIBRARIES:
SET(DEAL_II_LIBRARIES "debug;${DEAL_II_LIBRARIES_DEBUG};optimized;${DEAL_II_LIBRARIES_RELEASE}")


#
# Information about CUDA configuration
#

SET(DEAL_II_CUDA_TOOLKIT_ROOT_DIR "")
SET(DEAL_II_CUDA_COMPILER "")

# used for all cuda targets:
SET(DEAL_II_CUDA_FLAGS "")

# _additionally_ used for debug targets:
SET(DEAL_II_CUDA_FLAGS_DEBUG "")

# _additionally_ used for release targets:
SET(DEAL_II_CUDA_FLAGS_RELEASE "")


#
# Information about library targets and feature configuration
#

# The library targets file:
SET(DEAL_II_TARGET_CONFIG "${DEAL_II_PATH}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake")

# The Debug target:
SET(DEAL_II_TARGET_DEBUG "deal_II.g")

# The Release target:
SET(DEAL_II_TARGET_RELEASE "deal_II")

# Full list of targets with "debug" and "optimized" keywords for easy use with TARGET_LINK_LIBRARIES:
SET(DEAL_II_TARGET "debug;${DEAL_II_TARGET_DEBUG};optimized;${DEAL_II_TARGET_RELEASE}")

# The feature configuration file:
SET(DEAL_II_FEATURE_CONFIG "${DEAL_II_PATH}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}FeatureConfig.cmake")


#
# Feature configuration:
#

SET(DEAL_II_WITH_CXX11 ON)
SET(DEAL_II_WITH_CXX14 ON)
SET(DEAL_II_WITH_CXX17 FALSE)
SET(DEAL_II_WITH_THREADS ON)
SET(DEAL_II_WITH_64BIT_INDICES OFF)
SET(DEAL_II_WITH_ADOLC OFF)
SET(DEAL_II_WITH_ARBORX OFF)
SET(DEAL_II_WITH_ARPACK ON)
SET(DEAL_II_ARPACK_WITH_PARPACK FALSE)
SET(DEAL_II_WITH_ASSIMP ON)
SET(DEAL_II_WITH_BOOST ON)
SET(DEAL_II_BOOST_VERSION "1.71.0")
SET(DEAL_II_WITH_COMPLEX_VALUES ON)
SET(DEAL_II_WITH_CUDA OFF)
SET(DEAL_II_WITH_GINKGO OFF)
SET(DEAL_II_WITH_GMSH ON)
SET(DEAL_II_GMSH_WITH_API TRUE)
SET(DEAL_II_WITH_GSL ON)
SET(DEAL_II_GSL_VERSION "2.5")
SET(DEAL_II_WITH_HDF5 OFF)
SET(DEAL_II_WITH_KOKKOS OFF)
SET(DEAL_II_WITH_LAPACK ON)
SET(DEAL_II_LAPACK_WITH_MKL OFF)
SET(DEAL_II_WITH_METIS ON)
SET(DEAL_II_METIS_VERSION "5.1.0")
SET(DEAL_II_WITH_MPI OFF)
SET(DEAL_II_WITH_MUPARSER ON)
SET(DEAL_II_MUPARSER_VERSION "2.2.6")
SET(DEAL_II_WITH_OPENCASCADE ON)
SET(DEAL_II_OPENCASCADE_VERSION "7.3.0")
SET(DEAL_II_WITH_P4EST OFF)
SET(DEAL_II_WITH_PETSC OFF)
SET(DEAL_II_WITH_SCALAPACK OFF)
SET(DEAL_II_WITH_SLEPC OFF)
SET(DEAL_II_WITH_SUNDIALS ON)
SET(DEAL_II_SUNDIALS_VERSION "3.1.2")
SET(DEAL_II_SUNDIALS_WITH_IDAS TRUE)
SET(DEAL_II_WITH_SYMENGINE OFF)
SET(DEAL_II_WITH_TASKFLOW ON)
SET(DEAL_II_WITH_TBB ON)
SET(DEAL_II_TBB_VERSION "2020.1")
SET(DEAL_II_WITH_TRILINOS OFF)
SET(DEAL_II_WITH_UMFPACK ON)
SET(DEAL_II_UMFPACK_VERSION "5.7.8")
SET(DEAL_II_WITH_ZLIB ON)
SET(DEAL_II_ZLIB_VERSION "1.2.11")
