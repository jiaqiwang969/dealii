# Install script for directory: /workspaces/dealii/cmake/config

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DebugRelease")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/deal.II/macros" TYPE FILE MESSAGE_LAZY FILES
    "/workspaces/dealii/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake"
    "/workspaces/dealii/cmake/macros/macro_deal_ii_invoke_autopilot.cmake"
    "/workspaces/dealii/cmake/macros/macro_deal_ii_setup_target.cmake"
    "/workspaces/dealii/cmake/macros/macro_deal_ii_query_git_information.cmake"
    "/workspaces/dealii/cmake/macros/macro_deal_ii_add_test.cmake"
    "/workspaces/dealii/cmake/macros/macro_deal_ii_pickup_tests.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/deal.II" TYPE FILE MESSAGE_LAZY FILES
    "/workspaces/dealii/cmake/config/deal.IIConfig.cmake"
    "/workspaces/dealii/cmake/config/deal.IIConfigVersion.cmake"
    "/workspaces/dealii/lib/cmake/deal.II/deal.IIFeatureConfig.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE MESSAGE_LAZY FILES "/workspaces/dealii/cmake/config/deal.II_debug.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE MESSAGE_LAZY FILES "/workspaces/dealii/cmake/config/deal.II_release.pc")
endif()

