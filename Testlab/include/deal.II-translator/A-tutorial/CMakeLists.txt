include/deal.II-translator/A-tutorial/CMakeLists.txt
## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II Authors
##
## [This file is part of the deal.II library.
##
## [The deal.II library is free software; you can use it, redistribute
## ＃＃它，和/或修改它的条款的GNU较小的通用
## 自由软件基金会发布的GNU小通用许可证；或者
## 许可证的2.1版，或者（你可以选择）任何后来的版本。
## 许可证的全文可以在LICENSE.md文件中找到，该文件位于
## ＃＃＃deal.II的最高级别目录。
##
## ---------------------------------------------------------------------



#收集所有的教程程序的目录名称
FILE(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )

# 同时收集所有代码库项目的名称。为了
# 这样做，找到所有 "作者 "文件，然后去掉这些路径的最后两层
# 这些路径的最后两层。
#
# 由于不清楚的原因，glob返回这些文件为
# "/a/b/c/name/doc/author"，所以要确保我们在第二步中吃掉
# 第二步中的双斜线
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
FILE(GLOB _code_gallery_names
     "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}" )

#
# 定义教程的目标。它取决于
# 文件 tutorial.h 通过下面的目标建立，以及
# 以及我们从教程中创建的各种文件
# 下面的目录。这些依赖性被添加到
# 在各自的目标下面。
#
# 这个文件使用 DEAL_II_STEPS 变量，设置在
# .../CMakeLists.txt.
#

ADD_CUSTOM_TARGET(tutorial)

#
# 描述如何构建 tutorial.h。
#
# 首先收集这一切所依赖的各种文件。这
# 对于教程来说比较简单，因为目录名称
# 是结构化的。对于代码库，我们已经
# 确定了相关子目录的名称列表
#（可能是空的），所以我们可以创建文件列表
# 以下目标所依赖的文件列表，不需要GLOB，而是
# 只用一个 FOREACH 循环。
#
# 对于教程程序，step.pl脚本将读取
# tooltip, kind, 和 builds-on 文件。对于代码库
# 程序，"种类 "是显而易见的，我们只需要其他的
# 两个类别。
#

file(GLOB _deal_ii_steps_tooltip
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/tooltip
  )
FOREACH(_gallery ${_code_gallery_names})
  LIST(APPEND _deal_ii_code_gallery_tooltip
       ${_gallery}/doc/tooltip)
ENDFOREACH()
file(GLOB _deal_ii_steps_kind
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/kind
  )
文件(GLOB _deal_ii_steps_buildson
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/builds-on
  )
FOREACH(_gallery ${_code_gallery_names})
  LIST(APPEND _deal_ii_code_gallery_buildson
       ${_gallery}/doc/builds-on)
ENDFOREACH()


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  命令 ${perl_executable}
  许可证
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps}
    ${_code_gallery_names}。
    > ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  依附于
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps_tooltip}.
    ${_deal_ii_steps_kind}
    ${_deal_ii_steps_buildson} ${_deal_ii_steps_buildon}.
    ${_deal_ii_code_gallery_tooltip}
    ${_deal_ii_code_gallery_buildson}。
  )
ADD_CUSTOM_TARGET(build_tutorial_h
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
    评论
      "构建tutorial.h")
ADD_DEPENDENCIES(tutorial build_tutorial_h)


#
# 准备好生成文档的步骤
#

FOREACH(_step ${_deal_ii_steps})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)

  IF(EXISTS "${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cu" )
    SET(_FILE_EXTENSION "cu")
  ELSE()
    设置(_FILE_EXTENSION "cc")
  ENDIF()

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
    命令 ${perl_executable}.
    许可证
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      < ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
    隶属于
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}。
    脚本
    )

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    命令 ${perl_executable}.
    许可证
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${_step} ${CMAKE_SOURCE_DIR}
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    工作目录
      ${cmake_current_binary_dir}
    隶属于
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/intro2toc
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/create_anchors
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2doxygen
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}。
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/intro.dox
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/results.dox
    )

  ADD_CUSTOM_TARGET(tutorial_${_step})
    隶属于
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
      评论
        "为教程程序<${_step}>构建doxygen输入文件"
    )
  ADD_DEPENDENCIES(tutorial tutorial_${_step})
ENDFOREACH()
## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II Authors
##
##这个文件是deal.II库的一部分。
##
## [The deal.II library is free software; you can use it, redistribute
## ＃＃它，和/或修改它的条款的GNU较小的通用
## 自由软件基金会发布的GNU小通用许可证；或者
## 许可证的2.1版，或者（你可以选择）任何后来的版本。
## 许可证的全文可以在LICENSE.md文件中找到，该文件位于
## ＃＃＃deal.II的最高级别目录。
##
## ---------------------------------------------------------------------



#收集所有的教程程序的目录名称
FILE(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )

# 同时收集所有代码库项目的名称。为了
# 这样做，找到所有 "作者 "文件，然后去掉这些路径的最后两层
# 这些路径的最后两层。
#
# 由于不清楚的原因，glob返回这些文件为
# "/a/b/c/name/doc/author"，所以要确保我们在第二步中吃掉
# 第二步中的双斜线
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
FILE(GLOB _code_gallery_names
     "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}" )

#
# 定义教程的目标。它取决于
# 文件 tutorial.h 通过下面的目标建立，以及
# 以及我们从教程中创建的各种文件
# 下面的目录。这些依赖性被添加到
# 在各自的目标下面。
#
# 这个文件使用 DEAL_II_STEPS 变量，设置在
# .../CMakeLists.txt.
#

ADD_CUSTOM_TARGET(tutorial)

#
# 描述如何构建 tutorial.h。
#
# 首先收集这一切所依赖的各种文件。这
# 对于教程来说比较简单，因为目录名称
# 是结构化的。对于代码库，我们已经
# 确定了相关子目录的名称列表
#（可能是空的），所以我们可以创建文件列表
# 以下目标所依赖的文件列表，不需要GLOB，而是
# 只用一个 FOREACH 循环。
#
# 对于教程程序，step.pl脚本将读取
# tooltip, kind, 和 builds-on 文件。对于代码库
# 程序，"种类 "是显而易见的，我们只需要其他的
# 两个类别。
#

file(GLOB _deal_ii_steps_tooltip
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/tooltip
  )
FOREACH(_gallery ${_code_gallery_names})
  LIST(APPEND _deal_ii_code_gallery_tooltip
       ${_gallery}/doc/tooltip)
ENDFOREACH()
file(GLOB _deal_ii_steps_kind
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/kind
  )
文件(GLOB _deal_ii_steps_buildson
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/builds-on
  )
FOREACH(_gallery ${_code_gallery_names})
  LIST(APPEND _deal_ii_code_gallery_buildson
       ${_gallery}/doc/builds-on)
ENDFOREACH()


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  命令 ${perl_executable}
  许可证
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps}
    ${_code_gallery_names}。
    > ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  依附于
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps_tooltip}.
    ${_deal_ii_steps_kind}
    ${_deal_ii_steps_buildson} ${_deal_ii_steps_buildon}.
    ${_deal_ii_code_gallery_tooltip}
    ${_deal_ii_code_gallery_buildson}。
  )
ADD_CUSTOM_TARGET(build_tutorial_h
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
    评论
      "构建tutorial.h")
ADD_DEPENDENCIES(tutorial build_tutorial_h)


#
# 准备好生成文档的步骤
#

FOREACH(_step ${_deal_ii_steps})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)

  IF(EXISTS "${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cu" )
    SET(_FILE_EXTENSION "cu")
  ELSE()
    设置(_FILE_EXTENSION "cc")
  ENDIF()

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
    命令 ${perl_executable}.
    许可证
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      < ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
    隶属于
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}。
    脚本
    )

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    命令 ${perl_executable}.
    许可证
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${_step} ${CMAKE_SOURCE_DIR}
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    工作目录
      ${cmake_current_binary_dir}
    隶属于
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/intro2toc
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/create_anchors
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2doxygen
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.${_FILE_EXTENSION}。
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/intro.dox
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/results.dox
    )

  ADD_CUSTOM_TARGET(tutorial_${_step})
    隶属于
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_FILE_EXTENSION}。
      评论
        "为教程程序<${_step}>构建doxygen输入文件"
    )
  ADD_DEPENDENCIES(tutorial tutorial_${_step})
ENDFOREACH()


