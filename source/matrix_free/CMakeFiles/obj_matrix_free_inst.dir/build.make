# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /workspaces/dealii

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workspaces/dealii

# Utility rule file for obj_matrix_free_inst.

# Include the progress variables for this target.
include source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/progress.make

source/matrix_free/CMakeFiles/obj_matrix_free_inst: source/matrix_free/evaluation_template_factory.inst
source/matrix_free/CMakeFiles/obj_matrix_free_inst: source/matrix_free/fe_point_evaluation.inst
source/matrix_free/CMakeFiles/obj_matrix_free_inst: source/matrix_free/mapping_info.inst
source/matrix_free/CMakeFiles/obj_matrix_free_inst: source/matrix_free/matrix_free.inst
source/matrix_free/CMakeFiles/obj_matrix_free_inst: source/matrix_free/shape_info.inst


source/matrix_free/evaluation_template_factory.inst: bin/expand_instantiations
source/matrix_free/evaluation_template_factory.inst: share/deal.II/template-arguments
source/matrix_free/evaluation_template_factory.inst: source/matrix_free/evaluation_template_factory.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating evaluation_template_factory.inst"
	cd /workspaces/dealii/source/matrix_free && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/matrix_free/evaluation_template_factory.inst.in > /workspaces/dealii/source/matrix_free/evaluation_template_factory.inst.tmp
	cd /workspaces/dealii/source/matrix_free && /usr/bin/cmake -E rename /workspaces/dealii/source/matrix_free/evaluation_template_factory.inst.tmp /workspaces/dealii/source/matrix_free/evaluation_template_factory.inst

source/matrix_free/fe_point_evaluation.inst: bin/expand_instantiations
source/matrix_free/fe_point_evaluation.inst: share/deal.II/template-arguments
source/matrix_free/fe_point_evaluation.inst: source/matrix_free/fe_point_evaluation.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating fe_point_evaluation.inst"
	cd /workspaces/dealii/source/matrix_free && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/matrix_free/fe_point_evaluation.inst.in > /workspaces/dealii/source/matrix_free/fe_point_evaluation.inst.tmp
	cd /workspaces/dealii/source/matrix_free && /usr/bin/cmake -E rename /workspaces/dealii/source/matrix_free/fe_point_evaluation.inst.tmp /workspaces/dealii/source/matrix_free/fe_point_evaluation.inst

source/matrix_free/mapping_info.inst: bin/expand_instantiations
source/matrix_free/mapping_info.inst: share/deal.II/template-arguments
source/matrix_free/mapping_info.inst: source/matrix_free/mapping_info.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating mapping_info.inst"
	cd /workspaces/dealii/source/matrix_free && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/matrix_free/mapping_info.inst.in > /workspaces/dealii/source/matrix_free/mapping_info.inst.tmp
	cd /workspaces/dealii/source/matrix_free && /usr/bin/cmake -E rename /workspaces/dealii/source/matrix_free/mapping_info.inst.tmp /workspaces/dealii/source/matrix_free/mapping_info.inst

source/matrix_free/matrix_free.inst: bin/expand_instantiations
source/matrix_free/matrix_free.inst: share/deal.II/template-arguments
source/matrix_free/matrix_free.inst: source/matrix_free/matrix_free.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating matrix_free.inst"
	cd /workspaces/dealii/source/matrix_free && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/matrix_free/matrix_free.inst.in > /workspaces/dealii/source/matrix_free/matrix_free.inst.tmp
	cd /workspaces/dealii/source/matrix_free && /usr/bin/cmake -E rename /workspaces/dealii/source/matrix_free/matrix_free.inst.tmp /workspaces/dealii/source/matrix_free/matrix_free.inst

source/matrix_free/shape_info.inst: bin/expand_instantiations
source/matrix_free/shape_info.inst: share/deal.II/template-arguments
source/matrix_free/shape_info.inst: source/matrix_free/shape_info.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Generating shape_info.inst"
	cd /workspaces/dealii/source/matrix_free && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/matrix_free/shape_info.inst.in > /workspaces/dealii/source/matrix_free/shape_info.inst.tmp
	cd /workspaces/dealii/source/matrix_free && /usr/bin/cmake -E rename /workspaces/dealii/source/matrix_free/shape_info.inst.tmp /workspaces/dealii/source/matrix_free/shape_info.inst

obj_matrix_free_inst: source/matrix_free/CMakeFiles/obj_matrix_free_inst
obj_matrix_free_inst: source/matrix_free/evaluation_template_factory.inst
obj_matrix_free_inst: source/matrix_free/fe_point_evaluation.inst
obj_matrix_free_inst: source/matrix_free/mapping_info.inst
obj_matrix_free_inst: source/matrix_free/matrix_free.inst
obj_matrix_free_inst: source/matrix_free/shape_info.inst
obj_matrix_free_inst: source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/build.make

.PHONY : obj_matrix_free_inst

# Rule to build all files generated by this target.
source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/build: obj_matrix_free_inst

.PHONY : source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/build

source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/clean:
	cd /workspaces/dealii/source/matrix_free && $(CMAKE_COMMAND) -P CMakeFiles/obj_matrix_free_inst.dir/cmake_clean.cmake
.PHONY : source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/clean

source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/depend:
	cd /workspaces/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/dealii /workspaces/dealii/source/matrix_free /workspaces/dealii /workspaces/dealii/source/matrix_free /workspaces/dealii/source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/matrix_free/CMakeFiles/obj_matrix_free_inst.dir/depend

