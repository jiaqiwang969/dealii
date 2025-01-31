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

# Utility rule file for obj_opencascade_inst.

# Include the progress variables for this target.
include source/opencascade/CMakeFiles/obj_opencascade_inst.dir/progress.make

source/opencascade/CMakeFiles/obj_opencascade_inst: source/opencascade/manifold_lib.inst
source/opencascade/CMakeFiles/obj_opencascade_inst: source/opencascade/utilities.inst


source/opencascade/manifold_lib.inst: bin/expand_instantiations
source/opencascade/manifold_lib.inst: share/deal.II/template-arguments
source/opencascade/manifold_lib.inst: source/opencascade/manifold_lib.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating manifold_lib.inst"
	cd /workspaces/dealii/source/opencascade && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/opencascade/manifold_lib.inst.in > /workspaces/dealii/source/opencascade/manifold_lib.inst.tmp
	cd /workspaces/dealii/source/opencascade && /usr/bin/cmake -E rename /workspaces/dealii/source/opencascade/manifold_lib.inst.tmp /workspaces/dealii/source/opencascade/manifold_lib.inst

source/opencascade/utilities.inst: bin/expand_instantiations
source/opencascade/utilities.inst: share/deal.II/template-arguments
source/opencascade/utilities.inst: source/opencascade/utilities.inst.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating utilities.inst"
	cd /workspaces/dealii/source/opencascade && ../../bin/expand_instantiations /workspaces/dealii/share/deal.II/template-arguments < /workspaces/dealii/source/opencascade/utilities.inst.in > /workspaces/dealii/source/opencascade/utilities.inst.tmp
	cd /workspaces/dealii/source/opencascade && /usr/bin/cmake -E rename /workspaces/dealii/source/opencascade/utilities.inst.tmp /workspaces/dealii/source/opencascade/utilities.inst

obj_opencascade_inst: source/opencascade/CMakeFiles/obj_opencascade_inst
obj_opencascade_inst: source/opencascade/manifold_lib.inst
obj_opencascade_inst: source/opencascade/utilities.inst
obj_opencascade_inst: source/opencascade/CMakeFiles/obj_opencascade_inst.dir/build.make

.PHONY : obj_opencascade_inst

# Rule to build all files generated by this target.
source/opencascade/CMakeFiles/obj_opencascade_inst.dir/build: obj_opencascade_inst

.PHONY : source/opencascade/CMakeFiles/obj_opencascade_inst.dir/build

source/opencascade/CMakeFiles/obj_opencascade_inst.dir/clean:
	cd /workspaces/dealii/source/opencascade && $(CMAKE_COMMAND) -P CMakeFiles/obj_opencascade_inst.dir/cmake_clean.cmake
.PHONY : source/opencascade/CMakeFiles/obj_opencascade_inst.dir/clean

source/opencascade/CMakeFiles/obj_opencascade_inst.dir/depend:
	cd /workspaces/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/dealii /workspaces/dealii/source/opencascade /workspaces/dealii /workspaces/dealii/source/opencascade /workspaces/dealii/source/opencascade/CMakeFiles/obj_opencascade_inst.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/opencascade/CMakeFiles/obj_opencascade_inst.dir/depend

