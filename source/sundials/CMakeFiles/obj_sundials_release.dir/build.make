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

# Include any dependencies generated for this target.
include source/sundials/CMakeFiles/obj_sundials_release.dir/depend.make

# Include the progress variables for this target.
include source/sundials/CMakeFiles/obj_sundials_release.dir/progress.make

# Include the compile flags for this target's objects.
include source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make

source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.o: source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make
source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.o: source/sundials/arkode.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.o"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obj_sundials_release.dir/arkode.cc.o -c /workspaces/dealii/source/sundials/arkode.cc

source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obj_sundials_release.dir/arkode.cc.i"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/source/sundials/arkode.cc > CMakeFiles/obj_sundials_release.dir/arkode.cc.i

source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obj_sundials_release.dir/arkode.cc.s"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/source/sundials/arkode.cc -o CMakeFiles/obj_sundials_release.dir/arkode.cc.s

source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.o: source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make
source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.o: source/sundials/ida.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.o"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obj_sundials_release.dir/ida.cc.o -c /workspaces/dealii/source/sundials/ida.cc

source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obj_sundials_release.dir/ida.cc.i"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/source/sundials/ida.cc > CMakeFiles/obj_sundials_release.dir/ida.cc.i

source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obj_sundials_release.dir/ida.cc.s"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/source/sundials/ida.cc -o CMakeFiles/obj_sundials_release.dir/ida.cc.s

source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.o: source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make
source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.o: source/sundials/kinsol.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.o"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obj_sundials_release.dir/kinsol.cc.o -c /workspaces/dealii/source/sundials/kinsol.cc

source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obj_sundials_release.dir/kinsol.cc.i"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/source/sundials/kinsol.cc > CMakeFiles/obj_sundials_release.dir/kinsol.cc.i

source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obj_sundials_release.dir/kinsol.cc.s"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/source/sundials/kinsol.cc -o CMakeFiles/obj_sundials_release.dir/kinsol.cc.s

source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.o: source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make
source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.o: source/sundials/n_vector.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.o"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obj_sundials_release.dir/n_vector.cc.o -c /workspaces/dealii/source/sundials/n_vector.cc

source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obj_sundials_release.dir/n_vector.cc.i"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/source/sundials/n_vector.cc > CMakeFiles/obj_sundials_release.dir/n_vector.cc.i

source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obj_sundials_release.dir/n_vector.cc.s"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/source/sundials/n_vector.cc -o CMakeFiles/obj_sundials_release.dir/n_vector.cc.s

source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.o: source/sundials/CMakeFiles/obj_sundials_release.dir/flags.make
source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.o: source/sundials/sunlinsol_wrapper.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.o"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.o -c /workspaces/dealii/source/sundials/sunlinsol_wrapper.cc

source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.i"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/dealii/source/sundials/sunlinsol_wrapper.cc > CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.i

source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.s"
	cd /workspaces/dealii/source/sundials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/dealii/source/sundials/sunlinsol_wrapper.cc -o CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.s

obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/arkode.cc.o
obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/ida.cc.o
obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/kinsol.cc.o
obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/n_vector.cc.o
obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/sunlinsol_wrapper.cc.o
obj_sundials_release: source/sundials/CMakeFiles/obj_sundials_release.dir/build.make

.PHONY : obj_sundials_release

# Rule to build all files generated by this target.
source/sundials/CMakeFiles/obj_sundials_release.dir/build: obj_sundials_release

.PHONY : source/sundials/CMakeFiles/obj_sundials_release.dir/build

source/sundials/CMakeFiles/obj_sundials_release.dir/clean:
	cd /workspaces/dealii/source/sundials && $(CMAKE_COMMAND) -P CMakeFiles/obj_sundials_release.dir/cmake_clean.cmake
.PHONY : source/sundials/CMakeFiles/obj_sundials_release.dir/clean

source/sundials/CMakeFiles/obj_sundials_release.dir/depend:
	cd /workspaces/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/dealii /workspaces/dealii/source/sundials /workspaces/dealii /workspaces/dealii/source/sundials /workspaces/dealii/source/sundials/CMakeFiles/obj_sundials_release.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : source/sundials/CMakeFiles/obj_sundials_release.dir/depend

