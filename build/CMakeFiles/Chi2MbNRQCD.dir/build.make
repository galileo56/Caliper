# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/vicent/GitHub/Caliper

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/vicent/GitHub/Caliper/build

# Include any dependencies generated for this target.
include CMakeFiles/Chi2MbNRQCD.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Chi2MbNRQCD.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Chi2MbNRQCD.dir/flags.make

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o: CMakeFiles/Chi2MbNRQCD.dir/flags.make
CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o: ../src/Chi2MbNRQCD.F90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/vicent/GitHub/Caliper/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/vicent/GitHub/Caliper/src/Chi2MbNRQCD.F90 -o CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.i"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/vicent/GitHub/Caliper/src/Chi2MbNRQCD.F90 > CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.i

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.s"
	/usr/local/bin/gfortran $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/vicent/GitHub/Caliper/src/Chi2MbNRQCD.F90 -o CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.s

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.requires:

.PHONY : CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.requires

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.provides: CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.requires
	$(MAKE) -f CMakeFiles/Chi2MbNRQCD.dir/build.make CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.provides.build
.PHONY : CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.provides

CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.provides.build: CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o


# Object files for target Chi2MbNRQCD
Chi2MbNRQCD_OBJECTS = \
"CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o"

# External object files for target Chi2MbNRQCD
Chi2MbNRQCD_EXTERNAL_OBJECTS =

Chi2MbNRQCD: CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o
Chi2MbNRQCD: CMakeFiles/Chi2MbNRQCD.dir/build.make
Chi2MbNRQCD: lib/liblibCaliper.a
Chi2MbNRQCD: libsubdir_mods.a
Chi2MbNRQCD: CMakeFiles/Chi2MbNRQCD.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/vicent/GitHub/Caliper/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable Chi2MbNRQCD"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Chi2MbNRQCD.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Chi2MbNRQCD.dir/build: Chi2MbNRQCD

.PHONY : CMakeFiles/Chi2MbNRQCD.dir/build

CMakeFiles/Chi2MbNRQCD.dir/requires: CMakeFiles/Chi2MbNRQCD.dir/src/Chi2MbNRQCD.F90.o.requires

.PHONY : CMakeFiles/Chi2MbNRQCD.dir/requires

CMakeFiles/Chi2MbNRQCD.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Chi2MbNRQCD.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Chi2MbNRQCD.dir/clean

CMakeFiles/Chi2MbNRQCD.dir/depend:
	cd /Users/vicent/GitHub/Caliper/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/vicent/GitHub/Caliper /Users/vicent/GitHub/Caliper /Users/vicent/GitHub/Caliper/build /Users/vicent/GitHub/Caliper/build /Users/vicent/GitHub/Caliper/build/CMakeFiles/Chi2MbNRQCD.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Chi2MbNRQCD.dir/depend

