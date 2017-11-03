cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(Caliper C CXX Fortran)

###############
## REQUIREMENTS & SETUP

## check whether this is a linux system
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	set(OperatingSystem "Linux")
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
	set(OperatingSystem "Mac OS X")
endif()
message(STATUS "Operating system: ${OperatingSystem}")

## set correct RPath for external (dynamic) libraries
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

## find mathematica on the system
## Sascha Kratky's FindMathematica cmake module: https://github.com/sakra/FindMathematica/ [accessed 2015-03-06 01:00 AM]
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/tools")
find_package(Mathematica REQUIRED)
message(STATUS "Mathematica Version found: " ${Mathematica_VERSION})

###############
## COMPILER FLAGS

## additional flags for c++ compiler - for c++ and clang++ (should also take effect for the C compiler)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(addition_cxx "-std=c++11")
## elseif ("" ... "") ## additional flags for other compilers
endif()

## additional flags for the fortran compiler - for gfortran
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(addition_fortran "-Waliasing -Wampersand -Wconversion -Wc-binding-type -Wintrinsics-std -Wtabs -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Winteger-division -Wreal-q-constant -Wunused")
endif()

## let the flag changes take effect
set(CMAKE_CXX_FLAGS "${addition_cxx}" CACHE STRING "Flags used by the compiler during all build types." FORCE)
set(CMAKE_Fortran_FLAGS "${addition_fortran}" CACHE STRING "Flags used by the compiler during all build types." FORCE)

###############
## LIBS

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(lib)
include_directories(include)

include_directories(${Library_MODDIR})
FILE(GLOB libSources lib/*.F90)
add_library(subdir_mods ${libSources})

## include the information from ./libC/CmakeLists.txt
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libC)

###############
## EXE & LINKING & INSTALL

message(STATUS "Checking for mathlink library")
if(Mathematica_MathLink_LIBRARY)
	message(STATUS "Checking for mathlink library -- done")
	include_directories(${Mathematica_MathLink_INCLUDE_DIR})
	find_library(M_LIB m)
	find_library(DL_LIB dl)
	message(STATUS "mathlink version: " ${Mathematica_MathLink_VERSION} )

	if("${OperatingSystem}" STREQUAL "Linux")
		find_library(UUID_LIB uuid)
		find_library(RT_LIB rt)
	endif()

	## build mathlink
	execute_process(COMMAND ${Mathematica_MathLink_MPREP_EXECUTABLE} -o ${CMAKE_CURRENT_SOURCE_DIR}/src/Caliper.tm.c ${CMAKE_CURRENT_SOURCE_DIR}/src/Caliper.tm OUTPUT_QUIET)

	add_executable(Caliper src/Caliper.tm.c)
	target_link_libraries(Caliper libCaliper ${Mathematica_MathLink_LIBRARY} ${M_LIB} ${DL_LIB} ${UUID_LIB} ${RT_LIB})

	else(Mathematica_MathLink_LIBRARY)
	message(STATUS "Checking for mathlink library -- not found, proceeding without mathlink")
endif()

add_executable(hello-world ${CMAKE_CURRENT_SOURCE_DIR}/src/hello-world.cpp)
target_link_libraries(hello-world libradreturn libCaliper)

	add_executable(Chi2MbNRQCD ${CMAKE_CURRENT_SOURCE_DIR}/src/Chi2MbNRQCD.F90)
	target_link_libraries(Chi2MbNRQCD libCaliper subdir_mods)

		add_executable(MbFinderNRQCD ${CMAKE_CURRENT_SOURCE_DIR}/src/MbFinderNRQCD.F90)
		target_link_libraries(MbFinderNRQCD libCaliper subdir_mods)

###############
## INSTALL

if(Mathematica_MathLink_LIBRARY)
	install(TARGETS Caliper DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

	if("${OperatingSystem}" STREQUAL "Mac OS X")
		install(CODE "execute_process(COMMAND install_name_tool -change @executable_path/../Frameworks/mathlink.framework/Versions/${Mathematica_MathLink_VERSION}/mathlink ${Mathematica_MathLink_LIBRARY}/mathlink ${CMAKE_CURRENT_SOURCE_DIR}/bin/Caliper OUTPUT_QUIET)")
		install(CODE "execute_process(COMMAND install_name_tool -change @executable_path/../Frameworks/mathlink.framework/Versions/4.36/mathlink ${Mathematica_MathLink_LIBRARY}/mathlink ${CMAKE_CURRENT_SOURCE_DIR}/bin/Caliper OUTPUT_QUIET)")
	endif()
endif()

install(TARGETS hello-world DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
install(TARGETS Chi2MbNRQCD DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
install(TARGETS MbFinderNRQCD DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
