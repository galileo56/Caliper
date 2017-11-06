# Install script for directory: /Users/vmateu/GitHub/Caliper

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/Caliper")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/Caliper")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Caliper" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Caliper")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Caliper")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND install_name_tool -change @executable_path/../Frameworks/mathlink.framework/Versions/4.38/mathlink /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions/mathlink.framework/mathlink /Users/vmateu/GitHub/Caliper/bin/Caliper OUTPUT_QUIET)
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND install_name_tool -change @executable_path/../Frameworks/mathlink.framework/Versions/4.36/mathlink /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions/mathlink.framework/mathlink /Users/vmateu/GitHub/Caliper/bin/Caliper OUTPUT_QUIET)
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/hello-world")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/hello-world")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/hello-world" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/hello-world")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/hello-world")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/Chi2MbNRQCD")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/Chi2MbNRQCD")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2MbNRQCD" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2MbNRQCD")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2MbNRQCD")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCD")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/MbFinderNRQCD")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCD" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCD")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCD")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/NRQCDGrid")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/NRQCDGrid")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/NRQCDGrid" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/NRQCDGrid")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/NRQCDGrid")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/Chi2NRQCDAnalyzer")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/Chi2NRQCDAnalyzer")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2NRQCDAnalyzer" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2NRQCDAnalyzer")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/Chi2NRQCDAnalyzer")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCDAnalyzer")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/vmateu/GitHub/Caliper/bin" TYPE EXECUTABLE FILES "/Users/vmateu/GitHub/Caliper/build/MbFinderNRQCDAnalyzer")
  if(EXISTS "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCDAnalyzer" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCDAnalyzer")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/vmateu/GitHub/Caliper/bin/MbFinderNRQCDAnalyzer")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/vmateu/GitHub/Caliper/build/lib/cmake_install.cmake")
  include("/Users/vmateu/GitHub/Caliper/build/libC/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/vmateu/GitHub/Caliper/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
