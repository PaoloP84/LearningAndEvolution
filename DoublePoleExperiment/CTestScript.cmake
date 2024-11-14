# This script is run by CTest on test machines. It simply performs a compilation
# of FARSA on the host and submits results to CDash. This script doesn't
# perform the update of the repository as that is done externally. When using
# this script, pass the host as argument (either Linux, Windows or MacOSX). For
# example:
# 	ctest -S /mnt/sharedRepos/farsa/farsa/CTestScript.cmake,Linux

# The build name is the same on all sites
SET( CTEST_BUILD_NAME "CartPoleExperiment" )
# The build configuration to use is always Release
SET( CTEST_BUILD_CONFIGURATION "Release" )
# Not using memcheck nor code coverage, for the moment
SET( WITH_MEMCHECK FALSE )
SET( WITH_COVERAGE FALSE )
# Setting the site name
SITE_NAME( CTEST_SITE )

# Setting properties that are site-dependent
IF( CTEST_SCRIPT_ARG STREQUAL "Linux" )
	# Setting source and binary dir
	SET( CTEST_SOURCE_DIRECTORY "/mnt/sharedRepos/farsa/farsaPlugins/${CTEST_BUILD_NAME}" )
	SET( CTEST_BINARY_DIRECTORY "/home/nursery/build-${CTEST_BUILD_NAME}" )

	# Setting the generator to use
	SET( CTEST_CMAKE_GENERATOR "Unix Makefiles" )
	# Setting the build options
	SET( CTEST_BUILD_OPTIONS "\"-DCMAKE_CXX_FLAGS=-Wall -Wextra\" \"-DCMAKE_C_FLAGS=-Wall -Wextra\"" )
ELSEIF( CTEST_SCRIPT_ARG STREQUAL "Windows" )
	# Setting source and binary dir
	SET( CTEST_SOURCE_DIRECTORY "E:\\farsa\\farsaPlugins\\${CTEST_BUILD_NAME}" )
	SET( CTEST_BINARY_DIRECTORY "C:\\farsa\\build-${CTEST_BUILD_NAME}" )

	# Setting the generator to use
	SET( CTEST_CMAKE_GENERATOR "Visual Studio 10" )
	# Setting the build options
	SET( CTEST_BUILD_OPTIONS "" )
ELSEIF( CTEST_SCRIPT_ARG STREQUAL "MacOSX" )
	MESSAGE( FATAL_ERROR "TODO MacOSX" )
ELSE( CTEST_SCRIPT_ARG STREQUAL "Linux" )
	MESSAGE( FATAL_ERROR "Unknown system" )
ENDIF( CTEST_SCRIPT_ARG STREQUAL "Linux" )

# Clearing the build directory
CTEST_EMPTY_BINARY_DIRECTORY( ${CTEST_BINARY_DIRECTORY} )

# Generating the configure command
SET( CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION} -DBUILD_TESTING:BOOL=ON ${CTEST_BUILD_OPTIONS} \"-G${CTEST_CMAKE_GENERATOR}\" \"${CTEST_SOURCE_DIRECTORY}\"")

# The steps to execute. We skip the update (done externally)
CTEST_START( "Nightly" )
CTEST_CONFIGURE()
CTEST_BUILD()
CTEST_TEST()
IF ( WITH_COVERAGE AND CTEST_COVERAGE_COMMAND )
	CTEST_COVERAGE()
ENDIF ( WITH_COVERAGE AND CTEST_COVERAGE_COMMAND )
IF ( WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND )
	CTEST_MEMCHECK()
ENDIF ( WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND )
CTEST_SUBMIT()
