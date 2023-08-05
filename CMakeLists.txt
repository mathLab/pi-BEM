##
# CMake script for large project, with one executable
# and a library of objects, optionally
# depending on deal2lkit.
#
# The structure of the directory is assumed to be:
# ./source: containing implementations and one file with "int main()"
# ./include: containing all class declarations
# ./tests: containing pairs of files: test.cc/test.output to test
#
# If you set the environemnt variable D2K_DIR,
# everything will work out of the box
##

# Set the name of the project and target
# If your application follows the structure above, you don't need to
# specify anything else.
SET(TARGET bem_fma)
SET(_main source/main.cc)

# Set the _main variable to empty if you don't want an executable
# but only a library
#SET(_main "")

############################################################
# Normally you shouldn't need to change anything below.
############################################################
# Declare all source files the target consists of:
file(GLOB _files source/*cc)

# Don't compile the main file into the library
IF(NOT "${_main}" STREQUAL "")
    LIST(REMOVE_ITEM _files
      ${CMAKE_CURRENT_SOURCE_DIR}/${_main}
    )
ENDIF()


CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(deal.II 9.5 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING 
       "Choose the type of build, options are: Debug, Release or DebugRelease"
       )
ELSE()
   SET(CMAKE_BUILD_TYPE ${DEAL_II_BUILD_TYPE} CACHE STRING
       "Choose the type of build, options are: Debug, Release or DebugRelease"
       )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()


PROJECT(${TARGET})

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set(CMAKE_MACOSX_RPATH 1)
endif()

OPTION(BEM_ENABLE_TESTING 
"Enable deal.II style tests in this application." ON)


# We one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
    # Postfix to use everywhere
    SET(_p "${${_build_type}_postfix}")
    # Only build this type, if deal.II was compiled with it.
    IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
       DEAL_II_BUILD_TYPE MATCHES "${_build_type}")
	MESSAGE("-- Found ${_build_type} version of deal.II.")

	SET(_lib "${TARGET}-lib${_p}")
	MESSAGE("-- Configuring library ${_lib}")

	STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

	# Add all files, except the main file, to the library
	ADD_LIBRARY(${_lib} SHARED ${_files})
  	FOREACH(_dim 2 3)
		SET(_exe "${TARGET}_${_dim}d${${_build_type}_postfix}")
		MESSAGE("-- Configuring executable ${_exe}")
		ADD_EXECUTABLE(${_exe} ${_main})
		target_compile_definitions(${_exe} PUBLIC DEAL_II_DIMENSION=${_dim})
		TARGET_LINK_LIBRARIES(${_exe} ${_lib})
		DEAL_II_SETUP_TARGET(${_exe} ${_BUILD_TYPE})
    	ENDFOREACH()



	DEAL_II_SETUP_TARGET(${_lib} ${_BUILD_TYPE})
  SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})
  ENDIF()
ENDFOREACH()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests AND BEM_ENABLE_TESTING)
  ADD_SUBDIRECTORY(tests)
  ENABLE_TESTING()
endif()


# Doxygen
if(EXISTS ${CMAKE_SOURCE_DIR}/docs/Doxyfile)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/docs/Doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
    endif(DOXYGEN_FOUND)
endif()


message(STATUS "DEAL_II_INCLUDE_DIRECTORIES=${deal.II_INCLUDE_DIRECTORIES}")

ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./scripts/indent
)
