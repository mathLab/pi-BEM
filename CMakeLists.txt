##
# CMake script for large project, with one executable
# and a library
##

# Set the name of the project and target:
SET(TARGET bem)

# Specify which file contains the main() function
# Set it to empty if you don't want an executable
SET(_main source/main.cc)
#SET(_main "")

file(GLOB _files source/*cc)

# Declare all source files the target consists of:
SET(TARGET_SRC
    ${_files}
  )

IF(NOT "${_main}" STREQUAL "") 
    LIST(REMOVE_ITEM _files
      ${CMAKE_CURRENT_SOURCE_DIR}/${_main}
    )
ENDIF()

INCLUDE_DIRECTORIES(include)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.3 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

INCLUDE_DIRECTORIES(${deal.II_INCLUDE_DIRS}/deal.II)

PROJECT(${TARGET})

IF(DEAL_II_BUILD_TYPE MATCHES "Release")
    ADD_LIBRARY(${TARGET}-lib SHARED ${_files})
    IF(NOT "${_main}" STREQUAL "") 
	ADD_EXECUTABLE(${TARGET} ${_main})
	TARGET_LINK_LIBRARIES(${TARGET} ${TARGET}-lib)
	DEAL_II_SETUP_TARGET(${TARGET} RELEASE)
    ENDIF()
    SET(TEST_LIBRARIES_RELEASE ${TARGET}-lib)
    DEAL_II_SETUP_TARGET(${TARGET}-lib RELEASE)
ENDIF()

IF(DEAL_II_BUILD_TYPE MATCHES "Debug")
    ADD_LIBRARY(${TARGET}-lib.g SHARED ${_files})
    IF(NOT "${_main}" STREQUAL "") 
	ADD_EXECUTABLE(${TARGET}.g ${_main})
	TARGET_LINK_LIBRARIES(${TARGET}.g ${TARGET}-lib.g)
	DEAL_II_SETUP_TARGET(${TARGET}.g DEBUG)
    ENDIF()
    SET(TEST_LIBRARIES_DEBUG ${TARGET}-lib.g)
    DEAL_II_SETUP_TARGET(${TARGET}-lib.g DEBUG)
ENDIF()

ADD_SUBDIRECTORY(tests)
ENABLE_TESTING()
