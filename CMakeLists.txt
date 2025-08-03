# require a certain cmake version
cmake_minimum_required(VERSION 3.4)

# declare the project name
project(gammacomboCore)

# set project version
set(${PROJECT_NAME}_VERSION_MAJOR 1)
set(${PROJECT_NAME}_VERSION_MINOR 0)
set(${PROJECT_NAME}_VERSION_PATCH 0)

# install destinations can be passed via the command line: cmake
# -DLIBRARY_OUTPUT_DIRECTORY:PATH=.... -DINCLUDE_OUTPUT_DIRECTORY:PATH=....
# otherwise, set default install destinations
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
  set(LIBRARY_OUTPUT_DIRECTORY ./lib)
endif()
if(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
  set(INCLUDE_OUTPUT_DIRECTORY ./include)
endif()

# include cmake macros and settings
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(settings)
include(macros)
include(combiners)

# use ccache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# find FreeType
find_package(Freetype REQUIRED)
message("FREETYPE_LIBRARIES: " ${FREETYPE_LIBRARIES})
# MESSAGE("FREETYPE_LIBRARY_DIR: " ${FREETYPE_LIBRARY})

# Find ROOT ROOT 6.24.06 ensures that the compiler is at least C++17
find_package(ROOT 6.24.06 REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O3 ${ROOT_CXX_FLAGS}")
message("ROOT_INCLUDE_DIR: " ${ROOT_INCLUDE_DIR})
message("ROOT_LIBRARIES: " ${ROOT_LIBRARIES})
message("ROOT_GUI_LIBRARIES: " ${ROOT_GUI_LIBRARIES})
message("ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

# find Boost
set(Boost_NO_BOOST_CMAKE
    ON
    CACHE BOOL "Enable fix for FindBoost.cmake")
find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Boost is found at: ${Boost_INCLUDE_DIR}")
  include_directories(${Boost_INCLUDE_DIR})
endif()

# find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
  # add combiner modules to Doxygen
  set(COMBINER_MODULES_DOXYGEN_SOURCES)
  foreach(dir ${COMBINER_MODULES})
    set(COMBINER_MODULES_DOXYGEN_SOURCES
        "${COMBINER_MODULES_DOXYGEN_SOURCES}
            INPUT += ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include
            INPUT += ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/src")
  endforeach()
  # configure Doxygen config file
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  # add target
  add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
else()
  message(
    STATUS
      "Doxygen not installed or found -- API documentation can't be generated!")
endif(DOXYGEN_FOUND)

# ##############################################################################
#
# set includes
#
# ##############################################################################

# Set core include directory
set(CORE_INC_DIR ./core/include/)

# include directories
include_directories(BEFORE ${CORE_INC_DIR} ${ROOT_INCLUDE_DIR})
# declare the ROOT library path
link_directories(${ROOT_LIBRARY_DIR} /usr/local/lib)

# define the sources of the core library
aux_source_directory(./core/src CORE_LIB_SOURCES)

# ##############################################################################
#
# generate ROOT dictionaries
#
# ##############################################################################

set(CORE_DICTIONARY_SOURCES
    RooBinned2DBicubicBase.h
    RooCrossCorPdf.h
    # RooHistInterpol.h
    RooHistPdfAngleVar.h
    RooHistPdfVar.h
    RooPoly3Var.h
    RooPoly4Var.h
    RooSlimFitResult.h
    RooMultiPdf.h)
# add a make rule for the dictionary
root_generate_dictionary(${PROJECT_NAME}Dict ${CORE_DICTIONARY_SOURCES} LINKDEF
                         coreLinkDef.h)
# add a target calling the make rule for the dictionary
add_custom_target(${PROJECT_NAME}Dict DEPENDS ${PROJECT_NAME}Dict.cxx
                                              ${PROJECT_NAME}Dict.h)
# add the dictionary to the list of library sources
set(CORE_LIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Dict.cxx
                     ${CORE_LIB_SOURCES})

# ##############################################################################
#
# build core library
#
# ##############################################################################

add_shared_library(${PROJECT_NAME}Components ${CORE_LIB_SOURCES})
# add the ROOT dictionary to the list of dependencies
add_dependencies(${PROJECT_NAME}Components ${PROJECT_NAME}Dict)
# add ROOT libraries to the linker command
target_link_libraries(${PROJECT_NAME}Components ${ROOT_LIBRARIES}
                      ${FREETYPE_LIBRARIES})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Dict_rdict.pcm
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib)

# ##############################################################################
#
# subdirectorys for the combiners
#
# ##############################################################################

foreach(dir ${COMBINER_MODULES})
  add_subdirectory(${dir})
endforeach()
