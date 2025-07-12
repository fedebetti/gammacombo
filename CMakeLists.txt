cmake_minimum_required(VERSION 3.30)

project(GammaCombo VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#[[
Install destinations can be passed via the command line:

  cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=<lib-out-dir> -DINCLUDE_OUTPUT_DIRECTORY:PATH=<inc-out-dir>

otherwise, set default install destinations.
]]
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
  set(LIBRARY_OUTPUT_DIRECTORY ./lib)
endif()
if(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
  set(INCLUDE_OUTPUT_DIRECTORY ./include)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(settings)
include(macros)
include(combiners)

# Use ccache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

find_package(Freetype REQUIRED)
message("FREETYPE_LIBRARIES: " ${FREETYPE_LIBRARIES})

find_package(ROOT 6.34 REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O3 ${ROOT_CXX_FLAGS}")
message("ROOT_INCLUDE_DIR: " ${ROOT_INCLUDE_DIR})
message("ROOT_LIBRARIES: " ${ROOT_LIBRARIES})
message("ROOT_GUI_LIBRARIES: " ${ROOT_GUI_LIBRARIES})
message("ROOT_LIBRARY_DIR: " ${ROOT_LIBRARY_DIR})
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

set(Boost_NO_BOOST_CMAKE
    ON
    CACHE BOOL "Enable fix for FindBoost.cmake")
find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Boost is found at: ${Boost_INCLUDE_DIR}")
  include_directories(${Boost_INCLUDE_DIR})
endif()

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
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in
                 ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  add_custom_target(
    doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
else()
  message("Doxygen not found -- API documentation can't be generated!")
endif(DOXYGEN_FOUND)

# -----------------------------------------------------------------------------
# Set includes
# -----------------------------------------------------------------------------

set(CORE_INC_DIR ./core/include/)
include_directories(BEFORE ${CORE_INC_DIR} ${ROOT_INCLUDE_DIR})
# Declare the ROOT library path
link_directories(${ROOT_LIBRARY_DIR} /usr/local/lib)
configure_file(${CORE_INC_DIR}/VersionConfig.h.in VersionConfig.h)

# Define the sources of the core library
aux_source_directory(./core/src CORE_LIB_SOURCES)

# -----------------------------------------------------------------------------
# Generate ROOT dictionaries
# -----------------------------------------------------------------------------

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
# Add a make rule for the dictionary
root_generate_dictionary(${PROJECT_NAME}Dict ${CORE_DICTIONARY_SOURCES} LINKDEF
                         coreLinkDef.h)
# Add a target calling the make rule for the dictionary
add_custom_target(${PROJECT_NAME}Dict DEPENDS ${PROJECT_NAME}Dict.cxx
                                              ${PROJECT_NAME}Dict.h)
# Add the dictionary to the list of library sources
set(CORE_LIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Dict.cxx
                     ${CORE_LIB_SOURCES})

# -----------------------------------------------------------------------------
# Build core library
# -----------------------------------------------------------------------------

add_shared_library(${PROJECT_NAME}Components ${CORE_LIB_SOURCES})
# Add the ROOT dictionary to the list of dependencies
add_dependencies(${PROJECT_NAME}Components ${PROJECT_NAME}Dict)
# Add ROOT libraries to the linker command
target_link_libraries(${PROJECT_NAME}Components ${ROOT_LIBRARIES}
                      ${FREETYPE_LIBRARIES})
# For ConfigVersion.h
target_include_directories(${PROJECT_NAME}Components
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Dict_rdict.pcm
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/lib)

# -----------------------------------------------------------------------------
# Subdirectorys for the combiners
# -----------------------------------------------------------------------------

set(OUTPUT_FOLDERS
    root
    plots
    plots/C
    plots/cl
    plots/dot
    plots/eps
    plots/par
    plots/pdf
    plots/png
    plots/root
    plots/scanner)

foreach(dir ${COMBINER_MODULES})
  add_subdirectory(${dir})
endforeach()
