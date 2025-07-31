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

# https://root.cern/manual/integrate_root_into_my_cmake_project/
# https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
find_package(ROOT 6.34 CONFIG REQUIRED)
message("ROOT is found at: ${ROOT_INCLUDE_DIRS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O3 ${ROOT_CXX_FLAGS}")
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

find_package(Boost REQUIRED)
message("Boost is found at: ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})

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
include_directories(BEFORE ${CORE_INC_DIR})
configure_file(${CORE_INC_DIR}/VersionConfig.h.in VersionConfig.h)

# Define the sources of the core library
set(CORE_LIB_SOURCES
    ./core/src/BatchScriptWriter.cpp
    ./core/src/CLInterval.cpp
    ./core/src/CLIntervalMaker.cpp
    ./core/src/CLIntervalPrinter.cpp
    ./core/src/ColorBuilder.cpp
    ./core/src/Combiner.cpp
    ./core/src/ConfidenceContours.cpp
    ./core/src/Contour.cpp
    ./core/src/ControlPlots.cpp
    ./core/src/FileNameBuilder.cpp
    ./core/src/FitResultCache.cpp
    ./core/src/Fitter.cpp
    ./core/src/GammaComboEngine.cpp
    ./core/src/Graphviz.cpp
    ./core/src/LatexMaker.cpp
    ./core/src/MethodAbsScan.cpp
    ./core/src/MethodBergerBoosScan.cpp
    ./core/src/MethodCoverageScan.cpp
    ./core/src/MethodDatasetsPluginScan.cpp
    ./core/src/MethodDatasetsProbScan.cpp
    ./core/src/MethodPluginScan.cpp
    ./core/src/MethodProbScan.cpp
    ./core/src/OneMinusClPlot2d.cpp
    ./core/src/OneMinusClPlotAbs.cpp
    ./core/src/OneMinusClPlot.cpp
    ./core/src/OptParser.cpp
    ./core/src/ParameterCache.cpp
    ./core/src/Parameter.cpp
    ./core/src/ParameterEvolutionPlotter.cpp
    ./core/src/ParametersAbs.cpp
    ./core/src/PDF_Abs.cpp
    ./core/src/PDF_Datasets.cpp
    ./core/src/ProgressBar.cpp
    ./core/src/PullPlotter.cpp
    ./core/src/PValueCorrection.cpp
    ./core/src/RooBinned2DBicubicBase.cpp
    ./core/src/RooCrossCorPdf.cpp
    ./core/src/RooHistPdfAngleVar.cpp
    ./core/src/RooHistPdfVar.cpp
    ./core/src/RooMultiPdf.cpp
    ./core/src/RooPoly3Var.cpp
    ./core/src/RooPoly4Var.cpp
    ./core/src/RooSlimFitResult.cpp
    ./core/src/Rounder.cpp
    ./core/src/SharedArray.cpp
    ./core/src/ToyTree.cpp
    ./core/src/UtilsConfig.cpp
    ./core/src/Utils.cpp)

# -----------------------------------------------------------------------------
# Build core library
# -----------------------------------------------------------------------------

set(ROOT_REQUIRED_LIBS
    ROOT::Core
    ROOT::Core
    ROOT::Graf
    ROOT::Gui
    ROOT::Hist
    ROOT::RIO
    ROOT::RooFit
    ROOT::Tree)

set(CORE_LIB "${PROJECT_NAME}Core")
add_library(${CORE_LIB} SHARED ${CORE_LIB_SOURCES})
target_link_libraries(${CORE_LIB} ${FREETYPE_LIBRARIES} ${ROOT_REQUIRED_LIBS})
# For ConfigVersion.h
target_include_directories(${CORE_LIB} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# -----------------------------------------------------------------------------
# Add ROOT dictionaries to the core library
# -----------------------------------------------------------------------------

# See
# https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html#dictionary-generation
# But now RootNewMacros has been renamed RootMacros and is automatically
# included with find_package(ROOT)
set(CORE_DICTIONARY_SOURCES
    RooBinned2DBicubicBase.h
    RooCrossCorPdf.h
    RooHistPdfAngleVar.h
    RooHistPdfVar.h
    RooPoly3Var.h
    RooPoly4Var.h
    RooSlimFitResult.h
    RooMultiPdf.h)
root_generate_dictionary(G__${CORE_LIB} ${CORE_DICTIONARY_SOURCES} MODULE
                         ${CORE_LIB} LINKDEF core/include/coreLinkDef.h)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${CORE_LIB}_rdict.pcm
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
