# The combiner name should the the name of the subdirectory of the combiner.
set(COMBINER_NAME tutorial)

# Define all executables of this combiner. They have to be located in the `main`
# subdirectory, and end with `.cpp`.
set(COMBINER_EXECUTABLES
    ${COMBINER_NAME} cartesian tutorial_dataset
    tutorial_dataset_build_workspace tutorial_dataset_multipdf
    tutorial_dataset_multipdf_build_workspace)

# Configure if this combiner has custom ROOT objects
set(HAS_CUSTOMROOTOBJECTS FALSE)

# Define the custom ROOT objects
if(HAS_CUSTOMROOTOBJECTS)
  set(COMBINER_DICTIONARY_LINKDEF
      LinkDef.h # <-- change to name of LinkDef file
  )
  set(COMBINER_DICTIONARY_SOURCES Header.h # <-- add header files here
  )
endif()

# -----------------------------------------------------------------------------
# Set directory names and includes
# -----------------------------------------------------------------------------

set(COMBINER_INCLUDE_DIR ./include)
set(COMBINER_SOURCE_DIR ./src)
set(COMBINER_MAIN_DIR ./main)

# include directories
include_directories(BEFORE ${ROOT_INCLUDE_DIR} ${CORE_INC_DIR}
                    ${COMBINER_INCLUDE_DIR})

# -----------------------------------------------------------------------------
# Generate ROOT dictionaries
# -----------------------------------------------------------------------------

if(HAS_CUSTOMROOTOBJECTS)
  # Add a make rule for the dictionary
  root_generate_dictionary(${COMBINER_NAME}Dict ${COMBINER_DICTIONARY_SOURCES}
                           LINKDEF ${COMBINER_DICTIONARY_LINKDEF})
  # Add a target calling the make rule for the dictionary
  add_custom_target(${COMBINER_NAME}Dict DEPENDS ${COMBINER_NAME}Dict.cxx
                                                 ${COMBINER_NAME}Dict.h)
endif()

# -----------------------------------------------------------------------------
# Build the library
# -----------------------------------------------------------------------------

aux_source_directory(${COMBINER_SOURCE_DIR} COMBINER_LIB_SOURCES)
if(HAS_CUSTOMROOTOBJECTS)
  set(COMBINER_LIB_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/${COMBINER_NAME}Dict.cxx
                           ${COMBINER_LIB_SOURCES})
endif()
add_shared_library(${COMBINER_NAME}Components ${COMBINER_LIB_SOURCES})
if(HAS_CUSTOMROOTOBJECTS)
  add_dependencies(${COMBINER_NAME}Components ${COMBINER_NAME}Dict)
endif()
target_link_libraries(${COMBINER_NAME}Components ${ROOT_LIBRARIES}
                      ${PROJECT_NAME}Components)

# -----------------------------------------------------------------------------
# Build the executables
# -----------------------------------------------------------------------------

set(COMBINER_LIBS ${COMBINER_NAME}Components ${PROJECT_NAME}Components)
foreach(exec ${COMBINER_EXECUTABLES})
  add_executable(${exec} ${COMBINER_MAIN_DIR}/${exec}.cpp)
  target_link_libraries(${exec} ${COMBINER_LIBS})
endforeach()

# -----------------------------------------------------------------------------
# Install the binaries from the build directory back into the project
# subdirectory
# -----------------------------------------------------------------------------

foreach(exec ${COMBINER_EXECUTABLES})
  install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/${exec}
          DESTINATION ${COMBINER_NAME}/bin/)
endforeach()

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

install(DIRECTORY DESTINATION ${COMBINER_NAME}/root)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/C)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/cl)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/dot)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/eps)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/par)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/pdf)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/png)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/root)
install(DIRECTORY DESTINATION ${COMBINER_NAME}/plots/scanner)
