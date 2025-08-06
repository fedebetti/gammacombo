# collection of helpful macros for cmake created by J. Engels, DESY and Ch.
# Rosemann, DESY

# create symbolic lib target for calling library targets
macro(ADD_SHARED_LIBRARY _name)
  add_library(${_name} SHARED ${ARGN})

  # change lib_target properties
  set_target_properties(
    ${_name}
    PROPERTIES # create *nix style library versions + symbolic links
               VERSION ${${PROJECT_NAME}_VERSION} SOVERSION
                                                  ${${PROJECT_NAME}_SOVERSION})
endmacro(ADD_SHARED_LIBRARY)

# in order to include cmake projects into other projects config files must exist
# helper macro for generating project configuration file
macro(GENERATE_PACKAGE_CONFIGURATION_FILES)

  foreach(arg ${ARGN})
    if(${arg} MATCHES "Config.cmake")
      if(EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in")
        configure_file("${PROJECT_SOURCE_DIR}/cmake/${arg}.in"
                       "${PROJECT_BINARY_DIR}/${arg}" @ONLY)
        install(FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION .)
      endif()
    endif()

    if(${arg} MATCHES "ConfigVersion.cmake")
      # version configuration file
      if(EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in")
        configure_file("${PROJECT_SOURCE_DIR}/cmake/${arg}.in"
                       "${PROJECT_BINARY_DIR}/${arg}" @ONLY)
        install(FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION .)
      endif(EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in")
    endif()

    if(${arg} MATCHES "LibDeps.cmake")
      export_library_dependencies("${arg}")
      install(FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION lib/cmake)
    endif()

  endforeach()

endmacro(GENERATE_PACKAGE_CONFIGURATION_FILES)
