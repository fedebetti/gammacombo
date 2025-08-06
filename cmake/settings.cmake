# general settings for a cmake project less clutter in the real CMakeLists.txt
# file

# library *nix style versioning
set(${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}")
set(${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_VERSION_PATCH}")

# define output directories
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")
set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib")
mark_as_advanced(EXECUTABLE_OUTPUT_PATH)
mark_as_advanced(LIBRARY_OUTPUT_PATH)

# what happens when `make install` is invoked: set default install prefix to
# project root directory instead of the cmake default /usr/local
if(CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
endif()

# write this variable to cache
set(CMAKE_INSTALL_PREFIX
    "${CMAKE_INSTALL_PREFIX}"
    CACHE PATH "Where to install ${PROJECT_NAME}" FORCE)
