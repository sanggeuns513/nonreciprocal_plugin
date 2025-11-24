cmake_minimum_required(VERSION 3.15...3.27 FATAL_ERROR)
project(nonreciprocal_plugin LANGUAGES C CXX)

find_package(HOOMD REQUIRED)

# Version Check
if (DEFINED HOOMD_VERSION)
    message(STATUS "Found HOOMD ${HOOMD_VERSION}: ${HOOMD_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}")

    if (HOOMD_VERSION VERSION_LESS "4.0.0")
        message(FATAL_ERROR "HOOMD >= 4.0.0 is required.")
    elseif (HOOMD_VERSION VERSION_GREATER "6.0.0")
        message(WARNING "This plugin has not been tested with HOOMD-blue version 6 or higher, and may not function as expected with newer versions.")
    endif()
endif()

# Force installation to the HOOMD installation location.
set(CMAKE_INSTALL_PREFIX ${HOOMD_INSTALL_PREFIX} CACHE PATH "Installation prefix" FORCE)

# Enable compiler warnings on gcc and clang.
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wno-sign-conversion -Wno-unknown-pragmas -Wno-deprecated-declarations -Wno-unused-result")
endif()

add_subdirectory(src)