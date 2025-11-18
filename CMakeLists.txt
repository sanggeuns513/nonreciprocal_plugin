cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
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

set(SOURCES nonreciprocal_plugin.cc)

if (ENABLE_HIP)
list(APPEND SOURCES nonreciprocal_pluginGPU.cc)
endif()

set(CU_SOURCES nonreciprocal_pluginGPU.cu)
if (ENABLE_HIP)
set(_cuda_sources ${CU_SOURCES})
set_source_files_properties(${CU_SOURCES} PROPERTIES LANUAGE ${HOOMD_DEVICE_LANGUAGE})
endif (ENABLE_HIP)

set(python_files
    __init__.py
    force.py
    version.py
)

add_library(_${PROJECT_NAME} SHARED ${SOURCES} ${_cuda_sources})
set_target_properties(_${PROJECT_NAME} PROPERTIES PREFIX "" SUFFIX ".so" CXX_STANDARD 17)
if (APPLE)
    set_target_properties(_${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()
target_link_libraries(_${PROJECT_NAME} PUBLIC HOOMD::_hoomd PUBLIC HOOMD::_md)
target_include_directories(_${PROJECT_NAME} PRIVATE ${HOOMD_INCLUDE_DIRS})

# Install the library
install(TARGETS _${PROJECT_NAME} LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${PROJECT_NAME})
install(FILES ${python_files} DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${PROJECT_NAME})