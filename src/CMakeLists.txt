set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "default")

# Force a default install config (multi-config generators)
if(NOT CMAKE_INSTALL_CONFIG_NAME)
    set(CMAKE_INSTALL_CONFIG_NAME Release)
endif()

set(COMPONENT_NAME nonreciprocal_plugin)

set(SOURCES module-nonreciprocal_plugin.cc nonreciprocal_plugin.cc)

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

add_library(_${COMPONENT_NAME} SHARED ${SOURCES} ${_cuda_sources})
set_target_properties(_${COMPONENT_NAME} PROPERTIES PREFIX "" SUFFIX ".so" CXX_STANDARD 17)
if (APPLE)
    set_target_properties(_${COMPONENT_NAME} PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()
target_link_libraries(_${COMPONENT_NAME} PUBLIC HOOMD::_hoomd PUBLIC HOOMD::_md)
target_include_directories(_${COMPONENT_NAME} PRIVATE ${HOOMD_INCLUDE_DIRS})

# Install the library
install(TARGETS _${COMPONENT_NAME} LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${COMPONENT_NAME})
install(FILES ${python_files} DESTINATION ${PYTHON_SITE_INSTALL_DIR}/${COMPONENT_NAME})