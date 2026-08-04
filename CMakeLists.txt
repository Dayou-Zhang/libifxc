# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.16)

project(IFXC VERSION 0.4.0 LANGUAGES C)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

add_library(ifxc
  src/ifxc_api.c
  src/ifxc_derivatives.c
  src/ifxc_errors.c
  src/ifxc_eval.c
  src/ifxc_metadata.c
  src/ifxc_ml25_backend.c
  src/ifxc_ml25_combined.c
  src/ifxc_ml26_backend.c
  src/ifxc_ml26_combined.c
)

add_library(IFXC::ifxc ALIAS ifxc)

target_include_directories(ifxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_features(ifxc PUBLIC c_std_99)

set_target_properties(ifxc PROPERTIES
  OUTPUT_NAME ifxc
  VERSION ${PROJECT_VERSION}
  SOVERSION 0
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN YES
  POSITION_INDEPENDENT_CODE YES
)

if(WIN32 AND BUILD_SHARED_LIBS)
  target_compile_definitions(ifxc PRIVATE IFXC_BUILD_DLL)
endif()

if(UNIX AND NOT APPLE)
  target_link_libraries(ifxc PRIVATE m)
endif()

install(TARGETS ifxc
  EXPORT IFXCTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES include/ifxc.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES LICENSE NOTICE
  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/libifxc
)

install(EXPORT IFXCTargets
  NAMESPACE IFXC::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IFXC
)

configure_package_config_file(
  cmake/IFXCConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/IFXCConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IFXC
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/IFXCConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/IFXCConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/IFXCConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IFXC
)

configure_file(
  cmake/libifxc.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/libifxc.pc
  @ONLY
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libifxc.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
