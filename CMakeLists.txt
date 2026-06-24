cmake_minimum_required(VERSION 3.16)

project(IFXC VERSION 0.1.0 LANGUAGES C)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build the test suite" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(BASH_EXECUTABLE bash)

set(IFXC_GENERATED_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/generated/ifxc_mgga_xc_ml25_features.c)
set(IFXC_GENERATED_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/src/generated/ifxc_mgga_xc_ml25_features.h)
set(IFXC_MAPLE_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/maple/ifxc_mgga_xc_ml25_features.mpl)
set(IFXC_GENERATE_ARGS
  ${Python3_EXECUTABLE}
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/maple2c_ifxc.py
  --source ${IFXC_MAPLE_SOURCE}
  --feature-def ${CMAKE_CURRENT_SOURCE_DIR}/src/features/ifxc_ml25_features.def
  --feature-set ml25
  --max-deriv-order 4
  --backend-source ${CMAKE_CURRENT_SOURCE_DIR}/src/ifxc_ml25_backend.c
  --output-c ${IFXC_GENERATED_SOURCE}
  --output-h ${IFXC_GENERATED_HEADER}
)

add_custom_command(
  OUTPUT ${IFXC_GENERATED_SOURCE} ${IFXC_GENERATED_HEADER}
  COMMAND ${IFXC_GENERATE_ARGS}
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/maple2c_ifxc.py
    ${IFXC_MAPLE_SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/features/ifxc_ml25_features.def
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

add_custom_target(generate
  COMMAND ${IFXC_GENERATE_ARGS}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

add_library(ifxc
  src/ifxc_api.c
  src/ifxc_errors.c
  src/ifxc_eval.c
  src/ifxc_metadata.c
  ${IFXC_GENERATED_SOURCE}
  src/libxc_ref/src/mgga_xc_ml25.c
)

add_library(IFXC::ifxc ALIAS ifxc)

target_include_directories(ifxc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/libxc_ref/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/generated
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
  enable_testing()

  add_executable(test_api
    tests/test_api.c
  )
  target_link_libraries(test_api PRIVATE ifxc)

  add_test(NAME api COMMAND test_api)

  if(BASH_EXECUTABLE AND BUILD_SHARED_LIBS)
    add_test(
      NAME public_symbols
      COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_public_symbols.sh $<TARGET_FILE:ifxc>
    )
  endif()

  if(BASH_EXECUTABLE)
    add_test(
      NAME generated_freshness
      COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_generated.sh ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
endif()
