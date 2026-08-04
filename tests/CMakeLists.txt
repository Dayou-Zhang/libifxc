# SPDX-License-Identifier: MPL-2.0

find_program(IFXC_BASH_EXECUTABLE bash)
find_package(Python3 QUIET COMPONENTS Interpreter)

add_executable(test_api test_api.c)
target_link_libraries(test_api PRIVATE ifxc)
add_test(NAME api COMMAND test_api)

if(IFXC_BASH_EXECUTABLE AND BUILD_SHARED_LIBS)
  add_test(
    NAME public_symbols
    COMMAND ${IFXC_BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_public_symbols.sh $<TARGET_FILE:ifxc>
  )
endif()

if(Python3_Interpreter_FOUND)
  add_custom_target(check-generated
    COMMAND
      ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/scripts/check_ml25_manifest.py
      --source ${PROJECT_SOURCE_DIR}/maple/if_mgga/mgga_xc_ml25.mpl
      --feature-def ${PROJECT_SOURCE_DIR}/src/features/ifxc_ml25_features.def
    COMMAND
      ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/scripts/check_ml26_manifest.py
      --source ${PROJECT_SOURCE_DIR}/maple/if_mgga/mgga_xc_ml26.mpl
      --feature-def ${PROJECT_SOURCE_DIR}/src/features/ifxc_ml26_features.def
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    VERBATIM
  )

  if(IFXC_BASH_EXECUTABLE)
    add_test(
      NAME generated_manifests
      COMMAND
        ${IFXC_BASH_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/scripts/check_generated.sh
        ${Python3_EXECUTABLE}
    )
  endif()
else()
  message(STATUS
    "Python 3 not found: check-generated and generated_manifests are unavailable"
  )
endif()
