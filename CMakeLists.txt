cmake_minimum_required(VERSION 3.20)
project(RockSolver VERSION 0.1.0 LANGUAGES CXX)

option(ROCKSOLVER_BUILD_TESTS "Build RockSolver unit tests" ON)
option(ROCKSOLVER_BUILD_WASM "Build the Emscripten WebAssembly adapter" OFF)

add_library(rocksolver_core STATIC src/rocksolver.cpp src/smrm_material.cpp src/mohr_coulomb.cpp src/seepage.cpp)
target_include_directories(rocksolver_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(rocksolver_core PUBLIC cxx_std_20)
set_target_properties(rocksolver_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(MSVC)
  target_compile_options(rocksolver_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(rocksolver_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(rocksolver_native app/main.cpp)
target_link_libraries(rocksolver_native PRIVATE rocksolver_core)

if(ROCKSOLVER_BUILD_TESTS)
  enable_testing()
  add_executable(rocksolver_tests tests/test_main.cpp)
  target_link_libraries(rocksolver_tests PRIVATE rocksolver_core)
  add_executable(rocksolver_stage5_smrm_biot_tests tests/test_stage5_smrm_biot.cpp)
  target_link_libraries(rocksolver_stage5_smrm_biot_tests PRIVATE rocksolver_core)
  if(MSVC)
    target_compile_options(rocksolver_stage5_smrm_biot_tests PRIVATE /UNDEBUG)
  else()
    target_compile_options(rocksolver_stage5_smrm_biot_tests PRIVATE -UNDEBUG)
  endif()
  add_test(NAME rocksolver_tests COMMAND rocksolver_tests)
  add_test(NAME rocksolver_stage5_smrm_biot_tests COMMAND rocksolver_stage5_smrm_biot_tests)
endif()

if(ROCKSOLVER_BUILD_WASM)
  if(NOT EMSCRIPTEN)
    message(FATAL_ERROR "ROCKSOLVER_BUILD_WASM requires Emscripten/emcmake")
  endif()
  # Surface model validation errors as JSON through the existing try/catch adapter.
  target_compile_options(rocksolver_core PRIVATE -fexceptions)
  add_executable(rocksolver_wasm wasm/rocksolver_wasm.cpp)
  target_compile_options(rocksolver_wasm PRIVATE -fexceptions)
  target_link_libraries(rocksolver_wasm PRIVATE rocksolver_core)
  target_link_options(rocksolver_wasm PRIVATE
    "-sDISABLE_EXCEPTION_CATCHING=0"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sENVIRONMENT=web,worker,node"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMAXIMUM_MEMORY=2147483648"
    "-sEXPORTED_FUNCTIONS=['_malloc','_free','_rocksolver_version','_rocksolver_run_preset','_rocksolver_run_config','_rocksolver_run_imported','_rocksolver_run_steady_seepage','_rocksolver_run_baseline','_rocksolver_run_initial_state','_rocksolver_resume','_rocksolver_water_stage','_rocksolver_reduction_trial','_rocksolver_result_ptr','_rocksolver_result_size','_rocksolver_checkpoint_ptr','_rocksolver_checkpoint_size']"
    "-sEXPORTED_RUNTIME_METHODS=['UTF8ToString','stringToUTF8','lengthBytesUTF8']"
  )
  set_target_properties(rocksolver_wasm PROPERTIES SUFFIX ".js")
endif()
