cmake_minimum_required(VERSION 2.9)

set(PLUGINNAME "grplinst")
project(${PLUGINNAME} CXX)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)


# general options
set(CMAKE_CXX_FLAGS "-std=c++11")
set(crpropaSwigFile "crpropa.i")
set(swigModeFlag "")

# use SWIG's builtin
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON) # Segmentation errors if false
if(ENABLE_SWIG_BUILTIN)
  set(crpropaSwigFile "crpropa-builtin.i")
  set(swigModeFlag "-builtin")
endif(ENABLE_SWIG_BUILTIN)

# NumPy
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
  if(numpyIncludePath)
    MESSAGE(STATUS "Found NumPy headers in: " ${numpyIncludePath})
    SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
    list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
    include_directories(${numpyIncludePath})
  else(numpyIncludePath)
    MESSAGE(STATUS "NumPy not found.")
  endif(numpyIncludePath)

# get information on CRPropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE crpropaSwigPath)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE crpropaInstallPrefix)

find_path(CRPROPA_SWIG_INTERFACE_FILES ${crpropaSwigFile} HINTS share/crpropa share/crpropa/swig_interface python/ ${crpropaSwigPath})
find_path(CRPROPA_INCLUDE_PATH CRPropa.h HINTS crpropa/ include/ include/crpropa ${crpropaInstallPrefix}/include/)
find_library(CRPROPA_LIBRARY NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${crpropaInstallPrefix}/lib/)

message(STATUS "CRPropa install prefix: ${crpropaInstallPrefix}")
message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa include path : ${CRPROPA_INCLUDE_PATH}")
message(STATUS "CRPropa library: ${CRPROPA_LIBRARY}")
# message(STATUS "Python executable: ${PYTHON_EXECUTABLE}")
# message(STATUS "Python include path: ${PYTHON_INCLUDE_DIRS}")
# message(STATUS "Python library: ${PYTHON_LIBRARIES}")


# SWIG plugin
file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  COMMAND swig ${swigModeFlag} -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDE_PATH} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/${PLUGINNAME}.i
)
include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../include)


# build interface and source in one library
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Flow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Medium.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaInstability.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Simulation.cc
  )
target_link_libraries(${PLUGINNAME} ${CRPROPA_LIBRARY} ${PYTHON_LIBRARIES})
set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"  
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


# # ----------------------------------------------------------------------------
# # Tests
# # ----------------------------------------------------------------------------
# option(ENABLE_TESTING "Build tests and enable test target" ON)
# if(ENABLE_TESTING)
#   set(GTEST_INCLUDE_PATH ${CRPROPA_INSTALL_PREFIX}/../libs/gtest/include)
#   include_directories(${GTEST_INCLUDE_PATH})
#   link_directories(${CRPROPA_INSTALL_PREFIX}/libs/gtest)
#   message(STATUS "Google testing suite (from CRPropa): ${CRPROPA_INSTALL_PREFIX}/libs/gtest/include")
#   if(APPLE)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
#   endif(APPLE)

#   enable_testing()

#   # add_executable(testVector test/testVector.cpp)
#   # target_link_libraries(testVector ${CRPROPA_LIB} gtest gtest_main pthread)
#   # add_test(testVector testVector)

# endif(ENABLE_TESTING)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
install(TARGETS grplinst DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

