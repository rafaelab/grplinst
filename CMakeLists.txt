cmake_minimum_required(VERSION 3.0)

# Change name of the library to something meaningful!
set(PLUGINNAME "grplinst")
project(${PLUGINNAME} CXX)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)

# get information on crpropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/findCRPropa.py" swig_interface OUTPUT_VARIABLE crpropaSwigPath)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/findCRPropa.py" install_prefix OUTPUT_VARIABLE crpropaInstallPrefix)
message(STATUS "CRPropa includes : ${crpropaInstallPrefix}")

# use SWIG's builtin
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON) # Segmentation errors if false
if(ENABLE_SWIG_BUILTIN)
  set(crpropaSwigFile "crpropa-builtin.i")
  set(swigModeFlag "-builtin")
endif(ENABLE_SWIG_BUILTIN)

find_path(CRPROPA_SWIG_INTERFACE_FILES ${crpropaSwigFile} HINTS share/crpropa share/crpropa/swig_interface python/ ${crpropaSwigPath})
find_path(CRPROPA_INCLUDE_PATH CRPropa.h HINTS crpropa/ include/ include/crpropa ${crpropaInstallPrefix}/include/)
find_library(CRPROPA_LIBRARY NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${crpropaInstallPrefix}/lib/)

message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa includes : ${CRPROPA_INCLUDE_PATH}")
message(STATUS "CRPropa lib: ${CRPROPA_LIB}")

# swig plugin
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true )
add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx COMMAND swig  ${swigModeFlag} -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDE_PATH} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGINNAME}.i)
include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# build interface and source in one library
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGINNAME}.cc
)
target_link_libraries(${PLUGINNAME} ${CRPROPA_LIB} ${PYTHON_LIBRARIES})
set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(
sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )


# Uncomment these lines for installation of the plugin
# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
# install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")


