cmake_minimum_required(VERSION 3.14)
project(grplinst CXX)

message(STATUS "Building grplinst...")


# set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# compilation flags
set(CMAKE_CXX_FLAGS "-std=c++11")

# enable extra warnings on debug builds
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wextra")


# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")


# find CRPropa
include(cmake/FindCRPropa.cmake)


# assemble includes (using own Eigen version)
include_directories(include 
	${CRPropa_kiss_INCLUDE_DIR}
	${CRPropa_HepID_INCLUDE_DIR}
	${CRPropa_INCLUDE_DIR}  
	${CMAKE_CURRENT_BINARY_DIR}/include 
	${CMAKE_CURRENT_SOURCE_DIR}/../include
	)


# source files for compilation
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Flow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Medium.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaInstability.cc
	)
file(GLOB GRPLINST_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/grplinst/*.h")


# decide on whether to use Python or not
option(ENABLE_PYTHON "Create python library via SWIG" On)
message(STATUS "Build Python bindings via SWIG: ${ENABLE_PYTHON}")


# compile with Python support
if(ENABLE_PYTHON)
	# find Python
	find_package(Python COMPONENTS Interpreter Development NumPy)

	if(Python_FOUND)
		message(STATUS "Python version: ${Python_VERSION}")
		message(STATUS "Python executable: ${Python_EXECUTABLE}")
		message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python library: ${Python_LIBRARY}")
		message(STATUS "Python site-packges: ${Python_SITE_PACKAGES}")
	else(Python_FOUND)
		set(Python_HELPER_EXECUTABLE "Python_EXECUTABLE (path to python executable)")
		set(Python_HELPER_INCLUDE_DIRS "Python_INCLUDE_DIRS (directory where Python.h is located)")
		set(Python_HELPER_LIBRARY "Python_LIBRARY (path to file libpython.so or libpython.dylib [in OSX]; sometimes contains version before extension)")
		message(STATUS 	"Python: NOT FOUND. Try setting the relevant flags manually: \n\t${Python_HELPER_EXECUTABLE} \n\t${Python_HELPER_LIBRARY} \n\t${Python_HELPER_INCLUDE_DIRS}")
	endif(Python_FOUND)


	if(Python_NumPy_FOUND)
		set(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
		include_directories(${Python_NumPy_INCLUDE_DIRS})
		message(STATUS "NumPy headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
	else()
		message(STATUS "NumPy: NOT FOUND. The code might work just fine without NumPy, but some features in Python might not be straightfoward to use.")
	endif(Python_NumPy_FOUND)

	# find SWIG
	find_package(SWIG REQUIRED)


	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif()

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
	else(SWIG_FOUND)
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)

	# # add Python headers
	# include_directories(include ${Python_INCLUDE_DIRS})


	# SWIG plugin
	set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/grplinst_wrap.cxx PROPERTIES GENERATED True)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/grplinst_wrap.cxx
		COMMAND swig ${SWIG_MODE_FLAG} -c++ -python ${PY3} 
		-I${CRPropa_SWIG_INTERFACE_FILE}
		-I${CMAKE_SOURCE_DIR}/include 
		-I${Python_INCLUDE_DIRS}
		-I{Eigen_INCLUDE_PATH}
		-I${CRPropa_SWIG_PATH}  
		-I${CRPropa_INCLUDE_DIR} ${CRPropa_SWIG_DEFINES} 
		-dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/grplinst_wrap.cxx 
		-outdir ${CMAKE_CURRENT_BINARY_DIR} 
		${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst.i
		DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst.i ${GRPLINST_HEADERS}
	)
	include_directories(include 
		${CMAKE_CURRENT_BINARY_DIR}/include 
		${CRPropa_INCLUDE_PATH} 
		${Python_INCLUDE_DIRS} 
		${CMAKE_CURRENT_SOURCE_DIR}/../include
	)
	
	# build interface and source in one library
	add_library(grplinst MODULE ${CMAKE_CURRENT_BINARY_DIR}/grplinst_wrap.cxx ${SOURCE_FILES})
	target_link_libraries(grplinst ${CRPropa_LIBRARY} ${Python_LIBRARIES})

# compile code to work in C++ only, without Python support
else()

	# build interface and source in one library
	add_library(grplinst MODULE ${SOURCE_FILES})
	target_link_libraries(grplinst ${CRPropa_LIBRARY})
endif()

set_target_properties(grplinst PROPERTIES PREFIX "")
set_target_properties(grplinst PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(grplinst PROPERTIES OUTPUT_NAME "_grplinst")


# option(ENABLE_TESTING "Build tests and enable test target" On)
# if(ENABLE_TESTING)
# 	set(GTEST_INCLUDE_PATH ${CRPropa_INSTALL_PREFIX}/../libs/gtest/include)
# 	include_directories(${GTEST_INCLUDE_PATH})
# 	link_directories(${CRPropa_INSTALL_PREFIX}/libs/gtest)
# 	message(STATUS "Google testing suite (from CRPropa): ${CRPropa_INSTALL_PREFIX}/libs/gtest/include")
# 	enable_testing()
# endif(ENABLE_TESTING)

# set(grplinst_DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data)


# Install
add_definitions(-Dgrplinst_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS grplinst DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst.i DESTINATION swig_interface)
