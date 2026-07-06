cmake_minimum_required(VERSION 3.14)


# set the version to whatever you want
project(grplinst 
	VERSION 2.0.0 
	LANGUAGES Fortran C CXX
)


###############################################################################
#                   Compiler- and system-specific options                     #
###############################################################################

# C++ should be at least 17
set(CMAKE_CXX_STANDARD 17)

# set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

# set other needed compiler flags:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")

# enable extra warnings on debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wextra")

# find cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# compiler-specific flags for Linux and macOS
if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
	message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# RPATH handling (Linux and macOS)
set(CMAKE_MACOSX_RPATH On)
set(CMAKE_SKIP_BUILD_RPATH False)
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)
if(APPLE)
	set(originToken "@loader_path")
else()
	set(originToken "$ORIGIN")
endif(APPLE)

# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")

# search next to the binary and in the plugin's own install lib directory
set(CMAKE_INSTALL_RPATH "${originToken};${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "Install RPATH: ${CMAKE_INSTALL_RPATH} (+ linked library dirs)")



# general options
option(ENABLE_PYTHON "Create Python bindings with SWIG" On)
option(USE_OWN_CRPROPA "Use builtin CRPropa instead of in CRPROPA_PATH installed CRPropa" Off)
option(ENABLE_TESTING "Build tests and enable test target" On)




###############################################################################
#                   Compilation and installation of plugin                    #
###############################################################################

# finding Python
if(ENABLE_PYTHON)
	set(Python_FIND_VIRTUALENV FIRST)
	set(Python_FIND_STRATEGY LOCATION)
	find_package(Python 3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
	
	if(Python_FOUND)
		message(STATUS "Python version: ${Python_VERSION}")
		message(STATUS "Python executable: ${Python_EXECUTABLE}")
		message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python library: ${Python_LIBRARY}")
		message(STATUS "Python site-packages: ${Python_SITE_PACKAGES}")
	endif(Python_FOUND)

	if(Python_NumPy_FOUND)
		set(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND CRPropa_SWIG_DEFINES -DWITHNUMPY)
		include_directories(${Python_NumPy_INCLUDE_DIRS})
		message(STATUS "NumPy headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
	else()
		message(STATUS "NumPy: NOT FOUND.")
		message(STATUS "The code might work just fine without NumPy, but some features in Python might not be straightfoward to use.")
	endif(Python_NumPy_FOUND)

	#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
	set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
	message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")


	# find SWIG
	find_package(SWIG REQUIRED)

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif(SWIG_VERSION VERSION_LESS 4.0)	

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
	else()
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)

endif(ENABLE_PYTHON)

set(grplinst_EXTRA_INCLUDES)
set(grplinst_EXTRA_LIBRARIES)
set(grplinst_EXTRA_SOURCES)

# source files for compilation
set(SOURCE_FILES
	${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Flow.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Medium.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaInstability.cc
)
file(GLOB grplinst_HEADERS 
	"${CMAKE_CURRENT_SOURCE_DIR}/include/grplinst/*.h"
)


# build crpropa and add resulting libcrpropa.so as library
if(USE_OWN_CRPROPA)
	if(NOT CRPropa_POPULATED)
		include(FetchContent)
		FetchContent_Declare(crpropa GIT_REPOSITORY "https://github.com/CRPropa/CRPropa3.git")
		list(APPEND CMAKE_MODULE_PATH "${CRPropa_SOURCE_DIR}/cmake")
		FetchContent_MakeAvailable(CRPropa)
	endif(NOT CRPropa_POPULATED)

	list(APPEND grplinst_EXTRA_INCLUDES ${CRPropa_SOURCE_DIR}/include)
	list(APPEND grplinst_EXTRA_INCLUDES ${CRPROPA_EXTRA_INCLUDES})
	list(APPEND grplinst_EXTRA_LIBRARIES crpropa)
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SOURCE_DIR}/python)
else()
	find_package(CRPropa REQUIRED)
	list(APPEND grplinst_EXTRA_INCLUDES ${CRPropa_INCLUDE_DIR})
	list(APPEND grplinst_EXTRA_LIBRARIES ${CRPropa_LIBRARY})
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SWIG_PATH})
endif(USE_OWN_CRPROPA)


# source files
add_library(grplinst SHARED
	${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Flow.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Medium.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaInstability.cc
	${grplinst_EXTRA_SOURCES}
)

# headers
target_include_directories(grplinst PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${grplinst_EXTRA_INCLUDES}
)

# link crpropa and other needed libraries
target_link_libraries(grplinst ${grplinst_EXTRA_LIBRARIES})

# Python interface, this part should work basically fully automated.
if(ENABLE_PYTHON)
	if(Python_FOUND)
		find_package(SWIG 3.0 REQUIRED)
		include(UseSWIG) 
		add_definitions(-Dgrplinst_HAVE_PYTHON)
		list(APPEND grplinst_SWIG_DEFINES -Dgrplinst_HAVE_PYTHON)

		# print Python info in detail
		if(NOT USE_OWN_CRPROPA)
			message(STATUS "Python: Found!")
			message(STATUS "  version ${Python_VERSION}")
			message(STATUS "  executable: ${Python_EXECUTABLE}")
			message(STATUS "  libraries: ${Python_LIBRARIES}")
			message(STATUS "  headers: ${Python_INCLUDE_DIRS}")
			message(STATUS "  site packages: ${Python_SITELIB}")
			if(Python_Development_FOUND)
				message(STATUS "  development libraries: Found!")
			else()
				message(STATUS "  development libraries: NOT found!")
			endif(Python_Development_FOUND)

			# numpy is required, otherwise Python_FOUND = false
			message(STATUS "NumPy: Found!")
			message(STATUS "  headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
		endif(NOT USE_OWN_CRPROPA)

		#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
		set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
		message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")

		set(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND grplinst_SWIG_DEFINES -DWITHNUMPY)
		list(APPEND grplinst_SWIG_DEFINES -doxygen)

		# use SWIG's builtin only when explicitly requested; the non-builtin mode is more robust with CRPropa
		option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" Off)

		if(ENABLE_SWIG_BUILTIN)
			set(BUILTIN "-builtin")
			list(APPEND grplinst_SWIG_DEFINES "-builtin")
		else()
			set(BUILTIN "")
		endif(ENABLE_SWIG_BUILTIN)

		file(GLOB_RECURSE grplinst_SWIG_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/python/*.i)
		set(grplinst_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/grplinst")
		file(MAKE_DIRECTORY "${grplinst_PYTHON_PACKAGE_DIR}")
		file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst/" DESTINATION "${grplinst_PYTHON_PACKAGE_DIR}")

		set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst.i
			PROPERTIES
			CPLUSPLUS ON
			SWIG_MODULE_NAME "grplinst"
			SWIG_FLAGS "${grplinst_SWIG_DEFINES};-dirprot;${CMAKE_SWIG_FLAGS}"
			OUTPUT_DIR "${grplinst_PYTHON_PACKAGE_DIR}"
			INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include;${CRPropa_SWIG_INTERFACE_FILES};${grplinst_EXTRA_INCLUDES};"
		)
		
		swig_add_library(grplinst-swig TYPE SHARED LANGUAGE python
			SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst.i
		)
		
		target_include_directories(grplinst-swig PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${grplinst_EXTRA_INCLUDES}
			${Python_INCLUDE_DIRS}
			${Python_NumPy_INCLUDE_DIRS}
			${CRPropa_SWIG_INTERFACE_FILES}
		)
		
		target_link_libraries(grplinst-swig grplinst ${Python_LIBRARIES} ${Python_LIBRARY} Python::Module)
		set_target_properties(grplinst-swig PROPERTIES
			PREFIX ""
			OUTPUT_NAME "_grplinst"
			LIBRARY_OUTPUT_DIRECTORY "${grplinst_PYTHON_PACKAGE_DIR}"
			COMPILE_FLAGS  "-w"
			SUFFIX ".so"
		)

		if(APPLE)
			set_target_properties(grplinst-swig PROPERTIES INSTALL_RPATH "@loader_path;${CMAKE_INSTALL_PREFIX}/lib")
			set_property(TARGET grplinst-swig APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
		endif(APPLE)

	endif(Python_FOUND)

endif(ENABLE_PYTHON)


###############################################################################
#                                 Testing                                     #
###############################################################################

# the testing part, where tests for your code should be added
if(ENABLE_TESTING)
	enable_testing()

	# temporary workaround for newer clang versions due to its handling of unwinding
	#   see: https://github.com/google/googletest/issues/3062
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
		if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
			add_definitions(-DCRPROPA_TESTS_SKIP_EXCEPTIONS)
		endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

	# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
	#     Why is it not recommended use a pre-compiled copy of Google Test?)
	# Here are only the apple fixes for gtest
	if(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
	endif(APPLE)

	# gtest ships with CRPropa, so its headers and libraries follow CRPropa's paths:
	if(USE_OWN_CRPROPA)
		# gtest is built together with CRPropa (see FetchContent above); use its targets directly
		set(GTEST_INCLUDE_DIR ${CRPropa_SOURCE_DIR}/libs/gtest/include)
		set(GTEST_LIBRARIES gtest gtest_main)
	else()
		# headers sit under the CRPropa include dir (gtest/), libraries next to (install) or under lib/ (build tree) libcrpropa
		set(GTEST_INCLUDE_DIR ${CRPropa_INCLUDE_DIR})
		find_library(GTEST_LIBRARY NAMES gtest HINTS ${CRPropa_LIBRARY_DIR} ${CRPropa_LIBRARY_DIR}/lib)
		find_library(GTEST_MAIN_LIBRARY NAMES gtest_main HINTS ${CRPropa_LIBRARY_DIR} ${CRPropa_LIBRARY_DIR}/lib)
		set(GTEST_LIBRARIES ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
	endif(USE_OWN_CRPROPA)

	# C++ test executables — one per test suite
	foreach(TEST_NAME testFlow testMedium testPlasmaInstability)
		add_executable(${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_NAME}.cpp)
		target_include_directories(${TEST_NAME} PRIVATE
			${GTEST_INCLUDE_DIR}
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${grplinst_EXTRA_INCLUDES}
		)
		target_link_libraries(${TEST_NAME} grplinst ${GTEST_LIBRARIES} pthread)
		add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
	endforeach()

	if(ENABLE_PYTHON AND Python_FOUND)
		file(
			COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/
			DESTINATION  ${CMAKE_CURRENT_BINARY_DIR}
			FILES_MATCHING PATTERN "*.py"
		)

		if(USE_OWN_CRPROPA)
			set(CRPropa_PYTHON_TEST_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/crpropa")
			file(MAKE_DIRECTORY "${CRPropa_PYTHON_TEST_PACKAGE_DIR}")
			file(WRITE "${CRPropa_PYTHON_TEST_PACKAGE_DIR}/__init__.py" "from .crpropa import *\n")
			add_test(
				NAME moveCRPropaPlugin
				COMMAND ${CMAKE_COMMAND} -E copy
					"${CRPropa_BINARY_DIR}/crpropa.py"
					"$<TARGET_FILE:crpropa-swig>"
					"$<TARGET_FILE:crpropa>"
					"${CRPropa_PYTHON_TEST_PACKAGE_DIR}"
			)
		endif(USE_OWN_CRPROPA)

		# add python tests here, just copy paste the add_test statements and change test name and path
		add_test(NAME testPluginPython COMMAND ${Python_EXECUTABLE} testPlugin.py)
		if(USE_OWN_CRPROPA)
			set_tests_properties(testPluginPython PROPERTIES DEPENDS moveCRPropaPlugin)
		endif(USE_OWN_CRPROPA)
	endif(ENABLE_PYTHON AND Python_FOUND)

endif(ENABLE_TESTING)



###############################################################################
#                                Installation                                 #
###############################################################################

# specify installation instructions
install(TARGETS grplinst DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_PREFIX}/include" FILES_MATCHING PATTERN "*.h")
if(ENABLE_PYTHON AND Python_FOUND)
	install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/grplinst" DESTINATION "${Python_INSTALL_PACKAGE_DIR}")
	install(FILES "${grplinst_PYTHON_PACKAGE_DIR}/grplinst.py" DESTINATION "${Python_INSTALL_PACKAGE_DIR}/grplinst")
	install(TARGETS grplinst-swig LIBRARY DESTINATION "${Python_INSTALL_PACKAGE_DIR}/grplinst")
	install(FILES ${grplinst_SWIG_INPUTS} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/grplinst/swig_interface)
endif(ENABLE_PYTHON AND Python_FOUND)
