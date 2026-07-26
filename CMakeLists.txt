cmake_minimum_required(VERSION 3.14)


# set the version to whatever you want
project(grplinst 
	VERSION 2.0.0 
	LANGUAGES Fortran C CXX
)
cmake_minimum_required(VERSION 3.14)



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
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wextra")

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

# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")

# general options
option(ENABLE_PYTHON "Create Python bindings with SWIG" On)
option(USE_OWN_CRPROPA "Use builtin CRPropa instead of in CRPROPA_PATH installed CRPropa" Off)
option(ENABLE_TESTING "Build tests and enable test target" On)
option(ENABLE_COVERAGE "Add code for coverage analysis." Off)
option(ENABLE_PYBIND11_STUBGEN "Generate Pybind11 stubs" On)

if(APPLE)
	set(originToken "@loader_path")
else()
	set(originToken "$ORIGIN")
endif(APPLE)

set(CMAKE_INSTALL_RPATH "${originToken};${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "Install RPATH: ${CMAKE_INSTALL_RPATH} (+ linked library dirs)")


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
		add_definitions(-Dgrplinst_HAVE_PYTHON)
	endif(Python_FOUND)

	# find SWIG
	find_package(SWIG REQUIRED)
	include(UseSWIG) 

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif(SWIG_VERSION VERSION_LESS 4.0)	

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
		list(APPEND grplinst_SWIG_DEFINES -Dgrplinst_HAVE_PYTHON)
	else()
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)
	
	# find NumPy
	if(Python_NumPy_FOUND)
		include_directories(${Python_NumPy_INCLUDE_DIRS})
		message(STATUS "NumPy headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
	else()
		message(STATUS "NumPy: NOT FOUND.")
		message(STATUS "The code might work just fine without NumPy, but some features in Python might not be straightfoward to use.")
	endif(Python_NumPy_FOUND)

	#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
	set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
	message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")

endif(ENABLE_PYTHON)

# support coverage
if(ENABLE_COVERAGE)
	if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
		set(COVERAGE_LIBS "-coverage -lgcov")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
	endif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
endif(ENABLE_COVERAGE)


# build CRPropa and add resulting libcrpropa.so as library
if(USE_OWN_CRPROPA)
	if(NOT crpropa_POPULATED)
		include(FetchContent)
		FetchContent_Declare(
			crpropa
			GIT_REPOSITORY "https://github.com/CRPropa/CRPropa3.git"
			GIT_TAG master
		)
		list(APPEND CMAKE_MODULE_PATH "${CRPropa_SOURCE_DIR}/cmake")
		FetchContent_MakeAvailable(crpropa)
	endif(NOT crpropa_POPULATED)

	set(CRPropa_SOURCE_DIR "${crpropa_SOURCE_DIR}")
	set(CRPropa_BINARY_DIR "${crpropa_BINARY_DIR}")

	list(APPEND grplinst_EXTRA_INCLUDES ${CRPropa_SOURCE_DIR}/include)
	list(APPEND grplinst_EXTRA_INCLUDES ${CRPropa_BINARY_DIR}/include)
	list(APPEND grplinst_EXTRA_INCLUDES ${CRPROPA_EXTRA_INCLUDES})
	list(APPEND grplinst_EXTRA_LIBRARIES crpropa)
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SOURCE_DIR}/python)
else()
	find_package(CRPropa REQUIRED)
	list(APPEND grplinst_EXTRA_INCLUDES ${CRPropa_INCLUDE_DIR})
	list(APPEND grplinst_EXTRA_LIBRARIES ${CRPropa_LIBRARY})
	set(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SWIG_PATH})
endif(USE_OWN_CRPROPA)

add_library(grplinst SHARED
	${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Medium.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Flow.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PlasmaInstability.cc
	${grplinst_EXTRA_SOURCES}
)

# the plugin's own headers plus any extra includes (e.g. CRPropa's CMakeLists.txt)
target_include_directories(grplinst PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${grplinst_EXTRA_INCLUDES}
)

# link crpropa and other needed libraries
target_link_libraries(grplinst ${grplinst_EXTRA_LIBRARIES})

# Python interface
if(ENABLE_PYTHON)
	if(Python_FOUND)
		# use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
		option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On)

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

# generate tests if ENABLE_TESTING is set to true
if(ENABLE_TESTING)
	enable_testing()
	set(grplinst_TESTS_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "")  
	mark_as_advanced(grplinst_TESTS_PATH)

	# temporary workaround for newer clang versions due to its handling of unwinding
	#   see: https://github.com/google/googletest/issues/3062
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
		if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
			add_definitions(-Dgrplinst_TESTS_SKIP_EXCEPTIONS)
		endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

	if(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
	endif(APPLE)

	# Reuse CRPropa's own gtest if at all possible; only if that fails, look for the system's gtest (or download it)
	set(grplinst_GTEST_LIBRARIES)
	set(grplinst_GTEST_INCLUDE_DIRS)

	# install CRPropa together with plugin
	if(USE_OWN_CRPROPA)
		set(grplinst_GTEST_INCLUDE_DIRS "${CRPropa_SOURCE_DIR}/libs/gtest/include")
		set(grplinst_GTEST_LIBRARIES gtest gtest_main)
	else()
		# check if CRPropa's bundled gtest is available (it should be, if CRPropa was built with ENABLE_TESTING, for versions of CRPropa >= 3.3) 
		# Resolve GTest from CRPropa's bundle first, else fall back to system.
		find_library(gtestCRPropa     NAMES gtest      PATHS "${CRPropa_LIBRARY_DIR}" NO_DEFAULT_PATH)
		find_library(gtestCRPropaMain NAMES gtest_main PATHS "${CRPropa_LIBRARY_DIR}" NO_DEFAULT_PATH)
		if(IS_DIRECTORY "${CRPropa_INCLUDE_DIR}/gtest" AND gtestCRPropa AND gtestCRPropaMain)
			set(grplinst_GTEST_INCLUDE_DIRS "${CRPropa_INCLUDE_DIR}")
			set(grplinst_GTEST_LIBRARIES    "${gtestCRPropa}" "${gtestCRPropaMain}")
			message(STATUS "GTest: using CRPropa's bundled version, (${gtestCRPropa}, ${gtestCRPropaMain})")
		# otherwise try to find a system-installed gtest package (or download it if not found)
		else()
			find_package(GTest QUIET)
			if(TARGET GTest::gtest)
				set(grplinst_GTEST_LIBRARIES GTest::gtest)
				message(STATUS "GTest: using system's GTest package")
			elseif(TARGET GTest::GTest)
				set(grplinst_GTEST_LIBRARIES GTest::GTest)
				message(STATUS "GTest: using system's legacy GTest package.")
			else()
				include(FetchContent)
				FetchContent_Declare(
					googletest
					GIT_REPOSITORY https://github.com/google/googletest.git
					GIT_TAG v1.14.0 # v1.17.0 fails on some systems
				)
				FetchContent_MakeAvailable(googletest)
				set(grplinst_GTEST_LIBRARIES gtest)
				list(APPEND grplinst_GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include" "${gtest_SOURCE_DIR}")
				message(STATUS "Downloading GTests from Github.")
			endif()
		endif(IS_DIRECTORY "${CRPropa_INCLUDE_DIR}/gtest" AND gtestCRPropa AND gtestCRPropaMain)
	endif(USE_OWN_CRPROPA)

	# function to create C++ tests, which are compiled and linked against gtest and the plugin
	function(addTest testName)
		add_executable(${testName} ${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.cpp)
		target_include_directories(${testName} PRIVATE ${grplinst_GTEST_INCLUDE_DIRS})
		target_link_libraries(${testName} PRIVATE grplinst ${grplinst_GTEST_LIBRARIES} pthread ${COVERAGE_LIBS})
		add_test(NAME ${testName}Cpp COMMAND ${grplinst_TESTS_PATH}/${testName})
	endfunction()

	# C++ tests
	addTest(testMedium testMedium)
	addTest(testFlow testFlow)
	addTest(testPlasmaInstability testPlasmaInstability)


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
					"${crpropa_BINARY_DIR}/crpropa.py"
					"$<TARGET_FILE:crpropa-swig>"
					"$<TARGET_FILE:crpropa>"
					"${CRPropa_PYTHON_TEST_PACKAGE_DIR}"
			)
		endif(USE_OWN_CRPROPA)

		function(addTestPy testName)
			configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/${testName}.py ${CMAKE_CURRENT_BINARY_DIR}/${testName}.py COPYONLY)
			add_test(NAME ${testName}Py COMMAND ${Python_EXECUTABLE} ${grplinst_TESTS_PATH}/${testName}.py)
			if(USE_OWN_CRPROPA)
				set_tests_properties(${testName}Py PROPERTIES DEPENDS moveCRPropaPlugin)
			endif(USE_OWN_CRPROPA)
		endfunction()

		# the python tests
		addTestPy(testPython)

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


if(ENABLE_PYBIND11_STUBGEN AND ENABLE_PYTHON AND Python_FOUND)
	find_program(PYBIND11_STUBGEN_EXECUTABLE NAMES "pybind11-stubgen")
	if(PYBIND11_STUBGEN_EXECUTABLE)
		message(STATUS "Found pybind11-stubgen at ${PYBIND11_STUBGEN_EXECUTABLE}")
		install(CODE "execute_process(COMMAND ${PYBIND11_STUBGEN_EXECUTABLE} --ignore-all-errors -o ${Python_SITELIB} grplinst)")
		install(CODE "execute_process(COMMAND echo -- Generated python stubs over pybind11-stubgen for grplinst)")
	endif(PYBIND11_STUBGEN_EXECUTABLE)
endif(ENABLE_PYBIND11_STUBGEN AND ENABLE_PYTHON AND Python_FOUND)
