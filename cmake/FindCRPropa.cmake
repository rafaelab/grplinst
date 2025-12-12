# search for Python to help find CRPropa (if Python bindings are available)
find_package(Python COMPONENTS Interpreter Development)


# Determine which swig file to use, depending on CRPropa's installation
# For unknown reasons, if built-in is absent this gives segmentation errors in various systems
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On) 
if(ENABLE_SWIG_BUILTIN)
	set(CRPropa_SWIG_FILE "crpropa-builtin.i")
	set(SWIG_MODE_FLAG "-builtin")
else(ENABLE_SWIG_BUILTIN)
	set(CRPropa_SWIG_FILE "crpropa.i")
	set(SWIG_MODE_FLAG "")
endif(ENABLE_SWIG_BUILTIN)


# Find CRPropa headers
find_path(CRPropa_INCLUDE_DIR CRPropa.h 
	HINTS 
		${CRPropa_INSTALL_PREFIX}/../include
		${CRPropa_INSTALL_PREFIX}/include
		$ENV{CRPropa_DIR}/include
		$ENV{CRPropa_DIR}/build
		$ENV{CRPropa_DIR}/build/include
		crpropa 
		include 
		include/crpropa 
	)

# Find CRPropa associated SWIG files
find_path(CRPropa_SWIG_PATH crpropa.i
	HINTS 
		${CRPropa_INSTALL_PREFIX}/share/crpropa/swig_interface	
		${CRPropa_INSTALL_PREFIX}/../python
		$ENV{CRPropa_DIR}/python
		$ENV{CRPropa_DIR}/build/share/crpropa/swig_interface
		share/crpropa/swig_interface
	)

# Find CRPropa library
find_library(CRPropa_LIBRARY NAMES crpropa libcrpropa
	HINTS  
		${CRPropa_INSTALL_PREFIX}/lib
		$ENV{CRPropa_DIR}/build
		$ENV{CRPropa_DIR}/build/lib
		crpropa
		lib/crpropa 
		crpropa/lib
		)

# Find CRPropa's kiss library and headers
find_library(CRPropa_kiss_LIBRARY NAMES libkiss.a
	HINTS 
		${CRPropa_INSTALL_PREFIX}/libs/kiss
		$ENV{CRPropa_DIR}/build/libs/kiss
		lib/kiss
		lib
	)
find_path(CRPropa_kiss_INCLUDE_DIR kiss/logger.h
	HINTS 
		${CRPropa_INSTALL_PREFIX}/include
		$ENV{CRPropa_DIR}/build/include
		include
	)

# Find CRPropa's HepPID library
find_library(CRPropa_HepPID_LIBRARY NAMES libHepPID.a
	HINTS 
		${CRPropa_INSTALL_PREFIX}/libs/HepPID
		$ENV{CRPropa_DIR}/build/libs/HepPID
		lib/HepPID
		lib
	)
find_path(CRPropa_HepPID_INCLUDE_DIR HepPID/ParticleIDMethods.hh
	HINTS 
		${CRPropa_INSTALL_PREFIX}/include
		$ENV{CRPropa_DIR}/build/include
		include
	)

# # Find CRPropa's Eigen library
# find_path(CRPropa_Eigen_INCLUDE_DIR Eigen/Core
# 	HINTS 
# 		${CRPropa_INSTALL_PREFIX}/include
# 		$ENV{CRPropa_DIR}/build/include
# 		include
# 	)

# Define SWIG interface file
set(CRPropa_SWIG_INTERFACE_FILE "${CRPropa_SWIG_PATH}/${CRPropa_SWIG_FILE}")


# Determine whether CRPropa has really been found
if(NOT (${CRPropa_SWIG_INTERFACE_FILE} STREQUAL "") AND NOT (${CRPropa_INCLUDE_DIR} NOT STREQUAL "") AND NOT (${CRPropa_LIBRARY}  STREQUAL ""))
	set(CRPropa_FOUND True)
else()
	set(CRPropa_FOUND False)
endif()


# If CRPropa has not yet been found, try finding it via Python
if(NOT CRPropa_FOUND)
	if(Python_FOUND AND NOT CRPropa_SWIG_PATH)
		execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE CRPropa_SWIG_PATH)
		if(NOT (CRPropa_SWIG_PATH STREQUAL "") OR NOT CRPropa_SWIG_PATH)
			find_path(CRPropa_SWIG_PATH crpropa.i
				HINTS 
					share/crpropa/swig_interface
					${CRPropa_INSTALL_PREFIX}/share/crpropa/swig_interface	
					${CRPropa_INSTALL_PREFIX}/../python
					$ENV{CRPropa_DIR}/python
					$ENV{CRPropa_DIR}/build/share/crpropa/swig_interface
			)
		endif()
	endif(Python_FOUND AND NOT CRPropa_SWIG_PATH)
endif(NOT CRPropa_FOUND)
list(APPEND CMAKE_PREFIX_PATH ${CRPropa_INSTALL_PREFIX})

# If CRPropa not found, warn the user
if(NOT CRPropa_FOUND)
	message(STATUS "CRPropa could **NOT** be found!!!")
	return()
endif(NOT CRPropa_FOUND)



message(STATUS "CRPropa install prefix: ${CRPropa_INSTALL_PREFIX}")
message(STATUS "CRPropa SWIG interface file: ${CRPropa_SWIG_INTERFACE_FILE}")
message(STATUS "CRPropa include path: ${CRPropa_INCLUDE_DIR}")
message(STATUS "CRPropa library: ${CRPropa_LIBRARY}")
message(STATUS "CRPropa's kiss include path: ${CRPropa_kiss_INCLUDE_DIR}")
message(STATUS "CRPropa's kiss library: ${CRPropa_kiss_LIBRARY}")
message(STATUS "CRPropa's HepPID include path: ${CRPropa_HepPID_INCLUDE_DIR}")
message(STATUS "CRPropa's HepPID library: ${CRPropa_HepPID_LIBRARY}")
# message(STATUS "CRPropa's Eigen include path: ${CRPropa_Eigen_INCLUDE_DIR}")



