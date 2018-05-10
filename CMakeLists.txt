cmake_minimum_required(VERSION 2.8)

set(CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE x64)
set(CMAKE_GENERATOR_PLATFORM x64)

# Find MATLAB applications
if(WIN32)
file(GLOB MATLAB_APPS "C:/Program Files/MATLAB/R*")
set(MEXEXT ".mexw64")
endif()
if(APPLE)
file(GLOB MATLAB_APPS "/Applications/MATLAB_R*.app")
set(MEXEXT ".mexmaci64")
endif()
if(NOT MATLAB_APPS)
        message(WARNING "No MATLAB installation found, skipping building MATLAB bindings")
    return()
endif()
# Use most recent version
list(SORT MATLAB_APPS)
list(GET MATLAB_APPS -1 MATLAB_ROOT)

project(ojwul-mex)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()
message("Build type: " ${CMAKE_BUILD_TYPE})

# set locations for build products
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
endforeach()

# set some cmake options
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)
set(CMAKE_MACOSX_RPATH 1)

# header files
include_directories(
    "${MATLAB_ROOT}/extern/include"
)

# Libraries
set(LIB_SEARCH_PATHS  ${MATLAB_ROOT}/bin ${MATLAB_ROOT}/extern/lib ${MATLAB_ROOT}/sys/os)
set(LIB_SUFFIXES glnxa64 glnx86 win64 win64/microsoft win32/microsoft maci64)
find_library( MATLAB_MEX_LIBRARY NAMES libmex mex PATHS ${LIB_SEARCH_PATHS} PATH_SUFFIXES ${LIB_SUFFIXES})
find_library( MATLAB_MX_LIBRARY NAMES libmx mx PATHS ${LIB_SEARCH_PATHS} PATH_SUFFIXES ${LIB_SUFFIXES})
find_library( MATLAB_MAT_LIBRARY NAMES libmat mat PATHS ${LIB_SEARCH_PATHS} PATH_SUFFIXES ${LIB_SUFFIXES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -D_MATLAB_=901")

function(add_mex_function)
    string(FIND ${ARGV0} "/" LAST_SLASH REVERSE)
	if( ${LAST_SLASH} LESS 0)
		set(FUNCTION_DIR "")
		set(FUNCTION_NAME ${ARGV0})
	else()
		string(SUBSTRING ${ARGV0} 0 ${LAST_SLASH} FUNCTION_DIR)
		MATH(EXPR LAST_SLASH "${LAST_SLASH}+1")
		string(SUBSTRING ${ARGV0} ${LAST_SLASH} -1  FUNCTION_NAME)
	endif()
    add_library(${FUNCTION_NAME} SHARED ${ARGV0}.cpp)
    set_target_properties(${FUNCTION_NAME} PROPERTIES PREFIX "")
    set_target_properties(${FUNCTION_NAME} PROPERTIES SUFFIX ${MEXEXT})
	#set_target_properties(${FUNCTION_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/${FUNCTION_DIR})
	#set_target_properties(${FUNCTION_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/${FUNCTION_DIR})
    target_link_libraries(${FUNCTION_NAME} PRIVATE ${MATLAB_MX_LIBRARY} ${MATLAB_MEX_LIBRARY} ${MATLAB_MAT_LIBRARY})
	add_custom_command(TARGET ${FUNCTION_NAME}
		POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy 
			$<TARGET_FILE:${FUNCTION_NAME}> 
			${CMAKE_CURRENT_SOURCE_DIR}/${FUNCTION_DIR}/$<TARGET_FILE_NAME:${FUNCTION_NAME}>
	)
endfunction()

# specify individual mex functions here
add_mex_function(image/ojw_interp2)

# Call using:
# cmake ..