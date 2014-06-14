# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8.11)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(BundleEx)


#find_package(ITK REQUIRED)
#IF(ITK_FOUND)
#	INCLUDE(${USE_ITK_FILE})
#ENDIF(ITK_FOUND)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
   INCLUDE_DIRECTORIES($ENV{ITK_APPS_SRC_PATH}/QtITK )
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


#VTK

	FIND_PACKAGE(VTK REQUIRED)
	IF(VTK_FOUND)
	   INCLUDE(${VTK_USE_FILE})
	ELSE(VTK_FOUND)
	   MESSAGE(FATAL_ERROR "Cannot build without VTK.  Please set VTK_DIR.")
	ENDIF(VTK_FOUND)

#QT	
#find_package(Qt5Widgets)
	
	
#if (ITKVtkGlue_LOADED)
#  find_package(VTK REQUIRED)
#  include(${VTK_USE_FILE})
#else()
#  find_package(ItkVtkGlue REQUIRED)
#  include(${ItkVtkGlue_USE_FILE})
#  set(Glue ItkVtkGlue)
#endif()

ADD_DEFINITIONS(
    -std=c++11 # Or -std=c++0x
    # Other flags
)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#SET(BOOST_ROOT "C:\\Users\\p41123\\Documents\\internship\\programs\\boost_1_54_0\\boost_1_54_0")
#SET(BOOST_LIBRARYDIR "C:\\Users\\p41123\\Documents\\internship\\programs\\boost_1_54_0\\boost_1_54_0\\libs\\")

#SET(BOOST_ROOT "D:\\ABhattacharya\\Internship\\code\\dependencies\\boost_1_54_0")
#SET(BOOST_LIBRARYDIR "D:\\ABhattacharya\\Internship\\code\\dependencies\\boost_1_54_0\\libs\\")

#Updating the code to boost 1 55 0 (June 2 2014)
SET(BOOST_ROOT "D:\\ABhattacharya\\Internship2014\\Programs\\boost_1_55_0\\boost_1_55_0")
SET(BOOST_LIBRARYDIR "D:\\ABhattacharya\\Internship2014\\Programs\\boost_1_55_0\\boost_1_55_0\\libs\\")


MESSAGE(STATUS "** Search Boost root: ${BOOST_ROOT}")
find_package( Boost REQUIRED )
 MESSAGE(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
 MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
 MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
add_executable(xBundle xBundle_main.cxx xBundle_hessian.cxx xBundle_hessian.h
xBundle_input_param.h xBundle_sparseCoding.cxx xBundle_sparseCoding.h xBundle_fiber_info.h xBundle_fiber_info.cxx)

#qt5_use_modules(xBundle Widgets)
target_link_libraries(xBundle ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${Glue}  ${VTK_LIBRARIES} ${ITK_LIBRARIES})
