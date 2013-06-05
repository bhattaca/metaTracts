# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(BundleEx)


find_package(ITK REQUIRED)
IF(ITK_FOUND)
	INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)



add_executable(xBundle xBundle_main.cxx xBundle_hessian.cxx )

target_link_libraries(xBundle ${ITK_LIBRARIES})
