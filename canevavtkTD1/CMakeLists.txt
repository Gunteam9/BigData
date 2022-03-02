cmake_minimum_required(VERSION 2.8)
 


set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")



PROJECT(vtktp)



IF(APPLE)

SET( OperatingSystem "Mac os X")

find_package(VTK)



ENDIF (APPLE)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

SET(OperatingSystem "Linux")

find_package(VTK)


ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")


include(${VTK_USE_FILE})
 

file(GLOB SRC_LIST ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )



add_executable(${PROJECT_NAME} ${SRC_LIST} vtkMyCallback.cpp vtkMyCallback.h)

 
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})


enable_testing()
add_test(SimpleTest ${PROJECT_NAME})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS ${PROJECT_NAME})
