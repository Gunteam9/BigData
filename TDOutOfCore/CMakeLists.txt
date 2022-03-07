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

## Il faut mettre a jour la position de vos fichiers sn_resamp512 et eventuellement sn_resamp1024 dans le CMakeList.txt ici:
## set (MY_MESHES_PATH /home/etu/votre_chemin_ou_il_y_a_sn_resamp/)

set (MY_MESHES_PATH /share/etud)


configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

include_directories("${PROJECT_BINARY_DIR}")


include(${VTK_USE_FILE})
 

file(GLOB SRC_LIST ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )



add_executable(${PROJECT_NAME} ${SRC_LIST})

 
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})


enable_testing()
add_test(SimpleTest ${PROJECT_NAME})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS ${PROJECT_NAME})
