cmake_minimum_required(VERSION 2.8)
 


set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")



PROJECT(vtktp)



IF(APPLE)

SET( OperatingSystem "Mac os X")

find_package(VTK 7 REQUIRED)

set (MY_MESHES_PATH  /Users/emmanuelmelin/shared_folders/data_VTK/DEVOIR)


ENDIF (APPLE)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

SET(OperatingSystem "Linux")

find_package(VTK )
set (MY_MESHES_PATH /share/etud/mystere/)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")



MESSAGE(WARNING "System detected:${OperatingSystem}" )
MESSAGE(WARNING "Path to sn_resamp:${MY_MESHES_PATH}" )



configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

include_directories("${PROJECT_BINARY_DIR}")


include(${VTK_USE_FILE})
 
 
 
 
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

 
 

file(GLOB SRC_LIST ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp )



add_executable(${PROJECT_NAME} ${SRC_LIST})

 

target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${MPI_LIBRARIES})


if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()


add_custom_target(prun1
    COMMAND mpirun    -np 1 ${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(prun2
    COMMAND mpirun   -np 2 ${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(prun4
    COMMAND mpirun    -np 4 ${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
