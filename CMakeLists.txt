cmake_minimum_required(VERSION 2.12)

project(bingham_mac)

find_package(inmost REQUIRED)
if(inmost_FOUND)
	message("inmost found: include ${INMOST_INCLUDE_DIRS}")
	include_directories("${INMOST_INCLUDE_DIRS}")
endif()

add_executable(main main.cpp analytics.cpp)
target_link_libraries(main ${INMOST_LIBRARIES})

add_executable(cavity cavity.cpp analytics.cpp)
target_link_libraries(cavity ${INMOST_LIBRARIES})

add_executable(gmg_cavity cavity_gmg.cpp analytics.cpp)
target_link_libraries(gmg_cavity ${INMOST_LIBRARIES})
