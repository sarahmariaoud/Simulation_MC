cmake_minimum_required(VERSION 3.16.0)


file(READ ${CMAKE_CURRENT_SOURCE_DIR}/config.json MY_JSON_STRING)
string(JSON LOAD_OMP GET ${MY_JSON_STRING} openmp)
message("-- LOAD_OMP: ${LOAD_OMP}")

project(TP 	VERSION 1.0 
			DESCRIPTION "simulation code for multicomponent model"
			LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(TP.out src/main.cpp src/System.cpp)

target_compile_features(TP.out PRIVATE cxx_std_17)

find_package(OpenMP)

if(OpenMP_CXX_FOUND AND LOAD_OMP STREQUAL "true")
	if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
		    target_link_libraries(TP.out PUBLIC OpenMP::OpenMP_CXX)
	elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		    target_link_libraries(TP.out PUBLIC OpenMP::OpenMP_CXX stdc++fs)
	endif()
else()
	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		    target_link_libraries(TP.out PUBLIC stdc++fs)
	endif()
endif()



