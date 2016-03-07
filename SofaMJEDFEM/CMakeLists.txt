cmake_minimum_required(VERSION 2.8)

project(SofaMJEDFEM)
set(SOFAMJEDFEM_MAJOR_VERSION 1)
set(SOFAMJEDFEM_MINOR_VERSION 0)
set(SOFAMJEDFEM_VERSION ${SOFAMJEDFEM_MAJOR_VERSION}.${SOFAMJEDFEM_MINOR_VERSION})

find_package(SofaFramework REQUIRED)

set(HEADER_FILES

    initMJEDFEM.h 
    MJEDTetrahedralForceField.h 
    MJEDTetrahedralForceField.inl 
    fem/material/NeoHookeanIsotropicMJED.h 
    fem/material/HyperelasticMaterialMJED.h 
    fem/material/BoyceAndArrudaMJED.h 
    fem/material/STVenantKirchhoffMJED.h 
    fem/material/NeoHookeanMJED.h 
    fem/material/MooneyRivlinMJED.h 
    fem/material/VerondaWestmanMJED.h 
    fem/material/CostaMJED.h 
    fem/material/OgdenMJED.h    

	)

set(SOURCE_FILES
	initMJEDFEM.cpp 
	MJEDTetrahedralForceField.cpp
	)

set(README_FILE "SofaMJEDFEM.txt")

set(SOFAMJEDFEM_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${SOFAMJEDFEM_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_MJED_FEM")

target_link_libraries(${PROJECT_NAME} SofaBaseTopology)

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${SOFAMJEDFEM_INCLUDE_DIRS}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

## Install rules for the library and the headers; CMake package configurations files
sofa_create_package(${PROJECT_NAME} ${SOFAMJEDFEM_VERSION} ${PROJECT_NAME} ${PROJECT_NAME})

install(DIRECTORY sofa
        DESTINATION include/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.h" PATTERN "*hpp" PATTERN "*.inl"
)
