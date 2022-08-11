cmake_minimum_required(VERSION 3.12)
project(SofaMJEDFEM VERSION 21.12 LANGUAGES CXX)

find_package(Sofa.Framework REQUIRED)
find_package(Sofa.Component REQUIRED)
find_package(Sofa.GL REQUIRED)

set(SOFAMJEDFEM_SRC_DIR src/SofaMJEDFEM)

set(HEADER_FILES

    ${SOFAMJEDFEM_SRC_DIR}/config.h.in
    ${SOFAMJEDFEM_SRC_DIR}/MJEDTetrahedralForceField.h 
    ${SOFAMJEDFEM_SRC_DIR}/MJEDTetrahedralForceField.inl 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/NeoHookeanIsotropicMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/HyperelasticMaterialMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/BoyceAndArrudaMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/STVenantKirchhoffMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/NeoHookeanMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/MooneyRivlinMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/VerondaWestmanMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/CostaMJED.h 
    ${SOFAMJEDFEM_SRC_DIR}/fem/material/OgdenMJED.h    
	)


set(SOURCE_FILES

	${SOFAMJEDFEM_SRC_DIR}/initSofaMJEDFEM.cpp 
	${SOFAMJEDFEM_SRC_DIR}/MJEDTetrahedralForceField.cpp
	)


set(README_FILE "README.md")


# Create the plugin library
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILE})


# Link the plugin library to its dependency
target_link_libraries(${PROJECT_NAME} Sofa.Core Sofa.Type Sofa.DefaultType Sofa.Component.Topology Sofa.Component.StateContainer Sofa.Component.SolidMechanics.FEM.Hyperelastic Sofa.GL)


# Create package Config, Version & Target files
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )