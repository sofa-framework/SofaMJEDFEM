cmake_minimum_required(VERSION 2.8.12)

find_package(SofaFramework REQUIRED)
sofa_add_plugin(SofaMJEDFEM/ SofaMJEDFEM)

#sofa_install_git_version( "sofa-dev" ${CMAKE_CURRENT_SOURCE_DIR} )
