
CMAKE_MINIMUM_REQUIRED( VERSION 2.8.4 )

project (vq3demo) 
SET(PACKAGE_VERSION "2.20")
SET(PACKAGE_LICENSE "GPL")
SET(PROJECT_DESCRIPTION_SUMMARY "Demo for vq3")

SET(CMAKE_BUILD_TYPE Release)
SET(PROJECT_CFLAGS "-std=c++20 -O3")

 
###################################
# Installation common variables
###################################

SET(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${CMAKE_PROJECT_NAME}/" CACHE PATH "The subdirectory to the header prefix")

###################################
# Dependencies
###################################

find_package(PkgConfig)
pkg_check_modules(VQ3 vq3 REQUIRED)
pkg_check_modules(DEMO demo2d REQUIRED)

set(PKG_CONFIG_DEPENDS "vq3 demo2d")


###################################
# Pkg-config stuff
###################################

FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.pc
"
Name: ${CMAKE_PROJECT_NAME}
Description: ${PROJECT_DESCRIPTION_SUMMARY}
Version: ${PACKAGE_VERSION}
Requires: ${PKG_CONFIG_DEPENDS}
Libs: 
Cflags: -I${INCLUDE_INSTALL_DIR}
"
)

# Define the path to install the pc pkg-config file
SET(PKGCONFIG_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig/" CACHE STRING "Base directory for pkgconfig files")

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.pc
    	DESTINATION ${PKGCONFIG_INSTALL_DIR}
	COMPONENT devel)

###################################
#  Subdirectories
###################################

add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(examples)
