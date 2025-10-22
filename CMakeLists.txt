#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)

project(g4main)
set(CMAKE_CXX_STANDARD 11)

#----------------------------------------------------------------------------
# Find EXPAT package (ignore version requirement)
# You can optionally specify a custom EXPAT location with EXPAT_ROOT
find_package(EXPAT QUIET)  

#----------------------------------------------------------------------------
# Optionally build Geant4 with UI/Vis drivers
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
find_package(ROOT)
if(NOT ROOT_FOUND)
  message(STATUS "G4 Examples: ROOT package not found. --> g4main example disabled")  
  return()
else()
  message(STATUS "Found ROOT")  
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR}
    /opt/root/include/
)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(g4main g4main.cc ${sources} ${headers})
target_link_libraries(g4main ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy scripts to the build directory
set(g4main_SCRIPTS
    vis.mac
)

foreach(_script ${g4main_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS g4main DESTINATION bin)

