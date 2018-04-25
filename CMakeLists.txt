#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DICOM)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Add dicomReader subdirectory
#
if("$ENV{DICOM_USE_DCMTK}")
  find_package(DCMTK REQUIRED)
  add_definitions(-DG4_DCMTK)
  add_subdirectory(dicomReader)
  set(DICOM_READER_LIBRARY dicomReader)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/dicomReader/include 
                    ${Geant4_INCLUDE_DIR}
                    ${DCMTK_INCLUDE_DIRS})
                    
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# List any source specific properties here


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(DICOM DICOM.cc ${headers} ${sources})

target_link_libraries(DICOM ${Geant4_LIBRARIES} ${DICOM_READER_LIBRARY} ${DCMTK_LIBRARIES})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build DICOM. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# the macros
set(DICOM_MACROS
    run.mac vis.mac
    )
 
# original set of DICOM data
set(DICOM_SCRIPTS
    ColourMap.dat CT2Density.dat 
    Data.dat 
  )
  

foreach(_script ${DICOM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


foreach(_script ${DICOM_MACROS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

if("$ENV{DICOM_USE_DCMTK}")

set_source_files_properties(
     ${sources}
     PROPERTIES COMPILE_DEFINITIONS G4_DCMTK 

)
endif()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS DICOM DESTINATION bin)

