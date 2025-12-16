# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Geant4Garfield)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# ########################################################################################

# ---Find Geant4 package, activating all available UI and Vis drivers by default- You can
# set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui to build a batch
# mode only executable

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 COMPONENTS ui_all vis_all)
else()
  find_package(Geant4)
endif()

if(NOT Geant4_FOUND)
  message("Geant4 is not found. Example cannot be built.")
  return()
endif()

find_package(ROOT REQUIRED)
if(ROOT_FOUND)
  # message("ROOT find...")
  include(${ROOT_DIR}/ROOTUseFile.cmake)
  message(${ROOT_INCLUDE_DIRS})
  message(${ROOT_LIBRARIES})
endif()

# ---Setup Geant4 include directories and compile definitions--------------------
include(${Geant4_USE_FILE} ${ROOT_INCLUDE_DIRS})

# ---Locate sources and headers for this project---------------------------------
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

# ---Add the executable, and link it to the Geant4 libraries---------------------
add_executable(MUSIC MUSIC.cc ${sources} ${headers})
target_link_libraries(MUSIC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Garfield::Garfield)
target_include_directories(MUSIC
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# ---Copy all scripts to the build directory-------------------------------------
set(EXAMPLEGARFIELD_SCRIPTS
    icons.mac
    gui.mac
    physics.mac
    init_vis.mac
    run.mac
    vis.mac)

foreach(_script ${EXAMPLEGARFIELD_SCRIPTS})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_script}
                 ${CMAKE_CURRENT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

# ---Copy all files locally to the build directory-------------------------------
foreach(_file
GasFile/dielectrics.dat
GasFile/alpha.txt 
GasFile/gas.gas
GasFile/ar_90_ch4_10_760Torr.gas
GasFile/ar_90_ch4_10_O2_20ppm.gas
GasFile/ar_90_ch4_10_O2_20ppm_830mbar.gas
GasFile/ar_90_ch4_10_O2_100ppm_830mbar.gas
GasFile/ar_90_ch4_10_O2_500ppm.gas
GasFile/ar_90_ch4_10_O2_300ppm.gas
ComsolFile/F4MUSIC.STL
ComsolFile/pAnode1.txt 
ComsolFile/pCathode.txt
ComsolFile/wfield.txt
ComsolFile/mesh.mphtxt
)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
