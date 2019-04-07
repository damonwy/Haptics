CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Name of the project
PROJECT(Haptics)

# FOR LAB MACHINE
SET(DEF_DIR_GLM "D:\\Libraries\\glm-0.9.8.5\\glm")
SET(DEF_DIR_GLFW "D:\\Libraries\\glfw-3.2.1")
SET(DEF_DIR_GLEW "D:\\Libraries\\glew-2.1.0")
SET(DEF_DIR_EIGEN "D:\\Libraries\\eigen-eigen-6c3592226be0")
#SET(DEF_DIR_TETGEN "D:\\Libraries\\tetgen1.5.1-beta1")
SET(DEF_DIR_JSONHPP_INCLUDE "D:\\Libraries\\json-develop\\include")
SET(DEF_DIR_MOSEK "D:\\Libraries\\mosek_8.1_x64\\8\\tools\\platform\\win64x86")
#SET(DEF_DIR_ODEINT "D:\\Libraries\\odeint-v2-2.6_boost1.60\\include")
SET(DEF_DIR_BOOST "D:\\Libraries\\boost_1_65_1")
SET(DEF_DIR_LIBTORCH "D:\\Libraries\\libtorch-win-shared-with-deps-latest\\libtorch")
#SET(DEF_DIR_LIBIGL_INCLUDE "D:\\Libraries\\libigl-master\\libigl\\include")
#SET(DEF_DIR_PARDISO "D:\\Libraries\\PARDISO")
#SET(DEF_DIR_MKL "C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries_2018\\windows\\mkl")
#SET(DEF_DIR_JSONCPP "D:\\Libraries\\jsoncpp\\jsoncpp_x64-windows")
#SET(DEF_DIR_STB_INCLUDE "D:\\Libraries\\stb-master")
#SET(DEF_DIR_SUPERLU "D:\\Libraries\\superlu_5.2.1")
#SET(DEF_DIR_CMAES "D:\\Libraries\\CMA-ESpp-master")
#SET(DEF_DIR_LIBLBFGS "D:\\Libraries\\liblbfgs-master")
#SET(DEF_DIR_GUROBI "D:\\Libraries\\gurobi\\win64")
#SET(DEF_DIR_YIXUAN "D:\\Libraries\\yixuan-LBFGSpp-ad499fe")
#SET(DEF_DIR_DLIB "D:\\Libraries\\dlib-19.13")
#SET(DEF_DIR_QUADPROGPP "D:\\Libraries\\quadprogpp")
#SET(DEF_DIR_NLOPT "D:\\Libraries\\nlopt")

# FOR MAC
#SET(DEF_DIR_GLM "/Users/yingwang/Libraries/glm-0.9.9.3")
#SET(DEF_DIR_GLFW "/Users/yingwang/Libraries/glfw-master")
#SET(DEF_DIR_GLEW "/Users/yingwang/Libraries/glew-2.1.0")
#SET(DEF_DIR_EIGEN "/Users/yingwang/Libraries/eigen-eigen-08a1f23d4f46")
#SET(DEF_DIR_TETGEN "/Users/yingwang/Libraries/tetgen1.5.1")
#SET(DEF_DIR_MOSEK "/Users/yingwang/Libraries/mosek/8/tools/platform/osx64x86")
#SET(DEF_DIR_LIBIGL_INCLUDE "/Users/yingwang/Libraries/libigl-2.0.0/include")
#SET(DEF_DIR_PARDISO "/Users/yingwang/Libraries/pardiso")
#SET(DEF_DIR_MKL "/opt/intel/mkl")
#SET(DEF_DIR_JSONCPP "/Users/yingwang/Libraries/jsoncpp-0.10.7")
#SET(DEF_DIR_STB_INCLUDE "/Users/yingwang/Libraries/stb-master")
#SET(DEF_DIR_JSONHPP_INCLUDE "/Users/yingwang/Libraries/json-3.5.0/include/")

# Is this the solution?
# Override with `cmake -DSOL=ON ..`
OPTION(SOL "Solution" OFF)

# Use glob to get the list of all source files.
# We don't really need to include header and resource files to build, but it's
# nice to have them also show up in IDEs.
IF(${SOL})
  FILE(GLOB_RECURSE SOURCES "src0/*.cpp")
  FILE(GLOB_RECURSE HEADERS "src0/*.h")
ELSE()
  FILE(GLOB_RECURSE SOURCES "src/*.cpp")
  FILE(GLOB_RECURSE HEADERS "src/*.h")
ENDIF()
FILE(GLOB_RECURSE GLSL "resources/*.glsl")

# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

# Get JSON
SET(JSON_HPP_DIR "$ENV{JSON_HPP_DIR}")
IF(NOT JSON_HPP_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable JSON_HPP_DIR to the root directory of your JSON installation")
  IF(APPLE)
    MESSAGE(STATUS "Looking for JSONHPP in ${DEF_DIR_JSONHPP_INCLUDE}")
    IF(IS_DIRECTORY ${DEF_DIR_JSONHPP_INCLUDE})
      MESSAGE(STATUS "Found!")
      SET(JSON_HPP_DIR ${DEF_DIR_JSONHPP_INCLUDE})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSEIF(WIN32)
	MESSAGE(STATUS "Looking for JSONHPP in ${DEF_DIR_JSONHPP_INCLUDE}")
    IF(IS_DIRECTORY ${DEF_DIR_JSONHPP_INCLUDE})
      MESSAGE(STATUS "Found!")
      SET(JSON_HPP_DIR ${DEF_DIR_JSONHPP_INCLUDE})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${JSON_HPP_DIR})

# Get the GLM environment variable. Since GLM is a header-only library, we
# just need to add it to the include directory.
SET(GLM_INCLUDE_DIR "$ENV{GLM_INCLUDE_DIR}")
IF(NOT GLM_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLM_INCLUDE_DIR to the root directory of your GLM installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLM in ${DEF_DIR_GLM}")
    IF(IS_DIRECTORY ${DEF_DIR_GLM})
      MESSAGE(STATUS "Found!")
      SET(GLM_INCLUDE_DIR ${DEF_DIR_GLM})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Looking for GLM in ${DEF_DIR_GLM}")
    IF(IS_DIRECTORY ${DEF_DIR_GLM})
      MESSAGE(STATUS "Found!")
      SET(GLM_INCLUDE_DIR ${DEF_DIR_GLM})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLM_INCLUDE_DIR})

# Get the GLFW environment variable. There should be a CMakeLists.txt in the 
# specified directory.
SET(GLFW_DIR "$ENV{GLFW_DIR}")
IF(NOT GLFW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLFW_DIR to the root directory of your GLFW installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLFW in ${DEF_DIR_GLFW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLFW})
      MESSAGE(STATUS "Found!")
      SET(GLFW_DIR ${DEF_DIR_GLFW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Looking for GLFW in ${DEF_DIR_GLFW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLFW})
      MESSAGE(STATUS "Found!")
      SET(GLFW_DIR ${DEF_DIR_GLFW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ENDIF()
ENDIF()
OPTION(GLFW_BUILD_EXAMPLES "GLFW_BUILD_EXAMPLES" OFF)
OPTION(GLFW_BUILD_TESTS "GLFW_BUILD_TESTS" OFF)
OPTION(GLFW_BUILD_DOCS "GLFW_BUILD_DOCS" OFF)
IF(CMAKE_BUILD_TYPE MATCHES Release)
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/release)
ELSE()
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/debug)
ENDIF()
INCLUDE_DIRECTORIES(${GLFW_DIR}/include)
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} glfw ${GLFW_LIBRARIES})

# Get the GLEW environment variable.
SET(GLEW_DIR "$ENV{GLEW_DIR}")
IF(NOT GLEW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLEW_DIR to the root directory of your GLEW installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLEW in ${DEF_DIR_GLEW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLEW})
      MESSAGE(STATUS "Found!")
      SET(GLEW_DIR ${DEF_DIR_GLEW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
  MESSAGE(STATUS "Looking for GLEW in ${DEF_DIR_GLEW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLEW})
      MESSAGE(STATUS "Found!")
      SET(GLEW_DIR ${DEF_DIR_GLEW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()  
    ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLEW_DIR}/include)
IF(WIN32)
  # With prebuilt binaries
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_DIR}/lib/Release/x64/glew32s.lib)
ELSE()
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${GLEW_DIR}/lib/ libGLEW.a libGLEW.2.1.0.dylib)
ENDIF()

# Get the EIGEN environment variable. Since EIGEN is a header-only library, we
# just need to add it to the include directory.
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable EIGEN3_INCLUDE_DIR to the root directory of your EIGEN installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for EIGEN in ${DEF_DIR_EIGEN}")
    IF(IS_DIRECTORY ${DEF_DIR_EIGEN})
      MESSAGE(STATUS "Found!")
      SET(EIGEN3_INCLUDE_DIR ${DEF_DIR_EIGEN})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Looking for EIGEN in ${DEF_DIR_EIGEN}")
    IF(IS_DIRECTORY ${DEF_DIR_EIGEN})
      MESSAGE(STATUS "Found!")
      SET(EIGEN3_INCLUDE_DIR ${DEF_DIR_EIGEN})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})


# Get the MOSEK environment variable.
SET(MOSEK_DIR "$ENV{MOSEK_DIR}")
IF(NOT MOSEK_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable MOSEK_DIR to the root directory of your MOSEK installation")
  IF(WIN32)
	INCLUDE_DIRECTORIES(${DEF_DIR_MOSEK}/8/tools/platform/win64x86/h)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${DEF_DIR_MOSEK}/8/tools/platform/win64x86/bin/mosek64_8_1.lib)
  ELSEIF(APPLE)
    INCLUDE_DIRECTORIES(${DEF_DIR_MOSEK}/h)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${DEF_DIR_MOSEK}/bin/ libmosek64.8.1.dylib)
  ENDIF()
ELSE()
  IF(WIN32)
	INCLUDE_DIRECTORIES(${MOSEK_DIR}/8/tools/platform/win64x86/h)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MOSEK_DIR}/8/tools/platform/win64x86/bin/mosek64_8_1.lib)
  ELSEIF(APPLE)
    INCLUDE_DIRECTORIES(${DEF_DIR_MOSEK}/h)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${DEF_DIR_MOSEK}/bin/ libmosek64.8.1.dylib)
  ENDIF()
ENDIF()



SET(Torch_DIR ${DEF_DIR_LIBTORCH}/share/cmake/Torch)

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (default Debug)")
find_package(Torch REQUIRED)
if(NOT Torch_FOUND)
    message(FATAL_ERROR "Pytorch Not Found!")
endif(NOT Torch_FOUND)
message(STATUS "Pytorch status:")
message(STATUS "    libraries: ${TORCH_LIBRARIES}")
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${TORCH_LIBRARIES})
set_property(TARGET ${CMAKE_PROJECT_NAME} PROPERTY CXX_STANDARD 11)

# # Get the odeint environment variable. Since odeint is a header-only library, we
# # just need to add it to the include directory.
# SET(ODEINT_INCLUDE_DIR "$ENV{ODEINT_INCLUDE_DIR}")
# IF(NOT ODEINT_INCLUDE_DIR)
#   # The environment variable was not set
#   SET(ERR_MSG "Please point the environment variable ODEINT_INCLUDE_DIR to the root directory of your odeint installation.")
#   IF(WIN32)
#     # On Windows, try the default location
#     MESSAGE(STATUS "Looking for odeint in ${DEF_DIR_ODEINT}")
#     IF(IS_DIRECTORY ${DEF_DIR_ODEINT})
#       MESSAGE(STATUS "Found!")
#       SET(EIGEN3_INCLUDE_DIR ${DEF_DIR_ODEINT})
#     ELSE()
#       MESSAGE(FATAL_ERROR ${ERR_MSG})
#     ENDIF()
#   ELSE()
#     MESSAGE(FATAL_ERROR ${ERR_MSG})
#   ENDIF()
# ENDIF()
# INCLUDE_DIRECTORIES(${ODEINT_INCLUDE_DIR})

# Get the odeint environment variable. Since odeint is a header-only library, we
# just need to add it to the include directory.
SET(BOOST_INCLUDE_DIR "$ENV{BOOST_INCLUDE_DIR}")
IF(NOT BOOST_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable BOOST_INCLUDE_DIR to the root directory of your odeint installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for odeint in ${DEF_DIR_BOOST}")
    IF(IS_DIRECTORY ${DEF_DIR_BOOST})
      MESSAGE(STATUS "Found!")
      SET(BOOST_INCLUDE_DIR ${DEF_DIR_BOOST})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})




# OS specific options and libraries
IF(WIN32)
  # c++11 is enabled by default.
  # -Wall produces way too many warnings.
  # -pedantic is not supported.
  # Disable warning 4996.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} opengl32.lib)
ELSE()
  # Enable all pedantic warnings.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
  IF(APPLE)
    # Add required frameworks for GLFW.
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  ELSE()
    #Link the Linux OpenGL library
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "GL")
  ENDIF()
ENDIF()