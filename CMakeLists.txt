#cmake 版本限制
cmake_minimum_required(VERSION 3.0)
#工程名
project(PPPINS)

#编译模式 使用 IF(NOT ) ENDIF 放置重复设置
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

# 显示 编译模式信息
MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 -march=native")

# 检查c++11或者 c++0x 编译支持  Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

add_definitions(-DENAGLO -DENACMP -DENAGAL -DTRACE -D_DEBUG -D_CRT_SECURE_NO_DEPRECATE)

include_directories(
${PROJECT_SOURCE_DIR}/include
)

# yaml-cpp
find_package(yaml-cpp REQUIRED)

# 连接库文件
LINK_LIBRARIES(m)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(RTKLIB SHARED
src/options.c
src/rtkcmn.c

src/filepos.c
src/rinex.c

src/trop.c
src/ione.c
src/preceph.c
src/ephemeris.c

src/function.c
src/pntpos.c
src/ppppos.c
src/rtkpos.c

src/lambda.c
src/pppar.c

src/solution.c
src/download.c
src/geoid.c
src/tides.c
src/convkml.c
src/tdcpos.c

src/CGNSSConfig.cpp
src/CGNSSManage.cpp
)
if (CMAKE_SYSTEM_NAME MATCHES "Windows") 
target_link_libraries(RTKLIB Winmm.lib ${YAML_CPP_LIBRARIES}) 
else()
target_link_libraries(RTKLIB ${YAML_CPP_LIBRARIES})
endif()

add_library(PSINS SHARED
src/PSINS.cpp
src/CIMUConfig.cpp
)
target_link_libraries(PSINS ${YAML_CPP_LIBRARIES})       # YAML默认是静态库

# Build examples
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
add_executable(downobs 
Example/Tool/downobs.c
)
target_link_libraries(downobs RTKLIB)

add_executable(binfile 
Example/Tool/binfile.c
)
target_link_libraries(binfile RTKLIB)

add_executable(fcbpos
Example/Tool/fcbpos.c
)
target_link_libraries(fcbpos RTKLIB)

add_executable(difpos
Example/Tool/difpos.c
)
target_link_libraries(difpos RTKLIB)

add_executable(ionpos
Example/Tool/ionpos.c
)
target_link_libraries(ionpos RTKLIB)

add_executable(selsta
Example/Tool/selsta.c
)
target_link_libraries(selsta RTKLIB)

add_executable(ra
Example/Tool/ra.c
)
target_link_libraries(ra RTKLIB)

add_executable(qual
Example/Tool/qual.c
)
target_link_libraries(qual RTKLIB)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
add_executable(PPP 
Example/GNSS/PPP.c
)
target_link_libraries(PPP RTKLIB)

add_executable(RTK 
Example/GNSS/RTK.c
)
target_link_libraries(RTK RTKLIB)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
add_executable(GNSSINS 
Example/GNSS-INS/CKFApp.h
Example/GNSS-INS/CKFApp.cpp
Example/GNSS-INS/GNSSINS.cpp
)

target_link_libraries(GNSSINS RTKLIB PSINS)


