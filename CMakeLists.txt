#cmake 版本限制
cmake_minimum_required(VERSION 2.8)
#工程名
project(RTKLIB)

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
# 连接库文件
LINK_LIBRARIES(m)
# 运行实例
add_executable(${PROJECT_NAME} 
src/comfun.c
src/download.c
src/ephemeris.c
src/fcbfile.c
src/fcbpos.c
src/file2nav.c
src/geoid.c
src/ionex.c
src/ionoModel.c
src/lambda.c
src/main_nsat.c
src/main.c
src/options.c
src/pntpos7-7.c
src/postpos1.c
src/ppp_ar323.c
src/ppp_corr.c
src/ppp_out.c
src/ppp_pos.c
src/ppp_state.h
src/preceph7-7.c
src/rinex2.c
src/rtkcmn2.c
src/rtklib.h
src/rtkpos.c
src/solution.c
src/tides.c
src/tropx.c
)





