cmake_minimum_required(VERSION 3.0)

project(boost-odeint)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=5")

set(CMAKE_BUILD_TYPE debug)

# Find Boost
#find_package(Boost 1.74.0 COMPONENTS filesystem)
find_package(Boost REQUIRED)

if(Boost_FOUND)
  message(STATUS "Boost_INCLUDE_DIRS = '${Boost_INCLUDE_DIRS}'")
  message(STATUS "Boost_LIBRARIES = '${Boost_LIBRARIES}'")
else()
  message(STATUS "Boost was NOT found")
endif()

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/lib/tic-toc-timer/include
)

set(SRC
  #${CMAKE_SOURCE_DIR}/src/lorenz.cpp
)

set(timer_SRC
  ${CMAKE_SOURCE_DIR}/lib/tic-toc-timer/src/tic_toc_timer.cpp
)

add_library(tictoctimer ${timer_SRC})


add_executable(odeint src/lorenz.cpp)
target_link_libraries(odeint ${Boost_LIBRARIES} tictoctimer)


