CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(deal.II 8.5.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(myproject)
include_directories(include)

file(GLOB SOURCES "src/*.cpp")
ADD_EXECUTABLE(ns ${SOURCES})
DEAL_II_SETUP_TARGET(ns)
