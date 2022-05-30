find_package(ROOT)
find_package(PCL 1.2 REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/libs/inc)

add_library(libs plotter.cpp reader.cpp simarrays.cpp tools.cpp)