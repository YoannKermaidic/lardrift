find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(PCL 1.2 REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

include("${ROOT_USE_FILE}")

include_directories(${PCL_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/libs/inc)
link_directories(${PCL_LIBRARY_DIRS} ${PROJECT_SOURCE_DIR}/libs/src)

add_definitions(${PCL_DEFINITIONS})

add_executable (lardrift lardrift.cpp)
target_link_libraries (lardrift libs ${PCL_LIBRARIES} ${ROOT_LIBRARIES} nlohmann_json::nlohmann_json)
