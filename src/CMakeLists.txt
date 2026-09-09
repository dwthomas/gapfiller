cmake_minimum_required(VERSION 3.26)
set(CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
project(OpprotunisticSonar)
# Set debugging flags globally
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
# set(BOOST_ROOT /opt/boost-1.88)

# find_package(Boost REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
find_package(GeographicLib REQUIRED)

# add_subdirectory(trivis) 
# add_subdirectory(trivis_plus) 

# add_executable(mkvg mkvg.cpp)
# target_link_libraries(mkvg TrivisPlus)
# TARGET_LINK_LIBRARIES(mkvg ${Boost_LIBRARIES})

# add_executable(find_path find_path.cpp)
# target_link_libraries(find_path TrivisPlus) 
# TARGET_LINK_LIBRARIES(find_path ${Boost_LIBRARIES})

# add_executable(longest_path longest_path.cpp)
# TARGET_LINK_LIBRARIES(longest_path ${Boost_LIBRARIES})

find_package(PkgConfig REQUIRED)
pkg_check_modules(GDAL REQUIRED gdal)

find_package(s2 REQUIRED)

# add_executable(preprocess_maps preprocess_maps.cpp)
# target_include_directories(preprocess_maps PUBLIC ${GDAL_INCLUDE_DIRS})
# target_link_directories(preprocess_maps PUBLIC ${GDAL_LIBRARY_DIRS})
# target_link_libraries(preprocess_maps PUBLIC ${GDAL_LIBRARIES})
# TARGET_LINK_LIBRARIES(preprocess_maps PUBLIC ${Boost_LIBRARIES})

add_executable(local_search local_search.cpp dataset_load.cpp)
target_include_directories(local_search PUBLIC ${GDAL_INCLUDE_DIRS})
target_link_directories(local_search PUBLIC ${GDAL_LIBRARY_DIRS})
target_link_libraries(local_search PUBLIC ${GDAL_LIBRARIES})
TARGET_LINK_LIBRARIES(local_search PUBLIC ${Boost_LIBRARIES})

TARGET_LINK_LIBRARIES(local_search PUBLIC GeographicLib)

add_executable(local_search local_search.cpp dataset_load.cpp)
TARGET_LINK_LIBRARIES(local_search_s2 PUBLIC s2::s2)


add_executable(s2_test s2_test.cpp)
target_link_libraries(s2_test PRIVATE s2::s2)