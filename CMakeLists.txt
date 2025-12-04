cmake_minimum_required(VERSION 3.6)

project(loga)

set(Boost_NO_BOOST_CMAKE        ON  CACHE BOOL "Use FindBoost module, not BoostConfig")
set(Boost_USE_STATIC_LIBS       OFF)
set(Boost_USE_MULTITHREADED     ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

if (MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Armadillo CONFIG REQUIRED)
find_package(cereal REQUIRED)
find_package(igraph REQUIRED)

if(Armadillo_FOUND AND NOT TARGET Armadillo::Armadillo)
  add_library(Armadillo::Armadillo INTERFACE IMPORTED)
  set_target_properties(
    Armadillo::Armadillo
    PROPERTIES
      INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}"
  )
endif()

include_directories(includes)

set(SOURCES
  sources/index.cpp
  sources/face.cpp
  sources/collection.cpp
  sources/alignment.cpp
  sources/multi_alignment.cpp
  sources/segment.cpp
  sources/graph.cpp
  sources/path.cpp
  sources/zone.cpp
  sources/distance.cpp
  sources/group.cpp
  sources/bio_alignment.cpp
  sources/tokenized_collection.cpp
  sources/token.cpp
  sources/tokenized_group.cpp
  sources/tokenized_alignment.cpp
  sources/tokenized_distance.cpp
  sources/tokenized_multi_alignment.cpp
  sources/knn_cluster.cpp
  sources/community.cpp
  sources/automata.cpp
  sources/outliers.cpp
  sources/pattern_sequence.cpp
)

set(INCLUDES
  includes/loga/fwd.h
  includes/loga/face.h
  includes/loga/index.h
  includes/loga/collection.h
  includes/loga/alignment.h
  includes/loga/multi_alignment.h
  includes/loga/segment.h
  includes/loga/graph.h
  includes/loga/path.h
  includes/loga/zone.h
  includes/loga/package.h
  includes/loga/distance.h
  includes/loga/group.h
  includes/loga/bio_alignment.h
  includes/loga/token.h
  includes/loga/tokenized_collection.h
  includes/loga/tokenized_group.h
  includes/loga/tokenized_alignment.h
  includes/loga/tokenized_distance.h
  includes/loga/tokenized_multi_alignment.h
  includes/loga/colors.h
  includes/loga/knn_cluster.h
  includes/loga/community.h
  includes/loga/automata.h
  includes/loga/pattern_sequence.h
  includes/loga/outliers.h
)

add_library(loga_lib SHARED ${SOURCES} ${INCLUDES})
set_property(TARGET loga_lib PROPERTY CXX_STANDARD 20)
target_link_libraries(loga_lib PUBLIC Boost::headers igraph::igraph)
target_link_libraries(loga_lib PRIVATE cereal::cereal Armadillo::Armadillo)
set_target_properties(loga_lib PROPERTIES OUTPUT_NAME "loga")
if (WIN32)
    set_target_properties(loga_lib PROPERTIES ARCHIVE_OUTPUT_NAME "loga_lib")
endif()

add_executable(loga_exe sources/main.cpp)
set_property(TARGET loga_exe PROPERTY CXX_STANDARD 20)
target_link_libraries(loga_exe PRIVATE Boost::headers Boost::program_options Armadillo::Armadillo loga_lib)
set_property(TARGET loga_exe PROPERTY OUTPUT_NAME loga)

if(APPLE)
  set_target_properties(loga_exe PROPERTIES
      INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}"
      BUILD_RPATH   "@loader_path/../${CMAKE_INSTALL_LIBDIR}"
  )
elseif(UNIX)
  set_target_properties(loga_exe PROPERTIES
      INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
      BUILD_RPATH   "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
  )
endif()

include(GNUInstallDirs)

install(TARGETS loga_lib loga_exe
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES ${INCLUDES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/loga
)
