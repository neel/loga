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
if(DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg.cmake")
  find_package(Armadillo CONFIG REQUIRED)
else()
  find_package(Armadillo REQUIRED)
endif()
find_package(cereal REQUIRED)
find_package(igraph REQUIRED)

if(NOT TARGET Armadillo::Armadillo)
  if(TARGET armadillo)  # vcpkg / config package target
    add_library(Armadillo::Armadillo INTERFACE IMPORTED)
    set_property(TARGET Armadillo::Armadillo
      PROPERTY INTERFACE_LINK_LIBRARIES armadillo)
  elseif(DEFINED ARMADILLO_LIBRARIES OR DEFINED ARMADILLO_INCLUDE_DIRS) # FindArmadillo.cmake style
    add_library(Armadillo::Armadillo INTERFACE IMPORTED)
    set_target_properties(Armadillo::Armadillo PROPERTIES
      INTERFACE_LINK_LIBRARIES     "${ARMADILLO_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")
  else()
    message(FATAL_ERROR "Cannot construct Armadillo::Armadillo: neither 'armadillo' target nor ARMADILLO_LIBRARIES/INCLUDE_DIRS are set")
  endif()
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

include(GNUInstallDirs)

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

install(TARGETS loga_lib loga_exe
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set(LOGA_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
set(LOGA_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")

message(STATUS "LOGA_INSTALL_BINDIR: ${LOGA_INSTALL_BINDIR}")
message(STATUS "LOGA_INSTALL_LIBDIR: ${LOGA_INSTALL_LIBDIR}")

file(GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/loga_install_dirs_$<CONFIG>.cmake"
  CONTENT
    "set(LOGA_INSTALL_BINDIR \"${LOGA_INSTALL_BINDIR}\")\nset(LOGA_INSTALL_LIBDIR \"${LOGA_INSTALL_LIBDIR}\")\n"
)

file(GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/loga_runtime_deps_$<CONFIG>.cmake"
  CONTENT [==[
    include("${CMAKE_CURRENT_LIST_DIR}/loga_install_dirs_${CMAKE_INSTALL_CONFIG_NAME}.cmake")

    set(_exe "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${LOGA_INSTALL_BINDIR}/$<TARGET_FILE_NAME:loga_exe>")
    set(_lib "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${LOGA_INSTALL_LIBDIR}/$<TARGET_FILE_NAME:loga_lib>")

    get_filename_component(_exe "${_exe}" ABSOLUTE)
    get_filename_component(_lib "${_lib}" ABSOLUTE)

    message(STATUS "GET_RUNTIME_DEPENDENCIES scanning:")
    message(STATUS "  exe: ${_exe}")
    message(STATUS "  lib: ${_lib}")

    file(GET_RUNTIME_DEPENDENCIES
      EXECUTABLES "${_exe}"
      LIBRARIES   "${_lib}"
      RESOLVED_DEPENDENCIES_VAR   LOGA_RESOLVED
      UNRESOLVED_DEPENDENCIES_VAR LOGA_UNRESOLVED
      POST_EXCLUDE_REGEXES
         ".*/libc\\.so(\\..*)?$"
         ".*/libm\\.so(\\..*)?$"
         ".*/libpthread\\.so(\\..*)?$"
         ".*/libdl\\.so(\\..*)?$"
         ".*/libgcc_s\\.so(\\..*)?$"
         ".*/libstdc\\+\\+\\.so(\\..*)?$"
         ".*/librt\\.so(\\..*)?$"
         ".*/libutil\\.so(\\..*)?$"
         ".*/ld-linux-.*\\.so(\\..*)?$"
         ".*/libz\\.so(\\..*)?$"
         ".*/libxml2\\.so(\\..*)?$"
         ".*/libzstd\\.so(\\..*)?$"
         ".*/libudev\\.so(\\..*)?$"
         ".*/libevent.*\\.so(\\..*)?$"
    )

    message(STATUS "Loga runtime deps (resolved):")
    foreach(d IN LISTS LOGA_RESOLVED)
      message(STATUS "  ${d}")
    endforeach()

    if(LOGA_UNRESOLVED)
      message(WARNING "Loga runtime deps (unresolved): ${LOGA_UNRESOLVED}")
    endif()
]==]
)

file(GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/loga_runtime_deps_dispatch.cmake"
  CONTENT [==[
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME OR CMAKE_INSTALL_CONFIG_NAME STREQUAL "")
  set(CMAKE_INSTALL_CONFIG_NAME "Release")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/loga_runtime_deps_${CMAKE_INSTALL_CONFIG_NAME}.cmake")
]==]
)

# 4) Run the dispatcher at install time (no install(CODE))
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/loga_runtime_deps_dispatch.cmake")

install(FILES ${INCLUDES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/loga
)
