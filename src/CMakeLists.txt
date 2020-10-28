cmake_minimum_required(VERSION 3.0)
include(ExternalProject)

# set gcc and g++ compiler
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# set the project nameExternalProject_Add
project(HAPoG C)
enable_language(C)

# Add htslib
set(HTSLIB_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/htslib)
set(HTSLIB_INSTALL make install prefix=${HTSLIB_PREFIX}/src/htslib/)

ExternalProject_Add(htslib
  PREFIX ${HTSLIB_PREFIX}
  GIT_REPOSITORY "https://github.com/samtools/htslib"
  GIT_TAG "1.9"
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND autoheader && autoconf && ${HTSLIB_PREFIX}/src/htslib/configure --prefix=${HTSLIB_PREFIX}/src/htslib/
  BUILD_COMMAND make
  INSTALL_COMMAND ${HTSLIB_INSTALL}
  LOG_DOWNLOAD 0
  LOG_UPDATE 0
  LOG_CONFIGURE 0
  LOG_BUILD 0
  LOG_TEST 0
  LOG_INSTALL 0
)

ExternalProject_Get_Property(htslib install_dir)
include_directories(${HTSLIB_PREFIX}/src/htslib/include)

# Add the executable
set(CMAKE_C_FLAGS "-O3 -g")
add_executable(hapog polish_consensus.c alipile.c hash.c polished.c)

# Add htslib include dir
include_directories(${install_dir}/include)
add_dependencies(hapog htslib)

# Links
target_link_libraries(hapog m)
target_link_libraries(hapog ${HTSLIB_PREFIX}/src/htslib/lib/libhts.so.1.9)
target_link_libraries(hapog z)
target_link_libraries(hapog lzma)
target_link_libraries(hapog bz2)
