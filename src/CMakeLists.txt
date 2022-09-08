cmake_minimum_required(VERSION 3.0)
include(ExternalProject)

# Set gcc and g++ compiler
set(CMAKE_C_COMPILER clang)
#set(CMAKE_CXX_COMPILER g++)

# Set the project name
project(HAPoG C)
enable_language(C)

# Add htslib
if(NOT "$ENV{HTSLIB_ROOT}" STREQUAL "")
    message("\nWill attempt to find HTSlib at the provided path $ENV{HTSLIB_ROOT}")
    include(find_htslib.cmake)
    if("${htslib_found}" STREQUAL "TRUE")
        message("Found htslib at ${htslib_INCLUDE_DIR}")
        include_directories("${htslib_INCLUDE_DIR}/include")
        link_directories(${htslib_INCLUDE_DIR}/lib)
    else()
        message(FATAL_ERROR "Could not find the required files at the given $ENV{HTSLIB_ROOT} path.")
    endif()
else()
    message("\nNo htslib path was provided, trying to find it on system")
    include(find_htslib.cmake)
    if("${htslib_found}" STREQUAL "TRUE")
        message("Found htslib at ${htslib_INCLUDE_DIR}")
        include_directories("${htslib_INCLUDE_DIR}")
        link_directories(${htslib_LIBRARY})
    else()
        message("htslib was not found on system, attempting to download and compile it")
        set(HTSLIB_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/htslib)
        set(HTSLIB_INSTALL make install prefix=${HTSLIB_PREFIX}/src/htslib/)

        ExternalProject_Add(htslib
            PREFIX ${HTSLIB_PREFIX}
            GIT_REPOSITORY "https://github.com/samtools/htslib"
            GIT_TAG "1.11"
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
        include_directories(${HTSLIB_PREFIX}/src/htslib)
        link_directories(${HTSLIB_PREFIX}/src/htslib)
    endif()
endif()


# Add the executable
set(CMAKE_C_FLAGS "-O3 -g")
add_executable(hapog polish_consensus.c alipile.c hash.c polished.c)

# Links
target_link_libraries(hapog pthread)
target_link_libraries(hapog curl)
target_link_libraries(hapog crypto)
target_link_libraries(hapog m)
target_link_libraries(hapog hts)
target_link_libraries(hapog z)
target_link_libraries(hapog lzma)
target_link_libraries(hapog bz2)
