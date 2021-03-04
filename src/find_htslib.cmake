# - Try to find htslib
# Once done, this will define
#
#  htslib_found - system has htslib
#  htslib_INCLUDE_DIR - the htslib include directories
#  htslib_LIBRARY - link these to use htslib

set(HTSLIB_SEARCH_DIRS
        ${HTSLIB_SEARCH_DIRS}
        $ENV{HTLSIB_ROOT}
        )

set(_htslib_ver_path "htslib-${htslib_FIND_VERSION}")
include(LibFindMacros.cmake)

# Dependencies
#libfind_package(HTSlib)

# Include dir
find_path(htslib_INCLUDE_DIR
        NAMES include/htslib/sam.h
        PATHS ${HTSLIB_SEARCH_DIRS}
        HINTS ENV HTSLIB_ROOT
        )

# Finally the library itself
find_library(htslib_LIBRARY
        NAMES hts libhts.a hts.a
        PATHS ${htslib_INCLUDE_DIR} ${HTSLIB_SEARCH_DIRS}
        NO_DEFAULT_PATH
        PATH_SUFFIXES lib lib64 ${_htslib_ver_path}
        HINTS ENV HTSLIB_ROOT
        )

if(NOT "${htslib_INCLUDE_DIR}" STREQUAL "htslib_INCLUDE_DIR-NOTFOUND" AND NOT "${htslib_LIBRARY}" STREQUAL "htslib_LIBRARY-NOTFOUND")
        set(htslib_found "TRUE")
else()
        set(htslib_found "FALSE")
endif()
message(${htslib_INCLUDE_DIR})
message(${htslib_LIBRARY})
message("htslib_found ${htslib_found}\n")
