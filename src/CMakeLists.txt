cmake_minimum_required(VERSION 2.8)

# set gcc and g++ compiler
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# set the project name
project(HAPoG C)
enable_language(C)

# add compile options
set(CMAKE_C_FLAGS "-O3 -g")

# Add htslib
set(HTSLIB ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

# create a custom target called build_scintilla that is part of ALL
# and will run each time you type make 
add_custom_target(build_htslib ALL 
  COMMAND autoheader
  COMMAND autoconf
  COMMAND ./configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/htslib
  COMMAND ${CMAKE_MAKE_PROGRAM}
  COMMAND ${CMAKE_MAKE_PROGRAM} install
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib/
  COMMENT "Original htslib makefile target")

# now create an imported static target
add_library(htslib STATIC IMPORTED)
# Import target "scintilla" for configuration ""
set_property(TARGET htslib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(htslib PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${HTSLIB}")

# now you can use scintilla as if it were a regular cmake built target in your project
add_dependencies(htslib build_htslib)



# add the executable
add_executable(hapog polish_consensus.c alipile.c hash.c polished.c)

target_include_directories(hapog PRIVATE ${HAPoG_SOURCE_DIR}/htslib)
INCLUDE_DIRECTORIES(${HAPoG_SOURCE_DIR}/htslib/htslib)

target_link_libraries(hapog m)
target_link_libraries(hapog htslib)
target_link_libraries(hapog z)
target_link_libraries(hapog lzma)
target_link_libraries(hapog bz2)
