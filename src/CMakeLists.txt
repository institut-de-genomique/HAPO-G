cmake_minimum_required(VERSION 2.8)

# set gcc and g++ compiler
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# set the project name
project(yap C)
enable_language(C)

# add compile options
set(CMAKE_C_FLAGS "-O3 -g")

# Add htslib
link_directories(${yap_SOURCE_DIR}/../samtools/htslib-1.9)

# add the executable
add_executable(yap polish_consensus.c alipile.c hash.c polished.c)

target_include_directories(yap PRIVATE /env/ig/soft/ig/gcc-8.3.0/el6-x86_64-generic/include/c++/8 /env/ig/soft/ig/gcc-8.3.0/el6-x86_64-generic/include/c++/8/x86_64-redhat-linux-gnu)
target_include_directories(yap PRIVATE ${yap_SOURCE_DIR}/../samtools/htslib-1.9)

target_link_libraries(yap m)
target_link_libraries(yap hts)
target_link_libraries(yap z)
target_link_libraries(yap lzma)
target_link_libraries(yap bz2)
