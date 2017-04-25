project(sswlib LANGUAGES CXX C)

# project paths
set(ssw_RootDir     ${PROJECT_SOURCE_DIR})
set(ssw_SourceDir   ${ssw_RootDir}/src)

if(NOT ssw_OutputDir)
    set(ssw_OutputDir ${ssw_RootDir})
endif()

set(ssw_LibDir  ${ssw_OutputDir}/lib)

set(ssw_SOURCE_FILES
        ${ssw_SourceDir}/ssw.c
        ${ssw_SourceDir}/ssw_cpp.cpp
        )

add_library(sswlib ${ssw_SOURCE_FILES})
