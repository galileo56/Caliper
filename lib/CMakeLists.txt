FILE(GLOB libSources *.F90)
add_library(libCaliper STATIC ${libSources})
install(TARGETS libCaliper DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../lib)