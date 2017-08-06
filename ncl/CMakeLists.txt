file(GLOB ncl_INC "*.h")
file(GLOB ncl_SRC "*.cpp")

include_directories(${CMAKE_SOURCE_DIR})
add_library(ncl_shared SHARED ${ncl_SRC})
add_library(ncl_static STATIC ${ncl_SRC})

set_target_properties(ncl_shared PROPERTIES OUTPUT_NAME ncl)
set_target_properties(ncl_static PROPERTIES OUTPUT_NAME ncl)

install(TARGETS ncl_shared DESTINATION "lib/ncl")
install(TARGETS ncl_static DESTINATION "lib/ncl")
install(FILES ${ncl_INC} DESTINATION "include/ncl")

