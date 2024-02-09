#TODO: this whole thing needs to be made properly, this is just a tmp for development
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

set(libname sirius_jl)
#add_dependencies(${libname} sirius)


#TODO: hardcoded PATH, this needs to change
list(APPEND CMAKE_PREFIX_PATH /home/bussya/.julia/artifacts/b28969bfe0d4c551302f58ecbf5a91930aa66489)
find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
#TODO: no idea what this does
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")                                                   
message("Lib is in : ${CMAKE_INSTALL_RPATH}")


add_library(${libname} SHARED sirius_jl.cpp)
add_dependencies(${libname} sirius)

target_link_libraries(${libname} sirius)
target_link_libraries(${libname} JlCxx::cxxwrap_julia)

install(TARGETS ${libname} LIBRARY DESTINATION lib/julia_wrapper)
