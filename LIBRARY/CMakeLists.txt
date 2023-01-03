add_library(XRaySimulator STATIC)

list(APPEND XRLIB_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/src/utilities.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/src/photons.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/src/scene.cpp
    )

list(APPEND XRLIB_INCLUDES
     ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

message(${CMAKE_CURRENT_SOURCE_DIR})

target_sources(XRaySimulator PUBLIC ${XRLIB_SOURCES})
target_include_directories(XRaySimulator PUBLIC ${XRLIB_INCLUDES})

target_link_libraries(XRaySimulator compiler_flags)