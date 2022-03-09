add_executable(test)

target_sources(test
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

target_include_directories(test
PRIVATE 
    ${PROJECT_SOURCE_DIR}/include/matoy
)