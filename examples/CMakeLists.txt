add_executable(bench bench.cpp)
target_link_libraries(bench PRIVATE tc)

add_executable(path path.cpp)
target_link_libraries(path PRIVATE tc)

add_executable(group group.cpp)
target_link_libraries(group PRIVATE tc)

add_executable(complex complex.cpp)
target_link_libraries(complex PRIVATE tc)
