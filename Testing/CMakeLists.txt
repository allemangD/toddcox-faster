cmake_policy(SET CMP0110 NEW)

add_executable(schlafli schlafli.cpp)
target_link_libraries(schlafli PRIVATE tc)

add_executable(special special.cpp)
target_link_libraries(special PRIVATE tc)

add_test(NAME "schalfi [6]" COMMAND schlafli "6" "" "12")
add_test(NAME "schalfi [6] [0]" COMMAND schlafli "6" "0" "6")
add_test(NAME "schalfi [6] [1]" COMMAND schlafli "6" "1" "6")

add_test(NAME "schalfi [5 3 3] []" COMMAND schlafli "5 3 3" "" "14400")
add_test(NAME "schalfi [5 3 3] [0 1]" COMMAND schlafli "5 3 3" "0 1" "1440")
add_test(NAME "schalfi [5 3 3] [1 2 3]" COMMAND schlafli "5 3 3" "1 2 3" "600")

add_test(NAME "special E6" COMMAND special "E6" "" "51840")
add_test(NAME "special E7" COMMAND special "E7" "" "2903040")
#add_test(NAME "special E8" COMMAND special "E8" "" "696729600")  # too big.
add_test(NAME "special B6" COMMAND special "B6" "" "46080")
add_test(NAME "special B7" COMMAND special "B7" "" "645120")
add_test(NAME "special B8" COMMAND special "B8" "" "10321920")
