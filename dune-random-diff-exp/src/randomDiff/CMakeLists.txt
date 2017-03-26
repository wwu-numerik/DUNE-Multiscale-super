add_subdirectory("Solvers")
add_subdirectory("LocalOperators")

file(GLOB LOCAL_SRC "*.hh" "*.cc")
set(SRC ${SRC} ${LOCAL_SRC} PARENT_SCOPE)
