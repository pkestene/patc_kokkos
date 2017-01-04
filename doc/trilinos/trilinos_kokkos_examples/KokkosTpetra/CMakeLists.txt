#
# This CMakeLists.txt is slightly adapted from
# https://github.com/trilinos/Trilinos_tutorial/wiki/CMakeLists_Kokkos03
#

INCLUDE_DIRECTORIES ( ./ ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

MACRO(TUTO_ADD_TEST TEST_NAME)
  ADD_EXECUTABLE(KokkosTpetra_${TEST_NAME} KokkosTpetra_${TEST_NAME}.cpp)
  TARGET_LINK_LIBRARIES(KokkosTpetra_${TEST_NAME}
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES})
ENDMACRO()

TUTO_ADD_TEST(01_init)
TUTO_ADD_TEST(02_problem)
TUTO_ADD_TEST(03_fill)
TUTO_ADD_TEST(04_tpetra)
TUTO_ADD_TEST(05_solve)
