set(BINARY ${CMAKE_PROJECT_NAME})

add_executable(
        ${BINARY} main.cpp
        PDAG.cpp
        PDAG.h
        ScorerInterface.h
        Operators.h
        Operators.cpp
        set_ops.h
        set_ops.cpp
        EdgeQueueSet.cpp
        EdgeQueueSet.h
        Search.cpp
        Search.h
        
        BICScoreloader.cpp
        BICScorer.cpp
        BICScorer.h
        # BICScorer.cpp / BICScorer.h are NOT shipped here: download them from opensource
        # you may drop them into this directory, complete BICScoreloader.cpp, and re-enable the matrix-library
        # see the README section "BIC scorer"

        CircularBuffer.h
        utils.h
)

target_link_libraries(${BINARY} PRIVATE cnpy)
target_link_libraries(${BINARY} PRIVATE spdlog::spdlog)



find_package(Boost REQUIRED)
target_link_libraries(${BINARY} PRIVATE Boost::boost)

find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(${BINARY} PRIVATE ${BLAS_LIBRARIES})

# ---------------------------------------------------------------------------
# The BIC scorer sources you add (BICScorer.{h,cpp}) additionally require a
# matrix library (e.g., Eigen) that is intentionally NOT
# configured here. After placing the scorer sources and completing
# BICScoreloader.cpp, make that library available to this build YOURSELF --
# e.g. add its find_package()/link below, or pass its include path at configure
# time. See the README section "BIC scorer".
# ---------------------------------------------------------------------------