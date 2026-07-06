// BICScoreloader.cpp
// =============================================================================
//  >>> YOU MUST COMPLETE THIS FILE BEFORE BUILDING. <<<
//
//  This is the only bridge between the input .npy file and the Gaussian BIC
//  scorer. Its implementation relies on the BIC scorer from XGES, which we do
//  not redistribute. Everything you need is available, unmodified, from XGES:
//
//      https://github.com/ANazaret/XGES   (directory: src-cpp/)
//
//  ----------------------------------------------------------------------------
//  Steps to complete this file:
//
//   1. Download `BICScorer.h` and `BICScorer.cpp` from XGES (src-cpp/) and drop
//      them into this directory (breaking-bad/src/), unmodified.
//
//   2. At the top of this file, add the same includes and the small `.npy`
//      loader that XGES uses in `src-cpp/main.cpp`:
//         - `#include "BICScorer.h"` and `#include "cnpy/cnpy.h"`
//         - the row-major matrix typedef and the `load_npy<T>()` helper
//           (copy them verbatim from XGES `src-cpp/main.cpp` — do NOT rewrite).
//
//   3. Fill in the body of `make_bic_scorer_from_npy()` below so that it:
//         - loads `npy_path` into a matrix with `load_npy<double>()`,
//         - sets `n_variables` = number of columns, `n_samples` = number of rows,
//         - returns `std::make_unique<BICScorer>(matrix, alpha)`.
//
//  Once BICScorer.{h,cpp} are in place and the body below is filled in, the
//  project builds and runs normally. See the README section "BIC scorer".
// =============================================================================

#include <memory>
#include <stdexcept>
#include <string>

#include "ScorerInterface.h"

std::unique_ptr<ScorerInterface>
make_bic_scorer_from_npy(const std::string &npy_path, double alpha,
                         int &n_variables, int &n_samples) {
    // TODO: implement following steps 2-3 in the header comment above.
    (void) npy_path;
    (void) alpha;
    (void) n_variables;
    (void) n_samples;
    throw std::runtime_error(
            "make_bic_scorer_from_npy() is not implemented. Complete "
            "src/BICScoreloader.cpp (see the instructions at the top of the file "
            "and the README section \"BIC scorer\").");
}
