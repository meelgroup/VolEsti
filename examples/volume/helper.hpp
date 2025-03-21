#ifndef VOLESTI_HELPER_HPP
#define VOLESTI_HELPER_HPP

#define set_int cdd_set_int
extern "C" {
#include "cdd/setoper.h"
}
#undef set_int
#include "cdd/cdd.h"

int get_walk_len(int dim, int walk_len) {
  int walk_len_act = 10 + dim / 10;

  if (walk_len > 0) {
    std::cout << "c [volesti] Using walk length: " << walk_len
              << " calculated walk length: " << walk_len_act << "\n";
  } else {
    walk_len = walk_len_act;
    std::cout << "c [volesti] Using walk length: " << walk_len_act << "\n";
  }

  return walk_len;
}

void check_polytope_consistency(dd_MatrixPtr M) {
  if (!M || M->rowsize == 0 || M->colsize == 0) {
    throw std::runtime_error("Empty or null matrix");
  }
  for (int i = 0; i < M->rowsize; i++) {
    for (int j = 0; j < M->colsize; j++) {
      double val = dd_get_d(M->matrix[i][j]);
      if (!std::isfinite(val)) {
        throw std::runtime_error("Invalid entry in matrix");
      }
    }
  }

  // std::cout << "Matrix rows: " << M->rowsize
  //       << ", cols: " << M->colsize << std::endl;
  // for (int i = 0; i < M->rowsize; i++) {
  //   for (int j = 0; j < M->colsize; j++) {
  //     double val = dd_get_d(M->matrix[i][j]);
  //     std::cout << "M(" << i << "," << j << ") = " << val << std::endl;
  //   }
  // }

  // Uncomment if you wish to run an LP check:
  // dd_ErrorType errLP, lperr;
  // dd_LPPtr lp = dd_Matrix2LP(M, &errLP);
  // if (errLP != dd_NoError) throw std::runtime_error("Cannot create LP");
  // dd_boolean solveRet = dd_LPSolve(lp, dd_DualSimplex, &lperr);
  // if (lperr != dd_NoError || !solveRet) throw std::runtime_error("LP error");
  // dd_ErrorType copyErr;
  // dd_LPSolutionPtr lpSolution = dd_CopyLPSolution(lp, &copyErr);
  // if (copyErr != dd_NoError) throw std::runtime_error("LP copy error");
  // dd_FreeLPSolution(lpSolution);
  // dd_FreeLPData(lp);
}

#endif