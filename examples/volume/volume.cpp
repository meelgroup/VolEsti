// VolEsti ( volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Conributed and/or modified by Suraj Choubey

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "known_polytope_generators.h"
#include "random_walks/random_walks.hpp"
#include <vector>

#include "misc.h"
#include "vpolytope.h"
#include "hpolytope.h"
#include <fstream>
#include <iostream>
#include <limits>

#include "volume_cooling_balls.hpp"
#include "volume_cooling_gaussians.hpp"

typedef double NT;
typedef Cartesian<NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef VPolytope<Point> VPOLYTOPE;
typedef HPolytope<Point> HPOLYTOPE;

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

NT calculateVolumesHP(HPOLYTOPE &HP, std::string algo, int walk_len = 0) {
  // Setup parameters for calculating volume

	walk_len = get_walk_len(HP.dimension(), walk_len);

  NT e = 0.1;
  NT volume = 0;


  if (algo == "coolingball") {
    std::cout << "c [volesti] Using Cooling Balls method\n";
    volume = volume_cooling_balls<BallWalk, RNGType, HPOLYTOPE>(HP, e, walk_len)
                 .second;

  } else if (algo == "coolinggauss") {
    std::cout << "c [volesti] Using Cooling Gaussians method\n";
    volume =
        volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, e, walk_len);
  } else {
    std::cerr << "Unknown algorithm: " << algo << std::endl;
    return -1;
  }

	return volume;

}

NT calculateVolumesVP(VPOLYTOPE &VP, std::string algo, int walk_len = 0) {
  // Setup parameters for calculating volume

	walk_len = get_walk_len(VP.dimension(), walk_len);

  NT e = 0.1;
  NT volume = 0;


  if (algo == "coolingball") {
    std::cout << "c [volesti] Using Cooling Balls method\n";
    volume = volume_cooling_balls<BallWalk, RNGType, VPOLYTOPE>(VP, e, walk_len)
                 .second;

  } else if (algo == "coolinggauss") {
    std::cout << "c [volesti] Using Cooling Gaussians method\n";
    volume =
        volume_cooling_gaussians<GaussianBallWalk, RNGType>(VP, e, walk_len);
  } else {
    std::cerr << "Unknown algorithm: " << algo << std::endl;
    return -1;
  }

	return volume;

}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input_file> [--algorithm "
                 "<coolinggauss|coolingball|seqball>] [--walklen <length>]"
              << std::endl;
    return 1;
  }
  std::string fileName(argv[1]);
  std::string algorithm = "coolinggauss";
  int walk_len = 10;
	NT volume = 0;

	bool vpolytope = false;

	for (int i = 2; i < argc; ++i) {
		if (std::string(argv[i]) == "--algorithm" && i + 1 < argc) {
			algorithm = argv[++i];
			std::cout << "c [volesti] Using algorithm: " << algorithm << "\n";
		} else if (std::string(argv[i]) == "--walklen" && i + 1 < argc) {
			walk_len = std::stoi(argv[++i]);
		} else if (std::string(argv[i]) == "--vpoly") {
			vpolytope = true;
			std::cout << "c [volesti] Expecting VPolytope\n";
		}
	}
  std::ifstream inp;
  std::vector<std::vector<NT>> Pin;
  inp.open(fileName, std::ifstream::in);
  read_pointset(inp, Pin);

	if (vpolytope) {
		VPOLYTOPE VP2(Pin);
		volume = calculateVolumesVP(VP2, algorithm, walk_len);
	} else {
		HPOLYTOPE HP2(Pin);
		volume = calculateVolumesHP(HP2, algorithm, walk_len);
	}

	  if (volume < 0)
    std::cout << "c vol 0\n";
  else
    std::cout << "c vol " << volume << std::endl;

  return 0;
}
