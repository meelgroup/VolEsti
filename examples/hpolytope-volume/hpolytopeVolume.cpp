// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Contributed and/or modified by Vaibhav Thakkar

// Licensed under GNU LGPL.3, see LICENCE file

#include "Eigen/Eigen"
#include <vector>
#include <string>
#include "cartesian_geom/cartesian_kernel.h"
#include "hpolytope.h"
#include "known_polytope_generators.h"

#include "random_walks/random_walks.hpp"

#include "volume_sequence_of_balls.hpp"
#include "volume_cooling_gaussians.hpp"
#include "volume_cooling_balls.hpp"

#include <iostream>
#include <fstream>
#include "misc.h"

typedef double NT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
typedef HPolytope <Point> HPOLYTOPE;
void calculateVolumes(HPOLYTOPE &HP, std::string algo, int walk_len = 0) {
	// Setup parameters for calculating volume
	int walk_len_act = 10 + HP.dimension()/10;

	if (walk_len > 0) {
		std::cout<<"c [volesti] Using walk length: "<<walk_len<< " calculated walk length: "<<walk_len_act<<"\n";
	}
	else {
		walk_len = walk_len_act;
		std::cout<<"c [volesti] Using walk length: "<<walk_len_act<<"\n";
	}

	NT e=0.1;

	NT volume = 0;
	if (algo == "coolingball") {
		std::cout<<"c [volesti] Using Cooling Balls method\n";
		volume = volume_cooling_balls<BallWalk, RNGType, HPOLYTOPE>(HP, e, walk_len).second;
	} else if (algo == "seqball") {
		std::cout<<"c [volesti] Using Sequence of Balls method\n";
		volume = volume_sequence_of_balls<BallWalk, RNGType>(HP, e, walk_len);
	} else if (algo == "coolinggauss") {
		std::cout<<"c [volesti] Using Cooling Gaussians method\n";
		volume = volume_cooling_gaussians<GaussianBallWalk, RNGType>(HP, e, walk_len);
	} else {
		std::cerr << "Unknown algorithm: " << algo << std::endl;
		return;
	}

	if(volume < 0)
		std::cout<<"c vol 0\n";
	else
		std::cout<<"c vol "<<volume<< std::endl;

}


int main(int argc, char* argv[]) {
	// Generating a 4-dimensional cube centered at origin
	// HPOLYTOPE HP1 = generate_cube<HPOLYTOPE>(4, false);
	// std::cout<<"Polytope HP1: \n";
	// HP1.print();
	// std::cout<<"\n";

	// std::cout<<"Volume of HP1: \n";
	// calculateVolumes(HP1);


	// Reading a polytope from ine file
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <input_file> [--algorithm <coolinggauss|coolingball|seqball>] [--walklen <length>]" << std::endl;
		return 1;
	}
	std::string fileName(argv[1]);
	std::string algorithm = "coolinggauss";
	int walk_len = 10;

	for (int i = 2; i < argc; ++i) {
		if (std::string(argv[i]) == "--algorithm" && i + 1 < argc) {
			algorithm = argv[++i];
		} else if (std::string(argv[i]) == "--walklen" && i + 1 < argc) {
			walk_len = std::stoi(argv[++i]);
		}
	}
	// std::cout<<"Reading input from file..."<<std::endl;
	std::ifstream inp;
	std::vector<std::vector<NT> > Pin;
	inp.open(fileName, std::ifstream::in);
	read_pointset(inp,Pin);

	HPOLYTOPE HP2(Pin);
	// std::cout<<"Polytope HP2: \n";
	// HP2.print();
	// std::cout<<"\n";

	// std::cout<<"Volume of HP2: \n";
	calculateVolumes(HP2, algorithm, walk_len);

	return 0;
}
