/*
 * main.cpp
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#include "obstacle.h"
#include "sb_mpc.h"
#include "ship_model.h"

#include "Eigen/Dense"
#include "iostream"

int main(){

	// Guidance parameters
	double u_d = 9.34019;
	double psi_d = 1.71925;

	// Offsets
	double u_os;
	double psi_os;

	simulationBasedMpc *sb_mpc = new simulationBasedMpc();

	Eigen::Matrix<double, 6, 1> asv_state;
	asv_state << 0, 0, -0.0959975, 6.88277, 0, 0;

	Eigen::Matrix<double, 5, 9> obst_states;
	obst_states << 110.359,  146.154,  4.71239,        3,        0,       10,       10,       10,       10,
				 -826.865,  714.124,        0,        0,        0,       10,       10,       10,       10,
				 3444.27, -1950.35,  5.72293,  13.2727,        0,       10,       10,       10,       10,
				 -614.694,  735.825,  1.24093,  3.29244,        0,       10,       10,       10,       10,
				 -522.324,  6325.05,  2.27242,        0,        0,       10,       10,       10,       10;

	sb_mpc->getBestControlOffset(u_os, psi_os, u_d, psi_d, asv_state, obst_states);

	std::cout << "u_ os : " << u_os << std::endl;
	std::cout << "psi_os : " << psi_os << std::endl;

	return 0;
};
