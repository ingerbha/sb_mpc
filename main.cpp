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
#include <ctime>
#include "iostream"
#include <stdio.h>

int main(){
	// Guidance parameters
	double u_d = 3;
	double psi_d = 0.0;

	// ASV state
//	double x_s = 0.0;
//	double y_s = 0.0;
//	double psi_s = 0.3;
//	double u_s = 3.0;
//	double v_s = 0.0;
//	double r_s = 0.0;
	Eigen::Matrix<double, 6, 1> asv_state;
//	asv_state <<  x_s, y_s, psi_s, u_s, v_s, r_s;
    asv_state << 0, 0, 0.255222, 0.871141, 0, 0;


	// Obstacle states
//	double x_o = 100.0;
//	double y_o1 = 0.0; double y_o2 = -50.0;
//	double psi_o = -3.13;
//	double u_o = 3.0;
//	double v_o = 0.0;
//	double A = 5;
//	double B = 4;
//	double C = 1;
//	double D = 2;
	Eigen::Matrix<double, 5, 9> obst_states;
//	obst_states << x_o, y_o1, psi_o, u_o, v_o, A, B, C, D,
//				   x_o, y_o2, psi_o, u_o, v_o, A, B, C, D;
 	obst_states << 110.359,  146.154,  4.71239,        3,        0,       10,       10,       10,       10,
				 -826.865,  714.124,        0,        0,        0,       10,       10,       10,       10,
				 3444.27, -1950.35,  5.72293,  13.2727,        0,       10,       10,       10,       10,
				 -614.694,  735.825,  1.24093,  3.29244,        0,       10,       10,       10,       10,
				 -522.324,  6325.05,  2.27242,        0,        0,       10,       10,       10,       10;


	double u_os, psi_os;
	simulationBasedMpc *sb_mpc = new simulationBasedMpc();

	sb_mpc->getBestControlOffset(u_os, psi_os, u_d, psi_d, asv_state, obst_states);

	std::cout << "u_ os : " << u_os << "\npsi_os : " << psi_os << std::endl;

	return 0;

//	// Guidance parameters
//	double u_d = 3;
//	double psi_d = 0.0;
//
//	// ASV state
//	double x_s = 0.0;
//	double y_s = 0.0;
//	double psi_s = 0.3;
//	double u_s = 3.0;
//	double v_s = 0.0;
//	double r_s = 0.0;
//	Eigen::Matrix<double, 6, 1> asv_state;
//	asv_state <<  x_s, y_s, psi_s, u_s, v_s, r_s;
//
//	// Obstacle states
//	double x_o = 100.0;
//	double y_o1 = 0.0; double y_o2 = -50.0;
//	double psi_o = -3.13;
//	double u_o = 3.0;
//	double v_o = 0.0;
//	double A = 5;
//	double B = 4;
//	double C = 1;
//	double D = 2;
//	Eigen::Matrix<double, 2, 9> obst_states;
//	obst_states << x_o, y_o1, psi_o, u_o, v_o, A, B, C, D,
//				   x_o, y_o2, psi_o, u_o, v_o, A, B, C, D;
//
//
//	double u_os, psi_os;
//	simulationBasedMpc *sb_mpc = new simulationBasedMpc();
//
//	sb_mpc->getBestControlOffset(u_os, psi_os, u_d, psi_d, asv_state, obst_states);
//
//	std::cout << "u_ os : " << u_os << "\npsi_os : " << psi_os << std::endl;
//
//	return 0;
};
