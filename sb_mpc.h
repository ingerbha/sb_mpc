/*
 * sb_mpc.h
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#ifndef SB_MPC_H_
#define SB_MPC_H_

#include <Eigen/Dense>
#include "ship_model.h"
#include "obstacle.h"

#include <vector>

class simulationBasedMpc
{
public:
	/// Constructor
	simulationBasedMpc();

	/// Destructor
	~simulationBasedMpc();

	void getBestControlOffset(double &u_d_best, double &psi_d_best, double u_d, double psi_d, Eigen::Matrix<double,6,1> asv_state, Eigen::Matrix<double,Eigen::Dynamic,9> obst_states);

private:

	void eulerIntegration(double u_d, double psi_d);

	double costFunction(double P_ca, double Chi_ca, int k);

	double deltaP(double P_ca);

	double deltaChi(double Chi_ca);

	void rot2d(double yaw, Eigen::Vector2d &res);

	// Simulation Parameters
	double T_;
	double DT_;

	// Tuning Parameters
	double P_;
	double Q_;
	double D_CLOSE_;
	double D_SAFE_;
	double K_COLL_;
	double PHI_AH_;
	double PHI_OT_;
	double PHI_HO_;
	double PHI_CR_;
	double KAPPA_;
	double K_P_;
	double K_CHI_;
	double K_DP_;
	double K_DCHI_SB_;
	double K_DCHI_P_;

	// State Variables
	double Chi_ca_last_;
	double P_ca_last_;

	double cost_;

	Eigen::Matrix<double,13,1> Chi_ca_;
	Eigen::Vector4d P_ca_;

	shipModel *asv;
	std::vector<obstacle*> obst_vect;

};

#endif /* SB_MPC_H_ */
