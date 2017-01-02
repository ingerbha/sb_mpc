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

	/**
	 *  @brief Sets the offset pair (u_os_best, psi_os_best) with the lowst cost
	 *
	 * @param u_os_best The reference parameter to store the best speed offset.
	 * @param psi_os_best The reference parameter to store the best heading offset.
	 * @param u_d The nominal speed
	 * @param psi_d The nominal heading
	 * @param asv_state The state of the asv: x, y, psi, u, v, r.
	 * @param obst_states The states of the obstacles : x, y, u, v, A, B, C, D. (A, B, C, D - size from AIS)
	 */
	void getBestControlOffset(double &u_os_best, double &psi_os_best, double u_d, double psi_d, const Eigen::Matrix<double,6,1>& asv_state, const Eigen::Matrix<double,-1,9>& obst_states);

	/**
	 * @brief  Returns the simulation time (prediction horizon) [sec].
	 */
	double getT();
	/**
	 * @brief Returns the time step for the simulation [sec].
	 */
	double getDt();
	/**
	 * @brief Returns the wheight on time to evaluation instant
	 */
	double getP();
	/**
	 * @brief Returns the wheight on distance at evaluation instant
	 */
	double getQ();
	/**
	 * @brief Returns the distance where COLREGS are said to apply [m].
	 */
	double getDClose();
	/**
	 * @brief Returns the minimal distance which is considered as safe [m].
	 */
	double getDSafe();
	/**
	 * @brief Returns the collision cost
	 */
	double getKColl();
	/**
	 * @brief Returns the angle within which an obstacle is said to be ahead
	 * [deg].
	 */
	double getPhiAH();
	/**
	 * @brief Returns the angle outside of which an obstacle will be said to
	 * be overtaking, if the speed of the obstacle is larger than the ship's
	 * own speed
	 */
	double getPhiOT();
	/**
	 * @brief Returns the angle whitin which an obstacle is said to be head
	 * on [deg].
	 */
	double getPhiHO();
	/**
	 * @brief Returns the angle outside of which an obstacle is said to be
	 * crossing, if it is on the starboard side, heading towards the ship
	 * and not overtaking the ship.
	 */
	double getPhiCR();
	/**
	 * @brief Returns the cost of not complying with the COLREGS.
	 */
	double getKappa();
	/**
	 * @brief Returns the cost of deviating from the nominal speed.
	 */
	double getKP();
	/**
	 * @brief Returns the cost of changing the speed offset.
	 */
	double getKdP();
	/**
	 * @brief Returns the cost of deviating from the nominal heading
	 */
	double getKChi();
	/**
	 * @brief Returns the cost of changing the heading offset to starboard.
	 */
	double getKdChiSB();
	/**
	 * @brief Returns the cost of changing the heading offset to port.
	 */
	double getKdChiP();
	/**
	 * @brief Returns the possible offsets to the nominal heading [deg].
	 */
	Eigen::VectorXd getChiCA();
	/**
	 * @brief Returns the possible offsets to the nominal course, should be
	 * in the range [-1,1].
	 */
	Eigen::VectorXd getPCA();
	/**
	 * @brief Sets the prediction horizon [sec].
	 */
	void setT(double T);
	/**
	 * @brief Sets the simulation step time [sec].
	 */
	void setDt(double dt);

	void setP(double p);
	void setQ(double q);
	void setDClose(double d_close);
	void setDSafe(double d_safe);
	void setKColl(double k_coll);
	void setPhiAH(double phi_AH);
	void setPhiOT(double phi_OT);
	void setPhiHO(double phi_HO);
	void setPhiCR(double phi_CR);
	void setKappa(double kappa);
	void setKP(double K_P);
	void setKdP(double K_dP);
	void setKChi(double K_Chi);
	void setKdChiSB(double K_dChi_SB);
	void setKdChiP(double K_dChi_P);
	void setChiCA(Eigen::VectorXd Chi_ca);
	void setPCA(Eigen::VectorXd P_ca);


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

	Eigen::VectorXd Chi_ca_;
	Eigen::VectorXd P_ca_;

	shipModel *asv;
	std::vector<obstacle*> obst_vect;

};

#endif /* SB_MPC_H_ */
