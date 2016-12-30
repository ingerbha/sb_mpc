/*
 * ship_model.h
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#ifndef SHIP_MODEL_H_
#define SHIP_MODEL_H_

#include <vector>
#include <cmath>
#include <Eigen/Dense>

class shipModel
{
	public:

	/// Constructor
	shipModel(double T, double dt);

	/// Destructor
	~shipModel();

	void eulersMethod(const Eigen::Matrix<double,6,1>& state, double u_d, double psi_d);

	Eigen::VectorXd getX();
	Eigen::VectorXd getY();
	Eigen::VectorXd getPsi();
	Eigen::VectorXd getU();
	Eigen::VectorXd getV();
	Eigen::VectorXd getR();

	double getA();
	double getB();
	double getC();
	double getD();

	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd psi;
	Eigen::VectorXd u;
	Eigen::VectorXd v;
	Eigen::VectorXd r;

	double A, B, C, D;

	private:

	/// Assures that the numerical difference is at most PI
	double normalize_angle_diff(double angle, double angle_ref);

	/// Assures that angle is between [-PI, PI)
	double normalize_angle(double angle);

	void updateCtrlInput(double u_d, double psi_d, int i);

	Eigen::Vector3d tau;
	Eigen::Matrix3d Minv;
	Eigen::Vector3d Cvv;
	Eigen::Vector3d Dvv;

	// Model Parameters
	double rudder_d;
	double M; 	// [kg]
	double I_z; // [kg/m2]

	// Added mass terms
	double X_udot;
	double Y_vdot;
	double Y_rdot;
	double N_vdot;
	double N_rdot;

	// Linear damping terms [X_u, Y_v, Y_r, N_v, N_r]
	double X_u;
	double Y_v;
	double Y_r;
	double N_v;
	double N_r;

	// Nonlinear damping terms [X_|u|u, Y_|v|v, N_|r|r, X_uuu, Y_vvv, N_rrr]
	double X_uu;
	double Y_vv;
	double N_rr;
	double X_uuu;
	double Y_vvv;
	double N_rrr;

	//Force limits
	double Fx_min;
	double Fx_max;
	double Fy_min;
	double Fy_max;

	// Simulation parameters
	double DT_;
	double T_;
	const int n_samp_;

	// Controller parameters
	double Kp_u;
	double Kp_psi;
	double Kd_psi;
	double Kp_r;

};

#endif /* SHIP_MODEL_H_ */
