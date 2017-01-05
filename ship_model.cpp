/*
 * ship_model.cpp
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#include "ship_model.h"

shipModel::shipModel(double T, double dt)
: n_samp_(T/dt)
{
	x.resize(n_samp_);
	y.resize(n_samp_);
	psi.resize(n_samp_);
	u.resize(n_samp_);
	v.resize(n_samp_);
	r.resize(n_samp_);

	tau = Eigen::Vector3d::Zero();

	// Simulation parameters
	T_ = T;
	DT_ = dt;

	// Model parameters
	rudder_d = 4.0; // distance from rudder to CG
	A_ = 5; // [m]  in reality the length is 14,5 m.
	B_ = 5; // [m]
	C_ = 1.5; // [m]
	D_ = 1.5; // [m]
	calculate_position_offsets();
	M = 3980.0; // [kg]
	I_z = 19703.0; // [kg/m2]

	// Added M terms
	X_udot = 0.0;
	Y_vdot = 0.0;
	Y_rdot = 0.0;
	N_vdot = 0.0;
	N_rdot = 0.0;

	// Linear damping terms [X_u, Y_v, Y_r, N_v, N_r]
	X_u	= -50.0;
	Y_v = -200.0;
	Y_r = 0.0;
	N_v = 0.0;
	N_r = -1281;//-3224.0;

	// Nonlinear damping terms [X_|u|u, Y_|v|v, N_|r|r, X_uuu, Y_vvv, N_rrr]
	X_uu = -135.0;
	Y_vv = -2000.0;
	N_rr = 0.0;
	X_uuu = 0.0;
	Y_vvv = 0.0;
	N_rrr = -3224.0;

	Eigen::Matrix3d Mtot;
	Mtot << M - X_udot, 0, 0,
	        0, M-Y_vdot, -Y_rdot,
	        0, -Y_rdot, I_z-N_rdot;
	Minv = Mtot.inverse();

	//Force limits
	Fx_min = -6550.0;
	Fx_max = 13100.0;
	Fy_min = -645.0;
	Fy_max = 645.0;

	// Controller parameters
	Kp_u = 1.0;
	Kp_psi = 5.0;
	Kd_psi = 1.0;
	Kp_r = 8.0;

}


shipModel::~shipModel(){
}


Eigen::VectorXd shipModel::getX(){
	return x;
}

Eigen::VectorXd shipModel::getY(){
	return y;
}

Eigen::VectorXd shipModel::getPsi(){
	return psi;
}

Eigen::VectorXd shipModel::getU(){
	return u;
}

Eigen::VectorXd shipModel::getV(){
	return v;
}

Eigen::VectorXd shipModel::getR(){
	return r;
}

double shipModel::getA(){
	return A_;
}

double shipModel::getB(){
	return B_;
}

double shipModel::getC(){
	return C_;
}

double shipModel::getD(){
	return D_;
}

void shipModel::setB(double B){
	B_ = B;
}

void shipModel::calculate_position_offsets(){
	os_x = A_-B_;
	os_y = D_-C_;
}

void shipModel::eulersMethod(const Eigen::Matrix<double,6,1>& state, double u_d, double psi_d)
{

	psi(0) = normalize_angle(state(2));
	x(0) = state(0) + os_x*cos(psi(0)) - os_y*sin(psi(0));
	y(0) = state(1) + os_x*sin(psi(0)) + os_y*cos(psi(0));
	u(0) = state(3);
	v(0) = state(4);
	r(0) = state(5);


	double t = 0;
	Eigen::Vector3d temp;
	double r11, r12, r21, r22; // rotation matrix elements

	for (int i = 0; i < n_samp_-1; i++){

		t += DT_;

		psi_d = normalize_angle_diff(psi_d, psi(i));

		r11 = cos(psi(i));
		r12 = -sin(psi(i));
		r21 = sin(psi(i));
		r22 = cos(psi(i));

		// Calculate coriolis and dampening matrices according to Fossen, 2011 or Stenersen, 2014.
		Cvv(0) = (-M*v(i) + Y_vdot*v(i) + Y_rdot*r(i)) * r(i);
		Cvv(1) = ( M*u(i) - X_udot*u(i)) * r(i);
		Cvv(2) = (( M*v(i) - Y_vdot*v(i) - Y_rdot*r(i) ) * u(i) +
		            ( -M*u(i) + X_udot*u(i)) * v(i));

		Dvv(0) = - (X_u + X_uu*fabs(u(i)) + X_uuu*u(i)*u(i)) * u(i);
		Dvv(1) = - ((Y_v*v(i) + Y_r*r(i)) +
		              (Y_vv*fabs(v(i))*v(i) + Y_vvv*v(i)*v(i)*v(i)));
		Dvv(2) = - ((N_v*v(i) + N_r*r(i)) +
		              (N_rr*fabs(r(i))*r(i) + N_rrr*r(i)*r(i)*r(i)));

		this->updateCtrlInput(u_d, psi_d, i);

		// Integrate system

		x(i+1) = x(i) + DT_*(r11*u(i) + r12*v(i));
		y(i+1) = y(i) + DT_*(r21*u(i) + r22*v(i));
		psi(i+1) = psi(i) + DT_*r(i);
		temp = Minv * (tau - Cvv - Dvv);
		u(i+1) = u(i) + DT_*temp(0);
		v(i+1) = v(i) + DT_*temp(1);
		r(i+1) = r(i) + DT_*temp(2);

		// Keep yaw within [-PI,PI)
		psi(i+1) = normalize_angle(psi(i+1));

	}
}


void shipModel::updateCtrlInput(double u_d, double psi_d, int i){
	double Fx = Cvv[0] + Dvv[0] + Kp_u*M*(u_d - u(i));
	double Fy = 0.0;

    Fy = (Kp_psi * I_z ) * ((psi_d - psi(i)) - Kd_psi*r(i));
    Fy *= 1.0 / rudder_d;

	// Saturate
	if (Fx < Fx_min)
	  Fx = Fx_min;
	if (Fx > Fx_max)
	  Fx = Fx_max;

	if (Fy < Fy_min)
	  Fy = Fy_min;
	if (Fy > Fy_max)
	  Fy = Fy_max;

	tau[0] = Fx;
	tau[1] = Fy;
	tau[2] = rudder_d * Fy;
}


double shipModel::normalize_angle(double angle){

	if( isinf(angle)) return angle;

	while(angle <= -M_PI) angle += 2*M_PI;

	while (angle > M_PI) angle -= 2*M_PI;

	return angle;
}


double shipModel::normalize_angle_diff(double angle, double angle_ref){
	double new_angle;
	double diff = angle_ref - angle;

	if (isinf(angle) || isinf(angle_ref)) return angle;

	// Get angle within 2*PI of angle_ref
	if (diff > 0){
		new_angle = angle +(diff - fmod(diff, 2*M_PI));
	}else{
		new_angle = angle + (diff + fmod(-diff, 2*M_PI));
	}

	// Get angle on side closest to angle_ref
	diff = angle_ref - new_angle;
	if (diff > M_PI){
		new_angle += 2*M_PI;
	}else if (diff < -M_PI){
		new_angle -= 2*M_PI;
	}
	return new_angle;
}


