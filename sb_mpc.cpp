/*
 * sb_mpc.cpp
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#include "sb_mpc.h"
#include "ship_model.h"

#include <vector>

static const double DEG2RAD = M_PI/180.0f;

simulationBasedMpc::simulationBasedMpc(){
	T_ = 50.0;//100.0;
	DT_ = 1.0;//0.05;

	P_ = 1.0;
	Q_ = 4.0;
	D_CLOSE_ = 100.0;
	D_SAFE_ = 20.0;
	K_COLL_ = 0.05;
	PHI_AH_ = 15.0;
	PHI_OT_ = 68.5;
	PHI_HO_ = 22.5;
	PHI_CR_ = 15.0;
	KAPPA_ = 20;
	K_P_ = 2;
	K_CHI_ = 1.3;
	K_DP_ = 1;
	K_DCHI_SB_ = 0.35;
	K_DCHI_P_ = 1;

	P_ca_last_ = 1;
	Chi_ca_last_ = 0;

	cost_ = INFINITY;

	Chi_ca_ << -90.0,-75.0,-60.0,-45.0,-30.0,-15.0,0.0,15.0,30.0,45.0,60.0,75.0,90.0;
	Chi_ca_ *= DEG2RAD;
	P_ca_ << -1.0, 0.0, 0.5, 1.0;

	asv = new shipModel(T_,DT_);

}


simulationBasedMpc::~simulationBasedMpc(){
}


void simulationBasedMpc::getBestControlOffset(double &u_os_best, double &psi_os_best, double u_d, double psi_d, const Eigen::Matrix<double,6,1>& asv_state, const Eigen::Matrix<double,-1,9>& obst_states){
	double cost = INFINITY;
	double cost_i = 0;
	double cost_k;
	int n_obst;

	if (obst_states.rows() == 0){
		u_os_best = 1;
		psi_os_best = 0;
		P_ca_last_ = 1;
		Chi_ca_last_ = 0;
		return;
	}else{
		for (int i = 0; i < obst_states.rows(); i++){
			obstacle *obst = new obstacle(obst_states.row(i), T_, DT_);
			obst_vect.push_back(obst);
		}
		n_obst = obst_vect.size();
	}

	for (int i = 0; i < Chi_ca_.rows(); i++){
		for (int j = 0; j < P_ca_.rows(); j++){

			asv->eulersMethod(asv_state, u_d*P_ca_[j], psi_d+ Chi_ca_[i]);

			cost_i = -1;
			for (int k = 0; k < n_obst; k++){

				cost_k = costFunction(P_ca_[j], Chi_ca_[i], k);
				if (cost_k > cost_i){
					cost_i = cost_k;	// Maximizing cost associated with this scenario
				}
			}

			if (cost_i < cost){
				cost = cost_i; 			// Minimizing the overall cost
				u_os_best = P_ca_[j];
				psi_os_best = Chi_ca_[i];
			}
		}
	}

	for (int k = 0; k < n_obst; k++){
		delete(obst_vect[k]);
	}
	obst_vect.clear();

	P_ca_last_ = u_os_best;
	Chi_ca_last_ = psi_os_best;
}


double simulationBasedMpc::costFunction(double P_ca, double Chi_ca, int k){
	double dist, phi, psi_o, psi_rel, R, C, k_coll, d_safe_i;
	Eigen::Vector2d d, los,v_o, v_s;
	bool mu, OT, SB, HO, CR;
	// TODO: Adjust radius front/back according to AIS data
	double combined_radius = 10; //asv->radius + obstacles_vect[k]->radius_;
	double d_safe = D_SAFE_ + combined_radius;
	double d_close = D_CLOSE_ + combined_radius;
	double H0 = 0;
	double H1 = 0;
	double H2 = 0;
	double cost = 0;
	double t = 0;
	double t0 = 0;
	int n_samp = T_/DT_;

	for (int i = 0; i < n_samp-1; i++){

		t += DT_;

		d(0) = obst_vect[k]->x_(i) - asv->x(i);
		d(1) = obst_vect[k]->y_[i] - asv->y[i];
		dist = d.norm();

		R = 0;
		C = 0;
		mu = 0;

		if (dist < d_close){

			v_o(0) = obst_vect[k]->u_[i];
			v_o(1) = obst_vect[k]->v_[i];
			rot2d(obst_vect[k]->psi_,v_o);

			v_s(0) = asv->u[i];
			v_s(1) = asv->v[i];
			rot2d(asv->psi[i],v_s);

			psi_o = obst_vect[k]->psi_;
			while(psi_o <= -M_PI) phi += 2*M_PI;
			while (psi_o > M_PI) phi -= 2*M_PI;

			phi = atan2(d(1),d(0)) - asv->psi[i];
			while(phi <= -M_PI) phi += 2*M_PI;
			while (phi > M_PI) phi -= 2*M_PI;

			psi_rel = psi_o - asv->psi[i];
			while(psi_rel < -M_PI) psi_rel += 2*M_PI;
			while(psi_rel > M_PI) psi_rel -= 2*M_PI;

			if (phi < 0 && psi_rel > 0){
				d_safe_i = 0.001*d_safe;
			}else{
				d_safe_i = d_safe;
			}


			if (dist < d_safe_i){
				R = (1/pow(fabs(t-t0),P_))*pow(d_safe/dist,Q_);
				k_coll = K_COLL_*combined_radius;
				C = k_coll*pow((v_s-v_o).norm(),2);
			}

			los = d/dist;

			// Overtaken by obstacle
			OT = v_s.dot(v_o) > cos(PHI_OT_*DEG2RAD)*v_s.norm()*v_o.norm()
					&& v_s.norm() < v_o.norm();
			// Obstacle on starboard side
			SB = phi < 0;
			// Obstacle Head-on
			HO = v_o.norm() > 0.05
					&& v_s.dot(v_o) < -cos(PHI_HO_*DEG2RAD)*v_s.norm()*v_o.norm()
					&& v_s.dot(los) > cos(PHI_AH_*DEG2RAD)*v_s.norm();
			// Crossing situation
			CR = v_s.dot(v_o) < cos(PHI_CR_*DEG2RAD)*v_s.norm()*v_o.norm()
					&& ((SB && psi_rel > 0 ));

			mu = ( SB && HO ) || ( CR && !OT);

		}

		H0 = C*R + KAPPA_*mu;

		if (H0 > H1){
			H1 = H0;  // Maximizing the cost with regards to time
		}
	}

	H2 = K_P_*(1-P_ca) + K_CHI_*pow(Chi_ca,2) + deltaP(P_ca) + deltaChi(Chi_ca);
	cost =  H1 + H2;

	return cost;
}

double simulationBasedMpc::deltaP(double P_ca){
	return K_DP_*std::abs(P_ca_last_ - P_ca);
}


double simulationBasedMpc::deltaChi(double Chi_ca){
	double dChi = Chi_ca - Chi_ca_last_;
	if (dChi > 0){
		return K_DCHI_SB_*pow(dChi,2);
	}else if (dChi < 0){
		return K_DCHI_P_*pow(dChi,2);
	}else{
		return 0;
	}
}


void simulationBasedMpc::rot2d(double yaw, Eigen::Vector2d &res){
	Eigen::Matrix2d R;
	R << cos(yaw), -sin(yaw),
		 sin(yaw), cos(yaw);
	res = R*res;
}
