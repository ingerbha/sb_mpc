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
static const double RAD2DEG = 180.0f/M_PI;

simulationBasedMpc::simulationBasedMpc(){
	T_ = 400.0;
	DT_ = 0.1;
	method = LinearPrediction;

	P_ = 1.0;
	Q_ = 4.0;
	D_CLOSE_ = 200.0;
	D_SAFE_ = 40.0;
	K_COLL_ = 0.1;
	PHI_AH_ = 15.0;
	PHI_OT_ = 68.5;
	PHI_HO_ = 22.5;
	PHI_CR_ = 15.0;
	KAPPA_ = 3.0;
	K_P_ = 1.4;
	K_CHI_ = 1.3;
	K_DP_ = 1.3;
	K_DCHI_SB_ = 0.9;
	K_DCHI_P_ = 1.2;

	P_ca_last_ = 1;
	Chi_ca_last_ = 0;

	cost_ = INFINITY;

	Chi_ca_.resize(13);
	Chi_ca_ << -90.0,-75.0,-60.0,-45.0,-30.0,-15.0,0.0,15.0,30.0,45.0,60.0,75.0,90.0;
	Chi_ca_ *= DEG2RAD;
	P_ca_.resize(4);
	P_ca_ << -1.0, 0.0, 0.5, 1.0;

	asv = new shipModel(T_,DT_);

}


simulationBasedMpc::~simulationBasedMpc(){
}

double simulationBasedMpc::getT(){
	return T_;
}

double simulationBasedMpc::getDt(){
	return DT_;
}

double simulationBasedMpc::getP(){
	return P_;
}

double simulationBasedMpc::getQ(){
	return Q_;
}

double simulationBasedMpc::getDClose(){
	return D_CLOSE_;
}

double simulationBasedMpc::getDSafe(){
	return D_SAFE_;
}

double simulationBasedMpc::getKColl(){
	return K_COLL_;
}
double simulationBasedMpc::getPhiAH(){
	return PHI_AH_*RAD2DEG;
}

double simulationBasedMpc::getPhiOT(){
	return PHI_OT_*RAD2DEG;
}

double simulationBasedMpc::getPhiHO(){
	return PHI_HO_*RAD2DEG;
}

double simulationBasedMpc::getPhiCR(){
	return PHI_CR_*RAD2DEG;
}

double simulationBasedMpc::getKappa(){
	return KAPPA_;
}

double simulationBasedMpc::getKP(){
	return K_P_;
}

double simulationBasedMpc::getKdP(){
	return K_DP_;
}

double simulationBasedMpc::getKChi(){
	return K_CHI_;
}

double simulationBasedMpc::getKdChiSB(){
	return K_DCHI_SB_;
}

double simulationBasedMpc::getKdChiP(){
	return K_DCHI_P_;
}

Eigen::VectorXd simulationBasedMpc::getChiCA(){
	return Chi_ca_*RAD2DEG;
}

Eigen::VectorXd simulationBasedMpc::getPCA(){
	return P_ca_;
}

std::string simulationBasedMpc::getMethod(){
	std::string returnValue;
	switch (method){
		case EulerFirstOrder 	: returnValue = "EulerFirstOrder"; break;
		case LinearPrediction 	: returnValue = "LinearPrediction"; break;
		default 				: returnValue = "Failed";
	}
	return returnValue;
}

void simulationBasedMpc::setMethod(int i){
	switch (i){
	case 0 : method = EulerFirstOrder; break;
	case 1 : method = LinearPrediction; break;
	}
}

// Todo: Add validity checks for the set functions
void simulationBasedMpc::setT(double T){
	T_ = T;
}

void simulationBasedMpc::setDt(double dt){
	DT_ = dt;
}

void simulationBasedMpc::setP(double p){
	P_ = p;
}

void simulationBasedMpc::setQ(double q){
	Q_ = q;
}

void simulationBasedMpc::setDClose(double d_close){
	D_CLOSE_ = d_close;
}

void simulationBasedMpc::setDSafe(double d_safe){
	D_SAFE_ = d_safe;
}

void simulationBasedMpc::setKColl(double k_coll){
	K_COLL_ = k_coll;
}

void simulationBasedMpc::setPhiAH(double phi_AH){
	PHI_AH_ = phi_AH*DEG2RAD;
}

void simulationBasedMpc::setPhiOT(double phi_OT){
	PHI_OT_ = phi_OT*DEG2RAD;
}

void simulationBasedMpc::setPhiHO(double phi_HO){
	PHI_HO_ = phi_HO*DEG2RAD;
}

void simulationBasedMpc::setPhiCR(double phi_CR){
	PHI_CR_ = phi_CR*DEG2RAD;
}

void simulationBasedMpc::setKappa(double kappa){
	KAPPA_ = kappa;
}

void simulationBasedMpc::setKP(double K_P){
	K_P_ = K_P;
}

void simulationBasedMpc::setKdP(double K_dP){
	K_DP_ = K_dP;
}

void simulationBasedMpc::setKChi(double K_Chi){
	K_CHI_ = K_Chi;
}

void simulationBasedMpc::setKdChiSB(double K_dChi_SB){
	K_DCHI_SB_ = K_dChi_SB;
}

void simulationBasedMpc::setKdChiP(double K_dChi_P){
	K_DCHI_P_ = K_dChi_P;
}

void simulationBasedMpc::setChiCA(Eigen::VectorXd Chi_ca){
	Chi_ca_.resize(Chi_ca.size());
	Chi_ca_ = Chi_ca*DEG2RAD;
}

void simulationBasedMpc::setPCA(Eigen::VectorXd P_ca){
	P_ca_.resize(P_ca.size());
	P_ca_ = P_ca;
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

	for (int i = 0; i < Chi_ca_.size(); i++){
		for (int j = 0; j < P_ca_.size(); j++){

			switch(method){
			case EulerFirstOrder : asv->eulersMethod(asv_state, u_d*P_ca_[j], psi_d + Chi_ca_[i]);
			case LinearPrediction : asv->linearPrediction(asv_state, u_d*P_ca_[j], psi_d + Chi_ca_[i]);
			}


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
	double dist, phi, phi_o, psi_o, psi_rel, R, C, k_coll, d_safe_i;
	Eigen::Vector2d d, los, los_inv, v_o, v_s;
	bool mu, OT, SB, HO, CR;
	double combined_radius = asv->getL() + obst_vect[k]->getL();
	double d_safe = D_SAFE_;
	double d_close = D_CLOSE_;
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
			while(psi_o <= -M_PI) psi_o += 2*M_PI;
			while (psi_o > M_PI) psi_o -= 2*M_PI;

			phi = atan2(d(1),d(0)) - asv->psi[i];
			while(phi <= -M_PI) phi += 2*M_PI;
			while (phi > M_PI) phi -= 2*M_PI;

			psi_rel = psi_o - asv->psi[i];
			while(psi_rel < -M_PI) psi_rel += 2*M_PI;
			while(psi_rel > M_PI) psi_rel -= 2*M_PI;

			los = d/dist;
			los_inv = -d/dist;

			// Calculating d_safe
			if (phi < PHI_AH_){//v_s.dot(los) > cos(PHI_AH_*DEG2RAD)*v_s.norm()){ // obst ahead
				d_safe_i = d_safe + asv->getL()/2;
			}else if (phi > PHI_OT_){//v_s.dot(los) > cos(PHI_OT_*DEG2RAD)*v_s.norm()){ // obst behind
				d_safe_i = 0.5*d_safe + asv->getL()/2;
			}else{
				d_safe_i = d_safe + asv->getW()/2;
			}

			phi_o = atan2(-d(1),-d(0)) - obst_vect[k]->psi_;
			while(phi_o <= -M_PI) phi_o += 2*M_PI;
			while (phi_o > M_PI) phi_o -= 2*M_PI;

			if (phi_o < PHI_AH_){//v_o.dot(los_inv) > cos(PHI_AH_*DEG2RAD)*v_o.norm()){ // ship ahead
				d_safe_i += d_safe + obst_vect[k]->getL()/2;
			}else if(phi_o > PHI_OT_){//v_o.dot(los_inv) > cos(PHI_OT_*DEG2RAD)*v_o.norm()){ // ship behind
				d_safe_i += 0.5*d_safe + obst_vect[k]->getL()/2;
			}else{
				d_safe_i += d_safe + obst_vect[k]->getW()/2;
			}

			if (v_s.dot(v_o) > cos(PHI_OT_*DEG2RAD)*v_s.norm()*v_o.norm() && v_s.norm() > v_o.norm()){
				d_safe_i = d_safe + asv->getL()/2 + obst_vect[k]->getL()/2;
			}


			if (dist < d_safe_i){
				R = (1/pow(fabs(t-t0),P_))*pow(d_safe/dist,Q_);
				k_coll = K_COLL_*asv->getL()*obst_vect[k]->getL();
				C = k_coll*pow((v_s-v_o).norm(),2);
			}



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
