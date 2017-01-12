/*
 * obstacle.h
 *
 *  Created on: Dec 22, 2016
 *      Author: ingerbha
 */

#ifndef OBSTACLE_H_
#define OBSTACLE_H_

#include <vector>
#include <Eigen/Dense>

class obstacle
{
	public:

	/// Constructor
	obstacle(const Eigen::Matrix<double,9,1>& state, double T, double dt);

	/// Destructor
	~obstacle();

	Eigen::VectorXd getX();
	Eigen::VectorXd getY();
	Eigen::VectorXd getU();
	Eigen::VectorXd getV();
	double getPsi();
	double getA();
	double getB();
	double getC();
	double getD();
	double getL();
	double getW();

	Eigen::VectorXd x_;
	Eigen::VectorXd y_;
	Eigen::VectorXd u_;
	Eigen::VectorXd v_;

	double psi_;
	double A_, B_, C_, D_, l, w;
	double os_x, os_y;

	private:

	void calculatePosOffsets();

	void calculateTrajectory();

	const int n_samp_;
	double T_;
	double dt_;

	//Elements of the rotation matrix
	double r11_, r12_, r21_, r22_;

};

#endif /* OBSTACLE_H_ */
