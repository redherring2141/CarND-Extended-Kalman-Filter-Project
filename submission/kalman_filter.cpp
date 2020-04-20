#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
//Reference: Lesson24-8
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  /**
   * TODO: predict the state
   */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
	VectorXd z_pred = H_*x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_*Ht;
	MatrixXd S = H_*PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt*Si;

	//New estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
	//Recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	VectorXd h = VectorXd(3);
//	assert(px>0.0001);
//	assert(py>0.0001);
	if(px<0.0001 || py<0.0001)
		return;
	float rho = sqrt(px*px + py*py);
	float phi = atan2(py, px);
	float rho_dot = (px*vx + py*vy)/rho;
//	assert(phi!=0);
//	assert(rho_dot!=0);
	if(phi==0 || rho_dot==0)
		return;
	h << rho, phi, rho_dot;

	VectorXd y = z-h;
	phi = y(1);

	if(phi > M_PI)
	{
		phi = phi - M_PI*2.0;
	}
	else if (phi < -M_PI)
	{
		phi = phi + M_PI_2 * 2.0;
	}
	y(1) = phi;

	MatrixXd Hj = Tools_.CalculateJacobian(x_);
	MatrixXd Ht = Hj.transpose();
	MatrixXd PHt = P_*Ht;
	MatrixXd S = Hj*PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt*Si;

	//New estimate
	x_ = x_ + K*y;

	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K*Hj) * P_;
}
