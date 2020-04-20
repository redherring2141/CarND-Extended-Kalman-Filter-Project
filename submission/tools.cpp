#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	//Check the validity of the following inputs
	assert(estimations.size()!=0);//The estimation vector size should not be zero
	assert(estimations.size() == ground_truth.size());//The estimation vector size should be equal to the ground truth vector size

	VectorXd diff;
	//Accumulate squared residuals
	for(int i=0; i<estimations.size(); ++i)
	{
		diff = estimations[i] - ground_truth[i];
		diff = diff.array() * diff.array();
		rmse += diff;
	}

	//Calculate the mean
	rmse = rmse/estimations.size();
	//Calculate the squared root
	rmse = rmse.array().sqrt();

	//Return the result
	return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3,4);
	//Recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//Pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//Check division by zero
	if(fabs(c1) < 0.0001)
	{
		cout << "CalculateJacobian() - Error - Division by Zero" << endl;
		return Hj;
	}

	//Compute the Jacobian matrix
	Hj <<	 (px/c2), (py/c2), 0, 0,
			-(py/c1), (px/c1), 0, 0,
			py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
