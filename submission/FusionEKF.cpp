#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  H_laser_ <<	1, 0, 0, 0,
				0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  ekf_.F_ = MatrixXd(4,4);
  ekf_.H_ = H_laser_;

  aaT_ = MatrixXd(2,2);
  aaT_ << 9, 0,
		 0, 9;

  G_ = MatrixXd(4,2);
  G_ << 0, 0,
		0, 0,
		0, 0,
		0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}


//Reference: Lesson24-14, Lesson24-10
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;


	//State covariance matrix P_
	ekf_.P_ = MatrixXd(4,4);
	ekf_.P_ <<  1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1000, 0,
				0, 0, 0, 1000;

	//The initial transition matrix F_
	ekf_.F_ = MatrixXd(4,4);
	ekf_.F_ <<  1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;

	
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	{
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
		float rho = measurement_pack.raw_measurements_[0];
		float phi = measurement_pack.raw_measurements_[1];
		float x = rho*cos(phi);
		float y = rho*sin(phi);
		ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
	{
      // TODO: Initialize state.
		float x = measurement_pack.raw_measurements_[0];
		float y = measurement_pack.raw_measurements_[1];
		ekf_.x_ << x, y, 0, 0;
    }

	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  float dt2 = dt*dt;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
			 0, 1, 0, dt,
			 0, 0, 1,  0,
			 0, 0, 0,  1;

  G_ <<	dt2 / 2, 0,
		0, dt2 /2,
		dt, 0,
		0, dt;

  ekf_.Q_ = G_ * aaT_ * G_.transpose();
  ekf_.Predict();


  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
	  VectorXd z = VectorXd(3);
	  float rho = measurement_pack.raw_measurements_[0];
	  float phi = measurement_pack.raw_measurements_[1];
	  float rho_dot = measurement_pack.raw_measurements_[2];
	  z << rho, phi, rho_dot;

	  ekf_.R_ = R_radar_;
	  z += ekf_.R_.diagonal();

	  ekf_.UpdateEKF(z);
  }
  else
  {
    // TODO: Laser updates
	  VectorXd z = VectorXd(2);
	  float x = measurement_pack.raw_measurements_[0];
	  float y = measurement_pack.raw_measurements_[1];
	  z << x,y;

	  ekf_.R_ = R_laser_;
	  z += ekf_.R_.diagonal();
	  ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
