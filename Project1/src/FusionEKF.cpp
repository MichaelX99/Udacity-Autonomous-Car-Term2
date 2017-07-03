#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  float noise_ax_ = 9.0;
  float noise_ay_ = 9.0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    // first measurement
    ekf_.x_ = VectorXd(4);

    float x, y;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd polar(3);
      float rho = measurement_pack.raw_measurements_[0]; // range
	    float theta = measurement_pack.raw_measurements_[1]; // bearing
	    float rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
      polar << rho, theta, rho_dot;


      VectorXd cart(4);
      cart = tools.polar_to_cart(polar);
      x = cart(0);
      y = cart(1);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      x = measurement_pack.raw_measurements_[0];
	    y = measurement_pack.raw_measurements_[1];
    }
    ekf_.x_ << x, y, 0.0, 0.0;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.P_ << 1, 0, 0, 0,
			         0, 1, 0, 0,
			         0, 0, 1000, 0,
			         0, 0, 0, 1000;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   //dt_ = 10000.0 * (measurement_pack.timestamp_ - previous_timestamp_);
   dt_ = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	 previous_timestamp_ = measurement_pack.timestamp_;

	 dt_2_ = dt_ * dt_;
	 dt_3_ = dt_2_ * dt_;
	 dt_4_ = dt_3_ * dt_;

   //Modify the F matrix so that the time is integrated
   ekf_.F_ << 1, 0, dt_, 0,
			        0, 1, 0, dt_,
			        0, 0, 1, 0,
			        0, 0, 0, 1;

	//set the process covariance matrix Q
	ekf_.Q_ <<  dt_4_/4*noise_ax_, 0, dt_3_/2*noise_ax_, 0,
			        0, dt_4_/4*noise_ay_, 0, dt_3_/2*noise_ay_,
			        dt_3_/2*noise_ax_, 0, dt_2_*noise_ax_, 0,
			        0, dt_3_/2*noise_ay_, 0, dt_2_*noise_ay_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  //ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
  {
    // Laser updates
    ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  //ekf_.UpdateKF(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
