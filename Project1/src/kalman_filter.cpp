#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter()
{
  x_ = VectorXd(4);

  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

  Q_ = MatrixXd(4, 4);

  P_ = MatrixXd(4, 4);
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict()
{
  std::cout << F_ << std::endl;
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::UpdateKF(const VectorXd &z)
{
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;

  Update(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{

  VectorXd y = z - tools.cart_to_polar(x_);

  if (fabs(y(1)) > 3.14)
  {
    y(1) = atan2(sin(y(1)), cos(y(1)));
  }

  Update(y);
}

void KalmanFilter::Update(const VectorXd &y)
{
  MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
