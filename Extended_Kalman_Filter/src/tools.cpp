#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth)
{
  VectorXd rmse(4);
	rmse << 0,0,0,0;

  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

  //accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i)
  {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  MatrixXd Hj(3,4);
	//recover state parameters
	float x = x_state(0);
	float y = x_state(1);
	float x_dot = x_state(2);
	float y_dot = x_state(3);

  float x2 = x * x;
  float y2 = y * y;

  //pre-compute a set of terms to avoid repeated calculation
	float c1 = x2 + y2;

  float epsilon = .0001;
  if (fabs(c1) < epsilon)
  {
    c1 = epsilon;
  }

	float c2 = sqrt(c1);
  if (fabs(c2) < epsilon)
  {
    c2 = epsilon;
  }

	float c3 = (c1*c2);
  if (fabs(c3) < epsilon)
  {
    c3 = epsilon;
  }

  //compute the Jacobian matrix
  Hj(0,0) = x/c2;
  Hj(0,1) = y/c2;
  Hj(0,2) = 0;
  Hj(0,3) = 0;

  Hj(1,0) = -y/c1;
  Hj(1,1) = x/c1;
  Hj(1,2) = 0;
  Hj(1,3) = 0;

  Hj(2,0) = y*(x_dot*y - y_dot*x)/c3;
  Hj(2,1) = x*(x*y_dot - y*x_dot)/c3;
  Hj(2,2) = x/c2;
  Hj(2,3) = y/c2;

	return Hj;
}

VectorXd Tools::cart_to_polar(const VectorXd& cart)
{
  //printf("here\n");
  VectorXd polar(3);

  float x = cart(0);
  float y = cart(1);
  float x_dot = cart(2);
  float y_dot = cart(3);

  float epsilon = .0001;
  if (fabs(x) < epsilon) {
      x = epsilon;
  }
  if (fabs(x_dot) < epsilon) {
      x_dot = epsilon;
  }

  float rho = sqrt(x*x + y*y);

  if (fabs(rho) < epsilon)
  {
      rho = epsilon;
  }

  float theta = atan2(y, x);
  float rho_dot = (x*x_dot + y*y_dot) / rho;

  polar << rho, theta, rho_dot;

  return polar;
}

VectorXd Tools::polar_to_cart(const VectorXd& polar)
{
  VectorXd cart(4);

  float rho = polar(0);
  float theta = polar(1);
  float rho_dot = polar(2);

  float x = rho * cos(theta);
  float y = rho * sin(theta);
  float x_dot = rho_dot * cos(theta);
  float y_dot = rho_dot * sin(theta);

  cart << x, y, x_dot, y_dot;

  return cart;
}
