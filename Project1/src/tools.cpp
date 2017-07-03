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
	float c3 = (c1*c2);

  //compute the Jacobian matrix
	Hj << (x/c2), (y/c2), 0, 0,
		    -(y/c1), (x/c1), 0, 0,
		    y*(x_dot*y - y_dot*x)/c3, x*(x*y_dot - y*x_dot)/c3, x/c2, y/c2;

	return Hj;
}

VectorXd Tools::cart_to_polar(const VectorXd& cart)
{
  VectorXd polar(3);

  float x = cart(0);
  float y = cart(1);
  float x_dot = cart(2);
  float y_dot = cart(3);

  float rho = sqrt(x*x + y*y);

  /*if (fabs(rho) < .0001) {
      rho = .0001;
    }*/

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

float Tools::normalize_theta(float theta)
{
  float normalized = theta;

  /*if (theta > 0)
  {
    while (fabs(normalized) > 3.14)
    {
      normalized = normalized - 6.28;
    }
  }
  else if (theta < 0)
  {
    while (fabs(normalized) > 3.14)
    {
      normalized = normalized + 6.28;
    }
  }*/
  normalized = atan2(sin(theta), cos(theta));

  return normalized;
}
