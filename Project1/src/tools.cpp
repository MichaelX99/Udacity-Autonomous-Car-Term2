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
	double x = x_state(0);
	double y = x_state(1);
	double x_dot = x_state(2);
	double y_dot = x_state(3);

  double x2 = x * x;
  double y2 = y * y;

  //pre-compute a set of terms to avoid repeated calculation
	double c1 = x2 + y2;

  double epsilon = .0001;
  if (fabs(c1) < epsilon)
  {
    c1 = epsilon;
  }

	double c2 = sqrt(c1);
  if (fabs(c2) < epsilon)
  {
    c2 = epsilon;
  }

	double c3 = (c1*c2);
  if (fabs(c3) < epsilon)
  {
    c3 = epsilon;
  }

  //compute the Jacobian matrix
	Hj << (x/c2), (y/c2), 0, 0,
		    -(y/c1), (x/c1), 0, 0,
		    y*(x_dot*y - y_dot*x)/c3, x*(x*y_dot - y*x_dot)/c3, x/c2, y/c2;
  /*Hj(0,0) = x/c2;
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
  Hj(2,3) = y/c2;*/

	return Hj;
}

VectorXd Tools::cart_to_polar(const VectorXd& cart)
{
  VectorXd polar(3);

  double x = cart(0);
  double y = cart(1);
  double x_dot = cart(2);
  double y_dot = cart(3);

  double epsilon = .0001;
  if (fabs(x) < epsilon) {
      x = epsilon;
  }
  if (fabs(x_dot) < epsilon) {
      x_dot = epsilon;
  }

  double rho = sqrt(x*x + y*y);

  if (fabs(rho) < epsilon)
  {
      rho = epsilon;
  }

  double theta = atan2(y, x);
  double rho_dot = (x*x_dot + y*y_dot) / rho;

  polar << rho, theta, rho_dot;

  return polar;
}

VectorXd Tools::polar_to_cart(const VectorXd& polar)
{
  VectorXd cart(4);

  double rho = polar(0);
  double theta = polar(1);
  double rho_dot = polar(2);

  double x = rho * cos(theta);
  double y = rho * sin(theta);
  double x_dot = rho_dot * cos(theta);
  double y_dot = rho_dot * sin(theta);

  cart << x, y, x_dot, y_dot;

  return cart;
}
