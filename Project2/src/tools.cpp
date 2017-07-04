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

VectorXd Tools::ConvertRADAR(const VectorXd& z)
{
  VectorXd output(2);

  double x, y;

  x = output(0) * cos(output(1));
  y = output(0) * sin(output(1));

  output << x,y;

  return output;
}

double Tools::ComputeNIS()
{
  double NIS = 0;

  return NIS;
}
