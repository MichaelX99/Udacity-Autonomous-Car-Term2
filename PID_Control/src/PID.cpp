#include "PID.h"

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double in_Kp, double in_Ki, double in_Kd) {
  Kp = in_Kp;
  Ki = in_Ki;
  Kd = in_Kd;

  p_error = 0;
  i_error = 0;
  d_error = 0;
}

void PID::UpdateError(double cte) {
  d_error = cte - p_error;
  i_error += cte;
  p_error = cte;
}

double PID::TotalError() {
  double output;
  output = -Kp * p_error - Ki * i_error - Kd * d_error;

  return output;
}
