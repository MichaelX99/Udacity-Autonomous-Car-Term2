#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = .5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = .5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;


  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_x_;
  time_us_ = 0;
  Xsig_pred_ = MatrixXd(n_x_, n_aug_ * 2 + 1);

  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.5/(n_aug_+lambda_));
  weights_(0) = lambda_/(lambda_+n_aug_);

  n_radar_ = 3;
  n_lidar_ = 2;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  if (not is_initialized_)
  {
    x_(3) = 0;
    x_(4) = 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      x_.head(2) = tools.ConvertRADAR(meas_package.raw_measurements_);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    x_(2) = .1;

    P_ = MatrixXd::Identity(5, 5);

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  dt_ = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction();

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

  cout << "x = " << x_ << endl;
  cout << "P = " << P_ << endl;
  cout << "NIS = " << NIS_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt_ the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction()
{
  // Compute the Sigma Points of the current state distribution
  MatrixXd Xsig_aug = ComputeSigmaPoints();

  // Predict the Sigma points at the next time step with our nonlinear model
  PredictSigmaPoints(Xsig_aug);

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools.AngleNormalization(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

MatrixXd UKF::ComputeSigmaPoints()
{
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.col(0)  = x_aug;
  double coef = sqrt(lambda_+n_aug_);
  for (int i = 0; i < n_aug_; i++)
  {
    auto second( coef * A.col(i) );
    Xsig_aug.col(i+1) = x_aug + second;
    Xsig_aug.col(i+1+n_aug_) = x_aug - second;
  }

  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_aug)
{
  //predict sigma points forward in time
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double x = Xsig_aug(0,i);
    double y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double p_x, p_y;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
        p_x = x + v/yawd * ( sin (yaw + yawd*dt_) - sin(yaw));
        p_y = y + v/yawd * ( cos(yaw) - cos(yaw+yawd*dt_) );
    }
    else
    {
        p_x = x + v*dt_*cos(yaw);
        p_y = y + v*dt_*sin(yaw);
    }

    double p_v = v;
    double yaw_p = yaw + yawd*dt_;
    double yawd_p = yawd;

    //add noise
    p_x = p_x + 0.5 * nu_a * dt_ * dt_ * cos(yaw);
    p_y = p_y + 0.5 * nu_a * dt_ * dt_ * sin(yaw);
    p_v = p_v + nu_a*dt_;

    yaw_p = yaw_p + 0.5 * nu_yawdd * dt_ * dt_;
    yawd_p = yawd_p + nu_yawdd * dt_;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = p_x;
    Xsig_pred_(1,i) = p_y;
    Xsig_pred_(2,i) = p_v;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  // Move Sigma Points to Measurement space
  MatrixXd Zsig = SigmaToLIDAR();

  //mean predicted measurement
  VectorXd z_pred = MeasurementPredict(meas_package.sensor_type_, Zsig);

  //measurement covariance matrix S
  MatrixXd S = ComputeMeasurementCovariance(meas_package.sensor_type_, Zsig, z_pred);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_lidar_,n_lidar_);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = ComputeCrossCorrelation(meas_package.sensor_type_, Zsig, z_pred);

  MeasurementStateUpdate(meas_package.raw_measurements_, Tc, S, z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  // Move Sigma Points to Measurement space
  MatrixXd Zsig = SigmaToRADAR();

  //mean predicted measurement
  VectorXd z_pred = MeasurementPredict(meas_package.sensor_type_, Zsig);

  //measurement covariance matrix S
  MatrixXd S = ComputeMeasurementCovariance(meas_package.sensor_type_, Zsig, z_pred);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_radar_,n_radar_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = ComputeCrossCorrelation(meas_package.sensor_type_, Zsig, z_pred);

  MeasurementStateUpdate(meas_package.raw_measurements_, Tc, S, z_pred);
}

MatrixXd UKF::SigmaToRADAR()
{
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_radar_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points
    // extract values for better readibility
    double x = Xsig_pred_(0,i);
    double y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    if (fabs(x) < .001)
    {
      x = .001;
    }

    // measurement model
    Zsig(0,i) = sqrt(x*x + y*y);                        //r

    if (fabs(Zsig(0,i)) < .001)
    {
      Zsig(0,i) = .001;
    }

    Zsig(1,i) = atan2(y,x);                                 //phi
    Zsig(2,i) = (x*v1 + y*v2 ) / Zsig(0,i);   //r_dot
  }

  return Zsig;
}

MatrixXd UKF::SigmaToLIDAR()
{
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_lidar_, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }

  return Zsig;
}

VectorXd UKF::MeasurementPredict(MeasurementPackage::SensorType type, MatrixXd Zsig)
{
  //mean predicted measurement
  VectorXd z_pred;
  if (type == MeasurementPackage::LASER)
  {
    z_pred = VectorXd(n_lidar_);
  }
  else if (type == MeasurementPackage::RADAR)
  {
    z_pred = VectorXd(n_radar_);
  }
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  return z_pred;
}

MatrixXd UKF::ComputeMeasurementCovariance(MeasurementPackage::SensorType type, MatrixXd Zsig, VectorXd z_pred)
{
  MatrixXd S;
  if (type == MeasurementPackage::LASER)
  {
    S = MatrixXd(n_lidar_, n_lidar_);
  }
  else if (type == MeasurementPackage::RADAR)
  {
    S = MatrixXd(n_radar_, n_radar_);
  }
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools.AngleNormalization(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  return S;
}

MatrixXd UKF::ComputeCrossCorrelation(MeasurementPackage::SensorType type, MatrixXd Zsig, VectorXd z_pred)
{
  MatrixXd Tc;
  if (type == MeasurementPackage::LASER)
  {
    Tc = MatrixXd(n_x_, n_lidar_);
  }
  else if (type == MeasurementPackage::RADAR)
  {
    Tc = MatrixXd(n_x_, n_radar_);
  }
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools.AngleNormalization(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools.AngleNormalization(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  return Tc;
}

void UKF::MeasurementStateUpdate(VectorXd measurement, MatrixXd Tc, MatrixXd S, VectorXd z_pred)
{
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = measurement - z_pred;

  //angle normalization
  z_diff(1) = tools.AngleNormalization(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_ = z_diff.transpose() * S.inverse() * z_diff;
}
