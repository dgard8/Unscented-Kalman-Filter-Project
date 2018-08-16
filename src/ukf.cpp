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
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  // initial state vector
  x_ = VectorXd(5);
  
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;
  
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/10.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_aug_ = 7;
  n_x_ = 5;
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    x_ << 1, 1, .1, .1, .1;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      
      float x = rho * cos(theta);
      float y = rho * sin(theta);
      
      x_(0) = x;
      x_(1) = y;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }
    
    // TODO: probably change these values depending on the type of measurement
    P_ << .01, 0, 0, 0, 0,
    0, .01, 0, 0, 0,
    0, 0, .1, 0, 0,
    0, 0, 0, .1, 0,
    0, 0, 0, 0, .1;
    
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd sigmaPoints = Tools::calculateSigmaPoint(x_, P_, std_a_, std_yawdd_);
  Xsig_pred_ = predictSigmaPoints(sigmaPoints, delta_t);
  updateState(Xsig_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2; // lidar measures px and py
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  for (int i=0; i<2 * n_aug_ + 1; i++){
    VectorXd originalPoint = Xsig_pred_.col(i);
    float px = originalPoint(0);
    float py = originalPoint(1);
    
    Zsig.col(i) << px, py;
  }
  
  MatrixXd measurmentNoise = MatrixXd(n_z,n_z);
  measurmentNoise << std_laspx_*std_laspx_, 0,
  0, std_laspy_*std_laspy_;
  
  updateMeasurement(Zsig, meas_package, measurmentNoise);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3; // radar measures r, phi, and phi_dot
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //transform sigma points into measurement space
  for (int i=0; i<2 * n_aug_ + 1; i++){
    VectorXd originalPoint = Xsig_pred_.col(i);
    float px = originalPoint(0);
    float py = originalPoint(1);
    float v = originalPoint(2);
    float psi = originalPoint(3);
    
    float p = sqrt(px*px + py*py);
    float theta = atan2(py,px);
    float pdot = (px*cos(psi)*v + py*sin(psi)*v)/p;
    
    Zsig.col(i) << p, theta, pdot;
  }
  
  MatrixXd measurmentNoise = MatrixXd(n_z,n_z);
  measurmentNoise << std_radr_*std_radr_, 0, 0,
  0, std_radphi_*std_radphi_, 0,
  0, 0,std_radrd_*std_radrd_;
  
  updateMeasurement(Zsig, meas_package, measurmentNoise);
  
}

void UKF::updateMeasurement(MatrixXd Zsig, MeasurementPackage meas_package, MatrixXd measurmentNoise){
  int n_z = Zsig.col(0).size();
  
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S and cross correlation matrix
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    if (n_z == 3){
      z_diff(1) = Tools::normalizeAngle(z_diff(1));
    }
    x_diff(3) = Tools::normalizeAngle(x_diff(3));
    
    S = S + weights_(i) * z_diff * z_diff.transpose();
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  S = S + measurmentNoise;
  
  
  //get measurement
  VectorXd z = meas_package.raw_measurements_;
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  if (n_z == 3){
    z_diff(1)=Tools::normalizeAngle(z_diff(1));
  }
  float nis = z_diff.transpose() * S.inverse() * z_diff;
  nisValues.push_back(nis);
  
  x_ = x_ + K * (z_diff);
  P_ = P_ - K * S * K.transpose();

}

MatrixXd UKF::predictSigmaPoints(MatrixXd points, double delta_t){
  int n_aug = n_aug_;
  int n_x = n_x_;
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  
  for (int i=0; i<2 * n_aug + 1; i++){
    VectorXd original_point = points.col(i);
    float px = original_point(0);
    float py = original_point(1);
    float v = original_point(2);
    float psi = original_point(3);
    float psidot = original_point(4);
    float tangential_acceleration = original_point(5);
    float angle_acceleration = original_point(6);
    
    original_point = VectorXd(n_x);
    original_point << px, py, v, psi, psidot;
    
    VectorXd movement = VectorXd(n_x);
    if (fabs(psidot) < 0.001){
      movement(0) = v*cos(psi)*delta_t;
      movement(1) = v*sin(psi)*delta_t;
    }
    else {
      movement(0) = (v/psidot)*(sin(psi+psidot*delta_t)-sin(psi));
      movement(1) = (v/psidot)*(-cos(psi+psidot*delta_t)+cos(psi));
    }
    movement(2) = 0;
    movement(3) = psidot*delta_t;
    movement(4) = 0;
    
    VectorXd acceleration_noise = VectorXd(n_x);
    acceleration_noise(0) = 0.5 * delta_t * delta_t * cos(psi) * tangential_acceleration;
    acceleration_noise(1) = 0.5 * delta_t * delta_t * sin(psi) * tangential_acceleration;
    acceleration_noise(2) = delta_t * tangential_acceleration;
    acceleration_noise(3) = 0.5 * delta_t * delta_t * angle_acceleration;
    acceleration_noise(4) = delta_t * angle_acceleration;
    
    VectorXd new_point = original_point + movement + acceleration_noise;
    Xsig_pred.col(i) = new_point;
  }
  
  return Xsig_pred;
}

void UKF::updateState(MatrixXd predictedPoints){
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd x_diff = predictedPoints.col(i) - x_;
    
    x_diff(3) = Tools::normalizeAngle(x_diff(3));
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}
