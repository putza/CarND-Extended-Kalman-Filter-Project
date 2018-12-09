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
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

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

  /**
  TODO: (Completed)
    * (done) Finish initializing the FusionEKF.
    * (done) Set the process and measurement noises
  */

  H_laser_ << 1, 0 ,0 ,0,
              0, 1 , 0, 0;

  // Initial Transition Matrix
  ekf_.F_ = MatrixXd(4, 4); // Initializing the state prediction model
  ekf_.F_ << 1, 0, 1 , 0,
             0, 1, 0, 1,
             0, 0, 1, 0,    // 1 needs to be replaced by delta d
             0, 0, 0, 1;    // 1 needs to be replaced by delta d

  // State Covariance Matrix
  ekf_.P_ = MatrixXd(4, 4); // Initializing process covariance matrix
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0, // High penalty for initial velocity prediction vx
             0, 0, 0, 1000; // High penalty for initial velocity prediction vy

  //Tools tools; Initialized in header
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO: (completed)
      * (done) Initialize the state ekf_.x_ with the first measurement.
      * (done) Create the covariance matrix.
      * (done) Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0.001, 0.001, 0.001, 0.001;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      Note:
        * rho = measurement_pack.raw_measurements_[0];
        * phi = measurement_pack.raw_measurements_[1];
        * rho_dot = measurement_pack.raw_measurements_[2];
      */

        // Initial guess for position
        float rho = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        ekf_.x_(0) = rho*cos(phi);
        ekf_.x_(1) = rho*sin(phi);

        // Try a bad guess for the velocity. Seems to give slightly better results than initializing with 0
        //loat rho_dot = measurement_pack.raw_measurements_[2];
        //ekf_.x_(2) = rho_dot*cos(phi);
        //ekf_.x_(3) = rho_dot*sin(phi);

        // Alternative: Initialize with zero
        ekf_.x_(2) = 0;
        ekf_.x_(3) = 0;


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state. Taken from Lesson 25, coding example
      */
        ekf_.x_(0) = measurement_pack.raw_measurements_[0];
        ekf_.x_(1) = measurement_pack.raw_measurements_[1];

        ekf_.x_(2) = 0;
        ekf_.x_(3) = 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_; // Initialize timestamp with first measurement time.

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO: (completed)
     * (done) Update the state transition matrix F according to the new elapsed time.
        - (done) Time is measured in seconds.
     * (done) Update the process noise covariance matrix.
     * (done)Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     * Taken from Lesson 25, coding examples
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt/4;
  dt_3 = dt_3/2;

  // Replace placeholder values with the correct dt values in the motion model
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update Process Noise Covariance Matrix Q
  ekf_.Q_ = MatrixXd(4, 4);

  float noise_ax = 9.0;
  float noise_ay = 9.0;


  ekf_.Q_ <<  dt_4*noise_ax, 0, dt_3*noise_ax, 0,
             0, dt_4*noise_ay, 0, dt_3*noise_ay,
             dt_3*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO: (completed)
     * (done) Use the sensor type to perform the update step.
     * (done) Update the state and covariance matrices.
   */


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);;          // Assign Jacobian to ekf class
    ekf_.R_ = R_radar_;                                   // Assign R
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);   // Execute EKF update step


  } else {
    // Laser updates
    ekf_.H_ = H_laser_;                                   // Assign state->measurement projection matrix
    ekf_.R_ = R_laser_;                                   // Assign R
    ekf_.Update(measurement_pack.raw_measurements_);      // Execute classical Kalman update step

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
