#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO: (completed)
    * (done) predict the state
      Following coding example in lesson 25.14
  */

    x_ = F_*x_; // New state prediction based on motion model
    //MatrixXd Ft = F_.transpose();
    P_ = F_*P_*F_.transpose() + Q_; // New process covariance prediction  based on motion model
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO: (completed)
    * (done) update the state by using Kalman Filter equations
      Following coding example in lesson 25.14 and 25.15
  */

    VectorXd z_pred = H_ * x_;        // Project state prediction to measurement space
    VectorXd y = z - z_pred;          // Calculate residual to measurement
    MatrixXd Ht = H_.transpose();     // Prepare projection from measurement space to state space
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    x_ = x_ + (K * y);                // Update state with projected residual
    long x_size = x_.size();          // Update process coavriance matrix
    P_ = (MatrixXd::Identity(x_size, x_size) - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO: (dompleted)
    * (done) update the state by using Extended Kalman Filter equations
  */

    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];

    float rho = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    float rho_dot = (px*vx + py*vy)/rho;
    // Ignore radial velocity is radius is too small
    if (fabs(rho) < 0.0001) {
      rho_dot = 0;
    }

    // Calculate the exatc projection from state to measurement space: z_p =  h(x_p)
    VectorXd z_p(3);
    z_p << rho, phi, rho_dot;

    // Calculate residual y
    VectorXd y = z - z_p;

    // Adjust phi to be in [-pi,pi]
    if (y[1] > M_PI){
        y[1] -= 2*M_PI;
    } else if (y[1] < -M_PI){
        y[1] += 2*M_PI;
    }

    // Prepare projection matrixes from measurement space to state space
    // Identical to the standard Kalman filter
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_*P_*Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_*Ht;
    MatrixXd K = PHt*Si;

    // Project measurment space to state space
    x_ = x_ + (K*y);
    long x_size = x_.size();
    P_ = (MatrixXd::Identity(x_size, x_size) - K*H_)*P_;

}
