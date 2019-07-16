#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <iostream>

#include "bayesian_filter.h"
#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

#include "Eigen/Dense"

using std::cout;
using std::endl;

class Car : public Process {
 public:
  Car() {
    x_ = Vector(4);
    dt_ = 0.05;
    P_ = Matrix(4, 4);
    P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;
  }
  virtual ~Car() = default;
  void init(Vector z) { x_ = z; };
  void dt(double t) {  dt_ = t;  }

  virtual const Matrix F() {
    F_ = Matrix(4, 4);
    F_ << 1, 0, dt_, 0, 0, 1, 0, dt_, 0, 0, 1, 0, 0, 0, 0, 1;
    return F_;
  }
  virtual const Matrix P() { return P_; }
  virtual const Matrix Q() {
    double dt = dt_;
    double noise_ax = 9.0;
    double noise_ay = 9.0;

    double dt_2 = dt * dt;     // dt^2
    double dt_3 = dt_2 * dt;   // dt^3
    double dt_4 = dt_3 * dt;   // dt^4
    double dt_4_4 = dt_4 / 4;  // dt^4/4
    double dt_3_2 = dt_3 / 2;  // dt^3/2
    Q_ = Matrix(4, 4);
    Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0, 0, dt_4_4 * noise_ay, 0,
        dt_3_2 * noise_ay, dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0, 0,
        dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
    return Q_; }
};

class Laser : public Measurement {
 public:
  Laser() {
    R_ = Matrix(2, 2);
    R_ << 0.0225, 0, 0, 0.0225;
    H_ = Matrix(2, 4);
    H_ << 1, 0, 0, 0, 0, 1, 0, 0;
  }

  virtual ~Laser() = default;

  void measure(const Vector &z, const Vector &x_) { y_ = z - H_ * x_; }

  virtual Matrix H() { return H_; }
  virtual Matrix R() { return R_; }
  virtual Vector y() { return y_; }
};

class Radar : public Measurement {
 public:
  Radar() {
    R_ = Matrix(3, 3);
    R_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;
    H_ = Matrix(3, 4);
  }

  virtual ~Radar() = default;

  void measure(const Vector &z, const Vector &x_) {
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    double rho = sqrt(px * px + py * py);
    double theta = atan2(py, px);
    double rho_dot = (px * vx + py * vy) / rho;
    Vector h = Vector(3);
    h << rho, theta, rho_dot;
    Vector y = z - h;
    while (y(1) > M_PI || y(1) < -M_PI) {
      if (y(1) > M_PI) {
        y(1) -= M_PI;
      } else {
        y(1) += M_PI;
      }
    }

    y_ = y;

    H_ = Jacobian(x_);


  }

  Matrix Jacobian(const VectorXd& x_state){
    MatrixXd Hj(3, 4);

    if (x_state.size() != 4) {
      cout
          << "ERROR - CalculateJacobian () - The state vector must have size 4."
          << endl;
      return Hj;
    }
    // recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    // pre-compute a set of terms to avoid repeated calculation
    double c1 = px * px + py * py;
    double c2 = sqrt(c1);
    double c3 = (c1 * c2);

    // check division by zero
    if (fabs(c1) < 0.0001) {
      cout << "ERROR - CalculateJacobian () - Division by Zero" << endl;
      return Hj;
    }

    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;

  }

  virtual Matrix H() { return H_; }
  virtual Matrix R() { return R_; }
  virtual Vector y() { return y_; }
};

class FusionEKF {
 public:
  /**
   * Constructor.
   */
  FusionEKF();

  /**
   * Destructor.
   */
  virtual ~FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
   * Kalman Filter update and prediction math lives in here.
   */

  std::shared_ptr<Car> car_;
  std::shared_ptr<Laser> laser_;
  std::shared_ptr<Radar> radar_;
  KalmanFilter ekf_;

 private:
  // check whether the tracking toolbox was initialized or not (first
  // measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
};

#endif /* FusionEKF_H_ */
