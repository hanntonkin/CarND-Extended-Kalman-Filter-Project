#include "FusionEKF.h"


using std::cout;
using std::endl;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  car_ = std::make_shared<Car>();
  laser_ = std::make_shared<Laser>();
  radar_ = std::make_shared<Radar>();
  ekf_.init(car_);
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
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian
    coordinates.
    */
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "EKF : First measurement RADAR" << endl;
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];  // range
      double phi = measurement_pack.raw_measurements_[1];  // bearing
      double rho_dot =
          measurement_pack.raw_measurements_[2];  // velocity of rho
      // Coordinates convertion from polar to cartesian
      double x = rho * cos(phi);
      if (x < 0.0001) {
        x = 0.0001;
      }
      double y = rho * sin(phi);
      if (y < 0.0001) {
        y = 0.0001;
      }
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      auto measurement = VectorXd(4);
      measurement << x, y, vx, vy;
      car_->init(measurement);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // No velocity and coordinates are cartesian already.
      cout << "EKF : First measurement LASER" << endl;
      auto measurement = VectorXd(4);
      measurement << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1], 0, 0;
      car_->init(measurement);
    }

    // Saving first timestamp in seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // // State transition matrix update
  car_->dt(dt);
  ekf_.predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    // radar_->H_ = tools.CalculateJacobian(ekf_.x());
    radar_->measure(measurement_pack.raw_measurements_, ekf_.x());
    ekf_.update(radar_);

  } else {
    // Laser updates
    Vector z(2);
    z << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1];
    laser_->measure(z, ekf_.x());
    ekf_.update(laser_);
  }
}
