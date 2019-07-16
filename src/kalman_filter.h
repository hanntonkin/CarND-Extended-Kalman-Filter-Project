#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include <memory>

#include "Eigen/Dense"
#include "bayesian_filter.h"

class KalmanFilter : public BayesianFilter {
 protected:
  std::shared_ptr<Process> process_;

 public:
  void init(const std::shared_ptr<Process> &process) { process_ = process; }

  virtual void predict() {
    auto F = process_->F();
    auto P = process_->P();
    auto Q = process_->Q();

    process_->x_ = F * process_->x_;
    process_->P_ = F * P * F.transpose() + Q;
  }

  virtual void update(const std::shared_ptr<Measurement> &measurement) {
    auto P = process_->P();
    auto H = measurement->H();
    auto R = measurement->R();
    auto y = measurement->y();

    auto T = H * P * H.transpose() + R;
    auto K = P * H.transpose() * T.inverse();

    process_->x_ = process_->x_ + (K * y);
    process_->P_ = P - K * H * P;
  }

  Vector x() { return process_->x_; }
};

#endif /* KALMAN_FILTER_H_ */
