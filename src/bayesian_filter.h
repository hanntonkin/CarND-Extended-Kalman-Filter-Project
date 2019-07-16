#ifndef BAYESIAN_FILTER_H_
#define BAYESIAN_FILTER_H_

#include <iostream>
#include <memory>

#include "Eigen/Dense"

using std::cout;
using std::endl;

typedef double Scalar;
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

class Process {
 public:
  Scalar dt_;

  // state vector
  Vector x_;

  // state covariance matrix
  Matrix P_;

  // state transition matrix
  Matrix F_;

  // process covariance matrix
  Matrix Q_;

 public:
  virtual const Matrix F() = 0;
  virtual const Matrix P() = 0;
  virtual const Matrix Q() = 0;
};

class Measurement {
 public:
  // measurement
  Vector y_;

  // measurement matrix
  Matrix H_;

  // measurement covariance matrix
  Matrix R_;

 public:
  // Measurement(Matrix H): H_(H){};
  virtual Matrix H() = 0;
  virtual Matrix R() = 0;
  virtual Vector y() = 0;
};

class BayesianFilter {
 public:
  virtual void update(const std::shared_ptr<Measurement> &measurement) = 0;
  virtual void predict() = 0;
};

#endif
