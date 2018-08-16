#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);
  
  static MatrixXd calculateSigmaPoint(VectorXd x, MatrixXd P, double std_a, double std_yawdd);
  
  static MatrixXd calculateSigmaPoint(VectorXd x, MatrixXd P);

  static float normalizeAngle(float angle);
  
  static float percentageAbove(vector<float> values, float threshold);
};

#endif /* TOOLS_H_ */