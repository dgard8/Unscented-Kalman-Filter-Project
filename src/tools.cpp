#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd residuals(4);
    residuals << 0,0,0,0;
    
    
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return residuals;
    }
    
    for(int i=0; i < estimations.size(); ++i){
        VectorXd diff = estimations[i] - ground_truth[i];
        residuals += VectorXd(diff.array()*diff.array());
    }
    
    residuals /= estimations.size();
    residuals = residuals.array().sqrt();
    
    return residuals;
}

MatrixXd Tools::calculateSigmaPoint(VectorXd x, MatrixXd P, double std_a, double std_yawdd){
    int n_x = x.size();
    int n_aug = n_x + 2;
    
    VectorXd x_aug = VectorXd(n_aug);
    x_aug.fill(0.0);
    x_aug.head(n_x) = x;
    
    MatrixXd P_aug = MatrixXd(n_aug, n_aug);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x,n_x) = P;
    P_aug(n_x,n_x) = std_a*std_a;
    P_aug(n_x+1,n_x+1) = std_yawdd*std_yawdd;
    
    return calculateSigmaPoint(x_aug, P_aug);
}

MatrixXd Tools::calculateSigmaPoint(VectorXd x, MatrixXd P){
    int n_x = x.size();
    double lambda = 3 - n_x;
    
    MatrixXd sigmaPoints = MatrixXd(n_x, 2 * n_x + 1);
    MatrixXd A = P.llt().matrixL(); //square root of P
    
    sigmaPoints.col(0) = x;
    
    for (int i=1; i<(n_x)+1; i++){
        VectorXd offset = sqrt(lambda + n_x)*A.col(i-1);
        sigmaPoints.col(i) = x + offset;
        sigmaPoints.col(i+n_x) = x - offset;
    }
    
    return sigmaPoints;
}

float Tools::normalizeAngle(float angle){
    while (angle> M_PI) angle-=2.*M_PI;
    while (angle<-M_PI) angle+=2.*M_PI;
    return angle;
}

float Tools::percentageAbove(vector<float> values, float threshold){
    int count=0;
    for (int i=0; i<values.size(); i++){
        if (values[i] > threshold)
            count++;
    }
    return (float)count/values.size();
}