#include "kalman_filter.h"
#include "tools.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {
}

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

    //cout << "Predict Kalmanfilter: " << endl;
    //cout << "x_ = " << x_ << endl;
    //cout << "P_ = " << P_ << endl;
	//cout << "F_ = " << P_ << endl;

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    //cout << "Update Kalmanfilter: " << endl;
    //cout << "x_ = " << x_ << endl;
    //cout << "z = " << z << endl;
    //cout << "P_ = " << P_ << endl;
	//cout << "H_laser = " << H_ << endl;
	//cout << "R_laser = " << R_ << endl;

	//Formulas for Kalmanfilter
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	 TODO:
	 * update the state by using Extended Kalman Filter equations
	 */

    //cout << "Update Kalmanfilter EKF: " << endl;
    //cout << "x_ = " << x_ << endl;
    //cout << "z = " << z << endl;
    //cout << "P_ = " << P_ << endl;
	//cout << "H_radar = " << H_ << endl;
    //cout << "h(x) = " << h(x_) << endl;
	//cout << "R_radar = " << R_ << endl;

	// Instead of having initialize H_ with H_laser in EKF use Hj jacobian before calling make sure Hj is calculated correctly
	// Instead of using H*x for EKF we are using h(x)
	// Everything else is like with a Kalmanfilter
	VectorXd y = z - h(x_);
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

//calculate h(x) according to formula
VectorXd KalmanFilter::h(const VectorXd &x) {

	// extract position and velocity
	float px = x(0);
	float py = x(1);
	float vx = x(2);
	float vy = x(3);

	// make sure no division by zero is happening
	// Calculation according to formula for EKF
	float rho = sqrt(px*px + py*py);
    float theta = 0;
    float rho_dot = 0;
	if(rho>0.001){
        theta = atan(py/px);
        rho_dot = (px*vx + py*vy) / rho;
	}

	VectorXd hx = VectorXd(3);
	hx << rho, theta, rho_dot;

	return hx;
}
