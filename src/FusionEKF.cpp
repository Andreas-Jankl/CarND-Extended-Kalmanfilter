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

  // Initialize Matrices
  is_initialized_ = false;

  previous_timestamp_ = 0;

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  R_laser_ << 0.0225, 0,
              0, 0.0225;

  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;


  Q_ = MatrixXd(4, 4);
  F_ = MatrixXd(4, 4);
  P_ = MatrixXd(4, 4);
  x_ = VectorXd(4);
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
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    cout << "Initialise EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

            //Initialize the states for a first measurement with a radar
            //This means converting coordinates from polar to cartesian
            float rho = measurement_pack.raw_measurements_[0]; // Range - radial distance from origin
			float phi = measurement_pack.raw_measurements_[1]; // Bearing - angle between rho and x
			float rho_dot = measurement_pack.raw_measurements_[2]; // Radial Velocity - change of p (range rate)
			float x = rho * cos(phi);
			float y = rho * sin(phi);
			float vx = rho_dot * cos(phi);
			float vy = rho_dot * sin(phi);

			x_ << x, y, vx , vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

            //Initialize the states for a first measurement with a laser
            //Just take the position but and leave the speeds to zero
            x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    }

    //Set the state vektor in the Kalmanfilter
    ekf_.x_=x_;

    //Set the process covariance matrix in the Kalmanfilter
    P_ << 1.0, 0, 0, 0,
          0, 1.0, 0, 0,
          0, 0, 1.0, 0,
          0, 0, 0, 1.0;

    ekf_.P_ = P_;

    //cout << "x_initialize = " << ekf_.x_ << endl;
    //cout << "P_initialize = " << ekf_.P_ << endl;

    //Set the previous timestamp as the initial one. and set the filter to be intialized
    previous_timestamp_ = measurement_pack.timestamp_;
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

   	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

    cout << "Timestamp dt = " << dt << endl;

    //The predict step is only called in case the measurement is newer than the last one.
    //Else when we got two measurements at the exact same time we'll just call Predict once and call the according measurement update for the second measurement directly
    if(dt!=0){

    cout << "Predict" << endl;

	//set the transition matrix F_ with current dt
	F_ << 1, 0, dt, 0,
		  0, 1, 0, dt,
		  0, 0, 1, 0,
		  0, 0, 0, 1;

    //set F_ in the Kalmanfilte where its used
	ekf_.F_ = F_;

	//set the acceleration noise components according to hint above
	double noise_ax = 9.0;
	double noise_ay = 9.0;

	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//calculate and set the process covariance matrix Q
	Q_ <<   dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


	ekf_.Q_ = Q_;

	//Predict step in the Kalmanfilter
	ekf_.Predict();

    }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        // Radar update. Prepare measurement
		VectorXd z = VectorXd(3);
		z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];

        // Calculate the Jacobian and init the Kalmanfilter for EKF Update.
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);

		// If there as a calculated Jacobian call the update step
		if (!Hj_.isZero(0)){
			ekf_.UpdateEKF(z);}
  } else {
        // Laser updates. Prepare measurement
        VectorXd z = VectorXd(2);
        z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];

        // Init the Kalmanfilter for standard Update and call Update
        ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
        ekf_.Update(z);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
