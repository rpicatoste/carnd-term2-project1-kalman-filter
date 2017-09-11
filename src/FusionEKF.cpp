#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <assert.h> 

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() 
{
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

	// Measurement matrix - laser
	H_laser_ << 1,0,0,0,
				0,1,0,0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{

	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (!is_initialized_) {
		FirstMeasurement(measurement_pack);	
		return;
	}
	 
	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	
	// If the simulator is restarted, the delta time will be negative, and we need to consider the measurement
	// as the initialization one.
	if(dt < 0.0){
		FirstMeasurement(measurement_pack);	
		return;
	}
	
	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/
	// Update the state transition matrix F according to the new elapsed time.
    //  - Time is measured in seconds.
	this->ekf_.F_(0, 2) = dt;
	this->ekf_.F_(1, 3) = dt;
	
	// Update the process noise covariance matrix.
	// * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	const float noise_ax = 9, noise_ay = 9;
	this->ekf_.Q_ <<  	dt_4/4*noise_ax, 0, 				dt_3/2*noise_ax,	0,
						0,				 dt_4/4*noise_ay, 	0, 					dt_3/2*noise_ay,
						dt_3/2*noise_ax, 0, 				dt_2*noise_ax, 		0,
						0, 				 dt_3/2*noise_ay, 	0, 					dt_2*noise_ay;

	this->ekf_.Predict();

	/*****************************************************************************
	 *  Update
	 ****************************************************************************/

	// Use the sensor type to perform the update step.
	// Update the state and covariance matrices.
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates		
		
		this->Hj_ = this->ekf_.CalculateJacobian( this->ekf_.x_ );
		this->ekf_.R_ = this->R_radar_;
		this->ekf_.H_ = this->Hj_;
		
		std::cout << "Radar update! = " << measurement_pack.raw_measurements_.transpose() << std::endl;
		
		this->ekf_.UpdateEKF( measurement_pack.raw_measurements_ );
		
	} else {
		// Laser updates		
		this->ekf_.R_ = this->R_laser_;
		this->ekf_.H_ = this->H_laser_;
		
		cout << "Laser update! = " << measurement_pack.raw_measurements_.transpose() << endl;
		
		this->ekf_.Update( measurement_pack.raw_measurements_ );
	}

	previous_timestamp_ = measurement_pack.timestamp_;
	
	// Print the output
	cout << "dt = " << dt << endl;
	cout << "x_ = " << ekf_.x_.transpose() << endl;
	cout << endl;
}



void FusionEKF::FirstMeasurement(const MeasurementPackage &measurement_pack) 
{
	/**
	 * Initialize the state ekf_.x_ with the first measurement.
	 * Create the covariance matrix.
	 * Remember: you'll need to convert radar from polar to cartesian coordinates.
	 */
	
	MatrixXd P_in(4,4); 
	P_in << 1, 0, 0,    0,
			0, 1, 0,    0,
			0, 0, 1000, 0,
			0, 0, 0,    1000;
				
	MatrixXd F_in(4,4);
	F_in << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;
			
	MatrixXd Q_in(4,4);
	Q_in << 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0;
		
	// Initialize state.
	VectorXd x_in(4); 
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Convert radar from polar to cartesian coordinates and initialize state.
		x_in = this->ekf_.ConvertRadarToCartesian( measurement_pack.raw_measurements_ );
		 
		ekf_.Init(x_in, P_in, F_in, this->Hj_, this->R_radar_, Q_in);
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		
		x_in << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;

		ekf_.Init(x_in, P_in, F_in, this->H_laser_, this->R_laser_, Q_in);
	}
	else {
		// Wrong sensor
		assert(0);
	}

	previous_timestamp_ = measurement_pack.timestamp_;
	
	is_initialized_ = true;
	
	cout << "Kalman initialized, x = " << ekf_.x_.transpose() << endl;
	
	std::cout << "ekf_.F_:\n " << this->ekf_.F_ << std::endl; 
	std::cout << "ekf_.P_:\n " << this->ekf_.P_ << std::endl; 
	std::cout << "ekf_.Q_:\n " << this->ekf_.Q_ << std::endl; 
}