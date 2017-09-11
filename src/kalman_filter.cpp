#include "kalman_filter.h"
#include <iostream>

#define PI (3.141592653589793)

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
{
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() 
{
	// Predict the state
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
	
	// update the state by using Kalman Filter equations
	VectorXd z_predicted = this->H_ * this->x_;
	VectorXd y = z - z_predicted;
	MatrixXd Ht = this->H_.transpose();
	MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = this->P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	this->x_ = this->x_ + (K * y);
	long x_size = this->x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	this->P_ = (I - K * this->H_) * this->P_;
	
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
	// update the state by using Extended Kalman Filter equations
	VectorXd z_predicted = this->ConvertCartesianToRadar( this->x_ );
		
	VectorXd y = z - z_predicted;
	y(1) = this->LimitAngle(y(1));
	
	MatrixXd Ht = this->H_.transpose();
	MatrixXd S = this->H_ * this->P_ * Ht + this->R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = this->P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	this->x_ = this->x_ + (K * y);
	long x_size = this->x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	this->P_ = (I - K * this->H_) * this->P_;

}


VectorXd KalmanFilter::ConvertCartesianToRadar(const VectorXd &cartesian)
{
	VectorXd radar_measurement(3);
	
	float px = cartesian(0);
	float py = cartesian(1);
	float vx = cartesian(2);
	float vy = cartesian(3);
	
	radar_measurement(0) = sqrt( px*px + py*py );
	radar_measurement(1) = atan2( py, px );
	
	if( fabs(radar_measurement(0)) < 1.0e-6 ){
		radar_measurement(2) = 0.0;
		std::cout << "ConvertCartesianToRadar: Near zero division (by " << radar_measurement(0) << ")! - x = " << cartesian.transpose() << std::endl;
	}
	else{
		radar_measurement(2) = (px*vx + py*vy) / radar_measurement(0);		
	}
	
	return radar_measurement;
}

VectorXd KalmanFilter::ConvertRadarToCartesian(const VectorXd &radar_measurement)
{
	// This function will calculate the cartesian coordinates from the radar coordinate.
	// However, the speed information in radar coordinates is incomplete, and only 
	// the position will be returned (the cartesian speeds will be 0).
	VectorXd cartesian(4);
	
	float rho = radar_measurement(0);
	float phi = radar_measurement(1);
	
	cartesian(0) = rho * cos(phi);
	cartesian(1) = rho * sin(phi);
	cartesian(2) = 0.0;
	cartesian(3) = 0.0;
	
	return cartesian;
}

float KalmanFilter::LimitAngle(float angle)
{
	// Limit the angular error to +/- PI
	while( angle > PI ){
		angle -= 2*PI;
	}
	while( angle < -PI ){
		angle += 2*PI;
	}	
	return angle;
}

