/**
* Reference Udactity
**/

#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 2.1;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	/**
	TODO:

	Complete the initialization. See ukf.h for other member properties.

	Hint: one or more values initialized above might be wildly off...
	*/
	
	is_initialized_ = false;
	use_laser_ = true;
	use_radar_ = true;
	
	n_x_ = 5;
	n_aug_ = n_x_ + 2;
	lambda_ = 3 - n_aug_;
	
	std_radr_ = 0.3;
	std_radphi_ = 0.3;
	std_radrd_ = 0.3;
	
	previous_timestamp_ = 0;
	
	// initial covariance matrix
	P_ = MatrixXd::Zero(5, 5);
	//P_ << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1;
	P_ = MatrixXd::Identity(5,5);
	
	NIS_radar_ = 0.0;
	NIS_laser_ = 0.0;
	
	Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
	
	//set vector for weights
	weights_ = VectorXd(2*n_aug_+1);
	double weight_0 = lambda_/(lambda_+n_aug_);
	weights_(0) = weight_0;
	for (int i=1; i<2*n_aug_+1; i++) {  
		double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = weight;
	}
	
	x_ << 0, 0, 0, 0, 0;
	
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
	
	R_radar_ = MatrixXd(3,3);
	R_radar_ <<    std_radr_*std_radr_, 0, 0,
	0, std_radphi_*std_radphi_, 0,
	0, 0,std_radrd_*std_radrd_;
		
	R_laser_ = MatrixXd(2, 2);
	R_laser_ << pow(std_laspx_, 2), 0, 0, pow(std_laspy_, 2);	
	

	
	
}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


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
		// first measurement
		//cout << "EKF: " << endl;
		
	  
		float px = 0.0;
		float py = 0.0;
		

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			
			float ro = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			float ro_dot = meas_package.raw_measurements_[2];
		
			double px = ro * cos(phi);
			double py = ro * sin(phi);
			double vx = ro_dot * cos(phi);
			double vy = ro_dot * sin(phi);
			
			//check for valid data
			if ( fabs(px) < 0.001 && fabs(py) < 0.001 )  
			{
				px = 0.001;    
				py = 0.001;   
			}
			//x_ << px, py, sqrt(vx * vx + vy * vy), 0, 0;
			x_ << px, py, abs(ro_dot), 0, 0;
			
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			x_[0] = meas_package.raw_measurements_[0];
			x_[1] = meas_package.raw_measurements_[1];
			x_[2] = 0; 
			x_[3] = 0; 
		}
	
		//previous_timestamp_ = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		previous_timestamp_ = meas_package.timestamp_;
		return;
	}
	
	
	double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	while (delta_t > 0.1) {
		const double dt = 0.05;
		Prediction(dt);
		delta_t -= dt;
	}	
	

	// Prediction
	Prediction(delta_t);

	// Update
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		//std::cout << "rlidar" << std::endl;
		UpdateRadar(meas_package);
	} else {
		//std::cout << "radar" << std::endl;
		UpdateLidar(meas_package);
		
		
	}

	previous_timestamp_ = meas_package.timestamp_;
	
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/
		
	//std::cout << "prediciton" << std::endl;
	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_	;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();
	
	GenerateSigmaPoints(delta_t, L);
	
	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	
	

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

	
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//normalize
		NormalizeAngle(x_diff(3));

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}
		
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
		
	//::cout << "lidar" << std::endl;
		
	// create matrix for sigma points in measurement space
	int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
		// measurement model
		Zsig(0, i) = Xsig_pred_(0, i); // px
		Zsig(1, i) = Xsig_pred_(1, i); // py
	}

	// mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	// create cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// add measurement noise covariance matrix
	
	S = S + R_laser_;

	// update
	VectorXd z = VectorXd::Zero(n_z);
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

	// Kalman gain K;
	MatrixXd K = MatrixXd::Zero(n_x_, n_z);
	K = Tc * S.inverse();

	// residual
	VectorXd z_diff = VectorXd::Zero(n_z);
	z_diff = z - z_pred;

	// new estimate state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

	// calculate the lidar NIS
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
		

}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use radar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the radar NIS.
	*/


	//std::cout << "radar" << std::endl;
	//set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);


	//sigma points to measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// measurement model
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);
		const double v  = Xsig_pred_(2,i);
		const double yaw = Xsig_pred_(3,i);
		const double v1 = cos(yaw)*v;
		const double v2 = sin(yaw)*v;
		
		//check for valid data
		if(fabs(px) < 0.0001 && fabs(py) < 0.0001)
		{
			px = 0.0001;
			py = 0.0001;
		}
		
			
		Zsig(0,i) = sqrt(px * px + py * py);                        
		Zsig(1,i) = atan2(py,px);                                 
		Zsig(2,i) = (px * v1 + py * v2 ) / sqrt(px * px + py * py);   
	}
	
	// mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
	Tc.fill(0.0);

	//measurement covariance matrix S
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//normalize
		NormalizeAngle(z_diff(1));

		S = S + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//normalize
		NormalizeAngle(x_diff(3));

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	
	S = S + R_radar_;
	
	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z = VectorXd(n_z);
	z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;

	//normalize
	NormalizeAngle(z_diff(1));

	// new estimate state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();
		
	// calculate the radar NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
		
}

void UKF::GenerateSigmaPoints(double delta_t, MatrixXd L)
{
	
	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;
	
	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	
	//create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
	}
		
	//predict sigma points
	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw = Xsig_aug(3,i);
		double yawd = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}

}

void UKF::NormalizeAngle(double& phi)
{
	phi = atan2(sin(phi), cos(phi));
}
