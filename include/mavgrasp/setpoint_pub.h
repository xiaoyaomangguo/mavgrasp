#ifndef SETPOINT_PUB_H
#define SETPOINT_PUB_H

#include <iostream>
#include <fstream>          // file I/O support
#include <string>
#include <vector>
#include <cstdlib>          // support for exit()
#include <cmath>

#include <eigen3/Eigen/Eigen>

#include "ros/ros.h"
#include "geometry_msgs/PoseStamped.h"
#include "ros/time.h"

using namespace std;
using namespace Eigen;

const double AllowError = 0.1;

class SetPoint
{
public:
    SetPoint();
    virtual ~SetPoint();

    void poseCallBack(const geometry_msgs::PoseStampedConstPtr& msg);
    void start();            //make the mav fly to the starting points

    void diffFlatness(double t);
    void generatePosPub();
    void generatePosPub(Vector3d&  r_q);
    void generateAttPub(geometry_msgs::PoseStamped local_att);

private:
    ros::NodeHandle node_;
    ros::Time       time_now;
    ros::Time       time_prev;
    double          time_after_start;

    ros::Subscriber local_pos_sub_;
    ros::Publisher  pos_sp_pub_;
    ros::Publisher  att_sp_pub_;

    //differential flatness output
    Vector3d    r_q;
    double      beta;
    double      theta;
    double      u1;
    double      u3;
    //


    // the coefs
    MatrixXd    co_x_;
    MatrixXd    co_z_;
    MatrixXd    co_beta_;
    MatrixXd    co_yaw_;
    MatrixXd    co_d1x_;        //first Derivative
    MatrixXd    co_d1z_;
    MatrixXd    co_d1beta_;
    MatrixXd    co_d1yaw_;
    MatrixXd    co_d2x_;        //Second Derivative
    MatrixXd    co_d2z_;
    MatrixXd    co_d2beta_;
    MatrixXd    co_d2yaw_;
    MatrixXd    co_d3x_;
    MatrixXd    co_d3z_;
    MatrixXd    co_d3beta_;
    MatrixXd    co_d3yaw_;
    MatrixXd    co_d4x_;
    MatrixXd    co_d4z_;
    MatrixXd    co_d4beta_;
    MatrixXd    co_d4yaw_;
    MatrixXd    waypoints_;     //5*4
    VectorXd    timepoints_;    //5d

    //
    void converseQuatToEuler(int flag, double roll, double pitch, double yaw, geometry_msgs::Quaternion q);
    //read file to get coefficients
    ifstream infile;
    void readCoef();
    void readCoeftwp();     //read coefs of timepoints and waypoints
    void readCoefd(const char *filename, int n_rows, int n_cols, MatrixXd& mx, MatrixXd &mz, MatrixXd &mb, MatrixXd &my );
};

#endif // SETPOINT_PUB_H


