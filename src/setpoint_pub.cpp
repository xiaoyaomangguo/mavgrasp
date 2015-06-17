
#include "../include/mavgrasp/setpoint_pub.h"
#include <geometry_msgs/PoseStamped.h>

using namespace std;

const char coef_f0[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/coefd0.txt";
const char coef_f1[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/coefd1.txt";
const char coef_f2[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/coefd2.txt";
const char coef_f3[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/coefd3.txt";
const char coef_f4[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/coefd4.txt";
const char coef_twp[]= "/home/congleetea/catkin_ws/src/mavgrasp/src/config/Twaypoints.txt";
const int  n_segment = 4;
const int  n_order   = 7;
const int m_q   = 2;      //2kg
const int m_g   = 0.2;    //0.2kg
const int m_s   = (m_g + m_q);
const int L     = 0.1;      //the length of grasper
const int G     = 9.8;
const int J_g   = 1;
const int J_q   = 1;

SetPoint::SetPoint():
    node_("~"),
    time_after_start(0),
    time_prev(0),
    start_grasp(false),
    end_grasp(false)
{
    start();
}

SetPoint::~SetPoint(){}

void SetPoint::poseCallBack(const geometry_msgs::PoseStampedConstPtr &msg)
{

    if( fabs(msg->pose.position.x - waypoints_(0,0)) < AllowError
        && fabs(msg->pose.position.y)                < AllowError
        && fabs(msg->pose.position.z - waypoints_(0,1)) < AllowError  // not the first waypoint, continue to the first point
        && !start_grasp)
    {
        start_grasp = true;  //begin to grasp
        ROS_INFO("Begin to grasp!!!");
    }
    if(start_grasp == false && end_grasp == false) //not to the first point, do not to grasp
        {
            ROS_INFO("Waiting for grasping!!!");
            generatePosPub(0);
        }
        if(start_grasp) // begin to grasp
        {
            time_now = ros::Time::now();
            time_after_start += (time_now.toNSec() - time_prev.toNSec())/1e+9;
            diffFlatness(time_after_start);
            generatePosPub(r_q);
            generateAttPub(*msg);
            time_prev = time_now;
            if(time_after_start > timepoints_(4)) //if end grasping
            {
                end_grasp == true;
                start_grasp == false;
            }
        }

        if(end_grasp)//end to grasp
        {
            ROS_INFO("Grasp End!!!");
            generatePosPub(5);    //stop in the end point
        }
}

void SetPoint::start()
{
    readCoef();
    pos_sp_pub_ = node_.advertise<geometry_msgs::PoseStamped>("/mavros/setpoint_position/local",100);
    att_sp_pub_ = node_.advertise<geometry_msgs::PoseStamped>("/mavros/setpoint_attitude/attitude",100);
    local_pos_sub_ = node_.subscribe("/mavros/local_position/local",100,&SetPoint::poseCallBack, this);  //callback for class methods should add the third parameter point "this"
}

void SetPoint::diffFlatness(double t)
{
    Vector3d r_g;
    Vector3d dr_q;
    Vector3d dr_g;
    Vector3d d2r_q;
    Vector3d d2r_g;
    Vector3d d3r_q;
    Vector3d d3r_g;
    Vector3d d4r_q;
    Vector3d d4r_g;
    Vector3d r_s;
    Vector3d dr_s;
    Vector3d d2r_s;
    Vector3d d3r_s;
    Vector3d d4r_s;

    VectorXd timestamp(8);
    int n_segment;
    timestamp << 1 , t , pow(t,2) , pow(t,3) ,pow(t,4) , pow(t,5) , pow(t,6) , pow(t,7);
    if(t >= timepoints_(0) && t < timepoints_(1)) n_segment = 0;
    if(t >= timepoints_(1) && t < timepoints_(2)) n_segment = 1;
    if(t >= timepoints_(2) && t < timepoints_(3)) n_segment = 2;
    if(t >= timepoints_(3) && t < timepoints_(4)) n_segment = 3;
    double x = co_x_.row(n_segment) * timestamp;
    double z = co_z_.row(n_segment) * timestamp;
    beta = co_beta_.row(n_segment) * timestamp;

    double dx = co_d1x_.row(n_segment) * timestamp.segment(0, 7);
    double dz = co_d1z_.row(n_segment) * timestamp.segment(0, 7);
    double dbeta = co_d1beta_.row(n_segment) * timestamp.segment(0,7);

    double d2x = co_d2x_.row(n_segment) * timestamp.segment(0, 6);
    double d2z = co_d2z_.row(n_segment) * timestamp.segment(0, 6);
    double d2beta = co_d2beta_.row(n_segment) * timestamp.segment(0,6);

    double d3x = co_d3x_.row(n_segment) * timestamp.segment(0, 5);
    double d3z = co_d3z_.row(n_segment) * timestamp.segment(0, 5);
    double d3beta = co_d3beta_.row(n_segment) * timestamp.segment(0,5);

    double d4x = co_d4x_.row(n_segment) * timestamp.segment(0, 4);
    double d4z = co_d4z_.row(n_segment) * timestamp.segment(0, 4);
    double d4beta = co_d4beta_.row(n_segment) * timestamp.segment(0,4);

    Vector3d cs;
    cs << cos(beta),
          0 ,
          -sin(beta);
    Vector3d dcs;
    dcs << -sin(beta) * dbeta,
          0 ,
          -cos(beta) * dbeta;
    Vector3d d2cs;
    dcs << -cos(beta)*pow(dbeta,2) - sin(beta) * d2beta,
            0 ,
            sin(beta)*pow(dbeta,2) - cos(beta) * d2beta;
    Vector3d d3cs;
    d2cs << sin(beta) * pow(dbeta,3) - 3*cos(beta)*dbeta*d2beta - sin(beta)*d3beta,
            0,
            cos(beta) * pow(dbeta,3) + 3*sin(beta)*dbeta*d2beta - cos(beta)*d3beta;
    Vector3d d4cs;
    d4cs << cos(beta)*pow(dbeta,4)+4*sin(beta)*pow(dbeta,2)*d2beta-3*cos(beta)*pow(d2beta,2)-4*cos(beta)*dbeta*d3beta-sin(beta)*d4beta,
            0,
            -sin(beta)*pow(dbeta,4)+4*cos(beta)*pow(dbeta,2)*d2beta+3*sin(beta)*pow(d2beta,2)+4*sin(beta)*dbeta*d3beta-cos(beta)*d4beta;

    r_q << x, 0, z;
    r_g = r_q + L * cs;
    r_s = (m_q * r_q + m_g * r_g) / m_s;
    dr_q << dx, 0, dz;
    dr_g = dr_q + L * dcs;
    dr_s = (m_q * dr_q + m_g * dr_g) / m_s;
    d2r_q << d2x, 0, d2z;
    d2r_g = d2r_q + L * d2cs;
    d2r_s = (m_q * d2r_q + m_g * d2r_g) / m_s;
    d3r_q << d3x, 0, d3z;
    d3r_g = d3r_q + L * d3cs;
    d3r_s = (m_q * d3r_q + m_g * d3r_g) / m_s;
    d4r_q << d4x, 0, d4z;
    d4r_g = d4r_q + L * d4cs;
    d3r_s = (m_q * d3r_q + m_g * d3r_g) / m_s;

    Vector3d e1, e2, e3, b1, b2, b3;
    e1 << 1,0,0;
    e2 << 0,1,0;
    e3 << 0,0,1;

    b2 = e2;
    b3 = (d2r_s + G * e3)/(d2r_s + G * e3).norm();
    b1 = b2.cross(b3);
    double e3b3_cos = e3.dot(b3);
    double e3b3_sin = (e3.cross(b3)).norm();
    theta = atan2(e3b3_sin, e3b3_cos);
    u1 = m_s * (d2r_s + G * e3).norm();
    double du1 = m_s * b3.dot(d3r_s);
    double dtheta = m_s * (b1.dot(d3r_s));
    double d2theta = (m_s * b1.dot(d4r_s) - 2 * du1 * dtheta);
    double d2u1 = b3.dot(m_s * d4r_s) + dtheta * dtheta * u1;
    double tao = J_g * d2beta - L * m_g * (d2r_g(0) * sin(beta) + (d2r_g(2) + G) * cos(beta));

    u3 = d2theta * J_q + tao;
    //output: r_q, beta, theta, u1, u3
}

void SetPoint::generatePosPub(int n_point)
{
    geometry_msgs::PoseStamped msg;
    msg.pose.position.x = waypoints_(n_point,0);
    msg.pose.position.y = 0;
    msg.pose.position.z = waypoints_(n_point,1);
    pos_sp_pub_.publish(msg);
}

void SetPoint::generatePosPub(Vector3d &r_q)
{
    geometry_msgs::PoseStamped msg;
    msg.pose.position.x = r_q(0);
    msg.pose.position.y = r_q(1);
    msg.pose.position.z = r_q(2);
    pos_sp_pub_.publish(msg);
}

void SetPoint::generateAttPub(geometry_msgs::PoseStamped local_att)
{
    double yaw_const = 0.0;
    converseQuatToEuler(0, NULL,NULL,yaw_const, local_att.pose.orientation);
    geometry_msgs::PoseStamped msg;
    converseQuatToEuler(1, 0.0, theta, yaw_const, msg.pose.orientation);
    att_sp_pub_.publish(msg);
}


void SetPoint::converseQuatToEuler(int flag, double roll, double pitch, double yaw, geometry_msgs::Quaternion q)
{
    if(flag == 1)       //euler to q
    {
        double cosPhi_2 = cos(roll / 2.0);
        double sinPhi_2 = sin(roll / 2.0);
        double cosTheta_2 = cos(pitch / 2.0);
        double sinTheta_2 = sin(pitch / 2.0);
        double cosPsi_2 = cos(yaw / 2.0);
        double sinPsi_2 = sin(yaw / 2.0);

        q.w = cosPhi_2 * cosTheta_2 * cosPsi_2 + sinPhi_2 * sinTheta_2 * sinPsi_2;
        q.x = sinPhi_2 * cosTheta_2 * cosPsi_2 - cosPhi_2 * sinTheta_2 * sinPsi_2;
        q.y = cosPhi_2 * sinTheta_2 * cosPsi_2 + sinPhi_2 * cosTheta_2 * sinPsi_2;
        q.z = cosPhi_2 * cosTheta_2 * sinPsi_2 - sinPhi_2 * sinTheta_2 * cosPsi_2;
    }
    else if(flag == 0)  //converse quaternion to euler
    {
        yaw = atan2(2.0 * (q.w * q.x + q.y * q.z), 1.0 - 2.0 * (q.x * q.x + q.y * q.y));
        pitch = asin(2.0 * (q.w * q.y - q.z * q.x));
        roll = atan2(2.0 * (q.w * q.z + q.x * q.y), 1.0f - 2.0f * (q.y * q.y + q.z * q.z));
    }
}

void SetPoint::readCoef()
{
    readCoeftwp();
    readCoefd(coef_f0, n_segment,n_order+1,co_x_, co_z_, co_beta_, co_yaw_);
    readCoefd(coef_f1, n_segment,n_order,co_d1x_, co_d1z_, co_d1beta_, co_d1yaw_);
    readCoefd(coef_f2, n_segment,n_order-1,co_d2x_, co_d2z_, co_d2beta_, co_d2yaw_);
    readCoefd(coef_f3, n_segment,n_order-2,co_d3x_, co_d3z_, co_d3beta_, co_d3yaw_);
    readCoefd(coef_f4, n_segment,n_order-3,co_d4x_, co_d4z_, co_d4beta_, co_d4yaw_);
}

void SetPoint::readCoeftwp()
{
    MatrixXd    waypoints(5,4);
    VectorXd    timepoints(5);
    double  val;
    int     count = 0;
    infile.open(coef_twp);
    if(!infile.is_open())
    {
        cout << "Unable to open file:"<< coef_twp << "Program terminating.\n" << endl;
        exit(EXIT_FAILURE);
    }
    infile >> val;
    int ti=0;
    while(infile.good())
    {
        if(count % 5 == 0){ //this is timepoints
            timepoints(ti++) = val;
        }
        else {

            waypoints(ti-1,(count-ti)%4) = val;
        }
        ++count;
        infile >> val;
    }
    if(infile.eof())
        cout << "End of file reached.\n";
    infile.close();
    cout << "Read " << count << " waypoints numbers totally.\n";
    timepoints_ = timepoints;
    waypoints_ = waypoints;
}

void SetPoint::readCoefd(const char *filename, int n_rows, int n_cols, MatrixXd &mx, MatrixXd &mz, MatrixXd &mb, MatrixXd &my)
{
    MatrixXd    co_x(n_rows,n_cols);  //here 4 sgement *  order
    MatrixXd    co_z(n_rows,n_cols);
    MatrixXd    co_beta(n_rows,n_cols);
    MatrixXd    co_yaw(n_rows,n_cols);
    double  val;
    infile.open(filename);
    if(!infile.is_open()){
        cout << "Unable to open file:"<< filename << "Program terminating.\n" << endl;
        exit(EXIT_FAILURE);
    }
    int     count = 0;
    infile >> val;
    while(infile.good())
    {
        if(count >= 0 && count < n_rows*n_cols){
            co_x(count / n_cols, count % n_cols) = val;
        }
        if(count >= n_rows*n_cols && count < 2*n_rows*n_cols){
            co_z((count-n_rows*n_cols) / n_cols, count % n_cols) = val;
        }
        if(count >= 2*n_rows*n_cols && count < 3*n_rows*n_cols){
            co_beta((count-2*n_rows*n_cols) / n_cols, count % n_cols) = val;
        }
        if(count >= 3*n_rows*n_cols && count < 4*n_rows*n_cols){
            co_yaw((count-3*n_rows*n_cols) / n_cols,count % n_cols) = val;
        }
        ++count;
        infile >> val;
    }
    if(infile.eof())  cout << "End of file reached.\n";
    infile.close();
    cout << "Read " << count << " coef numbers totally.\n";
    mx = co_x;
    mz = co_z;
    mb = co_beta;
    my = co_yaw;
    cout << "The coefs of " << 8 - n_cols << "-th order = \n" << mx << endl;
    cout << "The coefs of " << 8 - n_cols << "-th = \n" << mz << endl;
    cout << "The coefs of " << 8 - n_cols << "-th = \n" << mb << endl;
    cout << "The coefs of " << 8 - n_cols << "-th = \n" << my << endl;
}
