
#include "../include/mavgrasp/setpoint_pub.h"

using namespace std;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "grasp_setpoint_node");
    SetPoint sp;
    ros::spin();
    return 0;
}
