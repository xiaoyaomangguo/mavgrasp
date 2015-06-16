cmake_minimum_required(VERSION 2.8.3)
project(mavgrasp)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  cmake_modules
)
find_package(Eigen REQUIRED)


catkin_package(
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)

add_executable(grasp_setpoint_node src/setpoint_pub.cpp
                             src/main.cpp)
add_dependencies(grasp_setpoint_node mavgrasp_generate_messages_cpp)
target_link_libraries(grasp_setpoint_node
  ${catkin_LIBRARIES}
)

