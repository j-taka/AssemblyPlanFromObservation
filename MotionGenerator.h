// MotionGenerator.h

#pragma once

#include "TaskAnalyzer.h"

class MotionGenerator
{
private:
	double pos_range;
	double angle_range;
public:
	// constructor
	MotionGenerator() : pos_range(10.0), angle_range(M_PI / 6) {}
	// 
	void RandomPoseGeneration(TaskAnalyzer::Pose &pose, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &c_state) const;
private:
	Eigen::Matrix3d RandomPose() const;
	Eigen::Matrix3d RRisTwo(const TaskAnalyzer::ContactState &c_state) const;
	Eigen::Matrix3d RRisOne(const TaskAnalyzer::ContactState &c_state) const;
};