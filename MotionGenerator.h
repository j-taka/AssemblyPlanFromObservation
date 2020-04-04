// MotionGenerator.h

#pragma once

#include "TaskAnalyzer.h"
#include "VisionErrorCorrector.h"

class MotionGenerator
{
private:
	double pos_range;
	double angle_range;
	double INTERPOLATIONANGLE;
	double INTERPOLATIONLENGTH;
public:
	// constructor
	MotionGenerator() : pos_range(10.0), angle_range(M_PI / 6), INTERPOLATIONANGLE(M_PI / 180.0), INTERPOLATIONLENGTH(1.0) {}
	// 
	void RandomPoseGeneration(TaskAnalyzer::Pose &pose, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &c_state) const;
	//
	void CalculateTrajectory(std::vector<TaskAnalyzer::Pose> &trajectory, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;

	void SetRangeInTranslation(double _pos_range) {
		pos_range = _pos_range;
	}
private:
	Eigen::Matrix3d RandomPose() const;
	Eigen::Matrix3d RRisTwo(const TaskAnalyzer::ContactState &c_state) const;
	Eigen::Matrix3d RRisOne(const TaskAnalyzer::ContactState &c_state) const;

	void SetTranslationParameter(Eigen::Vector3d &direction, double &distance, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::Pose &pose, const TaskAnalyzer::ContactState &end) const;
	void SetRotationParameter0To1(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;
	void SetRotationParameter0To2(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;
	void SetRotationParameterSR1(std::vector<Eigen::Vector3d> &axis, std::vector<double> &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;
	void SetRotationParameterER3(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;

	void InterpolationForTranslation(std::vector<TaskAnalyzer::Pose> &trajectory, const Eigen::Vector3d &direction, const double distance, const TaskAnalyzer::Pose &pose) const;
	void InterpolationForRotation(std::vector<TaskAnalyzer::Pose> &trajectory, const Eigen::Vector3d &axis, const double angle, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &c_state) const;

	void MoveALittle(std::vector<TaskAnalyzer::Pose> &trajectory, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &org_start, const TaskAnalyzer::ContactState &end) const;

	void AxisToMatrix0To1(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;
	void AxisToMatrix0To2(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;
	void AxisToMatrixSR1(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const;

	static void MatrixToRZYX(Eigen::Vector3d &zyx, const Eigen::Matrix3d &mat);
	static void MatrixToRZNyz(Eigen::Vector2d &zn, const Eigen::Vector2d &nyz, const Eigen::Matrix3d &mat);
};