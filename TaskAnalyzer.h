// TaskAnalyzer.h

#pragma once

#include "ContactCalculator.h"
#include "Shape.h"
#include "InfinitesimulDisplacement.h"
#include "DisplacementIndex.h"
#include "VisionErrorCorrector.h"

class TaskAnalyzer
{
public:
	typedef Eigen::Matrix<double, 3, 4> Pose;
public:
	struct ContactState
	{
		std::vector<Pose> poses; 
		std::vector<ContactElementBase> contact_elements;
		InfinitesimulDisplacement disp;
		DisplacementIndex index;
		ContactState(const ContactCalculator::ContactState &c_state) {
			contact_elements.resize(c_state.size());
			for (size_t i(0); i < c_state.size(); ++i) {
				contact_elements[i] = c_state[i];
			}
		}
		bool operator==(const ContactState &src) const {
			if (contact_elements.size() != src.contact_elements.size()) {
				return false;
			}
			for (size_t i(0); i < src.contact_elements.size(); ++i) {
				if (std::find(contact_elements.begin(), contact_elements.end(), src.contact_elements[i]) == contact_elements.end()) {
					return false;
				}
			}
			return true;
		}
	};
private:
	std::vector<ContactState> contact_states;
public:
	TaskAnalyzer(){}
	//
	size_t NumberOfStates() const {
		return contact_states.size();
	}
	void AppendContactState(const std::vector<Shape> &objects, const ContactCalculator::ContactState &c_state);
	
	void SetPoseOfFixedObject(const Eigen::Matrix3d &R, const Eigen::Vector3d &t);

	void Analyze(Shape &moving_object, Shape &fixed_object);
	//
	static void RecalculateContactState(ContactState &c_state, Shape &moving_object, Shape &fixed_object);
	// access result
	const ContactState& operator[](size_t src) const {
		return contact_states[src];
	}
	static void SetContactElements(VisionErrorCorrector::ContactState &vc_cs, const ContactState &c_state);
	friend std::ostream& operator<<(std::ostream &os, const TaskAnalyzer &src);
private:
	static void FurtherAnalysis(ContactState &c_state, Shape &moving_object, Shape &fixed_object);
	static Pose GetPose(const Shape &object) {
		Pose dest;
		dest.block(0, 0, 3, 3) = object.Rot();
		dest.block(0, 3, 3, 1) = object.Trans();
		return dest;
	}
	static void SetPose(Shape &object, const Pose &pose){
		Eigen::Matrix3d R = pose.block(0, 0, 3, 3);
		Eigen::Vector3d t = pose.block(0, 3, 3, 1);
		object.SetTransformation(R, t);
	}
};