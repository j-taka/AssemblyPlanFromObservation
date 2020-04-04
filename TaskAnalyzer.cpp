// TaskAnalyzer.cpp

#include "TaskAnalyzer.h"
#include "ContactStateForDisplacement.h"

void TaskAnalyzer::AppendContactState(const std::vector<Shape> &objects, const ContactCalculator::ContactState &c_state)
{
	ContactState _c_state(c_state);
	std::vector<ContactState>::iterator it;
	if ((it = std::find(contact_states.begin(), contact_states.end(), c_state)) == contact_states.end()) {
		_c_state.poses.push_back(GetPose(objects[0]));
		_c_state.poses.push_back(GetPose(objects[1]));
		contact_states.push_back(_c_state);
	}
	else {
		it->poses.push_back(GetPose(objects[0]));
		it->poses.push_back(GetPose(objects[1]));
	}
}

void TaskAnalyzer::Analyze(Shape &moving_object, Shape &fixed_object)
{
	for (size_t i(0); i < contact_states.size(); ++i) {
		RecalculateContactState(contact_states[i], moving_object, fixed_object);
	}
}

// when changing the pose
void TaskAnalyzer::RecalculateContactState(ContactState &c_state, Shape &moving_object, Shape &fixed_object)
{
	SetPose(moving_object, c_state.poses[0]);
	SetPose(fixed_object, c_state.poses[1]);

	ContactCalculator cc;
	cc.SetThreshold(0.01);
#if 1 
	cc.SetContact(c_state.contact_elements);
#else
	cc.Calc(moving_object, fixed_object);
	ContactCalculator::PrintContact(cc.GetContact());
	std::cout << std::endl;
	cc.DetailedAnalysis(moving_object, fixed_object);
	ContactCalculator::PrintContact(cc.GetContact());
#endif
	// 
	ContactStateForDisplacement cs;
	cc.Convert(cs, moving_object, fixed_object);
	// 
	c_state.disp.Calculate(0, cs);
	c_state.index.Calculate(c_state.disp);
	if (c_state.index.RotationRestrictedDOF() == 1) {
		// special case
		switch (c_state.index.SetPossibleAxis(cs)) {
		case -1:
			std::cerr << "Cannot handle this case in the current implementation" << std::endl;
			exit(-1);
		case 1:
			FurtherAnalysis(c_state, moving_object, fixed_object);
			break;
		default:
			break;
		}
	}
}

void TaskAnalyzer::FurtherAnalysis(ContactState &c_state, Shape &moving_object, Shape &fixed_object)
{
	// move small
	const double small_angle = M_PI / 36.0; // 5 degrees
	const Eigen::Vector3d axis1 = c_state.index.PossibleAxis().row(0).transpose();
	Eigen::Matrix3d m1;
	m1 = Eigen::AngleAxisd(small_angle, axis1);
	const Eigen::Vector3d axis2 = c_state.index.PossibleAxis().row(1).transpose();
	Eigen::Matrix3d m2;
	m2 = Eigen::AngleAxisd(small_angle, axis2);
	// check two types
	Pose pose;
	pose.block(0, 3, 3, 1) = c_state.poses[0].block(0, 3, 3, 1);
	moving_object.SetTransformation(m1 * m2 * pose.block(0, 0, 3, 3), pose.block(0, 3, 3, 1));
	fixed_object.SetTransformation(c_state.poses[1].block(0, 0, 3, 3), c_state.poses[1].block(0, 3, 3, 1));
	VisionErrorCorrector::ContactState vc_cs;
	SetContactElements(vc_cs, c_state);
	VisionErrorCorrector vc;
	vc.Translation(moving_object, fixed_object, vc_cs);
	const double error1 = vc.CalculateMaximumError(moving_object, fixed_object, vc_cs);
	moving_object.SetTransformation(m2 * m1 * pose.block(0, 0, 3, 3), pose.block(0, 3, 3, 1));
	vc.Translation(moving_object, fixed_object, vc_cs);
	const double error2 = vc.CalculateMaximumError(moving_object, fixed_object, vc_cs);
	if (error1 > error2) {
		c_state.index.SwapPossibleAxis();
	}
}

void TaskAnalyzer::SetContactElements(VisionErrorCorrector::ContactState &vc_cs, const TaskAnalyzer::ContactState &c_state)
{
	vc_cs.resize(c_state.contact_elements.size());
	for (size_t i(0); i < c_state.contact_elements.size(); ++i) {
		vc_cs[i] = c_state.contact_elements[i];
	}
}

void TaskAnalyzer::SetPoseOfFixedObject(const Eigen::Matrix3d &R, const Eigen::Vector3d &t)
{
	Eigen::Matrix4d fixM = Eigen::Matrix4d::Identity();
	fixM.block(0, 0, 3, 3) = R;
	fixM.block(0, 3, 3, 1) = t;
	for (size_t i(0); i < contact_states.size(); ++i) {
		Eigen::Matrix4d fixM_org = Eigen::Matrix4d::Identity();
		fixM_org.block(0, 0, 3, 4) = contact_states[i].poses[1];
		Eigen::Matrix4d d = fixM * fixM_org.inverse();
		Eigen::Matrix4d movM_org = Eigen::Matrix4d::Identity();
		movM_org.block(0, 0, 3, 4) = contact_states[i].poses[0];
		Eigen::Matrix4d movM = d * movM_org;
		contact_states[i].poses[0] = movM.block(0, 0, 3, 4);
		contact_states[i].poses[1] = fixM.block(0, 0, 3, 4);
	}
}

std::ostream& operator<<(std::ostream &os, const TaskAnalyzer &src)
{
	os << "Number of states: " << src.contact_states.size() << std::endl;
	for (size_t i(0); i < src.contact_states.size(); ++i) {
		os << "State: " << i + 1 << std::endl;
		if (src.contact_states[i].contact_elements.empty()) {
			os << "No contact" << std::endl;
		}
		else {
			for (size_t j(0); j < src.contact_states[i].contact_elements.size(); ++j) {
				os << src.contact_states[i].contact_elements[j] << std::endl;
			}
		}
		os << (src.contact_states[i].disp.isSingular() ? "Singular" : "Non singular") << std::endl;
		os << src.contact_states[i].index << std::endl;
	}
	return os;
}

