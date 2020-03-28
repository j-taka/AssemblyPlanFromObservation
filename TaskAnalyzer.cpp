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
	ContactCalculator cc;
	cc.SetThreshold(0.01);
	for (size_t i(0); i < contact_states.size(); ++i) {
		SetPose(moving_object, contact_states[i].poses[0]);
		SetPose(fixed_object, contact_states[i].poses[1]);
		cc.SetContact(contact_states[i].contact_elements);
		// 
		ContactStateForDisplacement cs;
		cc.Convert(cs, moving_object, fixed_object);
		// 
		contact_states[i].disp.Calculate(0, cs);
		contact_states[i].index.Calculate(contact_states[i].disp);
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

