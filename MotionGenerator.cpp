#include "MotionGenerator.h"
#include "s2rand.h"
#include "VisionErrorCorrector.h"

// only consider keeping the current contact state. there is the possibily to collide each other.
// asssume that edges and faces are infnitely bigger. 
void MotionGenerator::RandomPoseGeneration(TaskAnalyzer::Pose &pose, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &c_state) const
{
	// first orientation
	switch (c_state.index.RotationRestrictedDOF()) {
	case 0: // completele free
		pose.block(0, 0, 3, 3) = RandomPose() * c_state.poses[0].block(0, 0, 3, 3);
		break;
	case 1:
		pose.block(0, 0, 3, 3) = RRisOne(c_state) * c_state.poses[0].block(0, 0, 3, 3);
		break;
	case 2:
		pose.block(0, 0, 3, 3) = RRisTwo(c_state) * c_state.poses[0].block(0, 0, 3, 3);
		break;
	case 3: // just copy
		pose.block(0, 0, 3, 3) = c_state.poses[0].block(0, 0, 3, 3);
		break;
	}
	pose.block(0, 3, 3, 1) = c_state.poses[0].block(0, 3, 3, 1);
	moving_object.SetTransformation(pose.block(0, 0, 3, 3), pose.block(0, 3, 3, 1));
	fixed_object.SetTransformation(c_state.poses[1].block(0, 0, 3, 3), c_state.poses[1].block(0, 3, 3, 1));
	VisionErrorCorrector::ContactState vc_cs;
	vc_cs.resize(c_state.contact_elements.size());
	for (size_t i(0); i < c_state.contact_elements.size(); ++i) {
		vc_cs[i] = c_state.contact_elements[i];
	}
	VisionErrorCorrector vc;
	vc.Translation(moving_object, fixed_object, vc_cs);
	pose.block(0, 3, 3, 1) = moving_object.Trans();
	// next translation
	if (c_state.index.TranslactionRestrictedDOF() == 0) {
		pose.block(0, 3, 3, 1) += Eigen::Vector3d::Random() * pos_range;
	}
	else if (c_state.index.TranslactionRestrictedDOF() != 3){
		// since changing the orientation, it is necessary to calculate the possible translation directions
		ContactCalculator cc;
		cc.SetThreshold(0.01);
		cc.SetContact(c_state.contact_elements);
		// 
		ContactStateForDisplacement cs;
		cc.Convert(cs, moving_object, fixed_object);
		// 
		InfinitesimulDisplacement tmp_disp;
		tmp_disp.Calculate(0, cs);
		DisplacementIndex tmp_index;
		Eigen::MatrixXd tr_disp;
		tmp_index.CalculatePossibleTranslationDirections(tr_disp, tmp_disp);

		const int r = tr_disp.rows();
		pose.block(0, 3, 3, 1) += pos_range * (tr_disp.transpose() * Eigen::VectorXd::Random(r));
	}
	return;
}

Eigen::Matrix3d MotionGenerator::RRisTwo(const TaskAnalyzer::ContactState &c_state) const
{
	const Eigen::Vector3d axis = c_state.index.PossibleAxis().row(0).transpose();
	const double angle = (rand() / (double)RAND_MAX * 2.0 - 1.0) * angle_range;
	Eigen::Matrix3d res;
	res = Eigen::AngleAxisd(angle, axis);
	return res;
}

Eigen::Matrix3d MotionGenerator::RRisOne(const TaskAnalyzer::ContactState &c_state) const
{
	const Eigen::Vector3d axis1 = c_state.index.PossibleAxis().row(0).transpose();
	const double angle1 = (rand() / (double)RAND_MAX * 2.0 - 1.0) * angle_range;
	Eigen::Matrix3d m1;
	m1 = Eigen::AngleAxisd(angle1, axis1);
	const Eigen::Vector3d axis2 = c_state.index.PossibleAxis().row(1).transpose();
	const double angle2 = (rand() / (double)RAND_MAX * 2.0 - 1.0) * angle_range;
	Eigen::Matrix3d m2;
	m2 = Eigen::AngleAxisd(angle2, axis2);
	return m1 * m2;
}

Eigen::Matrix3d MotionGenerator::RandomPose() const
{
	double _axis[3];
	s2rand(_axis);
	const Eigen::Vector3d axis(_axis[0], _axis[1], _axis[2]);
	double angle = rand() / (double)RAND_MAX * angle_range;
	Eigen::Matrix3d res;
	res = Eigen::AngleAxisd(angle, axis);
	return res;
}