#include "MotionGenerator.h"
#include "s2rand.h"

static const double TRIFUNC_SMALL = 1.0e-3;

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
	TaskAnalyzer::SetContactElements(vc_cs, c_state);
	VisionErrorCorrector vc;
#if 1
	vc.Translation(moving_object, fixed_object, vc_cs);
	pose.block(0, 3, 3, 1) = moving_object.Trans();
#else
	vc.Calc(moving_object, fixed_object, vc_cs);
	pose.block(0, 0, 3, 3) = moving_object.Rot();
	pose.block(0, 3, 3, 1) = moving_object.Trans();
#endif
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

		const int r = static_cast<int>(tr_disp.rows());
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
	double angle = rand() / (double)RAND_MAX * angle_range * 1.5;
	Eigen::Matrix3d res;
	res = Eigen::AngleAxisd(angle, axis);
	return res;
}

void MotionGenerator::CalculateTrajectory(std::vector<TaskAnalyzer::Pose> &trajectory, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &org_start, const TaskAnalyzer::ContactState &end) const
{
	// calculate current admissible displacement
	TaskAnalyzer::ContactState cur_start(org_start);
	cur_start.poses[0] = trajectory.back();
	TaskAnalyzer::RecalculateContactState(cur_start, moving_object, fixed_object);
	TaskAnalyzer::Pose pose;
	if (cur_start.index.RotationRestrictedDOF() < end.index.RotationRestrictedDOF()) {
		if (cur_start.index.RotationRestrictedDOF() == 1) { // special case
			// 1 -> 2, 1 -> 3
			std::vector<Eigen::Vector3d> axis;
			std::vector<double> angle;
			SetRotationParameterSR1(axis, angle, cur_start, end);
			for (size_t i(0); i < angle.size(); ++i) {
				// std::cout << axis[i].transpose() << " " << angle[i] << std::endl;
				InterpolationForRotation(trajectory, axis[i], angle[i], moving_object, fixed_object, cur_start);
			}
		}
		else if (end.index.RotationRestrictedDOF() == 1) {
			// 0 -> 1
			Eigen::Vector3d axis;
			double angle; 
			SetRotationParameter0To1(axis, angle, cur_start, end);
			InterpolationForRotation(trajectory, axis, angle, moving_object, fixed_object, cur_start);
		}
		else if (end.index.RotationRestrictedDOF() == 3) { 
			// 0 -> 3, 2 -> 3
			Eigen::Vector3d axis;
			double angle;
			SetRotationParameterER3(axis, angle, cur_start, end);
			InterpolationForRotation(trajectory, axis, angle, moving_object, fixed_object, cur_start);
		}
		else {
			// 0 -> 2
			Eigen::Vector3d axis;
			double angle;
			SetRotationParameter0To2(axis, angle, cur_start, end);
			InterpolationForRotation(trajectory, axis, angle, moving_object, fixed_object, cur_start);
		}
	}
	pose = trajectory.back();
	if (cur_start.index.TranslactionRestrictedDOF() < end.index.TranslactionRestrictedDOF()) {
		Eigen::Vector3d direction;
		double distance;
		SetTranslationParameter(direction, distance, moving_object, fixed_object, pose, end);
		InterpolationForTranslation(trajectory, direction, distance, pose);
	}
	if (cur_start.disp.isSingular() && cur_start.index.TranslactionRestrictedDOF() >= end.index.TranslactionRestrictedDOF() &&
		cur_start.index.RotationRestrictedDOF() >= end.index.RotationRestrictedDOF()) {
		MoveALittle(trajectory, moving_object, fixed_object, org_start, end);
	}
}

void MotionGenerator::SetTranslationParameter(Eigen::Vector3d &direction, double &distance, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::Pose &pose, const TaskAnalyzer::ContactState &end) const
{
	moving_object.SetTransformation(pose.block(0, 0, 3, 3), pose.block(0, 3, 3, 1));
	fixed_object.SetTransformation(end.poses[1].block(0, 0, 3, 3), end.poses[1].block(0, 3, 3, 1));
	VisionErrorCorrector::ContactState vc_cs;
	TaskAnalyzer::SetContactElements(vc_cs, end);
	VisionErrorCorrector vc;
	vc.Translation(moving_object, fixed_object, vc_cs);
	direction = moving_object.Trans() - pose.block(0, 3, 3, 1);
	if (direction.norm() < 1.0e-4) {
		distance = 0;
		direction = Eigen::Vector3d::UnitZ(); // no meaning
	}
	else {
		distance = direction.norm();
		direction /= distance;
	}
}

void MotionGenerator::SetRotationParameterER3(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	const Eigen::Matrix3d S = start.poses[0].block(0, 0, 3, 3);
	const Eigen::Matrix3d E = end.poses[0].block(0, 0, 3, 3);
	const Eigen::Matrix3d ESI = E * S.transpose();
	const Eigen::AngleAxisd res(ESI);
	if (res.angle() < 0) {
		axis = -res.axis();
		angle = -res.angle();
	}
	else {
		axis = res.axis();
		angle = res.angle();
	}
}

// 1 -> 2, 1 -> 3
void MotionGenerator::SetRotationParameterSR1(std::vector<Eigen::Vector3d> &axis, std::vector<double> &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	const Eigen::Matrix3d SM = start.poses[0].block(0, 0, 3, 3);
	const Eigen::Matrix3d EM = end.poses[0].block(0, 0, 3, 3);
	Eigen::Matrix3d AM;
	AxisToMatrixSR1(AM, axis, start, end);
#if 0
	std::cout << "axis" << std::endl;
	std::cout << axis[0].transpose() << std::endl;
	std::cout << axis[1].transpose() << std::endl;
#endif
	// set rotation axes to x, y, z-axes
	const Eigen::Matrix3d SSM = AM * SM;
	const Eigen::Matrix3d EEM = AM * EM;
	const Eigen::Matrix3d SEM = EEM * SSM.transpose(); // inverse
	// how much rotates around each axis
	Eigen::Vector2d zn;
	Eigen::Vector3d t_axis = AM * axis[0];
	MatrixToRZNyz(zn, t_axis.block(1, 0, 2, 1), SEM);
#if 0
	// debug
	Eigen::Matrix3d R;
	R = Eigen::AngleAxisd(zn[0], Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(zn[1], t_axis);
	std::cout << SEM << std::endl;
	std::cout << R << std::endl << std::endl;
#endif
	if (end.index.RotationRestrictedDOF() == 2) {
		// std::cout << "target object: " << start.index.TargetObjectID() << std::endl;
		if (start.index.TargetObjectID() == 0) {
			axis.erase(axis.begin());
			angle.resize(1);
			angle[0] = zn[0];
		}
		else {
			axis.erase(axis.begin() + 1);
			angle.resize(1);
			angle[0] = zn[1];
		}
	}
	else {
		angle.resize(2);
		angle[0] = zn[1];
		angle[1] = zn[0];
	}
	for (size_t i(0); i < angle.size(); ++i) {
		if (angle[i] < 0) {
			angle[i] = -angle[i];
			axis[i] = -axis[i];
		}
	}
}

void MotionGenerator::SetRotationParameter0To1(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	const Eigen::Matrix3d SM = start.poses[0].block(0, 0, 3, 3);
	const Eigen::Matrix3d EM = end.poses[0].block(0, 0, 3, 3);
	std::vector<Eigen::Vector3d> axes;
	Eigen::Matrix3d AM;
	AxisToMatrix0To1(AM, axes, start, end);
	// set rotation axes to x, y, z-axes
	const Eigen::Matrix3d SSM = AM * SM;
	const Eigen::Matrix3d EEM = AM * EM;
	const Eigen::Matrix3d SEM = EEM * SSM.transpose(); // inverse
	// how much rotates around each axis
	Eigen::Vector3d zyx;
	MatrixToRZYX(zyx, SEM);
	Eigen::Matrix3d RM = Eigen::Matrix3d::Identity();
	for (size_t i(0); i < axes.size(); ++i) {
		Eigen::Matrix3d tM;
		tM = Eigen::AngleAxisd(zyx[2 - i], axes[i]);
		RM = tM * RM;
	}
	Eigen::AngleAxisd dest(RM);
	if (dest.angle() < 0) {
		angle = -dest.angle();
		axis = -dest.axis();
	}
	else {
		angle = dest.angle();
		axis = dest.axis();
	}
}

void MotionGenerator::SetRotationParameter0To2(Eigen::Vector3d &axis, double &angle, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	const Eigen::Matrix3d SM = start.poses[0].block(0, 0, 3, 3);
	const Eigen::Matrix3d EM = end.poses[0].block(0, 0, 3, 3);
	std::vector<Eigen::Vector3d> axes;
	Eigen::Matrix3d AM;
	AxisToMatrix0To2(AM, axes, start, end);
	// set rotation axes to x, y, z-axes
	const Eigen::Matrix3d SSM = AM * SM;
	const Eigen::Matrix3d EEM = AM * EM;
	const Eigen::Matrix3d SEM = EEM * SSM.transpose(); // inverse
	// how much rotates around each axis
	Eigen::Vector3d zyx;
	MatrixToRZYX(zyx, SEM);
	Eigen::Matrix3d RM = Eigen::Matrix3d::Identity();
	for (size_t i(0); i < axes.size(); ++i) {
		Eigen::Matrix3d tM;
		tM = Eigen::AngleAxisd(zyx[2 - i], axes[i]);
		RM = tM * RM;
	}
	Eigen::AngleAxisd dest(RM);
	if (dest.angle() < 0) {
		angle = -dest.angle();
		axis = -dest.axis();
	}
	else {
		angle = dest.angle();
		axis = dest.axis();
	}
}

void MotionGenerator::InterpolationForTranslation(std::vector<TaskAnalyzer::Pose> &trajectory, const Eigen::Vector3d &direction, const double distance, const TaskAnalyzer::Pose &pose) const
{
	assert(distance >= 0);
	int iNum = static_cast<int>(distance / INTERPOLATIONLENGTH) + 1;
	for (int i(0); i < iNum; ++i) {
		double tmp_distance = distance * (i + 1) / iNum;
		TaskAnalyzer::Pose tmp_pose = pose;
		tmp_pose.block(0, 3, 3, 1) += tmp_distance * direction;
		trajectory.push_back(tmp_pose);
	}
}

void MotionGenerator::InterpolationForRotation(std::vector<TaskAnalyzer::Pose> &trajectory, const Eigen::Vector3d &axis, const double angle, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &c_state) const
{
	assert(angle >= 0);
	Eigen::Matrix3d initR = moving_object.Rot();
	Eigen::Vector3d initT = moving_object.Trans();
	int iNum = static_cast<int>(angle / INTERPOLATIONANGLE) + 1;
	for (int i(0); i < iNum; ++i) {
		double tmp_angle = angle * (i + 1) / iNum;
		Eigen::Matrix3d dR;
		dR = Eigen::AngleAxisd(tmp_angle, axis);
		moving_object.SetTransformation(dR * initR, initT);
		VisionErrorCorrector::ContactState vc_cs;
		TaskAnalyzer::SetContactElements(vc_cs, c_state);
		VisionErrorCorrector vc;
		vc.Translation(moving_object, fixed_object, vc_cs);
		TaskAnalyzer::Pose pose;
		pose.block(0, 0, 3, 3) = moving_object.Rot();
		pose.block(0, 3, 3, 1) = moving_object.Trans();
		trajectory.push_back(pose);
	}
}

// follow the demonstration
void MotionGenerator::MoveALittle(std::vector<TaskAnalyzer::Pose> &trajectory, Shape &moving_object, Shape &fixed_object, const TaskAnalyzer::ContactState &org_start, const TaskAnalyzer::ContactState &end) const
{
	Eigen::Vector3d dt = end.poses[0].block(0, 3, 3, 1) - org_start.poses[0].block(0, 3, 3, 1);
	dt *= INTERPOLATIONLENGTH / dt.norm();
	// std::cout << dt.transpose() << std::endl;
	moving_object.SetTransformation(trajectory.back().block(0, 0, 3, 3), trajectory.back().block(0, 3, 3, 1) + dt);
	VisionErrorCorrector::ContactState vc_cs;
	TaskAnalyzer::SetContactElements(vc_cs, end);
	VisionErrorCorrector vc;
	vc.Translation(moving_object, fixed_object, vc_cs);
	TaskAnalyzer::Pose pose;
	pose.block(0, 0, 3, 3) = moving_object.Rot();
	pose.block(0, 3, 3, 1) = moving_object.Trans();
	trajectory.push_back(pose);
}


// sN == 0, eN == 1
void MotionGenerator::AxisToMatrix0To1(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	assert(start.index.RotationRestrictedDOF() == 0);
	assert(end.index.RotationRestrictedDOF() == 1);
	axes.resize(1);
	const Eigen::Vector3d e1 = end.index.PossibleAxis().row(1).transpose();
	const Eigen::Vector3d e2 = end.index.PossibleAxis().row(0).transpose();
	// std::cout << e1.dot(e2) << std::endl;
	// std::cout << end.index.PossibleAxis() << std::endl;
	if (fabs(e1.dot(e2)) > TRIFUNC_SMALL) {
		std::cerr << "Cannot handle this case in the current implementation." << std::endl;
			exit(-1);
	}
	axes[0] = e2.cross(e1);
	axes[0].normalize();
	AM.row(0) = axes[0];
	AM.row(1) = e2.transpose();
	AM.row(2) = e1.transpose();
}

// sN == 0, eN == 2
void MotionGenerator::AxisToMatrix0To2(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	assert(start.index.RotationRestrictedDOF() == 0);
	assert(end.index.RotationRestrictedDOF() == 2);
	axes.resize(2);
	const Eigen::Vector3d e1 = end.index.PossibleAxis().row(0).transpose();
	double min_dot = fabs(e1.dot(start.index.PossibleAxis().row(0).transpose()));
	int min_ID(0);
	for (int i(1); i < 3; ++i) {
		if (min_dot > fabs(e1.dot(start.index.PossibleAxis().row(i).transpose()))) {
			min_dot = fabs(e1.dot(start.index.PossibleAxis().row(i).transpose()));
			min_ID = 1;
		}		
	}
	axes[0] = start.index.PossibleAxis().row(min_ID).transpose();
	axes[0] -= e1.dot(axes[0]) * e1;
	axes[0].normalize();
	axes[1] = axes[0].cross(e1);
	AM.row(0) = axes[0].transpose();
	AM.row(1) = axes[1].transpose();
	AM.row(2) = e1.transpose();
}

void MotionGenerator::AxisToMatrixSR1(Eigen::Matrix3d &AM, std::vector<Eigen::Vector3d> &axes, const TaskAnalyzer::ContactState &start, const TaskAnalyzer::ContactState &end) const
{
	assert(start.index.RotationRestrictedDOF() == 1);
	axes.resize(2);
	axes[0] = start.index.PossibleAxis().row(1).transpose();
	axes[1] = start.index.PossibleAxis().row(0).transpose(); // z
	Eigen::Vector3d y = axes[0] - axes[0].dot(axes[1]) * axes[1];
	y.normalize();
	Eigen::Vector3d x = y.cross(axes[1]);
	AM.row(0) = x;
	AM.row(1) = y;
	AM.row(2) = axes[1].transpose();
}

void MotionGenerator::MatrixToRZYX(Eigen::Vector3d &zyx, const Eigen::Matrix3d &mat)
{
	double n1 = mat(0, 0);
	double n2 = mat(1, 0);
	double n3 = mat(2, 0);
	double o1 = mat(0, 1);
	double o2 = mat(1, 1);
	double a1 = mat(0, 2);
	double a2 = mat(1, 2);

	zyx[0] = atan2(n2, n1);
	zyx[1] = atan2(-n3, (cos(zyx[0]) * n1) + (sin(zyx[0]) * n2));
	zyx[2] = atan2((sin(zyx[0]) * a1) - (cos(zyx[0]) * a2),
		(cos(zyx[0]) * o2) - (sin(zyx[0]) * o1));
}

void MotionGenerator::MatrixToRZNyz(Eigen::Vector2d &zn, const Eigen::Vector2d &nyz, const Eigen::Matrix3d &mat)
{
	const double ny2 = nyz[0] * nyz[0];
	zn[1] = atan2(-mat(2, 0) / nyz[0], (mat(2, 2) + ny2 - 1) / ny2);
	const double ca = cos(zn[1]);
	const double sa = sin(zn[1]);
	if (fabs(sa) < TRIFUNC_SMALL){
		zn[0] = atan2(mat(1, 0), mat(0, 0));
	}
	else {
		zn[0] = atan2(-mat(1, 0) - mat(0, 1), mat(1, 1) - mat(0, 0));
	}
}