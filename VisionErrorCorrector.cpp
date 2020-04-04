// VisionErrorCorrector.cpp

#include "VisionErrorCorrector.h"
#include "ContactCalculator.h"
#include <Eigen/Dense>

const double _NEARLY_ZERO = 1.0e-6;

void VisionErrorCorrector::FixFixedObject(PoseListFileHandler &plHandler)
{
	Eigen::Matrix3d fixR;
	Eigen::Vector3d fixT;
	plHandler.GetTransformation(fixR, fixT, 0, 1);
	for (size_t i(1); i < plHandler.length(); ++i) {
		Eigen::Matrix3d R;
		Eigen::Vector3d t;
		plHandler.GetTransformation(R, t, i, 1);
		const Eigen::Matrix3d dR = fixR * R.transpose();
		const Eigen::Vector3d dt = -fixR * R.transpose() * t + fixT;
		Eigen::Matrix3d R2;
		Eigen::Vector3d t2;
		plHandler.GetTransformation(R2, t2, i, 0);
		const Eigen::Matrix3d newR = dR * R2;
		const Eigen::Vector3d newT = dR * t2 + dt;
		plHandler.SetTransformation(newR, newT, i, 0);
		plHandler.SetTransformation(fixR, fixT, i, 1);
	}
}

bool VisionErrorCorrector::Calc(Shape &moving_object, Shape &fixed_object, ContactState &c_state)
{
	if (c_state.empty()) {
		return true; // nothing to do 
	}
	const Eigen::Matrix3d initR = moving_object.Rot();
	const Eigen::Vector3d initT = moving_object.Trans();
	for (int i(0); i < _MAX_NUMBER_OF_REMOVE_ELEMENT; ++i) {
		if (CalcEach(moving_object, fixed_object, c_state)) {
			return true;
		}
		// remove furthest contact element
		double max_dist = fabs(c_state[0].GetDistance());
		size_t ID = 0;
		for (size_t j(1); j < c_state.size(); ++j) {
			if (max_dist < fabs(c_state[j].GetDistance())) {
				max_dist = fabs(c_state[j].GetDistance());
				ID = j;
			}
		}
		std::cerr << "Remove " << c_state[ID] << std::endl;
		c_state.erase(c_state.begin() + ID);
		moving_object.SetTransformation(initR, initT);
	}
	return false;
}

void VisionErrorCorrector::Translation(Shape &moving_object, Shape &fixed_object, const ContactState &c_state)
{
	if (c_state.empty()) {
		return; // nothing to do 
	}
	// center 
	CalcCenter(moving_object, fixed_object, c_state);
	// error correction only using translation
	for (size_t i(0); i < c_state.size(); ++i) {
		SetEquation(moving_object, fixed_object, c_state[i]);
	}
	Eigen::Vector3d dt;
	SolveEquation(dt);
	moving_object.SetTransformation(moving_object.Rot(), moving_object.Trans() + dt);
	// std::cout << CalculateMaximumError(moving_object, fixed_object, c_state) << std::endl;
}


bool VisionErrorCorrector::CalcEach(Shape &moving_object, Shape &fixed_object, const ContactState &c_state)
{
#ifdef _VC_DEBUG_MODE
	convergence.clear();
	Eigen::Matrix<double, 3, 4> pose;
	pose.block(0, 0, 3, 3) = moving_object.Rot();
	pose.block(0, 3, 3, 1) = moving_object.Trans();
	convergence.push_back(pose);
#endif
	// 
	ReduceSecondOrderEffect(moving_object, fixed_object, c_state);
#ifdef _VC_DEBUG_MODE
	pose.block(0, 0, 3, 3) = moving_object.Rot();
	pose.block(0, 3, 3, 1) = moving_object.Trans();
	convergence.push_back(pose);
#endif
	for (int i(0); i < _MAX_LOOP; ++i) {
		// center 
		CalcCenter(moving_object, fixed_object, c_state);
		double error(0);
		for (size_t i(0); i < c_state.size(); ++i) {
			error = std::max(error, fabs(SetEquation(moving_object, fixed_object, c_state[i])));
		}
		if (verbose) {
			std::cerr << "Loop " << i + 1 << " : " << error << std::endl;
		}
		if (error < tolerance) {
			eq.clear();
			return true;
		}
		Eigen::Matrix3d dR;
		Eigen::Vector3d dt;
		// Solve Equation
#if 0
		if (i == 0) {
			std::ofstream ofs("debug_matrix.txt");
			for (size_t i(0); i < eq.size(); ++i) {
				ofs << eq[i].elm[0] << " " << eq[i].elm[1] << " " << eq[i].elm[2] << " " << eq[i].elm[3] << " "
					<< eq[i].elm[4] << " " << eq[i].elm[5] << " : " << eq[i].elm[6] << std::endl;
			}
			ofs.close();
		}
#endif
		SolveEquation(dR, dt);
		// increment
		Eigen::Matrix3d R2 = dR * moving_object.Rot();
		Eigen::Vector3d t2 = dR * (moving_object.Trans() - center) + center + dt;
		// avoid overly translating
		if ((t2 - moving_object.Trans()).norm() > _MIN_TRANS_DIST) {
			t2 = moving_object.Trans();
		}
		moving_object.SetTransformation(R2, t2);
#ifdef _VC_DEBUG_MODE
		pose.block(0, 0, 3, 3) = R2;
		pose.block(0, 3, 3, 1) = t2;
		convergence.push_back(pose);
#endif
		// error correction only using translation
		for (size_t i(0); i < c_state.size(); ++i) {
			SetEquation(moving_object, fixed_object, c_state[i]);
		}
		SolveEquation(dt);
		if (dt.norm() <= _MIN_TRANS_DIST) {
			moving_object.SetTransformation(moving_object.Rot(), moving_object.Trans() + dt);
#ifdef _VC_DEBUG_MODE
			pose.block(0, 3, 3, 1) = moving_object.Trans();
			convergence.push_back(pose);
#endif
		}
	}
	return false;
}

void VisionErrorCorrector::CalcCenter(const Shape &moving_object, const Shape &fixed_object, const ContactState &c_state)
{
	center = Eigen::Vector3d::Zero();
	for (size_t i(0); i < c_state.size(); ++i) {
		switch (c_state[i].ContactType()) {
		case ContactElementBase::_VF_CONTACT:
			if (c_state[i].FirstElement().first == 0) {
				center += moving_object.V(c_state[i].FirstElement().second);
			}
			else {
				center += fixed_object.V(c_state[i].FirstElement().second);
			}
			break;
		case ContactElementBase::_EE_CONTACT:
			center += 0.25 * (moving_object.VonE(c_state[i].FirstElement().second, 0) + moving_object.VonE(c_state[i].FirstElement().second, 1)
				+ fixed_object.VonE(c_state[i].SecondElement().second, 0) + fixed_object.VonE(c_state[i].SecondElement().second, 1));
			break;
		default:
			std::cerr << "Including bugs..." << std::endl;
			exit(-1);
		}
	}
	center /= static_cast<double>(c_state.size());
}

double VisionErrorCorrector::CalculateMaximumError(const Shape &moving_object, const Shape &fixed_object, const ContactState &c_state)
{
	double error(0);
	for (size_t i(0); i < c_state.size(); ++i) {
		error = std::max(error, CalcError(moving_object, fixed_object, c_state[i]));
	}
	return error;
}

double VisionErrorCorrector::SetEquation(const Shape &moving_object, const Shape &fixed_object, const ContactElementForErrorCorrection &c_element)
{
	switch (c_element.ContactType()) {
	case ContactElementBase::_VF_CONTACT:
		if (c_element.FirstElement().first == 0) {
			return SetEquationVF(moving_object, fixed_object, c_element.FirstElement().second, c_element.SecondElement().second, true);
		}
		else {
			return SetEquationVF(fixed_object, moving_object, c_element.FirstElement().second, c_element.SecondElement().second, false);
		}
	case ContactElementBase::_EE_CONTACT:
		assert(c_element.FirstElement().first == 0);
		return SetEquationEE(moving_object, fixed_object, c_element.FirstElement().second, c_element.SecondElement().second);
	default:
		return 0;
	}
}

// v-f
double VisionErrorCorrector::SetEquationVF(const Shape &v_shape, const Shape &f_shape, const size_t vID, const size_t fID, bool v_move)
{
	Each tmp;
	const Eigen::Vector3d vertex_pos = v_shape.V(vID) - center;
	const Eigen::Vector3d norm = f_shape.OuterNormal(fID);
	const Eigen::Vector3d face_pos = f_shape.VonF(fID, 0) - center;
	const Eigen::Vector3d temp_3d_vec = vertex_pos.cross(norm);
	if (v_move) {
		for (int j(0); j < 3; ++j) { tmp.elm[j] = -temp_3d_vec[j]; }
		for (int j(0); j < 3; ++j) { tmp.elm[j + 3] = -norm[j]; }
	}
	else {
		for (int j(0); j < 3; ++j) { tmp.elm[j] = temp_3d_vec[j]; }
		for (int j(0); j < 3; ++j) { tmp.elm[j + 3] = norm[j]; }
	}
	tmp.elm[6] = norm.dot(vertex_pos - face_pos);
	eq.push_back(tmp);
	return fabs(tmp.elm[6]);
}

double VisionErrorCorrector::SetEquationEE(const Shape &e1_shape, const Shape &e2_shape, const size_t e1ID, const size_t e2ID, const double _PARALLEL)
{
	const Eigen::Vector3d edge1_pos = e1_shape.VonE(e1ID, 0) - center;
	const Eigen::Vector3d edge2_pos = e2_shape.VonE(e2ID, 0) - center;
	const Eigen::Vector3d edge1_dic = e1_shape.VonE(e1ID, 1) - e1_shape.VonE(e1ID, 0);
	const Eigen::Vector3d edge2_dic = e2_shape.VonE(e2ID, 1) - e2_shape.VonE(e2ID, 0);
	Eigen::Vector3d dP = edge2_pos - edge1_pos;
	Eigen::Vector3d temp_3d_vec = edge1_dic.cross(edge2_dic);
	// check outside
	if (temp_3d_vec.dot(e1_shape.OutsideDirectionOfE(e1ID)) > 0) {
		temp_3d_vec = -temp_3d_vec;
	}
	const double length = temp_3d_vec.norm();
	// parallel
	if ((length / edge1_dic.norm() / edge2_dic.norm()) < _PARALLEL) {
		return 0;
	}
	dP /= length;
	const Eigen::Vector3d norm = temp_3d_vec / length;
	Eigen::Vector3d P = dP + norm * norm.dot(dP);
	P = P * (1 / length);
	temp_3d_vec = edge1_pos.cross(norm) + (edge2_dic.cross(dP)).cross(edge1_dic);
	Each tmp;
	for (int j(0); j < 3; j++) { tmp.elm[j] = temp_3d_vec[j]; }
	for (int j(0); j < 3; j++) { tmp.elm[j + 3] = norm[j]; }
	tmp.elm[6] = norm.dot(edge2_pos - edge1_pos);
	eq.push_back(tmp);
	return fabs(tmp.elm[6]);
}

void VisionErrorCorrector::SolveEquation(Eigen::Matrix3d &dR, Eigen::Vector3d &dt) 
{
	Eigen::MatrixXd A(eq.size(), 6);
	Eigen::VectorXd b(eq.size());
	for (size_t i(0); i < eq.size(); ++i) {
		for (int j(0); j < 6; ++j) {
			A(i, j) = eq[i].elm[j];
		}
		b[i] = eq[i].elm[6];
	}
	eq.clear();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd x = PseudoInverse(svd) * b;
	Eigen::Vector3d axis(x[0], x[1], x[2]);
	double angle = axis.norm();
	if (angle < _NEARLY_ZERO) {
		axis = Eigen::Vector3d::UnitZ();
angle = 0;
	}
	else {
	axis /= angle;
	}
	dR = Eigen::AngleAxisd(angle, axis);
	dt = Eigen::Vector3d(x[3], x[4], x[5]);
}

void VisionErrorCorrector::SolveEquation(Eigen::Vector3d &dt)
{
	// only use translation part
	Eigen::MatrixXd A(eq.size(), 3);
	Eigen::VectorXd b(eq.size());
	for (size_t i(0); i < eq.size(); ++i) {
		for (int j(0); j < 3; ++j) {
			A(i, j) = eq[i].elm[j + 3];
		}
		b[i] = eq[i].elm[6];
	}
	eq.clear();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	// 
	Eigen::VectorXd x = PseudoInverse(svd) *b;
	dt = Eigen::Vector3d(x[0], x[1], x[2]);
	// std::cout << (A * x - b).transpose() << std::endl;
}

double VisionErrorCorrector::CalcError(const Shape &moving_object, const Shape &fixed_object, const ContactElementForErrorCorrection &c_element) const
{
	switch (c_element.ContactType()) {
	case ContactElementBase::_VF_CONTACT:
		if (c_element.FirstElement().first == 0) {
			return CalcErrorVF(moving_object, fixed_object, c_element.FirstElement().second, c_element.SecondElement().second);
		}
		else {
			return CalcErrorVF(fixed_object, moving_object, c_element.FirstElement().second, c_element.SecondElement().second);
		}
	case ContactElementBase::_EE_CONTACT:
		assert(c_element.FirstElement().first == 0);
		return CalcErrorEE(moving_object, fixed_object, c_element.FirstElement().second, c_element.SecondElement().second);
	default:
		return 0;
	}
}

// v-f
double VisionErrorCorrector::CalcErrorVF(const Shape &v_shape, const Shape &f_shape, const size_t vID, const size_t fID) const
{
	const Eigen::Vector3d vertex_pos = v_shape.V(vID) - center;
	const Eigen::Vector3d norm = f_shape.OuterNormal(fID);
	const Eigen::Vector3d face_pos = f_shape.VonF(fID, 0) - center;
	return fabs(norm.dot(vertex_pos - face_pos));
}

// e-e 
double VisionErrorCorrector::CalcErrorEE(const Shape &e1_shape, const Shape &e2_shape, const size_t e1ID, const size_t e2ID, const double _PARALLEL) const
{
	const Eigen::Vector3d edge1_pos = e1_shape.VonE(e1ID, 0) - center;
	const Eigen::Vector3d edge2_pos = e2_shape.VonE(e2ID, 0) - center;
	const Eigen::Vector3d edge1_dic = e1_shape.VonE(e1ID, 1) - e1_shape.VonE(e1ID, 0);
	const Eigen::Vector3d edge2_dic = e2_shape.VonE(e2ID, 1) - e2_shape.VonE(e2ID, 0);
	const Eigen::Vector3d temp_3d_vec = edge1_dic.cross(edge2_dic);
	const double length = temp_3d_vec.norm();
	// parallel
	if ((length / edge1_dic.norm() / edge2_dic.norm()) < _PARALLEL) {
		return 0;
	}
	const Eigen::Vector3d norm = temp_3d_vec / length;
	return fabs(norm.dot(edge2_pos - edge1_pos));
}

void VisionErrorCorrector::ReduceSecondOrderEffect(Shape &moving_object, Shape &fixed_object, const ContactState &c_state)
{
	// v-f
	if (c_state.size() < 2) {
		return;
	}
	for (size_t i(0); i < c_state.size() - 1; ++i) {
		for (size_t j(i + 1); j < c_state.size(); ++j) {
			if (c_state[i].ContactType() == c_state[j].ContactType() && c_state[i].FirstElement().first == c_state[j].FirstElement().first) {
				switch (c_state[i].ContactType()) {
				case ContactElementBase::_VF_CONTACT:
#if 0
					if (c_state[i].FirstElement().first == 0) {
						ReduceSecondOrderEffectVF(moving_object, fixed_object, c_state[i].FirstElement().second, c_state[i].SecondElement().second, c_state[j].FirstElement().second, c_state[j].SecondElement().second, true);
					}
					else {
						ReduceSecondOrderEffectVF(fixed_object, moving_object, c_state[i].FirstElement().second, c_state[i].SecondElement().second, c_state[j].FirstElement().second, c_state[j].SecondElement().second, false);
					}
#endif
					break;
				case ContactElementBase::_EE_CONTACT:
					assert(c_state[i].FirstElement().first == 0);
					if (c_state[i].FirstElement().second == c_state[j].FirstElement().second) {
						ReduceSecondOrderEffectEE(moving_object, fixed_object, c_state[i].FirstElement().second, c_state[i].SecondElement().second, c_state[j].SecondElement().second, true);
					}
					else if (c_state[i].SecondElement().second == c_state[j].SecondElement().second){
						ReduceSecondOrderEffectEE(fixed_object, moving_object, c_state[i].SecondElement().second, c_state[i].FirstElement().second, c_state[j].FirstElement().second, false);
					}
					break;
				}
			}
		}
	}
	// adjust translation
	for (size_t i(0); i < c_state.size(); ++i) {
		SetEquation(moving_object, fixed_object, c_state[i]);
	}
#if 0
	std::ofstream ofs("debug_matrix.txt");
	for (size_t i(0); i < eq.size(); ++i) {
		ofs << eq[i].elm[0] << " " << eq[i].elm[1] << " " << eq[i].elm[2] << " " << eq[i].elm[3] << " "
			<< eq[i].elm[4] << " " << eq[i].elm[5] << " : " << eq[i].elm[6] << std::endl;
	}
	ofs.close();
#endif
	Eigen::Vector3d dt;
	SolveEquation(dt);
	if (dt.norm() <= _MIN_TRANS_DIST) {
		moving_object.SetTransformation(moving_object.Rot(), moving_object.Trans() + dt);
	}
}

void VisionErrorCorrector::ReduceSecondOrderEffectVF(Shape &v_shape, Shape &f_shape, const size_t v1, const size_t f1, const size_t v2, const size_t f2, bool v_move)
{
	// if faces parallel?
	if (fabs(f_shape.OuterNormal(f1).dot(f_shape.OuterNormal(f2))) > 1 - _NEARLY_ZERO) {
		const double dist = fabs(f_shape.OuterNormal(f1).dot(f_shape.VonF(f1, 0) - f_shape.VonF(f2, 0)));
		// v-v distance
		Eigen::Vector3d dv = v_shape.V(v1) - v_shape.V(v2);
		if (fabs(dv.norm() - dist) < _NEARLY_ZERO) {
			dv /= dv.norm();
			if (dv.dot(f_shape.OuterNormal(f1)) < 0) {
				dv = -dv;
			}
			if (v_move) {
				v_shape.Align(dv, f_shape.OuterNormal(f1));
			}
			else {
				f_shape.Align(f_shape.OuterNormal(f1), dv);
			}
#if 0
			// check
			std::cout << f_shape.OuterNormal(f1).transpose() << " - " << (v_shape.V(v1) - v_shape.V(v2)).normalized().transpose() << std::endl;
#endif
		}
	}
}

void VisionErrorCorrector::ReduceSecondOrderEffectEE(Shape &e1_shape, Shape &e2_shape, const size_t e1, const size_t e21, const size_t e22, bool e1_move)
{
	// if edges of e2_shape are parallel?
	Eigen::Vector3d e21_dirc = e2_shape.VonE(e21, 1) - e2_shape.VonE(e21, 0);
	e21_dirc /= e21_dirc.norm();
	Eigen::Vector3d e22_dirc = e2_shape.VonE(e22, 1) - e2_shape.VonE(e22, 0);
	e22_dirc /= e22_dirc.norm();
	if (fabs(e21_dirc.dot(e22_dirc)) > 1 - _NEARLY_ZERO) {
		Eigen::Vector3d ortho_dirc = OrthogonalDirection(e2_shape, e21, e22);
		Eigen::Vector3d e1_dirc = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
		// std::cout << ortho_dirc.dot(e21_dirc) << " " << ortho_dirc.dot(e22_dirc) << std::endl;
		if (fabs(e1_dirc.norm() - ortho_dirc.norm()) < _NEARLY_ZERO) {
			e1_dirc /= e1_dirc.norm();
			ortho_dirc /= ortho_dirc.norm();
			if (e1_dirc.dot(ortho_dirc) < 0) {
				e1_dirc =- e1_dirc;
			}
			if (e1_move) {
				e1_shape.Align(e1_dirc, ortho_dirc);
			}
			else {
				e2_shape.Align(ortho_dirc, e1_dirc);
			}
			// std::cout << OrthogonalDirection(e2_shape, e21, e22).normalized().transpose() << " - " << (e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0)).normalized().transpose() << std::endl;
		}
	}
}

Eigen::Vector3d VisionErrorCorrector::OrthogonalDirection(const Shape &e_shape, const size_t e1, const size_t e2)
{
	Eigen::Vector3d e_dirc = e_shape.VonE(e1, 1) - e_shape.VonE(e1, 0);
	e_dirc /= e_dirc.norm();
	const Eigen::Vector3d d = e_shape.VonE(e2, 0) - e_shape.VonE(e1, 0);
	const double s = e_dirc.dot(d);
	return e_dirc * s - d;
}

Eigen::MatrixXd VisionErrorCorrector::PseudoInverse(const Eigen::JacobiSVD<Eigen::MatrixXd> &svd)
{
	const double sig_max = svd.singularValues()[0];
	const double sig_thresh = sig_max * 1.0e-4; // _NEARLY_ZERO;
	Eigen::VectorXd inv_sig = svd.singularValues();
	// std::cout << inv_sig.transpose() << std::endl;
	for (int i(0); i < inv_sig.size(); ++i) {
		if (inv_sig[i] < sig_thresh) {
			inv_sig[i] = 0;
		}
		else {
			inv_sig[i] = 1.0 / inv_sig[i];
		}
	}
	return svd.matrixV() * inv_sig.asDiagonal() * svd.matrixU().transpose();
}