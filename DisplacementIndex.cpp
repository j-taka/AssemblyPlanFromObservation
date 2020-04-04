// DisplacementIndex.cpp

#include "DisplacementIndex.h"
#include "PCCCalculator.h"

void DisplacementIndex::Calculate(const InfinitesimulDisplacement &disp)
{
	// 
	TranslationPart(disp);
	//
	RotationAndRigidPart(disp);
	// singular
	SingularPart(disp);
}

void DisplacementIndex::CalculatePossibleTranslationDirections(Eigen::MatrixXd &_tr_disp, const InfinitesimulDisplacement &disp)
{
	Eigen::MatrixXd r_transM;
	disp.SetRestrictedTranslationMatrix(r_transM);
	// std::cout << r_transM << std::endl;
	CalculateMaintainingDisplacementInSingular(_tr_disp, r_transM);
}

void DisplacementIndex::TranslationPart(const InfinitesimulDisplacement &disp)
{
	// set matrix
	Eigen::MatrixXd transM;
	disp.SetNonSingularTranslationMatrix(transM);
	if (transM.rows() == 0) {
		// set index
		t_m = 3;
		t_d = 0;
		t_c = 0;
	}
	else {
		PCCCalculator pcc;
		Eigen::MatrixXd tm_disp;
		Eigen::MatrixXd td_disp;
		// std::cout << transM << std::endl;
		pcc.Dual(tm_disp, td_disp, transM);
		// set index
		t_m = static_cast<int>(tm_disp.rows());
		t_d = pcc.DDOF();
		t_c = 3 - t_m - t_d;
	}
}

void DisplacementIndex::RotationAndRigidPart(const InfinitesimulDisplacement &disp)
{
	// set matrix
	Eigen::MatrixXd allM;
	disp.SetNonSingularMatrix(allM);
	if (allM.rows() == 0) {
		r_m = 3;
		r_c = 0;
		r_d = 0;
	}
	else {
		PCCCalculator all_pcc;
		Eigen::MatrixXd am_disp;
		Eigen::MatrixXd ad_disp;
		all_pcc.Dual(am_disp, ad_disp, allM);
		// set index
		r_m = static_cast<int>(am_disp.rows() - t_m);
		const int a_md = static_cast<int>(am_disp.rows() + all_pcc.DDOF());
		const int a_c = 6 - a_md;
		r_c = a_c - t_c;
		r_d = (3 - r_m - r_c);
	}
#if 0
	// axis
	PCCCalculator axis_pcc;
	Eigen::MatrixXd rm_disp;
	Eigen::MatrixXd rd_disp;
	std::cout << am_disp.block(0, 0, am_disp.rows(), 3) << std::endl;
	std::cout << ad_disp.block(0, 0, ad_disp.rows(), 3) << std::endl;
	axis_pcc.Dual(rm_disp, rd_disp, am_disp.block(0, 0, am_disp.rows(), 3), ad_disp.block(0, 0, ad_disp.rows(), 3));
#endif
}

void DisplacementIndex::SingularPart(const InfinitesimulDisplacement &disp)
{
	Eigen::MatrixXd r_transM;
	disp.SetRestrictedTranslationMatrix(r_transM);
	// std::cout << r_transM << std::endl;
	CalculateMaintainingDisplacementInSingular(tr_disp, r_transM);
	Eigen::MatrixXd r_allM;
	disp.SetRestrictedMatrix(r_allM);
	// std::cout << r_allM << std::endl;
	CalculateMaintainingDisplacementInSingular(ar_disp, r_allM);
	// std::cout << ar_disp << std::endl;
	// index
	t_r = 3 - static_cast<int>(tr_disp.rows());
	r_r = 6 - static_cast<int>(ar_disp.rows()) - t_r;
	// 
	CalculateRestrictionOfRotationAxis();
}

void DisplacementIndex::CalculateMaintainingDisplacementInSingular(Eigen::MatrixXd &disp, const Eigen::MatrixXd &M) const
{

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M.transpose() * M);
	// already sorted
	int dim(0);
	for (int i(0); i < eig.eigenvalues().size(); ++i) {
		if (eig.eigenvalues()[i] < EIGEPS) {
			dim++;
		}
		else {
			break;
		}
	}
	disp = eig.eigenvectors().block(0, 0, M.cols(), dim).transpose();
}

void DisplacementIndex::CalculateRestrictionOfRotationAxis()
{
	switch (r_r) {
	case 0:
		rr_disp = Eigen::MatrixXd::Identity(3, 3);
		break;
	case 1:
		rr_disp = AxisRestrictionTwo(); // should consider the second order effect 
		break;
	case 2:
		rr_disp = AxisRestrictionOne();
		break;
	case 3:
		rr_disp = Eigen::MatrixXd(0, 3); // nothing
		break;
	}
}

Eigen::MatrixXd DisplacementIndex::AxisRestrictionOne() const
{
	double max_val(0);
	int max_ID(0);
	for (int i(0); i < ar_disp.rows(); ++i) {
		Eigen::VectorXd tmp = ar_disp.row(i).transpose();
		tmp.normalize();
		if (tmp.block(0, 0, 3, 1).norm() > max_val) {
			max_val = tmp.block(0, 0, 3, 1).norm();
			max_ID = i;
		}
	}
	Eigen::MatrixXd res = ar_disp.block(max_ID, 0, 1, 3);
	res /= res.norm();
	return res;
}

// use eigen
Eigen::MatrixXd DisplacementIndex::AxisRestrictionTwo() const
{
	Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
	for (int i(0); i < ar_disp.rows(); ++i) {
		Eigen::VectorXd tmp = ar_disp.row(i).transpose();
		// std::cout << tmp.transpose() << std::endl;
		tmp.normalize();
		// std::cout << tmp.block(0, 0, 3, 1).norm() << std::endl;
		if (tmp.block(0, 0, 3, 1).norm() > DIC_THRESH) {
			Eigen::Vector3d tmp2 = tmp.block(0, 0, 3, 1);
			tmp2.normalize();
			cov += tmp2 * tmp2.transpose();
		}
	}
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
	// std::cout << eig.eigenvalues().transpose() << std::endl;
	return eig.eigenvectors().block(0, 1, 3, 2).transpose();
}

int DisplacementIndex::SetPossibleAxis(const ContactStateForDisplacement &c_state)
{
	bool check(false);
	int f_objectID(0);
	// std::cout << c_state << std::endl;
	// pattern 1 (multiple vf's)
	if (CheckPattern1(c_state, 0)) { return 0; }
	if (CheckPattern1(c_state, 1)) { return 0; }
	// pattern 2 
	if (CheckPattern2(c_state)) { return 1; }
	return -1;
}

// more than one vf's, no fv, contact positions on the line
bool DisplacementIndex::CheckPattern1(const ContactStateForDisplacement &c_state, size_t _target_objectID)
{
	size_t count(0);
	Eigen::Vector3d outer_normal;
	if (!isOnlyVFsORFVs(count, outer_normal, c_state.elements, _target_objectID)) {
		return false;
	}
	for (size_t i(0); i < c_state.singular_elements.size(); ++i) {
		if (!isOnlyVFsORFVs(count, outer_normal, c_state.singular_elements[i], _target_objectID)) {
			return false;
		}
	}
	if (count == 0) {
		return false;
	}
	if (!isOuterNormalSame(outer_normal, c_state.elements)) {
		return false;
	}
	for (size_t i(0); i < c_state.singular_elements.size(); ++i) {
		if (!isOuterNormalSame(outer_normal, c_state.elements)) {
			return false;
		}
	}
	Eigen::Vector3d direction;
	if (!arePointsOnTheLine(direction, c_state)) {
		return false;
	}
	// set
	target_objectID = _target_objectID;
	if (_target_objectID == 0) {
		rr_disp.row(0) = outer_normal.transpose();
		rr_disp.row(1) = direction.transpose();
	}
	else {
		rr_disp.row(0) = direction.transpose();
		rr_disp.row(1) = outer_normal.transpose();
	}
	return true;
}

bool DisplacementIndex::isOnlyVFsORFVs(size_t &count, Eigen::Vector3d &outer_normal, const std::vector<ContactElement> &c_elements, size_t target_objectID) const
{
	for (size_t i(0); i < c_elements.size(); ++i) {
		if (c_elements[i].ContactType() == ContactElementBase::_VF_CONTACT) {
			if (c_elements[i].FirstElement().first == target_objectID) {
				outer_normal = c_elements[i].OuterNormal();
				count++;
			}
			else {
				return false;
			}
		}
	}
	return true;
}

bool DisplacementIndex::isOuterNormalSame(const Eigen::Vector3d &outer_normal, const std::vector<ContactElement> &c_elements) const
{
	for (size_t i(0); i < c_elements.size(); ++i) {
		if ((outer_normal.cross(c_elements[i].OuterNormal())).norm() > DIC_THRESH) {
			return false;
		}
	}
	return true;
}

bool DisplacementIndex::arePointsOnTheLine(Eigen::Vector3d &direction, const ContactStateForDisplacement &c_state) const
{
	Eigen::Vector3d ave = Eigen::Vector3d::Zero();
	size_t count(0);
	for (size_t i(0); i < c_state.elements.size(); ++i) {
		ave += c_state.elements[i].ContactPosition();
		count++;
	}
	for (size_t i(0); i < c_state.singular_elements.size(); ++i) {
		for (size_t j(0); j < c_state.singular_elements[i].size(); ++j) {
			ave += c_state.singular_elements[i][j].ContactPosition();
			count++;
		}
	}
	ave /= static_cast<double>(count);
	Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
	for (size_t i(0); i < c_state.elements.size(); ++i) {
		const Eigen::Vector3d tmp = c_state.elements[i].ContactPosition() - ave;
		cov += tmp * tmp.transpose();
	}
	for (size_t i(0); i < c_state.singular_elements.size(); ++i) {
		for (size_t j(0); j < c_state.singular_elements[i].size(); ++j) {
			const Eigen::Vector3d tmp = c_state.singular_elements[i][j].ContactPosition() - ave;
			cov += tmp * tmp.transpose();
		}
	}
	cov /= static_cast<double>(count);
	// Eigen decomposition
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
	// already sorted
	// std::cout << eig.eigenvalues().transpose() << std::endl;
	const double thresh = eig.eigenvalues()[2] * EIGEPS;
	if (fabs(eig.eigenvalues()[0]) > thresh || fabs(eig.eigenvalues()[1]) > thresh) {
		return false;
	}
	direction = eig.eigenvectors().block(0, 2, 3, 1);
	return true;
}

bool DisplacementIndex::CheckPattern2(const ContactStateForDisplacement &c_state)
{
	const Eigen::Vector3d outer_normal = GetOuterNormal(c_state);
	if (!isOuterNormalSame(outer_normal, c_state.elements)) {
		return false;
	}
	for (size_t i(0); i < c_state.singular_elements.size(); ++i) {
		if (!isOuterNormalSame(outer_normal, c_state.elements)) {
			return false;
		}
	}
	Eigen::Vector3d direction;
	if (!arePointsOnTheLine(direction, c_state)) {
		return false;
	}
	// return only candidate
	target_objectID = 0; 
	rr_disp.row(0) = outer_normal.transpose();
	rr_disp.row(1) = direction.transpose();
	return true;
}

Eigen::Vector3d DisplacementIndex::GetOuterNormal(const ContactStateForDisplacement &c_state) const
{
	if (c_state.elements.empty()) {
		return c_state.singular_elements[0][0].OuterNormal();
	}
	else {
		return c_state.elements[0].OuterNormal();
	}
}

std::ostream& operator<<(std::ostream &os, const DisplacementIndex &src)
{
	os << "Translation" << std::endl;
	os << "Maintain DOF: " << src.t_m << std::endl;
	os << "Detaching DOF: " << src.t_d << std::endl;
	os << "Constraint DOF: " << src.t_c << std::endl;
	os << "Rectricted DOF: " << src.t_r << std::endl;
	os << "Rotation" << std::endl;
	os << "Maintain DOF: " << src.r_m << std::endl;
	os << "Detaching DOF: " << src.r_d << std::endl;
	os << "Constraint DOF: " << src.r_c << std::endl;
	os << "Restricted DOF: " << src.r_r << std::endl;
	return os;
}
