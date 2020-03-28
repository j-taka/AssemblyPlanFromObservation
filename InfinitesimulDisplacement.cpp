// InfinitesimulDisplacement.cpp

#include "InfinitesimulDisplacement.h"

const double _NEARLY_ZERO = 1.0e-6;

void InfinitesimulDisplacement::Calculate(size_t objectID, const ContactStateForDisplacement &src) 
{
	non_singular.clear();
	singular.clear();
	for (size_t i(0); i < src.elements.size(); ++i) {
		Calculate(objectID, src.elements[i]);
	}
	// singular
	for (size_t i(0); i < src.singular_elements.size(); ++i) {
		Calculate(objectID, src.singular_elements[i]);
	}

}

void InfinitesimulDisplacement::Calculate(size_t objectID, const ContactElement &src) 
{
	if (src.FirstElement().first == objectID) {
		non_singular.push_back(Coefficient(src.ContactPosition(), src.OuterNormal()));
	}
	if (src.SecondElement().first == objectID) {
		non_singular.push_back(Coefficient(src.ContactPosition(), -src.OuterNormal()));
	}
}

void InfinitesimulDisplacement::Calculate(size_t objectID, const std::vector<ContactElement> &src)
{
	ScrewVector tmp;
	for (size_t i(0); i < src.size(); ++i) {
		if (src[i].FirstElement().first == objectID) {
			const Screw s = Coefficient(src[i].ContactPosition(), src[i].OuterNormal());
			if (!isSimilar(s, tmp)) {
				tmp.push_back(s);
			}
		}
		if (src[i].SecondElement().first == objectID) {
			const Screw s = Coefficient(src[i].ContactPosition(), -src[i].OuterNormal());
			if (!isSimilar(s, tmp)) {
				tmp.push_back(s);
			}
		}
	}
	if (tmp.empty()) {
		return;
	}
	if (tmp.size() == 1) {
		non_singular.push_back(tmp[0]);
	}
	else {
		singular.push_back(tmp);
	}
}

bool InfinitesimulDisplacement::isSimilar(const Screw &s, const ScrewVector &sv)
{
	for (size_t i(0); i < sv.size(); ++i) {
		if (s.dot(sv[i]) / s.norm() / sv[i].norm() > 1 - _NEARLY_ZERO) {
			return true;
		}
	}
	return false;
}

bool InfinitesimulDisplacement::isFeasible(const Screw &disp) const
{
	for (size_t i(0); i < non_singular.size(); ++i) {
		if (non_singular[i].dot(disp) < -_NEARLY_ZERO) {
			return false;
		}
	}
	for (size_t i(0); i < singular.size(); ++i) {
		bool check(false);
		for (size_t j(0); j < singular[i].size(); ++j) {
			if (singular[i][j].dot(disp) >= -_NEARLY_ZERO) {
				check = true;
				break;
			}
		}
		if (!check) {
			return false;
		}
	}
	return true;
}

InfinitesimulDisplacement::Screw InfinitesimulDisplacement::Coefficient(const Eigen::Vector3d &pos, const Eigen::Vector3d &out_norm)
{
	Screw res;
	res.block(0, 0, 3, 1) = pos.cross(out_norm);
	res.block(3, 0, 3, 1) = out_norm;
	return res;
}

InfinitesimulDisplacement::Screw InfinitesimulDisplacement::Translation(const Eigen::Vector3d &dic)
{
	Screw res;
	res.block(0, 0, 3, 1) = Eigen::Vector3d::Zero(); 
	res.block(3, 0, 3, 1) = dic;
	return res;
}

InfinitesimulDisplacement::Screw InfinitesimulDisplacement::Rotation(const Eigen::Vector3d &axis, const Eigen::Vector3d &center, double trans_rot_ratio)
{
	Screw res;
	res.block(0, 0, 3, 1) = axis;
	res.block(3, 0, 3, 1) = center.cross(axis) + trans_rot_ratio * axis;
	return res;
}

std::ostream& operator<<(std::ostream &os, const InfinitesimulDisplacement &src)
{
	os << "Non-singular part" << std::endl;
	if (src.non_singular.empty()) {
		os << "Nothing" << std::endl;
	}
	else {
		for (size_t i(0); i < src.non_singular.size(); ++i) {
			os << src.non_singular[i].transpose() << std::endl;
		}
	}
	if (!src.singular.empty()) {
		os << "Singular part" << std::endl;
		for (size_t i(0); i < src.singular.size(); ++i) {
			os << "Part " << i + 1 << std::endl;
			for (size_t j(0); j < src.singular[i].size(); ++j) {
				os << src.singular[i][j] << std::endl;
			}
		}
	}
	return os;
}

void InfinitesimulDisplacement::SetNonSingularTranslationMatrix(Eigen::MatrixXd &dest) const
{
	dest = Eigen::MatrixXd(non_singular.size(), 3);
	for (int r(0); r < dest.rows(); ++r) {
		dest.row(r) = non_singular[r].block(3, 0, 3, 1).transpose();
	}
}

void InfinitesimulDisplacement::SetNonSingularMatrix(Eigen::MatrixXd &dest) const
{
	dest = Eigen::MatrixXd(non_singular.size(), 6);
	for (int r(0); r < dest.rows(); ++r) {
		dest.row(r) = non_singular[r].transpose();
	}
}

void InfinitesimulDisplacement::SetRestrictedTranslationMatrix(Eigen::MatrixXd &dest) const
{
	// count 
	size_t r = non_singular.size();
	for (size_t i(0); i < singular.size(); ++i) {
		r += singular[i].size();
	}
	dest = Eigen::MatrixXd(r, 3);
	r = 0;
	for (size_t i(0); i < non_singular.size(); ++i, ++r) {
		dest.row(r) = non_singular[i].block(3, 0, 3, 1).transpose();
	}
	for (size_t i(0); i < singular.size(); ++i) {
		for (size_t j(0); j < singular[i].size(); ++j, ++r) {
			dest.row(r) = singular[i][j].block(3, 0, 3, 1).transpose();
		}
	}
}

void InfinitesimulDisplacement::SetRestrictedMatrix(Eigen::MatrixXd &dest) const
{
	// count 
	size_t r = non_singular.size();
	for (size_t i(0); i < singular.size(); ++i) {
		r += singular[i].size();
	}
	dest = Eigen::MatrixXd(r, 6);
	r = 0;
	for (size_t i(0); i < non_singular.size(); ++i, ++r) {
		dest.row(r) = non_singular[i].transpose();
	}
	for (size_t i(0); i < singular.size(); ++i) {
		for (size_t j(0); j < singular[i].size(); ++j, ++r) {
			dest.row(r) = singular[i][j].transpose();
		}
	}
}
