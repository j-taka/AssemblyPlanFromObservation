// InfinitesimulDisplacement.h
#pragma once

#include "ContactStateForDisplacement.h"

class InfinitesimulDisplacement
{
public:
	typedef Eigen::Matrix<double, 6, 1> Screw;
	typedef std::vector<Screw, Eigen::aligned_allocator<Screw> > ScrewVector;
protected:
	ScrewVector non_singular;
	std::vector<ScrewVector> singular;
public:
	// constructor
	InfinitesimulDisplacement(){}
	// copy constructor
	InfinitesimulDisplacement(const InfinitesimulDisplacement &src) {
		non_singular = src.non_singular;
		singular = src.singular;
	}
	// convert
	void Calculate(size_t objectID, const ContactStateForDisplacement &src);
	bool isFeasible(const Screw &disp) const;
	bool isSingular() const {
		return (!singular.empty());
	}
	// 
	void SetNonSingularTranslationMatrix(Eigen::MatrixXd &dest) const;
	void SetNonSingularMatrix(Eigen::MatrixXd &dest) const;
	// 
	void SetRestrictedTranslationMatrix(Eigen::MatrixXd &dest) const;
	void SetRestrictedMatrix(Eigen::MatrixXd &dest) const;
	//
	static Screw Translation(const Eigen::Vector3d &dic);
	static Screw Rotation(const Eigen::Vector3d &axis, const Eigen::Vector3d &center, double trans_rot_ratio = 0.0);

	friend std::ostream& operator<<(std::ostream &os, const InfinitesimulDisplacement &src);
private:
	void Calculate(size_t objectID, const ContactElement &src);
	void Calculate(size_t objectID, const std::vector<ContactElement> &src);
	static Screw Coefficient(const Eigen::Vector3d &pos, const Eigen::Vector3d &out_norm);
	static bool InfinitesimulDisplacement::isSimilar(const Screw &s, const ScrewVector &sv);
};