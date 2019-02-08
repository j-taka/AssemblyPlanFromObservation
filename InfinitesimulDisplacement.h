// InfinitesimulDisplacement.h
#pragma once

#include "ContactState.h"

class InfinitesimulDisplacement
{
public:
	typedef Eigen::Matrix<double, 6, 1> Screw;
	typedef std::vector<Screw, Eigen::aligned_allocator<Screw> > ScrewVector;
private:
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
	void Calculate(size_t objectID, const ContactState &src) {
		non_singular.clear();
		singular.clear();
		for (size_t i(0); i < src.size(); ++i) {
			Calculate(objectID, src[i]);
		}
	}
	bool isFeasible(const Screw &disp) const;
	//
	static Screw Translation(const Eigen::Vector3d &dic);
	static Screw Rotation(const Eigen::Vector3d &axis, const Eigen::Vector3d &center, double trans_rot_ratio = 0.0);

	friend std::ostream& operator<<(std::ostream &os, const InfinitesimulDisplacement &src);
private:
	void Calculate(size_t objectID, const ContactElement &src);
	static Screw Coefficient(const Eigen::Vector3d &pos, const Eigen::Vector3d &out_norm);
};