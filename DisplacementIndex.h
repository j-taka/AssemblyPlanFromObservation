// DisplacementIndex.h
#pragma once

#include "InfinitesimulDisplacement.h"

class DisplacementIndex
{
private:
	int t_m;
	int t_d;
	int t_c;
	int r_m;
	int r_d; // do not distinguish type I and II now.
	int r_c;
	// for singular
	int t_r;
	int r_r;
	double EIGEPS;
	double DIC_THRESH;
	// restricted information
	Eigen::MatrixXd tr_disp;
	Eigen::MatrixXd rr_disp;
	Eigen::MatrixXd ar_disp;
	size_t target_objectID; // to handle spacial case
public:
	// constructor
	DisplacementIndex() : EIGEPS(1.0e-6), DIC_THRESH(1.0e-2) {}
	// calculate 
	void Calculate(const InfinitesimulDisplacement &disp);
	//
	void CalculatePossibleTranslationDirections(Eigen::MatrixXd &_tr_disp, const InfinitesimulDisplacement &disp);

	int SetPossibleAxis(const ContactStateForDisplacement &c_state);
	void SwapPossibleAxis() {
		target_objectID = 1;
		Eigen::Vector3d swap = rr_disp.row(0).transpose();
		rr_disp.row(0) = rr_disp.row(1);
		rr_disp.row(1) = swap.transpose();
	}
	//
	int TranslactionRestrictedDOF() const {
		return t_r;
	}
	int RotationRestrictedDOF() const {
		return r_r;
	}
	const Eigen::MatrixXd& PossibleAxis() const {
		return rr_disp;
	}
	const Eigen::MatrixXd& PossibleTranslationDirection() const {
		return tr_disp;
	}
	// handle the spacial case
	size_t TargetObjectID() const {
		return target_objectID;
	}
	friend std::ostream& operator<<(std::ostream &os, const DisplacementIndex &src);
private:
	void TranslationPart(const InfinitesimulDisplacement &disp);
	void RotationAndRigidPart(const InfinitesimulDisplacement &disp);
	void SingularPart(const InfinitesimulDisplacement &disp);
	void CalculateMaintainingDisplacementInSingular(Eigen::MatrixXd &disp, const Eigen::MatrixXd &M) const;
	void CalculateRestrictionOfRotationAxis();
	Eigen::MatrixXd AxisRestrictionOne() const;
	Eigen::MatrixXd AxisRestrictionTwo() const;

	bool CheckPattern1(const ContactStateForDisplacement &c_state, size_t _target_objectID);
	bool CheckPattern2(const ContactStateForDisplacement &c_state);
	bool isOnlyVFsORFVs(size_t &count, Eigen::Vector3d &outer_normal, const std::vector<ContactElement> &c_elements, size_t target_objectID) const;
	bool isOuterNormalSame(const Eigen::Vector3d &outer_normal, const std::vector<ContactElement> &c_elements) const;
	bool arePointsOnTheLine(Eigen::Vector3d &direction, const ContactStateForDisplacement &c_state) const;
	Eigen::Vector3d GetOuterNormal(const ContactStateForDisplacement &c_state) const;
};