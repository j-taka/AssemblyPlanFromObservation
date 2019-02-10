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
public:
	// constructor
	DisplacementIndex() : EIGEPS(1.0e-6) {}
	// calculate 
	void Calculate(const InfinitesimulDisplacement &disp);
	//
	friend std::ostream& operator<<(std::ostream &os, const DisplacementIndex &src);
private:
	void TranslationPart(const InfinitesimulDisplacement &disp);
	void RotationAndRigidPart(const InfinitesimulDisplacement &disp);
	void SingularPart(const InfinitesimulDisplacement &disp);
	void CalculateMaintainingDisplacementInSingular(Eigen::MatrixXd &disp, const Eigen::MatrixXd &M) const;
};