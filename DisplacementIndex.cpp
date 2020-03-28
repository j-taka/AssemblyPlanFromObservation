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
		t_m = tm_disp.rows();
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
		r_m = am_disp.rows() - t_m;
		const int a_md = am_disp.rows() + all_pcc.DDOF();
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
	if (!disp.isSingular()) {
		t_r = t_d + t_c;
		r_r = r_d + r_c;
	}
	else {
		Eigen::MatrixXd r_transM;
		disp.SetRestrictedTranslationMatrix(r_transM);
		// std::cout << r_transM << std::endl;
		Eigen::MatrixXd tr_disp;
		CalculateMaintainingDisplacementInSingular(tr_disp, r_transM);
		Eigen::MatrixXd r_allM;
		disp.SetRestrictedMatrix(r_allM);
		// std::cout << r_allM << std::endl;
		Eigen::MatrixXd ar_disp;
		CalculateMaintainingDisplacementInSingular(ar_disp, r_allM);
		// std::cout << ar_disp << std::endl;
		// index
		t_r = 3 - tr_disp.rows();
		r_r = 6 - ar_disp.rows() - t_r;
	}
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
