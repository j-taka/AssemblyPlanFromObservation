// APO_Recognition_sample.cpp

#include "InfinitesimulDisplacement.h"
#include "DisplacementIndex.h"

// parameters
const size_t O1 = 0;
const size_t O2 = 1;
const std::pair<size_t, size_t> O1_V1(0, 0);
const std::pair<size_t, size_t> O1_V2(0, 1);
const std::pair<size_t, size_t> O2_F1(1, 0);

static std::string answer(bool src)
{
	return (src ? "yes" : "no");
}

int main(int argc, char **argv)
{
	ContactState c_state;
	c_state.AddContact(ContactElement::VFContact(O1_V1, O2_F1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1)));
	c_state.AddContact(ContactElement::VFContact(O1_V2, O2_F1, Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 1)));
	// print
	std::cout << "Contact state:" << std::endl;
	std::cout << c_state << std::endl;
	//
	InfinitesimulDisplacement i_disp;
	i_disp.Calculate(O1, c_state);
	std::cout << "Feasible displacement:" << std::endl;
	std::cout << i_disp << std::endl;
	// check
	std::cout << "Translation (0, 0, 1): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Translation(Eigen::Vector3d::UnitZ()))) << std::endl;
	std::cout << "Translation (0, 0, -1): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Translation(-Eigen::Vector3d::UnitZ()))) << std::endl;
	std::cout << "Translation (0, 1, 0): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Translation(-Eigen::Vector3d::UnitY()))) << std::endl;
	std::cout << "Rotion X around (0, 0, 0): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Rotation(Eigen::Vector3d::UnitX(), Eigen::Vector3d::Zero()))) << std::endl;
	std::cout << "Rotion X around (0, 1, 0): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Rotation(Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY()))) << std::endl;
	std::cout << "Rotion X around (0, -1, 0): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Rotation(Eigen::Vector3d::UnitX(), -Eigen::Vector3d::UnitY()))) << std::endl;
	std::cout << "Rotion Z around (0, 0, 0): " << answer(i_disp.isFeasible(InfinitesimulDisplacement::Rotation(Eigen::Vector3d::UnitZ(), Eigen::Vector3d::Zero()))) << std::endl;
	//
	DisplacementIndex ind;
	ind.Calculate(i_disp);
	std::cout << "Indices of displacement" << std::endl;
	std::cout << ind << std::endl;
	return 0;
}