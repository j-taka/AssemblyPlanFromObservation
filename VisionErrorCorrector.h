// VisionErrorCorrector.h

#pragma once 

#include "PoseListFileHandler.h"
#include "Shape.h"
#include "ContactElement.h"

// #define _VC_DEBUG_MODE

class VisionErrorCorrector
{
public:
	typedef std::vector<ContactElementForErrorCorrection> ContactState;

private:
	struct Each {
		double elm[7];
	};
	std::vector<Each> eq;
	double tolerance;
	Eigen::Vector3d center;
	int _MAX_LOOP;
	int _MAX_NUMBER_OF_REMOVE_ELEMENT;
	double _MIN_TRANS_DIST;

	bool verbose;
#ifdef _VC_DEBUG_MODE
	std::vector<Eigen::Matrix<double, 3, 4> > convergence;
#endif
public:
	// constructor
	VisionErrorCorrector() : tolerance(1.0e-5), verbose(false), _MAX_LOOP(40), _MAX_NUMBER_OF_REMOVE_ELEMENT(2), _MIN_TRANS_DIST(10) {}
	bool Calc(Shape &moving_object, Shape &fixed_object, ContactState &c_state);
	double CalculateMaximumError(const Shape &moving_object, const Shape &fixed_object, const ContactState &c_state);

	void SetVerbose(bool src) {
		verbose = src;
	}
#ifdef _VC_DEBUG_MODE
	const std::vector<Eigen::Matrix<double, 3, 4> >& Convergence() const {
		return convergence;
	}
#endif
	static void FixFixedObject(PoseListFileHandler &plHandler);

private:
	bool CalcEach(Shape &moving_object, Shape &fixed_object, const ContactState &c_state);

	void CalcCenter(const Shape &moving_object, const Shape &fixed_object, const ContactState &c_state);

	double SetEquation(const Shape &moving_object, const Shape &fixed_object, const ContactElementForErrorCorrection &c_element);
	double SetEquationVF(const Shape &v_shape, const Shape &f_shape, const size_t vID, const size_t fID, bool v_move);
	double SetEquationEE(const Shape &e1_shape, const Shape &e2_shape, const size_t e1ID, const size_t e2ID, const double _PARALLEL = 0.05);

	void SolveEquation(Eigen::Matrix3d &dR, Eigen::Vector3d &dt);
	void SolveEquation(Eigen::Vector3d &dt); // only translation

	double CalcError(const Shape &moving_object, const Shape &fixed_object, const ContactElementForErrorCorrection &c_element) const;
	double CalcErrorVF(const Shape &v_shape, const Shape &f_shape, const size_t vID, const size_t fID) const;
	double CalcErrorEE(const Shape &e1_shape, const Shape &e2_shape, const size_t e1ID, const size_t e2ID, const double _PARALLEL = 0.05) const;

	void ReduceSecondOrderEffect(Shape &moving_object, Shape &fixed_object, const ContactState &c_state);
	void ReduceSecondOrderEffectVF(Shape &v_shape, Shape &f_shape, const size_t v1, const size_t f1, const size_t v2, const size_t f2, bool v_move);
	void ReduceSecondOrderEffectEE(Shape &e1_shape, Shape &e2_shape, const size_t e1, const size_t e21, const size_t e22, bool e1_move);

	static Eigen::Vector3d OrthogonalDirection(const Shape &e_shape, const size_t e1, const size_t e2);
	static Eigen::MatrixXd VisionErrorCorrector::PseudoInverse(const Eigen::JacobiSVD<Eigen::MatrixXd> &svd);
};