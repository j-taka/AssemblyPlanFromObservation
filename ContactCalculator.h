// ContactCalculator.h

#pragma once

#include "Shape.h"
#include "ContactStateForDisplacement.h"

class ContactCalculator
{
	// friend class TaskAnalyzer;
public:
	typedef std::vector<ContactElementForErrorCorrection> ContactState;
private:
	ContactState contact_elements;
	double threshold;
	struct ContactInfo {
		ContactElementForErrorCorrection elm;
		Eigen::Vector3d pos;
		bool parallel;
	};
public:
	// constructor
	ContactCalculator() : threshold(0.01) {}
	void Calc(const Shape &moving_object, const Shape &fixed_object);
	void DetailedAnalysis(const Shape &moving_object, const Shape &fixed_object);
	void Convert(ContactStateForDisplacement &cs, const Shape &moving_object, const Shape &fixed_object);
	void SetThreshold(double _threshold) {
		threshold = _threshold;
	}
	const ContactState& GetContact() const {
		return contact_elements;
	}
	void SetContact(const std::vector<ContactElementBase> &src) {
		contact_elements.resize(src.size());
		for (size_t i(0); i < src.size(); ++i) {
			contact_elements[i] = src[i];
		}
	}
	static void PrintContact(const ContactState &contacts);
private:
	static double DistanceFromVertexToFace(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape);
	double DistanceFromEdgeToEdge(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);
	double DistanceBetweenParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);
	double DistanceBetweenNonParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);

	void RemoveImpossibleContact(const Shape &moving_object, const Shape &fixed_object);
	bool isImpossibleContact(size_t vID, const Shape &v_shape, const Eigen::Vector3d &outer_normal, const double ang_th = 1.0e-5) const;
#if 0
	bool isImpossibleContactVF(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape, const double ang_th = 1.0e-5) const;
	bool isImpossibleContactEE(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape, const double ang_th = 1.0e-5) const;
#endif

	void SearchSingularContact(ContactStateForDisplacement &cs, const Shape &moving_object, const Shape &fixed_object);
	void CalculateContactInformation(std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const;

	bool ContactPositionEE(Eigen::Vector3d &pos, size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape, const double th = 1.0e-6) const;
	void SearchContactWiththeSampContactPosition(std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info) const;
	void LocalSingularContact(ContactStateForDisplacement &cs, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object);
	bool FindVVContact(const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object);
	int FindVEContact(const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object);
	static void Division(std::vector<size_t> &vf, std::vector<size_t> &fv, std::vector<size_t> &ee, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info);
	int FindVEContactPart(const std::vector<size_t> &vf, const std::vector<size_t> &ee, const std::vector<ContactInfo> &c_info, const Shape &v_shape, const Shape &e_shape, bool v_move);

	void AppendContacts(ContactStateForDisplacement &sc, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const;
	void AppendRelatedContacts(ContactStateForDisplacement &sc, const ContactElementForErrorCorrection &c_elm, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const;
	static void AppendContactEE(ContactStateForDisplacement &sc, const ContactInfo &c_info, const Shape &moving_object, const Shape &fixed_object);
		
#if 0	
	static double PointNearFace(Eigen::Vector3d &pos, const ShapeParser::Face &f, const Information &f_info);
#endif
};