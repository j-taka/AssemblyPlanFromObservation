// ContactCalculator.h

#pragma once

#include "Shape.h"
#include "ContactElement.h"

class ContactCalculator
{
private:
	std::vector<ContactElementForErrorCorrection> contact_elements;
	double threshold;
public:
	// constructor
	ContactCalculator() : threshold(0.01) {}
	void Calc(const Shape &moving_object, const Shape &fixed_object);
	void SetThreshold(double _threshold) {
		threshold = _threshold;
	}
	const std::vector<ContactElementForErrorCorrection> GetContact() const {
		return contact_elements;
	}
private:
	static double DistanceFromVertexToFace(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape);
	double DistanceFromEdgeToEdge(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);
	double DistanceBetweenParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);
	double DistanceBetweenNonParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape);
#if 0	
	static double PointNearFace(Eigen::Vector3d &pos, const ShapeParser::Face &f, const Information &f_info);
#endif
};