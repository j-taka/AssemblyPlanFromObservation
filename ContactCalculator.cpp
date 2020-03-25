// ContactCalculator.cpp

#include "ContactCalculator.h"

// const double _NEARLY_ZERO(1.0e-6);
const double _PARALLEL(0.05);

// only thresholding 
void ContactCalculator::Calc(const Shape &moving_object, const Shape &fixed_object)
{
	contact_elements.clear();
	// v-f contacts
	for (size_t i(0); i < moving_object.v_size(); ++i) {
		if (!moving_object.isVConvex(i)) {
			continue;
		}
		for (size_t j(0); j < fixed_object.f_size(); ++j) {
			const std::pair<size_t, size_t> O1(0, i);
			const std::pair<size_t, size_t> O2(1, j);
			const double dist = DistanceFromVertexToFace(i, j, moving_object, fixed_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_VF_CONTACT, dist));
			}
		}
	}
	// e-e contact
	for (size_t i(0); i < moving_object.e_size(); ++i) {
		if (!moving_object.isEConvex(i)) {
			continue;
		}
		for (size_t j(0); j < fixed_object.e_size(); ++j) {
			if (!fixed_object.isEConvex(j)) {
				continue;
			}
			const std::pair<size_t, size_t> O1(0, i);
			const std::pair<size_t, size_t> O2(1, j);
			const double dist = DistanceFromEdgeToEdge(i, j, moving_object, fixed_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_EE_CONTACT, dist));
			}
		}
	}
	// f-v contact
	for (size_t i(0); i < fixed_object.v_size(); ++i) {
		if (!fixed_object.isVConvex(i)) {
			continue;
		}
		for (size_t j(0); j < moving_object.f_size(); ++j) {
			const std::pair<size_t, size_t> O1(1, i);
			const std::pair<size_t, size_t> O2(0, j);
			const double dist = DistanceFromVertexToFace(i, j, fixed_object, moving_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_VF_CONTACT, dist));
			}
		}
	}
}

double ContactCalculator::DistanceFromVertexToFace(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape)
{
	const Eigen::Vector3d diff = v_shape.V(vID) - f_shape.VonF(fID, 0);
	double length = f_shape.OuterNormal(fID).dot(diff);
	// projection onto the face without boundary
	Eigen::Vector3d proj_v = v_shape.V(vID) - f_shape.OuterNormal(fID) * length;
	if (f_shape.GetVFRelation(proj_v, fID) == Shape::_OUTSIDE) {
		f_shape.PointNearFace(proj_v, fID);
		return (proj_v - v_shape.V(vID)).norm();
	}
	else {
		return length;
	}
}

double ContactCalculator::DistanceFromEdgeToEdge(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{
	// MyVector edge1, edge2;
	//double det;
	Eigen::Vector3d edge1 = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
	edge1.normalize();	
	Eigen::Vector3d edge2 = e2_shape.VonE(e2, 1) - e2_shape.VonE(e2, 0);
	edge2.normalize();
	const double det = edge1.dot(edge2); // cos
	if (fabs(det) >= 1 - _PARALLEL) { 
		// paralle 
		return DistanceBetweenParallelEdges(e1, e2, e1_shape, e2_shape);
	}
	else {
		// not parallel
		return DistanceBetweenNonParallelEdges(e1, e2, e1_shape, e2_shape);
	}
}


// distance between two parallel edges
double ContactCalculator::DistanceBetweenParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{	
	double length[4];
	const Eigen::Vector3d e11 = e1_shape.VonE(e1, 0);
	const Eigen::Vector3d e12 = e1_shape.VonE(e1, 1);
	const Eigen::Vector3d e21 = e2_shape.VonE(e2, 0);
	const Eigen::Vector3d e22 = e2_shape.VonE(e2, 1);
	Eigen::Vector3d tmp;
	length[0] = Shape::DistanceFromVertexToEdge(e11, e21, e22, tmp);
	length[1] = Shape::DistanceFromVertexToEdge(e12, e21, e22, tmp);
	length[2] = Shape::DistanceFromVertexToEdge(e21, e11, e12, tmp);
	length[3] = Shape::DistanceFromVertexToEdge(e22, e11, e12, tmp);
	double min_len = length[0];
	// find minimum
	for (int i(1); i < 4; ++i) {
		if (min_len > length[i]) {
			min_len = length[i];
		}
	}
	return min_len;
}

// 
double ContactCalculator::DistanceBetweenNonParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{
	const Eigen::Vector3d e1d = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
	const Eigen::Vector3d e2d = e2_shape.VonE(e2, 1) - e2_shape.VonE(e2, 0);
	const Eigen::Vector3d e12 = e2_shape.VonE(e2, 0) - e1_shape.VonE(e1, 0);
	// solve Equation
	Eigen::Matrix2d A;
	A(0, 0) = e1d.dot(e1d);
	A(1, 0) = e1d.dot(e2d);
	A(0, 1) = -A(1, 0);
	A(1, 1) = -e2d.dot(e2d);
	Eigen::Vector2d b;
	b[0] = e1d.dot(e12);
	b[1] = e2d.dot(e12);
	Eigen::Vector2d x = A.inverse() * b;
	if (x[0] < 0) { x[0] = 0; }
	if (x[0] > 1) { x[0] = 1; }
	if (x[1] < 0) { x[1] = 0; }
	if (x[1] > 1) { x[1] = 1; }
	// calculate distance
	return (e12 + e2d * x[1] - e1d * x[0]).norm();
}

