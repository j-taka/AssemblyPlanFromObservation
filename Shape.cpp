// shape.cpp

#include "Shape.h"

void Shape::SetTransformation(const Eigen::Matrix3d &_R, const Eigen::Vector3d &_t)
{
	R = _R;
	t = _t;
	RemoveSmallerValues(R);
	trans_pos.resize(vertices.size());
	for (size_t i(0); i < vertices.size(); ++i) {
		trans_pos[i] = R * vertices[i].pos + t;
	}
	trans_edge_outside.resize(edges.size());
	for (size_t i(0); i < edges.size(); ++i) {
		trans_edge_outside[i] = R * edges[i].outer_direction;
	}
	trans_norm.resize(faces.size());
	for (size_t i(0); i < faces.size(); ++i) {
		trans_norm[i] = R * faces[i].outer_normal;
	}
}

void Shape::Align(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, double _NEARLY_ZERO)
{
	Eigen::Vector3d axis = v1.cross(v2);
	const double angle = asin(axis.norm());
	if (axis.norm() < _NEARLY_ZERO) {
		return;
	}
	axis /= axis.norm();
	Eigen::Matrix3d dR;
	dR = Eigen::AngleAxisd(angle, axis);
	SetTransformation(dR * R, t);
}

void Shape::RemoveSmallerValues(Eigen::Matrix3d &R, double smaller_value)
{
	for (int r(0); r < 3; ++r) {
		for (int c(0); c < 3; ++c) {
			if (fabs(R(r, c)) < smaller_value) {
				R(r, c) = 0.0;
			}
		}
	}
}

Shape::VFRelation Shape::GetVFRelation(const Eigen::Vector3d &proj_v, size_t fID, double _NEARLY_ZERO) const
{
	double sum_ang(0);
	size_t memo = faces[fID].boundary_vertices.size() - 1;
	for (size_t i(0); i < faces[fID].boundary_vertices.size(); ++i) {
		const Eigen::Vector3d p1 = VonF(fID, memo);
		const Eigen::Vector3d p2 = VonF(fID, i);
		if (isOnEdge(proj_v, p1, p2)) {
			return _ON_EDGE;
		}
		sum_ang += AngleBetweenVectors(proj_v, p1, p2, trans_norm[fID]);
		memo = i;
	}
	if (fabs(sum_ang) >= M_PI) {
		return _INSIDE; // considering small errors
	}
	else {
		return _OUTSIDE;
	}
}

bool Shape::isOnEdge(const Eigen::Vector3d &proj_v, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, double _NEARLY_ZERO)
{
	const Eigen::Vector3d dirc = p2 - p1;
	const double tmp1 = proj_v.dot(dirc) - p1.dot(dirc);
	const double tmp2 = dirc.squaredNorm();
	const double t = tmp1 / tmp2;
	const Eigen::Vector3d on_edge = p1 + t * dirc;
	const Eigen::Vector3d diff = proj_v - on_edge;
	return (diff.norm() <= _NEARLY_ZERO && (t >= 0 && t <= 1));
}

static double proper_trifunc(double src)
{
	if (src < -1) {
		return -1;
	}
	if (src > 1) {
		return 1;
	}
	return src;
}

double Shape::AngleBetweenVectors(const Eigen::Vector3d &pos, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &normal)
{
	const Eigen::Vector3d va = v2 - pos;
	const Eigen::Vector3d vb = v1 - pos;
	const Eigen::Vector3d vc = vb.cross(va);
#if 0
	const double quot = vc.norm() / (va.norm() * vb.norm()); // sin of the angle
	const double ang_alpha = asin(proper_trifunc(quot));
	const double dot1 = va.dot(vb);
	const double dot2 = normal.dot(vc);
	const double ang = (dot2 < 0 ? -1 : 1);
	return (dot1 >= 0 ? ang * ang_alpha : ang * M_PI - ang_alpha);
#endif
	const double len = va.norm() * vb.norm();
	const double y = vc.dot(normal) / len;
	const double x = va.dot(vb) / len;
	return atan2(y, x);
}

double Shape::PointNearFace(Eigen::Vector3d &pos, size_t fID) const
{
	//int i, memo;
	//double length, templen;
	//MyVector ev1, ev2, ans, temp;
	Eigen::Vector3d ans(pos);
	double length = DBL_MAX; // start from the bigger value
	size_t memo = faces[fID].boundary_vertices.size() - 1;
	for (size_t i(0); i < faces[fID].boundary_vertices.size(); ++i) {
		const Eigen::Vector3d ev1 = VonF(fID, memo);
		const Eigen::Vector3d ev2 = VonF(fID, i);
		Eigen::Vector3d temppos;
		const double templen = DistanceFromVertexToEdge(pos, ev1, ev2, temppos);
		if (length > templen) {
			ans = temppos;
			length = templen;
		}
		memo = i;
	}
	pos = ans;
	return length;
}

double Shape::DistanceFromVertexToEdge(const Eigen::Vector3d &pos, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2, Eigen::Vector3d &ans)
{
	//MyVector dic, ep;
	//double t, temp, temp2;
	const Eigen::Vector3d ep = pos - e1;
	const Eigen::Vector3d dic = e2 - e1;
	const double temp = dic.dot(ep);
	const double temp2 = dic.squaredNorm();
	double t = temp / temp2;
	if (t < 0) {
		t = 0;
	}
	if (t > 1) {
		t = 1;
	}
	ans = e1 + dic * t;
	return (ans - pos).norm();
}

std::ostream& operator<<(std::ostream &os, const Shape &src)
{
	os << "Number of vertices: " << src.vertices.size() << std::endl;
	for (size_t i(0); i < src.vertices.size(); ++i) {
		os << "V" << i << " : " << src.trans_pos[i].transpose() << " " << (src.vertices[i].convex ? "(convex)" : "(not convex)") << std::endl;
	}
	os << "Number of edges: " << src.edges.size() << std::endl;
	for (size_t i(0); i < src.edges.size(); ++i) {
		os << "E" << i << " : (" << src.edges[i].v1 << " - " << src.edges[i].v2 << ") " << (src.edges[i].convex ? "(convex)" : "(not convex)") << std::endl;
	}
	os << "Number of faces: " << src.faces.size() << std::endl;
	for (size_t i(0); i < src.faces.size(); ++i) {
		os << "F" << i << " : (" << src.faces[i].boundary_vertices[0];

		for (size_t j(1); j < src.faces[i].boundary_vertices.size(); ++j) {
			os << " - " << src.faces[i].boundary_vertices[j];
		}
		os << ") (outer normal; " << src.trans_norm[i].transpose() << ")" << std::endl;
	}
	return os;
}

