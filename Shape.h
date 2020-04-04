// Shape.h

#pragma once

#include "ShapeParser.h"

class Shape : private ShapeParser
{
public:
	enum VFRelation { _INSIDE, _ON_EDGE, _OUTSIDE };
private:
	Eigen::Matrix3d R;
	Eigen::Vector3d t;
	std::vector<Eigen::Vector3d> trans_pos;
	std::vector<Eigen::Vector3d> trans_edge_outside;
	std::vector<Eigen::Vector3d> trans_norm;
public:
	const Eigen::Matrix3d& Rot() const {
		return R;
	}
	const Eigen::Vector3d& Trans() const {
		return t;
	}
	size_t v_size() const {
		return ShapeParser::v_size();
	}
	size_t e_size() const {
		return ShapeParser::e_size();
	}
	size_t f_size() const {
		return ShapeParser::f_size();
	}
	bool isVConvex(size_t src) const {
		return vertices[src].convex;
	}
	bool isEConvex(size_t src) const {
		return edges[src].convex;
	}
	void Set(const TopoDS_Shape &src) {
		ShapeParser::Set(src);
	}
	void SetTransformation(const Eigen::Matrix3d &_R, const Eigen::Vector3d &_t);
	const Eigen::Vector3d& V(size_t vID) const {
		return trans_pos[vID];
	}
	const std::vector<size_t>& AdjacentEdges(size_t vID) const {
		return vertices[vID].adjacent_edges;
	}
	void Align(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, double _NEARLY_ZERO = 1.0e-6);
	const size_t& vIDonE(size_t eID, size_t order) const {
		if (order == 0) {
			return edges[eID].v1;
		}
		else {
			return edges[eID].v2;
		}
	}
	const Eigen::Vector3d& VonE(size_t eID, size_t order) const {
		if (order == 0) {
			return trans_pos[edges[eID].v1];
		}
		else {
			return trans_pos[edges[eID].v2];
		}
	}
	const size_t& fIDAdjacentToE(size_t eID, size_t order) const {
		return edges[eID].adjacent_faces[order];
	}
	size_t NumberOfVerticesOnF(size_t fID) const {
		return faces[fID].boundary_vertices.size();
	}
 	const Eigen::Vector3d& VonF(size_t fID, size_t order) const {
		return trans_pos[faces[fID].boundary_vertices[order]];
	}
	const Eigen::Vector3d& OuterNormal(size_t fID) const {
		return trans_norm[fID];
	}
	const Eigen::Vector3d &OutsideDirectionOfE(size_t eID) const {
		return trans_edge_outside[eID];
	}
	VFRelation GetVFRelation(const Eigen::Vector3d &proj_v, size_t fID, double _NEARLY_ZERO = 1.0e-6) const;
	double Shape::PointNearFace(Eigen::Vector3d &pos, size_t fID) const;
	static double DistanceFromVertexToEdge(const Eigen::Vector3d &pos, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2, Eigen::Vector3d &ans);
	friend std::ostream& operator<<(std::ostream &os, const Shape &src);
private:
	static double AngleBetweenVectors(const Eigen::Vector3d &pos, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &normal);
	// 
	static bool isOnEdge(const Eigen::Vector3d &proj_v, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, double _NEARLY_ZERO = 1.0e-6);

	static void RemoveSmallerValues(Eigen::Matrix3d &R, double smaller_value = 1.0e-6);
};