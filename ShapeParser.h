// ShapeParser.h

#pragma once
#include <vector>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <Eigen/Dense>

class ShapeParser
{
public:
	struct Vertex {
		Eigen::Vector3d pos;
		std::vector<size_t> adjacent_edges;
		bool convex;
	};
	struct Edge {
		size_t v1, v2;
		Eigen::Vector3d outer_direction;
		std::vector<size_t> adjacent_faces;
		bool convex;
	};
	struct Face {
		std::vector<size_t> boundary_vertices;
		Eigen::Vector3d outer_normal;
	};
protected:
	std::vector<TopoDS_Vertex> aTopoV;
	std::vector<TopoDS_Edge> aTopoE;
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<Face> faces;
public:
	void Set(const TopoDS_Shape &shape, bool detailed = true);
	size_t v_size() const {
		return vertices.size();
	}
	const Vertex& V(size_t src) const {
		return vertices[src];
	}
	size_t e_size() const {
		return edges.size();
	}
	const Edge& E(size_t src) const {
		return edges[src];
	}
	size_t f_size() const {
		return faces.size();
	}
	const Face& F(size_t src) const {
		return faces[src];
	}
	friend std::ostream& operator<<(std::ostream &os, const ShapeParser &src);
private:
	void VertexIdentification(const TopoDS_Shape &shape);
	void EdgeIdentification(const TopoDS_Shape &shape);
	void FaceIdentification(const TopoDS_Shape &shape);
	void SetEdgeConvexity(const TopoDS_Shape &shape);
	void SetVertexConvexity();
	size_t GetVertexID(const TopoDS_Vertex &vertex) const;
	size_t GetEdgeID(const TopoDS_Edge &edge) const;

	void AdjacentEdgesOnVertices();
};