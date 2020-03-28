// ShapeParser.cpp

#include "ShapeParser.h"
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp_Face.hxx>

void ShapeParser::Set(const TopoDS_Shape &shape, bool detailed)
{
	// vertex
	VertexIdentification(shape);
	// edge
	EdgeIdentification(shape);
	// face
	FaceIdentification(shape);
	// convex
	SetEdgeConvexity(shape);
	SetVertexConvexity();
	// adjacency
	AdjacentEdgesOnVertices();
	// free memory
	aTopoV.clear();
	aTopoE.clear();
	return;
}

size_t ShapeParser::GetVertexID(const TopoDS_Vertex &vertex) const
{
	for (size_t i(0); i < aTopoV.size(); ++i) {
		if (vertex.IsSame(aTopoV[i])) {
			return i;
		}
	}
	return aTopoV.size();
}

size_t ShapeParser::GetEdgeID(const TopoDS_Edge &edge) const
{
	for (size_t i(0); i < aTopoE.size(); ++i) {
		if (edge.IsSame(aTopoE[i])) {
			return i;
		}
	}
	return aTopoE.size();
}

void ShapeParser::VertexIdentification(const TopoDS_Shape &shape)
{
	vertices.clear();
	for (TopExp_Explorer ex(shape, TopAbs_VERTEX); ex.More(); ex.Next()) {
		TopoDS_Vertex vertex = TopoDS::Vertex(ex.Current());
		// find...
		if (GetVertexID(vertex) == vertices.size()) {
			aTopoV.push_back(vertex);
			Vertex v;
			v.pos = Eigen::Vector3d(BRep_Tool::Pnt(vertex).X(), BRep_Tool::Pnt(vertex).Y(), BRep_Tool::Pnt(vertex).Z());
			vertices.push_back(v);
		}
	}
#if 0
	std::cout << vertices.size() << std::endl;
	for (size_t i(0); i < vertices.size(); ++i) {
		std::cout << vertices[i].transpose() << std::endl;
	}
#endif
}

void ShapeParser::EdgeIdentification(const TopoDS_Shape &shape)
{
	edges.clear();
	for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next()) {
		TopoDS_Edge edge = TopoDS::Edge(ex.Current());
		// find...
		size_t i(0);
		if (GetEdgeID(edge) == edges.size()){
			aTopoE.push_back(edge);
			ShapeParser::Edge e;
			e.v1 = GetVertexID(TopExp::FirstVertex(edge));
			e.v2 = GetVertexID(TopExp::LastVertex(edge));
			edges.push_back(e);
		}
	}
#if 0
	std::cout << edges.size() << " " << aTopoE.size() << std::endl;
	for (size_t i(0); i < edges.size(); ++i) {
		std::cout << edges[i].v1 << " - " << edges[i].v2 << " " << std::endl;
	}
#endif
}

void ShapeParser::FaceIdentification(const TopoDS_Shape &shape)
{
	faces.clear();
	for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
		TopoDS_Face face = TopoDS::Face(ex.Current());
		TopExp_Explorer ex1(face, TopAbs_EDGE);
		TopoDS_Edge edge = TopoDS::Edge(ex1.Current());
		const size_t v1 = GetVertexID(TopExp::FirstVertex(edge));
		const size_t v2 = GetVertexID(TopExp::LastVertex(edge));
		ex1.Next();
		edge = TopoDS::Edge(ex1.Current());
		const size_t v3 = GetVertexID(TopExp::FirstVertex(edge));
		const size_t v4 = GetVertexID(TopExp::LastVertex(edge));
		size_t tmp;
		Face f;
		if (v1 == v3 || v1 == v4) {
			f.boundary_vertices.push_back(v2);
			tmp = v1;
		}
		else {
			f.boundary_vertices.push_back(v1);
			tmp = v2;
		}
		for (; ex1.More(); ex1.Next()) {
			TopoDS_Edge edge = TopoDS::Edge(ex1.Current());
			const size_t v5 = GetVertexID(TopExp::FirstVertex(edge));
			const size_t v6 = GetVertexID(TopExp::LastVertex(edge));
			f.boundary_vertices.push_back(tmp);
			if (tmp == v5) {
				tmp = v6;
			}
			else {
				tmp = v5;
			}
		}
		// calculate outer normal
		BRepGProp_Face prop(face);
		Standard_Real u1, u2, s1, s2;
		prop.Bounds(u1, u2, s1, s2);
		double u = (u1 + u2) / 2;
		double s = (s1 + s2) / 2;
		gp_Vec vec;
		gp_Pnt pnt;
		prop.Normal(u, s, pnt, vec);
		f.outer_normal = Eigen::Vector3d(vec.X(), vec.Y(), vec.Z());
		faces.push_back(f);
	}
#if 0
	// debug
	std::cout << faces.size() << std::endl;
	for (size_t i(0); i < faces.size(); ++i) {
		for (size_t j(0); j < faces[i].boundary_vertices.size(); ++j) {
			std::cout << faces[i].boundary_vertices[j] << " ";
		}
		std::cout << std::endl << faces[i].outer_normal.transpose() << std::endl;
	}
#endif
}

void ShapeParser::SetEdgeConvexity(const TopoDS_Shape &shape)
{
	// find adjacent faces
	for (size_t i(0); i < edges.size(); ++i) {
		edges[i].adjacent_faces.resize(2, faces.size());
	}
	size_t i(0);
	for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next(), ++i) {
		TopoDS_Face face = TopoDS::Face(ex.Current());
		size_t j(0);
		for (TopExp_Explorer ex1(face, TopAbs_EDGE); ex1.More(); ex1.Next(), ++j) {
			TopoDS_Edge edge = TopoDS::Edge(ex1.Current());
			const size_t eID = GetEdgeID(edge);
#if 0
			std::cout << faces[i].boundary_vertices[j] << "-" << faces[i].boundary_vertices[(j + 1) % faces[i].boundary_vertices.size()] << " ";
			std::cout << edges[eID].v1 << "-" << edges[eID].v2 << std::endl;
#endif
			if (edges[eID].v1 == faces[i].boundary_vertices[j]) {
				edges[eID].adjacent_faces[0] = i;
			}
			else {
				edges[eID].adjacent_faces[1] = i;
			}
		}
	}
	for (size_t i(0); i < edges.size(); ++i) {
		const Eigen::Vector3d ed = vertices[edges[i].v1].pos - vertices[edges[i].v2].pos;
		const Eigen::Vector3d n1 = faces[edges[i].adjacent_faces[0]].outer_normal;
		const Eigen::Vector3d n2 = faces[edges[i].adjacent_faces[1]].outer_normal;
		edges[i].convex = (ed.dot(n1.cross(n2)) < 0);
		edges[i].outer_direction = n1 + n2;
		edges[i].outer_direction.normalize();
	}
#if 0
	// debug
	for (size_t i(0); i < edges.size(); ++i) {
		std::cout << edges[i].v1 << " - " << edges[i].v2 << " " << edges[i].convex << std::endl;
	}
#endif
}

void ShapeParser::SetVertexConvexity()
{
	for (size_t i(0); i < vertices.size(); ++i) {
		vertices[i].convex = true;
	}
	for (size_t i(0); i < edges.size(); ++i) {
		if (!edges[i].convex) {
			vertices[edges[i].v1].convex = false;
			vertices[edges[i].v2].convex = false;
		}
	}
#if 0
	// debug
	for (size_t i(0); i < vertices.size(); ++i) {
		std::cout << vertices[i].pos.transpose() << " " << vertices[i].convex << std::endl;
	}
#endif
}

void ShapeParser::AdjacentEdgesOnVertices()
{
	for (size_t i(0); i < vertices.size(); ++i) {
		vertices[i].adjacent_edges.clear();
	}
	for (size_t i(0); i < edges.size(); ++i) {
		vertices[edges[i].v1].adjacent_edges.push_back(i);
		vertices[edges[i].v2].adjacent_edges.push_back(i);
	}
}

std::ostream& operator<<(std::ostream &os, const ShapeParser &src)
{
	os << "Number of vertices: " << src.vertices.size() << std::endl;
	for (size_t i(0); i < src.vertices.size(); ++i) {
		os << "V" << i << " : " << src.vertices[i].pos.transpose() << " " << (src.vertices[i].convex ? "(convex)" : "(not convex)") << std::endl;
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
		os << ") (outer normal; " << src.faces[i].outer_normal.transpose() << ")" << std::endl;
	}
	return os;
}

