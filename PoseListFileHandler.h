// PoseListFileHandler
#pragma once

#include <vector>

#include <TopoDS.hxx>
#include <gp_TrsfForm.hxx>
#include "ShapeParser.h"

class PoseListFileHandler
{
public:
	static const int _FILENOTOPEN = -1;
	static const int _FORMATERROR = -2;
	static const int _LOADMODELERR = -3;
	static const int _FAILTOSAVE = -4;

private:
	std::vector<std::string> object_names;
	std::vector<TopoDS_Shape> aTopoObjects;

	size_t time;
	typedef std::vector<gp_Trsf> PosesInEachFrame;
	std::vector<PosesInEachFrame> poses;

	// error 
	std::string errored_model_file;
public:
	int Load(const std::string &filename);
	int Save(const std::string &filename);
	std::string ErroredModelFile() const {
		return errored_model_file;
	}
	size_t length() const {
		return poses.size();
	}
	size_t NumberOfObjects() const {
		return aTopoObjects.size();
	}
	void ObjectsAtTimeT(std::vector<TopoDS_Shape> &objects, const size_t t) const;
	const TopoDS_Shape& Object(size_t src) const {
		return aTopoObjects[src];
	}
	void GetTransformation(Eigen::Matrix3d &R, Eigen::Vector3d &t, const size_t time, const size_t objectID) const;
	void SetTransformation(const Eigen::Matrix3d &R, const Eigen::Vector3d &t, const size_t time, const size_t objectID);
	void EraseIthPose(const size_t time) {
		poses.erase(poses.begin() + time);
	}
	// for debugging 
	void SetPoses(const std::vector<Eigen::Matrix<double, 3, 4> > &poses, const Eigen::Matrix3d &R, const Eigen::Vector3d &t);
private: 
	int LoadModels();
};