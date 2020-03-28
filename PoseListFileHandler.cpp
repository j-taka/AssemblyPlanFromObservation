// PoseListFileHandler.cpp

#include "PoseListFileHandler.h"
#include <fstream>

#include <STEPControl_Reader.hxx>
#include <gp_Quaternion.hxx>
#include <NCollection_Mat4.hxx>

#include <BRepBuilderAPI_Transform.hxx>

int PoseListFileHandler::Load(const std::string &file)
{
	std::ifstream ifs(file);
	if (!ifs) {
		return _FILENOTOPEN;
	}
	try {
		std::string dummy;
		size_t num;
		ifs >> dummy >> dummy >> dummy >> num;
		object_names.resize(num);
		for (size_t i(0); i < object_names.size(); ++i) {
			ifs >> object_names[i];
		}
		ifs >> dummy >> dummy >> dummy >> num;
		poses.resize(num);
		for (size_t i(0); i < poses.size(); ++i) {
			poses[i].resize(object_names.size());
			for (size_t j(0); j < poses[i].size(); ++j) {
				Standard_Real x, y, z, w;
				// location
				ifs >> dummy >> x >> y >> z;
				gp_Vec v(x, y, z);
				// vector-angle 
				ifs >> x >> y >> z >> w;
				gp_Quaternion q;
				q.SetVectorAndAngle(gp_Vec(x, y, z), w);
				poses[i][j].SetTransformation(q, v);
			}
		}
		ifs.close();
		return LoadModels();
	}
	catch (...) {
		object_names.clear();
		poses.clear();
		return _FORMATERROR;
	}
}

int PoseListFileHandler::Save(const std::string &file)
{
	std::ofstream ofs(file);
	if (!ofs) {
		return _FILENOTOPEN;
	}
	try {
		ofs << "No of Models: " << object_names.size() << std::endl;
		for (size_t i(0); i < object_names.size(); ++i) {
			ofs << object_names[i] << std::endl;
		}
		ofs << std::endl;
		ofs << "No of Scenes: " << poses.size() << std::endl;
		ofs << std::setprecision(12); 
		for (size_t i(0); i < poses.size(); ++i) {
			for (size_t j(0); j < poses[i].size(); ++j) {
				ofs << (j == 0 ? "2 " : "0 ");
				// location
				ofs << poses[i][j].TranslationPart().X() << " "
					<< poses[i][j].TranslationPart().Y() << " "
					<< poses[i][j].TranslationPart().Z() << " ";
				// axis angle
				gp_Vec axis;
				Standard_Real angle;
				poses[i][j].GetRotation().GetVectorAndAngle(axis, angle);
				ofs << axis.X() << " " << axis.Y() << " " << axis.Z() << " " << angle << std::endl;
			}
			ofs << std::endl;
		}
		ofs.close();
		return 0;
	}
	catch (...) {
		return _FAILTOSAVE;
	}
}

int PoseListFileHandler::LoadModels()
{
	aTopoObjects.resize(object_names.size());
	for (size_t i(0); i < object_names.size(); ++i) {
		STEPControl_Reader reader;
		IFSelect_ReturnStatus stat = reader.ReadFile(object_names[i].c_str());
		if (stat != IFSelect_ReturnStatus::IFSelect_RetDone) {
			errored_model_file = object_names[i];
			aTopoObjects.clear();
			return _LOADMODELERR;
		}
		reader.TransferRoot();
		aTopoObjects[i] = reader.Shape();
	}
	return 0;
}

void PoseListFileHandler::ObjectsAtTimeT(std::vector<TopoDS_Shape> &objects, const size_t t) const
{
	objects.resize(aTopoObjects.size());
	for (size_t i(0); i < aTopoObjects.size(); ++i) {
		// transform
		BRepBuilderAPI_Transform transform(poses[t][i]);
		transform.Perform(aTopoObjects[i]);
		objects[i] = transform.Shape();
	}
}

void PoseListFileHandler::GetTransformation(Eigen::Matrix3d &R, Eigen::Vector3d &t, const size_t time, const size_t objectID) const
{
	NCollection_Mat4<double> mat;
	poses[time][objectID].GetMat4(mat);
	for (int r(0); r < 3; ++r) {
		for (int c(0); c < 3; ++c) {
			R(r, c) = mat.GetValue(r, c);
		}
		t(r) = mat.GetValue(r, 3);
	}
#if 0
	for (size_t r(0); r < 4; ++r) {
		for (size_t c(0); c < 4; ++c) {
			std::cout << mat.GetValue(r, c) << " ";
		}
		std::cout << std::endl;
	}
#endif
}


void PoseListFileHandler::SetTransformation(const Eigen::Matrix3d &R, const Eigen::Vector3d &t, const size_t time, const size_t objectID)
{
	poses[time][objectID].SetValues(
		R(0, 0), R(0, 1), R(0, 2), t(0),
		R(1, 0), R(1, 1), R(1, 2), t(1),
		R(2, 0), R(2, 1), R(2, 2), t(2));
#if 0
	for (size_t r(0); r < 4; ++r) {
		for (size_t c(0); c < 4; ++c) {
			std::cout << mat.GetValue(r, c) << " ";
		}
		std::cout << std::endl;
	}
#endif
}

void PoseListFileHandler::SetPoses(const std::vector<Eigen::Matrix<double, 3, 4> > &_poses, const Eigen::Matrix3d &R, const Eigen::Vector3d &t)
{
	poses.resize(_poses.size());
	for (size_t i(0); i < _poses.size(); ++i) {
		poses[i].resize(2);
		SetTransformation(_poses[i].block(0, 0, 3, 3), _poses[i].block(0, 3, 3, 1), i, 0);
		SetTransformation(R, t, i, 1);
	}
}

