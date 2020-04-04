// analyze_task.cpp
#include <iostream>
#include "PoseListFileHandler.h"

#include "ContactCalculator.h"
#include "ContactElement.h"

#include "TaskAnalyzer.h"

int main(int argc, char **argv)
{
	if (argc != 3) {
		std::cerr << argv[0] << " (output pose list file) (input pose list file)" << std::endl;
		return -1;
	}
	PoseListFileHandler plHandler;
	int err;
	if ((err = plHandler.Load(argv[2])) != 0) {
		switch (err) {
		case PoseListFileHandler::_FILENOTOPEN:
			std::cerr << "Cannot open: " << argv[2] << std::endl;
			return -1;
		case PoseListFileHandler::_FORMATERROR:
			std::cerr << "Format error: " << argv[2] << std::endl;
			return -1;
		default:
			std::cerr << "Cannot open: " << plHandler.ErroredModelFile() << std::endl;
			std::cerr << err << std::endl;
			return -1;
		}
	}
	// parse shape information
	std::vector<Shape> objects(plHandler.NumberOfObjects());
	for (size_t i(0); i < objects.size(); ++i) {
		objects[i].Set(plHandler.Object(i));
	}
	// Set Contact Calculator
	ContactCalculator cc;
	cc.SetThreshold(0.01); // for task analyzer
	TaskAnalyzer ta;
	for (int i(0), time(0); i < plHandler.length(); ++i, ++time) {
		for (size_t j(0); j < objects.size(); ++j) {
			Eigen::Matrix3d R;
			Eigen::Vector3d t;
			plHandler.GetTransformation(R, t, i, j);
			objects[j].SetTransformation(R, t);
		}
		cc.Calc(objects[0], objects[1]);
		cc.DetailedAnalysis(objects[0], objects[1]);
		ta.AppendContactState(objects, cc.GetContact());
	}
	ta.Analyze(objects[0], objects[1]);
	std::cout << ta << std::endl;
	// save
	PoseListFileHandler save(plHandler);
	std::vector<Eigen::Matrix<double, 3, 4> > poses(ta.NumberOfStates());
	for (size_t i(0); i < ta.NumberOfStates(); ++i) {
		poses[i] = ta[i].poses[0];
	}
	save.SetPoses(poses, objects[1].Rot(), objects[1].Trans());
	if (save.Save(argv[1]) != 0) {
		std::cerr << "Cannot save: " << argv[1] << std::endl;
	}
	return 0;
}