// random_pose_generation.cpp
// this may be useful for randomized path planning

#include "MotionGenerator.h"
#include "PoseListFileHandler.h"

static const int NUMBER_OF_SAMPLE = 100;

int main(int argc, char **argv)
{
	int target_state_ID;
	int res;
	if (argc != 4 || (res = sscanf(argv[2], "%d", &target_state_ID)) == 0) {
		std::cerr << argv[0] << " (output pose list) (target state ID) (input pose list)" << std::endl;
		return -1;
	}
	PoseListFileHandler plHandler;
	int err;
	if ((err = plHandler.Load(argv[3])) != 0) {
		switch (err) {
		case PoseListFileHandler::_FILENOTOPEN:
			std::cerr << "Cannot open: " << argv[3] << std::endl;
			return -1;
		case PoseListFileHandler::_FORMATERROR:
			std::cerr << "Format error: " << argv[3] << std::endl;
			return -1;
		default:
			std::cerr << "Cannot open: " << plHandler.ErroredModelFile() << std::endl;
			std::cerr << err << std::endl;
			return -1;
		}
	}
	if (target_state_ID < 0 || target_state_ID >= plHandler.length()) {
		std::cout << "Out of range: " << target_state_ID << " (0 to " << plHandler.length() - 1 << ")" << std::endl;
		return -1;
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
	// first analyze
	for (size_t j(0); j < objects.size(); ++j) {
		Eigen::Matrix3d R;
		Eigen::Vector3d t;
		plHandler.GetTransformation(R, t, target_state_ID, j);
		objects[j].SetTransformation(R, t);
	}
	cc.Calc(objects[0], objects[1]);
	cc.DetailedAnalysis(objects[0], objects[1]);
	// 
	ta.AppendContactState(objects, cc.GetContact());
	ta.Analyze(objects[0], objects[1]);
	std::cout << ta << std::endl;
	// 
	MotionGenerator mg;
	std::vector<TaskAnalyzer::Pose> poses(NUMBER_OF_SAMPLE);
	for (int i(0); i < NUMBER_OF_SAMPLE; ++i) {
		mg.RandomPoseGeneration(poses[i], objects[0], objects[1], ta[0]);
	}
	// save
	PoseListFileHandler save(plHandler);
	save.SetPoses(poses, objects[1].Rot(), objects[1].Trans());
	if (save.Save(argv[1]) != 0) {
		std::cerr << "Cannot save: " << argv[1] << std::endl;
		return -1;
	}
	return 0;
}