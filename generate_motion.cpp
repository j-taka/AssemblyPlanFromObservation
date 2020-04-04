// generate_motion.cpp

#include "MotionGenerator.h"
#include <boost/program_options.hpp>

int main(int argc, char **argv)
{
	// options
	namespace po = boost::program_options;
	po::options_description opt("Option");
	opt.add_options()
		("random,r", "add random noise to the initial pose of a moving object")
		("help,h", "show help");
	po::options_description hidden("Hidden option");
	hidden.add_options()
		("input-file", boost::program_options::value<std::vector<std::string> >(), "input file");

	po::options_description cmdline_options;
	cmdline_options.add(opt).add(hidden);

	po::positional_options_description p;
	p.add("input-file", -1);
	// parse
	try {
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
		po::notify(vm);
		if (vm.count("help") == 1 || (vm.count("input-file") == 0 && vm["input-file"].as<std::vector<std::string> >().size() != 3)) {
			std::cerr << argv[0] << " (output pose list) (input pose list) (init pose list)" << std::endl;
			std::cerr << opt << std::endl;
			return -1;
		}
		PoseListFileHandler plHandler;
		int err;
		const std::string plFile = vm["input-file"].as<std::vector<std::string> >()[1];
		if ((err = plHandler.Load(plFile)) != 0) {
			switch (err) {
			case PoseListFileHandler::_FILENOTOPEN:
				std::cerr << "Cannot open: " << plFile << std::endl;
				return -1;
			case PoseListFileHandler::_FORMATERROR:
				std::cerr << "Format error: " << plFile << std::endl;
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
		// first analyze
		for (size_t i(0); i < plHandler.length(); ++i) {
			for (size_t j(0); j < objects.size(); ++j) {
				Eigen::Matrix3d R;
				Eigen::Vector3d t;
				plHandler.GetTransformation(R, t, i, j);
				objects[j].SetTransformation(R, t);
			}
			cc.Calc(objects[0], objects[1]);
			cc.DetailedAnalysis(objects[0], objects[1]);
			// 
			ta.AppendContactState(objects, cc.GetContact());
		}
		MotionGenerator mg;
		PoseListFileHandler for_init;
		const std::string initFile = vm["input-file"].as<std::vector<std::string> >()[2];
		if ((err = for_init.Load(initFile)) != 0) {
			switch (err) {
			case PoseListFileHandler::_FILENOTOPEN:
				std::cerr << "Cannot open: " << initFile << std::endl;
				return -1;
			case PoseListFileHandler::_FORMATERROR:
				std::cerr << "Format error: " << initFile << std::endl;
				return -1;
			default:
				std::cerr << "Cannot open: " << for_init.ErroredModelFile() << std::endl;
				std::cerr << err << std::endl;
				return -1;
			}
		}
		Eigen::Matrix3d R;
		Eigen::Vector3d t;
		for_init.GetTransformation(R, t, 0, 0);
		TaskAnalyzer::Pose pose;
		pose.block(0, 0, 3, 3) = R;
		pose.block(0, 3, 3, 1) = t;
		// set init pose
		if (vm.count("random") == 1) {
			srand((unsigned int)time(NULL));
			TaskAnalyzer::ContactState init(ta[0]);
			init.poses[0] = pose;
			TaskAnalyzer::RecalculateContactState(init, objects[0], objects[1]);
			mg.SetRangeInTranslation(1.0);
			mg.RandomPoseGeneration(pose, objects[0], objects[1], init);
		}
		// motion generation
		for_init.GetTransformation(R, t, 0, 1);
		ta.SetPoseOfFixedObject(R, t);
		ta.Analyze(objects[0], objects[1]);
		std::cout << ta << std::endl;
		std::vector<TaskAnalyzer::Pose> trajectory;
		trajectory.push_back(pose);
		for (size_t i(1); i < ta.NumberOfStates(); ++i) {
			std::cout << i << " -> " << i + 1 << std::endl;
			mg.CalculateTrajectory(trajectory, objects[0], objects[1], ta[i - 1], ta[i]);
		}
		// save
		PoseListFileHandler save(plHandler);
		save.SetPoses(trajectory, objects[1].Rot(), objects[1].Trans());
		std::string saveFile = vm["input-file"].as<std::vector<std::string> >()[0];
		if (save.Save(saveFile) != 0) {
			std::cerr << "Cannot save: " << saveFile << std::endl;
			return -1;
		}
		return 0;
	}
	catch (std::exception &e) {
		std::cout << e.what() << std::endl;
		return -1;
	}
}