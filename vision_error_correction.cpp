// VisionErrorCorrection
#include <iostream>
#include "PoseListFileHandler.h"

#include "ContactCalculator.h"
#include "ContactElement.h"

#include "VisionErrorCorrector.h"

#undef _VC_DEBUG_MODE
size_t debug_target = 15;

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
#if 0
	for (size_t i(0); i < objects.size(); ++i) {
		std::cout << "Object " << i + 1 << std::endl;
		objects[i].SetTransformation(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
		std::cout << objects[i] << std::endl;
	}
	return 0;
#endif
	// Set Contact Calculator
	ContactCalculator cc;
	cc.SetThreshold(5.0); // for error correction
	VisionErrorCorrector vc;
	VisionErrorCorrector::FixFixedObject(plHandler);
#ifndef _VC_DEBUG_MODE
	for (int i(0), time(0); i < plHandler.length(); ++i, ++time) {
#else
	{
		size_t i(debug_target);
#endif
		for (size_t j(0); j < objects.size(); ++j) {
			Eigen::Matrix3d R;
			Eigen::Vector3d t;
			plHandler.GetTransformation(R, t, i, j);
			objects[j].SetTransformation(R, t);
		}
		cc.Calc(objects[0], objects[1]);
		std::cout << "Time: " << time + 1 << std::endl;
		// vc.SetVerbose(true);
		std::vector<ContactElementForErrorCorrection> c_state = cc.GetContact();
		// ContactCalculator::PrintContact(c_state);
		const bool success = vc.Calc(objects[0], objects[1], c_state);
		std::cout << "Error: " << vc.CalculateMaximumError(objects[0], objects[1], c_state) << std::endl;
		if (success) {
			plHandler.SetTransformation(objects[0].Rot(), objects[0].Trans(), i, 0);
		}
		else {
			plHandler.EraseIthPose(i);
			i--;
		}
#ifdef _VC_DEBUG_MODE
		PoseListFileHandler for_debug(plHandler);
		for_debug.SetPoses(vc.Convergence(), objects[1].Rot(), objects[1].Trans());
		for_debug.Save("debug.pose.list");
		return 0;
/*
#else
		if (i == debug_target) {
			PoseListFileHandler for_debug(plHandler);
			for_debug.SetPoses(vc.Convergence(), objects[1].Rot(), objects[1].Trans());
			for_debug.Save("debug.pose.list");
			return 0;
		}
*/
#endif
	}
	// save
	if (plHandler.Save(argv[1]) != 0) {
		std::cerr << "Cannot save: " << argv[1] << std::endl;
		return -1;
	}
	return 0;
}