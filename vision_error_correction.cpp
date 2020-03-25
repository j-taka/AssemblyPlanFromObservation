// VisionErrorCorrection
#include <iostream>
#include "PoseListFileHandler.h"

#include "ContactCalculator.h"
#include "ContactElement.h"

#include "VisionErrorCorrector.h"

#undef _VC_DEBUG_MODE
size_t debug_target = 15;

static void PrintContact(const std::vector<ContactElementForErrorCorrection> &contacts)
{
	if (contacts.empty()) {
		std::cout << "No contact" << std::endl;
	}
	else {
		for (size_t i(0); i < contacts.size(); ++i) {
			std::cout << contacts[i] << std::endl;
		}
	}
}

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
#ifndef _VC_DEBUG_MODE
	for (size_t i(0); i < plHandler.length(); ++i) {
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
		std::cout << "Time: " << i + 1 << std::endl;
		PrintContact(cc.GetContact());
#ifdef _VC_DEBUG_MODE
		vc.SetVerbose(true);
#endif
		vc.Calc(objects[0], objects[1], cc.GetContact());
		std::cout << "Error: " << vc.CalculateMaximumError(objects[0], objects[1], cc.GetContact()) << std::endl;
		plHandler.SetTransformation(objects[0].Rot(), objects[0].Trans(), i, 0);
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