PROJECT(AssemblyPlanFromObservation)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

#eigen
set(CMAKE_MODULE_PATH $ENV{Eigen3_DIR}/cmake)
set(EIGEN3_INCLUDE_DIR $ENV{Eigen3_DIR})
FIND_PACKAGE(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#qhull
set(QHULL_PATH "C:/Program Files/qhull")
include_directories(${QHULL_PATH}/include)
link_directories(${QHULL_PATH}/lib)

#PCC
set(PCC_PATH ../PolyhedralConvexCones)
include_directories(${PCC_PATH})

# tbb
set(TBB_PATH "C:/OpenCASCADE-7.3.0-vc14-64/tbb_2017.0.100")
link_directories(${TBB_PATH}/lib/intel64/vc10)

# occ
set(OpenCASCADE_DIR C:/OpenCASCADE-7.3.0-vc14-64/opencascade-7.3.0/cmake)
find_package(OpenCASCADE REQUIRED)
include_directories(${OpenCASCADE_INCLUDE_DIR})
link_directories(${OpenCASCADE_LIBRARY_DIR})

ADD_EXECUTABLE(
	APO_Recognition_sample
	APO_Recognition_sample.cpp
	ContactElement.cpp
	ContactElement.h
	ContactStateForDisplacement.cpp
	ContactStateForDisplacement.h
	InfinitesimulDisplacement.cpp
	InfinitesimulDisplacement.h
	DisplacementIndex.cpp
	DisplacementIndex.h
	${PCC_PATH}/PCCCalculator.cpp
	${PCC_PATH}/PCCCalculator.h
	${PCC_PATH}/Combination.cpp
	${PCC_PATH}/Combination.h
)

ADD_EXECUTABLE(
	VisionErrorCorrection
	vision_error_correction.cpp
	PoseListFileHandler.cpp
	PoseListFileHandler.h
	ShapeParser.cpp
	ShapeParser.h
	Shape.cpp
	Shape.h
	ContactCalculator.cpp
	ContactCalculator.h
	ContactElement.cpp
	ContactElement.h
	VisionErrorCorrector.cpp
	VisionErrorCorrector.h
)

ADD_EXECUTABLE(
	TaskAnalyzer
	analyze_task.cpp
	PoseListFileHandler.cpp
	PoseListFileHandler.h
	ShapeParser.cpp
	ShapeParser.h
	Shape.cpp
	Shape.h
	ContactCalculator.cpp
	ContactCalculator.h
	ContactElement.cpp
	ContactElement.h
	TaskAnalyzer.cpp
	TaskAnalyzer.h
	ContactStateForDisplacement.cpp
	ContactStateForDisplacement.h
	InfinitesimulDisplacement.cpp
	InfinitesimulDisplacement.h
	DisplacementIndex.cpp
	DisplacementIndex.h
	${PCC_PATH}/PCCCalculator.cpp
	${PCC_PATH}/PCCCalculator.h
	${PCC_PATH}/Combination.cpp
	${PCC_PATH}/Combination.h
)

TARGET_LINK_LIBRARIES(APO_Recognition_sample qhullstatic_r.lib)	
TARGET_LINK_LIBRARIES(VisionErrorCorrection ${OpenCASCADE_ModelingAlgorithms_LIBRARIES} TKSTEP.lib TKSTEP209.lib TKSTEPAttr.lib TKSTEPBase.lib TKShHealing.lib TKXSBase.lib)
TARGET_LINK_LIBRARIES(TaskAnalyzer ${OpenCASCADE_ModelingAlgorithms_LIBRARIES} TKSTEP.lib TKSTEP209.lib TKSTEPAttr.lib TKSTEPBase.lib TKShHealing.lib TKXSBase.lib qhullstatic_r.lib)