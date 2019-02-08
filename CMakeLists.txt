PROJECT(AssemblyPlanFromObservation)

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

#eigen
set(CMAKE_MODULE_PATH $ENV{Eigen3_DIR}/cmake)
set(EIGEN3_INCLUDE_DIR $ENV{Eigen3_DIR})
FIND_PACKAGE(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#qhull
set(QHULL_PATH $ENV{QHULL_ROOT})
include_directories(${QHULL_PATH}/include)
link_directories(${QHULL_PATH}/lib)

#PCC
set(PCC_PATH ../PolyhedralConvexCones)
include_directories(${PCC_PATH})

ADD_EXECUTABLE(
	APO_Recognition_sample
	APO_Recognition_sample.cpp
	ContactElement.cpp
	ContactElement.h
	ContactState.cpp
	ContactState.h
	InfinitesimulDisplacement.cpp
	InfinitesimulDisplacement.h
)

TARGET_LINK_LIBRARIES(APO_Recognition_sample)	
