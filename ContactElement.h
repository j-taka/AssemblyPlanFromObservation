// ContactElement.h
#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Dense>

class ContactElement
{
public:
	enum Type {
		_VF_CONTACT, 
		_EE_CONTACT, 
		_VE_CONTACT, 
		_VV_CONTACT 
	};
	typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > Normals;
private:
	std::pair<size_t, size_t> ID1;
	std::pair<size_t, size_t> ID2;
	Type type;
	std::vector<double> parameters;
public:
	// constructor (no action)
	ContactElement(){}
	// copy constructor
	ContactElement(const ContactElement &src) {
		ID1 = src.ID1;
		ID2 = src.ID2;
		type = src.type;
		parameters = src.parameters;
	}
	size_t FirstObject() const {
		return ID1.first;
	}
	size_t SecondObject() const {
		return ID2.first;
	}
	Type Type() const {
		return type;
	}
	bool isSingular() const {
		return (type == _VE_CONTACT || type == _VV_CONTACT);
	}
	Eigen::Vector3d ContactPosition() const {
		return Eigen::Vector3d(parameters[0], parameters[1], parameters[2]);
	}
	Eigen::Vector3d OuterNormal() const {
		switch (type) {
		case _VF_CONTACT:
			return Eigen::Vector3d(parameters[3], parameters[4], parameters[5]);
		case _EE_CONTACT:
		{
			const Eigen::Vector3d e1(parameters[3], parameters[4], parameters[5]);
			const Eigen::Vector3d e2(parameters[6], parameters[7], parameters[8]);
			Eigen::Vector3d outer_norm = e1.cross(e2);
			outer_norm.normalize();
			return outer_norm;
		}
		default:
			std::cerr << "Do not use in the type" << std::endl << *this << std::endl;
			exit(-1); // do not use in the other type
		}
	}
	Eigen::Vector3d OuterNormal(int i) const {
		assert(type == _VE_CONTACT);
		return Eigen::Vector3d(parameters[i * 3], parameters[i * 3 + 1], parameters[i * 3 + 2]);
	}
	void GetNormalInformation(Eigen::MatrixXd &n1, Eigen::MatrixXd &n2) const;
	// non-singular contact element
	static ContactElement VFContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &f_norm);
	static ContactElement EEContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &e1_direction, const Eigen::Vector3d &e2_direction);
	// singular contact	element
	static ContactElement VEContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &norm1, const Eigen::Vector3d &norm2);
	static ContactElement VVContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Normals &norm1, const Normals &norm2);
	friend std::ostream& operator<<(std::ostream &os, const ContactElement &src);
};