// ContactElement.cpp

#include "ContactElement.h"

std::ostream& operator<<(std::ostream &ofs, const ContactElementForErrorCorrection &src)
{
	switch (src.type) {
	case ContactElementBase::_VF_CONTACT:
		ofs << "v-f contact:";
		break;
	case ContactElementBase::_EE_CONTACT:
		ofs << "e-e contact: ";
		break;
	case ContactElementBase::_VE_CONTACT:
		ofs << "v-e contact: ";
		break;
	case ContactElementBase::_VV_CONTACT:
		ofs << "v-v contact: ";
		break;
	}
	ofs << "(" << src.ID1.first << "," << src.ID1.second << ") - (" << src.ID2.first << "," << src.ID2.second << ") (distance: " << src.distance << ")";
	return ofs;
}

ContactElement ContactElement::VFContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &f_norm)
{
	ContactElement res;
	res.type = _VF_CONTACT;
	res.ID1 = ID1;
	res.ID2 = ID2;
	res.parameters.resize(6);
	res.parameters[0] = contact_position[0];
	res.parameters[1] = contact_position[1];
	res.parameters[2] = contact_position[2];
	res.parameters[3] = f_norm[0];
	res.parameters[4] = f_norm[1];
	res.parameters[5] = f_norm[2];
	return res;
}

ContactElement ContactElement::EEContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &e1_direction, const Eigen::Vector3d &e2_direction)
{
	if (ID1.first > ID2.first) {
		return EEContact(ID2, ID1, contact_position, e2_direction, e1_direction);
	}
	ContactElement res;
	res.type = _EE_CONTACT;
	res.ID1 = ID1;
	res.ID2 = ID2;
	res.parameters.resize(9);
	res.parameters[0] = contact_position[0];
	res.parameters[1] = contact_position[1];
	res.parameters[2] = contact_position[2];
	res.parameters[3] = e1_direction[0];
	res.parameters[4] = e1_direction[1];
	res.parameters[5] = e1_direction[2];
	res.parameters[6] = e2_direction[0];
	res.parameters[7] = e2_direction[1];
	res.parameters[8] = e2_direction[2];
	return res;
}


// singular contact	element
ContactElement ContactElement::VEContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Eigen::Vector3d &norm1, const Eigen::Vector3d &norm2)
{
	ContactElement res;
	res.type = _VE_CONTACT;
	res.ID1 = ID1;
	res.ID2 = ID2;
	res.parameters.resize(9);
	res.parameters[0] = contact_position[0];
	res.parameters[1] = contact_position[1];
	res.parameters[2] = contact_position[2];
	res.parameters[3] = norm1[0];
	res.parameters[4] = norm1[1];
	res.parameters[5] = norm1[2];
	res.parameters[6] = norm2[0];
	res.parameters[7] = norm2[1];
	res.parameters[8] = norm2[2];
	return res;
}

ContactElement ContactElement::VVContact(const std::pair<size_t, size_t> &ID1, const std::pair<size_t, size_t> &ID2, const Eigen::Vector3d &contact_position, const Normals &norm1, const Normals &norm2)
{
	ContactElement res;
	res.type = _VV_CONTACT;
	res.ID1 = ID1;
	res.ID2 = ID2;
	res.parameters.resize(5 + norm1.size() * 3 + norm2.size() * 3);
	res.parameters[0] = contact_position[0];
	res.parameters[1] = contact_position[1];
	res.parameters[2] = contact_position[2];
	// memorize nubmers
	res.parameters[3] = (double) norm1.size();
	res.parameters[4] = (double) norm2.size();
	size_t offset = 5;
	for (size_t i(0); i < norm1.size(); ++i) {
		res.parameters[offset] = norm1[i][0]; offset++;
		res.parameters[offset] = norm1[i][1]; offset++;
		res.parameters[offset] = norm1[i][2]; offset++;
	}
	for (size_t i(0); i < norm2.size(); ++i) {
		res.parameters[offset] = norm2[i][0]; offset++;
		res.parameters[offset] = norm2[i][1]; offset++;
		res.parameters[offset] = norm2[i][2]; offset++;
	}
	return res;
}

void ContactElement::GetNormalInformation(Eigen::MatrixXd &n1, Eigen::MatrixXd &n2) const
{
	assert(type == _VV_CONTACT);
	n1 = Eigen::MatrixXd((int) parameters[3], 3);
	n2 = Eigen::MatrixXd((int) parameters[4], 3);
	// copy
	size_t offset = 5;
	for (int r(0); r < n1.rows(); ++r) {
		n1(r, 0) = parameters[offset]; offset++;
		n1(r, 1) = parameters[offset]; offset++;
		n1(r, 2) = parameters[offset]; offset++;
	}
	for (int r(0); r < n2.rows(); ++r) {
		n2(r, 0) = parameters[offset]; offset++;
		n2(r, 1) = parameters[offset]; offset++;
		n2(r, 2) = parameters[offset]; offset++;
	}
}

std::ostream& operator<<(std::ostream &os, const ContactElement &src)
{
	switch (src.type) {
	case ContactElement::_VF_CONTACT:
		os << "Vertex " << src.ID1.second + 1 << " of Object " << src.ID1.first + 1 << std::endl;
		os << "Face " << src.ID2.second + 1 << " of Object " << src.ID2.first + 1 << std::endl;
		os << "Contact position: (" << src.parameters[0] << ", " << src.parameters[1] << ", " << src.parameters[2] << ")" << std::endl;
		os << "Face outer normal: (" << src.parameters[3] << ", " << src.parameters[4] << ", " << src.parameters[5] << ")" << std::endl;
		break;
	case ContactElement::_EE_CONTACT:
		os << "Edge " << src.ID1.second + 1 << " of Object " << src.ID1.first + 1 << std::endl;
		os << "Edge " << src.ID2.second + 1 << " of Object " << src.ID2.first + 1 << std::endl;
		os << "Contact position: (" << src.parameters[0] << ", " << src.parameters[1] << ", " << src.parameters[2] << ")" << std::endl;
		os << "Direction of Edge of Object " << src.ID1.first + 1 << ": (" << src.parameters[3] << ", " << src.parameters[4] << ", " << src.parameters[5] << ")" << std::endl;
		os << "Direction of Edge of Object " << src.ID2.first + 1 << ": (" << src.parameters[6] << ", " << src.parameters[7] << ", " << src.parameters[8] << ")" << std::endl;
		break;
	case ContactElement::_VE_CONTACT:
		os << "Vertex " << src.ID1.second + 1 << " of Object " << src.ID1.first + 1 << std::endl;
		os << "Edge " << src.ID2.second + 1 << " of Object " << src.ID2.first + 1 << std::endl;
		os << "Contact position: (" << src.parameters[0] << ", " << src.parameters[1] << ", " << src.parameters[2] << ")" << std::endl;
		os << "Outer normal of one face adjacent to Edge " << src.ID2.second + 1 << ": (" << src.parameters[3] << ", " << src.parameters[4] << ", " << src.parameters[5] << ")" << std::endl;
		os << "Outer normal of the other face adjacent to Edge " << src.ID2.second + 1 << ": (" << src.parameters[6] << ", " << src.parameters[7] << ", " << src.parameters[8] << ")" << std::endl;
		break;
	case ContactElement::_VV_CONTACT:
	{
		os << "Vertex " << src.ID1.second + 1 << " of Object " << src.ID1.first + 1 << std::endl;
		os << "Vertex " << src.ID2.second + 1 << " of Object " << src.ID2.first + 1 << std::endl;
		os << "Contact position: (" << src.parameters[0] << ", " << src.parameters[1] << ", " << src.parameters[2] << ")" << std::endl;
		os << (int) src.parameters[3] << " faces adjacent to vertex of Object " << src.ID1.first + 1 << std::endl;
		size_t offset = 5;
		for (int i(0); i < (int)src.parameters[3]; ++i, offset += 3) {
			os << "(" << src.parameters[offset] << ", " << src.parameters[offset + 1] << ", " << src.parameters[offset + 2] << ")" << std::endl;
		}
		os << (int)src.parameters[4] << " faces adjacent to vertex of Object " << src.ID2.first + 1 << std::endl;
		for (int i(0); i < (int)src.parameters[4]; ++i, offset += 3) {
			os << "(" << src.parameters[offset] << ", " << src.parameters[offset + 1] << ", " << src.parameters[offset + 2] << ")" << std::endl;
		}
		break;
	}
	}
	return os;
}
