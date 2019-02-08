// ContactState.cpp

#include "ContactState.h"

std::ostream& operator<<(std::ostream &os, const ContactState &src)
{
	if (src.elements.empty()) {
		os << "No contact" << std::endl;
	}
	else {
		for (size_t i(0); i < src.elements.size(); ++i) {
			os << "Element " << i + 1 << std::endl << src.elements[i] << std::endl;
		}
	}
	return os;
}