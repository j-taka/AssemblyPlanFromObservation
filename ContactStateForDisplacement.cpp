// ContactStateForDisplacement.cpp

#include "ContactStateForDisplacement.h"

std::ostream& operator<<(std::ostream &os, const ContactStateForDisplacement &src)
{
	if (src.elements.empty() && src.singular_elements.empty()) {
		os << "No contact" << std::endl;
	}
	else {
		if (!src.elements.empty()) {
			os << "Non singular:" << std::endl;
			for (size_t i(0); i < src.elements.size(); ++i) {
				os << "Element " << i + 1 << std::endl << src.elements[i] << std::endl;
			}
		}
		if (!src.singular_elements.empty()) {
			for (size_t i(0); i < src.singular_elements.size(); ++i) {
				os << "Singular " << i + 1 << std::endl;
				for (size_t j(0); j < src.singular_elements[i].size(); ++j) {
					os << "Element " << j + 1 << std::endl << src.singular_elements[i][j] << std::endl;
				}
			}
		}
	}
	return os;
}