// ContactStateForDisplacement.h
#pragma once

#include "ContactElement.h"
#include <cassert>

class ContactStateForDisplacement
{
	friend class InfinitesimulDisplacement;
private:
	std::vector<ContactElement> elements;
	std::vector<std::vector<ContactElement> > singular_elements;
public:
	// constructor
	ContactStateForDisplacement() {}
	// copy constructor
	ContactStateForDisplacement(const ContactStateForDisplacement &src) {
		elements = src.elements;
		singular_elements = src.singular_elements;
	}
	void AddContact(const ContactElement &src) {
		elements.push_back(src);
	}
	void AddSingularContact(const ContactStateForDisplacement &src) {
		singular_elements.push_back(src.elements);
	}
	void Clear() {
		elements.clear();
		singular_elements.clear();
	}
#if 0
	size_t size() const {
		return elements.size();
	}
	const ContactElement& operator[](size_t src) const {
		assert(src < elements.size());
		return elements[src];
	}
#endif
	friend std::ostream& operator<<(std::ostream &os, const ContactStateForDisplacement &src);
};