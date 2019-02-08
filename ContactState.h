// ContactState.h
#pragma once

#include "ContactElement.h"
#include <cassert>

class ContactState
{
private:
	std::vector<ContactElement> elements;
public:
	// constructor
	ContactState() {}
	// copy constructor
	ContactState(const ContactState &src) {
		elements = src.elements;
	}
	void AddContact(const ContactElement &src) {
		elements.push_back(src);
	}
	void Clear() {
		elements.clear();
	}
	size_t size() const {
		return elements.size();
	}
	const ContactElement& operator[](size_t src) const {
		assert(src < elements.size());
		return elements[src];
	}
	friend std::ostream& operator<<(std::ostream &os, const ContactState &src);
};