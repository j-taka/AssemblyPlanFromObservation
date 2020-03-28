// ContactCalculator.cpp

#include "ContactCalculator.h"

// const double _NEARLY_ZERO(1.0e-6);
const double _PARALLEL(0.05);

// only thresholding 
void ContactCalculator::Calc(const Shape &moving_object, const Shape &fixed_object)
{
	contact_elements.clear();
	// v-f contacts
	for (size_t i(0); i < moving_object.v_size(); ++i) {
		if (!moving_object.isVConvex(i)) {
			continue;
		}
		for (size_t j(0); j < fixed_object.f_size(); ++j) {
			const std::pair<size_t, size_t> O1(0, i);
			const std::pair<size_t, size_t> O2(1, j);
			const double dist = DistanceFromVertexToFace(i, j, moving_object, fixed_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_VF_CONTACT, dist));
			}
		}
	}
	// e-e contact
	for (size_t i(0); i < moving_object.e_size(); ++i) {
		if (!moving_object.isEConvex(i)) {
			continue;
		}
		for (size_t j(0); j < fixed_object.e_size(); ++j) {
			if (!fixed_object.isEConvex(j)) {
				continue;
			}
			const std::pair<size_t, size_t> O1(0, i);
			const std::pair<size_t, size_t> O2(1, j);
			const double dist = DistanceFromEdgeToEdge(i, j, moving_object, fixed_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_EE_CONTACT, dist));
			}
		}
	}
	// f-v contact
	for (size_t i(0); i < fixed_object.v_size(); ++i) {
		if (!fixed_object.isVConvex(i)) {
			continue;
		}
		for (size_t j(0); j < moving_object.f_size(); ++j) {
			const std::pair<size_t, size_t> O1(1, i);
			const std::pair<size_t, size_t> O2(0, j);
			const double dist = DistanceFromVertexToFace(i, j, fixed_object, moving_object);
			if (fabs(dist) < threshold) {
				contact_elements.push_back(ContactElementForErrorCorrection(O1, O2, ContactElementBase::_VF_CONTACT, dist));
			}
		}
	}
}

void ContactCalculator::DetailedAnalysis(const Shape &moving_object, const Shape &fixed_object)
{
	//
	RemoveImpossibleContact(moving_object, fixed_object);
}

void ContactCalculator::Convert(ContactStateForDisplacement &cs, const Shape &moving_object, const Shape &fixed_object)
{
	SearchSingularContact(cs, moving_object, fixed_object);
}

double ContactCalculator::DistanceFromVertexToFace(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape)
{
	const Eigen::Vector3d diff = v_shape.V(vID) - f_shape.VonF(fID, 0);
	double length = f_shape.OuterNormal(fID).dot(diff);
	// projection onto the face without boundary
	Eigen::Vector3d proj_v = v_shape.V(vID) - f_shape.OuterNormal(fID) * length;
	if (f_shape.GetVFRelation(proj_v, fID) == Shape::_OUTSIDE) {
		f_shape.PointNearFace(proj_v, fID);
		return (proj_v - v_shape.V(vID)).norm();
	}
	else {
		return length;
	}
}

double ContactCalculator::DistanceFromEdgeToEdge(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{
	// MyVector edge1, edge2;
	//double det;
	Eigen::Vector3d edge1 = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
	edge1.normalize();	
	Eigen::Vector3d edge2 = e2_shape.VonE(e2, 1) - e2_shape.VonE(e2, 0);
	edge2.normalize();
	const double det = edge1.dot(edge2); // cos
	if (fabs(det) >= 1 - _PARALLEL) { 
		// paralle 
		return DistanceBetweenParallelEdges(e1, e2, e1_shape, e2_shape);
	}
	else {
		// not parallel
		return DistanceBetweenNonParallelEdges(e1, e2, e1_shape, e2_shape);
	}
}


// distance between two parallel edges
double ContactCalculator::DistanceBetweenParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{	
	double length[4];
	const Eigen::Vector3d e11 = e1_shape.VonE(e1, 0);
	const Eigen::Vector3d e12 = e1_shape.VonE(e1, 1);
	const Eigen::Vector3d e21 = e2_shape.VonE(e2, 0);
	const Eigen::Vector3d e22 = e2_shape.VonE(e2, 1);
	Eigen::Vector3d tmp;
	length[0] = Shape::DistanceFromVertexToEdge(e11, e21, e22, tmp);
	length[1] = Shape::DistanceFromVertexToEdge(e12, e21, e22, tmp);
	length[2] = Shape::DistanceFromVertexToEdge(e21, e11, e12, tmp);
	length[3] = Shape::DistanceFromVertexToEdge(e22, e11, e12, tmp);
	double min_len = length[0];
	// find minimum
	for (int i(1); i < 4; ++i) {
		if (min_len > length[i]) {
			min_len = length[i];
		}
	}
	return min_len;
}

// 
double ContactCalculator::DistanceBetweenNonParallelEdges(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape)
{
	const Eigen::Vector3d e1d = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
	const Eigen::Vector3d e2d = e2_shape.VonE(e2, 1) - e2_shape.VonE(e2, 0);
	const Eigen::Vector3d e12 = e2_shape.VonE(e2, 0) - e1_shape.VonE(e1, 0);
	// solve Equation
	Eigen::Matrix2d A;
	A(0, 0) = e1d.dot(e1d);
	A(1, 0) = e1d.dot(e2d);
	A(0, 1) = -A(1, 0);
	A(1, 1) = -e2d.dot(e2d);
	Eigen::Vector2d b;
	b[0] = e1d.dot(e12);
	b[1] = e2d.dot(e12);
	Eigen::Vector2d x = A.inverse() * b;
	if (x[0] < 0) { x[0] = 0; }
	if (x[0] > 1) { x[0] = 1; }
	if (x[1] < 0) { x[1] = 0; }
	if (x[1] > 1) { x[1] = 1; }
	// calculate distance
	return (e12 + e2d * x[1] - e1d * x[0]).norm();
}

// remove impossible contact elements due to the penetration
void ContactCalculator::RemoveImpossibleContact(const Shape &moving_object, const Shape &fixed_object)
{
	for (int i(0); i < contact_elements.size(); ++i) {
		switch (contact_elements[i].ContactType()) {
		case ContactElementBase::_VF_CONTACT:
			if (contact_elements[i].FirstElement().first == 0) {
				if (isImpossibleContactVF(contact_elements[i].FirstElement().second, contact_elements[i].SecondElement().second, moving_object, fixed_object)) {
					contact_elements.erase(contact_elements.begin() + i);
					--i;
				}
			}
			else {
				if (isImpossibleContactVF(contact_elements[i].FirstElement().second, contact_elements[i].SecondElement().second, fixed_object, moving_object)) {
					contact_elements.erase(contact_elements.begin() + i);
					--i;
				}
			}
			break;
		case ContactElementBase::_EE_CONTACT:
			assert(contact_elements[i].FirstElement().first == 0);
			if (isImpossibleContactEE(contact_elements[i].FirstElement().second, contact_elements[i].SecondElement().second, moving_object, fixed_object) ||
				isImpossibleContactEE(contact_elements[i].SecondElement().second, contact_elements[i].FirstElement().second, fixed_object, moving_object)) {
				contact_elements.erase(contact_elements.begin() + i);
				--i;
			}
			break;
		default:
			std::cerr << "Include some bugs" << std::endl;
			exit(-1);
		}
	}
}

static double proper_trifunc(double src)
{
	if (src < -1) {
		return -1;
	}
	if (src > 1) {
		return 1;
	}
	return src;
}

bool ContactCalculator::isImpossibleContactVF(size_t vID, size_t fID, const Shape &v_shape, const Shape &f_shape, const double ang_th) const
{
	for (size_t i(0); i < v_shape.AdjacentEdges(vID).size(); ++i) {
		const size_t eID = v_shape.AdjacentEdges(vID)[i];
		const Eigen::Vector3d temp = (v_shape.vIDonE(eID, 0) == vID ? v_shape.VonE(eID, 1) - v_shape.VonE(eID, 0) : v_shape.VonE(eID, 0) - v_shape.VonE(eID, 1));
		double prod = f_shape.OuterNormal(fID).dot(temp);
		prod /= temp.norm();
		const double ang = acos(proper_trifunc(prod));
		if (ang > M_PI / 2.0 + ang_th) {
			return true; // penetrate
		}
	}
	return false;
}

// check if an edge penetrates to a face
bool ContactCalculator::isImpossibleContactEE(size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape, const double ang_th) const
{
	const Eigen::Vector3d norm1 = e2_shape.OuterNormal(e2_shape.fIDAdjacentToE(e2, 0));
	const Eigen::Vector3d norm2 = e2_shape.OuterNormal(e2_shape.fIDAdjacentToE(e2, 1));
	const Eigen::Vector3d edge_dic = (e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0)).normalized();
	const double dot1 = norm1.dot(edge_dic);
	const double dot2 = norm2.dot(edge_dic);
	if (fabs(dot1) > ang_th && fabs(dot2) > ang_th && dot1 * dot2 < 0) {
		return false; // not penetrate
	}
	const double ang1 = asin(proper_trifunc(fabs(dot1)));
	const double ang2 = asin(proper_trifunc(fabs(dot2)));
	const Eigen::Vector3d norm3 = e1_shape.OuterNormal(e1_shape.fIDAdjacentToE(e1, 0));
	const Eigen::Vector3d norm4 = e1_shape.OuterNormal(e1_shape.fIDAdjacentToE(e1, 1));
	if (ang1 < ang_th) {
		const double ang3 = asin(proper_trifunc(norm3.dot(norm1)));
		const double ang4 = asin(proper_trifunc(norm4.dot(norm1)));
		if ((ang3 < ang_th) && (ang4 < ang_th)) {
			return false; // not penetrate considering rounding
		}
	}
	if (ang2 < ang_th) {
		const double ang3 = asin(proper_trifunc(norm3.dot(norm2)));
		const double ang4 = asin(proper_trifunc(norm4.dot(norm2)));
		if ((ang3 < ang_th) && (ang4 < ang_th)) {
			return false; // not penetrate considering rounding
		}
	}
	return true;
}

// find singular contact
void ContactCalculator::SearchSingularContact(ContactStateForDisplacement &cs, const Shape &moving_object, const Shape &fixed_object)
{
	// 
	std::vector<ContactInfo> c_info;
	CalculateContactInformation(c_info, moving_object, fixed_object);
	contact_elements.clear(); // update
	for (;!c_info.empty();) {
		std::vector<size_t> cIDs;
		SearchContactWiththeSampContactPosition(cIDs, c_info);
		if (cIDs.size() == 1) {
			AppendContacts(cs, cIDs, c_info, moving_object, fixed_object);
			contact_elements.push_back(c_info[cIDs[0]].elm);
			// remove
			c_info.erase(c_info.begin() + cIDs[0]);
		}
		else {
			// check singular contact
			LocalSingularContact(cs, cIDs, c_info, moving_object, fixed_object);
			// remove
			for (size_t i(0); i < cIDs.size(); ++i) {
				size_t target = cIDs[cIDs.size() - i - 1];
				c_info.erase(c_info.begin() + target);
			}
		}
	}
}

void ContactCalculator::AppendContacts(ContactStateForDisplacement &sc, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const
{ 
	// first append
	for (size_t i(0); i < cIDs.size(); ++i) {
		const ContactInfo *c = &(c_info[cIDs[i]]);
		if (c->parallel) {
			continue;
		}
		switch (c->elm.ContactType()) {
		case ContactElementBase::_VF_CONTACT: {
			const Shape *t_obj = (c->elm.FirstElement().first == 0 ? &fixed_object : &moving_object);
			sc.AddContact(ContactElement::VFContact(c->elm.FirstElement(), c->elm.SecondElement(), c->pos, t_obj->OuterNormal(c->elm.SecondElement().second)));
			break;
		}
		case ContactElementBase::_EE_CONTACT:
			AppendContactEE(sc, *c, moving_object, fixed_object);
			break;
		default:
			std::cerr << "include bugs ... " << std::endl;
			exit(-1);
		}
	}
}

void ContactCalculator::AppendContactEE(ContactStateForDisplacement &sc, const ContactInfo &c_info, const Shape &moving_object, const Shape &fixed_object)
{
	assert(c_info.elm.FirstElement().first == 0);
	const size_t e1 = c_info.elm.FirstElement().second;
	const size_t e2 = c_info.elm.SecondElement().second;
	const Eigen::Vector3d e1_dic = (moving_object.VonE(e1, 1) - moving_object.VonE(e1, 0)).normalized();
	const Eigen::Vector3d e2_dic = (fixed_object.VonE(e2, 1) - fixed_object.VonE(e2, 0)).normalized();
	const Eigen::Vector3d out_dic = fixed_object.OutsideDirectionOfE(e2);
	if ((e1_dic.cross(e2_dic)).dot(out_dic) > 0) {
		sc.AddContact(ContactElement::EEContact(c_info.elm.FirstElement(), c_info.elm.SecondElement(), c_info.pos, e1_dic, e2_dic));
	}
	else {
		sc.AddContact(ContactElement::EEContact(c_info.elm.FirstElement(), c_info.elm.SecondElement(), c_info.pos, e1_dic, -e2_dic));
	}
}

void ContactCalculator::CalculateContactInformation(std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const
{
	c_info.resize(contact_elements.size());
	for (size_t i(0); i < contact_elements.size(); ++i) {
		c_info[i].elm = contact_elements[i];
		switch (contact_elements[i].ContactType()) {
		case ContactElementBase::_VF_CONTACT:
			if (contact_elements[i].FirstElement().first == 0) {
				c_info[i].pos = moving_object.V(contact_elements[i].FirstElement().second);
			}
			else {
				c_info[i].pos = fixed_object.V(contact_elements[i].FirstElement().second);
			}
			c_info[i].parallel = false;
			break;
		case ContactElementBase::_EE_CONTACT:
			assert(contact_elements[i].FirstElement().first == 0);
			c_info[i].parallel = ContactPositionEE(c_info[i].pos, contact_elements[i].FirstElement().second, contact_elements[i].SecondElement().second, moving_object, fixed_object);
			break;
		}
	}
}

bool ContactCalculator::ContactPositionEE(Eigen::Vector3d &pos, size_t e1, size_t e2, const Shape &e1_shape, const Shape &e2_shape, const double th) const
{
	const Eigen::Vector3d e1d = e1_shape.VonE(e1, 1) - e1_shape.VonE(e1, 0);
	const Eigen::Vector3d e2d = e2_shape.VonE(e2, 1) - e2_shape.VonE(e2, 0);
	const Eigen::Vector3d B = e2_shape.VonE(e2, 0) - e1_shape.VonE(e1, 0);
	Eigen::Matrix2d mat;
	mat(0, 0) = e1d.dot(e1d);
	mat(0, 1) = mat(1, 0) = e1d.dot(e2d);
	mat(1, 1) = e2d.dot(e2d);
	if (fabs(mat.determinant()) < th) {
		return true; // parallel
	}
	Eigen::Vector2d BB;
	BB[0] = B.dot(e1d);
	BB[1] = B.dot(e2d);
	Eigen::Vector2d D = mat.inverse() * BB;
	pos = e1_shape.VonE(e1, 0) + e1d * D[0];
	return false; 
}

void ContactCalculator::SearchContactWiththeSampContactPosition(std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info) const
{
	cIDs.push_back(0);
	if (c_info[0].parallel) {
		return; // parallel -> no contact position
	}
	for (size_t i(1); i < c_info.size(); ++i) {
		if (c_info[i].parallel) {
			continue;
		}
		if ((c_info[i].pos - c_info[0].pos).norm() < threshold) {
			cIDs.push_back(i);
		}
	}
}

void ContactCalculator::LocalSingularContact(ContactStateForDisplacement &cs, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object)
{
	// v-v contact
	size_t num;
	if (FindVVContact(cIDs, c_info, moving_object, fixed_object)) {
		ContactStateForDisplacement singular;
		AppendContacts(singular, cIDs, c_info, moving_object, fixed_object);
		cs.AddSingularContact(singular);		
	}
	else if ((num = FindVEContact(cIDs, c_info, moving_object, fixed_object)) != 0) {
		for (size_t i(0); i < num; ++i) {
			ContactStateForDisplacement singular;
			size_t tar = contact_elements.size() - num + i;
			AppendRelatedContacts(singular, contact_elements[tar], cIDs, c_info, moving_object, fixed_object);
			cs.AddSingularContact(singular);
		}
	}
	else {
		// non singular
		AppendContacts(cs, cIDs, c_info, moving_object, fixed_object);
		for (size_t i(0); i < cIDs.size(); ++i) {
			contact_elements.push_back(c_info[cIDs[i]].elm);
		}
	}
}

bool ContactCalculator::FindVVContact(const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object)
{
	// find from v-f andf-v contact elements
	size_t vID_vf(moving_object.v_size()), vID_fv(fixed_object.v_size());
	for (size_t i(0); i < cIDs.size(); ++i) {
		if (c_info[cIDs[i]].elm.ContactType() == ContactElementBase::_VF_CONTACT) {
			if (c_info[cIDs[i]].elm.FirstElement().first == 0) {
				vID_vf = c_info[cIDs[i]].elm.FirstElement().second;
			}
			else {
				vID_fv = c_info[cIDs[i]].elm.FirstElement().second;
			}
		}
	}
	if (vID_vf != moving_object.v_size() && vID_fv != fixed_object.v_size()) {
		const std::pair<size_t, size_t> first(0, vID_vf);
		const std::pair<size_t, size_t> second(1, vID_fv);
		contact_elements.push_back(ContactElementForErrorCorrection(first, second, ContactElementBase::_VV_CONTACT, 0.0));
		return true;
	}
	else {
		return false;
	}
}

int ContactCalculator::FindVEContact(const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object)
{
	std::vector<size_t> vf, fv, ee;
	Division(vf, fv, ee, cIDs, c_info);
	if (ee.empty()) { return false; }
	//
	const int num1 = FindVEContactPart(vf, ee, c_info, moving_object, fixed_object, true);
	const int num2 = FindVEContactPart(fv, ee, c_info, fixed_object, moving_object, false);

	return num1 + num2;
}

int ContactCalculator::FindVEContactPart(const std::vector<size_t> &vf, const std::vector<size_t> &ee, const std::vector<ContactInfo> &c_info, const Shape &v_shape, const Shape &e_shape, bool v_move)
{
	int num(0);
	for (size_t i(0); i < vf.size(); ++i) {
		const size_t vID = c_info[vf[i]].elm.FirstElement().second;
		if (v_shape.isVConvex(vID)) {
			for (size_t j(0); j < ee.size(); ++j) {
				const size_t eID = (v_move ? c_info[ee[j]].elm.SecondElement().second : c_info[ee[j]].elm.FirstElement().second);
				if (e_shape.isEConvex(eID)) {
					if (e_shape.fIDAdjacentToE(eID, 0) == c_info[vf[i]].elm.SecondElement().second ||
						e_shape.fIDAdjacentToE(eID, 1) == c_info[vf[i]].elm.SecondElement().second) {
						size_t objID = (v_move ? 0 : 1);
						ContactElementForErrorCorrection tmp(std::pair<size_t, size_t>(objID, vID), std::pair<size_t, size_t>(1 - objID, eID), ContactElementBase::_VE_CONTACT, 0.0);
						if (std::find(contact_elements.begin(), contact_elements.end(), tmp) == contact_elements.end()) {
							contact_elements.push_back(tmp);
							num++;
						}
					}
				}
			}
		}
	}
	return num;
}

void ContactCalculator::Division(std::vector<size_t> &vf, std::vector<size_t> &fv, std::vector<size_t> &ee, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info)
{
	for (size_t i(0); i < cIDs.size(); ++i) {
		switch (c_info[cIDs[i]].elm.ContactType()) {
		case ContactElementBase::_VF_CONTACT:
			if (c_info[cIDs[i]].elm.FirstElement().first == 0) {
				vf.push_back(cIDs[i]);
			}
			else {
				fv.push_back(cIDs[i]);
			}
			break;
		case ContactElementBase::_EE_CONTACT:
			assert(c_info[cIDs[i]].elm.FirstElement().first == 0);
			ee.push_back(cIDs[i]);
			break;
		default:
			std::cerr << "include bugs..." << std::endl;
			exit(-1);
		}
	}
}

void ContactCalculator::PrintContact(const ContactState &contacts)
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

void ContactCalculator::AppendRelatedContacts(ContactStateForDisplacement &sc, const ContactElementForErrorCorrection &c_elm, const std::vector<size_t> &cIDs, const std::vector<ContactInfo> &c_info, const Shape &moving_object, const Shape &fixed_object) const
{
	// candiates
	ContactState candidates;
	// std::cout << c_elm << std::endl;
	// v-f
	const Shape *f_obj = (c_elm.FirstElement().first == 0 ? &fixed_object : &moving_object);
	for (int i(0); i < 2; ++i) {
		const std::pair<size_t, size_t> s(c_elm.SecondElement().first, f_obj->fIDAdjacentToE(c_elm.SecondElement().second, i));
		candidates.push_back(ContactElementForErrorCorrection(c_elm.FirstElement(), s, ContactElementBase::_VF_CONTACT, 0.0));
	}
	// e-e
	const size_t vID = c_elm.FirstElement().second;
	if (c_elm.FirstElement().first == 0) {
		for (size_t i(0); i < moving_object.AdjacentEdges(vID).size(); ++i) {
			const std::pair<size_t, size_t> f(0, moving_object.AdjacentEdges(vID)[i]);
			candidates.push_back(ContactElementForErrorCorrection(f, c_elm.SecondElement(), ContactElementBase::_EE_CONTACT, 0.0));
		}
	}
	else {
		for (size_t i(0); i < fixed_object.AdjacentEdges(vID).size(); ++i) {
			const std::pair<size_t, size_t> s(1, fixed_object.AdjacentEdges(vID)[i]);
			candidates.push_back(ContactElementForErrorCorrection(c_elm.SecondElement(), s, ContactElementBase::_EE_CONTACT, 0.0));
		}
	}
#if 0
	std::cout << "candidate" << std::endl;
	for (size_t i(0); i < candidates.size(); ++i) {
		std::cout << candidates[i] << std::endl;
	}
#endif
	// find
	for (size_t i(0); i < cIDs.size(); ++i) {
		const ContactInfo *c = &(c_info[cIDs[i]]);
		// std::cout << c->elm;
		if (std::find(candidates.begin(), candidates.end(), c->elm) != candidates.end()) {
			// std::cout << " found" << std::endl;
			if (c->elm.ContactType() == ContactElementBase::_VF_CONTACT) {
				Eigen::Vector3d out_norm = (c->elm.FirstElement().first == 0 ? fixed_object : moving_object).OuterNormal(c->elm.SecondElement().second);
				sc.AddContact(ContactElement::VFContact(c->elm.FirstElement(), c->elm.SecondElement(), c->pos, out_norm));
			}
			else if (c->elm.ContactType() == ContactElementBase::_EE_CONTACT) {
				AppendContactEE(sc, *c, moving_object, fixed_object);
			}
		}
#if 0
		else {
			std::cout << " not found" << std::endl;
		}
#endif
	}
}
