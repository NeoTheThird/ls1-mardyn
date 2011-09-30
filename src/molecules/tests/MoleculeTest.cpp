/*
 * MoleculeTest.cpp
 *
 * @Date: 04.02.2011
 * @Author: eckhardw
 */

#include "MoleculeTest.h"
#include "molecules/BasicMolecule.h"
#include "molecules/CachingMolecule.h"
#include "molecules/MoleculeTypes.h"

TEST_SUITE_REGISTRATION(MoleculeTest);

MoleculeTest::MoleculeTest() {

}

MoleculeTest::~MoleculeTest() {
}


void MoleculeTest::testIsLessThan() {
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,0,0,0,0,false);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);
	Molecule b(0, 0, 2.0, 2.0, 2.0,0,0,0,0,0,0,0,0,0,0, &components);

	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(2, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));

	a.setr(2, 2.0);
	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(1, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));

	a.setr(1, 2.0);
	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(0, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));
}

void MoleculeTest::testRotationZero() {
	Component component;
	component.addLJcenter(1, 0, 0, 1, 1, 1, 5, false);
	component.addCharge(0, 2, 0, 1.5, 1.5);
	component.addDipole(0, 0, 3, 0.5, 0.5, 0.5, 0.5);
	component.addQuadrupole(3, 3, 3, 0.7, 0.7, 0.7, 0.7);

	std::vector<Component> components;
	components.push_back(component);

	BasicMolecule::setComponentsVector(&components);

	// this corresponds to a molecule for which the body-fixed coordinate system
	// is aligned with the global one, i.e. rotation is 0 degree.
	BasicMolecule basicMolecule(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &components);
	CachingMolecule cachingMolecule(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &components);
	cachingMolecule.upd_cache();

	double position[3];
	basicMolecule.ljcenter_d(0, position);
	ASSERT_DOUBLES_EQUAL(1.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	cachingMolecule.ljcenter_d(0, position);
	ASSERT_DOUBLES_EQUAL(1.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	basicMolecule.charge_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(2.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	cachingMolecule.charge_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(2.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	basicMolecule.dipole_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	cachingMolecule.dipole_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	basicMolecule.quadrupole_d(0, position);
	ASSERT_DOUBLES_EQUAL(3.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	cachingMolecule.quadrupole_d(0, position);
	ASSERT_DOUBLES_EQUAL(3.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);
}


void MoleculeTest::testRotationPI() {
	Component component;
	component.addLJcenter(1, 0, 0, 1, 1, 1, 5, false);
	component.addCharge(0, 2, 0, 1.5, 1.5);
	component.addDipole(0, 0, 3, 0.5, 0.5, 0.5, 0.5);
	component.addQuadrupole(3, 3, 3, 0.7, 0.7, 0.7, 0.7);

	std::vector<Component> components;
	components.push_back(component);

	BasicMolecule::setComponentsVector(&components);

	// this corresponds to a molecule for which the body-fixed coordinate system
	// is rotated by pi.
	BasicMolecule basicMolecule(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, &components);
	CachingMolecule cachingMolecule(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, &components);
	cachingMolecule.upd_cache();

	double position[3];
	basicMolecule.ljcenter_d(0, position);
	ASSERT_DOUBLES_EQUAL(-1.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	cachingMolecule.ljcenter_d(0, position);
	ASSERT_DOUBLES_EQUAL(-1.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	basicMolecule.charge_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(-2.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	cachingMolecule.charge_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(-2.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[2], 10e-7);

	basicMolecule.dipole_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	cachingMolecule.dipole_d(0, position);
	ASSERT_DOUBLES_EQUAL(0.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(0.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	basicMolecule.quadrupole_d(0, position);
	ASSERT_DOUBLES_EQUAL(-3.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(-3.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);

	cachingMolecule.quadrupole_d(0, position);
	ASSERT_DOUBLES_EQUAL(-3.0, position[0], 10e-7);
	ASSERT_DOUBLES_EQUAL(-3.0, position[1], 10e-7);
	ASSERT_DOUBLES_EQUAL(3.0, position[2], 10e-7);
}

void MoleculeTest::testForceMoment() {
	Component component;
	component.addLJcenter(1, 0, 0, 1, 1, 1, 5, false);
	component.addCharge(0, 2, 0, 1.5, 1.5);
	component.addDipole(0, 0, 3, 0.5, 0.5, 0.5, 0.5);
	component.addQuadrupole(3, 3, 3, 0.7, 0.7, 0.7, 0.7);

	std::vector<Component> components;
	components.push_back(component);

	BasicMolecule::setComponentsVector(&components);

	CachingMolecule cachingMolecule(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &components);
	cachingMolecule.upd_cache();
	BasicMolecule basicMolecule(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &components);

	/*******  LJ add and sub ******/

	double force[] = { 0.3, 0.4, 0.7};
	cachingMolecule.Fljcenteradd(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fljcenteradd(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	cachingMolecule.Fljcentersub(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fljcentersub(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	/*** Charge add and sub ***/
	cachingMolecule.Fchargeadd(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fchargeadd(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	cachingMolecule.Fchargesub(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fchargesub(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	/*** Dipole add and sub ***/
	cachingMolecule.Fdipoleadd(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fdipoleadd(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	cachingMolecule.Fdipolesub(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fdipolesub(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	/*** Quadrupole add and sub ***/
	cachingMolecule.Fquadrupoleadd(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fquadrupoleadd(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();

	cachingMolecule.Fquadrupolesub(0, force);
	cachingMolecule.calcFM();
	basicMolecule.Fquadrupolesub(0, force);

	for (int i = 0; i < 3; i++) {
		std::stringstream stream;
		stream << "force component i=" << i;
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.M(i), basicMolecule.M(i), 10e-7);
		ASSERT_DOUBLES_EQUAL_MSG(stream.str(), cachingMolecule.F(i), basicMolecule.F(i), 10e-7);
	}

	cachingMolecule.clearFM();
	basicMolecule.clearFM();
}
