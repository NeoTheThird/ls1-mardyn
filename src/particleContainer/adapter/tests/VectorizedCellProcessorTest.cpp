/*
 * VectorizedCellProcessorTest.h
 *
 * @Date: 18.08.2014
 * @Author: tchipevn
 */

#include "VectorizedCellProcessorTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"

TEST_SUITE_REGISTRATION(VectorizedCellProcessorTest);

VectorizedCellProcessorTest::VectorizedCellProcessorTest() {
#if VCP_VEC_TYPE==VCP_NOVEC
	test_log->info() << "VectorizedCellProcessorTest: testing no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	test_log->info() << "VectorizedCellProcessorTest: testing SSE intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	test_log->info() << "VectorizedCellProcessorTest: testing AVX intrinsics." << std::endl;
#endif
}

VectorizedCellProcessorTest::~VectorizedCellProcessorTest() {
}

void VectorizedCellProcessorTest::testForcePotentialCalculationU0() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )

	double forces[4][3] = { { -24, -24, 0 },
	                        {  24, -24, 0 },
	                        { -24,  24, 0 },
	                        {  24,  24, 0 }};

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestU0.inp", 1.1);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	VectorizedCellProcessor cellProcessor( *_domain, 1.1, 1.1);
	container->traverseCells(cellProcessor);

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		m->calcFM();
	}

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), forces[m->id()-1][i], m->F(i), 1e-8);
		}
	}

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(96, _domain->getLocalVirial(), 1e-8);

	delete container;
}

void VectorizedCellProcessorTest::testForcePotentialCalculationF0() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}


	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )
	// -> potential U = -1 per particle pair.
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestF0.inp", 1.3);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	VectorizedCellProcessor cellProcessor( *_domain, 1.3, 1.3);
	container->traverseCells(cellProcessor);

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		m->calcFM();
	}

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << " F[i]=" << m->F(i) << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), 0.0, m->F(i), 1e-7);
		}
	}

	ASSERT_DOUBLES_EQUAL(-4, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-6);

	delete container;
}

void VectorizedCellProcessorTest::testLennardJonesVectorization() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	double ScenarioCutoff = 35.0;
	double Tolerance = 1e-12; // goes through up until 1e-16. Leave it at 1e-12 to be on the safe side
	const char filename[] = {"VectorizationLennardJones.inp"};

	ParticleContainer* container_1 = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, ScenarioCutoff);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), Tolerance);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), Tolerance);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	LegacyCellProcessor cellProcessor( ScenarioCutoff, ScenarioCutoff, ScenarioCutoff, &forceAdapter);
	container_1->traverseCells(cellProcessor);

	for (Molecule* m = container_1->begin(); m != container_1->end(); m = container_1->next()) {
		m->calcFM();
	}

	// store potential and virial for verification
	double legacy_u_pot = _domain ->getLocalUpot();
	double legacy_virial = _domain ->getLocalVirial();


	//-----------------------------------------------------------------------------//
	// Clear any global/static variables between applying the two cell-processors, //
	// but do not delete container_1                                               //
	tearDown(); setUp();
	//-----------------------------------------------------------------------------//


	ParticleContainer* container_2 = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, ScenarioCutoff);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), Tolerance);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), Tolerance);

	VectorizedCellProcessor vectorized_cell_proc(*_domain, ScenarioCutoff, ScenarioCutoff);

	container_2->traverseCells(vectorized_cell_proc);

	for (Molecule* m = container_2->begin(); m != container_2->end(); m = container_2->next()) {
		m->calcFM();
	}

	double vectorized_u_pot = _domain ->getLocalUpot();
	double vectorized_virial = _domain ->getLocalVirial();

	// Traverse both lists simultaneously, advancing both iterators together
	// and assert that the force on the same molecule within both lists is the same:
	Molecule* m_1 = container_1->begin();
	for (Molecule* m_2 = container_2->begin(); m_2 != container_2->end(); m_2 = container_2->next()) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m_2->id() << " index i="<< i << std::endl;
			// check force
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), m_1->F(i), m_2->F(i), Tolerance);
			// check torque
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), m_1->M(i), m_2->M(i), Tolerance);
		}
		// advance molecule of first container
		m_1 = container_1->next();
	}


	// Assert that macroscopic quantities are the same
	ASSERT_DOUBLES_EQUAL(legacy_u_pot,  vectorized_u_pot, Tolerance);
	ASSERT_DOUBLES_EQUAL(legacy_virial, vectorized_virial, Tolerance);

	delete container_1;
	delete container_2;
}

void VectorizedCellProcessorTest::testElectrostaticVectorization(const char* filename, double ScenarioCutoff) {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info()
				<< "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs()
				<< std::endl;
		return;
	}

	// NoVec breaks at 1e-14
	// SSE3 breaks at 1e-13
	// AVX breaks at 1e-14
	// probably architecture dependent. set at 1e-11 to be on the safe side
	// also on other architectures
	double Tolerance = 1e-11;



	ParticleContainer* container_1 = initializeFromFile(
			ParticleContainerFactory::LinkedCell, filename,
			ScenarioCutoff);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), Tolerance);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), Tolerance);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	LegacyCellProcessor cellProcessor(ScenarioCutoff, ScenarioCutoff, ScenarioCutoff, &forceAdapter);
	container_1->traverseCells(cellProcessor);

	for (Molecule* m = container_1->begin(); m != container_1->end(); m = container_1->next()) {
		m->calcFM();
	}

	// store potential and virial for verification
	double legacy_u_pot = _domain->getLocalUpot();
	double legacy_virial = _domain->getLocalVirial();

	//-----------------------------------------------------------------------------//
	// Clear any global/static variables between applying the two cell-processors, //
	// but do not delete container_1                                               //
	tearDown();
	setUp();
	//-----------------------------------------------------------------------------//

	ParticleContainer* container_2 = initializeFromFile(
			ParticleContainerFactory::LinkedCell, filename,
			ScenarioCutoff);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), Tolerance);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), Tolerance);

	VectorizedCellProcessor vectorized_cell_proc(*_domain, ScenarioCutoff,
			ScenarioCutoff);

	container_2->traverseCells(vectorized_cell_proc);

	for (Molecule* m = container_2->begin(); m != container_2->end(); m = container_2->next()) {
		m->calcFM();
	}

	double vectorized_u_pot = _domain->getLocalUpot();
	double vectorized_virial = _domain->getLocalVirial();

	// Traverse both lists simultaneously, advancing both iterators together
	// and assert that the force on the same molecule within both lists is the same:
	Molecule* m_1 = container_1->begin();
	for (Molecule* m_2 = container_2->begin(); m_2 != container_2->end(); m_2 = container_2->next()) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << filename << std::endl;
			str << "Molecule id=" << m_2->id() << " index i=" << i << std::endl;
			// check force
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), m_1->F(i), m_2->F(i), Tolerance);
			// check torque
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), m_1->M(i), m_2->M(i), Tolerance);
		}
		// advance molecule of first container
		m_1 = container_1->next();
	}

	// Assert that macroscopic quantities are the same
	ASSERT_DOUBLES_EQUAL(legacy_u_pot, vectorized_u_pot, Tolerance);
	ASSERT_DOUBLES_EQUAL(legacy_virial, vectorized_virial, Tolerance);

	//-----------------------------------------------------------------------------//
	// Clear any global/static variables between applying the two cell-processors, //
	delete container_1;
	delete container_2;
	tearDown();
	setUp();
	//-----------------------------------------------------------------------------//

}

void VectorizedCellProcessorTest::testChargeChargeVectorization() {
	const char* filename = "VectorizationCharge.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testChargeDipoleVectorization() {
	const char* filename = "VectorizationChargeDipole.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testChargeQuadrupoleVectorization() {
	const char* filename = "VectorizationChargeQuadrupole.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testDipoleDipoleVectorization() {
	const char* filename = "VectorizationDipole.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testDipoleQuadrupoleVectorization() {
	const char* filename = "VectorizationDipoleQuadrupole.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testQuadrupoleQuadrupoleVectorization() {
	const char* filename = "VectorizationQuadrupole.inp";
	testElectrostaticVectorization(filename, 35.0);
}

void VectorizedCellProcessorTest::testWaterVectorization() {
	const char* filename = "VectorizationWater.inp";
	testElectrostaticVectorization(filename, 6.16);
}

