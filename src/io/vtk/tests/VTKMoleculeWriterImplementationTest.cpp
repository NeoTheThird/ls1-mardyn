/*
 * VTKMoleculeWriterImplementationTest.cpp
 *
 * @Date: 25.08.2010
 * @Author: eckhardw
 */

#include "VTKMoleculeWriterImplementationTest.h"
#include "io/vtk/VTKMoleculeWriterImplementation.h"
#include "utils/FileUtils.h"
#include "molecules/MoleculeTypes.h"
#include <vector>
#ifdef PARALLEL
#include <mpi.h>
#endif

TEST_SUITE_REGISTRATION(VTKMoleculeWriterImplementationTest);

VTKMoleculeWriterImplementationTest::VTKMoleculeWriterImplementationTest() {
}

VTKMoleculeWriterImplementationTest::~VTKMoleculeWriterImplementationTest() {
}

void VTKMoleculeWriterImplementationTest::testInitialization() {
	VTKMoleculeWriterImplementation writer(0);

	ASSERT_EQUAL(writer.isVTKFileInitialized(), false);
	writer.initializeVTKFile();
	ASSERT_EQUAL(writer.isVTKFileInitialized(), true);

	ASSERT_EQUAL(writer.isParallelVTKFileInitialized(), false);
	std::vector<std::string> v;
	writer.initializeParallelVTKFile(v);
	ASSERT_EQUAL(writer.isParallelVTKFileInitialized(), true);
}

void VTKMoleculeWriterImplementationTest::testWriteVTKFile() {
#ifdef PARALLEL
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank != 0) {
		return;
	}
#endif
	VTKMoleculeWriterImplementation writer(0);

	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,0,0,0,0,false);
	components.push_back(dummyComponent);
	Molecule dummyMolecule(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &components);

	std::vector<std::string> v;
	writer.initializeVTKFile();
	writer.initializeParallelVTKFile(v);

	writer.plotMolecule(dummyMolecule);
	writer.plotMolecule(dummyMolecule);
	ASSERT_EQUAL(writer.getNumMoleculesPlotted(), 2u);

	writer.writeVTKFile("VTKMoleculeWriter00.vtu");
	ASSERT_TRUE(fileExists("VTKMoleculeWriter00.vtu"));
	writer.writeParallelVTKFile("VTKMoleculeWriter00.pvtu");
	ASSERT_TRUE(fileExists("VTKMoleculeWriter00.pvtu"));

	removeFile("VTKMoleculeWriter00.vtu");
	removeFile("VTKMoleculeWriter00.pvtu");
}
