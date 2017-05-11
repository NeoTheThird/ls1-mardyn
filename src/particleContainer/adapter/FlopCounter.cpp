/**
 * \file
 * \brief A CellProcessor that produces Flop information.
 * \author Johannes Heckl
 */

#include "FlopCounter.h"
#include "WrapOpenMP.h"

#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

FlopCounter::_Counts::_Counts():
	_moleculeDistances(0.), _distanceMultiplier(8){
	// 3 sub + 3 square + 2 add


//inverse R squared is one, because only 1/(R^2) has to be calculated, while R^2 already is calculated.

	// Kernel: 15 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale) + 3 (virial tensor)
	// Macro: 4 = 2 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 0 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_LJ, "Lennard-Jones", 15, 4, 12, 2);

	// Kernel: 10 = 1 (inverse R squared) + 1 (square root) + 2 (compute scale) + 3 (apply scale) + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 0 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE, "Charge", 10, 2, 12, 2);

	// Kernel: 34 = 1 (inverse R squared) + 1 (square root) + 29 + 3 (virial tensor)
	// Macro: 3 = 1 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 3 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE_DIPOLE, "Charge-Dipole", 34, 3, 15, 2);

	// Kernel: 101 = 1 (inverse R squared) + 1 (square root) + 96 + 3 (virial tensor)
	// Macro: 5 = 3 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 1 (RF)
	initPotCounter(I_DIPOLE, "Dipole", 101, 5, 18, 3);

	// Kernel: 52 = 1 (inverse R squared) + 1 (square root) + 47 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 3 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_CHARGE_QUADRUPOLE, "Charge-Quadrupole", 52, 2, 15, 2);

	// Kernel: 121 = 1 (inverse R squared) + 1 (square root) + 116 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_DIPOLE_QUADRUPOLE, "Dipole-Quadrupole", 121, 2, 18, 2);

	// Kernel: 131 = 1 (inverse R squared) + 1 (square root) + 126 + 3 (virial tensor)
	// Macro: 2 = 0 (upot) + 2 (virial)
	// sum Forces, Virials and Torques: 6 (forces) + 6 (virials) + 6 (torques)
	// sum Macro: 2 (upot + virial) + 0 (RF)
	initPotCounter(I_QUADRUPOLE, "Quadrupole", 131, 2, 18, 2);
}

void FlopCounter::_PotentialCounts::collCommAppend() {
	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommAppendDouble(_numKernelCalls);
	domainDecomp.collCommAppendDouble(_numMacroCalls);
}

void FlopCounter::_PotentialCounts::collCommGet() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_numKernelCalls = domainDecomp.collCommGetDouble();
	_numMacroCalls = domainDecomp.collCommGetDouble();
}

void FlopCounter::_Counts::allReduce() {
	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommInit(15);

	domainDecomp.collCommAppendDouble(_moleculeDistances);

	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		_potCounts[i].collCommAppend();//adds 2 values each
	}

	domainDecomp.collCommAllreduceSum();
	_moleculeDistances = domainDecomp.collCommGetDouble();

	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		_potCounts[i].collCommGet();
	}

	domainDecomp.collCommFinalize();
}


double FlopCounter::_Counts::process() const {
	using std::endl;
	using Log::global_log;

//	global_log->info() << " Molecule distance: " << getMoleculeDistanceFlops() << endl;
//	global_log->info() << " Center distance: " << getCenterDistanceFlops() << endl;
	for (int i = 0; i < NUM_POTENTIALS; ++i) {
		std::string str = _potCounts[i].printNameKernelAndMacroFlops();
//		if(str.length() > 0) global_log->info() << str;
	}
//	global_log->info() << " Force/Torque Sum: " << getForceTorqueSumFlops() << endl;
//	global_log->info() << " Macroscopic value sum: " << getMacroValueSumFlops() << endl;
//	global_log->info() << "Total: " << getTotalFlops() << endl;

	return getTotalFlops();
}

FlopCounter::FlopCounter(double cutoffRadius, double LJCutoffRadius) : CellProcessor(cutoffRadius, LJCutoffRadius),
		_currentCounts(), _totalFlopCount(0.), _myFlopCount(0.), _synchronized(true){

	const int numThreads = mardyn_get_max_threads();

	global_log->info() << "FlopCounter: allocate data for " << numThreads << " threads." << std::endl;

	_threadData.resize(numThreads);

	#if defined (_OPENMP)
	#pragma omp parallel
	#endif
	{
		_Counts * myown = new _Counts();
		const int myid = mardyn_get_thread_num();
		_threadData[myid] = myown;
	}
}

void FlopCounter::initTraversal() {
	#if defined (_OPENMP)
	#pragma omp master
	#endif
	{
		_currentCounts.clear();
	}
}

void FlopCounter::endTraversal() {
	_Counts sum_omp_values;
	sum_omp_values.clear();

	#if defined (_OPENMP)
	#pragma omp declare reduction(myadd : _Counts : omp_out.addCounts(omp_in))
	#pragma omp parallel reduction(myadd : sum_omp_values)
	#endif
	{
		const int myid = mardyn_get_thread_num();
		_Counts & myown = *_threadData[myid];
		sum_omp_values.addCounts(myown);
		myown.clear();
	}

	_currentCounts = sum_omp_values;

	_myFlopCount = _currentCounts.getTotalFlops();
	if(_synchronized){
		_currentCounts.allReduce();
	}

//	Log::global_log->info() << "FLOP counts in force calculation for this iteration:" << std::endl;
	_totalFlopCount = _currentCounts.process();

//	_totalCounts.addCounts(_currentCounts);
//	Log::global_log->info()
//			<< "Accumulated FLOP counts in force calculation:" << std::endl;
//	_totalFlopCount = _totalCounts.process();
}

class CellPairPolicy_ {
public:

	inline static size_t InitJ(const size_t /*i*/) {
		return 0;
	}
};
class SingleCellPolicy_ {
public:

	inline static size_t InitJ(const size_t i) {
		return i+1;
	}
};

void FlopCounter::processCell(ParticleCell & c) {
#ifndef MARDYN_WR
	FullParticleCell & full_c = downcastReferenceFull(c);
	CellDataSoA& soa = full_c.getCellDataSoA();
#else
	ParticleCell_WR & wr_c = downcastReferenceWR(c);
	CellDataSoA_WR& soa = wr_c.getCellDataSoA();
#endif

	if (c.isHaloCell() or soa._mol_num < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	_calculatePairs<SingleCellPolicy_, CalculateMacroscopic>(soa, soa);
}

void FlopCounter::processCellPair(ParticleCell & c1, ParticleCell & c2) {
	mardyn_assert(&c1 != &c2);
#ifndef MARDYN_WR
	FullParticleCell & full_c1 = downcastReferenceFull(c1);
	FullParticleCell & full_c2 = downcastReferenceFull(c2);
	const CellDataSoA& soa1 = full_c1.getCellDataSoA();
	const CellDataSoA& soa2 = full_c2.getCellDataSoA();
#else
	ParticleCell_WR & wr_c1 = downcastReferenceWR(c1);
	ParticleCell_WR & wr_c2 = downcastReferenceWR(c2);
	const CellDataSoA_WR& soa1 = wr_c1.getCellDataSoA();
	const CellDataSoA_WR& soa2 = wr_c2.getCellDataSoA();
#endif


	const bool c1Halo = c1.isHaloCell();
	const bool c2Halo = c2.isHaloCell();

	// this variable determines whether
	// _calcPairs(soa1, soa2) or _calcPairs(soa2, soa1)
	// is more efficient
	const bool calc_soa1_soa2 = (soa1._mol_num <= soa2._mol_num);

	// if one cell is empty, or both cells are Halo, skip
	if (soa1._mol_num == 0 or soa2._mol_num == 0 or (c1Halo and c2Halo)) {
		return;
	}

	// Macroscopic conditions:
	// if none of the cells is halo, then compute
	// if one of them is halo:
	// 		if full_c1-index < full_c2-index, then compute
	// 		else, then don't compute
	// This saves the Molecule::isLessThan checks
	// and works similar to the "Half-Shell" scheme

	const bool ApplyCutoff = true;

	if ((not c1Halo and not c2Halo) or						// no cell is halo or
			(c1.getCellIndex() < c2.getCellIndex())) 		// one of them is halo, but full_c1.index < full_c2.index
	{
		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic>(soa2, soa1);
		}

	} else {
		mardyn_assert(c1Halo != c2Halo);							// one of them is halo and
		mardyn_assert(not (c1.getCellIndex() < c2.getCellIndex()));// full_c1.index not < full_c2.index

		const bool CalculateMacroscopic = false;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_, CalculateMacroscopic>(soa2, soa1);
		}
	}
}



template<class ForcePolicy, bool CalculateMacroscopic>
void FlopCounter::_calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_pos.zBegin();
	const int * const soa1_mol_ljc_num = soa1._mol_ljc_num;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;
	const int * const soa1_mol_dipoles_num = soa1._mol_dipoles_num;
	const int * const soa1_mol_quadrupoles_num = soa1._mol_quadrupoles_num;

	const vcp_real_calc * const soa2_mol_pos_x = soa2._mol_pos.xBegin();
	const vcp_real_calc * const soa2_mol_pos_y = soa2._mol_pos.yBegin();
	const vcp_real_calc * const soa2_mol_pos_z = soa2._mol_pos.zBegin();
	const int * const soa2_mol_ljc_num = soa2._mol_ljc_num;
	const int * const soa2_mol_charges_num = soa2._mol_charges_num;
	const int * const soa2_mol_dipoles_num = soa2._mol_dipoles_num;
	const int * const soa2_mol_quadrupoles_num = soa2._mol_quadrupoles_num;

	const size_t end_i = soa1._mol_num;
	const size_t end_j = soa2._mol_num;

	unsigned long int i_lj = 0, i_charge=0, i_charge_dipole=0, i_dipole=0, i_charge_quadrupole=0, i_dipole_quadrupole=0, i_quadrupole=0;
	unsigned long int i_mm = 0;

	for (size_t i = 0 ; i < end_i ; ++i) {
		const double m1_x = soa1_mol_pos_x[i];
		const double m1_y = soa1_mol_pos_y[i];
		const double m1_z = soa1_mol_pos_z[i];

		const int numLJcenters_i 	= soa1_mol_ljc_num[i];
		const int numCharges_i 		= soa1_mol_charges_num[i];
		const int numDipoles_i 		= soa1_mol_dipoles_num[i];
		const int numQuadrupoles_i 	= soa1_mol_quadrupoles_num[i];

		#if defined(_OPENMP)
		#pragma omp simd reduction(+ : i_lj, i_charge, i_charge_dipole, i_dipole, i_charge_quadrupole, i_dipole_quadrupole, i_quadrupole, i_mm) \
		aligned(soa1_mol_pos_x, soa1_mol_pos_y, soa1_mol_pos_z, soa1_mol_ljc_num, soa1_mol_charges_num, soa1_mol_dipoles_num, soa1_mol_quadrupoles_num, \
				soa2_mol_pos_x, soa2_mol_pos_y, soa2_mol_pos_z, soa2_mol_ljc_num, soa2_mol_charges_num, soa2_mol_dipoles_num, soa2_mol_quadrupoles_num: 64)
		#endif
		for (size_t j = ForcePolicy::InitJ(i); j < end_j ; ++j) {
			const double m2_x = soa2_mol_pos_x[j];
			const double m2_y = soa2_mol_pos_y[j];
			const double m2_z = soa2_mol_pos_z[j];

			const int numLJcenters_j 	= soa2_mol_ljc_num[j];
			const int numCharges_j 		= soa2_mol_charges_num[j];
			const int numDipoles_j 		= soa2_mol_dipoles_num[j];
			const int numQuadrupoles_j 	= soa2_mol_quadrupoles_num[j];

			const double d_x = m1_x - m2_x;
			const double d_y = m1_y - m2_y;
			const double d_z = m1_z - m2_z;
			const double dist2 = d_x * d_x + d_y * d_y + d_z * d_z;

			++ i_mm;

			if (dist2 < _LJCutoffRadiusSquare) {
				i_lj += numLJcenters_i * numLJcenters_j;
			}

			if(dist2 < _cutoffRadiusSquare) {
				i_charge += numCharges_i * numCharges_j;
				i_charge_dipole += numCharges_i * numDipoles_j + numDipoles_i * numCharges_j;
				i_dipole += numDipoles_i * numDipoles_j;
				i_charge_quadrupole += numCharges_i * numQuadrupoles_j + numQuadrupoles_i * numCharges_j;
				i_dipole_quadrupole += numDipoles_i * numQuadrupoles_j + numQuadrupoles_i * numDipoles_j;
				i_quadrupole += numQuadrupoles_i * numQuadrupoles_j;

			}
		}
	}

	const int tid = mardyn_get_thread_num();
	_Counts &my_threadData = *_threadData[tid];

	my_threadData._moleculeDistances += i_mm;
	my_threadData.addKernelAndMacro(I_LJ, i_lj, CalculateMacroscopic);

	my_threadData.addKernelAndMacro(I_CHARGE, i_charge, CalculateMacroscopic);

	my_threadData.addKernelAndMacro(I_CHARGE_DIPOLE, i_charge_dipole, CalculateMacroscopic);
	my_threadData.addKernelAndMacro(I_DIPOLE, i_dipole, CalculateMacroscopic);

	my_threadData.addKernelAndMacro(I_CHARGE_QUADRUPOLE, i_charge_quadrupole, CalculateMacroscopic);
	my_threadData.addKernelAndMacro(I_DIPOLE_QUADRUPOLE, i_dipole_quadrupole, CalculateMacroscopic);
	my_threadData.addKernelAndMacro(I_QUADRUPOLE, i_quadrupole, CalculateMacroscopic);
}

template<class ForcePolicy, bool CalculateMacroscopic>
void FlopCounter::_calculatePairs(const CellDataSoA_WR & soa1, const CellDataSoA_WR & soa2) {
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_r.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_r.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_r.zBegin();

	const vcp_real_calc * const soa2_mol_pos_x = soa2._mol_r.xBegin();
	const vcp_real_calc * const soa2_mol_pos_y = soa2._mol_r.yBegin();
	const vcp_real_calc * const soa2_mol_pos_z = soa2._mol_r.zBegin();

	const size_t end_i = soa1._mol_num;
	const size_t end_j = soa2._mol_num;

	unsigned long int i_lj = 0;
	unsigned long int i_mm = 0;

	for (size_t i = 0 ; i < end_i ; ++i) {
		const double m1_x = soa1_mol_pos_x[i];
		const double m1_y = soa1_mol_pos_y[i];
		const double m1_z = soa1_mol_pos_z[i];

		const int numLJcenters_i 	= 1;

		#if defined(_OPENMP)
		#pragma omp simd reduction(+ : i_lj, i_mm) \
		aligned(soa1_mol_pos_x, soa1_mol_pos_y, soa1_mol_pos_z, \
				soa2_mol_pos_x, soa2_mol_pos_y, soa2_mol_pos_z: 64)
		#endif
		for (size_t j = ForcePolicy::InitJ(i); j < end_j ; ++j) {
			const double m2_x = soa2_mol_pos_x[j];
			const double m2_y = soa2_mol_pos_y[j];
			const double m2_z = soa2_mol_pos_z[j];

			const int numLJcenters_j 	= 1;

			const double d_x = m1_x - m2_x;
			const double d_y = m1_y - m2_y;
			const double d_z = m1_z - m2_z;
			const double dist2 = d_x * d_x + d_y * d_y + d_z * d_z;

			++ i_mm;

			if (dist2 < _LJCutoffRadiusSquare) {
				i_lj += numLJcenters_i * numLJcenters_j;
			}
		}
	}

	const int tid = mardyn_get_thread_num();
	_Counts &my_threadData = *_threadData[tid];

	my_threadData._moleculeDistances += i_mm;
	my_threadData.addKernelAndMacro(I_LJ, i_lj, CalculateMacroscopic);
}
