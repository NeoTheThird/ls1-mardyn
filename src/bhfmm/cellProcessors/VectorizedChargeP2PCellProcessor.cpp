/**
 * \file
 * \brief VectorizedChargeP2PCellProcessor.cpp
 */

#include "particleContainer/adapter/CellDataSoA.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"
#include <algorithm>
#include "particleContainer/adapter/vectorization/MaskGatherChooser.h"
#include "VectorizedChargeP2PCellProcessor.h"

using namespace Log;
using namespace std;
namespace bhfmm {

VectorizedChargeP2PCellProcessor::VectorizedChargeP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		_domain(domain),
		// maybe move the following to somewhere else:
		_upotXpoles(0.0), _virial(0.0){

#if VCP_VEC_TYPE==VCP_NOVEC
	global_log->info() << "VectorizedChargeP2PCellProcessor: using no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	global_log->info() << "VectorizedChargeP2PCellProcessor: using SSE3 intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	global_log->info() << "VectorizedChargeP2PCellProcessor: using AVX2 intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_KNC) || (VCP_VEC_TYPE==VCP_VEC_KNC_GATHER)
	global_log->info() << "VectorizedChargeP2PCellProcessor: using KNC intrinsics." << std::endl;
#elif (VCP_VEC_TYPE==VCP_VEC_KNL) || (VCP_VEC_TYPE==VCP_VEC_KNL_GATHER)
	global_log->info() << "VectorizedChargeP2PCellProcessor: using KNL intrinsics." << std::endl;
#endif

	// initialize thread data
	_numThreads = mardyn_get_max_threads();
	global_log->info() << "VectorizedChargeP2PCellProcessor: allocate data for " << _numThreads << " threads." << std::endl;
	_threadData.resize(_numThreads);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		VCP2PCPThreadData * myown = new VCP2PCPThreadData();
		const int myid = mardyn_get_thread_num();
		_threadData[myid] = myown;
	} // end pragma omp parallel

#ifdef ENABLE_MPI
	_timer.set_sync(false);
#endif
}

VectorizedChargeP2PCellProcessor :: ~VectorizedChargeP2PCellProcessor () {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = mardyn_get_thread_num();
		delete _threadData[myid];
	}
}

void VectorizedChargeP2PCellProcessor::printTimers() {
	std::cout << "FMM: Time spent in Charge P2P " << _timer.get_etime() << std::endl;
}


void VectorizedChargeP2PCellProcessor::initTraversal() {
	_timer.start();

	#if defined(_OPENMP)
	#pragma omp master
	#endif
	{
		_upotXpoles = 0.0;
		_virial = 0.0;
	} // end pragma omp master

}


void VectorizedChargeP2PCellProcessor::endTraversal() {
	double currentVirial = _domain.getLocalVirial();
	double currentUpot = _domain.getLocalUpot();
	double glob_upotXpoles = 0.0;
	double glob_virial = 0.0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:glob_upotXpoles, glob_virial)
	#endif
	{
		const int tid = mardyn_get_thread_num();

		// reduce vectors and clear local variable
		vcp_real_calc thread_upotXpoles = 0.0, thread_virial = 0.0;

		load_hSum_Store_Clear(&thread_upotXpoles, _threadData[tid]->_upotXpolesV);
		load_hSum_Store_Clear(&thread_virial, _threadData[tid]->_virialV);

		// add to global sum
		glob_upotXpoles += thread_upotXpoles;
		glob_virial += thread_virial;
	} // end pragma omp parallel reduction

	_upotXpoles = glob_upotXpoles;
	_virial = glob_virial;
	_domain.setLocalVirial(currentVirial + _virial);
	_domain.setLocalUpot(currentUpot + _upotXpoles);
	_timer.stop();
}


void VectorizedChargeP2PCellProcessor::preprocessCell(ParticleCellPointers & c) {
	// as pre new integration of Caches in SoAs, 
	// this function work as before, as it builds secondary SoAs

	// Determine the total number of centers.
	size_t numMolecules = c.getMoleculeCount();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;
	
	for (auto m = c.moleculesBegin();  m != c.moleculesEnd(); ++m) {
		nCharges += (*m)->numCharges();
	}

	// Construct the SoA.
	CellDataSoA & soa = c.getCellDataSoA();
	soa.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	ComponentList components = *(_simulation.getEnsemble()->getComponents());

	vcp_real_calc* const soa_charges_m_r_x = soa.charges_m_r_xBegin();
	vcp_real_calc* const soa_charges_m_r_y = soa.charges_m_r_yBegin();
	vcp_real_calc* const soa_charges_m_r_z = soa.charges_m_r_zBegin();
	vcp_real_calc* const soa_charges_r_x = soa.charges_r_xBegin();
	vcp_real_calc* const soa_charges_r_y = soa.charges_r_yBegin();
	vcp_real_calc* const soa_charges_r_z = soa.charges_r_zBegin();
	vcp_real_calc* const soa_charges_f_x = soa.charges_f_xBegin();
	vcp_real_calc* const soa_charges_f_y = soa.charges_f_yBegin();
	vcp_real_calc* const soa_charges_f_z = soa.charges_f_zBegin();
	vcp_real_calc* const soa_charges_V_x = soa.charges_V_xBegin();
	vcp_real_calc* const soa_charges_V_y = soa.charges_V_yBegin();
	vcp_real_calc* const soa_charges_V_z = soa.charges_V_zBegin();

	size_t iCharges = 0;
	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < numMolecules; ++i) {
		Molecule& mol = c.moleculesAt(i);

		const size_t mol_charges_num = mol.numCharges();
		const vcp_real_calc mol_pos_x = mol.r(0);
		const vcp_real_calc mol_pos_y = mol.r(1);
		const vcp_real_calc mol_pos_z = mol.r(2);

		soa._mol_pos.x(i) = mol_pos_x;
		soa._mol_pos.y(i) = mol_pos_y;
		soa._mol_pos.z(i) = mol_pos_z;
		soa._mol_charges_num[i] = mol_charges_num;

		for (size_t j = 0; j < mol_charges_num; ++j, ++iCharges)
		{
			soa_charges_m_r_x[iCharges] = mol_pos_x;
			soa_charges_m_r_y[iCharges] = mol_pos_y;
			soa_charges_m_r_z[iCharges] = mol_pos_z;
			soa_charges_r_x[iCharges] = mol.charge_d(j)[0] + mol_pos_x;
			soa_charges_r_y[iCharges] = mol.charge_d(j)[1] + mol_pos_y;
			soa_charges_r_z[iCharges] = mol.charge_d(j)[2] + mol_pos_z;
			soa_charges_f_x[iCharges] = 0.0;
			soa_charges_f_y[iCharges] = 0.0;
			soa_charges_f_z[iCharges] = 0.0;
			soa_charges_V_x[iCharges] = 0.0;
			soa_charges_V_y[iCharges] = 0.0;
			soa_charges_V_z[iCharges] = 0.0;
			//soa._charges_dist_lookup[iCharges] = 0.0;
			// Get the charge
			soa._charges_q[iCharges] = components[mol.componentid()].charge(j).q();
		}
	}
}


void VectorizedChargeP2PCellProcessor::postprocessCell(ParticleCellPointers & c) {
	// as pre new integration of Caches in SoAs, 
	// this function work as before, as it builds secondary SoAs
	using std::isnan; // C++11 required
	CellDataSoA& soa = c.getCellDataSoA();

	vcp_real_calc* const soa_charges_f_x = soa.charges_f_xBegin();
	vcp_real_calc* const soa_charges_f_y = soa.charges_f_yBegin();
	vcp_real_calc* const soa_charges_f_z = soa.charges_f_zBegin();
	vcp_real_calc* const soa_charges_V_x = soa.charges_V_xBegin();
	vcp_real_calc* const soa_charges_V_y = soa.charges_V_yBegin();
	vcp_real_calc* const soa_charges_V_z = soa.charges_V_zBegin();

	// For each molecule iterate over all its centers.
	size_t iCharges = 0;
	for (auto m = c.moleculesBegin(); m != c.moleculesEnd(); ++m) {

		const size_t mol_charges_num = (*m)->numCharges();

		for (size_t i = 0; i < mol_charges_num; ++i, ++iCharges) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = static_cast<double>(soa_charges_f_x[iCharges]);
			f[1] = static_cast<double>(soa_charges_f_y[iCharges]);
			f[2] = static_cast<double>(soa_charges_f_z[iCharges]);
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			(*m)->Fchargeadd(i, f);

			// Store the resulting virial in the molecule.
			double V[3];
			V[0] = static_cast<double>(soa_charges_V_x[iCharges]*0.5);
			V[1] = static_cast<double>(soa_charges_V_y[iCharges]*0.5);
			V[2] = static_cast<double>(soa_charges_V_z[iCharges]*0.5);
			assert(!isnan(V[0]));
			assert(!isnan(V[1]));
			assert(!isnan(V[2]));
			(*m)->Viadd(V);
		}
	}
}



	//const DoubleVec minus_one = DoubleVec::set1(-1.0); //currently not used, would produce warning
	const RealCalcVec zero = RealCalcVec::zero();
	const RealCalcVec one = RealCalcVec::set1(1.0);
	const RealCalcVec two = RealCalcVec::set1(2.0);
	const RealCalcVec three = RealCalcVec::set1(3.0);
	const RealCalcVec four = RealCalcVec::set1(4.0);
	const RealCalcVec five = RealCalcVec::set1(5.0);
	const RealCalcVec six = RealCalcVec::set1(6.0);
	const RealCalcVec ten = RealCalcVec::set1(10.0);
	const RealCalcVec _05 = RealCalcVec::set1(0.5);
	const RealCalcVec _075 = RealCalcVec::set1(0.75);
	const RealCalcVec _1pt5 = RealCalcVec::set1(1.5);
	const RealCalcVec _15 = RealCalcVec::set1(15.0);

	template<bool calculateMacroscopic>
	inline void VectorizedChargeP2PCellProcessor :: _loopBodyCharge(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& qii,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			const RealCalcVec& qjj,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealCalcVec& V_x, RealCalcVec& V_y, RealCalcVec& V_z,
			RealCalcVec& sum_upotXpoles, RealCalcVec& sum_virial,
			const MaskVec& forceMask)
	{
		const RealCalcVec c_dx = r1_x - r2_x;
		const RealCalcVec c_dy = r1_y - r2_y;
		const RealCalcVec c_dz = r1_z - r2_z;//fma not possible since they will be reused...

		const RealCalcVec c_dr2 = RealCalcVec::scal_prod(c_dx, c_dy, c_dz, c_dx, c_dy, c_dz);

		const RealCalcVec c_dr2_inv_unmasked = one / c_dr2;
		const RealCalcVec c_dr2_inv = RealCalcVec::apply_mask(c_dr2_inv_unmasked, forceMask);//masked
	    const RealCalcVec c_dr_inv = RealCalcVec::sqrt(c_dr2_inv);//masked

		const RealCalcVec q1q2per4pie0 = qii * qjj;
		const RealCalcVec upot = q1q2per4pie0 * c_dr_inv;//masked
		const RealCalcVec fac = upot * c_dr2_inv;//masked

		f_x = c_dx * fac;
		f_y = c_dy * fac;
		f_z = c_dz * fac;
		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		V_x = m_dx * f_x;
		V_y = m_dy * f_y;
		V_z = m_dz * f_z;
		// Check if we have to add the macroscopic values up
		if (calculateMacroscopic) {
			sum_upotXpoles = sum_upotXpoles + upot;
			sum_virial = sum_virial + V_x + V_y + V_z;//DoubleVec::scal_prod(m_dx, m_dy, m_dz, f_x, f_y, f_z);
		}
	}

template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
void VectorizedChargeP2PCellProcessor::_calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2) {
	const int tid = mardyn_get_thread_num();
	VCP2PCPThreadData &my_threadData = *_threadData[tid];

	// initialize dist lookups
	soa2.initDistLookupPointersSingle(my_threadData._centers_dist_lookup,
			my_threadData._charges_dist_lookup, soa2._charges_num);

	// Pointer for molecules
	const vcp_real_calc * const soa1_mol_pos_x = soa1._mol_pos.xBegin();
	const vcp_real_calc * const soa1_mol_pos_y = soa1._mol_pos.yBegin();
	const vcp_real_calc * const soa1_mol_pos_z = soa1._mol_pos.zBegin();

	// Pointer for charges
	const vcp_real_calc * const soa1_charges_r_x = soa1.charges_r_xBegin();
	const vcp_real_calc * const soa1_charges_r_y = soa1.charges_r_yBegin();
	const vcp_real_calc * const soa1_charges_r_z = soa1.charges_r_zBegin();
	      vcp_real_calc * const soa1_charges_f_x = soa1.charges_f_xBegin();
	      vcp_real_calc * const soa1_charges_f_y = soa1.charges_f_yBegin();
	      vcp_real_calc * const soa1_charges_f_z = soa1.charges_f_zBegin();
	      vcp_real_calc * const soa1_charges_V_x = soa1.charges_V_xBegin();
	      vcp_real_calc * const soa1_charges_V_y = soa1.charges_V_yBegin();
	      vcp_real_calc * const soa1_charges_V_z = soa1.charges_V_zBegin();
	const vcp_real_calc * const soa1_charges_q = soa1._charges_q;
	const int * const soa1_mol_charges_num = soa1._mol_charges_num;

	const vcp_real_calc * const soa2_charges_m_r_x = soa2.charges_m_r_xBegin();
	const vcp_real_calc * const soa2_charges_m_r_y = soa2.charges_m_r_yBegin();
	const vcp_real_calc * const soa2_charges_m_r_z = soa2.charges_m_r_zBegin();
	const vcp_real_calc * const soa2_charges_r_x   = soa2.charges_r_xBegin();
	const vcp_real_calc * const soa2_charges_r_y   = soa2.charges_r_yBegin();
	const vcp_real_calc * const soa2_charges_r_z   = soa2.charges_r_zBegin();
	      vcp_real_calc * const soa2_charges_f_x   = soa2.charges_f_xBegin();
	      vcp_real_calc * const soa2_charges_f_y   = soa2.charges_f_yBegin();
	      vcp_real_calc * const soa2_charges_f_z   = soa2.charges_f_zBegin();
	      vcp_real_calc * const soa2_charges_V_x   = soa2.charges_V_xBegin();
	      vcp_real_calc * const soa2_charges_V_y   = soa2.charges_V_yBegin();
	      vcp_real_calc * const soa2_charges_V_z   = soa2.charges_V_zBegin();
	const vcp_real_calc * const soa2_charges_q = soa2._charges_q;

	vcp_lookupOrMask_single* const soa2_charges_dist_lookup = my_threadData._charges_dist_lookup;


	RealCalcVec sum_upotXpoles = RealCalcVec::zero();
	RealCalcVec sum_virial = RealCalcVec::zero();

	const RealCalcVec cutoffRadiusSquare = RealCalcVec::set1(_cutoffRadiusSquare);

	/*
	 *  Here different end values for the loops are defined. For loops, which do not vectorize over the last (possibly "uneven") amount of indices, the normal values are computed. These mark the end of the vectorized part.
	 *  The longloop values mark the end of the vectorized part, if vectorization is performed for all elements. For these, various conditions have to be fulfilled, to be sure, that no NaN values are stored and no segfaults arise:
	 *  * arrays have to be long enough
	 *  * arrays have to be filled with something
	 *  * _ljc_id has to be set to existing values for each index
	 *  All of these conditions are fulfilled by setting the non existing values within CellDataSoA.h and AlignedArray.h to zero.
	 */
	const size_t end_charges_j = vcp_floor_to_vec_size(soa2._charges_num);
	const size_t end_charges_j_longloop = vcp_ceil_to_vec_size(soa2._charges_num);//this is ceil _charges_num, VCP_VEC_SIZE
	size_t i_charge_idx = 0;

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._mol_num; ++i) {//over the molecules
		const RealCalcVec m1_r_x = RealCalcVec::broadcast(soa1_mol_pos_x + i);
		const RealCalcVec m1_r_y = RealCalcVec::broadcast(soa1_mol_pos_y + i);
		const RealCalcVec m1_r_z = RealCalcVec::broadcast(soa1_mol_pos_z + i);
		// Iterate over centers of second cell
		const countertype32 compute_molecule_charges = calcDistLookup<ForcePolicy, MaskGatherChooser>(i_charge_idx, soa2._charges_num, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);

		size_t end_charges_loop = MaskGatherChooser::getEndloop(end_charges_j_longloop, compute_molecule_charges);

		// Computation of site interactions with charges

		if (compute_molecule_charges==0) {
			i_charge_idx += soa1_mol_charges_num[i];
		}
		else {
			// Computation of charge-charge interactions

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_charges_num[i]; local_i++) {

				const RealCalcVec q1 = RealCalcVec::broadcast(soa1_charges_q + i_charge_idx + local_i);
				const RealCalcVec r1_x = RealCalcVec::broadcast(soa1_charges_r_x + i_charge_idx + local_i);
				const RealCalcVec r1_y = RealCalcVec::broadcast(soa1_charges_r_y + i_charge_idx + local_i);
				const RealCalcVec r1_z = RealCalcVec::broadcast(soa1_charges_r_z + i_charge_idx + local_i);

				RealCalcVec sum_f1_x = RealCalcVec::zero();
				RealCalcVec sum_f1_y = RealCalcVec::zero();
				RealCalcVec sum_f1_z = RealCalcVec::zero();

				RealCalcVec sum_V1_x = RealCalcVec::zero();
				RealCalcVec sum_V1_y = RealCalcVec::zero();
				RealCalcVec sum_V1_z = RealCalcVec::zero();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ2(i_charge_idx + local_i);

				for (; j < end_charges_loop; j += VCP_VEC_SIZE) {
					const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMask(soa2_charges_dist_lookup, j);
					// Check if we have to calculate anything for at least one of the pairs
					if (MaskGatherChooser::computeLoop(lookupORforceMask)) {
						const RealCalcVec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);

						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealCalcVec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								MaskGatherChooser::getForceMask(lookupORforceMask));



						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask);
						vcp_simd_load_sub_store<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask);

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask);
						vcp_simd_load_add_store<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask);
					}
				}
#if VCP_VEC_TYPE == VCP_VEC_KNC_GATHER or VCP_VEC_TYPE == VCP_VEC_KNL_GATHER
				if(MaskGatherChooser::hasRemainder()){//remainder computations, that's not an if, but a constant branch... compiler is wise.
					const __mmask8 remainderM = MaskGatherChooser::getRemainder(compute_molecule_charges);
					if(remainderM != 0x00){
						const vcp_lookupOrMask_vec lookupORforceMask = MaskGatherChooser::loadLookupOrForceMaskRemainder(soa2_charges_dist_lookup, j, remainderM);

						const RealCalcVec q2 = MaskGatherChooser::load(soa2_charges_q, j, lookupORforceMask);
						const RealCalcVec r2_x = MaskGatherChooser::load(soa2_charges_r_x, j, lookupORforceMask);
						const RealCalcVec r2_y = MaskGatherChooser::load(soa2_charges_r_y, j, lookupORforceMask);
						const RealCalcVec r2_z = MaskGatherChooser::load(soa2_charges_r_z, j, lookupORforceMask);

						const RealCalcVec m2_r_x = MaskGatherChooser::load(soa2_charges_m_r_x, j, lookupORforceMask);
						const RealCalcVec m2_r_y = MaskGatherChooser::load(soa2_charges_m_r_y, j, lookupORforceMask);
						const RealCalcVec m2_r_z = MaskGatherChooser::load(soa2_charges_m_r_z, j, lookupORforceMask);

						RealCalcVec f_x, f_y, f_z;
						RealCalcVec Vx, Vy, Vz;

						_loopBodyCharge<CalculateMacroscopic>(
								m1_r_x, m1_r_y, m1_r_z,	r1_x, r1_y, r1_z, q1,
								m2_r_x, m2_r_y, m2_r_z, r2_x, r2_y, r2_z, q2,
								f_x, f_y, f_z,
								Vx, Vy, Vz,
								sum_upotXpoles, sum_virial,
								remainderM);



						sum_f1_x = sum_f1_x + f_x;
						sum_f1_y = sum_f1_y + f_y;
						sum_f1_z = sum_f1_z + f_z;

						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_x, j, f_x, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_y, j, f_y, lookupORforceMask, remainderM);
						vcp_simd_load_sub_store_masked<MaskGatherChooser>(soa2_charges_f_z, j, f_z, lookupORforceMask, remainderM);

						sum_V1_x = sum_V1_x + Vx;
						sum_V1_y = sum_V1_y + Vy;
						sum_V1_z = sum_V1_z + Vz;

						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_x, j, Vx, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_y, j, Vy, lookupORforceMask, remainderM);
						vcp_simd_load_add_store_masked<MaskGatherChooser>(soa2_charges_V_z, j, Vz, lookupORforceMask, remainderM);
					}
				}
#endif

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx + local_i, sum_f1_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx + local_i, sum_f1_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx + local_i, sum_f1_z);
				// Add old virial and summed calculated virials for center 1
				hSum_Add_Store(soa1_charges_V_x + i_charge_idx + local_i, sum_V1_x);
				hSum_Add_Store(soa1_charges_V_y + i_charge_idx + local_i, sum_V1_y);
				hSum_Add_Store(soa1_charges_V_z + i_charge_idx + local_i, sum_V1_z);

			}
			i_charge_idx += soa1_mol_charges_num[i];
		}
	}

	hSum_Add_Store(my_threadData._upotXpolesV, sum_upotXpoles);
	hSum_Add_Store(my_threadData._virialV, sum_virial);
} // void VectorizedChargeP2PCellProcessor::_calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2)

void VectorizedChargeP2PCellProcessor::processCell(ParticleCellPointers & c) {
	CellDataSoA& soa = c.getCellDataSoA();
	if (c.isHaloCell() or soa._mol_num < 2) {
		return;
	}
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = false;
	_calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa, soa);
}

void VectorizedChargeP2PCellProcessor::processCellPair(ParticleCellPointers & c1, ParticleCellPointers & c2) {
	assert(&c1 != &c2);
	const CellDataSoA& soa1 = c1.getCellDataSoA();
	const CellDataSoA& soa2 = c2.getCellDataSoA();
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
	// 		if c1-index < c2-index, then compute
	// 		else, then don't compute
	// This saves the Molecule::isLessThan checks
	// and works similar to the "Half-Shell" scheme

	const bool ApplyCutoff = false;

	if ((not c1Halo and not c2Halo) or						// no cell is halo or
			(c1.getCellIndex() < c2.getCellIndex())) 		// one of them is halo, but c1.index < c2.index
	{
		const bool CalculateMacroscopic = true;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}

	} else {
		assert(c1Halo != c2Halo);							// one of them is halo and
		assert(not (c1.getCellIndex() < c2.getCellIndex()));// c1.index not < c2.index

		const bool CalculateMacroscopic = false;

		if (calc_soa1_soa2) {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa1, soa2);
		} else {
			_calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(soa2, soa1);
		}
	}
}

} // namespace bhfmm

