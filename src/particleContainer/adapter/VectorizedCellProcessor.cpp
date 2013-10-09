/**
 * \file
 * \brief VectorizedCellProcessor.cpp
 */

#include "VectorizedCellProcessor.h"
#include "CellDataSoA.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#include "Domain.h"
#include "utils/Logger.h"

#include <algorithm>

using namespace Log;

VectorizedCellProcessor::VectorizedCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius) :
		_domain(domain), _cutoffRadiusSquare(cutoffRadius * cutoffRadius), _LJcutoffRadiusSquare(LJcutoffRadius * LJcutoffRadius),
		_compIDs(), _eps_sig(), _shift6(), _upot6lj(0.0), _virial(0.0), _centers_dist_lookup(128) {
#if VLJCP_VEC_TYPE==VLJCP_NOVEC
	Log::global_log->info() << "VectorizedLJCellProcessor: no vectorization."
	<< std::endl;
#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
	Log::global_log->info() << "VectorizedLJCellProcessor: using SSE3."
			<< std::endl;
#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX
	Log::global_log->info() << "VectorizedLJCellProcessor: using AVX."
			<< std::endl;
#endif

	ComponentList components = _domain.getComponents();
	// Get the maximum Component ID.
	size_t maxID = 0;
	const ComponentList::const_iterator end = components.end();
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c)
		maxID = std::max(maxID, static_cast<size_t>(c->ID()));

	// Assign a center list start index for each component.
	_compIDs.resize(maxID + 1, 0);
	size_t centers = 0;
	for (ComponentList::const_iterator c = components.begin(); c != end; ++c) {
		_compIDs[c->ID()] = centers;
		centers += c->numLJcenters();
	}

	// One row for each LJ Center, one pair (epsilon*24, sigma^2) for each LJ Center in each row.
	_eps_sig.resize(centers, DoubleArray(centers * 2));
	_shift6.resize(centers, DoubleArray(centers));

	// Construct the parameter tables.
	for (size_t comp_i = 0; comp_i < components.size(); ++comp_i) {
		for (size_t comp_j = 0; comp_j < components.size(); ++comp_j) {
			ParaStrm & p = _domain.getComp2Params()(components[comp_i].ID(),
					components[comp_j].ID());
			p.reset_read();
			for (size_t center_i = 0;
					center_i < components[comp_i].numLJcenters(); ++center_i) {
				for (size_t center_j = 0;
						center_j < components[comp_j].numLJcenters();
						++center_j) {
					if ((components[comp_i].ID() == components[comp_j].ID())
							&& (components[comp_i].numTersoff() > 0
									|| components[comp_j].numTersoff() > 0)) {
						// No LJ interaction between solid atoms of the same component.
						_eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j)] = 0.0;
						_eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j) + 1] = 0.0;
						_shift6[_compIDs[comp_i] + center_i][_compIDs[comp_j] + center_j] = 0.0;
					} else {
						// Extract epsilon*24.0, sigma^2 and shift*6.0 from paramStreams.
						p >> _eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j)];
						p >> _eps_sig[_compIDs[comp_i] + center_i][2 * (_compIDs[comp_j] + center_j) + 1];
						p >> _shift6[_compIDs[comp_i] + center_i][_compIDs[comp_j] + center_j];
					}
				}
			}
		}
	}
}

VectorizedCellProcessor :: ~VectorizedCellProcessor () {
	for (size_t i = 0; i < _particleCellDataVector.size(); ++i) {
		delete _particleCellDataVector[i];
	}
	_particleCellDataVector.clear();
}


void VectorizedCellProcessor::initTraversal(const size_t numCells) {
	_virial = 0.0;
	_upot6lj = 0.0;
	_upotXpoles = 0.0;

	global_log->debug() << "VectorizedLJCellProcessor::initTraversal() to " << numCells << " cells." << std::endl;

	if (numCells > _particleCellDataVector.size()) {
//		_particleCellDataVector.resize(numCells);
		for (size_t i = _particleCellDataVector.size(); i < numCells; i++) {
			_particleCellDataVector.push_back(new CellDataSoA(64,64,64));
		}
		global_log->debug() << "resize CellDataSoA to " << numCells << " cells." << std::endl;
	}

}


void VectorizedCellProcessor::endTraversal() {
	_domain.setLocalVirial(_virial);
	_domain.setLocalUpot(_upot6lj / 6.0 + _upotXpoles);
}


void VectorizedCellProcessor::preprocessCell(ParticleCell & c) {
	assert(!c.getCellDataSoA());

	const MoleculeList & molecules = c.getParticlePointers();

	// Determine the total number of LJ centers.
	size_t numMolecules = molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += molecules[m]->numLJcenters();
		nCharges += molecules[m]->numCharges();
	}

	// Construct the SoA.
	assert(!_particleCellDataVector.empty());
	CellDataSoA* soaPtr = _particleCellDataVector.back();
	global_log->debug() << " _particleCellDataVector.size()=" << _particleCellDataVector.size() << " soaPtr=" << soaPtr << " nLJCenters=" << nLJCenters << std::endl;
	CellDataSoA & soa = *soaPtr;
	soa.setDistLookup(_centers_dist_lookup);
	soa.resize(numMolecules,nLJCenters,nCharges);
	c.setCellDataSoA(soaPtr);
	_particleCellDataVector.pop_back();

	// To get the charges
	ComponentList components = _domain.getComponents();

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	// For each molecule iterate over all its LJ centers and charges.
	for (size_t i = 0; i < molecules.size(); ++i) {
		const size_t mol_num_ljc = molecules[i]->numLJcenters();
		const size_t mol_num_charges = molecules[i]->numCharges();
		const double mol_pos_x = molecules[i]->r(0);
		const double mol_pos_y = molecules[i]->r(1);
		const double mol_pos_z = molecules[i]->r(2);

		soa._mol_pos_x[i] = mol_pos_x;
		soa._mol_pos_y[i] = mol_pos_y;
		soa._mol_pos_z[i] = mol_pos_z;
		soa._mol_num_ljc[i] = mol_num_ljc;
		soa._mol_num_charges[i] = mol_num_charges;

		for (size_t j = 0; j < mol_num_ljc; ++j, ++iLJCenters) {
			// Store a copy of the molecule position for each center, and the position of
			// each LJ center. Assign each LJ center its ID and set the force to 0.0.
			soa._ljc_m_r_x[iLJCenters] = mol_pos_x;
			soa._ljc_m_r_y[iLJCenters] = mol_pos_y;
			soa._ljc_m_r_z[iLJCenters] = mol_pos_z;
			soa._ljc_r_x[iLJCenters] = molecules[i]->ljcenter_d(j)[0] + mol_pos_x;
			soa._ljc_r_y[iLJCenters] = molecules[i]->ljcenter_d(j)[1] + mol_pos_y;
			soa._ljc_r_z[iLJCenters] = molecules[i]->ljcenter_d(j)[2] + mol_pos_z;
			soa._ljc_f_x[iLJCenters] = 0.0;
			soa._ljc_f_y[iLJCenters] = 0.0;
			soa._ljc_f_z[iLJCenters] = 0.0;
			soa._ljc_id[iLJCenters] = _compIDs[molecules[i]->componentid()] + j;
			soa._ljc_dist_lookup[iLJCenters] = 0.0;
		}

		for (size_t j = 0; j < mol_num_charges; ++j, ++iCharges)
		{
			soa._charges_m_r_x[iCharges] = mol_pos_x;
			soa._charges_m_r_y[iCharges] = mol_pos_y;
			soa._charges_m_r_z[iCharges] = mol_pos_z;
			soa._charges_r_x[iCharges] = molecules[i]->charge_d(j)[0] + mol_pos_x;
			soa._charges_r_y[iCharges] = molecules[i]->charge_d(j)[1] + mol_pos_y;
			soa._charges_r_z[iCharges] = molecules[i]->charge_d(j)[2] + mol_pos_z;
			soa._charges_f_x[iCharges] = 0.0;
			soa._charges_f_y[iCharges] = 0.0;
			soa._charges_f_z[iCharges] = 0.0;
			soa._charges_dist_lookup[iCharges] = 0.0;
			// Get the charge
			soa._charges_q[iCharges] = components[molecules[i]->componentid()].charge(j).q();
		}
	}
}


void VectorizedCellProcessor::postprocessCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	CellDataSoA& soa = *c.getCellDataSoA();

	MoleculeList & molecules = c.getParticlePointers();

	// For each molecule iterate over all its centers.
	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t numMols = molecules.size();
	for (size_t m = 0; m < numMols; ++m) {
		const size_t mol_num_ljc = molecules[m]->numLJcenters();
		const size_t mol_num_charges = molecules[m]->numCharges();

		for (size_t i = 0; i < mol_num_ljc; ++i, ++iLJCenters) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._ljc_f_x[iLJCenters];
			f[1] = soa._ljc_f_y[iLJCenters];
			f[2] = soa._ljc_f_z[iLJCenters];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fljcenteradd(i, f);
		}

		for (size_t i = 0; i < mol_num_charges; ++i, ++iCharges) {
			// Store the resulting force in the molecule.
			double f[3];
			f[0] = soa._charges_f_x[iCharges];
			f[1] = soa._charges_f_y[iCharges];
			f[2] = soa._charges_f_z[iCharges];
			assert(!isnan(f[0]));
			assert(!isnan(f[1]));
			assert(!isnan(f[2]));
			molecules[m]->Fchargeadd(i, f);
		}
	}
	// Delete the SoA.
	_particleCellDataVector.push_back(&soa);
	c.setCellDataSoA(0);
}

template<class ForcePolicy, class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovec (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		// Distance vector from center 2 to center 1. Together with the rest of the
		// calculation, we get the proper force from center 1 to center 2.
		// This is done because the calculation of the potential fits in nicely that way.
		const double c_dx = soa1._ljc_r_x[i] - soa2._ljc_r_x[j];
		const double c_dy = soa1._ljc_r_y[i] - soa2._ljc_r_y[j];
		const double c_dz = soa1._ljc_r_z[i] - soa2._ljc_r_z[j];

		const double c_r2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double r2_inv = 1.0 / c_r2;

		const double eps_24 = _eps_sig[soa1._ljc_id[i]][2 * soa2._ljc_id[j]];
		const double sig2 = _eps_sig[soa1._ljc_id[i]][2 * soa2._ljc_id[j] + 1];

		const double lj2 = sig2 * r2_inv;
		const double lj6 = lj2 * lj2 * lj2;
		const double lj12 = lj6 * lj6;
		const double lj12m6 = lj12 - lj6;

		const double scale = eps_24 * r2_inv * (lj12 + lj12m6);

		const double fx = c_dx * scale;
		const double fy = c_dy * scale;
		const double fz = c_dz * scale;

		const double m_dx = soa1._ljc_m_r_x[i] - soa2._ljc_m_r_x[j];
		const double m_dy = soa1._ljc_m_r_y[i] - soa2._ljc_m_r_y[j];
		const double m_dz = soa1._ljc_m_r_z[i] - soa2._ljc_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			_upot6lj += eps_24 * lj12m6 + _shift6[soa1._ljc_id[i]][soa2._ljc_id[j]];
			_virial += m_dx * fx + m_dy * fy + m_dz * fz;
		}
		// Add the force to center 1, and subtract it from center 2.
		soa1._ljc_f_x[i] += fx;
		soa1._ljc_f_y[i] += fy;
		soa1._ljc_f_z[i] += fz;

		soa2._ljc_f_x[j] -= fx;
		soa2._ljc_f_y[j] -= fy;
		soa2._ljc_f_z[j] -= fz;
	}
}  /* end of method VectorizedLJCellProcessor :: _loopBodyNovec */

template<class ForcePolicy, class MacroPolicy>
inline
void VectorizedCellProcessor :: _loopBodyNovecCharges (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask)
{
	// Check if we have to calculate anything for this pair.
	if (*forceMask) {
		// TODO: Fix assuming that 4pie0 = 1 !!!
		const double q1q2per4pie0 = soa1._charges_q[i] * soa2._charges_q[j];

		const double c_dx = soa1._charges_r_x[i] - soa2._charges_r_x[j];
		const double c_dy = soa1._charges_r_y[i] - soa2._charges_r_y[j];
		const double c_dz = soa1._charges_r_z[i] - soa2._charges_r_z[j];

		const double c_dr2 = c_dx * c_dx + c_dy * c_dy + c_dz * c_dz;
		const double c_dr2_inv = 1.0 / c_dr2;
		const double c_dr_inv = sqrt(c_dr2_inv);

		const double upot = q1q2per4pie0 * c_dr_inv;
		const double fac = upot * c_dr2_inv;

		const double f_x = c_dx * fac;
		const double f_y = c_dy * fac;
		const double f_z = c_dz * fac;

		const double m_dx = soa1._charges_m_r_x[i] - soa2._charges_m_r_x[j];
		const double m_dy = soa1._charges_m_r_y[i] - soa2._charges_m_r_y[j];
		const double m_dz = soa1._charges_m_r_z[i] - soa2._charges_m_r_z[j];

		// Check if we have to add the macroscopic values up for this pair.
		if (MacroPolicy :: MacroscopicValueCondition(m_dx, m_dy, m_dz)) {
			_upotXpoles += upot;
			_virial += m_dx * f_x + m_dy * f_y + m_dz * f_z;
		}

		// Add the force to center 1, and subtract it from center 2.
		soa1._charges_f_x[i] += f_x;
		soa1._charges_f_y[i] += f_y;
		soa1._charges_f_z[i] += f_z;

		soa2._charges_f_x[j] -= f_x;
		soa2._charges_f_y[j] -= f_y;
		soa2._charges_f_z[j] -= f_z;
	}
}

#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3

inline void hSum_Add_Store( double * const mem_addr, const __m128d & a ) {
	_mm_store_sd(
			mem_addr,
			_mm_add_sd(_mm_hadd_pd(a, a),
			_mm_load_sd(mem_addr)));
}

#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX

const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);

inline void hSum_Add_Store( double * const mem_addr, const __m256d & a ) {
	const __m256d a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
	const __m256d a_t2 = _mm256_hadd_pd(a, a_t1);
	const __m256d a_t3 = _mm256_hadd_pd(a_t2, a_t2);
	_mm256_maskstore_pd(
			mem_addr,
			memoryMask_first,
			_mm256_add_pd(
				a_t3,
				_mm256_maskload_pd(mem_addr, memoryMask_first)
			)

	);
}

#endif


template<class ForcePolicy, class MacroPolicy>
#if VLJCP_VEC_TYPE==VLJCP_NOVEC
	unsigned long
#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
	__m128d
#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX
	__m256d
#endif
inline VectorizedCellProcessor::calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
		, const __m128d & cutoffRadiusSquareD, size_t end_j, const __m128d m1_r_x, const __m128d m1_r_y, const __m128d m1_r_z
#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		, const __m256d & cutoffRadiusSquareD, size_t end_j, const __m256d m1_r_x, const __m256d m1_r_y, const __m256d m1_r_z
#endif
		) {

#if VLJCP_VEC_TYPE==VLJCP_NOVEC

	unsigned long compute_molecule = 0;

	for (size_t j = ForcePolicy :: InitJ(i_center_idx); j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];
		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		const signed long forceMask = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? (~0l) : 0l;
		compute_molecule |= forceMask;
		*(soa2_center_dist_lookup + j) = forceMask;
	}

	return compute_molecule;

#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3

	__m128d compute_molecule = _mm_setzero_pd();

	// Iterate over centers of second cell
	size_t j = ForcePolicy :: InitJ(i_center_idx);
	for (; j < end_j; j+=2) {
		const __m128d m2_r_x = _mm_load_pd(soa2_m_r_x + j);
		const __m128d m2_r_y = _mm_load_pd(soa2_m_r_y + j);
		const __m128d m2_r_z = _mm_load_pd(soa2_m_r_z + j);

		const __m128d m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
		const __m128d m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
		const __m128d m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

		const __m128d m_dx2 = _mm_mul_pd(m_dx, m_dx);
		const __m128d m_dy2 = _mm_mul_pd(m_dy, m_dy);
		const __m128d m_dz2 = _mm_mul_pd(m_dz, m_dz);

		const __m128d m_r2 = _mm_add_pd(_mm_add_pd(m_dx2, m_dy2), m_dz2);

		const __m128d forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD);
		_mm_store_pd(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = _mm_or_pd(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// TODO: Can't we do this nicer?

		const signed long forceMask_l = ForcePolicy :: Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		// this casting via void* is required for gcc
		const void* forceMask_tmp = reinterpret_cast<const void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* const>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const __m128d forceMask_128 = _mm_set1_pd(forceMask);
		compute_molecule = _mm_or_pd(compute_molecule, forceMask_128);
	}

	return compute_molecule;

#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX

	__m256d compute_molecule = _mm256_setzero_pd();

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	__m256d initJ_mask = ForcePolicy::InitJ_Mask(i_center_idx);
	for (; j < end_j; j+=4) {
		const __m256d m2_r_x = _mm256_load_pd(soa2_m_r_x + j);
		const __m256d m2_r_y = _mm256_load_pd(soa2_m_r_y + j);
		const __m256d m2_r_z = _mm256_load_pd(soa2_m_r_z + j);

		const __m256d m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
		const __m256d m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
		const __m256d m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

		const __m256d m_dx2 = _mm256_mul_pd(m_dx, m_dx);
		const __m256d m_dy2 = _mm256_mul_pd(m_dy, m_dy);
		const __m256d m_dz2 = _mm256_mul_pd(m_dz, m_dz);

		const __m256d m_r2 = _mm256_add_pd(_mm256_add_pd(m_dx2, m_dy2), m_dz2);

		const __m256d forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);
		_mm256_store_pd(soa2_center_dist_lookup + j, forceMask);
		compute_molecule = _mm256_or_pd(compute_molecule, forceMask);
	}

	// End iteration over centers with possible left over center
	for (; j < soa2_num_centers; ++j) {
		const double m_dx = soa1._mol_pos_x[i] - soa2_m_r_x[j];
		const double m_dy = soa1._mol_pos_y[i] - soa2_m_r_y[j];
		const double m_dz = soa1._mol_pos_z[i] - soa2_m_r_z[j];

		const double m_r2 = m_dx * m_dx + m_dy * m_dy + m_dz * m_dz;

		// TODO: Can't we do this nicer?

		signed long forceMask_l;
		if (ForcePolicy::DetectSingleCell()) {
			forceMask_l = (ForcePolicy::Condition(m_r2, cutoffRadiusSquare) && j > i_center_idx) ? ~0l : 0l;
		} else {
			forceMask_l = ForcePolicy::Condition(m_r2, cutoffRadiusSquare) ? ~0l : 0l;
		}

//			this casting via void* is required for gcc
		void* forceMask_tmp = reinterpret_cast<void*>(&forceMask_l);
		double forceMask = *reinterpret_cast<double const* const>(forceMask_tmp);

		*(soa2_center_dist_lookup + j) = forceMask;
		const __m256d forceMask_256 = _mm256_set1_pd(forceMask);
		compute_molecule = _mm256_or_pd(compute_molecule, forceMask_256);
	}

	return compute_molecule;

#endif

}

template<class ForcePolicy, class MacroPolicy>
void VectorizedCellProcessor::_calculatePairs(const CellDataSoA & soa1,
		const CellDataSoA & soa2) {
#if VLJCP_VEC_TYPE==VLJCP_NOVEC
	// For the unvectorized version, we only have to iterate over all pairs of
	// LJ centers and apply the unvectorized loop body.
	size_t i_ljc_idx = 0;
	size_t i_charge_idx = 0;
	for (size_t i = 0; i < soa1._num_molecules; ++i) {
		// Computation of LJ interaction
		unsigned long compute_molecule_lj = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_ljc_idx, soa2._num_ljcenters, _LJcutoffRadiusSquare,
				soa2._ljc_dist_lookup, soa2._ljc_m_r_x, soa2._ljc_m_r_y, soa2._ljc_m_r_z);
		unsigned long compute_molecule_charge  = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_charge_idx, soa2._num_charges, _cutoffRadiusSquare,
				soa2._charges_dist_lookup, soa2._charges_m_r_x, soa2._charges_m_r_y, soa2._charges_m_r_z);

		if (!compute_molecule_lj) {
			i_ljc_idx += soa1._mol_num_ljc[i];
		}
		else {
			for (int local_i = 0; local_i < soa1._mol_num_ljc[i]; local_i++ ) {
				for (size_t j = ForcePolicy :: InitJ(i_ljc_idx); j < soa2._num_ljcenters; ++j) {
					_loopBodyNovec<CellPairPolicy_, MacroPolicy>(soa1, i_ljc_idx, soa2, j, soa2._ljc_dist_lookup + j);
				}
				i_ljc_idx++;
			}
		}

		// Computation of Charge-Charge interaction

		if (!compute_molecule_charge) {
			i_charge_idx += soa1._mol_num_charges[i];
		}
		else {
			for (int local_i = 0; local_i < soa1._mol_num_charges[i]; local_i++ ) {
				for (size_t j = ForcePolicy :: InitJ(i_charge_idx); j < soa2._num_charges; ++j) {
					_loopBodyNovecCharges<CellPairPolicy_, MacroPolicy>(soa1, i_charge_idx, soa2, j, soa2._charges_dist_lookup + j);
				}
				i_charge_idx++;
			}
		}
	}


#else
	// Pointer for molecules
	const double * const soa1_mol_pos_x = soa1._mol_pos_x;
	const double * const soa1_mol_pos_y = soa1._mol_pos_y;
	const double * const soa1_mol_pos_z = soa1._mol_pos_z;

	// Pointer for LJ centers
	const double * const soa1_ljc_r_x = soa1._ljc_r_x;
	const double * const soa1_ljc_r_y = soa1._ljc_r_y;
	const double * const soa1_ljc_r_z = soa1._ljc_r_z;
	double * const soa1_ljc_f_x = soa1._ljc_f_x;
	double * const soa1_ljc_f_y = soa1._ljc_f_y;
	double * const soa1_ljc_f_z = soa1._ljc_f_z;
	const int * const soa1_mol_num_ljc = soa1._mol_num_ljc;
	const size_t * const soa1_ljc_id = soa1._ljc_id;

	const double * const soa2_ljc_m_r_x = soa2._ljc_m_r_x;
	const double * const soa2_ljc_m_r_y = soa2._ljc_m_r_y;
	const double * const soa2_ljc_m_r_z = soa2._ljc_m_r_z;
	const double * const soa2_ljc_r_x = soa2._ljc_r_x;
	const double * const soa2_ljc_r_y = soa2._ljc_r_y;
	const double * const soa2_ljc_r_z = soa2._ljc_r_z;
	double * const soa2_ljc_f_x = soa2._ljc_f_x;
	double * const soa2_ljc_f_y = soa2._ljc_f_y;
	double * const soa2_ljc_f_z = soa2._ljc_f_z;
	const size_t * const soa2_ljc_id = soa2._ljc_id;

	double* const soa2_ljc_dist_lookup = soa2._ljc_dist_lookup;

	// Pointer for charges
	const double * const soa1_charges_r_x = soa1._charges_r_x;
	const double * const soa1_charges_r_y = soa1._charges_r_y;
	const double * const soa1_charges_r_z = soa1._charges_r_z;
	double * const soa1_charges_f_x = soa1._charges_f_x;
	double * const soa1_charges_f_y = soa1._charges_f_y;
	double * const soa1_charges_f_z = soa1._charges_f_z;
	const double * const soa1_charges_q = soa1._charges_q;

	const double * const soa2_charges_m_r_x = soa2._charges_m_r_x;
	const double * const soa2_charges_m_r_y = soa2._charges_m_r_y;
	const double * const soa2_charges_m_r_z = soa2._charges_m_r_z;
	const double * const soa2_charges_r_x = soa2._charges_r_x;
	const double * const soa2_charges_r_y = soa2._charges_r_y;
	const double * const soa2_charges_r_z = soa2._charges_r_z;
	double * const soa2_charges_f_x = soa2._charges_f_x;
	double * const soa2_charges_f_y = soa2._charges_f_y;
	double * const soa2_charges_f_z = soa2._charges_f_z;
	const int * const soa1_mol_num_charges = soa1._mol_num_charges;
	const double * const soa2_charges_q = soa2._charges_q;

	double* const soa2_charges_dist_lookup = soa2._charges_dist_lookup;

#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3

	const __m128d one = _mm_set1_pd(1.0);

	__m128d sum_upot6lj = _mm_setzero_pd();
	__m128d sum_upotXpoles = _mm_setzero_pd();
	__m128d sum_virial = _mm_setzero_pd();

	const __m128d cutoffRadiusSquare = _mm_set1_pd(_cutoffRadiusSquare);
	const size_t end_charges_j = soa2._num_charges & (~1);
	size_t i_charge_idx = 0;

	const size_t end_ljc_j = soa2._num_ljcenters & (~1);
	const __m128d rc2 = _mm_set1_pd(_LJcutoffRadiusSquare);
	size_t i_ljc_idx = 0;

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._num_molecules; ++i) {
		const __m128d m1_r_x = _mm_loaddup_pd(soa1_mol_pos_x + i);
		const __m128d m1_r_y = _mm_loaddup_pd(soa1_mol_pos_y + i);
		const __m128d m1_r_z = _mm_loaddup_pd(soa1_mol_pos_z + i);


		const __m128d compute_molecule_ljc = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_ljc_idx, soa2._num_ljcenters, _LJcutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const __m128d compute_molecule_charges = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_charge_idx, soa2._num_charges, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);

		if (!_mm_movemask_pd(compute_molecule_ljc)) {
			i_ljc_idx += soa1_mol_num_ljc[i];
		}
		else {

			// actual force computation
			for (int local_i = 0; local_i < soa1_mol_num_ljc[i]; local_i++ ) {
				__m128d sum_fx1 = _mm_setzero_pd();
				__m128d sum_fy1 = _mm_setzero_pd();
				__m128d sum_fz1 = _mm_setzero_pd();
				const __m128d c_r_x1 = _mm_loaddup_pd(soa1_ljc_r_x + i_ljc_idx);
				const __m128d c_r_y1 = _mm_loaddup_pd(soa1_ljc_r_y + i_ljc_idx);
				const __m128d c_r_z1 = _mm_loaddup_pd(soa1_ljc_r_z + i_ljc_idx);
				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ(i_ljc_idx);
				for (; j < end_ljc_j; j += 2) {
					const __m128d forceMask = _mm_load_pd(soa2_ljc_dist_lookup + j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (_mm_movemask_pd(forceMask) > 0) {
						const __m128d c_r_x2 = _mm_load_pd(soa2_ljc_r_x + j);
						const __m128d c_dx = _mm_sub_pd(c_r_x1, c_r_x2);
						const __m128d c_r_y2 = _mm_load_pd(soa2_ljc_r_y + j);
						const __m128d c_dy = _mm_sub_pd(c_r_y1, c_r_y2);
						const __m128d c_r_z2 = _mm_load_pd(soa2_ljc_r_z + j);
						const __m128d c_dz = _mm_sub_pd(c_r_z1, c_r_z2);
						const __m128d c_dxdx = _mm_mul_pd(c_dx, c_dx);
						const __m128d c_dydy = _mm_mul_pd(c_dy, c_dy);
						const __m128d c_dzdz = _mm_mul_pd(c_dz, c_dz);
						const __m128d c_dxdx_dydy = _mm_add_pd(c_dxdx, c_dydy);
						const __m128d c_r2 = _mm_add_pd(c_dxdx_dydy, c_dzdz);
						const __m128d r2_inv_unmasked = _mm_div_pd(one, c_r2);
						const __m128d r2_inv = _mm_and_pd(r2_inv_unmasked, forceMask);
						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						const size_t id_j0 = soa2_ljc_id[j];
						const size_t id_j1 = soa2_ljc_id[j + 1];
						const __m128d e1s1 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j0);
						const __m128d e2s2 = _mm_load_pd(_eps_sig[id_i] + 2 * id_j1);
						const __m128d eps_24 = _mm_unpacklo_pd(e1s1, e2s2);
						const __m128d sig2 = _mm_unpackhi_pd(e1s1, e2s2);
						const __m128d lj2 = _mm_mul_pd(sig2, r2_inv);
						const __m128d lj4 = _mm_mul_pd(lj2, lj2);
						const __m128d lj6 = _mm_mul_pd(lj4, lj2);
						const __m128d lj12 = _mm_mul_pd(lj6, lj6);
						const __m128d lj12m6 = _mm_sub_pd(lj12, lj6);
						const __m128d eps24r2inv = _mm_mul_pd(eps_24, r2_inv);
						const __m128d lj12lj12m6 = _mm_add_pd(lj12, lj12m6);
						const __m128d scale = _mm_mul_pd(eps24r2inv, lj12lj12m6);
						const __m128d fx = _mm_mul_pd(c_dx, scale);
						const __m128d fy = _mm_mul_pd(c_dy, scale);
						const __m128d fz = _mm_mul_pd(c_dz, scale);

						const __m128d m_r_x2 = _mm_load_pd(soa2_ljc_m_r_x + j);
						const __m128d m_dx = _mm_sub_pd(m1_r_x, m_r_x2);
						const __m128d m_r_y2 = _mm_load_pd(soa2_ljc_m_r_y + j);
						const __m128d m_dy = _mm_sub_pd(m1_r_y, m_r_y2);
						const __m128d m_r_z2 = _mm_load_pd(soa2_ljc_m_r_z + j);
						const __m128d m_dz = _mm_sub_pd(m1_r_z, m_r_z2);

						const __m128d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
						// Only go on if at least 1 macroscopic value has to be calculated.
						if (_mm_movemask_pd(macroMask) > 0) {
							const __m128d sh1 = _mm_load_sd(_shift6[id_i] + id_j0);
							const __m128d sh2 = _mm_load_sd(_shift6[id_i] + id_j1);
							const __m128d shift6 = _mm_unpacklo_pd(sh1, sh2);
							const __m128d upot = _mm_mul_pd(eps_24, lj12m6);
							const __m128d upot_sh = _mm_add_pd(shift6, upot);
							const __m128d upot_masked = _mm_and_pd(upot_sh, macroMask);
							sum_upot6lj = _mm_add_pd(sum_upot6lj, upot_masked);
							const __m128d vir_x = _mm_mul_pd(m_dx, fx);
							const __m128d vir_y = _mm_mul_pd(m_dy, fy);
							const __m128d vir_z = _mm_mul_pd(m_dz, fz);
							const __m128d vir_xy = _mm_add_pd(vir_x, vir_y);
							const __m128d virial = _mm_add_pd(vir_xy, vir_z);
							const __m128d vir_masked = _mm_and_pd(virial, macroMask);
							sum_virial = _mm_add_pd(sum_virial, vir_masked);
						}
						const __m128d old_fx2 = _mm_load_pd(soa2_ljc_f_x + j);
						const __m128d new_fx2 = _mm_sub_pd(old_fx2, fx);
						_mm_store_pd(soa2_ljc_f_x + j, new_fx2);
						const __m128d old_fy2 = _mm_load_pd(soa2_ljc_f_y + j);
						const __m128d new_fy2 = _mm_sub_pd(old_fy2, fy);
						_mm_store_pd(soa2_ljc_f_y + j, new_fy2);
						const __m128d old_fz2 = _mm_load_pd(soa2_ljc_f_z + j);
						const __m128d new_fz2 = _mm_sub_pd(old_fz2, fz);
						_mm_store_pd(soa2_ljc_f_z + j, new_fz2);
						sum_fx1 = _mm_add_pd(sum_fx1, fx);
						sum_fy1 = _mm_add_pd(sum_fy1, fy);
						sum_fz1 = _mm_add_pd(sum_fz1, fz);
					}
				}

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				// Unvectorized calculation for leftover pairs.
				switch (soa2._num_ljcenters & 1) {
					case 1: {
						_loopBodyNovec<ForcePolicy, MacroPolicy>(soa1, i_ljc_idx, soa2, end_ljc_j, soa2_ljc_dist_lookup + j);
					}
					break;
				}

				i_ljc_idx++;
			}
		}

		// Continue with next molecule if no force has to be calculated
		if (!_mm_movemask_pd(compute_molecule_charges)) {
			i_charge_idx += soa1_mol_num_charges[i];
			continue;
		}
		else {
			// Force calculation

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_num_charges[i]; local_i++ ) {

				const __m128d c1_q = _mm_loaddup_pd(soa1_charges_q + i_charge_idx);
				const __m128d c1_r_x = _mm_loaddup_pd(soa1_charges_r_x + i_charge_idx);
				const __m128d c1_r_y = _mm_loaddup_pd(soa1_charges_r_y + i_charge_idx);
				const __m128d c1_r_z = _mm_loaddup_pd(soa1_charges_r_z + i_charge_idx);

				__m128d sum_c1_f_x = _mm_setzero_pd();
				__m128d sum_c1_f_y = _mm_setzero_pd();
				__m128d sum_c1_f_z = _mm_setzero_pd();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 2) {
					const __m128d forceMask = _mm_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm_movemask_pd(forceMask) > 0) {
						const __m128d c2_q = _mm_load_pd(soa2_charges_q + j);
						// TODO: Fix assuming that 4pie0 = 1 !!!
						const __m128d q1q2per4pie0 = _mm_mul_pd(c1_q, c2_q);

						const __m128d c2_r_x = _mm_load_pd(soa2_charges_r_x + j);
						const __m128d c2_r_y = _mm_load_pd(soa2_charges_r_y + j);
						const __m128d c2_r_z = _mm_load_pd(soa2_charges_r_z + j);

						const __m128d c_dx = _mm_sub_pd(c1_r_x, c2_r_x);
						const __m128d c_dy = _mm_sub_pd(c1_r_y, c2_r_y);
						const __m128d c_dz = _mm_sub_pd(c1_r_z, c2_r_z);

						const __m128d c_dx2 = _mm_mul_pd(c_dx, c_dx);
						const __m128d c_dy2 = _mm_mul_pd(c_dy, c_dy);
						const __m128d c_dz2 = _mm_mul_pd(c_dz, c_dz);

						const __m128d c_dr2 = _mm_add_pd(_mm_add_pd(c_dx2, c_dy2), c_dz2);
						const __m128d c_dr2_inv_unmasked = _mm_div_pd(one, c_dr2);
						const __m128d c_dr2_inv = _mm_and_pd(c_dr2_inv_unmasked, forceMask);
						const __m128d c_dr_inv = _mm_sqrt_pd(c_dr2_inv);

						const __m128d upot = _mm_mul_pd(q1q2per4pie0, c_dr_inv);
						const __m128d fac = _mm_mul_pd(upot, c_dr2_inv);

						const __m128d f_x = _mm_mul_pd(c_dx, fac);
						const __m128d f_y = _mm_mul_pd(c_dy, fac);
						const __m128d f_z = _mm_mul_pd(c_dz, fac);

						const __m128d m2_r_x = _mm_load_pd(soa2_charges_m_r_x + j);
						const __m128d m2_r_y = _mm_load_pd(soa2_charges_m_r_y + j);
						const __m128d m2_r_z = _mm_load_pd(soa2_charges_m_r_z + j);

						const __m128d m_dx = _mm_sub_pd(m1_r_x, m2_r_x);
						const __m128d m_dy = _mm_sub_pd(m1_r_y, m2_r_y);
						const __m128d m_dz = _mm_sub_pd(m1_r_z, m2_r_z);

						const __m128d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
						// Check if we have to add the macroscopic values up for at least one of this pairs
						if (_mm_movemask_pd(macroMask) > 0) {
							const __m128d upot_masked = _mm_and_pd(upot, macroMask);
							sum_upotXpoles = _mm_add_pd(sum_upotXpoles, upot_masked);

							const __m128d virial_x = _mm_mul_pd(m_dx, f_x);
							const __m128d virial_y = _mm_mul_pd(m_dy, f_y);
							const __m128d virial_z = _mm_mul_pd(m_dz, f_z);

							const __m128d virial = _mm_add_pd(_mm_add_pd(virial_x, virial_y), virial_z);
							const __m128d virial_masked = _mm_and_pd(virial, macroMask);
							sum_virial = _mm_add_pd(sum_virial, virial_masked);

						}

						sum_c1_f_x = _mm_add_pd(sum_c1_f_x, f_x);
						sum_c1_f_y = _mm_add_pd(sum_c1_f_y, f_y);
						sum_c1_f_z = _mm_add_pd(sum_c1_f_z, f_z);

						__m128d c2_f_x = _mm_load_pd(soa2_charges_f_x + j);
						__m128d c2_f_y = _mm_load_pd(soa2_charges_f_y + j);
						__m128d c2_f_z = _mm_load_pd(soa2_charges_f_z + j);

						c2_f_x = _mm_sub_pd(c2_f_x, f_x);
						c2_f_y = _mm_sub_pd(c2_f_y, f_y);
						c2_f_z = _mm_sub_pd(c2_f_z, f_z);

						_mm_store_pd(soa2_charges_f_x + j, c2_f_x);
						_mm_store_pd(soa2_charges_f_y + j, c2_f_y);
						_mm_store_pd(soa2_charges_f_z + j, c2_f_z);
					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx, sum_c1_f_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx, sum_c1_f_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx, sum_c1_f_z);


				// End iteration over centers with possible left over center
				for (; j < soa2._num_charges; ++j) {
					_loopBodyNovecCharges<ForcePolicy, MacroPolicy>(soa1, i_charge_idx, soa2, j, soa2_charges_dist_lookup + j);
				}

				i_charge_idx++;
			}
		}
	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);

#elif VLJCP_VEC_TYPE==VLJCP_VEC_AVX


	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);
	const __m256d one = _mm256_set1_pd(1.0);

	__m256d sum_upot6lj = _mm256_setzero_pd();
	__m256d sum_upotXpoles = _mm256_setzero_pd();
	__m256d sum_virial = _mm256_setzero_pd();

	const __m256d cutoffRadiusSquare = _mm256_set1_pd(_cutoffRadiusSquare);
	const size_t end_charges_j = soa2._num_charges & (~3);
	size_t i_charge_idx = 0;

	// Begin LJ interaction

	const size_t end_ljc_j = soa2._num_ljcenters & ~static_cast<size_t>(3);
	const __m256d rc2 = _mm256_set1_pd(_LJcutoffRadiusSquare);

	size_t i_ljc_idx = 0;

	// Iterate over each center in the first cell.
	for (size_t i = 0; i < soa1._num_molecules; ++i) {
		const __m256d m1_r_x = _mm256_broadcast_sd(soa1_mol_pos_x + i);
		const __m256d m1_r_y = _mm256_broadcast_sd(soa1_mol_pos_y + i);
		const __m256d m1_r_z = _mm256_broadcast_sd(soa1_mol_pos_z + i);

		// Iterate over centers of second cell
		const __m256d compute_molecule_ljc = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_ljc_idx, soa2._num_ljcenters, _LJcutoffRadiusSquare,
				soa2_ljc_dist_lookup, soa2_ljc_m_r_x, soa2_ljc_m_r_y, soa2_ljc_m_r_z,
				rc2, end_ljc_j, m1_r_x, m1_r_y, m1_r_z);
		const __m256d compute_molecule_charges = calcDistLookup<ForcePolicy, MacroPolicy>(soa1, i, i_charge_idx, soa2._num_charges, _cutoffRadiusSquare,
				soa2_charges_dist_lookup, soa2_charges_m_r_x, soa2_charges_m_r_y, soa2_charges_m_r_z,
				cutoffRadiusSquare,	end_charges_j, m1_r_x, m1_r_y, m1_r_z);

		if (!_mm256_movemask_pd(compute_molecule_ljc)) {
			i_ljc_idx += soa1_mol_num_ljc[i];
		}
		else {
			// actual force computation
			for (int local_i = 0; local_i < soa1_mol_num_ljc[i]; local_i++ ) {
				__m256d sum_fx1 = _mm256_setzero_pd();
				__m256d sum_fy1 = _mm256_setzero_pd();
				__m256d sum_fz1 = _mm256_setzero_pd();
				const __m256d c_r_x1 = _mm256_broadcast_sd(soa1_ljc_r_x + i_ljc_idx);
				const __m256d c_r_y1 = _mm256_broadcast_sd(soa1_ljc_r_y + i_ljc_idx);
				const __m256d c_r_z1 = _mm256_broadcast_sd(soa1_ljc_r_z + i_ljc_idx);
				// Iterate over each pair of centers in the second cell.
				size_t j = ForcePolicy::InitJ(i_ljc_idx);
				for (; j < end_ljc_j; j += 4) {
					const __m256d forceMask = _mm256_load_pd(soa2_ljc_dist_lookup + j);
					// Only go on if at least 1 of the forces has to be calculated.
					if (_mm256_movemask_pd(forceMask) > 0) {
						const __m256d c_r_x2 = _mm256_load_pd(soa2_ljc_r_x + j);
						const __m256d c_dx = _mm256_sub_pd(c_r_x1, c_r_x2);
						const __m256d c_r_y2 = _mm256_load_pd(soa2_ljc_r_y + j);
						const __m256d c_dy = _mm256_sub_pd(c_r_y1, c_r_y2);
						const __m256d c_r_z2 = _mm256_load_pd(soa2_ljc_r_z + j);
						const __m256d c_dz = _mm256_sub_pd(c_r_z1, c_r_z2);
						const __m256d c_dxdx = _mm256_mul_pd(c_dx, c_dx);
						const __m256d c_dydy = _mm256_mul_pd(c_dy, c_dy);
						const __m256d c_dzdz = _mm256_mul_pd(c_dz, c_dz);
						const __m256d c_dxdx_dydy = _mm256_add_pd(c_dxdx, c_dydy);
						const __m256d c_r2 = _mm256_add_pd(c_dxdx_dydy, c_dzdz);
						const __m256d r2_inv_unmasked = _mm256_div_pd(one, c_r2);
						const __m256d r2_inv = _mm256_and_pd(r2_inv_unmasked, forceMask);

						const size_t id_i = soa1_ljc_id[i_ljc_idx];
						const size_t id_j0 = soa2_ljc_id[j];
						const size_t id_j1 = soa2_ljc_id[j + 1];
						const size_t id_j2 = soa2_ljc_id[j + 2];
						const size_t id_j3 = soa2_ljc_id[j + 3];

						const __m256d e0s0 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j0, memoryMask_first_second);
						const __m256d e1s1 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j1, memoryMask_first_second);
						const __m256d e2s2 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j2, memoryMask_first_second);
						const __m256d e3s3 = _mm256_maskload_pd(_eps_sig[id_i] + 2 * id_j3, memoryMask_first_second);

						const __m256d e0e1 = _mm256_unpacklo_pd(e0s0, e1s1);
						const __m256d s0s1 = _mm256_unpackhi_pd(e0s0, e1s1);
						const __m256d e2e3 = _mm256_unpacklo_pd(e2s2, e3s3);
						const __m256d s2s3 = _mm256_unpackhi_pd(e2s2, e3s3);

						const __m256d eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
						const __m256d sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);

						const __m256d lj2 = _mm256_mul_pd(sig2, r2_inv);
						const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
						const __m256d lj6 = _mm256_mul_pd(lj4, lj2);
						const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
						const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);

						const __m256d eps24r2inv = _mm256_mul_pd(eps_24, r2_inv);
						const __m256d lj12lj12m6 = _mm256_add_pd(lj12, lj12m6);
						const __m256d scale = _mm256_mul_pd(eps24r2inv, lj12lj12m6);

						const __m256d fx = _mm256_mul_pd(c_dx, scale);
						const __m256d fy = _mm256_mul_pd(c_dy, scale);
						const __m256d fz = _mm256_mul_pd(c_dz, scale);

						const __m256d m_r_x2 = _mm256_load_pd(soa2_ljc_m_r_x + j);
						const __m256d m_dx = _mm256_sub_pd(m1_r_x, m_r_x2);
						const __m256d m_r_y2 = _mm256_load_pd(soa2_ljc_m_r_y + j);
						const __m256d m_dy = _mm256_sub_pd(m1_r_y, m_r_y2);
						const __m256d m_r_z2 = _mm256_load_pd(soa2_ljc_m_r_z + j);
						const __m256d m_dz = _mm256_sub_pd(m1_r_z, m_r_z2);

						const __m256d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);

						// Only go on if at least 1 macroscopic value has to be calculated.
						if (_mm256_movemask_pd(macroMask) > 0) {
							const __m256d sh0 = _mm256_maskload_pd(_shift6[id_i] + id_j0, memoryMask_first);
							const __m256d sh1 = _mm256_maskload_pd(_shift6[id_i] + id_j1, memoryMask_first);
							const __m256d sh2 = _mm256_maskload_pd(_shift6[id_i] + id_j2, memoryMask_first);
							const __m256d sh3 = _mm256_maskload_pd(_shift6[id_i] + id_j3, memoryMask_first);

							const __m256d sh0sh1 = _mm256_unpacklo_pd(sh0, sh1);
							const __m256d sh2sh3 = _mm256_unpacklo_pd(sh2, sh3);

							const __m256d shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);

							const __m256d upot = _mm256_mul_pd(eps_24, lj12m6);
							const __m256d upot_sh = _mm256_add_pd(shift6, upot);
							const __m256d upot_masked = _mm256_and_pd(upot_sh, macroMask);

							sum_upot6lj = _mm256_add_pd(sum_upot6lj, upot_masked);

							const __m256d vir_x = _mm256_mul_pd(m_dx, fx);
							const __m256d vir_y = _mm256_mul_pd(m_dy, fy);
							const __m256d vir_z = _mm256_mul_pd(m_dz, fz);

							const __m256d vir_xy = _mm256_add_pd(vir_x, vir_y);
							const __m256d virial = _mm256_add_pd(vir_xy, vir_z);
							const __m256d vir_masked = _mm256_and_pd(virial, macroMask);

							sum_virial = _mm256_add_pd(sum_virial, vir_masked);
						}
						const __m256d old_fx2 = _mm256_load_pd(soa2_ljc_f_x + j);
						const __m256d new_fx2 = _mm256_sub_pd(old_fx2, fx);
						_mm256_store_pd(soa2_ljc_f_x + j, new_fx2);
						const __m256d old_fy2 = _mm256_load_pd(soa2_ljc_f_y + j);
						const __m256d new_fy2 = _mm256_sub_pd(old_fy2, fy);
						_mm256_store_pd(soa2_ljc_f_y + j, new_fy2);
						const __m256d old_fz2 = _mm256_load_pd(soa2_ljc_f_z + j);
						const __m256d new_fz2 = _mm256_sub_pd(old_fz2, fz);
						_mm256_store_pd(soa2_ljc_f_z + j, new_fz2);
						sum_fx1 = _mm256_add_pd(sum_fx1, fx);
						sum_fy1 = _mm256_add_pd(sum_fy1, fy);
						sum_fz1 = _mm256_add_pd(sum_fz1, fz);
					}
				}

				hSum_Add_Store(soa1_ljc_f_x + i_ljc_idx, sum_fx1);
				hSum_Add_Store(soa1_ljc_f_y + i_ljc_idx, sum_fy1);
				hSum_Add_Store(soa1_ljc_f_z + i_ljc_idx, sum_fz1);

				// Unvectorized calculation for leftover pairs.
				for (; j < soa2._num_ljcenters; ++j) {
					_loopBodyNovec<ForcePolicy, MacroPolicy>(soa1, i_ljc_idx, soa2, j, soa2_ljc_dist_lookup + j);
				}

				i_ljc_idx++;
			}
		}

		// Continue with next molecule if no force has to be calculated
		if (!_mm256_movemask_pd(compute_molecule_charges)) {
			i_charge_idx += soa1_mol_num_charges[i];
		}
		else {
			// Force calculation

			// Iterate over centers of actual molecule
			for (int local_i = 0; local_i < soa1_mol_num_charges[i]; local_i++ ) {

				const __m256d c1_q = _mm256_broadcast_sd(soa1_charges_q + i_charge_idx);
				const __m256d c1_r_x = _mm256_broadcast_sd(soa1_charges_r_x + i_charge_idx);
				const __m256d c1_r_y = _mm256_broadcast_sd(soa1_charges_r_y + i_charge_idx);
				const __m256d c1_r_z = _mm256_broadcast_sd(soa1_charges_r_z + i_charge_idx);

				__m256d sum_c1_f_x = _mm256_setzero_pd();
				__m256d sum_c1_f_y = _mm256_setzero_pd();
				__m256d sum_c1_f_z = _mm256_setzero_pd();

				// Iterate over centers of second cell
				size_t j = ForcePolicy::InitJ(i_charge_idx);
				for (; j < end_charges_j; j += 4) {
					const __m256d forceMask = _mm256_load_pd(soa2_charges_dist_lookup + j);
					// Check if we have to calculate anything for at least one of the pairs
					if (_mm256_movemask_pd(forceMask) > 0) {

						const __m256d c2_q = _mm256_load_pd(soa2_charges_q + j);
						// TODO: Fix assuming that 4pie0 = 1 !!!
						const __m256d q1q2per4pie0 = _mm256_mul_pd(c1_q, c2_q);

						const __m256d c2_r_x = _mm256_load_pd(soa2_charges_r_x + j);
						const __m256d c2_r_y = _mm256_load_pd(soa2_charges_r_y + j);
						const __m256d c2_r_z = _mm256_load_pd(soa2_charges_r_z + j);

						const __m256d c_dx = _mm256_sub_pd(c1_r_x, c2_r_x);
						const __m256d c_dy = _mm256_sub_pd(c1_r_y, c2_r_y);
						const __m256d c_dz = _mm256_sub_pd(c1_r_z, c2_r_z);

						const __m256d c_dx2 = _mm256_mul_pd(c_dx, c_dx);
						const __m256d c_dy2 = _mm256_mul_pd(c_dy, c_dy);
						const __m256d c_dz2 = _mm256_mul_pd(c_dz, c_dz);

						const __m256d c_dr2 = _mm256_add_pd(_mm256_add_pd(c_dx2, c_dy2), c_dz2);
						const __m256d c_dr2_inv_unmasked = _mm256_div_pd(one, c_dr2);
						const __m256d c_dr2_inv = _mm256_and_pd(c_dr2_inv_unmasked, forceMask);
						const __m256d c_dr_inv = _mm256_sqrt_pd(c_dr2_inv);

						const __m256d upot = _mm256_mul_pd(q1q2per4pie0, c_dr_inv);
						const __m256d fac = _mm256_mul_pd(upot, c_dr2_inv);

						const __m256d f_x = _mm256_mul_pd(c_dx, fac);
						const __m256d f_y = _mm256_mul_pd(c_dy, fac);
						const __m256d f_z = _mm256_mul_pd(c_dz, fac);

						const __m256d m2_r_x = _mm256_load_pd(soa2_charges_m_r_x + j);
						const __m256d m2_r_y = _mm256_load_pd(soa2_charges_m_r_y + j);
						const __m256d m2_r_z = _mm256_load_pd(soa2_charges_m_r_z + j);

						const __m256d m_dx = _mm256_sub_pd(m1_r_x, m2_r_x);
						const __m256d m_dy = _mm256_sub_pd(m1_r_y, m2_r_y);
						const __m256d m_dz = _mm256_sub_pd(m1_r_z, m2_r_z);

						const __m256d macroMask = MacroPolicy::GetMacroMask(forceMask, m_dx, m_dy, m_dz);
						// Check if we have to add the macroscopic values up for at least one of this pairs
						if (_mm256_movemask_pd(macroMask) > 0) {
							const __m256d upot_masked = _mm256_and_pd(upot, macroMask);
							sum_upotXpoles = _mm256_add_pd(sum_upotXpoles, upot_masked);

							const __m256d virial_x = _mm256_mul_pd(m_dx, f_x);
							const __m256d virial_y = _mm256_mul_pd(m_dy, f_y);
							const __m256d virial_z = _mm256_mul_pd(m_dz, f_z);

							const __m256d virial = _mm256_add_pd(_mm256_add_pd(virial_x,virial_y),virial_z);
							const __m256d virial_masked = _mm256_and_pd(virial, macroMask);
							sum_virial = _mm256_add_pd(sum_virial, virial_masked);

						}

						sum_c1_f_x = _mm256_add_pd(sum_c1_f_x, f_x);
						sum_c1_f_y = _mm256_add_pd(sum_c1_f_y, f_y);
						sum_c1_f_z = _mm256_add_pd(sum_c1_f_z, f_z);

						__m256d c2_f_x = _mm256_load_pd(soa2_charges_f_x + j);
						__m256d c2_f_y = _mm256_load_pd(soa2_charges_f_y + j);
						__m256d c2_f_z = _mm256_load_pd(soa2_charges_f_z + j);

						c2_f_x = _mm256_sub_pd(c2_f_x, f_x);
						c2_f_y = _mm256_sub_pd(c2_f_y, f_y);
						c2_f_z = _mm256_sub_pd(c2_f_z, f_z);

						_mm256_store_pd(soa2_charges_f_x + j, c2_f_x);
						_mm256_store_pd(soa2_charges_f_y + j, c2_f_y);
						_mm256_store_pd(soa2_charges_f_z + j, c2_f_z);

					}
				}

				// Add old force and summed calculated forces for center 1
				hSum_Add_Store(soa1_charges_f_x + i_charge_idx, sum_c1_f_x);
				hSum_Add_Store(soa1_charges_f_y + i_charge_idx, sum_c1_f_y);
				hSum_Add_Store(soa1_charges_f_z + i_charge_idx, sum_c1_f_z);

				// End iteration over centers with possible left over center
				for (; j < soa2._num_charges; ++j) {
					_loopBodyNovecCharges<ForcePolicy, MacroPolicy>(soa1, i_charge_idx, soa2, j, soa2_charges_dist_lookup + j);
				}

				i_charge_idx++;
			}
		}
	}

	hSum_Add_Store(&_upot6lj, sum_upot6lj);
	hSum_Add_Store(&_upotXpoles, sum_upotXpoles);
	hSum_Add_Store(&_virial, sum_virial);

#endif
#endif
} // void LennardJonesCellHandler::CalculatePairs_(LJSoA & soa1, LJSoA & soa2)

void VectorizedCellProcessor::processCell(ParticleCell & c) {
	assert(c.getCellDataSoA());
	if (c.isHaloCell() || ((c.getCellDataSoA()->_num_ljcenters < 2) && (c.getCellDataSoA()->_num_charges < 2)))
		return;

	_calculatePairs<SingleCellPolicy_, AllMacroPolicy_>(*(c.getCellDataSoA()), *(c.getCellDataSoA()));
}

void VectorizedCellProcessor::processCellPair(ParticleCell & c1,
		ParticleCell & c2) {
	assert(&c1 != &c2);
	assert(c1.getCellDataSoA());
	assert(c2.getCellDataSoA());

	if (((c1.getCellDataSoA()->_num_ljcenters == 0) || (c2.getCellDataSoA()->_num_ljcenters == 0)) &&
			((c1.getCellDataSoA()->_num_charges == 0) || (c2.getCellDataSoA()->_num_charges == 0))) {
		return;
	}

	if (!(c1.isHaloCell() || c2.isHaloCell())) {
		_calculatePairs<CellPairPolicy_, AllMacroPolicy_>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
	} else if (c1.isHaloCell() == (!c2.isHaloCell())) {
		_calculatePairs<CellPairPolicy_, SomeMacroPolicy_>(*(c1.getCellDataSoA()), *(c2.getCellDataSoA()));
	} else {
		return;
	}
}
