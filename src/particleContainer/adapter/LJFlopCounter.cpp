/**
 * \file
 * \brief A CellProcessor that produces Flop information.
 * \author Johannes Heckl
 */

#include "LJFlopCounter.h"

#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

#define VLENGTH 4

LJFlopCounter::LJFlopCounter(double rc) : _rc2(rc * rc) {
	_totalCounts.clear();
}

void LJFlopCounter::initTraversal(const size_t numCells) {
	_currentCounts.clear();
}

void LJFlopCounter::endTraversal() {

	DomainDecompBase& domainDecomp =  global_simulation->domainDecomposition();
	domainDecomp.collCommInit(6);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_molDist);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_LJ);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_Charges);
	domainDecomp.collCommAppendDouble(_currentCounts.calc_Macro);
	domainDecomp.collCommAppendDouble(_currentCounts.vect_op);
	domainDecomp.collCommAppendDouble(_currentCounts.masked_elements);
	domainDecomp.collCommAllreduceSum();
	_currentCounts.calc_molDist = domainDecomp.collCommGetDouble();
	_currentCounts.calc_LJ = domainDecomp.collCommGetDouble();
	_currentCounts.calc_Charges = domainDecomp.collCommGetDouble();
	_currentCounts.calc_Macro = domainDecomp.collCommGetDouble();
	_currentCounts.vect_op = domainDecomp.collCommGetDouble();
	_currentCounts.masked_elements = domainDecomp.collCommGetDouble();
	domainDecomp.collCommFinalize();

	_totalCounts.addCounts(_currentCounts);

	const double cflMolDist = _currentCounts.calc_molDist * _flops_MolDist;
	const double cflCenterDist = (_currentCounts.calc_LJ + _currentCounts.calc_Charges) * _flops_CenterDist;
	const double cflLJKernel = _currentCounts.calc_LJ * _flops_LJKernel;
	const double cflChargesKernel = _currentCounts.calc_Charges * _flops_ChargesKernel;
	const double cflSum = (_currentCounts.calc_LJ + _currentCounts.calc_Charges) * _flops_ForcesSum;
	const double cflMacro = _currentCounts.calc_Macro * _flops_LJMacroValues;
	const double cflMacroSum = _currentCounts.calc_Macro * _flops_MacroSum;
	const double cflTotal = cflMolDist + cflCenterDist + cflLJKernel + cflChargesKernel + cflSum + cflMacro + cflMacroSum;

	Log::global_log->info()
			<< "FLOP counts in LJ force calculation for this iteration:"
			<< std::endl << " Molecule distance: " << cflMolDist
			<< " Center distance: " << cflCenterDist << " LJ Kernel: "
			<< cflLJKernel << " LJ Sum: " << cflSum << " Macroscopic values: "
			<< cflMacro << " Macroscopic value sum: " << cflMacroSum << std::endl
			<< "Current total FLOPS: " << cflTotal << std::endl;


	const double flMolDist = _totalCounts.calc_molDist * _flops_MolDist;
	const double flCenterDist = (_totalCounts.calc_LJ + _totalCounts.calc_Charges) * _flops_CenterDist;
	const double flLJKernel = _totalCounts.calc_LJ * _flops_LJKernel;
	const double flChargesKernel = _totalCounts.calc_Charges * _flops_ChargesKernel;
	const double flSum = (_totalCounts.calc_LJ + _totalCounts.calc_Charges) * _flops_ForcesSum;
	const double flMacro = _totalCounts.calc_Macro * _flops_LJMacroValues;
	const double flMacroSum = _totalCounts.calc_Macro * _flops_MacroSum;
	_totalFlopCount = flMolDist + flCenterDist + flLJKernel + flChargesKernel + flSum + flMacro + flMacroSum;

	Log::global_log->info()
			<< "Accumulated FLOP counts in force calculation for this iteration:"
			<< std::endl << " Molecule distance: " << cflMolDist
			<< " Center distance: " << cflCenterDist << " LJ Kernel: "
			<< cflLJKernel << " LJ Sum: " << cflSum << " Macroscopic values: "
			<< cflMacro << " Macroscopic value sum: " << cflMacroSum << std::endl
			<< "Accumulated total FLOPS: " << _totalFlopCount << std::endl;
	Log::global_log->info() << "(VLENGTH="<<VLENGTH<<", counting only LJ) #VectorOperations: " << _totalCounts.vect_op << " #el.s processed: " << (_totalCounts.vect_op * VLENGTH)
				<< " #el.s masked: " << (_totalCounts.masked_elements) << std::endl;
}

void LJFlopCounter::preprocessCell(ParticleCell & c) {
}

void LJFlopCounter::postprocessCell(ParticleCell & c) {
}

void LJFlopCounter::processCell(ParticleCell & c) {
  if(c.isHaloCell()) return; // don't count Halo cells

	const MoleculeList & molecules = c.getParticlePointers();
	if (molecules.size() > 1) {
		const MoleculeList::const_iterator end_i = --(molecules.end());
		const MoleculeList::const_iterator end_j = molecules.end();

		for (MoleculeList::const_iterator i = molecules.begin(); i != end_i;
				++i) {

			bool vcompute = false;
			int vpos = 0;
			int vcomputations = 0; // #vcomputations
			int mask = 0; // #masked_elements

			MoleculeList::const_iterator j = i;
			++j;
			for (; j != end_j; ++j) {

				// Have to compare the distance between 2 molecules.
				_currentCounts.calc_molDist += 1;

				const double d_x = (*i)->r(0) - (*j)->r(0);
				const double d_y = (*i)->r(1) - (*j)->r(1);
				const double d_z = (*i)->r(2) - (*j)->r(2);
				const double d2 = d_x * d_x + d_y * d_y + d_z * d_z;
				const size_t numLJcenters_j = (*j)->numLJcenters();
				if (d2 < _rc2) {
					const size_t numLJcenters_i = (*i)->numLJcenters();
					const size_t numCharges_i = (*i)->numCharges();
					const size_t numCharges_j = (*j)->numCharges();

					// Have to calculate the LJ force for each pair of centers.
					_currentCounts.calc_LJ += numLJcenters_i * numLJcenters_j;
					// Have to calculate the charge force for each pair of centers.
					_currentCounts.calc_Charges += numCharges_i * numCharges_j;
					// Have to calculate macroscopic values for each pair of centers.
					_currentCounts.calc_Macro += numLJcenters_i * numLJcenters_j + numCharges_i * numCharges_j;
				}

				for (size_t c = 0; c < numLJcenters_j; c++) {
					if (!vcompute && d2 < _rc2) {
						vcompute = true;
						vcomputations++;
						mask += vpos;
					} else if (vcompute && d2 > _rc2) {
						mask++;
					}

					vpos++;

					if (vpos == VLENGTH) {
						vpos = 0;
						vcompute = false;
					}
				}

			}
			vcomputations *= (*i)->numLJcenters();
			mask *= (*i)->numLJcenters();
			_currentCounts.vect_op += vcomputations;
			_currentCounts.masked_elements += mask;
		}
	}
}

void LJFlopCounter::processCellPair(ParticleCell & c1, ParticleCell & c2) {
  if(c1.isHaloCell() and c2.isHaloCell()) return; // don't count Halo cells

	const MoleculeList & molecules1 = c1.getParticlePointers();
	const MoleculeList & molecules2 = c2.getParticlePointers();
	if ((molecules1.size() > 0) && (molecules2.size() > 0)) {
		const MoleculeList::const_iterator end_i = molecules1.end();
		const MoleculeList::const_iterator end_j = molecules2.end();

		for (MoleculeList::const_iterator i = molecules1.begin(); i != end_i;
				++i) {

			bool vcompute = false;
			int vpos = 0;
			int vcomputations = 0; // #vcomputations
			int mask = 0; // #masked_elements

			for (MoleculeList::const_iterator j = molecules2.begin();
					j != end_j; ++j) {

				// Have to compare the distance between 2 molecules.
				_currentCounts.calc_molDist += 1;

				const double d_x = (*i)->r(0) - (*j)->r(0);
				const double d_y = (*i)->r(1) - (*j)->r(1);
				const double d_z = (*i)->r(2) - (*j)->r(2);
				const size_t numLJcenters_j = (*j)->numLJcenters();

				const double d2 = d_x * d_x + d_y * d_y + d_z * d_z;
				if (d2 < _rc2) {
					const size_t numLJcenters_i = (*i)->numLJcenters();
					const size_t numCharges_i = (*i)->numCharges();
					const size_t numCharges_j = (*j)->numCharges();

					// Have to calculate the LJ force for each pair of centers.
					_currentCounts.calc_LJ += numLJcenters_i * numLJcenters_j;
					// Have to calculate the charge force for each pair of centers.
					_currentCounts.calc_Charges += numCharges_i * numCharges_j;

					if ((c1.isHaloCell() == (!c2.isHaloCell()))
							&& ((*i)->isLessThan(**j))) {
						// Have to calculate macroscopic values for each pair of centers.
						_currentCounts.calc_Macro += numLJcenters_i * numLJcenters_j + numCharges_i * numCharges_j;
					}
				}

				for (size_t c = 0; c < numLJcenters_j; c++) {
					if (!vcompute && d2 < _rc2) {
						vcompute = true;
						vcomputations++;
						mask += vpos;
					} else if (vcompute && d2 > _rc2) {
						mask++;
					}

					vpos++;

					if (vpos == VLENGTH) {
						vpos = 0;
						vcompute = false;
					}
				}
			}

			vcomputations *= (*i)->numLJcenters();
			mask *= (*i)->numLJcenters();
			_currentCounts.vect_op += vcomputations;
			_currentCounts.masked_elements += mask;
		}
	}
}
