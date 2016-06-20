/*
 * UniformPseudoParticleContainer.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#include "UniformPseudoParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "bhfmm/utils/RotationParameterLookUp.h"
#include "particleContainer/ParticleContainer.h"

#include <algorithm>

namespace bhfmm {

#ifndef WIGNER
#define WIGNER 0 // 0: original, 1: Wigner rotation acceleration
#endif


#define IsOdd(x) ((x) & 1)
#define ToEven(x) ((x) & ~1)

UniformPseudoParticleContainer::UniformPseudoParticleContainer(
		double domainLength[3], double bBoxMin[3], double bBoxMax[3],
		double LJCellLength[3], unsigned LJSubdivisionFactor, int orderOfExpansions,
		bool periodic) :
		PseudoParticleContainer(orderOfExpansions), _leafContainer(0), _wellSep(1),
		_M2M_Wigner(4, WignerMatrix(orderOfExpansions, true)), _L2L_Wigner(4, WignerMatrix(orderOfExpansions, true)) {

	_periodicBC = periodic;

#ifdef ENABLE_MPI
	_timerProcessCells.set_sync(false);
	_timerAllreduce.set_sync(false);
	_timerAllreduce_me.set_sync(false);
	_timerCombineMpCell.set_sync(false);
	_timerGatherWellSepLo.set_sync(false);
	_timerPropagateCellLo.set_sync(false);
	_timerProcessFarField.set_sync(false);
#endif

	_leafContainer = new LeafNodesContainer(bBoxMin, bBoxMax, LJCellLength,
			LJSubdivisionFactor, periodic);

	double cellLength[3];

	for (int i = 0; i < 3; i++) {
		cellLength[i] = _leafContainer->getCellLength()[i];
		_cellLength[i] = cellLength[i];
	}

	_globalNumCellsPerDim = domainLength[0] / cellLength[0];
	_maxLevel = log2(_globalNumCellsPerDim);
	assert(_maxLevel == log2(domainLength[1] / cellLength[1]));
	assert(_maxLevel == log2(domainLength[2] / cellLength[2]));


	//allocate Multipole and Local particles
	int num_cells_in_level = 1;
	_mpCell.reserve(_maxLevel + 1);
	for (int n = 0; n <= _maxLevel; n++) {
		_mpCell.push_back(std::vector<MpCell>(num_cells_in_level, _maxOrd));
		num_cells_in_level *= 8;
	}
//	assert(
//			num_cells_in_level
//					== _globalNumCellsPerDim * _globalNumCellsPerDim
//							* _globalNumCellsPerDim);

	// initalize centers and radii
	num_cells_in_level = 1;
	int num_cells_in_level_one_dim = 1;
	Vector3<double> current_pos;
	Vector3<double> current_cell_length(domainLength);

	for (int n = 0; n <= _maxLevel; ++n) {
		for (int z = 0; z < num_cells_in_level_one_dim; ++z) {
			for (int y = 0; y < num_cells_in_level_one_dim; ++y) {
				for (int x = 0; x < num_cells_in_level_one_dim; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2];
					int cellIndex = ((z * num_cells_in_level_one_dim + y)
							* num_cells_in_level_one_dim) + x;
					_mpCell[n][cellIndex].multipole.setCenter(current_pos);
					_mpCell[n][cellIndex].multipole.setRadius(
							current_cell_length.L2Norm() * 0.5);

					_mpCell[n][cellIndex].local.setCenter(current_pos);
					_mpCell[n][cellIndex].local.setRadius(
							current_cell_length.L2Norm() * 0.5);
				}
			}
		}
		num_cells_in_level_one_dim *= 2;
		num_cells_in_level *= 8; //ToDo: is it needed anymore?
		current_cell_length = Vector3<double>(domainLength)
				* (1.0 / num_cells_in_level_one_dim);
	}

	_domain = global_simulation->getDomain();

	_globalNumCells = pow(_globalNumCellsPerDim, 3);
	_occVector = new int[_globalNumCells];
	std::fill(_occVector, _occVector + _globalNumCells, 0);

	_coeffVectorLength = 0;
//	_coeffVectorLength = _mpCell[0][0].multipole.get
	for (int j = 0; j <= _maxOrd; j++) {
		for (int k = 0; k <= j; k++) {
			_coeffVectorLength++;
		}
	}
	_coeffVectorLength *= _globalNumCells;
	_coeffVector = new double[_coeffVectorLength * 2];
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

	Log::global_log->info() << "UniformPseudoParticleContainer: coeffVectorLength="
			<< _coeffVectorLength << " Size of MPI Buffers is "
			<< (8 * (_coeffVectorLength * 2 + _globalNumCells)
					/ (1024.0 * 1024.0)) << " MB;" << std::endl;

	/* Initialize sin/cos(phi) lookup */
	double phi[4];
	phi[0] = -3*M_PI/4;
	phi[1] = -M_PI/4;
	phi[2] = 3*M_PI/4;
	phi[3] = M_PI/4;

	_CosSin = new double[(_maxOrd+1)*4];
	// We only use negative phi -> cos(-x)=cos(x), sin(-x)=-sin(x)
	for (int m = 0; m <= _maxOrd; ++m) {
		CosSin_ptr(0)[2*m] 		= cos(m*phi[0]);
		CosSin_ptr(0)[2*m + 1] 	= sin(m*phi[0]);
		CosSin_ptr(1)[2*m] 		= cos(m*phi[1]);
		CosSin_ptr(1)[2*m + 1] 	= sin(m*phi[1]);
	}

	/* Initialize M2M-Wigner matrices */
	double theta = atan(sqrt(2.0));
	_M2M_Wigner[0].evaluate(M_PI-theta);
	_M2M_Wigner[1].evaluate(theta);
	_M2M_Wigner[2].evaluate(-(M_PI-theta));
	_M2M_Wigner[3].evaluate(-theta);
	/* Initialize Lookup table */
	RotationParameterLookUp::tab = new RotationParameterLookUp(_maxOrd);
	RotationParameterLookUp::tab->initFromDirEval();

	// pre-multiply prefactors
	for (unsigned i = 0; i<_M2M_Wigner.size(); ++i) {
		WignerMatrix &W = _M2M_Wigner[i];
		for (int l = 0; l <= _maxOrd; ++l) {
			for (int m = 0; m <= l; ++m) {
				for (int k = -l; k<=l; ++k) {
					const double factor = RotationParameterLookUp::tab->acc_c(l,m,k);
					W.acc(l,m,k) *= factor ;
				}
			}
		}

	}

	_L2L_Wigner[0].evaluate(M_PI-theta);
	_L2L_Wigner[1].evaluate(theta);
	_L2L_Wigner[2].evaluate(-(M_PI-theta));
	_L2L_Wigner[3].evaluate(-theta);

	// pre-multiply prefactors
	for (unsigned i = 0; i<_L2L_Wigner.size(); ++i) {
		WignerMatrix &W = _L2L_Wigner[i];
		for (int l = 0; l <= _maxOrd; ++l) {
			for (int m = 0; m <= l; ++m) {
				for (int k = -l; k<=l; ++k) {
					const double factor = RotationParameterLookUp::tab->acc_c(l,k,m);
					W.acc(l,m,k) *= factor ;
				}
			}
		}

	}

	/* Initialize M2L-Wigner matrices */
	processTreeInitM2LWigner();

	delete RotationParameterLookUp::tab;

#ifdef FMM_FFT
	//FFT objects initialization
	FFTSettings::autoSetting(_maxOrd); //auto set FFT acceleration's settings to optimal values
	//FFTSettings::TFMANAGER_VERBOSE = false;
	FFTSettings::printCurrentOptions();
	if (FFTSettings::issetFFTAcceleration()) {
		_FFTAcceleration = FFTFactory::getFFTAccelerationAPI(_maxOrd);
		_FFT_TM = FFTFactory::getTransferFunctionManagerAPI(_maxOrd,
				_FFTAcceleration);
	} else {
		_FFTAcceleration = NULL;
		_FFT_TM = NULL;
	}
#endif  /* FMM_FFT */

}

UniformPseudoParticleContainer::~UniformPseudoParticleContainer() {
	delete _leafContainer;
	delete[] _coeffVector;
	delete[] _occVector;
	delete[] _CosSin;
	
#ifdef FMM_FFT
	if (FFTSettings::issetFFTAcceleration()) {
		delete _FFT_TM;
		delete _FFTAcceleration;
	}
#endif  /* FMM_FFT */
}

void UniformPseudoParticleContainer::build(ParticleContainer* pc) {
	_leafContainer->clearParticles();

	Molecule* tM;
	for(tM = pc->begin(); tM != pc->end(); tM = pc->next()) {
		_leafContainer->addParticle(*tM);
	}
}

void UniformPseudoParticleContainer::upwardPass(P2MCellProcessor* cp) {
	// P2M
	_leafContainer->traverseCellPairs(*cp);

	// M2M
	AllReduceMultipoleMoments();

	_timerCombineMpCell.start();

	int curCellsEdge=_globalNumCellsPerDim;
	double cellWid[3];

	for(int i=0; i <3; i++)	cellWid[i]=_cellLength[i];

	// when considering periodic boundary conditions, there is actually work up to level 1!
	for(int curLevel=_maxLevel-1; curLevel>=1; curLevel--){
		curCellsEdge /=2;

		for(int i=0; i <3; i++)	cellWid[i] *=2;
#if WIGNER==0
		CombineMpCell(cellWid, curCellsEdge, curLevel);
#else
		CombineMpCell_Wigner(cellWid, curCellsEdge, curLevel);
#endif
	}
	_timerCombineMpCell.stop();
}

void UniformPseudoParticleContainer::horizontalPass(
		VectorizedChargeP2PCellProcessor* cp) {
	// P2P
	_leafContainer->traverseCellPairs(*cp);

	// M2L
	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}

#ifndef FMM_FFT
#if WIGNER==0
		GatherWellSepLo(cellWid, curCellsEdge, curLevel);
#else
		GatherWellSepLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
#else
		GatherWellSepLo_FFT(cellWid, curCellsEdge, curLevel);
#endif  /* FMM_FFT */
		AllReduceLocalMoments(curCellsEdge, curLevel);
	}
}

void UniformPseudoParticleContainer::downwardPass(L2PCellProcessor* cp) {
	// L2L
	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /= 2;
		}

#if WIGNER==0
		PropagateCellLo(cellWid, curCellsEdge, curLevel);
#else
		PropagateCellLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
	}

	// L2P
	_leafContainer->traverseCellPairs(*cp);
}



void UniformPseudoParticleContainer::CombineMpCell(double */*cellWid*/, int& mpCells, int& curLevel){
	int iDir, m1=0, m1x, m1y, m1z, m2=0;
	int m2v[3] = {0, 0, 0};
	int mpCellsN=2*mpCells;


	for(m1z=0; m1z<mpCells; m1z++){
		for(m1y=0; m1y<mpCells; m1y++){
			for(m1x=0; m1x<mpCells; m1x++){
				m1=(m1z*mpCells + m1y)*mpCells + m1x;

				for(iDir=0; iDir<8; iDir++){

					m2v[0]=2*m1x;
					m2v[1]=2*m1y;
					m2v[2]=2*m1z;

					if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
					if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
					if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


					m2=(m2v[2]*mpCellsN + m2v[1])*mpCellsN + m2v[0];

					if(_mpCell[curLevel+1][m2].occ==0) continue;

					_mpCell[curLevel][m1].occ +=_mpCell[curLevel+1][m2].occ;

					_mpCell[curLevel][m1].multipole.addMultipoleParticle(_mpCell[curLevel+1][m2].multipole);
				} // iDir closed
			}// m1x closed
		}// m1y closed
	} // m1z closed
}

void UniformPseudoParticleContainer::CombineMpCell_Wigner(double *cellWid, int& mpCells, int& curLevel){
	int iDir, m1=0, m1x, m1y, m1z, m2=0;
	int m2v[3] = {0, 0, 0};
	int mpCellsN=2*mpCells;

	const double magnitude = sqrt(3.0)/4.0 *cellWid[0];

	for(m1z=0; m1z<mpCells; m1z++){
		for(m1y=0; m1y<mpCells; m1y++){
			for(m1x=0; m1x<mpCells; m1x++){
				m1=(m1z*mpCells + m1y)*mpCells + m1x;

				for(iDir=0; iDir<8; iDir++){

					m2v[0]=2*m1x;
					m2v[1]=2*m1y;
					m2v[2]=2*m1z;

					if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
					if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
					if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


					m2=(m2v[2]*mpCellsN + m2v[1])*mpCellsN + m2v[0];

					if(_mpCell[curLevel+1][m2].occ==0) continue;

					_mpCell[curLevel][m1].occ +=_mpCell[curLevel+1][m2].occ;

					/* get Wigner matrix */
					const int idxPhi = iDir%2;
					const int negatePhi = 1 - ((iDir%4)/2)*2;
					const int idxWig = iDir/4;

					_mpCell[curLevel][m1].multipole.addMultipoleParticle_Wigner(_mpCell[curLevel+1][m2].multipole, M2M_Wigner(idxWig),
							M2M_Wigner(idxWig+2), CosSin_ptr(idxPhi), negatePhi, magnitude);
				} // iDir closed
			}// m1x closed
		}// m1y closed
	} // m1z closed
}


#define HiLim(t) ToEven(m1v[t])+ 2*_wellSep+1
#define LoLim(t) ToEven(m1v[t])- 2*_wellSep

void UniformPseudoParticleContainer::GatherWellSepLo(double *cellWid, int mpCells, int& curLevel){
	_timerGatherWellSepLo.start();

	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			if (_periodicBC == false and (m2z < 0 or m2z >= mpCells)) {
				continue;
			}

			// to get periodic image
			m22z = (mpCells + m2z) % mpCells;
			periodicShift[2] = 0.0;
			if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
			if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				if (_periodicBC == false and (m2y < 0 or m2y >= mpCells)) {
					continue;
				}

				// to get periodic image
				m22y = (mpCells + m2y) % mpCells;

				periodicShift[1] = 0.0;
				if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
				if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					if (_periodicBC == false and (m2x < 0 or m2x >= mpCells)) {
						continue;
					}

					// to get periodic image
					m22x = (mpCells + m2x) % mpCells;

					periodicShift[0] = 0.0;
					if (m2x < 0)		periodicShift[0] = -mpCells * cellWid[0];
					if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
					//
					m2v[0] = m2x;

					if (abs(m2v[0] - m1v[0]) <= _wellSep &&
						abs(m2v[1] - m1v[1]) <= _wellSep &&
						abs(m2v[2] - m1v[2]) <= _wellSep)
						continue;
					m2 = (m22z * mpCells + m22y) * mpCells + m22x;

					if (_mpCell[curLevel][m2].occ == 0)
						continue;

					_mpCell[curLevel][m1].local.addMultipoleParticle(
							_mpCell[curLevel][m2].multipole, periodicShift);
				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed

	_timerGatherWellSepLo.stop();

} // GatherWellSepLo closed

void UniformPseudoParticleContainer::GatherWellSepLo_Wigner(double *cellWid, int mpCells, int& curLevel){
	_timerGatherWellSepLo.start();

	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			if (_periodicBC == false and (m2z < 0 or m2z >= mpCells)) {
				continue;
			}

			// to get periodic image
			m22z = (mpCells + m2z) % mpCells;
			periodicShift[2] = 0.0;
			if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
			if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				if (_periodicBC == false and (m2y < 0 or m2y >= mpCells)) {
					continue;
				}

				// to get periodic image
				m22y = (mpCells + m2y) % mpCells;

				periodicShift[1] = 0.0;
				if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
				if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					if (_periodicBC == false and (m2x < 0 or m2x >= mpCells)) {
						continue;
					}

					// to get periodic image
					m22x = (mpCells + m2x) % mpCells;

					periodicShift[0] = 0.0;
					if (m2x < 0)		periodicShift[0] = -mpCells * cellWid[0];
					if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
					//
					m2v[0] = m2x;

					if (abs(m2v[0] - m1v[0]) <= _wellSep &&
						abs(m2v[1] - m1v[1]) <= _wellSep &&
						abs(m2v[2] - m1v[2]) <= _wellSep)
						continue;
					m2 = (m22z * mpCells + m22y) * mpCells + m22x;

					if (_mpCell[curLevel][m2].occ == 0)
						continue;

					_mpCell[curLevel][m1].local.addMultipoleParticle_Wigner(
							_mpCell[curLevel][m2].multipole, periodicShift, cellWid, _M2L_Wigner);
				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed

	_timerGatherWellSepLo.stop();

}

#ifdef FMM_FFT
void UniformPseudoParticleContainer::GatherWellSepLo_FFT(double *cellWid, int mpCells, int& curLevel) {
	if (FFTSettings::USE_VECTORIZATION) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	} else {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	}
}
  

template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_template(
		double *cellWid, int mpCells, int& curLevel) {
	_timerGatherWellSepLo.start();

	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length = mpCells * mpCells * mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank = domainDecomp.getRank();
	_size = domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	//FFT param
	double radius;
	double base_unit = 2.0 / sqrt(3);
	int M2L_order;
	FFTDataContainer* tf;

	//Initialize FFT
	for (m1 = loop_min; m1 < loop_max; m1++) {
		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		radius = _mpCell[curLevel][m1].local.getRadius();
		FFTAccelerableExpansion& source =
				static_cast<bhfmm::SHMultipoleParticle&>(_mpCell[curLevel][m1].multipole).getExpansion();
		FFTAccelerableExpansion& target =
				static_cast<bhfmm::SHLocalParticle&>(_mpCell[curLevel][m1].local).getExpansion();
		_FFTAcceleration->FFT_initialize_Source(source, radius);
		_FFTAcceleration->FFT_initialize_Target(target);

	}

	//M2L in Fourier space
	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			// to get periodic image
			m22z = (mpCells + m2z) % mpCells;
			periodicShift[2] = 0.0;
			if (m2z < 0)
				periodicShift[2] = -mpCells * cellWid[2];
			if (m2z >= mpCells)
				periodicShift[2] = mpCells * cellWid[2];

			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				// to get periodic image
				m22y = (mpCells + m2y) % mpCells;

				periodicShift[1] = 0.0;
				if (m2y < 0)
					periodicShift[1] = -mpCells * cellWid[1];
				if (m2y >= mpCells)
					periodicShift[1] = mpCells * cellWid[1];

				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					// to get periodic image
					m22x = (mpCells + m2x) % mpCells;

					periodicShift[0] = 0.0;
					if (m2x < 0)
						periodicShift[0] = -mpCells * cellWid[0];
					if (m2x >= mpCells)
						periodicShift[0] = mpCells * cellWid[0];
					//
					m2v[0] = m2x;

					if (abs(m2v[0] - m1v[0]) <= _wellSep
							&& abs(m2v[1] - m1v[1]) <= _wellSep
							&& abs(m2v[2] - m1v[2]) <= _wellSep)
						continue;
					m2 = (m22z * mpCells + m22y) * mpCells + m22x;

					if (_mpCell[curLevel][m2].occ == 0)
						continue;

					if (UseM2L_2way) {
						if (m1 > m2)
							continue;
					}

					tf = _FFT_TM->getTransferFunction(m2x - m1v[0],
							m2y - m1v[1], m2z - m1v[2], base_unit, base_unit,
							base_unit);

					if (UseM2L_2way) {
						FFTAccelerableExpansion& source2 =
								static_cast<bhfmm::SHMultipoleParticle&>(_mpCell[curLevel][m2].multipole).getExpansion();
						FFTAccelerableExpansion& target2 =
								static_cast<bhfmm::SHLocalParticle&>(_mpCell[curLevel][m2].local).getExpansion();
						FFTAccelerableExpansion& source1 =
								static_cast<bhfmm::SHMultipoleParticle&>(_mpCell[curLevel][m1].multipole).getExpansion();
						FFTAccelerableExpansion& target1 =
								static_cast<bhfmm::SHLocalParticle&>(_mpCell[curLevel][m1].local).getExpansion();

						if (UseOrderReduction) {
							M2L_order = FFTOrderReduction::getM2LOrder(
									m2v[0] - m1v[0], m2v[1] - m1v[1],
									m2v[2] - m1v[2], _maxOrd);
							if (UseVectorization) {
								static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_2way_ORed_vec(
										source2, source1, target2, target1, tf,
										M2L_order);
							} else {
								static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_2way_ORed(
										source2, source1, target2, target1, tf,
										M2L_order);
							}
						} else {
							if (UseVectorization) {
								static_cast<FFTAccelerationAPI_2Way*>(_FFTAcceleration)->FFT_M2L_2way_vec(
										source2, source1, target2, target1, tf);
							} else {
								static_cast<FFTAccelerationAPI_2Way*>(_FFTAcceleration)->FFT_M2L_2way(
										source2, source1, target2, target1, tf);
							}
						}
					} else {
						FFTAccelerableExpansion& source =
								static_cast<bhfmm::SHMultipoleParticle&>(_mpCell[curLevel][m2].multipole).getExpansion();
						FFTAccelerableExpansion& target =
								static_cast<bhfmm::SHLocalParticle&>(_mpCell[curLevel][m1].local).getExpansion();

						if (UseOrderReduction) {
							M2L_order = FFTOrderReduction::getM2LOrder(
									m2v[0] - m1v[0], m2v[1] - m1v[1],
									m2v[2] - m1v[2], _maxOrd);
							if (UseVectorization) {
								static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
										source, target, tf, M2L_order);
							} else {
								static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction(
										source, target, tf, M2L_order);
							}
						} else {
							if (UseVectorization) {
								_FFTAcceleration->FFT_M2L_vec(source, target,
										tf);
							} else {
								_FFTAcceleration->FFT_M2L(source, target, tf);
							}
						}
					}

					if (!UseTFMemoization) {
						delete tf; //free useless memory
					}
					//_mpCell[curLevel][m1].local.addMultipoleParticle(_mpCell[curLevel][m2].multipole, periodicShift);
				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed

	//Finalize FFT
	for (m1 = loop_min; m1 < loop_max; m1++) {
		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		radius = _mpCell[curLevel][m1].local.getRadius();
		FFTAccelerableExpansion& target =
				static_cast<bhfmm::SHLocalParticle&>(_mpCell[curLevel][m1].local).getExpansion();
		_FFTAcceleration->FFT_finalize_Target(target, radius);
	}

	_timerGatherWellSepLo.stop();

} // GatherWellSepLo_FFT_template closed
#endif /* FMM_FFT */

void UniformPseudoParticleContainer::PropagateCellLo(double */*cellWid*/, int mpCells, int& curLevel){
	_timerPropagateCellLo.start();
	int m1v[3];
	int m2v[3];

	int iDir, m1, m1x, m1y, m1z, m2;

	int mpCellsN = 2*mpCells;

// TODO: parallelization is broken currently, but parallelizing L2L is not all that important
//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();


//	int _rank= domainDecomp.getRank();
//	int _size= domainDecomp.getNumProcs();

//	int loop_min = (int) ((long) (_rank + 0) * (long) (mpCells * mpCells * mpCells) / (long) _size);
//	int loop_max = (int) ((long) (_rank + 1) * (long) (mpCells * mpCells * mpCells) / (long) _size);
	int loop_min = 0;
	int loop_max = mpCells * mpCells * mpCells;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		m1x = m1v[0];
		m1y = m1v[1];
		m1z = m1v[2];

		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (iDir = 0; iDir < 8; iDir++) {
			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			_mpCell[curLevel][m1].local.actOnLocalParticle(
					_mpCell[curLevel + 1][m2].local);
		} // iDir
	}
	_timerPropagateCellLo.stop();
} // PropogateCellLo

void UniformPseudoParticleContainer::PropagateCellLo_Wigner(double *cellWid, int mpCells, int& curLevel){
	_timerPropagateCellLo.start();
	int m1v[3];
	int m2v[3];

	int iDir, m1, m1x, m1y, m1z, m2;

	int mpCellsN = 2*mpCells;
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int _rank= domainDecomp.getRank();
	int _size= domainDecomp.getNumProcs();

	int loop_min = (int) ((long) (_rank + 0) * (long) (mpCells * mpCells * mpCells) / (long) _size);
	int loop_max = (int) ((long) (_rank + 1) * (long) (mpCells * mpCells * mpCells) / (long) _size);

	const double magnitude = sqrt(3.0)/4.0 *cellWid[0];

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		m1x = m1v[0];
		m1y = m1v[1];
		m1z = m1v[2];

		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (iDir = 0; iDir < 8; iDir++) {
			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			/* get Wigner matrix */
			const int idxPhi = iDir%2;
			const int negatePhi = 1 - ((iDir%4)/2)*2;
			const int idxWig = iDir/4;

			_mpCell[curLevel][m1].local.actOnLocalParticle_Wigner(
					_mpCell[curLevel + 1][m2].local, L2L_Wigner(idxWig),
							L2L_Wigner(idxWig+2), CosSin_ptr(idxPhi), negatePhi, magnitude);

		} // iDir
	}
	_timerPropagateCellLo.stop();
} // PropogateCellLo


void UniformPseudoParticleContainer::processMultipole(ParticleCell& cell){
	int cellIndexV[3];
	for (int i = 0; i < 3; i++) {
		cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
	}

	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1])	* _globalNumCellsPerDim) + cellIndexV[0];

//	assert(cell.isInActiveWindow());

	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();


	int Occupied = 0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		++Occupied;
		Molecule& molecule1 = *currentCellParticles[i];
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const std::array<double,3> dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			double dr[3];

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}	// for k closed

			bhfmm::Vector3<double> site_pos_vec3(dr);
			_mpCell[_maxLevel][cellIndex].multipole.addSource(site_pos_vec3, chargei.q());

		}// for j closed
	} // current particle closed

	_mpCell[_maxLevel][cellIndex].occ = Occupied;
}

void UniformPseudoParticleContainer::processFarField(ParticleCell& cell) {
	int cellIndexV[3];
	for (int i = 0; i < 3; i++) {
		cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
	}

	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];

	bhfmm::SolidHarmonicsExpansion leLocal(_maxOrd);
	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();
	double u = 0;
	double uSum = 0.0;
	double f[3] = {0.0, 0.0, 0.0};
	Vector3<double>f_vec3;
	double virialSum=0.0;
	double P_xxSum=0.0;
	double P_yySum=0.0;
	double P_zzSum=0.0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		Molecule& molecule1 = *currentCellParticles[i];
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const std::array<double,3> dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			Vector3<double> dr;

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}       // for k closed

			_mpCell[_maxLevel][cellIndex].local.actOnTarget(dr,chargei.q(),u,f_vec3);
			f[0] = f_vec3[0];
			f[1] = f_vec3[1];
			f[2] = f_vec3[2];

			double virial = 0.0;
			for(int l=0; l<3; l++){
				virial +=-f[l]*dr[l];
			}
			P_xxSum +=0.5*-f[0]*dr[0];
			P_yySum +=0.5*-f[1]*dr[1];
			P_zzSum +=0.5*-f[2]*dr[2];
			molecule1.Fchargeadd(j, f);
			uSum +=0.5*chargei.q()*u;
			virialSum +=0.5*virial;
		}// for j closed
	} // current particle closed

//	_domain->addLocalUpot(uSum);
//	_domain->addLocalVirial(virialSum);
//	_domain->addLocalP_xx(P_xxSum);
//	_domain->addLocalP_yy(P_yySum);
//	_domain->addLocalP_zz(P_zzSum);
}

void UniformPseudoParticleContainer::clear() {
	for (int n = _maxLevel; n >= 1; n--) {
		int mpCells = pow(2, n);

		for (int m1z = 0; m1z < mpCells; m1z++) {
			for (int m1y = 0; m1y < mpCells; m1y++) {
				for (int m1x = 0; m1x < mpCells; m1x++) {
					int cellIndexNew = (m1z * mpCells + m1y) * mpCells + m1x;

					_mpCell[n][cellIndexNew].occ = 0;

					_mpCell[n][cellIndexNew].multipole.clear();
					_mpCell[n][cellIndexNew].local.clear();
				}
			}
		}
	}

	// clear the MPI buffers
#ifdef ENABLE_MPI
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength*2, 0.0);
	std::fill(_occVector, _occVector + _globalNumCells, 0);
#endif
}

void UniformPseudoParticleContainer::AllReduceMultipoleMoments() {
	_timerAllreduce.start();
#ifdef ENABLE_MPI

	int coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {
		const MpCell & currentCell = _mpCell[_maxLevel][cellIndex];

		// NOTE: coeffIndex modified in following call:
		currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

		assert(cellIndex < _globalNumCells);
		_occVector[cellIndex] = currentCell.occ;
	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _occVector, _globalNumCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {

		MpCell & currentCell = _mpCell[_maxLevel][cellIndex];

		currentCell.occ = _occVector[cellIndex];
		currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce.stop();
}

void UniformPseudoParticleContainer::AllReduceLocalMoments(int mpCells, int _curLevel) {
	_timerAllreduce_me.start();

#ifdef ENABLE_MPI

	const int _row_Length=pow(mpCells, 3);
	int coeffIndex = 0;

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		const MpCell & currentCell = _mpCell[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, coeffIndex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		MpCell & currentCell = _mpCell[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce_me.stop();
}



void UniformPseudoParticleContainer::processTree() {

	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}
#if WIGNER==0
		GatherWellSepLo(cellWid, curCellsEdge, curLevel);
#else
		GatherWellSepLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
		AllReduceLocalMoments(curCellsEdge, curLevel);

		if(curLevel<_maxLevel) {
#if WIGNER==0
			PropagateCellLo(cellWid, curCellsEdge, curLevel);
#else
			PropagateCellLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
		}
	}
}

void UniformPseudoParticleContainer::processTreeInitM2LWigner() {


	int curCellsEdge = 1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}
		GatherWellSepLoInitM2LWigner(cellWid, curCellsEdge, curLevel);

	}
}

void UniformPseudoParticleContainer::GatherWellSepLoInitM2LWigner(double *cellWid, int mpCells, int& curLevel) {
	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
//		if (_mpCell[curLevel][m1].occ == 0)
//			continue;

		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			// to get periodic image
			m22z = (mpCells + m2z) % mpCells;
			periodicShift[2] = 0.0;
			if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
			if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				// to get periodic image
				m22y = (mpCells + m2y) % mpCells;

				periodicShift[1] = 0.0;
				if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
				if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					// to get periodic image
					m22x = (mpCells + m2x) % mpCells;

					periodicShift[0] = 0.0;
					if (m2x < 0)		periodicShift[0] = -mpCells * cellWid[0];
					if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
					//
					m2v[0] = m2x;

					if (abs(m2v[0] - m1v[0]) <= _wellSep &&
						abs(m2v[1] - m1v[1]) <= _wellSep &&
						abs(m2v[2] - m1v[2]) <= _wellSep)
						continue;
					m2 = (m22z * mpCells + m22y) * mpCells + m22x;

					const SHMultipoleParticle& sh_multipole = _mpCell[curLevel][m2].multipole;
					const SHLocalParticle& sh_local = _mpCell[curLevel][m1].local;

					// compute periodically-shifted center
					Vector3<double> shifted_center = sh_multipole.getCenter() + periodicShift;

					// distance-vector FROM local TO multipole
					Vector3<double> r_target_to_source = shifted_center - sh_local.getCenter();

					/* compute index vector and add param to map */
					Vector3<int> idxVec(rint(r_target_to_source[0]/(cellWid[0])),
							rint(r_target_to_source[1]/(cellWid[1])),
							rint(r_target_to_source[2]/(cellWid[2])));

					// continue if element is already in container
					if (_M2L_Wigner.find(idxVec) != _M2L_Wigner.end() ) continue;

					double projXY_len = sqrt(
							r_target_to_source[0] * r_target_to_source[0]
									+ r_target_to_source[1] * r_target_to_source[1]);

					double phi = atan2(r_target_to_source[1], r_target_to_source[0]);
					double theta = atan2(projXY_len, r_target_to_source[2]);

					double* CosSin = new double[(_maxOrd+1)*2];
					for (int m = 0; m <= _maxOrd; ++m) {
						CosSin[2*m] 		= cos(m*phi);
						CosSin[2*m + 1] 	= sin(m*phi);
					}

					WignerMatrix W_pos(_maxOrd, true);
					WignerMatrix W_neg(_maxOrd, true);

					W_pos.evaluate(theta);
					W_neg.evaluate(-theta);

					// pre-multiply prefactors
					for (int l = 0; l <= _maxOrd; ++l) {
						for (int m = 0; m <= l; ++m) {
							for (int k = -l; k<=l; ++k) {
								W_pos.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,m,k);
								W_neg.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,k,m);
							}
						}
					}

					RotationParams& param = _M2L_Wigner[idxVec];// {W_pos, W_neg, CosSin};
					param.W[0] = W_pos;
					param.W[1] = W_neg;
					param.SinCos = CosSin;

				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed


}

void UniformPseudoParticleContainer::printTimers() {
	std::cout << "\t\t" << _timerAllreduce.get_etime()       		<< "\t\t" <<"s in Allreduce" << std::endl;
	std::cout << "\t\t" << _timerAllreduce_me.get_etime()			<< "\t\t" <<"s in Allreduce_me"<<std::endl;
	std::cout << "\t\t" << _timerCombineMpCell.get_etime()     		<< "\t\t" <<"s in CombineMpCell" << std::endl;
	std::cout << "\t\t" << _timerGatherWellSepLo.get_etime() 		<< "\t\t" <<"s in GatherWellSepLo" << std::endl;
	std::cout << "\t\t" << _timerPropagateCellLo.get_etime() 		<< "\t\t" <<"s in PropagateCellLo" << std::endl;
}


} /* namespace bhfmm */

