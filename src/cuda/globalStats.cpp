/*
 * globalStats.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: andreas
 */

#include "globalStats.h"

void GlobalStats::preInteractionCalculation() {
	_cellStats.resize( _linkedCells.getCells().size() );
	_cellStats.zeroDevice();

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
	_cellStatsLocks.resize( _linkedCells.getCells().size() );
	_cellStatsLocks.zeroDevice();
#endif
}

void GlobalStats::postInteractionCalculation() {
	const std::vector<CellStatsStorage> &cellStats = _cellStats.copyToHost();

	const std::vector<unsigned long> &innerCellIndices = _linkedCells.getInnerCellIndices();
	const std::vector<unsigned long> &boundaryCellIndices = _linkedCells.getBoundaryCellIndices();
	const std::vector<unsigned long> &haloCellIndices = _linkedCells.getHaloCellIndices();

	_potential = 0.0f;
	_virial = 0.0f;
	for( int i = 0 ; i < innerCellIndices.size() ; i++ ) {
		int innerCellIndex = innerCellIndices[i];
		_potential += cellStats[ innerCellIndex ].potential;
		_virial += cellStats[ innerCellIndex ].virial;
	}
	for( int i = 0 ; i < boundaryCellIndices.size() ; i++ ) {
		int boundaryCellIndex = boundaryCellIndices[ i ];
		_potential += cellStats[ boundaryCellIndex ].potential;
		_virial += cellStats[ boundaryCellIndex ].virial;
	}

	// to be used in conjunction with a different halo potential treatment
#if 0
	for( int i = 0 ; i < haloCellIndices.size() ; i++ ) {
		int haloCellIndex = haloCellIndices[ i ];
		if( _linkedCells.getCells()[haloCellIndex].haloOwnsInteractions() ) {
			_potential += cellStats[ haloCellIndex ].potential;
			_virial += cellStats[ haloCellIndex ].virial;
		}
		else {
			_potential -= cellStats[ haloCellIndex ].potential;
			_virial -= cellStats[ haloCellIndex ].virial;
		}
	}
#endif

	// every contribution is added twice so divide by 2
	_potential /= 2.0f;
	_virial /= 2.0f;

	// the CPU code has moleculeAtoB = A - B (which is obviously semantically wrong) and I use B - A
	// this only affects the virial sign
	_virial = -_virial;
}
