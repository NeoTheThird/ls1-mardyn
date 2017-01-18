/*
 * UniformPseudoParticleContainer.h
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#ifndef UNIFORMPSEUDOPARTICLECONTAINER_H_
#define UNIFORMPSEUDOPARTICLECONTAINER_H_

#include "PseudoParticleContainer.h"
#include "LeafNodesContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Timer.h"
#include "bhfmm/utils/WignerMatrix.h"
#include "bhfmm/utils/RotationParameter.h"

#ifdef FMM_FFT
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/FFTAccelerationAPI_extensions.h"
#include "bhfmm/fft/FFTSettings.h"
#include "bhfmm/fft/FFTFactory.h"
#include "bhfmm/fft/FFTOrderReduction.h"
#include "bhfmm/fft/TransferFunctionManagerAPI.h"
#endif /* FMM_FFT */

#include <vector>
#include <map>
#include "bhfmm/HaloBufferNoOverlap.h"
#include "bhfmm/HaloBufferOverlap.h"

class Domain;
class DomainDecompBase;

namespace bhfmm {

class UniformPseudoParticleContainer: public PseudoParticleContainer {
public:
	UniformPseudoParticleContainer(double domainLength[3], double bBoxMin[3], double bBoxMax[3],
			double LJCellLength[3], unsigned LJSubdivisionFactor, int orderOfExpansions,bool periodic = true);
	~UniformPseudoParticleContainer();

	void clear();
	void build(ParticleContainer* pc);
	void upwardPass(P2MCellProcessor * cp);
	void horizontalPass(VectorizedChargeP2PCellProcessor * cp);
	void downwardPass(L2PCellProcessor *cp);

	// P2M
	void processMultipole(ParticleCellPointers& cell);

	// L2P
	void processFarField(ParticleCellPointers& cell);

	// M2M, M2L, L2L
	void processTree();

	void printTimers();



private:
	LeafNodesContainer* _leafContainer;

	int _wellSep;
	int _maxLevel;
	int _globalLevel;
	//In the parallel version the octree is divided into two trees:
	//- A local subtree starting at _globalLevel + 1 which contains only the
	//ancestor nodes from the node that was assigned to the MPI process on
	//the _globalLevel
	//- A global tree which contains all the octree elements of the FMM tree
	//from level 0 to _globalLevel
	std::vector<std::vector<MpCell> > _mpCellGlobalTop;
	std::vector<std::vector<MpCell> > _mpCellLocal;
	double _cellLength[3];
	int _globalNumCellsPerDim;
	Domain* _domain;
	int _globalNumCells;
	int* _occVector; // array for MPI allgather

	int _coeffVectorLength;
	double* _coeffVector; // array for MPI allgather
	double* _coeffVector_me;
#ifdef ENABLE_MPI
	HaloBufferNoOverlap<double> * _multipoleRecBuffer, *_multipoleBuffer;
	HaloBufferOverlap<double> * _multipoleRecBufferOverlap, *_multipoleBufferOverlap;
	MPI_Request _allReduceRequest;

#endif
	bool _periodicBC;


	Vector3<int> _numProcessorsPerDim;
	Vector3<int> _numCellsOnGlobalLevel;
	Vector3<int> _processorPositionGlobalLevel;
	Vector3<double> _bBoxMin;
	std::vector<int> _neighbours;

	int _globalLevelNumCells;

#ifdef FMM_FFT
	TransferFunctionManagerAPI* _FFT_TM;
	FFTAccelerationAPI* _FFTAcceleration;
#endif  /* FMM_FFT */


	// M2M
	void CombineMpCell(double *cellWid, int& mpCells, int curLevel);

	// M2M
	void CombineMpCell_MPI(double *cellWid, Vector3<int> localMpCells, int curLevel, Vector3<int> offset);

	// M2L
	void GatherWellSepLo(double *cellWid, int mpCells, int curLevel);

	// M2L
	void GatherWellSepLo_MPI(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);
  

#ifdef FMM_FFT
	// M2L
	void GatherWellSepLo_FFT(double *cellWid, int mpCells, int& curLevel);

	void GatherWellSepLo_FFT_MPI(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);

	template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
	void GatherWellSepLo_FFT_template(double *cellWid, int mpCells, int& curLevel);

	template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
	void GatherWellSepLo_FFT_MPI_template(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);
#endif /* FMM_FFT */

	// L2L
	void PropagateCellLo(double *cellWid, int mpCells, int curLevel);

	// L2L
	void PropagateCellLo_MPI(double *cellWid, Vector3<int> localMpCells, int curLevel, Vector3<int> offset);

	// for parallelization
	void AllReduceMultipoleMoments();
	void AllReduceLocalMoments(int mpCells, int _curLevel);
	void AllReduceMultipoleMomentsLevelToTop(int mpCells, int _curLevel);
	void AllReduceMultipoleMomentsSetValues(int mpCells, int _curLevel);


	void getHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, double *buffer,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh);
	void setHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, double *bufferRec,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh);
	//for parallelization
	void communicateHalosNoOverlap();
	void communicateHalosOverlapStart();
	//has to be called after receive finished
	void communicateHalosOverlapSetHalos();
	void communicateHalos();
	void communicateHalosX();
	void communicateHalosY();
	void communicateHalosZ();
	void communicateHalosAlongAxis(double * lowerNeighbourBuffer, double * higherNeighbourBuffer,
			double * lowerNeighbourBufferRec, double * higherNeighbourBufferRec,
			int lowerNeighbour, int higherNeighbour, int haloSize
			);
	int busyWaiting();
	void initBusyWaiting(){
		_allReduceProcessed = 0;
		_halosProcessed = 0;
		_sendProcessed = 0;
	}
	int _allReduceProcessed;
	int _halosProcessed;
	int _sendProcessed;
	Timer _timerProcessCells;
	Timer _timerAllreduce;
	Timer _timerGatherWellSepLoGlobal;
	Timer _timerPropagateCellLoGlobal;
	Timer _timerCombineMpCellGlobal;
	Timer _timerCombineMpCellLokal;
	Timer _timerGatherWellSepLoLokal;
	Timer _timerPropagateCellLoLokal;
	Timer _timerProcessFarField;
	Timer _timerCommunicationHalos;
	Timer _timerHaloGather;
	Timer _timerBusyWaiting;
	Timer _timerFMMcomplete;

	Timer _timerGatherEvalM;
	Timer _timerGatherEvalLM;
	Timer _timerAllreduce_me;
#ifdef ENABLE_MPI

	MPI_Comm _comm;
#endif
	int _overlapComm;

};

} /* namespace bhfmm */

#endif /* UNIFORMPSEUDOPARTICLECONTAINER_H_ */
