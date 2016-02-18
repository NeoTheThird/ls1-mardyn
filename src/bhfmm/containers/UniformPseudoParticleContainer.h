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
#include <vector>
#include <map>

class Domain;
class DomainDecompBase;

namespace bhfmm {

class UniformPseudoParticleContainer: public PseudoParticleContainer {
public:
	UniformPseudoParticleContainer(double domainLength[3], double bBoxMin[3], double bBoxMax[3],
			double LJCellLength[3], unsigned LJSubdivisionFactor, int orderOfExpansions, bool peridic = true);
	~UniformPseudoParticleContainer();

	void clear();
	void build(ParticleContainer* pc);
	void upwardPass(P2MCellProcessor * cp);
	void horizontalPass(VectorizedChargeP2PCellProcessor * cp);
	void downwardPass(L2PCellProcessor *cp);

	// P2M
	void processMultipole(ParticleCell& cell);

	// L2P
	void processFarField(ParticleCell& cell);

	// M2M, M2L, L2L
	void processTree();

	void printTimers();

private:
	LeafNodesContainer* _leafContainer;

	int _wellSep;
	int _maxLevel;
	std::vector<std::vector<MpCell> > _mpCell;
	double _cellLength[3];
	int _globalNumCellsPerDim;
	Domain* _domain;
	int _globalNumCells;
	int* _occVector; // array for MPI allgather

	int _coeffVectorLength;
	double* _coeffVector; // array for MPI allgather
	double* _coeffVector_me;

	bool _periodicBC;

	// M2M
	void CombineMpCell(double *cellWid, int& mpCells, int& curLevel);

	// M2M
	void CombineMpCell_Wigner(double *cellWid, int& mpCells, int& curLevel);

	// M2L
	void GatherWellSepLo(double *cellWid, int mpCells, int& curLevel);

	// M2L
	void GatherWellSepLo_Wigner(double *cellWid, int mpCells, int& curLevel);

	// L2L
	void PropagateCellLo(double *cellWid, int mpCells, int& curLevel);

	// L2L
	void PropagateCellLo_Wigner(double *cellWid, int mpCells, int& curLevel);

	// for parallelization
	void AllReduceMultipoleMoments();
	void AllReduceLocalMoments(int mpCells, int _curLevel);

	// Lookup
	inline const WignerMatrix& M2M_Wigner(const int& idx) const {
		return _M2M_Wigner[idx];
	}

	// Lookup
	inline const WignerMatrix& L2L_Wigner(const int& idx) const {
		return _L2L_Wigner[idx];
	}

	inline double* CosSin_ptr(const int& idx) const {
		return _CosSin + idx*2*(_maxOrd+1);
	}

	std::vector<WignerMatrix> _M2M_Wigner;
	std::vector<WignerMatrix> _L2L_Wigner;
	double* _CosSin;

	std::map<Vector3<int>, RotationParams,  Vector3<int>::compare> _M2L_Wigner;

	// Initializing M2L Wigner matrices
	void processTreeInitM2LWigner();

	void GatherWellSepLoInitM2LWigner(double *cellWid, int mpCells, int& curLevel);


	Timer _timerProcessCells;
	Timer _timerAllreduce;
	Timer _timerCombineMpCell;
	Timer _timerGatherWellSepLo;
	Timer _timerPropagateCellLo;
	Timer _timerProcessFarField;

	Timer _timerGatherEvalM;
	Timer _timerGatherEvalLM;
	Timer _timerAllreduce_me;


};

} /* namespace bhfmm */

#endif /* UNIFORMPSEUDOPARTICLECONTAINER_H_ */