/*
 * PseudoParticleContainer.h
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#ifndef PSEUDOPARTICLECONTAINER_H_
#define PSEUDOPARTICLECONTAINER_H_

#include "bhfmm/pseudoParticles/SHMultipoleParticle.h"
#include "bhfmm/pseudoParticles/SHLocalParticle.h"
#include "particleContainer/ParticleCell.h"
#include "bhfmm/cellProcessors/P2MCellProcessor.h"
#include "bhfmm/cellProcessors/L2PCellProcessor.h"
#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"

class ParticleContainer;

namespace bhfmm {

class MpCell {
public:

	MpCell(int order) :
			occ(0), multipole(order), local(order) {
	}

	int occ;
	bhfmm::SHMultipoleParticle multipole;
	bhfmm::SHLocalParticle local;
};

class PseudoParticleContainer {
public:
	PseudoParticleContainer(int maxOrd) :
			_maxOrd(maxOrd) {
	}
	virtual ~PseudoParticleContainer() {
	}
	
	virtual void build(ParticleContainer* pc) = 0;
	virtual void upwardPass(P2MCellProcessor * cp) = 0;
	virtual void horizontalPass(VectorizedChargeP2PCellProcessor * cp) = 0;
	virtual void downwardPass(L2PCellProcessor *cp) = 0;

	// P2M
	virtual void processMultipole(ParticleCell& cell) = 0;
	// L2P
	virtual void processFarField(ParticleCell& cell) = 0;
	// M2M M2L L2L
	virtual void processTree() = 0;

	virtual void printTimers() = 0;

	virtual void clear() = 0;

protected:
	int _maxOrd;

};

} /* namespace bhfmm */

#endif /* PSEUDOPARTICLECONTAINER_H_ */
