/*
 * NeighbourCommunicationScheme.h
 *
 *  Created on: Sep 29, 2016
 *      Author: seckler
 */

#pragma once

#include "parallel/CommunicationPartner.h"

#include <vector>

class DomainDecompMPIBase;
class Domain;
class FullShell;
class HaloRegion;
class NeighbourCommunicationScheme {
public:
	/**
	 * Specifies the amount of sequential communication steps needed for the communication scheme.
	 * This is also the outer size of DomainDecompMPIBase::_neighbours
	 * @return
	 */
	unsigned int getCommDims() {
		return _commDimms;
	}
	NeighbourCommunicationScheme() = delete;
	NeighbourCommunicationScheme(unsigned int commDimms);
	virtual ~NeighbourCommunicationScheme();

	virtual void prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp) = 0;

	virtual void finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp) = 0;

	virtual void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) = 0;

	void setCoverWholeDomain(unsigned int d, bool covers) {
		_coversWholeDomain[d] = covers;
	}

	virtual void initCommunicationPartners(double cutoffRadius, Domain * domain, DomainDecompMPIBase* domainDecomp) = 0;

	virtual std::vector<int> get3StageNeighbourRanks() = 0;

	virtual std::vector<int> getFullShellNeighbourRanks() {
		std::vector<int> neighbourRanks;
		for (unsigned int i = 0; i < _fullShellNeighbours.size(); i++) {
			neighbourRanks.push_back(_fullShellNeighbours[i].getRank());
		}
		return neighbourRanks;
	}
protected:

	//! vector of neighbours. The first dimension should be of size getCommDims().
	std::vector<std::vector<CommunicationPartner>> _neighbours;

	//! flag, which tells whether a processor covers the whole domain along a dimension
	//! if true, we will use the methods provided by the base class for handling the
	//! respective dimension, instead of packing and unpacking messages to self
	bool _coversWholeDomain[3];

	unsigned int _commDimms;

	//! communication scheme (FullShell, EightShell, ...)
	FullShell* _commScheme;

	//! list of all neighbours (non-squeezed)
	std::vector<CommunicationPartner> _fullShellNeighbours;

};

class DirectNeighbourCommunicationScheme: public NeighbourCommunicationScheme {
public:
	DirectNeighbourCommunicationScheme() :
			NeighbourCommunicationScheme(1) {
	}
	virtual ~DirectNeighbourCommunicationScheme() {
	}
	virtual void initCommunicationPartners(double cutoffRadius, Domain * domain, DomainDecompMPIBase* domainDecomp)
			override;
	virtual std::vector<int> get3StageNeighbourRanks() override {
		std::vector<int> neighbourRanks;
		for (unsigned int i = 0; i < _neighbours[0].size(); i++) {
			if (_neighbours[0][i].isFaceCommunicator()) {
				neighbourRanks.push_back(_neighbours[0][i].getRank());
			}
		}
		return neighbourRanks;
	}

	virtual void prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

	virtual void finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

	virtual void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp);

protected:
	void finalizeExchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* /*domain*/, MessageType /*msgType*/,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp);
	void initExchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* /*domain*/, MessageType msgType,
			bool /*removeRecvDuplicates*/, DomainDecompMPIBase* domainDecomp);
};

class IndirectNeighbourCommunicationScheme: public NeighbourCommunicationScheme {
public:

	IndirectNeighbourCommunicationScheme() :
			NeighbourCommunicationScheme(3) {
	}
	virtual ~IndirectNeighbourCommunicationScheme() {
	}
	void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) override;
	virtual void initCommunicationPartners(double cutoffRadius, Domain * domain, DomainDecompMPIBase* domainDecomp)
			override;
	virtual std::vector<int> get3StageNeighbourRanks() override {
		std::vector<int> neighbourRanks;
		for (unsigned int i = 0; i < _fullShellNeighbours.size(); i++) {
			if (_fullShellNeighbours[i].isFaceCommunicator()) {
				neighbourRanks.push_back(_fullShellNeighbours[i].getRank());
			}
		}
		return neighbourRanks;
	}

	virtual void prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

	virtual void finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

protected:
	void initExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);

	void finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);
	void exchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);
	void convert1StageTo3StageNeighbours(const std::vector<CommunicationPartner>& commPartners,
			std::vector<std::vector<CommunicationPartner>>& neighbours, HaloRegion& ownRegion, double cutoffRadius);

};
