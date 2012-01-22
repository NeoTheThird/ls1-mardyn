/*
 * MemoryManager.h
 *
 * @Date: 30.09.2011
 * @Author: eckhardw
 */

#ifndef MEMORYMANAGER_H_
#define MEMORYMANAGER_H_

#include "molecules/MoleculeTypes.h"
#include "utils/DynamicArray.h"
#include <vector>
#include <list>

//#define VECTORIZE

#define PROFILE_MEMORY

#ifdef PROFILE_MEMORY
#include "utils/CProcessMemoryInformation.hpp"
#endif

class MemoryManager {

public:

	typedef float fp_memory_type;

	MemoryManager();

	virtual ~MemoryManager();

	HandlerMoleculeTypeArray* getMoleculeArray();

	void releaseMoleculeArray(HandlerMoleculeTypeArray* memory);

	/* === Methods based on the ring buffer used from HybridReorderedLinkedCells === */
	void initRingBuffer(size_t size);

	void deleteRingBuffer();

	HandlerMoleculeType* getMoleculeArray(size_t size);

	void releaseMoleculeArray(HandlerMoleculeType* memory, size_t size);
	/** End RB related methods */

	/**
	 * return aligned memory (aligned to 16, size is a multiple of 4*sizeof(fp_memory_type)
	 */
	fp_memory_type* getFPMemory(size_t size);

	void releaseFPMemory(fp_memory_type* memory);

	/**
	 * return a two-dimensional array of at least size x times y.
	 */
	fp_memory_type** getScratchMemory(size_t x, size_t y);


private:

	HandlerMoleculeType* _ringBuffer;
	size_t _rbSize;
	size_t _rbFirstAllocated;
	size_t _rbFirstFree;

	std::vector<HandlerMoleculeTypeArray*> _moleculeArrays;

	std::list<std::pair<size_t, fp_memory_type*> > _fpMemory;

	std::vector<std::pair<size_t, fp_memory_type*> > _freeFPMemory;

	size_t _sizeX;
	size_t _sizeY;
	fp_memory_type** _scratchMemory;

	void allocateScratchMemory(size_t x, size_t y);

	void deallocateScratchMemory(size_t x, size_t y);

};

#endif /* MEMORYMANAGER_H_ */
