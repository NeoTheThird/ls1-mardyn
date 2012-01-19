/*
 * MemoryManager.cpp
 *
 * @Date: 30.09.2011
 * @Author: eckhardw
 */

#include "MemoryManager.h"
#include "utils/Logger.h"
#include <malloc.h>
#include <xmmintrin.h>

using namespace std;
using namespace utils;

//#define REUSE_MEMORY

unsigned long dynamicArraySize = 0;

MemoryManager::MemoryManager() {
	_sizeX = 16;
	_sizeY = 16;
	allocateScratchMemory(16, 16);
}

MemoryManager::~MemoryManager() {
	unsigned long memUsage = 0;
	// sizeof _scratchMemory
	memUsage += _sizeX * _sizeY * sizeof(fp_memory_type);

	std::cout << "Number of Molecule arrays=" << _moleculeArrays.size() << endl;

	std::vector<HandlerMoleculeTypeArray*>::iterator it1 =  _moleculeArrays.begin();
	while (it1 != _moleculeArrays.end()) {
		memUsage += (*it1)->capacity() * sizeof (HandlerMoleculeType);
		it1++;
	}

	std::vector<std::pair<size_t, fp_memory_type*> >::iterator it2 = _freeFPMemory.begin();
	while (it2 != _freeFPMemory.end()) {
		memUsage += it2->first * sizeof(fp_memory_type);
		it2++;
	}
	std::cout << "Memory Manager allocated " << (memUsage / 1024) << " kBytes." << endl;

#ifdef REUSE_MEMORY
	size_t size = _moleculeArrays.size();
	for (size_t i = 0; i < size; i++) {
		delete _moleculeArrays.back();
		_moleculeArrays.pop_back();
	}
#endif

	deallocateScratchMemory(_sizeX, _sizeY);
}


HandlerMoleculeTypeArray* MemoryManager::getMoleculeArray() {
#ifdef REUSE_MEMORY
	if (_moleculeArrays.empty()) {
		return new HandlerMoleculeTypeArray(8, 1.2);
	} else {
		HandlerMoleculeTypeArray* last = _moleculeArrays.back();
		_moleculeArrays.pop_back();
		return last;
	}
#else
	return new HandlerMoleculeTypeArray(8);
#endif
}

void MemoryManager::releaseMoleculeArray(HandlerMoleculeTypeArray* memory) {
#ifdef REUSE_MEMORY
	memory->clear();
//	std::cout << "Clear memory of capacity " << memory->capacity() << std::endl;
//	std::cout << "sizeof MoleculeArrays: " << _moleculeArrays.capacity() * sizeof (HandlerMoleculeType) << " Bytes." << endl;
	_moleculeArrays.push_back(memory);
#else
	delete memory;
#endif
}


MemoryManager::fp_memory_type* MemoryManager::getFPMemory(size_t size) {

	size = (size / 4 + 1) * 4;
	assert(size % 4 == 0);

#ifdef REUSE_MEMORY
	if (_freeFPMemory.empty()) {
		fp_memory_type* newMemory = static_cast<fp_memory_type*> (_mm_malloc(size * sizeof(fp_memory_type), 16));
		_fpMemory.push_back(make_pair(size, newMemory));
		return newMemory;
	} else {
		std::pair<size_t, fp_memory_type*> memoryPair = _freeFPMemory.back();
		_freeFPMemory.pop_back();
		if (memoryPair.first < size) {
			_mm_free(memoryPair.second);
			memoryPair.second = static_cast<fp_memory_type*> (_mm_malloc(size * sizeof(fp_memory_type), 16));
			memoryPair.first = size;
		}
		_fpMemory.push_back(memoryPair);
		return memoryPair.second;
	}
#else
	Log::global_log->debug() << "MemoryManager::getFPMemory() size=" << size << std::endl;
	return static_cast<fp_memory_type*> (_mm_malloc(size * sizeof(fp_memory_type), 16));
#endif
}

void MemoryManager::releaseFPMemory(fp_memory_type* memory) {
#ifdef REUSE_MEMORY
	std::pair<size_t, fp_memory_type*> memoryPair = _fpMemory.front();
	_fpMemory.pop_front();
	assert(memoryPair.second == memory);
	_freeFPMemory.push_back(memoryPair);
#else
	_mm_free(memory);
#endif
}


MemoryManager::fp_memory_type** MemoryManager::getScratchMemory(size_t x, size_t y) {
	if(x > _sizeX) {
		deallocateScratchMemory(_sizeX, _sizeY);
		_sizeX = (x / 4 + 1 ) * 4;
		if (y > _sizeY) {
			_sizeY = (y / 4 + 1) *4;
		}

		allocateScratchMemory(_sizeX, _sizeY);
	} else if (y > _sizeY) {
		deallocateScratchMemory(_sizeX, _sizeY);
		_sizeY = (y / 4 + 1 ) * 4;
		allocateScratchMemory(_sizeX, _sizeY);
	}

	assert(_sizeX % 4 == 0);
	assert(_sizeY % 4 == 0);
	return _scratchMemory;
}


void MemoryManager::allocateScratchMemory(size_t x, size_t y) {
	_scratchMemory = new fp_memory_type*[x];

	for (size_t i = 0; i < x; i++) {
		_scratchMemory[i] = static_cast<fp_memory_type*> (_mm_malloc(y * sizeof(fp_memory_type), 16));
	}

	if (Log::global_log != NULL) {
		Log::global_log->debug() << "MemoryManager::allocateScratchMemory() Allocated memory x=" << x << " y=" << y << endl;
	}
}

void MemoryManager::deallocateScratchMemory(size_t x, size_t y) {
	for (size_t i = 0; i < x; i++) {
		_mm_free(_scratchMemory[i]);
	}
	Log::global_log->debug() << "MemoryManager::deallocateScratchMemory() Deleted memory x=" << x << " y=" << y << endl;
	delete[] _scratchMemory;
}
