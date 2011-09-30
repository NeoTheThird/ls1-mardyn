/*
 * MemoryManager.cpp
 *
 * @Date: 30.09.2011
 * @Author: eckhardw
 */

#include "MemoryManager.h"

using namespace std;
using namespace utils;

MemoryManager::MemoryManager() {

}

MemoryManager::~MemoryManager() {
#ifdef REUSE_MEMORY
	for (int i = 0; i < _moleculeArrays.size(); i++) {
		DynamicArray<HandlerMoleculeType, true, false>* array = _moleculeArrays.pop_back();
		delete array;
	}
#endif
}


DynamicArray<HandlerMoleculeType, true, false>* MemoryManager::getMoleculeArray() {
#ifdef REUSE_MEMORY
	if (_moleculeArrays.empty()) {
		return new DynamicArray<HandlerMoleculeType, true, false>();
	} else {
		return _moleculeArrays.pop_back();
	}
#else
	return new DynamicArray<HandlerMoleculeType, true, false>();
#endif
}

void MemoryManager::releaseMoleculeArray(DynamicArray<HandlerMoleculeType, true, false>* memory) {
#ifdef REUSE_MEMORY
	_moleculeArrays.push_back(memory);
#else
	delete memory;
#endif
}

