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

class MemoryManager {

public:

	MemoryManager();

	virtual ~MemoryManager();

	utils::DynamicArray<HandlerMoleculeType, true, false>* getMoleculeArray();

	void releaseMoleculeArray(utils::DynamicArray<HandlerMoleculeType, true, false>* memory);

private:

	std::vector<utils::DynamicArray<Molecule, true, false>*> _moleculeArrays;


};

#endif /* MEMORYMANAGER_H_ */
